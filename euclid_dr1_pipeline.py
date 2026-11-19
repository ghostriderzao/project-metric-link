
"""
euclid_dr1_pipeline.py
Euclid DR1-Foundation TAP/ADQL Extraction, Section 15.2/15.3 Filtering,
and ECLIPSE Single-Lens & Double-Source-Plane Lens (DSPL) Ratio Audit
Repository: ghostriderzao/project-metric-link
Archive DOI: https://doi.org/10.5281/zenodo.24182514
Constants: G_c = ±1.0857, |det M| = 1.1787, delta = 0.1787
"""

from pathlib import Path
import numpy as np
import pandas as pd
from astroquery.esa.euclid import Euclid

# Metric Synchronization Constants
G_C = 1.0857
DET_M = 1.1787
DELTA = 0.1787


def build_section15_filtered_adql(limit: int = 10000) -> str:
    """
    Constructs an ADQL query incorporating Euclid DR1-Foundation Explanatory
    Supplement Section 15.2 and 15.3 quality controls:
      - Section 15.2.1 (Slitless Spectroscopy): spe_rank = 0, spe_z_prob >= 0.95
      - Section 15.3.2 / 15.3.4 (T-PHOT & Worst-PSF bugs): Selects A-PHOT aperture
        fluxes over T-PHOT
      - Section 15.3.7 (Duplicate tile-edge sources): Filters spurious_flag = 0
        and deduplicates on object_id in post-processing
      - PHZ precision: Filters 90% credible interval width < 0.10 * (1 + z)
    """
    return f"""
    SELECT TOP {limit}
        mer.object_id,
        mer.ra,
        mer.dec,
        mer.vis_det,
        mer.spurious_flag,
        mer.flux_vis_1fwhm_aper,
        mer.fluxerr_vis_1fwhm_aper,
        mer.flux_y_1fwhm_aper,
        mer.flux_j_1fwhm_aper,
        mer.flux_h_1fwhm_aper,
        mer.kron_radius,
        spe.spe_z,
        spe.spe_z_prob,
        phz.phz_median,
        phz.phz_90_int1,
        phz.phz_90_int2
    FROM euclid_dr1_mer_catalogue AS mer
    LEFT JOIN euclid_dr1_spectro_z_spe_galaxy_candidates AS spe
        ON mer.object_id = spe.object_id
        AND spe.spe_rank = 0
    LEFT JOIN euclid_dr1_phz_photo_z AS phz
        ON mer.object_id = phz.object_id
    WHERE
        mer.vis_det = 1
        AND mer.spurious_flag = 0
        AND mer.flux_vis_1fwhm_aper > 0
        AND (
            (spe.spe_z IS NOT NULL AND spe.spe_z_prob >= 0.95)
            OR (
                phz.phz_median IS NOT NULL
                AND (phz.phz_90_int2 - phz.phz_90_int1) / (1.0 + phz.phz_median) < 0.10
            )
        )
    """


def fetch_and_clean_dr1_tap(limit: int = 10000) -> pd.DataFrame:
    """
    Executes the Section 15.2/15.3 ADQL query against the ESA Euclid Science Archive
    and enforces Section 15.3.7 deduplication and Section 15.2.1 blue-grism flags.
    """
    query = build_section15_filtered_adql(limit=limit)
    job = Euclid.launch_job_async(query)
    df = job.get_results().to_pandas()

    # Section 15.3.7: Deduplicate overlapping MER tile-edge entries by object_id
    df = df.drop_duplicates(subset=["object_id"], keep="first").copy()

    # Select best available redshift (spectroscopic priority over validated PHZ)
    has_spec = df["spe_z"].notna() & (df["spe_z_prob"] >= 0.95)
    df["z_best"] = np.where(has_spec, df["spe_z"], df["phz_median"])
    df["z_source"] = np.where(has_spec, "SPE", "PHZ")

    # Section 15.2.1: Flag redshifts where H-alpha (656.28 nm) falls inside the
    # uncalibrated NISP blue-grism window (926 - 1366 nm -> 0.411 < z < 1.081)
    df["blue_grism_halpha_flag"] = (df["z_source"] == "SPE") & df["z_best"].between(0.411, 1.081)

    return df


def _comoving_distance_mpc(z: np.ndarray, h0: float = 70.0, om: float = 0.3) -> np.ndarray:
    """Vectorized flat-LCDM line-of-sight comoving distance in Mpc."""
    c_kms = 299792.458
    ol = 1.0 - om
    z_arr = np.asarray(z, dtype=float)
    out = np.zeros_like(z_arr)
    for idx, zi in np.ndenumerate(z_arr):
        if zi > 0.0:
            grid = np.linspace(0.0, zi, 2000)
            ez = np.sqrt(om * (1.0 + grid) ** 3 + ol)
            out[idx] = (c_kms / h0) * np.trapezoid(1.0 / ez, grid)
    return out


def _angular_diameter_distance_mpc(
    z1: np.ndarray, z2: np.ndarray, h0: float = 70.0, om: float = 0.3
) -> np.ndarray:
    """Angular diameter distance D_A(z1, z2) in Mpc for a flat universe."""
    dc1 = _comoving_distance_mpc(z1, h0=h0, om=om)
    dc2 = _comoving_distance_mpc(z2, h0=h0, om=om)
    return (dc2 - dc1) / (1.0 + np.asarray(z2, dtype=float))


def evaluate_eclipse_sync_residuals(
    eclipse_df: pd.DataFrame,
    dr1_clean_df: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict]:
    """
    Tests ECLIPSE strong-lens Einstein radius residuals against G_c = 1.0857
    and |det M| = 1.1787. Optionally cross-matches with Section 15.2/15.3 cleaned
    MER/SPE/PHZ catalog rows on 'object_id'.
    """
    df = eclipse_df.copy()
    if dr1_clean_df is not None and "object_id" in df.columns and "object_id" in dr1_clean_df.columns:
        df = pd.merge(df, dr1_clean_df, on="object_id", how="inner", suffixes=("_eclipse", "_mer"))

    valid = (
        df["einstein_radius_eff_median_pdf"].notna()
        & df["log10_stellar_mass_median_pdf"].notna()
        & (df["einstein_radius_eff_median_pdf"] > 0.0)
    )
    df = df.loc[valid].copy()

    theta_obs = df["einstein_radius_eff_median_pdf"].to_numpy()
    log_m_star = df["log10_stellar_mass_median_pdf"].to_numpy()

    # Empirical log-linear baseline anchored at ECLIPSE-A median log10(M*/M_sun) = 11.39
    x = log_m_star - 11.39
    slope, intercept = np.polyfit(x, np.log10(theta_obs), 1)
    theta_pred_std = 10.0 ** (slope * x + intercept)
    theta_pred_gc = theta_pred_std * G_C

    df["theta_pred_std"] = theta_pred_std
    df["theta_pred_gc"] = theta_pred_gc
    df["resid_std_arcsec"] = np.abs(theta_obs - theta_pred_std)
    df["resid_gc_arcsec"] = np.abs(theta_obs - theta_pred_gc)
    df["gc_improves_fit"] = df["resid_gc_arcsec"] < df["resid_std_arcsec"]

    summary = {
        "single_lenses_evaluated": int(len(df)),
        "gc_improves_fit_count": int(df["gc_improves_fit"].sum()),
        "gc_improves_fit_pct": round(float(df["gc_improves_fit"].mean() * 100.0), 2) if len(df) else 0.0,
        "median_resid_std_arcsec": round(float(df["resid_std_arcsec"].median()), 4) if len(df) else np.nan,
        "median_resid_gc_arcsec": round(float(df["resid_gc_arcsec"].median()), 4) if len(df) else np.nan,
    }
    return df, summary


def audit_dspl_einstein_ratios(
    eclipse_df: pd.DataFrame,
    dspl_mass_df: pd.DataFrame,
    dr1_clean_df: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict]:
    """
    Isolates the 76 Double-Source-Plane Lenses (DSPLs) from the ECLIPSE DR1 catalogs
    and tests their observed outer-to-inner Einstein radius ratios (theta_E2 / theta_E1)
    against the standard cosmological distance ratio eta = (D_s1 * D_ds2) / (D_ds1 * D_s2)
    and the G_c = 1.0857 metric synchronization operator.
    """
    key_col = "object_id" if "object_id" in eclipse_df.columns else "id_str"

    if "is_dspl" in eclipse_df.columns:
        dspl_mask = eclipse_df["is_dspl"].astype(bool)
    elif "lens_type" in eclipse_df.columns:
        dspl_mask = eclipse_df["lens_type"].str.upper().str.contains("DSPL", na=False)
    else:
        raise KeyError("eclipse_df must contain either 'is_dspl' or 'lens_type'.")

    dspl_cat = eclipse_df.loc[dspl_mask & eclipse_df["grade"].isin(["A", "B"])].copy()
    dspl = pd.merge(dspl_cat, dspl_mass_df, on=key_col, how="inner", suffixes=("", "_mass"))

    if dr1_clean_df is not None and "object_id" in dr1_clean_df.columns:
        dspl = pd.merge(
            dspl,
            dr1_clean_df[["object_id", "z_best", "z_source", "blue_grism_halpha_flag"]],
            left_on=key_col,
            right_on="object_id",
            how="left",
        )
        if "z_lens" in dspl.columns:
            dspl["z_lens"] = dspl["z_best"].fillna(dspl["z_lens"])
        else:
            dspl["z_lens"] = dspl["z_best"]

    th1 = np.minimum(dspl["theta_e1_median_pdf"].to_numpy(), dspl["theta_e2_median_pdf"].to_numpy())
    th2 = np.maximum(dspl["theta_e1_median_pdf"].to_numpy(), dspl["theta_e2_median_pdf"].to_numpy())
    valid = np.isfinite(th1) & np.isfinite(th2) & (th1 > 0.0)
    dspl = dspl.loc[valid].copy()
    th1, th2 = th1[valid], th2[valid]

    ratio_obs = th2 / th1
    dspl["theta_e1_inner"] = th1
    dspl["theta_e2_outer"] = th2
    dspl["ratio_obs_e2_e1"] = ratio_obs
    dspl["delta_from_gc"] = np.abs(ratio_obs - G_C)
    dspl["delta_from_det_m"] = np.abs(ratio_obs - DET_M)

    summary = {
        "dspl_candidates_isolated": int(len(dspl_cat)),
        "dspl_modeled_pairs_evaluated": int(len(dspl)),
        "median_ratio_obs_e2_e1": round(float(np.median(ratio_obs)), 4) if len(dspl) else np.nan,
        "mean_ratio_obs_e2_e1": round(float(np.mean(ratio_obs)), 4) if len(dspl) else np.nan,
        "within_5pct_of_gc_count": int(np.sum(np.abs(ratio_obs / G_C - 1.0) <= 0.05)),
        "within_5pct_of_det_m_count": int(np.sum(np.abs(ratio_obs / DET_M - 1.0) <= 0.05)),
    }

    z_cols = {"z_lens", "z_source_1", "z_source_2"}
    if z_cols.issubset(dspl.columns):
        z_mask = (
            dspl["z_lens"].notna()
            & dspl["z_source_1"].notna()
            & dspl["z_source_2"].notna()
            & (dspl["z_source_2"] > dspl["z_source_1"])
            & (dspl["z_source_1"] > dspl["z_lens"])
        )
        if z_mask.any():
            zd = dspl.loc[z_mask, "z_lens"].to_numpy()
            zs1 = dspl.loc[z_mask, "z_source_1"].to_numpy()
            zs2 = dspl.loc[z_mask, "z_source_2"].to_numpy()

            d_s1 = _angular_diameter_distance_mpc(np.zeros_like(zs1), zs1)
            d_s2 = _angular_diameter_distance_mpc(np.zeros_like(zs2), zs2)
            d_ds1 = _angular_diameter_distance_mpc(zd, zs1)
            d_ds2 = _angular_diameter_distance_mpc(zd, zs2)

            eta_std = (d_s1 * d_ds2) / (d_ds1 * d_s2)
            eta_gc = eta_std * G_C

            r_sub = dspl.loc[z_mask, "ratio_obs_e2_e1"].to_numpy()
            err_std = np.abs(r_sub - eta_std)
            err_gc = np.abs(r_sub - eta_gc)

            dspl.loc[z_mask, "eta_cosmological_std"] = eta_std
            dspl.loc[z_mask, "eta_cosmological_gc"] = eta_gc
            dspl.loc[z_mask, "dspl_gc_improves_fit"] = err_gc < err_std

            summary["dspl_with_redshifts_count"] = int(z_mask.sum())
            summary["dspl_gc_improves_fit_count"] = int(np.sum(err_gc < err_std))
            summary["dspl_median_err_std"] = round(float(np.median(err_std)), 4)
            summary["dspl_median_err_gc"] = round(float(np.median(err_gc)), 4)

    return dspl, summary


def run_pipeline(
    eclipse_csv_path: str | None = None,
    dspl_csv_path: str | None = None,
    use_live_tap: bool = False,
    tap_limit: int = 10000,
    output_dir: str = "outputs",
) -> None:
    """
    Unified runner for single-lens residual testing and the 76-DSPL ratio audit.
    Exports CSV results when real catalog paths are supplied.
    """
    dr1_clean_df = fetch_and_clean_dr1_tap(limit=tap_limit) if use_live_tap else None

    if eclipse_csv_path and Path(eclipse_csv_path).exists():
        eclipse_df = pd.read_csv(eclipse_csv_path)
        dspl_mass_df = pd.read_csv(dspl_csv_path) if (dspl_csv_path and Path(dspl_csv_path).exists()) else eclipse_df
        mode = "LIVE ECLIPSE CATALOG"
    else:
        mode = "OFFLINE VERIFICATION (SYNTHETIC SCHEMA)"
        eclipse_df = pd.DataFrame({
            "object_id": [1001, 1002, 1003, 1004],
            "grade": ["A", "A", "B", "A"],
            "is_dspl": [True, True, False, False],
            "einstein_radius_eff_median_pdf": [0.85, 0.93, 0.78, 1.12],
            "log10_stellar_mass_median_pdf": [11.39, 11.42, 11.31, 11.55],
        })
        dspl_mass_df = pd.DataFrame({
            "object_id": [1001, 1002],
            "theta_e1_median_pdf": [0.72, 0.85],
            "theta_e2_median_pdf": [1.06, 0.93],
            "z_lens": [0.50, 0.45],
            "z_source_1": [1.20, 1.15],
            "z_source_2": [2.50, 1.85],
        })

    single_df, single_stats = evaluate_eclipse_sync_residuals(eclipse_df, dr1_clean_df=dr1_clean_df)
    dspl_df, dspl_stats = audit_dspl_einstein_ratios(eclipse_df, dspl_mass_df, dr1_clean_df=dr1_clean_df)

    print(f"=== {mode} ===")
    print("Single-Lens Audit:", single_stats)
    print("DSPL Ratio Audit: ", dspl_stats)

    if mode == "LIVE ECLIPSE CATALOG":
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        single_df.to_csv(out_path / "eclipse_single_lens_sync_audit.csv", index=False)
        dspl_df.to_csv(out_path / "eclipse_dspl_ratio_audit.csv", index=False)
        print(f"Saved audit CSVs to {out_path.resolve()}")


if __name__ == "__main__":
    # Run offline self-test by default; pass real ECLIPSE CSV paths and use_live_tap=True for full catalog run
    run_pipeline(use_live_tap=False)
