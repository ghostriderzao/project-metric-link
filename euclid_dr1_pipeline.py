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
        raise KeyError("eclipse_df must contain either 'is_dspl' or 'lens_type
