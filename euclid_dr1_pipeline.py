"""
Euclid DR1-Foundation ADQL TAP Extraction & Anomaly Filtering Pipeline
Author: Michael James Gross
Repository: ghostriderzao/project-metric-link
"""

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.constants import c
from astropy.cosmology import FlatLambdaCDM
from astroquery.utils.tap.core import TapPlus

from sync_operators import (
    G_C,
    rescale_distance_first_order,
    euclid_einstein_radius_offset,
)

H0_LCDM = 67.4 * u.km / u.s / u.Mpc
OM0_LCDM = 0.315
cosmo_lcdm = FlatLambdaCDM(H0=H0_LCDM, Om0=OM0_LCDM)
EUCLID_TAP_URL = "https://eas.esac.esa.int/tap-server/tap"


def fetch_clean_dr1_lenses(
    mer_table: str = "dr1.mer_catalogue_deep_survey",
    ra_deg: float | None = None,
    dec_deg: float | None = None,
    radius_deg: float = 0.5,
    row_limit: int = 500,
) -> pd.DataFrame:
    """
    Queries Euclid DR1-Foundation MER + SPE tables while stripping documented
    Section 15.3 release anomalies (Skypatch 40 astrometric warp, spurious
    detections, corrupted detector masks, and MER tile duplicates).
    """
    tap = TapPlus(url=EUCLID_TAP_URL)

    spatial_clause = ""
    if ra_deg is not None and dec_deg is not None:
        spatial_clause = f"""
        AND 1 = CONTAINS(
            POINT('ICRS', mer.right_ascension, mer.declination),
            CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_deg})
        )
        """

    adql_query = f"""
    SELECT TOP {row_limit}
        mer.object_id,
        mer.right_ascension AS ra,
        mer.declination AS dec,
        mer.tile_index,
        mer.spurious_prob,
        mer.det_quality_flag,
        mer.semimajor_axis,
        mer.ellipticity,
        ESDC_NULLIF(mer.flux_vis_sersic, 'NaN') AS flux_vis_sersic,
        spe.spe_z AS z_lens,
        spe.spe_z_err AS z_lens_err,
        spe.spe_velocity_dispersion AS sigma_v_kms,
        spe.spe_velocity_dispersion_err AS sigma_v_err_kms
    FROM {mer_table} AS mer
    INNER JOIN dr1.spectro_z_spe_galaxy_candidates AS spe
        ON mer.object_id = spe.object_id
    WHERE
        mer.spurious_prob < 0.7
        AND mer.det_quality_flag = 0
        AND mer.flag_vis = 0
        AND mer.flag_y = 0
        AND mer.flag_j = 0
        AND mer.flag_h = 0
        AND mer.skypatch_id != 40
        AND spe.spe_z BETWEEN 0.1 AND 1.5
        AND (spe.spe_z_err / (1.0 + spe.spe_z)) < 0.01
        AND spe.spe_velocity_dispersion > 150.0
        {spatial_clause}
    ORDER BY mer.object_id ASC
    """

    job = tap.launch_job_async(adql_query)
    results = job.get_results().to_pandas()
    return results.drop_duplicates(subset=["object_id"], keep="first").reset_index(drop=True)


def evaluate_dr1_lens_residuals(
    df: pd.DataFrame, z_source_assumed: float = 2.0
) -> pd.DataFrame:
    """
    Computes baseline LCDM SIS Einstein radii vs. G_c and euclid_einstein_radius_offset
    predictions on cleaned DR1-Foundation catalog rows.
    """
    z_l = df["z_lens"].to_numpy()
    sigma_v = (df["sigma_v_kms"].to_numpy() * u.km / u.s).to(u.m / u.s)

    d_s = cosmo_lcdm.angular_diameter_distance(z_source_assumed)
    d_ls = cosmo_lcdm.angular_diameter_distance_z1z2(z_l, z_source_assumed)

    theta_e_rad = 4.0 * np.pi * (sigma_v / c) ** 2 * (d_ls / d_s)
    theta_e_lcdm = (theta_e_rad.decompose().value * u.rad).to(u.arcsec).value

    frac_offsets = np.array([euclid_einstein_radius_offset(zl, z_source_assumed) for zl in z_l])
    theta_e_delta_model = theta_e_lcdm * (1.0 + frac_offsets)
    theta_e_gc_direct = rescale_distance_first_order(theta_e_lcdm)

    out = df.copy()
    out["z_source_ref"] = z_source_assumed
    out["theta_E_lcdm_arcsec"] = np.round(theta_e_lcdm, 5)
    out["frac_offset_delta"] = np.round(frac_offsets, 5)
    out["theta_E_delta_sync_arcsec"] = np.round(theta_e_delta_model, 5)
    out["theta_E_Gc_direct_arcsec"] = np.round(theta_e_gc_direct, 5)
    return out


if __name__ == "__main__":
    print("Querying Euclid DR1-Foundation TAP archive (EDF-N cutout)...")
    clean_df = fetch_clean_dr1_lenses(
        mer_table="dr1.mer_catalogue_deep_survey",
        ra_deg=269.73,
        dec_deg=66.02,
        radius_deg=0.5,
        row_limit=250,
    )
    residual_table = evaluate_dr1_lens_residuals(clean_df, z_source_assumed=2.0)
    residual_table.to_csv("euclid_dr1_metric_sync_residuals.csv", index=False)
    print(f"Clean rows extracted: {len(residual_table)}")
