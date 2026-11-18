"""
dspl_lensing_ratio.py
Metric Synchronization Operator for Double-Source-Plane Lenses (DSPL)
Target: Euclid DR1-Foundation DSPL Catalog (arXiv:2611.20166, 84 candidates)
Author: Michael James Gross
Archive: https://doi.org/10.5281/zenodo.24182514
"""

G_C = 1.0857
DET_M = 1.1787
DELTA = 0.1787


def fractional_einstein_offset(z_lens: float, z_source: float) -> float:
    """Fractional Einstein radius residual predicted by metric synchronization."""
    if z_source <= z_lens:
        raise ValueError("Source redshift must exceed lens redshift.")
    return DELTA * ((z_source - z_lens) / (1.0 + z_source))


def dspl_beta_ratio_shift(z_lens: float, z_source1: float, z_source2: float) -> dict:
    """
    Computes the differential shift in the inner-to-outer Einstein ring ratio
    (theta_E1 / theta_E2) for a double-source-plane lens system.
    """
    if not (z_lens < z_source1 < z_source2):
        raise ValueError("Require z_lens < z_source1 < z_source2.")

    offset_1 = fractional_einstein_offset(z_lens, z_source1)
    offset_2 = fractional_einstein_offset(z_lens, z_source2)

    beta_multiplier = (1.0 + offset_1) / (1.0 + offset_2)
    percent_shift = (beta_multiplier - 1.0) * 100.0

    return {
        "z_lens": z_lens,
        "z_source1": z_source1,
        "z_source2": z_source2,
        "inner_ring_offset_pct": round(offset_1 * 100.0, 4),
        "outer_ring_offset_pct": round(offset_2 * 100.0, 4),
        "beta_ratio_multiplier": round(beta_multiplier, 6),
        "beta_ratio_shift_pct": round(percent_shift, 4),
    }


if __name__ == "__main__":
    sample = dspl_beta_ratio_shift(z_lens=0.5, z_source1=1.2, z_source2=2.5)
    print("Euclid DR1 DSPL Metric Synchronization Benchmark:")
    for k, v in sample.items():
        print(f"  {k}: {v}")
