"""
Metric Synchronization Transformation Operators (Gc = ±1.0857)
Author: Michael James Gross
Reference: Metric Synchronization Across High-Redshift JWST AGN and Euclid DR1-Foundation
"""

import math

# Core Synchronization Constants
G_C = 1.0857
DET_M = -1.1787  # |det M| = G_c^2 = 1.17874449
ABS_DET_M = 1.1787
DELTA = 0.1787   # |det M| - 1


def rescale_time(delta_t_myr: float) -> float:
    """Temporal interval rescaling: Δt' = G_c * Δt (in Myr)."""
    return G_C * delta_t_myr


def rescale_distance_first_order(d_perp: float) -> float:
    """First-order proper spatial scaling: d'_perp = G_c * d_perp."""
    return G_C * d_perp


def rescale_distance_second_order(d_perp: float) -> float:
    """Second-order determinant coupling: d'_perp = |det M| * d_perp."""
    return ABS_DET_M * d_perp


def hubble_synchronized_distance(d_perp: float, h0_base: float = 70.0) -> tuple[float, float]:
    """Returns synchronized H'_0 = G_c * H_0 and contracted angular diameter distance."""
    h0_sync = G_C * h0_base
    d_sync = d_perp / G_C
    return h0_sync, d_sync


def salpeter_mass_boost(delta_t_myr: float, t_salp_myr: float = 45.0) -> float:
    """
    Computes the multiplicative black hole mass gain M(t') / M(t)
    from temporal synchronization at fixed radiative efficiency (default ε = 0.10).
    """
    efold_unscaled = delta_t_myr / t_salp_myr
    efold_sync = rescale_time(delta_t_myr) / t_salp_myr
    return math.exp(efold_sync - efold_unscaled)


def debias_virial_mass(m_bh_obs: float) -> float:
    """Removes (1 + delta) * G_c = 1.2797 metric projection inflation from virial M_BH."""
    projection_factor = (1.0 + DELTA) * G_C
    return m_bh_obs / projection_factor


def euclid_einstein_radius_offset(z_lens: float, z_source: float) -> float:
    """
    Predicts fractional Einstein radius residual Δθ_E / θ_E for Euclid ECLIPSE strong lenses.
    """
    return DELTA * ((z_source - z_lens) / (1.0 + z_source))


if __name__ == "__main__":
    print("=== JWST & Euclid DR1-Foundation Synchronization Benchmarks ===")
    print(f"450 Myr -> {rescale_time(450.0):.3f} Myr (Salpeter boost: {salpeter_mass_boost(450.0):.4f}x)")
    print(f"900 Myr -> {rescale_time(900.0):.3f} Myr")

    d_napo = 18100.0  # J1602+4228AB unscaled separation in ly
    h0_sync, d_hubble = hubble_synchronized_distance(d_napo)
    print(f"J1602+4228AB 1st-order: {rescale_distance_first_order(d_napo):.2f} ly")
    print(f"J1602+4228AB 2nd-order: {rescale_distance_second_order(d_napo):.2f} ly")
    print(f"J1602+4228AB Hubble-sync (H'_0={h0_sync:.2f}): {d_hubble:.2f} ly")

    print(f"Virial M_BH proper/obs ratio: {debias_virial_mass(1.0):.4f}")
    print(f"Euclid ECLIPSE Δθ_E/θ_E (z_l=0.5, z_s=2.0): +{euclid_einstein_radius_offset(0.5, 2.0)*100:.2f}%")
