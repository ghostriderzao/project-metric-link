"""
Universal Metric Synchronization Core Operators
Author: Michael James Gross
Repository: ghostriderzao/project-metric-link
Primary Archive: https://doi.org/10.5281/zenodo.24182514
"""

import math

# Fundamental Metric Synchronization Constants
G_C = 1.0857
G_C_NEG = -1.0857
G_C_RECIPROCAL = 1.0 / G_C          # ~0.92106
DET_M_ABS = round(G_C ** 2, 4)      # 1.1787 (Second-Order Metric Determinant)
DELTA = round(DET_M_ABS - 1.0, 4)   # 0.1787 (Dimensionless Metric Offset)


def rescale_temporal_interval(delta_t_obs: float) -> float:
    """
    Expands observed past-light-cone evolutionary intervals (Δt' = G_c * Δt).
    Example: 450 Myr (GHZ9, z = 10.4) -> 488.565 Myr.
    """
    return delta_t_obs * G_C


def salpeter_accretion_multiplier(radiative_efficiency: float = 0.10) -> float:
    """
    Computes the Salpeter mass-growth enhancement factor over the rescaled
    temporal interval (exp(0.857) = 2.3561x at ε = 0.10).
    """
    return math.exp(0.0857 / radiative_efficiency)


def debias_virial_mass(mass_obs: float) -> float:
    """
    Corrects the (1 + delta) * G_c = 1.2797 metric projection distortion in
    high-redshift virial black hole mass estimators (M_BH,proper = 0.7814 * M_BH,obs).
    """
    projection_factor = (1.0 + DELTA) * G_C
    return mass_obs / projection_factor


def rescale_distance_first_order(distance_lcdm: float) -> float:
    """
    Applies direct first-order metric scaling (G_c = 1.0857) to baseline distances
    or angular scales.
    """
    return distance_lcdm * G_C


def rescale_distance_second_order(distance_lcdm: float) -> float:
    """
    Applies second-order determinant coupling (|det M| = 1.1787) to spatial scales.
    """
    return distance_lcdm * DET_M_ABS


def euclid_einstein_radius_offset(z_lens: float, z_source: float = 2.0) -> float:
    """
    Computes the predicted differential Einstein radius offset in Euclid ECLIPSE
    strong lenses: Δθ_E / θ_E = 0.1787 * ((z_s - z_l) / (1 + z_s)).
    Returns +0.08935 (+8.94%) at z_l = 0.5, z_s = 2.0.
    """
    if z_source <= z_lens or z_lens <= 0.0:
        return 0.0
    return DELTA * ((z_source - z_lens) / (1.0 + z_source))


if __name__ == "__main__":
    print(f"G_c = ±{G_C} | |det M| = {DET_M_ABS} | delta = {DELTA}")
    print(f"GHZ9 Rescaled Interval (450 Myr): {rescale_temporal_interval(450.0):.3f} Myr")
    print(f"Salpeter Multiplier (eps=0.10): {salpeter_accretion_multiplier(0.10):.4f}x")
    print(f"Virial Mass De-bias Factor: {1.0 / ((1.0 + DELTA) * G_C):.4f}")
    print(f"J1602+4228AB 1st Order (18,100 ly): {rescale_distance_first_order(18100.0):.2f} ly")
    print(f"J1602+4228AB 2nd Order (18,100 ly): {rescale_distance_second_order(18100.0):.2f} ly")
    print(f"Euclid ECLIPSE Offset (z_l=0.5, z_s=2.0): +{euclid_einstein_radius_offset(0.5, 2.0) * 100:.2f}%")
