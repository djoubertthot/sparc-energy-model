"""
Estimation spectrale raffinée.
"""
import numpy as np

def spectral_estimate():
    # Énergie stellaire totale (10 Gyr, 10^10 M☉)
    E_stars = 1.2e54
    E_cmb = 1.2e54
    E_nu = 0.8e54

    # Répartition stellaire
    E_uv = 0.4 * E_stars
    E_nir = 0.3 * E_stars
    E_fir = 0.3 * E_stars

    # Facteurs
    f_uv = 30
    f_nir = 10
    f_fir = 2
    f_cmb = 1.5
    f_nu = 1.0

    E_eff = (E_uv * f_uv + E_nir * f_nir + E_fir * f_fir +
             E_cmb * f_cmb + E_nu * f_nu)

    E_needed = 1.8e57
    factor = E_needed / E_eff

    print(f"Énergie effective : {E_eff:.2e} J")
    print(f"Facteur manquant : {factor:.1f}")
    return factor

if __name__ == "__main__":
    spectral_estimate()
