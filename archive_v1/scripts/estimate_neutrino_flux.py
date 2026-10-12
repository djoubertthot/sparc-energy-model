"""
Estimation de l'énergie des neutrinos cosmologiques et autres sources.
"""
import numpy as np

def estimate_neutrino_energy():
    # Densité d'énergie des neutrinos cosmologiques ~ 0.5 × CMB (ordre de grandeur)
    u_cmb = 0.26 * 1.602e-19 * 1e6  # J/m³
    u_nu = 0.5 * u_cmb
    print(f"Densité énergie CMB : {u_cmb:.2e} J/m³")
    print(f"Densité énergie neutrinos : {u_nu:.2e} J/m³")
    print(f"Total CMB + ν : {u_cmb + u_nu:.2e} J/m³")
    return u_cmb + u_nu

if __name__ == "__main__":
    estimate_neutrino_energy()
