"""
Estimation du flux nécessaire pour m_res ~ M_bary
"""
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98847e30
KPC = 3.085677581e19
YR = 365.25 * 24 * 3600

def estimate_flux(M_msun=1e10, R_kpc=10, T_gyr=10, eta_eps=0.01):
    M = M_msun * MSUN
    R = R_kpc * KPC
    T = T_gyr * 1e9 * YR
    S = 4 * np.pi * R**2

    # m_res = M  ⇒  E_res = M c²
    # E_res ≈ eta_eps * (G M / R) * (Phi * S * T) / c²
    # ⇒ Phi ≈ (M c^4 * R) / (eta_eps * G M * S * T)
    # Simplification : Phi ≈ (c^4 * R) / (eta_eps * G * S * T)

    Phi = (C**4 * R) / (eta_eps * G * S * T)
    print(f"Flux nécessaire : {Phi:.2e} W/m²")
    print(f"Comparaison CMB (~3e-6 W/m²) : facteur {Phi / 3e-6:.2e}")
    return Phi

if __name__ == "__main__":
    estimate_flux()
    estimate_flux(eta_eps=0.1)
    estimate_flux(eta_eps=1.0)
