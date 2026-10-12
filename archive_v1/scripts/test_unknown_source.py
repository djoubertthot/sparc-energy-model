"""
Teste le modèle complet avec source externe inconnue.
"""
import numpy as np
import pandas as pd
from pathlib import Path

A0 = 1.2e-10
BETA = 1.4
YSTAR = 0.6
GLOBAL_NORM = 1e-3
F_UNKNOWN = 80

G = 6.67430e-11
MSUN = 1.98847e30
KPC = 3.085677581e19

def compute(df, f_galaxy=1.0):
    R_kpc = df['R'].values
    R = R_kpc * KPC
    Vgas = df['Vgas'].values * 1000.0
    Vdisk = df['Vdisk'].values * 1000.0
    Vbul = df['Vbul'].values * 1000.0 if 'Vbul' in df.columns else np.zeros_like(Vgas)

    Vbar2 = Vgas**2 + (YSTAR * Vdisk)**2 + (YSTAR * Vbul)**2
    Vbar2 = np.maximum(Vbar2, 1.0)
    a_bar = Vbar2 / np.maximum(R, 1e15)

    E_factor = (A0 / np.maximum(a_bar, 1e-13))**BETA * F_UNKNOWN * f_galaxy
    rho_eff = E_factor / np.max(E_factor)

    dR = np.gradient(R)
    dM = 4 * np.pi * R**2 * rho_eff * dR
    M_cum = np.cumsum(dM) * GLOBAL_NORM * MSUN
    Veq2 = G * M_cum / np.maximum(R, 1e15)
    Vmodel = np.sqrt(Vbar2 + np.maximum(Veq2, 0)) / 1000.0

    residual = df['Vobs'].values - Vmodel
    eVobs = df['eVobs'].values if 'eVobs' in df.columns else np.ones_like(residual)
    eVobs = np.maximum(eVobs, 0.1)
    weights = 1.0 / eVobs**2
    return np.sqrt(np.average(residual**2, weights=weights))

def run(rotmod_dir="data/rotmod/"):
    files = list(Path(rotmod_dir).glob("*_rotmod.dat"))
    if not files:
        print("Aucun fichier rotmod trouvé.")
        return

    wrms_list = []
    for f in files:
        try:
            data = np.loadtxt(f, comments='#')
            cols = ['R', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
            df = pd.DataFrame(data, columns=cols[:data.shape[1]])
            if len(df) < 5:
                continue
            wrms_list.append(compute(df))
        except Exception:
            continue

    if wrms_list:
        print(f"WRMS médian : {np.median(wrms_list):.1f} km/s  (N={len(wrms_list)})")
        pd.DataFrame({"WRMS": wrms_list}).to_csv("results/unknown_source_rms.csv", index=False)

if __name__ == "__main__":
    run()
