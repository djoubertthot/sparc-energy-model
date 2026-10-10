import numpy as np
import pandas as pd
from pathlib import Path

# Explicit formulation: m_eff = E / (γ c²) with γ ≈ 1
A0 = 1.2e-10
BETA = 1.4
E0 = 1.0e-24 * (3e8)**2  # equivalent to ρ0 = 1e-24
YSTAR = 0.6

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98847e30
KPC = 3.085677581e19

def load_rotmod(filepath):
    data = np.loadtxt(filepath, comments='#')
    n = data.shape[1]
    cols = ['R', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
    return pd.DataFrame(data, columns=cols[:n])

def compute_radial_model(df, Ystar=YSTAR, beta=BETA, e0=E0):
    R_kpc = df['R'].values
    R = R_kpc * KPC
    Vgas = df['Vgas'].values * 1000.0
    Vdisk = df['Vdisk'].values * 1000.0
    Vbul = df['Vbul'].values * 1000.0 if 'Vbul' in df.columns else np.zeros_like(Vgas)

    Vbar2 = Vgas**2 + (Ystar * Vdisk)**2 + (Ystar * Vbul)**2
    Vbar2 = np.maximum(Vbar2, 1.0)
    a_bar = Vbar2 / np.maximum(R, 1e15)

    # Explicit energy density
    E = e0 * (A0 / np.maximum(a_bar, 1e-13))**beta
    rho_eff = E / (C**2)

    dR = np.gradient(R)
    dM = 4 * np.pi * R**2 * rho_eff * dR
    M_cum = np.cumsum(dM)

    Veq2 = G * M_cum / np.maximum(R, 1e15)
    Vmodel = np.sqrt(Vbar2 + np.maximum(Veq2, 0)) / 1000.0

    out = df.copy()
    out['Vbar'] = np.sqrt(Vbar2) / 1000.0
    out['Vmodel'] = Vmodel
    out['residual'] = out['Vobs'] - out['Vmodel']
    return out

def process_directory(rotmod_dir, max_files=None):
    files = list(Path(rotmod_dir).glob('*_rotmod.dat'))
    if max_files:
        files = files[:max_files]
    results = []
    for f in files:
        try:
            df = load_rotmod(f)
            if len(df) < 5:
                continue
            modeled = compute_radial_model(df)
            rms = np.sqrt(np.nanmean(modeled['residual']**2))
            results.append({
                'Galaxy': f.stem.replace('_rotmod', ''),
                'N': len(df),
                'RMS': rms
            })
        except Exception:
            continue
    return pd.DataFrame(results)

if __name__ == '__main__':
    res = process_directory('data/rotmod/')
    print(res['RMS'].describe())
