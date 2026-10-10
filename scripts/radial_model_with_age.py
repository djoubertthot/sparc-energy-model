import numpy as np
import pandas as pd
from pathlib import Path

# Paramètres optimaux
A0 = 1.2e-10
BETA = 1.4
YSTAR = 0.6
RHO0_FACTOR = 1.0
GLOBAL_NORM = 1e-3

G = 6.67430e-11
MSUN = 1.98847e30
KPC = 3.085677581e19

def load_rotmod(filepath):
    data = np.loadtxt(filepath, comments='#')
    n = data.shape[1]
    cols = ['R', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
    return pd.DataFrame(data, columns=cols[:n])

def compute_radial_model(df, f_age=1.0, Ystar=YSTAR, norm=GLOBAL_NORM, weighted=True):
    R_kpc = df['R'].values
    R = R_kpc * KPC
    Vgas = df['Vgas'].values * 1000.0
    Vdisk = df['Vdisk'].values * 1000.0
    Vbul = df['Vbul'].values * 1000.0 if 'Vbul' in df.columns else np.zeros_like(Vgas)

    Vbar2 = Vgas**2 + (Ystar * Vdisk)**2 + (Ystar * Vbul)**2
    Vbar2 = np.maximum(Vbar2, 1.0)
    a_bar = Vbar2 / np.maximum(R, 1e15)

    log_local = BETA * np.log(A0 / np.maximum(a_bar, 1e-13))
    rho_eff = np.exp(log_local - np.max(log_local)) * f_age * RHO0_FACTOR

    dR = np.gradient(R)
    dM = 4 * np.pi * R**2 * rho_eff * dR
    M_cum = np.cumsum(dM) * norm * MSUN

    Veq2 = G * M_cum / np.maximum(R, 1e15)
    Vmodel = np.sqrt(Vbar2 + np.maximum(Veq2, 0)) / 1000.0

    residual = df['Vobs'].values - Vmodel
    if weighted and 'eVobs' in df.columns:
        eVobs = np.maximum(df['eVobs'].values, 0.1)
        weights = 1.0 / eVobs**2
        rms = np.sqrt(np.average(residual**2, weights=weights))
    else:
        rms = np.sqrt(np.nanmean(residual**2))

    out = df.copy()
    out['Vbar'] = np.sqrt(Vbar2) / 1000.0
    out['Vmodel'] = Vmodel
    out['residual'] = residual
    out['rms'] = rms
    return out

def process_directory(rotmod_dir, weighted=True):
    files = list(Path(rotmod_dir).glob('*_rotmod.dat'))
    results = []
    for f in files:
        try:
            name = f.stem.replace('_rotmod', '')
            df = load_rotmod(f)
            if len(df) < 5:
                continue
            modeled = compute_radial_model(df, weighted=weighted)
            rms = modeled['rms'].iloc[0]
            results.append({
                'Galaxy': name,
                'N': len(df),
                'RMS': rms,
                'Weighted': weighted
            })
        except Exception:
            continue
    return pd.DataFrame(results)

if __name__ == '__main__':
    res = process_directory("data/rotmod/", weighted=True)
    print(res["RMS"].describe())
    res.to_csv("results/rms_weighted.csv", index=False)
