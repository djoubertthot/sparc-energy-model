import numpy as np
import pandas as pd
from pathlib import Path

# Parameters
A0 = 1.2e-10
BETA = 1.4
YSTAR = 0.6
GLOBAL_NORM = 1e-3

G = 6.67430e-11
MSUN = 1.98847e30
KPC = 3.085677581e19

def load_rotmod(filepath):
    data = np.loadtxt(filepath, comments='#')
    n = data.shape[1]
    cols = ['R', 'Vobs', 'eVobs', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
    return pd.DataFrame(data, columns=cols[:n])

def compute_radial_model(df, Ystar=YSTAR, norm=GLOBAL_NORM):
    R_kpc = df['R'].values
    R = R_kpc * KPC
    Vgas = df['Vgas'].values * 1000.0
    Vdisk = df['Vdisk'].values * 1000.0
    Vbul = df['Vbul'].values * 1000.0 if 'Vbul' in df.columns else np.zeros_like(Vgas)

    Vbar2 = Vgas**2 + (Ystar * Vdisk)**2 + (Ystar * Vbul)**2
    Vbar2 = np.maximum(Vbar2, 1.0)
    a_bar = Vbar2 / np.maximum(R, 1e15)

    log_local = BETA * np.log(A0 / np.maximum(a_bar, 1e-13))
    rho_eff = np.exp(log_local - np.max(log_local))

    dR = np.gradient(R)
    dM = 4 * np.pi * R**2 * rho_eff * dR
    M_cum = np.cumsum(dM) * norm * MSUN

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
                'RMS': rms,
                'Vobs_last': modeled['Vobs'].iloc[-1],
                'Vmodel_last': modeled['Vmodel'].iloc[-1]
            })
        except Exception:
            continue
    return pd.DataFrame(results)

if __name__ == '__main__':
    res = process_directory('data/rotmod/')
    print(res['RMS'].describe())
    res.to_csv('results/rms_results.csv', index=False)
