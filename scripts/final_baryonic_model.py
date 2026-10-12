"""Final model of tag 1.0: pure baryonic, no energy term (archive_v1/docs/CONSOLIDATED_NOTES.md).

    V^2(R) = f_gas * Vgas|Vgas| + Y* * (Vdisk|Vdisk| + Vbul|Vbul|) + G * M_CGM(<R) / R

- Y* = 0.88 on disk and bulge                      (archive_v1/docs/CONSOLIDATED_NOTES.md)
- gas factor f_gas = 2.7, V_gas_eff = V_gas sqrt(f) (docs/GAS_FACTOR_IN_FORMULA.md)
- E0 = 0, no energy or unknown source              (docs/MODEL_OPTION_A.md)
- CGM: beta profile rho ~ [1 + (r/r_c)^2]^(-3 beta / 2), r_c = 6 kpc, beta = 0.68,
  M_CGM(R) = M_tot * f_beta(R / r_c, beta)          (docs/BETA_INTEGRATED.md)
- M_tot = M x (baryonic mass), M = 3.8               (docs/FULL_RETEST_ALL_SCALES.md)

The documents do not say which baryonic mass M multiplies nor at which radius
M_tot is reached, so three readings are implemented (norm argument):
  "last_std"   M_CGM(<R_last) = M * (0.5 L36 + 1.33 MHI)      standard SPARC baryons
  "last_model" M_CGM(<R_last) = M * (Y* L36 + f_gas 1.33 MHI)  the model's own baryons
  "r25_std"    M_CGM(<25 kpc) = M * (0.5 L36 + 1.33 MHI)       fixed radius, scripts/add_extended_cgm.py
"""
import numpy as np
import pandas as pd

G = 6.67430e-11
MSUN = 1.98847e30
KPC = 3.085677581e19

YSTAR = 0.88
F_GAS = 2.7
M_CGM = 3.8
R_C = 6.0      # kpc
BETA = 0.68


def load_rotmod(filepath):
    data = np.loadtxt(filepath, comments="#")
    cols = ["R", "Vobs", "eVobs", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
    return pd.DataFrame(data, columns=cols[: data.shape[1]])


def baryonic2(df, ystar=YSTAR, f_gas=F_GAS, ybul=None):
    """V_bar^2 in (m/s)^2, SPARC convention V|V|, factors multiplying the mass."""
    ybul = ystar if ybul is None else ybul
    gas = df["Vgas"].values * 1e3
    disk = df["Vdisk"].values * 1e3
    bul = df["Vbul"].values * 1e3 if "Vbul" in df.columns else np.zeros_like(gas)
    return f_gas * gas * np.abs(gas) + ystar * disk * np.abs(disk) + ybul * bul * np.abs(bul)


def beta_enclosed(r_kpc, r_c=R_C, beta=BETA, n=4000):
    """Unnormalised enclosed mass of a beta profile at radii r_kpc (array)."""
    rmax = max(float(np.max(r_kpc)), 1e-3)
    grid = np.linspace(0.0, rmax, n)
    integrand = grid**2 * (1.0 + (grid / r_c) ** 2) ** (-1.5 * beta)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(grid))])
    return np.interp(r_kpc, grid, cum)


def cgm_mass(df, l36_1e9, mhi_1e9, m_factor=M_CGM, r_c=R_C, beta=BETA, norm="last_std",
             ystar=YSTAR, f_gas=F_GAS):
    """Enclosed CGM mass in kg at the curve radii, under one of the three normalisations."""
    r = df["R"].values
    if norm == "last_model":
        m_ref = (ystar * l36_1e9 + f_gas * 1.33 * mhi_1e9) * 1e9 * MSUN
        r_norm = r[-1]
    elif norm == "r25_std":
        m_ref = (0.5 * l36_1e9 + 1.33 * mhi_1e9) * 1e9 * MSUN
        r_norm = 25.0
    else:  # last_std
        m_ref = (0.5 * l36_1e9 + 1.33 * mhi_1e9) * 1e9 * MSUN
        r_norm = r[-1]
    radii = np.append(r, r_norm)
    enc = beta_enclosed(radii, r_c, beta)
    scale = m_factor * m_ref / max(enc[-1], 1e-30)
    return enc[:-1] * scale


def compute(df, l36_1e9, mhi_1e9, ystar=YSTAR, f_gas=F_GAS, m_factor=M_CGM, r_c=R_C, beta=BETA,
            norm="last_std", ybul=None):
    R = df["R"].values * KPC
    vbar2 = baryonic2(df, ystar, f_gas, ybul)
    m_cgm = cgm_mass(df, l36_1e9, mhi_1e9, m_factor, r_c, beta, norm, ystar, f_gas) if m_factor > 0 else 0.0
    v2 = vbar2 + G * m_cgm / np.maximum(R, 1e15)
    vmodel = np.sqrt(np.maximum(v2, 0.0)) / 1e3
    out = df.copy()
    out["Vbar"] = np.sqrt(np.maximum(vbar2, 0.0)) / 1e3
    out["Vmodel"] = vmodel
    out["residual"] = out["Vobs"] - vmodel
    out["M_cgm_Msun"] = np.asarray(m_cgm) / MSUN if m_factor > 0 else 0.0
    return out


def rms(out, weighted=False):
    res = out["residual"].values
    if weighted:
        w = 1.0 / np.maximum(out["eVobs"].values, 0.1) ** 2
        return float(np.sqrt(np.average(res**2, weights=w)))
    return float(np.sqrt(np.mean(res**2)))
