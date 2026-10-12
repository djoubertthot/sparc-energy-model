"""Final documented model, built from the formula in docs/THEORY.md.

    E_total(r) = (E_known(r) + E_unknown) * f_amplification(r) * f_accumulation
    m_eff      = E_total / (gamma c^2)                      (docs/FORM_AND_E_MC2.md)
    f_amplification(r) = (a0 / a_bary(r))^beta              (THEORY.md, beta = 1.4)
    f_accumulation     = age proxy x size dependence        (docs/AGE_RD_PHYSICS.md)

Every number comes from the repository documents:
- reference galaxy 1e10 Msun, R = 10 kpc, 10 Gyr        (docs/CUMULATIVE_ENERGY.md)
- stellar radiation 1.2e54 J, trapping x5               (docs/REALISTIC_BUDGET.md)
- CMB 1.2e54 J x1.5, neutrinos 0.8e54 J x1              (scripts/spectral_estimate.py)
- CIB brings the missing factor from 110 to 80          (docs/REALISTIC_BUDGET.md)
- E_known + E_unknown = M_bary c^2 for the reference    (docs/MISSING_ENERGY_SUMMARY.md)
- unknown source uniform, cosmological-like             (docs/THEORY.md, COSMOLOGICAL_BACKGROUND.md)
- Y* = 0.5 for disk and bulge, a0 = 1.2e-10, beta = 1.4 (THEORY.md)

Two normalisations are documented:
- "budget": the literal energy inventory above, no free parameter
- "calibrated": the hybrid form where f_missing is chosen so that the uniform
  part equals the phenomenological rho0 (scripts/test_hybrid_rms.py, MISSING_FACTOR.md)

The equivalent density is integrated in spherical shells, as THEORY.md states.
"""
import numpy as np
import pandas as pd

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98847e30
PC = 3.085677581e16
KPC = 1e3 * PC
YR = 365.25 * 86400

A0 = 1.2e-10
BETA = 1.4
YSTAR = 0.5

# --- documented energy inventory, reference galaxy -------------------------
M_REF = 1e10 * MSUN
R_REF = 10 * KPC
V_REF = 4 / 3 * np.pi * R_REF**3
E_NEEDED = M_REF * C**2                       # 1.8e57 J, m_res = M_bary
F_TRAP_INT = 5.0                              # realistic trapping, internal
E_INT = 1.2e54 * F_TRAP_INT                   # stellar radiation over 10 Gyr
E_CMB = 1.2e54 * 1.5                          # CMB with collective factor 1.5
E_NU = 0.8e54 * 1.0                           # neutrinos, no trapping
F_MISSING = 80.0                              # remaining factor after trapping + CIB
E_KNOWN = E_NEEDED / F_MISSING                # 2.25e55 J
E_CIB = E_KNOWN - E_INT - E_CMB - E_NU        # what the CIB must add to reach 1/80
E_UNKNOWN = E_NEEDED - E_KNOWN                # 1.78e57 J, uniform

EPS_INT = E_INT / (M_REF * C**2)              # accumulated internal energy per unit stellar rest energy
RHO_EXT = (E_CMB + E_NU + E_CIB) / (V_REF * C**2)   # kg/m3, uniform
RHO_UNK = E_UNKNOWN / (V_REF * C**2)                # kg/m3, uniform
RHO_UNIFORM_BUDGET = RHO_EXT + RHO_UNK              # 1.6e-22 kg/m3
T_REF = 10e9 * YR


def load_rotmod(filepath):
    data = np.loadtxt(filepath, comments="#")
    cols = ["R", "Vobs", "eVobs", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
    return pd.DataFrame(data, columns=cols[: data.shape[1]])


def baryonic(df, ystar=YSTAR):
    """V_bar^2 in (m/s)^2 with the SPARC convention V|V| and Y* multiplying the mass."""
    gas = df["Vgas"].values * 1e3
    disk = df["Vdisk"].values * 1e3
    bul = df["Vbul"].values * 1e3 if "Vbul" in df.columns else np.zeros_like(gas)
    return gas * np.abs(gas) + ystar * disk * np.abs(disk) + ystar * bul * np.abs(bul)


def stellar_density(df, rdisk_kpc, ystar=YSTAR):
    """Local stellar mass density from the 3.6 um surface brightness, kg/m3.
    Thin disk of scale height h = 0.1 R_disk (the docs give no thickness)."""
    sb = df["SBdisk"].values + (df["SBbul"].values if "SBbul" in df.columns else 0.0)
    sigma = ystar * sb * MSUN / PC**2
    h = 0.1 * max(rdisk_kpc, 0.1) * KPC
    return sigma / (2 * h)


def compute_final_model(df, rdisk_kpc=2.0, f_acc=1.0, beta=BETA, ystar=YSTAR,
                        norm="budget", rho0=1e-24, include_internal=True, gamma=True):
    """Return a DataFrame with Vbar, Vmodel, residual and the density components.

    norm = "budget"     : uniform part = RHO_UNIFORM_BUDGET (literal inventory)
    norm = "calibrated" : uniform part = rho0 (hybrid, f_missing matched to the calibration)
    The internal known part keeps its documented ratio to the uniform part in both cases.
    """
    R = df["R"].values * KPC
    vbar2 = np.maximum(baryonic(df, ystar), 1.0)
    a_bar = vbar2 / np.maximum(R, 1e15)
    f_amp = (A0 / np.maximum(a_bar, 1e-13)) ** beta

    rho_uniform = RHO_UNIFORM_BUDGET if norm == "budget" else rho0
    scale = rho_uniform / RHO_UNIFORM_BUDGET
    rho_int = EPS_INT * stellar_density(df, rdisk_kpc, ystar) * scale if include_internal else 0.0
    rho_pre = rho_int + rho_uniform

    g = 1.0 / np.sqrt(1.0 - np.minimum(vbar2, 0.99 * C**2) / C**2) if gamma else 1.0
    rho_eff = rho_pre * f_amp * f_acc / g

    # spherical integration including the sphere inside the first point
    shell = 4 * np.pi * R**2 * rho_eff
    M = np.empty_like(R)
    M[0] = 4 / 3 * np.pi * R[0] ** 3 * rho_eff[0]
    M[1:] = M[0] + np.cumsum(0.5 * (shell[1:] + shell[:-1]) * np.diff(R))

    veq2 = G * M / np.maximum(R, 1e15)
    vmodel = np.sqrt(vbar2 + veq2) / 1e3

    out = df.copy()
    out["Vbar"] = np.sqrt(vbar2) / 1e3
    out["Vmodel"] = vmodel
    out["residual"] = out["Vobs"] - vmodel
    out["rho_pre"] = rho_pre if np.ndim(rho_pre) else np.full_like(R, rho_pre)
    out["f_amp"] = f_amp
    out["M_eff_Msun"] = M / MSUN
    return out


def rms(out, weighted=False):
    res = out["residual"].values
    if weighted:
        w = 1.0 / np.maximum(out["eVobs"].values, 0.1) ** 2
        return float(np.sqrt(np.average(res**2, weights=w)))
    return float(np.sqrt(np.mean(res**2)))


if __name__ == "__main__":
    print(f"documented inventory: E_known={E_KNOWN:.2e} J (int {E_INT:.1e}, CMB {E_CMB:.1e}, nu {E_NU:.1e}, CIB {E_CIB:.1e})")
    print(f"E_unknown={E_UNKNOWN:.2e} J ; uniform pre-amplification density = {RHO_UNIFORM_BUDGET:.2e} kg/m3 "
          f"(= {RHO_UNIFORM_BUDGET / 1e-24:.0f} x rho0=1e-24)")
    print(f"internal accumulated energy fraction eps_int = {EPS_INT:.1e} of stellar rest energy")
