# -*- coding: utf-8 -*-
"""
Is a0 tied to the expansion or constant?  Test on rotation curves at z = 0.5 - 4.6.

Two hypotheses from results/a0_evolution.md:
  H_A : a0(z) = a0(0) * H(z)/H0      (a0 ~ c H / 2 pi, expansion-tied)
  H_B : a0(z) = a0(0)                (a0 ~ c sqrt(Lambda/3) / 2 pi, vacuum-tied)
plus the empirical law of MUSE-DARK III (Ciocan et al. 2026): a0(z) = 1.0e-10 + 1.59e-10 z.

For every galaxy with a measured total acceleration g_obs = V_c^2/R and a baryonic acceleration
g_bar at the same radius, the law g_obs = nu(g_bar/a0) g_bar is fitted for r = a0(z)/a0(0).
Samples (data/external/highz/):
  RC100   Nestor Shachar et al. 2023, ApJ 944, 78 : 100 massive SFGs z = 0.6-2.5, f_DM(<Re), V_c(Re), Re (Table 3)
  MSA-3D  JWST/NIRSpec 2026 (arXiv:2606.27853) : 30 SFGs z = 0.5-1.7, log M* 9-11, f_DM(Re), V_rot(Re), sigma0, Re
  KMOS3D  Ubler et al. 2017, ApJ 842, 121 : 135 SFGs z = 0.6-2.6, v_circ,max and M_bar (sizes from van der Wel 2014)
  KROSS   Harrison et al. 2017, MNRAS 467, 1965 : z ~ 0.9, v_C at 2 R_1/2, M* (gas from Tacconi 2018)
  z~4.5   Roman-Oliveira et al. 2024 (A&A 684, A118; curves on Zenodo 10707348), Rizzo et al. 2021, Fraternali et al. 2021
Published constraints quoted in the report: MUSE-DARK III, Puech et al. 2010, Milgrom 2017, Lelli et al. 2023.
Outputs: results/a0_highz.md, results/a0_highz.csv, results/a0_highz.png
"""
import os
import re
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import i0e, i1e, k0e, k1e

warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
D = "data/external/highz"

A0 = 1.2e-10                      # m/s^2
G_KPC = 4.301e-6                  # kpc (km/s)^2 / Msun
KMS2_KPC_TO_SI = 1e6 / 3.0857e19  # (km/s)^2/kpc -> m/s^2
OM, H0 = 0.3, 70.0
SIG_INT = 0.08                    # intrinsic scatter allowed on log g_obs at fixed g_bar (dex)
LOGR_GRID = np.arange(-1.0, 1.41, 0.01)


def Hz(z):
    return np.sqrt(OM * (1 + z) ** 3 + 1 - OM)


def r_muse(z):
    return (1.0 + 1.59 * z) / 1.2


def nu(y):
    y = np.maximum(y, 1e-12)
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))


def logg_pred(logg_bar, r):
    gb = 10 ** logg_bar
    return np.log10(nu(gb / (r * A0)) * gb)


def t_universe(z):
    from scipy.integrate import quad
    f = lambda zz: 1.0 / ((1 + zz) * Hz(zz))  # noqa: E731
    return 977.8 / H0 * quad(f, z, np.inf)[0]  # Gyr


# ---------------- mass models (kpc, Msun, km/s) ----------------
def v2_expdisc(R, M, Rd):
    y = np.asarray(R, float) / (2 * Rd)
    return 2 * G_KPC * M / Rd * y ** 2 * (i0e(y) * k0e(y) - i1e(y) * k1e(y))


def v2_hernquist(R, M, Re):
    a = Re / 1.8153
    return G_KPC * M * R / (R + a) ** 2


def g_si(v2, R):
    return v2 / R * KMS2_KPC_TO_SI


# ---------------- generic likelihood on the law ----------------
def chi2_table(df):
    """df columns: z, logg_obs, s_obs, logg_bar, s_bar (independent errors, dex).
    Returns chi2[galaxy, logr grid], and chi2 for H_A, H_B, MUSE per galaxy."""
    lb, lo = df.logg_bar.values[:, None], df.logg_obs.values[:, None]
    so, sb = df.s_obs.values[:, None], df.s_bar.values[:, None]
    R = 10 ** LOGR_GRID[None, :]
    pred = logg_pred(lb, R)
    slope = (logg_pred(lb + 0.05, R) - logg_pred(lb - 0.05, R)) / 0.1
    var = so ** 2 + (slope * sb) ** 2 + SIG_INT ** 2
    chi = (lo - pred) ** 2 / var
    out = {}
    for name, rr in (("B", np.ones(len(df))), ("A", Hz(df.z.values)), ("M", r_muse(df.z.values))):
        rr = rr[:, None]
        p = logg_pred(lb, rr)
        sl = (logg_pred(lb + 0.05, rr) - logg_pred(lb - 0.05, rr)) / 0.1
        out[name] = ((lo - p) ** 2 / (so ** 2 + (sl * sb) ** 2 + SIG_INT ** 2))[:, 0]
    return chi, out


def fit_r(chi):
    c = chi.sum(0)
    i = int(np.argmin(c))
    ok = c <= c[i] + 1.0
    lo, hi = LOGR_GRID[ok].min(), LOGR_GRID[ok].max()
    return LOGR_GRID[i], lo, hi, c[i]


def summarize(df, label, zbins):
    rows = []
    for zlo, zhi in zbins:
        s = df[(df.z >= zlo) & (df.z < zhi)]
        if len(s) < 3:
            continue
        chi, out = chi2_table(s)
        lr, lo, hi, cmin = fit_r(chi)
        # systematics: M_bar +/- 0.2 dex (g_bar shifts, g_obs unchanged), V +/- 10 % (g_obs shifts 0.083 dex)
        sysm = []
        for dm in (+0.2, -0.2):
            s2 = s.copy(); s2["logg_bar"] = s2.logg_bar + dm
            sysm.append(fit_r(chi2_table(s2)[0])[0])
        sysv = []
        for dv in (+0.083, -0.083):
            s2 = s.copy(); s2["logg_obs"] = s2.logg_obs + dv
            sysv.append(fit_r(chi2_table(s2)[0])[0])
        zmed = s.z.median()
        # errors rescaled to chi2_min / N = 1 (the quoted errors are conservative): interval from Delta chi2 = chi2_min / N
        c = chi.sum(0); scale = max(cmin / len(s), 1e-3)
        okr = c <= cmin + scale
        lo_r, hi_r = LOGR_GRID[okr].min(), LOGR_GRID[okr].max()
        # uniform shift of the baryonic masses that would make H_A the best fit
        dms = np.arange(-0.6, 0.61, 0.02)
        rf = []
        for dm in dms:
            s2 = s.copy(); s2["logg_bar"] = s2.logg_bar + dm
            rf.append(fit_r(chi2_table(s2)[0])[0])
        rf = np.array(rf)
        k = int(np.argmin(np.abs(rf - np.log10(Hz(zmed)))))
        dm_A = dms[k] if np.abs(rf[k] - np.log10(Hz(zmed))) < 0.1 else np.nan
        rows.append(dict(sample=label, zlo=zlo, zhi=zhi, N=len(s), z_med=zmed, logr_lo_r=lo_r, logr_hi_r=hi_r, dM_for_A=dm_A,
                         chi2A_minus_B_rescaled=(out["A"].sum() - out["B"].sum()) / scale,
                         x_med=np.median(10 ** s.logg_bar / A0), nu_med=np.median(10 ** (s.logg_obs - s.logg_bar)),
                         logr=lr, logr_lo=lo, logr_hi=hi, r=10 ** lr, r_lo=10 ** lo, r_hi=10 ** hi,
                         rA=Hz(zmed), rM=r_muse(zmed),
                         chi2_B=out["B"].sum(), chi2_A=out["A"].sum(), chi2_M=out["M"].sum(), chi2_min=cmin,
                         sys_M_plus=sysm[0] - lr, sys_M_minus=sysm[1] - lr, sys_V_plus=sysv[0] - lr, sys_V_minus=sysv[1] - lr))
    return pd.DataFrame(rows)


# =====================================================================
# 1. RC100
# =====================================================================
def load_rc100():
    t = pd.read_csv(f"{D}/rc100_nestorshachar2023_table3.csv")
    f = t.fDM.clip(0.0, 0.95)
    ef = t.e_fDM.clip(0.02, None)
    gobs = g_si(t.Vc_Re ** 2, t.Re_kpc)
    s_obs = np.sqrt((2 * t.e_Vc_Re / t.Vc_Re / np.log(10)) ** 2 + (t.e_Re_kpc / t.Re_kpc / np.log(10)) ** 2)
    # authors' baryonic acceleration at Re: (1 - fDM) g_obs ; error from fDM (independent of the V, Re errors)
    gbar_auth = (1 - f) * gobs
    s_nu = ef / (1 - f) / np.log(10)
    # cross-check: mass model with their log M_bar, Re, M_bulge (Hernquist, Re_b = 1 kpc) and exponential disc
    Mb, Mbul = 10 ** t.logMbar, 10 ** t.logMbulge
    Mbul = np.minimum(Mbul, 0.9 * Mb)
    v2 = v2_expdisc(t.Re_kpc, Mb - Mbul, t.Re_kpc / 1.678) + v2_hernquist(t.Re_kpc, Mbul, 1.0)
    gbar_mm = g_si(v2, t.Re_kpc)
    d = pd.DataFrame(dict(sample="RC100", name=t.galaxy, z=t.z, R_kpc=t.Re_kpc,
                          logg_obs=np.log10(gobs), s_obs=s_obs,
                          logg_bar=np.log10(gbar_auth), s_bar=s_nu, fDM=f, e_fDM=ef,
                          logg_bar_mm=np.log10(gbar_mm)))
    # the fit is on nu = g_obs/g_bar: g_obs enters both sides, so the relevant observable error is on nu;
    # emulate by giving g_obs its error and g_bar the quadrature of both (conservative)
    d["s_bar"] = np.sqrt(s_nu ** 2 + s_obs ** 2)
    return d


# =====================================================================
# 2. MSA-3D (JWST)
# =====================================================================
def load_msa3d():
    t = pd.read_csv(f"{D}/msa3d_2026.csv")
    t["ID"] = t.ID.astype(str).str.replace(r"\D", "", regex=True)
    t = t.groupby("ID", as_index=False).first()
    t = t.dropna(subset=["z", "Vrot_Re", "Re_disk", "fDM"])
    f = t.fDM.clip(0.0, 0.95)
    ef = np.nanmean(np.c_[t.fDM_up, t.fDM_lo], axis=1)
    ef = np.where(np.isfinite(ef), ef, 0.15)
    ef = np.clip(ef, 0.03, None)
    vc2 = t.Vrot_Re ** 2 + 3.36 * t.sigma0 ** 2          # pressure support (Burkert et al. 2010) at Re = 1.678 Rd
    gobs = g_si(vc2, t.Re_disk)
    s_obs = np.sqrt((2 * 0.05 / np.log(10)) ** 2 + (0.15 / np.log(10)) ** 2)  # 10 % on V_c, 15 % on Re
    s_nu = ef / (1 - f) / np.log(10)
    d = pd.DataFrame(dict(sample="MSA-3D", name=t.ID, z=t.z, R_kpc=t.Re_disk,
                          logg_obs=np.log10(gobs), s_obs=s_obs,
                          logg_bar=np.log10((1 - f) * gobs), s_bar=np.sqrt(s_nu ** 2 + s_obs ** 2),
                          fDM=f, e_fDM=ef, logMstar=t.logMstar, Vrot=t.Vrot_Re, sigma0=t.sigma0))
    return d


# =====================================================================
# 3. KMOS3D (Ubler 2017): v_circ,max at ~2.2 Rd, M_bar; sizes from van der Wel et al. 2014 (late types)
# =====================================================================
def re_vdw(z, logM):
    zz = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75])
    la = np.array([0.86, 0.78, 0.70, 0.65, 0.55, 0.51])
    al = np.array([0.25, 0.22, 0.22, 0.23, 0.22, 0.18])
    return 10 ** (np.interp(z, zz, la) + np.interp(z, zz, al) * (logM - np.log10(5e10)))


def load_kmos3d():
    t = pd.read_csv(f"{D}/ubler2017_kmos3d_tfr.csv")
    Re = re_vdw(t.z, t["logM*"])
    Rd = Re / 1.678
    R = 2.2 * Rd
    Mb = 10 ** t.logMb
    gbar = g_si(v2_expdisc(R, Mb, Rd), R)
    gobs = g_si(t.Vcirc ** 2, R)
    d = pd.DataFrame(dict(sample="KMOS3D", name=t.Seq.astype(str), z=t.z, R_kpc=R,
                          logg_obs=np.log10(gobs), s_obs=np.sqrt((2 * 0.1 / np.log(10)) ** 2 + (0.16 / np.log(10)) ** 2),
                          logg_bar=np.log10(gbar), s_bar=np.sqrt(0.2 ** 2 + (0.16) ** 2),  # M_bar 0.2 dex, size scatter 0.16 dex
                          logMb=t.logMb, Vcirc=t.Vcirc, Re_kpc=Re))
    return d, t


# =====================================================================
# 4. KROSS (Harrison 2017): v_C at 2 R_1/2, M*, gas from Tacconi et al. 2018
# =====================================================================
def mu_mol_tacconi(z, logM, sfr):
    tz = np.array([t_universe(zz) for zz in np.atleast_1d(z)])
    logsfr_ms = (0.84 - 0.026 * tz) * logM - (6.51 - 0.11 * tz)        # Speagle et al. 2014
    dms = np.clip(np.log10(np.maximum(sfr, 1e-3)) - logsfr_ms, -1, 1)
    return 10 ** (0.12 - 3.62 * (np.log10(1 + z) - 0.66) ** 2 + 0.53 * dms - 0.35 * (logM - 10.7))


def load_kross():
    t = pd.read_csv(f"{D}/harrison2017_kross.csv")
    sel = (t.Qual == 1) & (t.fAGN == 0) & (t.fIRR == 0) & (t.thetaim > 25) & (t.f_thetaim == 0) & (t.fEXTRAP == 0)
    sel &= t.VC.notna() & (t.VC > 0) & t.sigma0.notna() & (t.Rim > 0) & t.KinType.astype(str).str.contains("RT")
    sel &= ((t.e_VC + t.E_VC) / 2 / t.VC < 0.3)
    t = t[sel].copy()
    logM = np.log10(t.Mass)
    Re = t.Rim
    Rd = Re / 1.678
    R = 2.0 * Re
    mu = mu_mol_tacconi(t.z.values, logM.values, t.SFR.values)
    Mb = t.Mass * (1 + mu)
    gbar = g_si(v2_expdisc(R, Mb, Rd), R)
    vc2 = t.VC ** 2 + 2 * t.sigma0 ** 2 * (R / Rd)
    gobs = g_si(vc2, R)
    gobs_nop = g_si(t.VC ** 2, R)
    ev = (t.e_VC + t.E_VC) / 2 / t.VC
    d = pd.DataFrame(dict(sample="KROSS", name=t.Name.astype(str).str.strip(), z=t.z, R_kpc=R,
                          logg_obs=np.log10(gobs), s_obs=np.sqrt((2 * np.maximum(ev, 0.1) / np.log(10)) ** 2 + (t.e_Rim / t.Rim / np.log(10)) ** 2),
                          logg_bar=np.log10(gbar), s_bar=np.sqrt(0.25 ** 2 + 0.1 ** 2),
                          logg_obs_nopress=np.log10(gobs_nop), logMstar=logM, mu_mol=mu, VC=t.VC, sigma0=t.sigma0))
    return d


# =====================================================================
# 5. z ~ 4.5 discs
# =====================================================================
def load_z45():
    rows = []
    # Roman-Oliveira et al. 2024: stars Hernquist (Re_*), gas exponential (Re_bar); outermost point of the Zenodo curves
    RO = [  # name, z, logM*, e+, e-, logMgas, e+, e-, Re_star, Re_bar, R_out, V_out, e_V
        ("BRI1335-0417", 4.41, 10.4, 0.1, 0.1, 9.8, 0.3, 0.5, 2.3, 2.9, 4.604, 162.2, 25.8),
        ("J081740", 4.54, 10.6, 0.2, 0.2, 10.5, 0.3, 0.5, 2.4, 2.8, 3.15, 227.2, 60.4),
        ("SGP38326-1", 4.28, 11.0, 0.2, 0.2, 11.3, 0.1, 0.3, 2.2, 2.8, 3.983, 530.6, 67.5),
        ("SGP38326-2", 4.28, 10.3, 0.4, 0.6, 11.1, 0.2, 0.3, 1.7, 2.0, 2.042, 419.2, 53.8),
    ]
    for n, z, lms, ep, em, lmg, gp, gm, res, reb, R, V, eV in RO:
        Ms, Mg = 10 ** lms, 10 ** lmg
        v2 = v2_hernquist(R, Ms, res) + v2_expdisc(R, Mg, reb / 1.678)
        sb = np.sqrt(((Ms * (ep + em) / 2) ** 2 + (Mg * (gp + gm) / 2) ** 2)) / (Ms + Mg)  # dex, mass-weighted
        rows.append(dict(sample="z4.5", name=n, z=z, R_kpc=R, logg_obs=np.log10(g_si(V ** 2, R)),
                         s_obs=np.sqrt((2 * eV / V / np.log(10)) ** 2 + (0.1 / np.log(10)) ** 2),
                         logg_bar=np.log10(g_si(v2, R)), s_bar=np.sqrt(sb ** 2 + 0.1 ** 2), group="curve"))
    # Rizzo et al. 2021 (lensed DSFGs, flat V, R_out not tabulated: 2.5 kpc assumed, +/- factor 1.5) and
    # Fraternali et al. 2021 (R_out = 3 scale lengths assumed): stars + gas as one exponential disc of scale length Rd
    RZ = [  # name, z, M*(1e10), Mgas(1e10), Rd_gas (kpc), V, eV, R_out
        ("SPT0113-46", 4.23, 6.0, 4.3, 1.5, 358, 3, 2.5), ("SPT0345-47", 4.29, 2.3, 1.7, 0.76, 280, 25, 2.5),
        ("SPT0441-46", 4.48, 1.8, 1.4, 0.63, 342, 4, 2.5), ("SPT2146-55", 4.57, 1.0, 1.2, 1.4, 194, 7, 2.5),
        ("SPT2132-58", 4.77, 1.9, 1.9, 2.1, 196, 14, 2.5),
        ("AzTEC/C159", 4.567, 4.5, 15.0, 1.3, 506, 110, 3.9), ("J1000+0234", 4.539, 4.0, 12.5, 1.5, 538, 40, 4.5),
    ]
    for n, z, Ms, Mg, Rd, V, eV, R in RZ:
        M = (Ms + Mg) * 1e10
        v2 = v2_expdisc(R, M, Rd)
        rows.append(dict(sample="z4.5", name=n, z=z, R_kpc=R, logg_obs=np.log10(g_si(V ** 2, R)),
                         s_obs=np.sqrt((2 * eV / V / np.log(10)) ** 2 + 0.18 ** 2),   # R_out uncertain by x1.5
                         logg_bar=np.log10(g_si(v2, R)), s_bar=0.3, group="assumed"))
    return pd.DataFrame(rows)


# =====================================================================
def main():
    rc = load_rc100()
    msa = load_msa3d()
    km, kmraw = load_kmos3d()
    kr = load_kross()
    z45 = load_z45()
    allg = pd.concat([rc, msa, km, kr, z45], ignore_index=True)
    allg["x"] = 10 ** allg.logg_bar / A0
    allg["nu_obs"] = 10 ** (allg.logg_obs - allg.logg_bar)
    allg["rA"] = Hz(allg.z)
    allg.to_csv("results/a0_highz.csv", index=False)

    L = []
    L.append("# a0 at high redshift: expansion-tied or constant? (2026-10-11)\n")
    L.append("Law g_obs = nu(g_bar/a0) g_bar with nu(y) = 1/(1 - exp(-sqrt(y))), a0(0) = 1.2e-10 m/s^2, fitted for r = a0(z)/a0(0) "
             "on every galaxy with a total and a baryonic acceleration at the same radius. H_A: r = H(z)/H0 (Omega_m = 0.3). "
             "H_B: r = 1. MUSE: r = (1.0 + 1.59 z)/1.2 (Ciocan et al. 2026). chi2 uses the quoted errors, the sensitivity of the law "
             f"to g_bar, and an intrinsic scatter of {SIG_INT} dex. 68 % intervals from Delta chi2 = 1. Systematics: shift of log r when "
             "all baryonic masses move by +/- 0.2 dex (sys_M) or all velocities by +/- 10 % (sys_V).\n")

    zb_rc = [(0.6, 1.0), (1.0, 1.7), (1.7, 2.6)]
    zb_msa = [(0.5, 1.1), (1.1, 1.7)]
    zb_km = [(0.6, 1.1), (1.1, 1.9), (1.9, 2.7)]
    zb_kr = [(0.6, 1.1)]
    zb_z45 = [(4.0, 5.0)]
    tabs = [summarize(rc, "RC100", zb_rc), summarize(msa, "MSA-3D", zb_msa), summarize(km, "KMOS3D", zb_km),
            summarize(kr, "KROSS", zb_kr), summarize(z45, "z4.5 all", zb_z45), summarize(z45[z45.group == "curve"], "z4.5 curves", zb_z45)]
    T = pd.concat(tabs, ignore_index=True)

    L.append("## Fit of r = a0(z)/a0(0) per sample and redshift bin\n")
    L.append("| sample | z range | N | z med | g_bar/a0 med | nu_obs med | r (68 %, quoted errors) | r (errors rescaled to chi2/N = 1) | r_A = H/H0 | chi2 H_B | chi2 H_A | chi2 MUSE | (chi2_A - chi2_B) rescaled | sys_M (+0.2/-0.2 dex) | sys_V (+10/-10 %) | mass shift making H_A best |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in T.iterrows():
        L.append(f"| {r['sample']} | {r.zlo:.1f}-{r.zhi:.1f} | {r.N} | {r.z_med:.2f} | {r.x_med:.2f} | {r.nu_med:.2f} | "
                 f"{r.r:.2f} ({r.r_lo:.2f}-{r.r_hi:.2f}) | {10**r.logr_lo_r:.2f}-{10**r.logr_hi_r:.2f} | {r.rA:.2f} | {r.chi2_B:.1f} | {r.chi2_A:.1f} | {r.chi2_M:.1f} | "
                 f"{r.chi2A_minus_B_rescaled:+.0f} | {r.sys_M_plus:+.2f}/{r.sys_M_minus:+.2f} | {r.sys_V_plus:+.2f}/{r.sys_V_minus:+.2f} | "
                 f"{'none within +/-0.6 dex' if np.isnan(r.dM_for_A) else f'{r.dM_for_A:+.2f} dex'} |")
    L.append("")
    L.append("Reading: r is the a0 ratio that the sample prefers at face value; 'mass shift making H_A best' is the uniform change of all "
             "baryonic masses that would make a0 ~ H(z) the best fit. Where it is -0.2 dex or less in size, a 0.2 dex error in the "
             "baryonic masses decides the test; where it is larger, H_A needs an error larger than the mass budget allows.\n")
    # chi2 totals for z >= 1.7 only (where H_A predicts r >= 3)
    hz = pd.concat([rc[rc.z >= 1.7], z45[z45.group == "curve"]], ignore_index=True)
    chi_h, out_h = chi2_table(hz)
    lr_h, lo_h, hi_h, cmin_h = fit_r(chi_h)
    L.append(f"## z >= 1.7 only (RC100 z > 1.7 + z~4.5 curves, N = {len(hz)})\n")
    L.append(f"- r = {10**lr_h:.2f} ({10**lo_h:.2f}-{10**hi_h:.2f}); chi2(H_B) = {out_h['B'].sum():.1f}, chi2(H_A) = {out_h['A'].sum():.1f}, "
             f"Delta = {out_h['A'].sum() - out_h['B'].sum():+.1f} with the quoted errors, {(out_h['A'].sum() - out_h['B'].sum()) * len(hz) / cmin_h:+.0f} with errors rescaled to chi2/N = 1.")
    for dm in (-0.1, -0.2, -0.3):
        h2 = hz.copy(); h2["logg_bar"] += dm
        c2, o2 = chi2_table(h2)
        L.append(f"- baryonic masses lowered by {dm:+.1f} dex: r = {10**fit_r(c2)[0]:.2f}, Delta chi2 (A - B) = {o2['A'].sum() - o2['B'].sum():+.1f}")
    L.append("")
    # MUSE-DARK: what their a0(z) becomes once the molecular gas of Tacconi et al. 2018 is added (deep regime: a0 ~ 1/M_bar)
    L.append("## MUSE-DARK III corrected for the molecular gas it omits (a0 ~ 1/M_bar in the low-acceleration regime)\n")
    L.append("| log M* | z | a0 published (1e-10) | M_mol/M* (Tacconi 2018, main sequence) | a0 with H2 (1e-10) | with H2 and HI (xGASS-like) |")
    L.append("|---|---|---|---|---|---|")
    for lm in (9.3, 9.7, 10.1):
        for z, a0p in ((0.45, 1.99), (0.9, 2.38), (1.3, 2.71)):
            mu = float(mu_mol_tacconi(np.array([z]), np.array([lm]), np.array([10 ** ((0.84 - 0.026 * t_universe(z)) * lm - (6.51 - 0.11 * t_universe(z)))]))[0])
            muhi = 10 ** (-0.4 - 0.7 * (lm - 9.75))
            L.append(f"| {lm:.1f} | {z:.2f} | {a0p:.2f} | {mu:.2f} | {a0p/(1+mu):.2f} | {a0p/(1+mu+muhi):.2f} |")
    L.append("")

    # combined chi2 over the samples with homogeneous modelling
    core = pd.concat([rc, msa, km, z45[z45.group == "curve"]], ignore_index=True)
    chi, out = chi2_table(core)
    lr, lo, hi, cmin = fit_r(chi)
    L.append("## Combined (RC100 + MSA-3D + KMOS3D + z~4.5 curves, N = %d)\n" % len(core))
    L.append(f"- single r for all z: r = {10**lr:.2f} ({10**lo:.2f}-{10**hi:.2f}); chi2 = {cmin:.1f}")
    L.append(f"- chi2(H_B, r = 1) = {out['B'].sum():.1f} ; chi2(H_A, r = H(z)/H0) = {out['A'].sum():.1f} ; chi2(MUSE law) = {out['M'].sum():.1f}")
    L.append(f"- Delta chi2 (H_A - H_B) = {out['A'].sum() - out['B'].sum():+.1f} ; (MUSE - H_B) = {out['M'].sum() - out['B'].sum():+.1f}")
    # power-law fit r = (H/H0)^n
    ns = np.arange(-1.0, 2.01, 0.02)
    cn = []
    for n in ns:
        rr = Hz(core.z.values) ** n
        p = logg_pred(core.logg_bar.values, rr)
        sl = (logg_pred(core.logg_bar.values + 0.05, rr) - logg_pred(core.logg_bar.values - 0.05, rr)) / 0.1
        cn.append((((core.logg_obs - p) ** 2) / (core.s_obs ** 2 + (sl * core.s_bar) ** 2 + SIG_INT ** 2)).sum())
    cn = np.array(cn); i = int(np.argmin(cn)); ok = cn <= cn[i] + 1
    L.append(f"- a0 ~ H(z)^n : n = {ns[i]:.2f} ({ns[ok].min():.2f} to {ns[ok].max():.2f}) ; H_A is n = 1, H_B is n = 0")
    # same with masses shifted
    for dm in (0.2, -0.2):
        c2 = core.copy(); c2["logg_bar"] += dm
        cn2 = []
        for n in ns:
            rr = Hz(c2.z.values) ** n
            p = logg_pred(c2.logg_bar.values, rr)
            sl = (logg_pred(c2.logg_bar.values + 0.05, rr) - logg_pred(c2.logg_bar.values - 0.05, rr)) / 0.1
            cn2.append((((c2.logg_obs - p) ** 2) / (c2.s_obs ** 2 + (sl * c2.s_bar) ** 2 + SIG_INT ** 2)).sum())
        cn2 = np.array(cn2); j = int(np.argmin(cn2)); ok2 = cn2 <= cn2[j] + 1
        L.append(f"- with all baryonic masses shifted by {dm:+.1f} dex: n = {ns[j]:.2f} ({ns[ok2].min():.2f} to {ns[ok2].max():.2f})")
    L.append("")

    # ----- checks
    L.append("## Checks\n")
    dd = rc.logg_bar - rc.logg_bar_mm
    L.append(f"- RC100: authors' g_bar(Re) = (1 - f_DM) V_c^2/Re versus an exponential disc + Hernquist bulge built from their "
             f"log M_bar, Re, M_bulge: median difference {dd.median():+.3f} dex, 16-84 %: {np.percentile(dd,16):+.2f} to {np.percentile(dd,84):+.2f}.")
    off = kmraw.logMb - (3.75 * np.log10(kmraw.Vcirc) + 2.18)
    L.append(f"- KMOS3D: offset of log M_bar from the Lelli et al. 2016 bTFR at v_circ,max, median z < 1.1: "
             f"{off[kmraw.z < 1.1].median():+.2f} dex, z > 1.9: {off[kmraw.z > 1.9].median():+.2f} dex (Ubler et al. 2017: -0.44 and -0.27).")
    L.append(f"- KROSS: {len(kr)} galaxies kept (quality 1, rotation-dominated, no AGN, i > 25, v_C not extrapolated, error < 30 %); "
             f"median molecular gas fraction M_mol/M* = {kr.mu_mol.median():.2f}; the pressure-support term raises g_obs by a median "
             f"{(kr.logg_obs - kr.logg_obs_nopress).median():.2f} dex.")
    for lab, s in (("RC100", rc), ("MSA-3D", msa), ("KMOS3D", km), ("KROSS", kr), ("z4.5", z45)):
        L.append(f"- {lab}: N = {len(s)}, g_bar/a0 16-50-84 %: {np.percentile(10**s.logg_bar/A0,16):.2f} / {np.percentile(10**s.logg_bar/A0,50):.2f} / "
                 f"{np.percentile(10**s.logg_bar/A0,84):.2f}; nu_obs median {np.median(10**(s.logg_obs-s.logg_bar)):.2f}; "
                 f"law predicts nu median {np.median(nu(10**s.logg_bar/A0)):.2f} (H_B) and {np.median(nu(10**s.logg_bar/A0/Hz(s.z))):.2f} (H_A).")
    L.append("")
    # per-galaxy z~4.5 table
    L.append("## z ~ 4.5 discs, galaxy by galaxy\n")
    L.append("| galaxy | z | R (kpc) | g_bar/a0 | nu_obs | nu H_B | nu H_A (r = %.1f) | note |" % Hz(4.5))
    L.append("|---|---|---|---|---|---|---|---|")
    for _, g in z45.iterrows():
        x = 10 ** g.logg_bar / A0
        L.append(f"| {g['name']} | {g.z:.2f} | {g.R_kpc:.1f} | {x:.1f} | {10**(g.logg_obs-g.logg_bar):.2f} +/- {np.sqrt(g.s_obs**2+g.s_bar**2)*np.log(10)*10**(g.logg_obs-g.logg_bar):.2f} | "
                 f"{nu(x):.2f} | {nu(x/Hz(g.z)):.2f} | {g.group} |")
    L.append("")
    L.append("## Published constraints used as they stand\n")
    L.append("- MUSE-DARK III (Ciocan et al. 2026, A&A, arXiv:2604.22613): 79 SFGs, 0.33 < z < 1.44, log M* > 8.8; RAR fitted on model rotation "
             "curves gives a0 = 2.38 (+0.12/-0.10) e-10 at z ~ 1 and a0(z) = (1.0 + 1.59 z) e-10, from 1.99e-10 in the lowest z bin to 2.71e-10 "
             "in the highest. Baryons = stars + a constant-surface-density HI disc; no molecular gas. The authors state that a0 = 1.2e-10 at all z "
             "would need stellar masses larger by +0.2 dex (lowest bins) to +0.45 dex (highest), i.e. gas fractions of 37 to 65 % of the baryons, "
             "which is what Tacconi et al. 2018 give for log M* = 9-10 at z = 1-1.4 before any HI. Their evolution is degenerate with the gas they do not include.")
    L.append("- Puech et al. 2010 (A&A 510, A68): bTFR at z ~ 0.6 with V_flat from 3D modelling, zero-point shift -0.05 +/- 0.08 dex "
             "(no evolution within 0.08 dex). H_A predicts -0.14 dex at z = 0.6 (r = 1.39), H_B 0.")
    L.append("- Milgrom 2017 (arXiv:1703.06110) on the six Genzel et al. 2017 galaxies: g(R_1/2) = 3-11 a0, observed phantom fractions "
             "0.07-0.21 match MOND with the local a0; a value ~4 a0 at z ~ 2 is 'all but excluded' (would give fractions ~0.5).")
    L.append("- Lelli et al. 2023 (A&A 672, A106): two CO discs at z = 1.47 and 2.24, flat to 8 kpc, fitted with the local a0 and the local "
             "interpolation function with mass-to-light ratios 10-30 % below the Newtonian baryons-only fit; consistency, not a measurement of a0.")
    L.append("- Ubler et al. 2017 (ApJ 842, 121): bTFR zero point -0.44 dex at z ~ 0.9 and -0.27 dex at z ~ 2.3 relative to Lelli et al. 2016, "
             "with v_circ,max at 2.2 Rd: the sign and the non-monotonic z dependence are the signature of v_circ,max > V_flat in compact "
             "high-acceleration discs, not of a0; the direct fit above uses the same data with sizes.")
    L.append("")
    L.append("Figure: results/a0_highz.png. Per-galaxy accelerations: results/a0_highz.csv")
    open("results/a0_highz.md", "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))

    # ---------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    zz = np.linspace(0, 5, 200)
    for a in ax:
        a.plot(zz, Hz(zz), "k-", label="H_A : a0 ~ H(z)")
        a.plot(zz, np.ones_like(zz), "k--", label="H_B : a0 constant")
        a.plot(zz, r_muse(zz), "k:", label="MUSE-DARK III law")
        a.set_xlabel("z"); a.set_ylabel("a0(z) / 1.2e-10")
        a.set_yscale("log")
    cols = {"RC100": "C0", "MSA-3D": "C1", "KMOS3D": "C2", "KROSS": "C3", "z4.5 curves": "C4", "z4.5 all": "C5"}
    seen = set()
    for _, r in T.iterrows():
        a = ax[0]
        lab = None if r["sample"] in seen else r["sample"]
        seen.add(r["sample"])
        a.errorbar(r.z_med, r.r, yerr=[[r.r - r.r_lo], [r.r_hi - r.r]], fmt="o", color=cols.get(r["sample"], "k"), capsize=3, label=lab)
        lo = 10 ** (r.logr + min(r.sys_M_plus, r.sys_M_minus, 0)); hi = 10 ** (r.logr + max(r.sys_M_plus, r.sys_M_minus, 0))
        a.plot([r.z_med, r.z_med], [lo, hi], "-", color=cols.get(r["sample"], "k"), alpha=0.35, lw=6)
    ax[0].errorbar(0.6, 10 ** 0.05, yerr=[[10 ** 0.05 - 10 ** -0.03], [10 ** 0.13 - 10 ** 0.05]], fmt="s", color="grey", capsize=3, label="Puech+10 bTFR z~0.6")
    ax[0].errorbar([0.45, 1.3], [1.99 / 1.2, 2.71 / 1.2], yerr=0.15, fmt="^", color="purple", capsize=3, label="MUSE-DARK bins (no H2)")
    ax[0].set_ylim(0.2, 12); ax[0].legend(fontsize=7, loc="upper left"); ax[0].set_title("a0(z)/a0(0) fitted per sample (thick bar: masses +/- 0.2 dex)")
    # right: per-galaxy nu_obs vs x for z > 1.7 with the two predictions
    a = ax[1]
    s = allg[(allg.z > 1.7) & (allg["sample"].isin(["RC100", "KMOS3D"]))]
    a.errorbar(s.x, s.nu_obs, yerr=s.nu_obs * np.sqrt(s.s_obs ** 2 + s.s_bar ** 2) * np.log(10), fmt=".", color="C0", alpha=0.5, label="RC100 + KMOS3D, z > 1.7")
    s = allg[allg["sample"] == "z4.5"]
    a.errorbar(s.x, s.nu_obs, yerr=s.nu_obs * np.sqrt(s.s_obs ** 2 + s.s_bar ** 2) * np.log(10), fmt="D", color="C4", alpha=0.8, label="z ~ 4.5")
    xx = np.logspace(-1, 2.5, 200)
    a.plot(xx, nu(xx), "k--", label="H_B (a0 constant)")
    a.plot(xx, nu(xx / Hz(2.2)), "k-", label="H_A at z = 2.2 (r = %.1f)" % Hz(2.2))
    a.plot(xx, nu(xx / Hz(4.5)), "k-.", label="H_A at z = 4.5 (r = %.1f)" % Hz(4.5))
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("g_bar / a0(0)"); a.set_ylabel("nu_obs = g_obs / g_bar")
    a.set_ylim(0.5, 10); a.set_xlim(0.1, 300); a.legend(fontsize=7); a.set_title("Mass discrepancy at the measured radius")
    a.set_yscale("log"); a.set_ylim(0.6, 8)
    plt.tight_layout(); plt.savefig("results/a0_highz.png", dpi=130)
    print("figure written")


if __name__ == "__main__":
    main()
