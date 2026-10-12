# -*- coding: utf-8 -*-
"""
a0 at high redshift, second pass: samples with measured cold gas (CO, [CI], [CII]) and the rotation-curve shape test.

Adds to scripts/a0_highz.py:
  NOEMA3D  Jolly et al. 2026 (arXiv:2604.18503/18504): 10 massive main-sequence discs at z = 1.1-1.6, CO gas masses, DysmalPy
           fits (log M_bary, f_DM(Re), Re_disk, Re_bulge, B/T)
  CRISTAL  Lee et al. 2025 (arXiv:2507.11600): 14 [CII] discs at z = 4.4-5.7, DysmalPy fits (log M_tot, f_DM(Re), Re, B/T)
  ALPAKA   Roman-Oliveira, Rizzo & Fraternali 2026 (arXiv:2601.03338): 3 discs (z = 0.56, 1.45, 2.10) with CO/[CI] curves to
           3.8-4.2 Re, JWST stellar masses and sizes, gas normalisation fitted
  z~4.5    Rizzo et al. 2021 (five lensed DSFGs, baryonic masses and asymptotic velocities), Fraternali et al. 2021 (two),
           Lelli et al. 2021 (ALESS 073.1), plus the four Roman-Oliveira et al. 2024 curves (shape test)
For the f_DM-based samples the observable is nu_obs = 1/(1 - f_DM) and g_bar comes from the authors' baryonic mass model
(exponential disc + Hernquist bulge with their masses and sizes); velocities only enter through the authors' fits.
Outputs: results/a0_highz_coldgas.md, results/a0_highz_coldgas.csv, results/a0_highz_coldgas.png
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a0_highz import (A0, D, LOGR_GRID, SIG_INT, Hz, g_si, load_kmos3d, load_kross, load_msa3d, load_rc100, load_z45,  # noqa: E402
                      logg_pred, nu, r_muse, v2_expdisc, v2_hernquist)

LN10 = np.log(10)


# ---------------------------------------------------------------- likelihood on nu_obs (f_DM samples)
def lognu_pred(logg_bar, r):
    return np.log10(nu(10 ** logg_bar / (r * A0)))


def chi2_nu(df):
    """df: z, lognu_obs, s_nu, logg_bar, s_bar. Returns chi2[galaxy, grid] and chi2 for H_B, H_A, MUSE."""
    lb, ln = df.logg_bar.values[:, None], df.lognu_obs.values[:, None]
    sn, sb = df.s_nu.values[:, None], df.s_bar.values[:, None]

    def chi(rr):
        p = lognu_pred(lb, rr)
        sl = (lognu_pred(lb + 0.05, rr) - lognu_pred(lb - 0.05, rr)) / 0.1
        return (ln - p) ** 2 / (sn ** 2 + (sl * sb) ** 2 + SIG_INT ** 2)

    grid = chi(10 ** LOGR_GRID[None, :])
    out = {"B": chi(np.ones((len(df), 1)))[:, 0], "A": chi(Hz(df.z.values)[:, None])[:, 0], "M": chi(r_muse(df.z.values)[:, None])[:, 0]}
    return grid, out


def chi2_rar(df):
    lb, lo = df.logg_bar.values[:, None], df.logg_obs.values[:, None]
    so, sb = df.s_obs.values[:, None], df.s_bar.values[:, None]

    def chi(rr):
        p = logg_pred(lb, rr)
        sl = (logg_pred(lb + 0.05, rr) - logg_pred(lb - 0.05, rr)) / 0.1
        return (lo - p) ** 2 / (so ** 2 + (sl * sb) ** 2 + SIG_INT ** 2)

    grid = chi(10 ** LOGR_GRID[None, :])
    out = {"B": chi(np.ones((len(df), 1)))[:, 0], "A": chi(Hz(df.z.values)[:, None])[:, 0], "M": chi(r_muse(df.z.values)[:, None])[:, 0]}
    return grid, out


def chi2_any(df):
    return chi2_nu(df) if df.kind.iloc[0] == "nu" else chi2_rar(df)


def fit_r(chi):
    c = chi.sum(0)
    i = int(np.argmin(c))
    ok = c <= c[i] + 1.0
    return LOGR_GRID[i], LOGR_GRID[ok].min(), LOGR_GRID[ok].max(), c[i]


def shift_mass(df, dm):
    """All baryonic masses wrong by dm dex: g_bar moves by dm; the measured V_c does not, so for the f_DM samples
    nu_obs = V_c^2 / V_bar^2 moves by -dm (the authors would have fitted a different f_DM)."""
    s2 = df.copy()
    s2["logg_bar"] = s2.logg_bar + dm
    if s2.kind.iloc[0] == "nu":
        s2["lognu_obs"] = s2.lognu_obs - dm
    return s2


def summarize(df, label, zbins):
    rows = []
    for zlo, zhi in zbins:
        s = df[(df.z >= zlo) & (df.z < zhi)]
        if len(s) < 2:
            continue
        chi, out = chi2_any(s)
        lr, lo, hi, cmin = fit_r(chi)
        scale = max(cmin / len(s), 1e-3)
        c = chi.sum(0)
        okr = c <= cmin + scale
        dms = np.arange(-0.6, 0.61, 0.02)
        rf = np.array([fit_r(chi2_any(shift_mass(s, dm))[0])[0] for dm in dms])
        zmed = s.z.median()
        k = int(np.argmin(np.abs(rf - np.log10(Hz(zmed)))))
        dm_A = dms[k] if np.abs(rf[k] - np.log10(Hz(zmed))) < 0.1 else np.nan
        x = 10 ** s.logg_bar / A0
        nuo = 10 ** s.lognu_obs if "lognu_obs" in s else 10 ** (s.logg_obs - s.logg_bar)
        rows.append(dict(sample=label, zlo=zlo, zhi=zhi, N=len(s), z_med=zmed, x_med=np.median(x), nu_med=np.median(nuo),
                         r=10 ** lr, r_lo=10 ** lo, r_hi=10 ** hi, r_lo_resc=10 ** LOGR_GRID[okr].min(), r_hi_resc=10 ** LOGR_GRID[okr].max(),
                         rA=Hz(zmed), chi2_B=out["B"].sum(), chi2_A=out["A"].sum(), chi2_M=out["M"].sum(), chi2_min=cmin,
                         dB=out["B"].sum() - cmin, dA=out["A"].sum() - cmin,
                         dAB_resc=(out["A"].sum() - out["B"].sum()) / scale,
                         sys_M_plus=fit_r(chi2_any(shift_mass(s, 0.2))[0])[0] - lr, sys_M_minus=fit_r(chi2_any(shift_mass(s, -0.2))[0])[0] - lr,
                         dM_for_A=dm_A))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- loaders
def inc_err_dex(inc_deg, dinc_deg=5.0):
    """error on log V^2 from the inclination: 2 cot(i) di / ln10"""
    return 2 * np.abs(1 / np.tan(np.radians(inc_deg))) * np.radians(dinc_deg) / LN10


def to_nu_frame(sample, name, z, R, gbar, s_bar, f, ef, extra=None):
    f = np.clip(np.asarray(f, float), 0.0, 0.95)
    ef = np.clip(np.asarray(ef, float), 0.02, None)
    d = pd.DataFrame(dict(sample=sample, name=name, z=z, R_kpc=R, kind="nu",
                          logg_bar=np.log10(gbar), s_bar=s_bar, fDM=f, e_fDM=ef,
                          lognu_obs=-np.log10(1 - f), s_nu=ef / (1 - f) / LN10))
    d["logg_obs"] = d.logg_bar + d.lognu_obs
    if extra is not None:
        for k, v in extra.items():
            d[k] = np.asarray(v)
    return d


def load_noema3d():
    t = pd.read_csv(f"{D}/noema3d_jolly2026_dysmalpy.csv")
    M = 10 ** t.logMbary
    Mb = t.BT * M
    v2 = v2_expdisc(t.Re_disk, M - Mb, t.Re_disk / 1.678) + v2_hernquist(t.Re_disk, Mb, t.Re_bulge)
    gbar = g_si(v2, t.Re_disk)
    s_bar = np.sqrt(((t.e_logMbary_up + t.e_logMbary_lo) / 2) ** 2 + 0.1 ** 2)
    ef = np.sqrt(((t.e_fDM_up + t.e_fDM_lo) / 2) ** 2 + (0.5 * inc_err_dex(t.inc) * LN10 * (1 - t.fDM)) ** 2)  # inclination moves V^2, hence f
    d = to_nu_frame("NOEMA3D", t.name, t.z, t.Re_disk, gbar, s_bar, t.fDM, ef,
                    extra=dict(Vc_Re=t.Vc_Re, sigma0=t.sigma0, inc=t.inc, logMstar=t.logMstar, logMmol=t.logMmol, logMbary=t.logMbary))
    # check: g_obs from V_c(Re) directly versus nu_obs g_bar
    d["logg_obs_V"] = np.log10(g_si(t.Vc_Re ** 2, t.Re_disk))
    return d


def load_cristal():
    t = pd.read_csv(f"{D}/cristal_lee2025_dysmalpy.csv")
    M = 10 ** t.logMtot
    Mb = t.BT * M
    v2 = v2_expdisc(t.Re_disk, M - Mb, t.Re_disk / 1.678) + v2_hernquist(t.Re_disk, Mb, 1.0)
    gbar = g_si(v2, t.Re_disk)
    s_bar = np.sqrt(((t.e_logMtot_up + t.e_logMtot_lo) / 2) ** 2 + ((t.e_Re_up + t.e_Re_lo) / 2 / t.Re_disk / LN10) ** 2)
    ef = (t.e_fDM_up + t.e_fDM_lo) / 2
    d = to_nu_frame("CRISTAL", "CRISTAL-" + t.ID.astype(str), t.z, t.Re_disk, gbar, s_bar, t.fDM, ef,
                    extra=dict(Vrot_Re=t.Vrot_Re, sigma0=t.sigma0, inc=t.inc, logMtot=t.logMtot, Rout_over_Re=t.Rout_over_Re))
    d["logg_obs_V"] = np.log10(g_si(t.Vrot_Re ** 2 + 3.36 * t.sigma0 ** 2, t.Re_disk))
    return d


def load_alpaka3():
    t = pd.read_csv(f"{D}/alpaka_romanoliveira2026.csv")
    R = t.nRE * t.Reff_gas
    Ms, Mg = 10 ** t.logMstar_fit, 10 ** t.logMgas
    Mb = t.BT * Ms
    Rd_gas = t.Reff_gas / 1.678
    v2 = v2_expdisc(R, Ms - Mb, t.Reff_stellar_disc / 1.678) + v2_hernquist(R, Mb, t.Reff_bulge) + v2_expdisc(R, Mg, Rd_gas)
    gbar = g_si(v2, R)
    vc2 = t.V_ext ** 2 + 2 * t.sigma_ext ** 2 * (R / Rd_gas)
    gobs = g_si(vc2, R)
    eV = (t.e_V_ext_up + t.e_V_ext_lo) / 2 / t.V_ext
    s_obs = np.sqrt((2 * eV / LN10) ** 2 + inc_err_dex(t.inc_ALMA, (t.e_inc_up + t.e_inc_lo) / 2) ** 2 + (0.1 / LN10) ** 2)
    sM = np.sqrt((Ms * (t.e_logMstar_up + t.e_logMstar_lo) / 2) ** 2 + (Mg * (t.e_logMgas_up + t.e_logMgas_lo) / 2) ** 2) / (Ms + Mg)
    d = pd.DataFrame(dict(sample="ALPAKA+JWST", name=t.name, z=t.z, R_kpc=R, kind="rar", logg_obs=np.log10(gobs), s_obs=s_obs,
                          logg_bar=np.log10(gbar), s_bar=np.sqrt(sM ** 2 + 0.05 ** 2), V_ext=t.V_ext, inc=t.inc_ALMA))
    return d


def load_z45_assumed():
    """Rizzo et al. 2021 (their M_bar and asymptotic V_t; outer radius not tabulated: 2.5 kpc assumed),
    Fraternali et al. 2021 (3 scale lengths assumed), Lelli et al. 2021 ALESS 073.1 (2.5 kpc)."""
    RZ = [  # name, z, M_bar(1e10), e, R_bar half-mass (kpc), V, eV, R_out
        ("SPT0113-46", 4.23, 11.0, 1.0, 2.2, 352, 8, 2.5), ("SPT0345-47", 4.29, 4.0, 0.3, 0.64, 200, 25, 2.5),
        ("SPT0441-46", 4.48, 3.2, 0.2, 0.35, 293, 7, 2.5), ("SPT2146-55", 4.57, 2.2, 0.2, 1.3, 176, 18, 2.5),
        ("SPT2132-58", 4.77, 3.9, 0.35, 2.5, 219, 22, 2.5),
        ("AzTEC/C159", 4.567, 19.5, 6.0, 1.3 * 1.678, 506, 110, 3.9), ("J1000+0234", 4.539, 16.5, 8.0, 1.5 * 1.678, 538, 40, 4.5),
        ("ALESS073.1", 4.755, 5.2, 1.5, 1.2, 309, 56, 2.5),
    ]
    rows = []
    for n, z, M, eM, Rh, V, eV, R in RZ:
        v2 = v2_expdisc(R, M * 1e10, Rh / 1.678)
        rows.append(dict(sample="z4.5 assumed", name=n, z=z, R_kpc=R, kind="rar", logg_obs=np.log10(g_si(V ** 2, R)),
                         s_obs=np.sqrt((2 * eV / V / LN10) ** 2 + 0.18 ** 2), logg_bar=np.log10(g_si(v2, R)),
                         s_bar=np.sqrt((eM / M / LN10) ** 2 + 0.3 ** 2)))   # 0.3 dex: alpha_[CII] and R_out
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- shape test on the Roman-Oliveira 2024 curves
def shape_test():
    cur = pd.read_csv(f"{D}/z45/circular_speed_romanoliveira2024.csv")
    RO = {"BRI1335-0417": (4.41, 10.4, 9.8, 2.3, 2.9), "J081740": (4.54, 10.6, 10.5, 2.4, 2.8),
          "SGP38326-1": (4.28, 11.0, 11.3, 2.2, 2.8), "SGP38326-2": (4.28, 10.3, 11.1, 1.7, 2.0)}
    dmg = np.arange(-0.8, 0.81, 0.02)
    out, prof_all = [], np.zeros(len(LOGR_GRID))
    for n, (z, lms, lmg, res, reb) in RO.items():
        c = cur[cur.galaxy == n]
        R, V, eV = c.R_kpc.values, c.Vcirc_kms.values, np.maximum(c.e_Vcirc_kms.values, 0.1 * c.Vcirc_kms.values)
        prof = np.zeros(len(LOGR_GRID))
        best = None
        for j, lr in enumerate(LOGR_GRID):
            r = 10 ** lr
            chis = []
            for dm in dmg:
                v2 = (v2_hernquist(R, 10 ** (lms + dm), res) + v2_expdisc(R, 10 ** (lmg + dm), reb / 1.678))
                gb = g_si(v2, R)
                vmod = np.sqrt(nu(gb / (r * A0)) * v2)
                chis.append((((V - vmod) / eV) ** 2).sum() + (dm / 0.3) ** 2)
            chis = np.array(chis)
            prof[j] = chis.min()
            if best is None or prof[j] < best[1]:
                best = (lr, prof[j], dmg[int(np.argmin(chis))])
        ok = prof <= prof.min() + 1
        rA = Hz(z)
        jA = int(np.argmin(np.abs(LOGR_GRID - np.log10(rA)))); jB = int(np.argmin(np.abs(LOGR_GRID)))
        out.append(dict(name=n, z=z, N=len(R), R_max=R.max(), r=10 ** best[0], r_lo=10 ** LOGR_GRID[ok].min(), r_hi=10 ** LOGR_GRID[ok].max(),
                        dlogM_best=best[2], chi2_min=best[1], chi2_B=prof[jB], chi2_A=prof[jA], rA=rA))
        prof_all += prof
    ok = prof_all <= prof_all.min() + 1
    comb = dict(r=10 ** LOGR_GRID[int(np.argmin(prof_all))], r_lo=10 ** LOGR_GRID[ok].min(), r_hi=10 ** LOGR_GRID[ok].max(),
                chi2_B=prof_all[int(np.argmin(np.abs(LOGR_GRID)))], chi2_A=prof_all[int(np.argmin(np.abs(LOGR_GRID - np.log10(Hz(4.4)))))])
    return pd.DataFrame(out), comb


# ---------------------------------------------------------------- main
def main():
    rc = load_rc100(); rc["kind"] = "nu"; rc["lognu_obs"] = -np.log10(1 - rc.fDM); rc["s_nu"] = rc.e_fDM / (1 - rc.fDM) / LN10
    rc["s_bar"] = np.sqrt(rc.s_obs ** 2 + 0.05 ** 2)   # authors' g_bar = (1 - f) V^2/Re: V and Re errors move g_bar
    msa = load_msa3d(); msa["kind"] = "nu"; msa["lognu_obs"] = -np.log10(1 - msa.fDM); msa["s_nu"] = msa.e_fDM / (1 - msa.fDM) / LN10
    msa["s_bar"] = np.sqrt(msa.s_obs ** 2 + 0.05 ** 2)
    km, _ = load_kmos3d(); km["kind"] = "rar"
    kr = load_kross(); kr["kind"] = "rar"
    z45 = load_z45(); z45["kind"] = "rar"
    z45c = z45[z45.group == "curve"].copy(); z45c["sample"] = "z4.5 curves"
    no = load_noema3d(); cr = load_cristal(); al = load_alpaka3(); za = load_z45_assumed()
    allg = pd.concat([rc, msa, km, kr, z45c, no, cr, al, za], ignore_index=True)
    allg["x"] = 10 ** allg.logg_bar / A0
    allg["rA"] = Hz(allg.z)
    allg.to_csv("results/a0_highz_coldgas.csv", index=False)

    T = pd.concat([summarize(rc, "RC100", [(0.6, 1.0), (1.0, 1.7), (1.7, 2.6)]),
                   summarize(msa, "MSA-3D", [(0.5, 1.7)]),
                   summarize(km, "KMOS3D", [(0.6, 1.1), (1.1, 1.9), (1.9, 2.7)]),
                   summarize(kr, "KROSS", [(0.6, 1.1)]),
                   summarize(no, "NOEMA3D (CO gas)", [(1.0, 1.7)]),
                   summarize(al, "ALPAKA+JWST (CO gas, 4 Re)", [(0.5, 1.0), (1.0, 1.7), (1.7, 2.6)]),
                   summarize(cr, "CRISTAL ([CII], z~5)", [(4.0, 6.0)]),
                   summarize(z45c, "z4.5 curves", [(4.0, 5.0)]),
                   summarize(za, "z4.5 assumed R_out", [(4.0, 5.0)])], ignore_index=True)
    sh, comb = shape_test()

    L = ["# a0 at high redshift, second pass: measured cold gas and curve shapes (2026-10-12)\n",
         "Same test as results/a0_highz.md (r = a0(z)/a0(0); H_A: r = H(z)/H0; H_B: r = 1). For the f_DM samples the observable is now "
         "nu_obs = 1/(1 - f_DM) against nu(g_bar/(r a0)) with g_bar from the authors' own baryonic masses and sizes (exponential disc + "
         "Hernquist bulge); velocities enter only through the authors' fits. New samples: NOEMA3D (CO gas measured), CRISTAL ([CII] discs at "
         "z = 4.4-5.7), ALPAKA+JWST (three discs with curves to 4 Re and fitted gas normalisation), the Rizzo et al. 2021 lensed discs with "
         "their own baryonic masses. Systematics as before (sys_M: all masses +/- 0.2 dex; 'mass shift for H_A': uniform change making a0 ~ H(z) "
         "the best fit).\n",
         "## r = a0(z)/a0(0) per sample\n",
         "| sample | z | N | g_bar/a0 med | nu_obs med | r (68 %) | r, errors rescaled | r_A | chi2 min | H_B: chi2 - min | H_A: chi2 - min | (A-B) rescaled | sys_M +0.2/-0.2 | mass shift for H_A |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in T.iterrows():
        L.append(f"| {r['sample']} | {r.zlo:.1f}-{r.zhi:.1f} | {r.N} | {r.x_med:.2f} | {r.nu_med:.2f} | {r.r:.2f} ({r.r_lo:.2f}-{r.r_hi:.2f}) | "
                 f"{r.r_lo_resc:.2f}-{r.r_hi_resc:.2f} | {r.rA:.2f} | {r.chi2_min:.1f} | {r.dB:+.1f} | {r.dA:+.1f} | {r.dAB_resc:+.0f} | "
                 f"{r.sys_M_plus:+.2f}/{r.sys_M_minus:+.2f} | {'none' if np.isnan(r.dM_for_A) else f'{r.dM_for_A:+.2f} dex'} |")
    L.append("")
    # combined high-z with measured gas
    hz = pd.concat([rc[rc.z >= 1.7], cr, al[al.z > 1.7], z45c], ignore_index=True)
    parts = {"RC100 z>1.7": rc[rc.z >= 1.7], "CRISTAL": cr, "ALPAKA ID13": al[al.z > 1.7], "z4.5 curves": z45c}
    L.append("## z >= 1.7, all samples (N = %d)\n" % len(hz))
    tot = {"B": 0.0, "A": 0.0, "M": 0.0}
    prof = np.zeros(len(LOGR_GRID))
    for k, s in parts.items():
        g, o = chi2_any(s)
        prof += g.sum(0)
        for key in tot:
            tot[key] += o[key].sum()
        L.append(f"- {k}: N = {len(s)}, chi2(H_B) = {o['B'].sum():.1f}, chi2(H_A) = {o['A'].sum():.1f}")
    i = int(np.argmin(prof)); ok = prof <= prof[i] + 1
    L.append(f"- combined: r = {10**LOGR_GRID[i]:.2f} ({10**LOGR_GRID[ok].min():.2f}-{10**LOGR_GRID[ok].max():.2f}); chi2(H_B) = {tot['B']:.1f}, "
             f"chi2(H_A) = {tot['A']:.1f}, chi2(MUSE) = {tot['M']:.1f}; Delta(A - B) = {tot['A'] - tot['B']:+.1f} with quoted errors")
    for dm in (-0.1, -0.2, -0.3):
        p2 = np.zeros(len(LOGR_GRID)); tb = ta = 0.0
        for k, s in parts.items():
            g, o = chi2_any(shift_mass(s, dm)); p2 += g.sum(0); tb += o["B"].sum(); ta += o["A"].sum()
        L.append(f"- all baryonic masses {dm:+.1f} dex: r = {10**LOGR_GRID[int(np.argmin(p2))]:.2f}, Delta(A - B) = {ta - tb:+.1f}")
    L.append("")
    L.append("## Shape test on the four z ~ 4.4 curves (Roman-Oliveira et al. 2024), baryonic mass free with a 0.3 dex prior\n")
    L.append("| galaxy | z | N rings | R_max (kpc) | r (68 %) | best dlog M | chi2 min | chi2 H_B | chi2 H_A | r_A |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for _, s in sh.iterrows():
        L.append(f"| {s['name']} | {s.z:.2f} | {s.N} | {s.R_max:.1f} | {s.r:.2f} ({s.r_lo:.2f}-{s.r_hi:.2f}) | {s.dlogM_best:+.2f} | {s.chi2_min:.1f} | {s.chi2_B:.1f} | {s.chi2_A:.1f} | {s.rA:.1f} |")
    L.append(f"\nCombined: r = {comb['r']:.2f} ({comb['r_lo']:.2f}-{comb['r_hi']:.2f}); chi2(H_B) = {comb['chi2_B']:.1f}, chi2(H_A at z = 4.4) = {comb['chi2_A']:.1f}. "
             "The mass is marginalised, so this is the curve shape (rise of nu with radius) speaking, not the amplitude.\n")
    L.append("## Checks\n")
    for lab, s in (("NOEMA3D", no), ("CRISTAL", cr)):
        dd = s.logg_obs_V - s.logg_obs
        L.append(f"- {lab}: g_obs from the authors' velocity at Re (V_c, or V_rot^2 + 3.36 sigma0^2) minus nu_obs g_bar(mass model): "
                 f"median {dd.median():+.2f} dex, 16-84 %: {np.percentile(dd,16):+.2f} to {np.percentile(dd,84):+.2f}.")
    L.append(f"- NOEMA3D inclinations: {', '.join(f'{n} {i:.0f} deg' for n, i in zip(no.name, no.inc))}; the three discs below 21 deg carry "
             "0.3-0.4 dex inclination errors on nu_obs.")
    L.append(f"- CRISTAL: sigma0 = {cr.sigma0.min():.0f}-{cr.sigma0.max():.0f} km/s against V_rot(Re) = {cr.Vrot_Re.min():.0f}-{cr.Vrot_Re.max():.0f}: "
             "the pressure-support term is comparable to the rotation term, and the f_DM values have 68 % intervals of 0.3-0.5.")
    L.append("- ALPAKA+JWST galaxy by galaxy (outer ring, 3.8-4.2 Re):")
    for _, g in al.iterrows():
        x = 10 ** g.logg_bar / A0; nuo = 10 ** (g.logg_obs - g.logg_bar)
        L.append(f"  - {g['name']} z = {g.z:.2f}, R = {g.R_kpc:.1f} kpc, i = {g.inc:.0f} deg, g_bar/a0 = {x:.2f}, nu_obs = {nuo:.2f} "
                 f"(+/- {np.sqrt(g.s_obs**2+g.s_bar**2)*LN10*nuo:.2f}), law H_B nu = {nu(x):.2f}, H_A nu = {nu(x/Hz(g.z)):.2f}")
    L.append("")
    L.append("Figure: results/a0_highz_coldgas.png. Per-galaxy values: results/a0_highz_coldgas.csv")
    open("results/a0_highz_coldgas.md", "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    zz = np.linspace(0, 6, 200)
    ax.plot(zz, Hz(zz), "k-", label="H_A : a0 ~ H(z)")
    ax.plot(zz, np.ones_like(zz), "k--", label="H_B : a0 constant")
    ax.plot(zz, r_muse(zz), "k:", label="MUSE-DARK III law")
    cols = {"RC100": "C0", "MSA-3D": "C1", "KMOS3D": "C2", "KROSS": "C3", "NOEMA3D (CO gas)": "C6", "ALPAKA+JWST (CO gas, 4 Re)": "C8",
            "CRISTAL ([CII], z~5)": "C9", "z4.5 curves": "C4", "z4.5 assumed R_out": "C5"}
    seen = set()
    for _, r in T.iterrows():
        lab = None if r["sample"] in seen else r["sample"]; seen.add(r["sample"])
        ax.errorbar(r.z_med, r.r, yerr=[[max(r.r - r.r_lo, 1e-3)], [max(r.r_hi - r.r, 1e-3)]], fmt="o", color=cols.get(r["sample"], "k"), capsize=3, label=lab)
        lo = 10 ** (np.log10(r.r) + min(r.sys_M_plus, r.sys_M_minus, 0)); hi = 10 ** (np.log10(r.r) + max(r.sys_M_plus, r.sys_M_minus, 0))
        ax.plot([r.z_med, r.z_med], [lo, hi], "-", color=cols.get(r["sample"], "k"), alpha=0.3, lw=6)
    ax.errorbar(4.4, comb["r"], yerr=[[max(comb["r"] - comb["r_lo"], 1e-3)], [max(comb["r_hi"] - comb["r"], 1e-3)]], fmt="*", color="purple", ms=12, capsize=3, label="shape test, 4 curves z~4.4")
    ax.set_yscale("log"); ax.set_ylim(0.1, 15); ax.set_xlabel("z"); ax.set_ylabel("a0(z) / 1.2e-10")
    ax.legend(fontsize=7, loc="upper left", ncol=2); ax.set_title("a0(z)/a0(0): fitted per sample (thick bar: all baryonic masses +/- 0.2 dex)")
    plt.tight_layout(); plt.savefig("results/a0_highz_coldgas.png", dpi=130)
    print("figure written")


if __name__ == "__main__":
    main()
