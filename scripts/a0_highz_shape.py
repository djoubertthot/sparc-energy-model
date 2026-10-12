# -*- coding: utf-8 -*-
"""
a0 at high redshift, third pass: the rotation-curve SHAPE test, which does not depend on the baryonic mass scale.

Samples:
  KURVS  Puglisi et al. 2023 (MNRAS 524, 2814; arXiv:2305.04382): 22 star-forming galaxies at z ~ 1.5 with H-alpha rotation
         curves to 4 Re (10-15 kpc). Table 4 gives the intrinsic rotation velocity at 3 and 6 disc scale radii (t = v6/v3, the
         outer slope of Tiley et al. 2019), Table 2 sigma0 and v/sigma0, Table 3 f_DM(<Re) for the ten rotation-dominated discs.
  z~4.4  the four Roman-Oliveira et al. 2024 [CII] curves (shape test of scripts/a0_highz_coldgas.py, repeated here).
Two uses of KURVS:
  (a) amplitude: f_DM(<Re) against the law, g_bar from M* + molecular gas (Tacconi et al. 2018), same test as RC100;
  (b) shape: the outer slope t, converted to circular velocities with the Burkert et al. 2010 pressure correction
      (V_c^2 = v^2 + 2 sigma0^2 R/Rd, i.e. + 6 sigma0^2 at 3 Rd and + 12 sigma0^2 at 6 Rd), against the slope the law predicts
      for an exponential disc of the galaxy's size, with the baryonic mass free (0.3 dex prior). H_A (a0 ~ H(z), r = 2.3 at
      z = 1.5) predicts outer slopes closer to flat or rising than H_B.
Outputs: results/a0_highz_shape.md, results/a0_highz_shape.png
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a0_highz import A0, D, LOGR_GRID, Hz, g_si, mu_mol_tacconi, nu, v2_expdisc  # noqa: E402
from a0_highz_coldgas import chi2_nu, fit_r, shape_test, shift_mass, to_nu_frame  # noqa: E402

LN10 = np.log(10)


def load_kurvs():
    t = pd.read_csv(f"{D}/kurvs_puglisi2023.csv")
    mu = mu_mol_tacconi(t.z.values, t.logMstar.values, t.SFR.values)
    t["logMbar"] = np.log10(10 ** t.logMstar * (1 + mu))
    t["mu_mol"] = mu
    t["Rd"] = t.Reff / 1.678
    return t


def kurvs_amplitude(t):
    s = t[t.fDM_Re.notna()].copy()
    M = 10 ** s.logMbar
    gbar = g_si(v2_expdisc(s.Reff, M, s.Rd), s.Reff)
    s_bar = np.sqrt(0.2 ** 2 + 0.1 ** 2 + (s.e_Reff / s.Reff / LN10) ** 2)   # M* 0.15, gas 0.15, size
    ef = np.sqrt(s.e_fDM_Re ** 2 + 0.05 ** 2)
    d = to_nu_frame("KURVS f_DM", "KURVS-" + s.ID.astype(str), s.z.values, s.Reff.values, gbar.values, s_bar.values, s.fDM_Re.values, ef.values,
                    extra=dict(vsig=s.vsig.values, sigma0=s.sigma0.values))
    return d


GAS_EXT = 1.0   # scale length of the gas disc in units of the stellar one (set per variant)
GAS_FRAC = None  # gas fraction of the baryons per galaxy (set when GAS_EXT != 1)


def model_slope(M, Rd, r, fgas=0.0):
    """V_c(6 Rd) / V_c(3 Rd) for the baryonic disc: stars (mass M (1 - fgas), scale length Rd) and gas (M fgas, GAS_EXT Rd),
    under the law with a0 -> r a0. Rd in kpc, M in Msun."""
    R3, R6 = 3 * Rd, 6 * Rd
    v23 = v2_expdisc(R3, M * (1 - fgas), Rd) + v2_expdisc(R3, M * fgas, GAS_EXT * Rd)
    v26 = v2_expdisc(R6, M * (1 - fgas), Rd) + v2_expdisc(R6, M * fgas, GAS_EXT * Rd)
    g3, g6 = g_si(v23, R3), g_si(v26, R6)
    return np.sqrt(nu(g6 / (r * A0)) * v26 / (nu(g3 / (r * A0)) * v23)), g3 / A0


def kurvs_shape(t, vsig_min, pressure=True, gas_ext=1.0):
    global GAS_EXT
    GAS_EXT = gas_ext
    s = t[(t.vsig >= vsig_min) & (t.v_3Rd > 20)].copy()
    rows = []
    dmg = np.arange(-0.8, 0.81, 0.02)
    prof_all = np.zeros(len(LOGR_GRID))
    chiA_all = chiB_all = 0.0
    for _, g in s.iterrows():
        if pressure:
            c3, c6 = g.v_3Rd ** 2 + 6 * g.sigma0 ** 2, g.v_6Rd ** 2 + 12 * g.sigma0 ** 2
            tobs = np.sqrt(c6 / c3)
            # error: v errors and sigma0 error (assume 5 km/s minimum on sigma0, 10 % systematic on the correction)
            dsig = max(g.e_sigma0, 5.0)
            dt_v = tobs * np.sqrt((g.v_6Rd * g.e_v_6Rd / c6) ** 2 + (g.v_3Rd * g.e_v_3Rd / c3) ** 2)
            dt_s = tobs * abs(12 * g.sigma0 / c6 - 6 * g.sigma0 / c3) * dsig
            dt_sys = tobs * 0.5 * abs(12 * g.sigma0 ** 2 / c6 - 6 * g.sigma0 ** 2 / c3) * 0.3   # 30 % on the pressure term
            et = np.sqrt(dt_v ** 2 + dt_s ** 2 + dt_sys ** 2 + 0.03 ** 2)
        else:
            tobs = g.v_6Rd / g.v_3Rd
            et = np.sqrt(tobs ** 2 * ((g.e_v_6Rd / g.v_6Rd) ** 2 + (g.e_v_3Rd / g.v_3Rd) ** 2) + 0.03 ** 2)
        fg = g.mu_mol / (1 + g.mu_mol) if gas_ext != 1.0 else 0.0
        prof = np.zeros(len(LOGR_GRID))
        for j, lr in enumerate(LOGR_GRID):
            tm = np.array([model_slope(10 ** (g.logMbar + dm), g.Rd, 10 ** lr, fg)[0] for dm in dmg])
            prof[j] = (((tobs - tm) / et) ** 2 + (dmg / 0.3) ** 2).min()
        ok = prof <= prof.min() + 1
        rA = Hz(g.z)
        jA, jB = int(np.argmin(np.abs(LOGR_GRID - np.log10(rA)))), int(np.argmin(np.abs(LOGR_GRID)))
        tB, x3 = model_slope(10 ** g.logMbar, g.Rd, 1.0, fg)
        tA, _ = model_slope(10 ** g.logMbar, g.Rd, rA, fg)
        rows.append(dict(ID=g.ID, z=g.z, vsig=g.vsig, sigma0=g.sigma0, x3=x3, t_obs=tobs, e_t=et, t_B=tB, t_A=tA,
                         r=10 ** LOGR_GRID[int(np.argmin(prof))], r_lo=10 ** LOGR_GRID[ok].min(), r_hi=10 ** LOGR_GRID[ok].max(),
                         chi2_min=prof.min(), chi2_B=prof[jB], chi2_A=prof[jA]))
        prof_all += prof; chiA_all += prof[jA]; chiB_all += prof[jB]
    ok = prof_all <= prof_all.min() + 1
    comb = dict(N=len(rows), r=10 ** LOGR_GRID[int(np.argmin(prof_all))], r_lo=10 ** LOGR_GRID[ok].min(), r_hi=10 ** LOGR_GRID[ok].max(),
                chi2_min=prof_all.min(), chi2_B=chiB_all, chi2_A=chiA_all)
    return pd.DataFrame(rows), comb, prof_all


def main():
    t = load_kurvs()
    L = ["# a0 at high redshift, third pass: the shape of the outer rotation curves (2026-10-12)\n",
         "r = a0(z)/a0(0); H_A: r = H(z)/H0 (2.3 at z = 1.5, 6.9 at z = 4.4); H_B: r = 1. The shape test compares the measured "
         "outer slope of the rotation curve with the slope the law predicts for the galaxy's own baryonic disc, the baryonic mass "
         "being free with a 0.3 dex prior: it is insensitive to alpha_CO, M/L and the gas fraction at the 0.3 dex level.\n"]

    # (a) KURVS amplitude with f_DM
    ka = kurvs_amplitude(t)
    chi, out = chi2_nu(ka)
    lr, lo, hi, cmin = fit_r(chi)
    L.append("## KURVS, amplitude: f_DM(<Re) of the ten rotation-dominated discs (z ~ 1.5)\n")
    L.append(f"- g_bar from M* + molecular gas (Tacconi et al. 2018, median M_mol/M* = {t.mu_mol.median():.2f}), exponential disc; "
             f"g_bar/a0 at Re: 16-50-84 % = {np.percentile(10**ka.logg_bar/A0,16):.2f} / {np.percentile(10**ka.logg_bar/A0,50):.2f} / {np.percentile(10**ka.logg_bar/A0,84):.2f}; "
             f"nu_obs median {np.median(10**ka.lognu_obs):.2f}, law H_B {np.median(nu(10**ka.logg_bar/A0)):.2f}, H_A {np.median(nu(10**ka.logg_bar/A0/Hz(ka.z))):.2f}.")
    L.append(f"- r = {10**lr:.2f} ({10**lo:.2f}-{10**hi:.2f}); chi2(H_B) = {out['B'].sum():.1f}, chi2(H_A) = {out['A'].sum():.1f}, min {cmin:.1f}; "
             f"all masses +0.2 dex: r = {10**fit_r(chi2_nu(shift_mass(ka, 0.2))[0])[0]:.2f}; -0.2 dex: r = {10**fit_r(chi2_nu(shift_mass(ka, -0.2))[0])[0]:.2f}. "
             "Same degeneracy as the other z ~ 1-1.5 samples: not discriminating.\n")

    # (b) KURVS shape
    L.append("## KURVS, shape: outer slope t = V_c(6 Rd) / V_c(3 Rd)\n")
    for vmin, label in ((3.0, "v/sigma0 >= 3, pressure term small"), (1.5, "v/sigma0 >= 1.5, rotation-dominated")):
        for press in (True, False):
            sh, comb, prof = kurvs_shape(t, vmin, pressure=press)
            L.append(f"### {label}, {'with' if press else 'without'} the Burkert et al. 2010 pressure correction (N = {comb['N']})\n")
            L.append("| ID | z | v/sigma0 | sigma0 | g_bar(3 Rd)/a0 | t observed | t law H_B | t law H_A | r (68 %) | chi2 min | chi2 H_B | chi2 H_A |")
            L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
            for _, s in sh.iterrows():
                L.append(f"| {s.ID} | {s.z:.2f} | {s.vsig:.1f} | {s.sigma0:.0f} | {s.x3:.2f} | {s.t_obs:.2f} +/- {s.e_t:.2f} | {s.t_B:.2f} | {s.t_A:.2f} | "
                         f"{s.r:.2f} ({s.r_lo:.2f}-{s.r_hi:.2f}) | {s.chi2_min:.1f} | {s.chi2_B:.1f} | {s.chi2_A:.1f} |")
            L.append(f"\ncombined: r = {comb['r']:.2f} ({comb['r_lo']:.2f}-{comb['r_hi']:.2f}); chi2 min {comb['chi2_min']:.1f}, "
                     f"H_B {comb['chi2_B']:.1f} (Delta {comb['chi2_B']-comb['chi2_min']:+.1f}), H_A {comb['chi2_A']:.1f} (Delta {comb['chi2_A']-comb['chi2_min']:+.1f})\n")
    # (b') gas disc twice as extended as the stars (half the baryons at z ~ 1.5 are molecular gas)
    L.append("### Same, no pressure correction, molecular gas in a disc twice as extended as the stars\n")
    for vmin in (3.0, 1.5):
        sh, comb, _ = kurvs_shape(t, vmin, pressure=False, gas_ext=2.0)
        L.append(f"- v/sigma0 >= {vmin:.1f} (N = {comb['N']}): t law H_B median {sh.t_B.median():.3f}, H_A {sh.t_A.median():.3f}, observed median {sh.t_obs.median():.2f}; "
                 f"r = {comb['r']:.2f} ({comb['r_lo']:.2f}-{comb['r_hi']:.2f}); H_B Delta chi2 {comb['chi2_B']-comb['chi2_min']:+.1f}, H_A {comb['chi2_A']-comb['chi2_min']:+.1f}")
    L.append("- deep-MOND limit of the slope for a single exponential disc: t = [M(<6Rd)/M(<3Rd)]^(1/4) = 1.05; the law cannot produce "
             "t > 1.05 whatever a0, so slopes above it measure the pressure correction, the extrapolation of the fitted curve beyond "
             "R_Halpha,max (3-4 Rd for most KURVS discs) or an extended gas disc, not a0.\n")
    # (c) z~4.4 curves (from a0_highz_coldgas)
    sh4, comb4 = shape_test()
    L.append("## z ~ 4.4, shape of the four Roman-Oliveira et al. 2024 curves (mass free, 0.3 dex prior)\n")
    L.append(f"- combined: r = {comb4['r']:.2f} ({comb4['r_lo']:.2f}-{comb4['r_hi']:.2f}); chi2(H_B) = {comb4['chi2_B']:.1f}, chi2(H_A) = {comb4['chi2_A']:.1f}, "
             f"min = {sh4.chi2_min.sum():.1f}: H_A at Delta chi2 = {comb4['chi2_A']-sh4.chi2_min.sum():+.1f}, H_B at {comb4['chi2_B']-sh4.chi2_min.sum():+.1f}.\n")
    # (d) overall shape verdict
    shk, combk, profk = kurvs_shape(t, 1.5, pressure=True)
    L.append("## Shape test, all (KURVS v/sigma0 >= 1.5 with pressure correction + four z ~ 4.4 curves)\n")
    L.append(f"- KURVS: H_A - min = {combk['chi2_A']-combk['chi2_min']:+.1f}, H_B - min = {combk['chi2_B']-combk['chi2_min']:+.1f} (N = {combk['N']})")
    L.append(f"- z ~ 4.4: H_A - min = {comb4['chi2_A']-sh4.chi2_min.sum():+.1f}, H_B - min = {comb4['chi2_B']-sh4.chi2_min.sum():+.1f} (N = 4, 17 rings)")
    dA = combk['chi2_A'] - combk['chi2_min'] + comb4['chi2_A'] - sh4.chi2_min.sum()
    dB = combk['chi2_B'] - combk['chi2_min'] + comb4['chi2_B'] - sh4.chi2_min.sum()
    L.append(f"- total: H_A at Delta chi2 = {dA:+.1f}, H_B at {dB:+.1f} above the best shape fit.\n")
    L.append("Caveats: the KURVS slopes are H-alpha rotation velocities with v/sigma0 = 1.5-4.4; the pressure correction adds 6-12 sigma0^2 "
             "to V_c^2 and is itself uncertain by tens of percent (30 % assumed here), which is why the v/sigma0 >= 3 subset and the "
             "no-correction variant are shown; beam smearing and extrapolation to 6 Rd are handled by the authors' models. The z ~ 4.4 "
             "curves are [CII] circular speeds after asymmetric-drift correction, 3-5 rings each.")
    L.append("\nFigure: results/a0_highz_shape.png")
    open("results/a0_highz_shape.md", "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    a = ax[0]
    shp, _, _ = kurvs_shape(t, 1.5, pressure=True)
    a.errorbar(shp.x3, shp.t_obs, yerr=shp.e_t, fmt="o", mfc="none", color="C3", label="KURVS v/sigma0 >= 1.5, pressure-corrected (6-12 sigma0^2 added)")
    shn, _, _ = kurvs_shape(t, 1.5, pressure=False)
    a.errorbar(shn.x3, shn.t_obs, yerr=shn.e_t, fmt="o", color="C3", label="KURVS v/sigma0 >= 1.5, rotation velocities as published")
    Ms = np.logspace(9, 12, 200); Rd = 3.0
    for r, ls, lab in ((1.0, "--", "law H_B (r = 1)"), (Hz(1.5), "-", "law H_A at z = 1.5 (r = 2.3)")):
        global GAS_EXT
        GAS_EXT = 1.0
        vals = np.array([model_slope(M, Rd, r) for M in Ms])
        a.plot(vals[:, 1], vals[:, 0], "k" + ls, label=lab + ", stars only")
        GAS_EXT = 2.0
        vals = np.array([model_slope(M, Rd, r, 0.5) for M in Ms])
        a.plot(vals[:, 1], vals[:, 0], "C0" + ls, label=lab + ", half the baryons in a gas disc twice as extended")
    a.axhline(1.05, color="grey", lw=0.8); a.text(0.02, 1.06, "deep-MOND limit, single disc", fontsize=7, color="grey")
    a.set_xscale("log"); a.set_xlabel("g_bar(3 Rd) / a0"); a.set_ylabel("t = V(6 Rd) / V(3 Rd)"); a.set_ylim(0.6, 1.6); a.legend(fontsize=6.5, loc="upper right")
    a.set_title("KURVS z ~ 1.5: outer slope versus the law")
    b = ax[1]
    for lab, pr, col in (("KURVS shape", profk - profk.min(), "C3"), ("z~4.4 curves shape", None, "C4")):
        if pr is not None:
            b.plot(10 ** LOGR_GRID, pr, color=col, label=lab)
    b.axvline(1.0, color="k", ls="--", label="H_B"); b.axvline(Hz(1.5), color="k", ls="-", label="H_A at z = 1.5"); b.axvline(Hz(4.4), color="k", ls="-.", label="H_A at z = 4.4")
    b.set_xscale("log"); b.set_xlabel("r = a0(z)/a0(0)"); b.set_ylabel("Delta chi2 (mass marginalised)"); b.set_ylim(0, 40); b.legend(fontsize=8)
    b.set_title("Profile likelihood of the shape test")
    plt.tight_layout(); plt.savefig("results/a0_highz_shape.png", dpi=130)
    print("figure written")


if __name__ == "__main__":
    main()
