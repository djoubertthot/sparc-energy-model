# -*- coding: utf-8 -*-
"""
a0 at high redshift, fourth pass: the rotation-curve shape test on every cold-gas curve that could be recovered.

For each galaxy the circular-velocity curve V_c(R) is compared with the law V_c^2 = nu(g_bar/(r a0)) g_bar R, where g_bar comes
from the galaxy's own baryonic structure (sizes, Sersic index, bulge fraction from the papers) and the TOTAL baryonic mass is
free with a 0.3 dex Gaussian prior around the published value. The profile likelihood in r = a0(z)/a0(0) therefore reflects the
SHAPE of the curve (how fast nu rises outward), not its amplitude, and is insensitive to alpha_CO, M/L and the gas fraction.

Curves:
  Rizzo et al. 2021 (five lensed [CII] discs, z = 4.2-4.8) and Rizzo et al. 2020 (SPT0418-47, z = 4.22): V_rot(R) from their
    multi-parameter fits (V_t, R_t, beta, xi), asymmetric drift V_A^2 = sigma^2 R (1/R_gas + 2/R_sigma) with their exponential
    sigma(R) and gas disc; stars as a deprojected Sersic sphere (M_star, R_e, n), gas as an exponential disc (R_gas).
  Roman-Oliveira et al. 2024 (four [CII] discs, z = 4.3-4.5): circular speeds from Zenodo (see a0_highz_coldgas.py).
  Lelli et al. 2023 (two CO discs, z = 1.47 and 2.24): ring velocities read off their Figure (data/external/highz/lelli2023_rings_byeye.csv),
    bulge + disc from their surface-brightness fits, B/T from their baryons-only mass models.
  ALPAKA (Rizzo et al. 2023, 19 CO/[CI] discs, z = 0.56-3.63): 3DBarolo ring velocities read off their Fig. 6
    (data/external/highz/alpaka_rizzo2023_rings_byeye.csv), exponential disc with the stellar size from the mass-size relation
    (van der Wel et al. 2014, late types) and a gas disc 1.5 times larger, pressure correction with the mean sigma.
Outputs: results/a0_highz_shape2.md, results/a0_highz_shape2.png
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.special import gamma as Gamma, gammainc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a0_highz import A0, D, G_KPC, LOGR_GRID, Hz, g_si, nu, re_vdw, v2_expdisc, v2_hernquist  # noqa: E402
from a0_highz_coldgas import shape_test as shape_test_ro24  # noqa: E402

LN10 = np.log(10)
DM_GRID = np.arange(-0.8, 0.81, 0.02)
PRIOR_DEX = 0.3


def kpc_per_arcsec(z):
    """Flat LCDM, Omega_m = 0.3, H0 = 70."""
    from scipy.integrate import quad
    dc = 299792.458 / 70.0 * quad(lambda zz: 1.0 / Hz(zz), 0, z)[0]   # Mpc
    return dc / (1 + z) * np.pi / 180 / 3600 * 1000.0                   # kpc per arcsec


def v2_sersic_sphere(R, M, Re, n):
    """Circular velocity squared of a deprojected Sersic sphere (Prugniel & Simien 1997)."""
    b = 2 * n - 1.0 / 3 + 0.009876 / n
    p = 1.0 - 0.6097 / n + 0.05463 / n ** 2
    a = n * (3 - p)
    Menc = M * gammainc(a, b * (np.asarray(R, float) / Re) ** (1.0 / n))
    return G_KPC * Menc / R


def law_v2(v2bar, R, r):
    gb = g_si(v2bar, R)
    return nu(gb / (r * A0)) * v2bar


def profile(R, Vobs, eV, v2bar_fun, logM0, label, z, size_factors=(1.0,)):
    """Profile likelihood in log r with the baryonic mass marginalised (Gaussian prior PRIOR_DEX around logM0),
    optionally also over a set of size scalings (flat prior) when the structure is uncertain."""
    prof = np.zeros(len(LOGR_GRID)); best_dm = np.zeros(len(LOGR_GRID))
    for j, lr in enumerate(LOGR_GRID):
        chis = []
        for dm in DM_GRID:
            best = np.inf
            for sf in size_factors:
                v2b = v2bar_fun(R, 10 ** (logM0 + dm), sf) if len(size_factors) > 1 else v2bar_fun(R, 10 ** (logM0 + dm))
                vmod = np.sqrt(law_v2(v2b, R, 10 ** lr))
                best = min(best, (((Vobs - vmod) / eV) ** 2).sum())
            chis.append(best + (dm / PRIOR_DEX) ** 2)
        chis = np.array(chis); k = int(np.argmin(chis))
        prof[j], best_dm[j] = chis[k], DM_GRID[k]
    i = int(np.argmin(prof)); ok = prof <= prof[i] + 1
    rA = Hz(z); jA = int(np.argmin(np.abs(LOGR_GRID - np.log10(rA)))); jB = int(np.argmin(np.abs(LOGR_GRID)))
    v2_0 = v2bar_fun(R, 10 ** logM0, 1.0) if len(size_factors) > 1 else v2bar_fun(R, 10 ** logM0)
    x_out = g_si(v2_0, R)[-1] / A0
    return dict(name=label, z=z, N=len(R), R_max=float(np.max(R)), x_out=float(x_out), r=10 ** LOGR_GRID[i],
                r_lo=10 ** LOGR_GRID[ok].min(), r_hi=10 ** LOGR_GRID[ok].max(), dlogM=best_dm[i], chi2_min=prof[i],
                chi2_B=prof[jB], chi2_A=prof[jA], rA=rA), prof


# ------------------------------------------------------------------ Rizzo 2020/2021 analytic curves
RIZZO = [  # name, z, V_t, eV_t, R_t, beta, xi, sigma0, e_sigma0, R_sigma, M_star(1e10), Re, n, M_gas(1e10), R_gas, R_max, form
    ("SPT0113-46", 4.23, 352, 8, 0.020, 0.79, 2.3, 145, 12, 1.5, 6.0, 1.4, 6.0, 4.3, 1.5, 3.4, "multi"),
    ("SPT0345-47", 4.29, 200, 25, 0.68, 1.2, 2.8, 123, 13, 2.2, 2.3, 0.30, 1.5, 1.7, 0.76, 3.1, "multi"),
    ("SPT0441-46", 4.48, 293, 7, 0.17, 1.1, 2.9, 151, 23, 0.68, 1.8, 0.11, 2.0, 1.4, 0.63, 2.0, "multi"),
    ("SPT2146-55", 4.57, 176, 18, 0.34, 0.68, 2.1, 76, 10, 2.0, 1.0, 0.39, 1.6, 1.2, 1.4, 2.7, "multi"),
    ("SPT2132-58", 4.77, 219, 22, 0.48, np.nan, np.nan, 50, 7, 2.6, 1.9, 1.4, 0.94, 1.9, 2.1, 4.1, "arctan"),
    ("SPT0418-47", 4.22, 245.1, 0.8, 0.14, 0.80, 2.0, 58.1, 1.8, 1.7, 1.2, 0.22, 2.2, 1.3, 0.65, 2.5, "multi"),
]


def rizzo_curve(p, R):
    name, z, Vt, eVt, Rt, beta, xi, s0, es0, Rs, Ms, Re, n, Mg, Rg, Rmax, form = p
    if form == "multi":
        vrot = Vt * (1 + Rt / R) ** beta / (1 + (Rt / R) ** xi) ** (1 / xi)
    else:
        vrot = 2 / np.pi * Vt * np.arctan(R / Rt)
    sig = s0 * np.exp(-R / Rs)
    va2 = sig ** 2 * R * (1 / Rg + 2 / Rs)
    vc = np.sqrt(vrot ** 2 + va2)
    # error: V_t fractional, sigma0 fractional on the drift term, 10 % floor
    evc = np.sqrt((vrot * eVt / Vt) ** 2 + (0.5 * va2 / vc * 2 * es0 / s0) ** 2 + (0.1 * vc) ** 2)
    return vc, evc


def run_rizzo():
    out, profs = [], {}
    for p in RIZZO:
        name, z, *_ = p
        Ms, Re, n, Mg, Rg, Rmax = p[10], p[11], p[12], p[13], p[14], p[15]
        R = np.linspace(0.4 * Rg, Rmax, 6)
        vc, evc = rizzo_curve(p, R)
        M0 = (Ms + Mg) * 1e10; fs = Ms / (Ms + Mg)

        def v2bar(RR, M, fs=fs, Re=Re, n=n, Rg=Rg):
            return v2_sersic_sphere(RR, M * fs, Re, n) + v2_expdisc(RR, M * (1 - fs), Rg)
        res, prof = profile(R, vc, evc, v2bar, np.log10(M0), name, z)
        res["sample"] = "Rizzo 2020/2021 lensed [CII]"
        out.append(res); profs[name] = prof
    return out, profs


# ------------------------------------------------------------------ Lelli 2023 rings
LELLI = {"zC-400569": dict(z=2.24, Re_bul_arcsec=0.12, Re_disk_arcsec=0.76, BT=0.20, Mbar=12.8e10),
         "zC-488879": dict(z=1.47, Re_bul_arcsec=0.35, Re_disk_arcsec=0.88, BT=0.54, Mbar=26.8e10)}


def run_lelli():
    t = pd.read_csv(f"{D}/lelli2023_rings_byeye.csv")
    out, profs = [], {}
    for name, g in LELLI.items():
        s = t[(t.galaxy == name) & (t.line.str.startswith("CO"))]
        k = kpc_per_arcsec(g["z"])
        R, V, eV = s.R_arcsec.values * k, s.Vrot.values, s.e_Vrot.values * np.sqrt(2)   # two CO lines of the same curve: inflate
        Reb, Red, BT = g["Re_bul_arcsec"] * k, g["Re_disk_arcsec"] * k, g["BT"]

        def v2bar(RR, M, Reb=Reb, Red=Red, BT=BT):
            return v2_hernquist(RR, M * BT, Reb) + v2_expdisc(RR, M * (1 - BT), Red / 1.678)
        res, prof = profile(R, V, eV, v2bar, np.log10(g["Mbar"]), name, g["z"])
        res["sample"] = "Lelli 2023 CO"
        out.append(res); profs[name] = prof
    return out, profs


# ------------------------------------------------------------------ ALPAKA rings
def run_alpaka():
    rings = pd.read_csv(f"{D}/alpaka_rizzo2023_rings_byeye.csv")
    props = {1: (10.04, 10), 2: (10.76, 15), 3: (10.51, 49), 6: (10.76, 42), 7: (10.88, 52), 8: (11.04, 44), 9: (10.59, 30),
             11: (10.80, 30), 12: (10.77, 30), 13: (11.08, 32), 15: (11.40, 99), 18: (11.08, 40), 19: (10.90, 42), 20: (11.51, 38),
             22: (10.85, 52), 23: (11.72, 24), 24: (11.0, 30), 25: (11.04, 48), 28: (11.11, 100)}   # log M* (Table 3), sigma_m (Table 5)
    out, profs = [], {}
    for gid, s in rings.groupby("ID"):
        z = float(s.z.iloc[0]); k = kpc_per_arcsec(z)
        R, V, eV = s.R_arcsec.values * k, s.Vrot.values.astype(float), s.e_Vrot.values.astype(float)
        logM, sig = props[int(gid)]
        Re = re_vdw(z, logM); Rd = Re / 1.678
        vc = np.sqrt(V ** 2 + 2 * sig ** 2 * R / (1.5 * Rd))    # pressure term with the gas disc (1.5 Rd)
        evc = np.sqrt(eV ** 2 + (0.5 * 2 * sig ** 2 * R / (1.5 * Rd) / vc) ** 2)
        M0 = 10 ** logM * 1.8   # gas fraction ~ 45 % (molecular, Tacconi et al. 2018 at these masses and z); free anyway

        def v2bar(RR, M, sf, Rd=Rd):   # size free: both discs scaled by sf (0.25 to 2 times the mass-size relation)
            return v2_expdisc(RR, 0.55 * M, sf * Rd) + v2_expdisc(RR, 0.45 * M, 1.5 * sf * Rd)
        res, prof = profile(R, vc, evc, v2bar, np.log10(M0), f"ALPAKA ID{int(gid)}", z, size_factors=(0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0))
        res["sample"] = "ALPAKA CO/[CI]"
        out.append(res); profs[f"ID{int(gid)}"] = prof
    return out, profs


# ------------------------------------------------------------------ Roman-Oliveira 2024 curves (same machinery)
RO24 = {"BRI1335-0417": (4.41, 10.4, 9.8, 2.3, 2.9), "J081740": (4.54, 10.6, 10.5, 2.4, 2.8),
        "SGP38326-1": (4.28, 11.0, 11.3, 2.2, 2.8), "SGP38326-2": (4.28, 10.3, 11.1, 1.7, 2.0)}


def run_ro24():
    cur = pd.read_csv(f"{D}/z45/circular_speed_romanoliveira2024.csv")
    out, profs = [], {}
    for name, (z, lms, lmg, res_, reb) in RO24.items():
        c = cur[cur.galaxy == name]
        R, V, eV = c.R_kpc.values, c.Vcirc_kms.values, np.maximum(c.e_Vcirc_kms.values, 0.1 * c.Vcirc_kms.values)
        fs = 10 ** lms / (10 ** lms + 10 ** lmg)

        def v2bar(RR, M, fs=fs, res_=res_, reb=reb):
            return v2_hernquist(RR, M * fs, res_) + v2_expdisc(RR, M * (1 - fs), reb / 1.678)
        res, prof = profile(R, V, eV, v2bar, np.log10(10 ** lms + 10 ** lmg), name, z)
        res["sample"] = "Roman-Oliveira 2024 [CII]"
        out.append(res); profs[name] = prof
    return out, profs


def main():
    rz, pz = run_rizzo(); ll, pl = run_lelli(); al, pa = run_alpaka(); ro, pr = run_ro24()
    T = pd.DataFrame(rz + ll + ro + al)
    sh4 = pd.DataFrame(ro)
    L = ["# a0 at high redshift, fourth pass: shape test on all recovered cold-gas curves (2026-10-12)\n",
         "Per galaxy: circular-velocity curve against the law with the galaxy's own baryonic structure and the total baryonic mass free "
         "(0.3 dex prior). r = a0(z)/a0(0); H_A: r = H(z)/H0; H_B: r = 1. 'x_out' = g_bar/a0 at the last radius with the published mass. "
         "Delta chi2 columns are relative to the best r. Curves: Rizzo et al. 2020/2021 (analytic from their fitted parameters, "
         "with asymmetric drift), Lelli et al. 2023 and ALPAKA (ring velocities read off the published figures, see data/external/highz/README.md), "
         "Roman-Oliveira et al. 2024 (Zenodo).\n",
         "| galaxy | sample | z | N | R_max (kpc) | x_out | r (68 %) | best dlog M | chi2 min | H_B - min | H_A - min | r_A |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, s in T.iterrows():
        L.append(f"| {s['name']} | {s['sample']} | {s.z:.2f} | {s.N} | {s.R_max:.1f} | {s.x_out:.1f} | {s.r:.2f} ({s.r_lo:.2f}-{s.r_hi:.2f}) | "
                 f"{s.dlogM:+.2f} | {s.chi2_min:.1f} | {s.chi2_B - s.chi2_min:+.1f} | {s.chi2_A - s.chi2_min:+.1f} | {s.rA:.1f} |")
    L.append("")
    L.append("## Sums of Delta chi2 (relative to each galaxy's best r)\n")
    groups = {"Rizzo 2020/2021 lensed [CII], z = 4.2-4.8": T[T['sample'].str.startswith("Rizzo")],
              "Roman-Oliveira 2024 [CII], z = 4.3-4.5": T[T['sample'].str.startswith("Roman")],
              "Lelli 2023 CO, z = 1.5 and 2.2": T[T['sample'].str.startswith("Lelli")],
              "ALPAKA CO/[CI], z = 0.6-3.6 (sizes free, rings read off a small figure)": T[T['sample'].str.startswith("ALPAKA")]}
    for k, g in groups.items():
        dA, dB = (g.chi2_A - g.chi2_min).sum(), (g.chi2_B - g.chi2_min).sum()
        L.append(f"- {k}: N = {len(g)}, chi2_min/N_points = {g.chi2_min.sum()/g.N.sum():.2f}, sum(H_B - min) = {dB:+.1f}, sum(H_A - min) = {dA:+.1f}")
    clean = T[~T['sample'].str.startswith("ALPAKA")]
    L.append(f"- clean set (Rizzo + Roman-Oliveira + Lelli, N = {len(clean)}, z = 1.5-4.8): sum(H_B - min) = {(clean.chi2_B - clean.chi2_min).sum():+.1f}, "
             f"sum(H_A - min) = {(clean.chi2_A - clean.chi2_min).sum():+.1f}, difference A - B = {(clean.chi2_A - clean.chi2_B).sum():+.1f}")
    hz = clean[clean.z >= 4.0]
    L.append(f"- clean set, z > 4 only (N = {len(hz)}, H_A needs r = 6.6-7.6): sum(H_B - min) = {(hz.chi2_B - hz.chi2_min).sum():+.1f}, "
             f"sum(H_A - min) = {(hz.chi2_A - hz.chi2_min).sum():+.1f}")
    # global fit of a single scaling r = (H/H0)^n across z, clean set and all
    ns = np.arange(-1.0, 2.01, 0.05)

    def nfit(profs, zs):
        cn = []
        for n in ns:
            tot = 0.0
            for prof, z in zip(profs, zs):
                j = int(np.argmin(np.abs(LOGR_GRID - n * np.log10(Hz(z))))); tot += prof[j]
            cn.append(tot)
        cn = np.array(cn); i = int(np.argmin(cn)); ok = cn <= cn[i] + 1
        return ns[i], ns[ok].min(), ns[ok].max(), cn[int(np.argmin(np.abs(ns)))] - cn[i], cn[int(np.argmin(np.abs(ns - 1)))] - cn[i]
    zs_clean = [p[1] for p in RIZZO] + [LELLI[k]["z"] for k in LELLI] + [RO24[k][0] for k in RO24]
    profs_clean = list(pz.values()) + list(pl.values()) + list(pr.values())
    n0, nlo, nhi, d0, d1 = nfit(profs_clean, zs_clean)
    L.append(f"- a0 ~ H(z)^n on the clean set ({len(profs_clean)} galaxies): n = {n0:.2f} ({nlo:.2f} to {nhi:.2f}); chi2(n = 0) - min = {d0:+.1f}, chi2(n = 1) - min = {d1:+.1f}")
    zs_al = [float(pd.read_csv(f"{D}/alpaka_rizzo2023_rings_byeye.csv").groupby("ID").z.first()[int(k[2:])]) for k in pa]
    n0, nlo, nhi, d0, d1 = nfit(profs_clean + list(pa.values()), zs_clean + zs_al)
    L.append(f"- same with ALPAKA added ({len(profs_clean) + len(pa)} galaxies): n = {n0:.2f} ({nlo:.2f} to {nhi:.2f}); chi2(n = 0) - min = {d0:+.1f}, chi2(n = 1) - min = {d1:+.1f}")
    L.append("")
    L.append("## Caveats\n")
    L.append("- Rizzo et al. curves are parametric fits (V_t, R_t, beta, xi) to lensed [CII] cubes at 20-300 pc resolution; six radii per "
             "galaxy are sampled from the smooth curve with a 10 % floor, so the chi2 values are indicative, the r profiles are what matters. "
             "Their asymmetric-drift term uses sigma0 = 50-150 km/s and is large for SPT0113-46 and SPT0441-46.")
    L.append("- ALPAKA rings were read off a 937 x 1020 pixel figure (2-5 rings per galaxy, 10-15 % errors), sizes from the mass-size relation: "
             "weak individual constraints, kept for completeness; x_out is large for most (Newtonian regime), so H_A and H_B differ little.")
    L.append("- Lelli et al. 2023 rings were read off their Figure 7 (CO lines only); structure from their HST fits.")
    L.append("- The law's maximum outer slope for a pure exponential disc is +5 % between 3 and 6 Rd; rising curves beyond that are not a0.")
    L.append("\nFigure: results/a0_highz_shape2.png")
    open("results/a0_highz_shape2.md", "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    a = ax[0]
    for name, prof in pz.items():
        a.plot(10 ** LOGR_GRID, prof - prof.min(), label=name, color="C4", alpha=0.8)
    for name, prof in pl.items():
        a.plot(10 ** LOGR_GRID, prof - prof.min(), label=name, color="C1", ls="--")
    for name, prof in pr.items():
        a.plot(10 ** LOGR_GRID, prof - prof.min(), label=name, color="C2", ls=":")
    tot = sum(list(pz.values()) + list(pl.values()) + list(pr.values()))
    a.plot(10 ** LOGR_GRID, tot - tot.min(), color="k", lw=2.5, label="clean set, sum (12 galaxies)")
    a.axvline(1, color="k", ls="--"); a.axvline(Hz(4.5), color="k", ls="-."); a.axvline(Hz(2.2), color="k", ls=":")
    a.text(1.05, 18, "H_B", fontsize=8); a.text(Hz(4.5) * 1.05, 18, "H_A z=4.5", fontsize=8); a.text(Hz(2.2) * 1.05, 16, "H_A z=2.2", fontsize=8)
    a.set_xscale("log"); a.set_ylim(0, 20); a.set_xlabel("r = a0(z)/a0(0)"); a.set_ylabel("Delta chi2, mass marginalised"); a.legend(fontsize=7)
    a.set_title("Shape profiles: Rizzo 2020/2021 (z~4.5), Roman-Oliveira 2024 (z~4.4), Lelli 2023 (z 1.5, 2.2)")
    b = ax[1]
    # example: SPT2132-58 curve with the two predictions at the published mass
    p = [q for q in RIZZO if q[0] == "SPT2132-58"][0]
    R = np.linspace(0.3, 4.1, 50); vc, evc = rizzo_curve(p, R)
    Ms, Re, n, Mg, Rg = p[10], p[11], p[12], p[13], p[14]
    v2b = v2_sersic_sphere(R, Ms * 1e10, Re, n) + v2_expdisc(R, Mg * 1e10, Rg)
    b.fill_between(R, vc - evc, vc + evc, color="C4", alpha=0.2); b.plot(R, vc, color="C4", label="SPT2132-58 V_c (Rizzo+21)")
    b.plot(R, np.sqrt(v2b), "k:", label="baryons, Newton")
    b.plot(R, np.sqrt(law_v2(v2b, R, 1.0)), "k--", label="law, a0 constant")
    b.plot(R, np.sqrt(law_v2(v2b, R, Hz(4.77))), "k-", label="law, a0 x %.1f (H_A)" % Hz(4.77))
    b.set_xlabel("R (kpc)"); b.set_ylabel("V (km/s)"); b.legend(fontsize=8); b.set_title("Example at z = 4.77, published baryonic mass")
    plt.tight_layout(); plt.savefig("results/a0_highz_shape2.png", dpi=130)
    print("figure written")


if __name__ == "__main__":
    main()
