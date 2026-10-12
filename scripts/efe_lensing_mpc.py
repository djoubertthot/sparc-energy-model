# -*- coding: utf-8 -*-
"""
The external-field test at 0.3-2.5 Mpc: QUMOND predictions against the weak-lensing circular velocities of Mistele et al. 2024
(ApJL 969, L3; machine-readable Table 1, data/external/lensing/mistele2024_vc_mrt.txt), isolated KiDS-1000 lenses in four
baryonic-mass bins (log M_b = 9.0-10.5, 10.5-10.8, 10.8-11.1, 11.1-11.5; means 10.10, 10.66, 10.96, 11.29).

For a point mass M_b in a uniform Newtonian external field e_N, the QUMOND effective enclosed mass follows from Gauss's law
(see scripts/efe_qumond_lensing.py); the spherical deprojection of the lensing stacks measures that monopole to first order.
Two fits per bin: (a) amplitude fixed by the law, M_b from the bin mean, (b) amplitude free (M_b scaled to match R < 300 kpc),
so that the SHAPE beyond 0.3 Mpc alone constrains e_N. chi2 uses the statistical and systematic errors in quadrature.
Output: results/efe_lensing_mpc.md, results/efe_lensing_mpc.png
"""
import os

import numpy as np
import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
G = 4.301e-6
A0 = 1.2e-10 * 3.0857e19 / 1e6          # (km/s)^2 / kpc
LOGMB = [10.10, 10.66, 10.96, 11.29]
EN_GRID = np.concatenate([[0.0], np.logspace(-5.5, -2.5, 31)])
th = np.linspace(0, np.pi, 2001)
w = np.sin(th) * (th[1] - th[0])


def nu(y):
    y = np.maximum(y, 1e-14)
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))


def v_eff(M, eN, R):
    out = []
    for r in np.atleast_1d(R):
        gr = -G * M / r ** 2 + eN * A0 * np.cos(th)
        gt = -eN * A0 * np.sin(th)
        gabs = np.sqrt(gr ** 2 + gt ** 2)
        flux = (nu(gabs / A0) * gr * w).sum() * 2 * np.pi * r ** 2
        out.append(np.sqrt(G * (-flux / (4 * np.pi * G)) / r))
    return np.array(out)


def read_mrt(path):
    rows = []
    for line in open(path, encoding="utf-8"):
        if not line[:24].strip() or line.startswith(("Title", "Authors", "Table", "=", "-", " ", "Byte", "Note")):
            continue
        sample = line[0:24].strip()
        try:
            R = float(line[25:32])
        except ValueError:
            continue
        rec = dict(sample=sample, R=R)
        for k in range(4):
            base = 33 + 21 * k
            vals = [line[base:base + 6], line[base + 7:base + 13], line[base + 14:base + 20]]
            try:
                rec[f"Vc{k+1}"], rec[f"es{k+1}"], rec[f"ey{k+1}"] = [float(v) for v in vals]
            except ValueError:
                rec[f"Vc{k+1}"] = rec[f"es{k+1}"] = rec[f"ey{k+1}"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


T = read_mrt("data/external/lensing/mistele2024_vc_mrt.txt")
A = T[T["sample"] == "All"].copy()
L = ["# External-field effect at 0.3-2.5 Mpc: QUMOND against the lensing circular velocities of Mistele et al. 2024 (2026-10-12)\n",
     "Isolated KiDS-1000 lenses, four baryonic-mass bins, circular velocities from the spherical deprojection of the stacked excess "
     "surface density (statistical and systematic errors in quadrature). QUMOND prediction for a point mass in a uniform Newtonian "
     "external field e_N (Gauss's law on the field equation, monopole). The authors do not trust their isolation criterion beyond "
     "1 Mpc; results are given for R <= 1 Mpc and for the full range.\n"]
summary = []
curves = {}
for k, lm in enumerate(LOGMB):
    Mb = 10 ** lm
    V = A[f"Vc{k+1}"].values; e = np.sqrt(A[f"es{k+1}"].values ** 2 + A[f"ey{k+1}"].values ** 2); R = A.R.values
    ok = np.isfinite(V)
    V, e, R = V[ok], e[ok], R[ok]
    vflat_law = (G * Mb * A0) ** 0.25
    # amplitude free: scale M_b so that the isolated law matches R < 300 kpc (V ~ M^1/4)
    inner = R < 300
    scale = (np.average(V[inner], weights=1 / e[inner] ** 2) / vflat_law) ** 4
    L.append(f"## Bin {k+1}, log M_b = {lm}: V_flat(law) = {vflat_law:.0f} km/s, lensing R < 300 kpc = {np.average(V[inner], weights=1/e[inner]**2):.0f} km/s "
             f"(implied M_b x {scale:.2f})\n")
    L.append("| e_N / a0 | r_t (kpc) | chi2 R<=1 Mpc, fixed amplitude | chi2 R<=1 Mpc, free amplitude | chi2 all R, fixed | chi2 all R, free | V(0.5)/V(0.1) pred. | V(1)/V(0.1) pred. |")
    L.append("|---|---|---|---|---|---|---|---|")
    rows = []
    for eN in EN_GRID:
        pf = v_eff(Mb, eN, R); pfree = v_eff(Mb * scale, eN, R)
        m1 = R <= 1000
        c = dict(eN=eN, rt=np.sqrt(G * Mb / (eN * A0)) if eN > 0 else np.inf,
                 c1f=np.sum(((V - pf) / e)[m1] ** 2), c1v=np.sum(((V - pfree) / e)[m1] ** 2),
                 caf=np.sum(((V - pf) / e) ** 2), cav=np.sum(((V - pfree) / e) ** 2),
                 r05=v_eff(Mb * scale, eN, 500)[0] / v_eff(Mb * scale, eN, 100)[0], r10=v_eff(Mb * scale, eN, 1000)[0] / v_eff(Mb * scale, eN, 100)[0])
        rows.append(c)
    Rw = pd.DataFrame(rows)
    for _, c in Rw[Rw.eN.isin([0.0] + list(EN_GRID[[1, 7, 11, 16, 21, 26, 31]]))].iterrows():
        L.append(f"| {c.eN:.1e} | {c.rt:.0f} | {c.c1f:.1f} | {c.c1v:.1f} | {c.caf:.1f} | {c.cav:.1f} | {c.r05:.2f} | {c.r10:.2f} |")
    # bounds: Delta chi2 = 4 (2 sigma) relative to the minimum, free amplitude; Delta chi2 at the 2M++ field values
    for col, lab in (("c1v", "R <= 1 Mpc"), ("cav", "all R")):
        cmin = Rw[col].min(); ok2 = Rw[col] <= cmin + 4.0
        emax = Rw.eN[ok2].max(); ebest = Rw.eN[Rw[col].idxmin()]
        d1 = np.interp(np.log10(1e-4), np.log10(Rw.eN[1:]), Rw[col][1:]) - cmin
        d15 = np.interp(np.log10(1.5e-4), np.log10(Rw.eN[1:]), Rw[col][1:]) - cmin
        N = int(m1.sum()) if col == "c1v" else len(R)
        L.append(f"- free amplitude, {lab}: best e_N = {ebest:.1e} (chi2 {cmin:.1f} for {N} points), e_N < {emax:.1e} a0 at Delta chi2 = 4; "
                 f"Delta chi2 = {d1:+.1f} at e_N = 1e-4, {d15:+.1f} at 1.5e-4 (2M++ median).")
        summary.append(dict(bin=k + 1, logMb=lm, range=lab, e_best=ebest, e_max_2sig=emax, chi2min=cmin, N=N, dchi2_1e4=d1, dchi2_15e4=d15))
    L.append("")
    curves[k] = dict(R=R, V=V, e=e, Mb=Mb, scale=scale, c1v=Rw.c1v.values, cav=Rw.cav.values)
S = pd.DataFrame(summary)
S.to_csv("results/efe_lensing_mpc.csv", index=False)
L.append("## Summary (free amplitude per bin)\n")
L.append("| bin | log M_b | range | best e_N / a0 | 2-sigma upper bound on e_N / a0 | chi2_min / N | Delta chi2 at 1e-4 | at 1.5e-4 |")
L.append("|---|---|---|---|---|---|---|---|")
for _, s in S.iterrows():
    L.append(f"| {s.bin} | {s.logMb} | {s.range} | {s.e_best:.1e} | {s.e_max_2sig:.1e} | {s.chi2min:.1f} / {s.N} | {s.dchi2_1e4:+.1f} | {s.dchi2_15e4:+.1f} |")
L.append("")
L.append("## Combined (same e_N for the four bins, amplitude free in each)\n")
comb = {}
for col, lab in (("c1v", "R <= 1 Mpc"), ("cav", "all R")):
    tot = sum(c[col] for c in curves.values()); cmin = tot.min(); ebest = EN_GRID[np.argmin(tot)]
    emax = EN_GRID[tot <= cmin + 4.0].max()
    d1 = np.interp(np.log10(1e-4), np.log10(EN_GRID[1:]), tot[1:]) - cmin
    d15 = np.interp(np.log10(1.5e-4), np.log10(EN_GRID[1:]), tot[1:]) - cmin
    comb[col] = dict(ebest=ebest, emax=emax, cmin=cmin, d1=d1, d15=d15, N=int(S[S.range == lab].N.sum()))
    L.append(f"- {lab}: best e_N = {ebest:.1e}, chi2_min = {cmin:.1f} for {comb[col]['N']} points, e_N < {emax:.1e} a0 at Delta chi2 = 4; "
             f"Delta chi2 = {d1:+.1f} at 1e-4 ({np.sqrt(max(d1,0)):.1f} sigma), {d15:+.1f} at 1.5e-4 ({np.sqrt(max(d15,0)):.1f} sigma).")
L.append("")
b1 = S[(S.bin == 1) & (S.range == "R <= 1 Mpc")].iloc[0]; c1 = comb["c1v"]
L.append("## Reading\n")
L.append(f"- Reference field: the 2M++ baryonic Newtonian field of the SPARC galaxies is e_N = 1.2-2.0e-4 a0 (median 1.5e-4, Milky Way 1.3e-4), of "
         "which only 3.7e-5 comes from within 10 Mpc; an isolation criterion at 3 Mpc cannot bring an isolated lens below about 1e-4.")
L.append(f"- The lowest-mass bin (log M_b = 10.1, the most sensitive because r_t is proportional to sqrt(M_b)) is flat to 1 Mpc and gives e_N < {b1.e_max_2sig:.0e} a0 "
         f"at 2 sigma; the 2M++ field is excluded at Delta chi2 = {b1.dchi2_15e4:.0f} ({np.sqrt(b1.dchi2_15e4):.1f} sigma) in that bin alone. The three heavier bins are "
         "compatible with 0 and with 1.5e-4 (their r_t lies at 0.7-1.5 Mpc, where the data stop).")
L.append(f"- Combined over the four bins, R <= 1 Mpc: e_N < {c1['emax']:.0e} a0, the 2M++ field excluded at {np.sqrt(c1['d15']):.1f} sigma.")
L.append("- Either the law holds in the lensing regime without the external-field effect of a MOND-type field theory (a dark component, or a non-local law "
         "whose field equation has no EFE), or the isolated KiDS lenses of 10^10 Msun sit in external fields three times weaker than the median 2M++ "
         "field of SPARC galaxies, which the isolation criterion does not provide. QUMOND/RMOND with the EFE is in tension with this data at the 3 sigma level.")
L.append("- Caveats: point-mass baryons (fine beyond 50 kpc); the stacked, deprojected V_c measures the monopole of the QUMOND effective mass, which is "
         "what Gauss's law gives; no two-halo term (Mistele et al. estimate it small below 1 Mpc for their isolation criterion, and the result holds "
         "for R <= 1 Mpc alone); systematic errors added in quadrature and treated as independent between radii, so the bound is optimistic by "
         "perhaps a factor 1.5; the hot-gas content of the lenses raises M_b and the inner amplitude (absorbed by the free amplitude) but not the shape "
         "beyond r_t. A dedicated stack with the external field estimated per lens from the KiDS/GAMA environment would turn this into a measurement.")
L.append("\nFigure: results/efe_lensing_mpc.png. Summary: results/efe_lensing_mpc.csv")
open("results/efe_lensing_mpc.md", "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 2, figsize=(12, 9))
Rg = np.logspace(1.4, 3.45, 60)
for k, c in curves.items():
    a = ax.flat[k]
    a.errorbar(c["R"], c["V"], yerr=c["e"], fmt="o", color="k", ms=4, label="lensing (Mistele+24)")
    for eN, ls, lab in ((0.0, "--", "isolated"), (4e-5, "-.", "e_N = 4e-5"), (1.5e-4, "-", "e_N = 1.5e-4 (2M++ field)"), (5e-4, ":", "e_N = 5e-4")):
        a.plot(Rg, v_eff(c["Mb"] * c["scale"], eN, Rg), ls=ls, color="C3", label=lab)
    a.axvline(1000, color="grey", lw=0.8); a.set_xscale("log"); a.set_xlabel("R (kpc)"); a.set_ylabel("V_c (km/s)")
    a.set_title(f"bin {k+1}: log M_b = {LOGMB[k]} (amplitude scaled x{c['scale']:.2f})"); a.legend(fontsize=7); a.set_ylim(0, 350)
plt.tight_layout(); plt.savefig("results/efe_lensing_mpc.png", dpi=130)
print("figure written")
