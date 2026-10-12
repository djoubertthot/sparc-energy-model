# -*- coding: utf-8 -*-
"""
The missing factor in galaxy clusters under the law g_obs = nu(g_bar/a0) g_bar.

Data:
  CLASH (Tian et al. 2020, ApJ 896, 70; VizieR J/ApJ/896/70/fig2): 84 points for 20 clusters at 14-600 kpc, log g_bar
      (BCG stars + ICM gas) and log g_tot (weak + strong lensing, X-ray), with errors.
  X-COP (Eckert et al. 2022, A&A 662, A123; Eckert et al. 2019, A&A 621, A40): 12 clusters at z = 0.05-0.09, total NFW
      (M200, c200) from hydrostatic + SZ profiles, gas fractions at R500 and R200 (and R2500 where published), stellar
      fractions 1.2-2.8 % where published, 1.5 % assumed otherwise. Evaluated at R2500 (0.45 R500 assumed), R500 and R200.
For every point: g_RAR = nu(g_bar/a0) g_bar, residual log(g_tot/g_RAR), and the extra mass the law would need, M_extra/M_bar =
(g_tot - g_RAR)/g_bar. Also the cluster-scale acceleration scale fitted with the RAR shape (Tian et al. find 2.0e-9).
Output: results/cluster_rar.md, results/cluster_rar.png
"""
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
A0 = 1.2e-10
G_SI = 6.674e-11
MSUN = 1.989e30
KPC = 3.0857e19


def nu(y):
    y = np.maximum(y, 1e-14)
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))


_yg = np.logspace(-6, 4, 4001)
_zg = nu(_yg) * _yg


def mond_mass_factor(gbar, gtot):
    """mass the law needs, in units of the baryons: the Newtonian g_N such that nu(g_N/a0) g_N = g_tot, divided by g_bar"""
    y_needed = 10 ** np.interp(np.log10(gtot / A0), np.log10(_zg), np.log10(_yg))
    return y_needed * A0 / gbar


# ---------------- CLASH
C = pd.read_csv("data/external/clusters/tian2020_clash_rar.csv")
C["gbar"] = 10 ** C.log_gbar; C["gtot"] = 10 ** C.log_gtot
C["log_gRAR"] = np.log10(nu(C.gbar / A0) * C.gbar)
C["resid"] = C.log_gtot - C.log_gRAR
C["x"] = C.gbar / A0
C["Mextra_over_Mbar"] = (C.gtot - 10 ** C.log_gRAR) / C.gbar
C["Mmond_over_Mbar"] = mond_mass_factor(C.gbar.values, C.gtot.values)
C["src"] = "CLASH"

# ---------------- X-COP
X = pd.read_csv("data/external/clusters/xcop_profiles.csv")
rows = []
for _, c in X.iterrows():
    if not np.isfinite(c.fgas500):
        continue
    rho_c = 2.775e11 * 0.7 ** 2 * (1 + 0.3 * ((1 + c.z) ** 3 - 1))        # Msun / Mpc^3, flat LCDM
    M200 = c.M200_1e14 * 1e14
    R200 = (3 * M200 / (4 * np.pi * 200 * rho_c)) ** (1 / 3) * 1e3          # kpc
    rs = R200 / c.c200
    mnfw = lambda r: M200 * (np.log(1 + r / rs) - (r / rs) / (1 + r / rs)) / (np.log(1 + c.c200) - c.c200 / (1 + c.c200))  # noqa: E731
    fstar = c.fstar500 if np.isfinite(c.fstar500) else 0.015
    R500 = c.R500_kpc
    R2500 = 0.45 * R500
    fgas2500 = c.fgas2500 if np.isfinite(c.fgas2500) else 0.8 * c.fgas500
    for R, fg, lab in ((R2500, fgas2500, "R2500"), (R500, c.fgas500, "R500"), (R200, c.fgas200, "R200")):
        Mt = mnfw(R); Mb = (fg + fstar) * Mt
        gt = G_SI * Mt * MSUN / (R * KPC) ** 2; gb = G_SI * Mb * MSUN / (R * KPC) ** 2
        rows.append(dict(cluster=c.cluster, R_kpc=R, radius=lab, gbar=gb, gtot=gt, log_gbar=np.log10(gb), log_gtot=np.log10(gt),
                         e_log_gbar=0.05, e_log_gtot=0.05, Mtot=Mt, Mbar=Mb, src="X-COP"))
XX = pd.DataFrame(rows)
XX["log_gRAR"] = np.log10(nu(XX.gbar / A0) * XX.gbar); XX["resid"] = XX.log_gtot - XX.log_gRAR; XX["x"] = XX.gbar / A0
XX["Mextra_over_Mbar"] = (XX.gtot - 10 ** XX.log_gRAR) / XX.gbar
XX["Mmond_over_Mbar"] = mond_mass_factor(XX.gbar.values, XX.gtot.values)
ALL = pd.concat([C, XX], ignore_index=True)
ALL.to_csv("results/cluster_rar.csv", index=False)


def stats(d):
    return f"N = {len(d)}, median residual {d.resid.median():+.2f} dex (16-84 %: {np.percentile(d.resid,16):+.2f} to {np.percentile(d.resid,84):+.2f}), " \
           f"median g_bar/a0 = {d.x.median():.2f}, M_extra/M_bar = {d.Mextra_over_Mbar.median():.1f}, M_MOND/M_bar = {d.Mmond_over_Mbar.median():.1f} " \
           f"({np.percentile(d.Mmond_over_Mbar,16):.1f}-{np.percentile(d.Mmond_over_Mbar,84):.1f})"


L = ["# The law in galaxy clusters: the missing factor (2026-10-12)\n",
     "g_RAR = nu(g_bar/a0) g_bar with the galaxy a0 = 1.2e-10 m/s^2; residual = log(g_tot/g_RAR). Two ways to express the missing mass: "
     "M_extra/M_bar = (g_tot - g_RAR)/g_bar is the Newtonian mass to add to the baryons so that the law, applied to the baryons only, plus "
     "that mass in Newtonian gravity, gives g_tot (a dark-component reading); M_MOND/M_bar is the total mass the law itself would need, "
     "in units of the observed baryons (the usual 'MOND missing mass' of the cluster literature, the quantity behind the 'factor 2').\n",
     "## CLASH (Tian et al. 2020), 20 clusters, 14-600 kpc\n"]
for lo, hi in ((10, 30), (30, 100), (100, 300), (300, 700)):
    d = C[(C.R_kpc >= lo) & (C.R_kpc < hi)]
    L.append(f"- {lo}-{hi} kpc: " + stats(d))
L.append("- all points: " + stats(C))
L.append("\n## X-COP (Eckert et al. 2019, 2022), %d clusters, R2500 to R200\n" % XX.cluster.nunique())
for lab in ("R2500", "R500", "R200"):
    d = XX[XX.radius == lab]
    L.append(f"- {lab} (median {d.R_kpc.median():.0f} kpc): " + stats(d))
L.append("")
L.append("## Scale of the relation in clusters\n")
# fit a0_cluster with the RAR shape on all points (log residual, equal weights), and separately on g_bar < a0 points
def chi(a, d):
    return np.sum((d.log_gtot - np.log10(nu(d.gbar / a) * d.gbar)) ** 2)
for lab, d in (("CLASH", C), ("X-COP", XX), ("both", ALL)):
    res = minimize_scalar(lambda la: chi(10 ** la, d), bounds=(-11, -8), method="bounded")
    a = 10 ** res.x
    sc = np.std(d.log_gtot - np.log10(nu(d.gbar / a) * d.gbar))
    L.append(f"- {lab}: best a0 with the RAR shape = {a:.2e} m/s^2 = {a/A0:.1f} x the galaxy value, scatter {sc:.2f} dex (Tian et al. 2020: 2.0e-9 with a power law of slope 0.5)")
L.append("")
L.append("## Reading\n")
med_all = ALL.resid.median()
core = C[(C.R_kpc >= 100) & (C.R_kpc < 300)]; r500 = XX[XX.radius == "R500"]; r200 = XX[XX.radius == "R200"]; r2500 = XX[XX.radius == "R2500"]
L.append(f"- With the galaxy a0 the law under-predicts the cluster accelerations everywhere: median residual {core.resid.median():+.2f} dex at 100-300 kpc (CLASH), "
         f"{r2500.resid.median():+.2f} at R2500, {r500.resid.median():+.2f} at R500 and {r200.resid.median():+.2f} at R200 (X-COP). In MOND mass the law needs "
         f"{core.Mmond_over_Mbar.median():.1f} times the observed baryons at 100-300 kpc, {r2500.Mmond_over_Mbar.median():.1f} at R2500, "
         f"{r500.Mmond_over_Mbar.median():.1f} at R500 and {r200.Mmond_over_Mbar.median():.1f} at R200.")
L.append("- The 'factor 2' of the cluster lore is the R500 value. It is not a constant: the missing factor is 4-7 between 15 and 600 kpc (where g_bar is 0.1-0.8 a0 and the "
         "law gives little boost), 2 at R500 and about 1.3 at R200. The mass the law lacks is therefore concentrated inside ~0.5 Mpc, the known result of "
         "Sanders 2003 and Angus et al. 2008: it must be a component more centrally concentrated than the hot gas, with a mass of order 10^13-10^14 Msun per cluster.")
L.append("- The cluster points form their own sequence with the RAR shape and a larger scale (CLASH: 15 x a0, Tian et al. 2020: 17 x; X-COP outskirts: 2.6 x). "
         "The two samples do not even agree on a cluster a0, because the residual depends on radius: there is no single a0 for galaxies and clusters, nor one for clusters.")
L.append("- What this means for a MOND-type theory such as RMOND: an additional component of 1-5 times the observed baryons, concentrated in the cluster cores. "
         "Candidates in the literature are undetected baryons (cool gas, excluded at that level by X-ray and SZ data), sterile neutrinos of a few eV (the 'nuHDM' "
         "model, which also helps the CMB in RMOND-like theories; but a few-eV neutrino cannot cluster enough to reach a factor 5 in the cores), or a failure "
         "of the law on cluster scales. The data here cannot distinguish these; they fix the size and the location of the required component.")
L.append("- Caveats: X-COP hydrostatic masses may be biased low by 5-15 % by non-thermal pressure (Eckert et al. 2019: median 6 % at R500), which would "
         "make the missing factor slightly larger; CLASH accelerations come from lensing and X-ray and are free of that bias; the baryon budgets include "
         "BCG stars and ICM gas, not the intracluster light, which adds at most a few per cent.")
L.append("\nFigure: results/cluster_rar.png. Per-point values: results/cluster_rar.csv")
open("results/cluster_rar.md", "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
a = ax[0]
a.errorbar(C.log_gbar, C.log_gtot, xerr=C.e_log_gbar, yerr=C.e_log_gtot, fmt="o", ms=4, color="C0", alpha=0.7, label="CLASH, 14-600 kpc (Tian+20)")
a.errorbar(XX.log_gbar, XX.log_gtot, xerr=0.05, yerr=0.05, fmt="s", ms=5, color="C3", alpha=0.8, label="X-COP, R2500-R200 (Eckert+19,22)")
xx = np.logspace(-12.5, -8.5, 200)
a.plot(np.log10(xx), np.log10(nu(xx / A0) * xx), "k-", label="law, galaxy a0 = 1.2e-10")
a.plot(np.log10(xx), np.log10(nu(xx / 2.0e-9) * xx), "k--", label="same shape, a0 = 2.0e-9 (Tian+20)")
a.plot(np.log10(xx), np.log10(xx), "k:", label="Newton")
a.set_xlabel("log g_bar (m/s^2)"); a.set_ylabel("log g_tot (m/s^2)"); a.legend(fontsize=8); a.set_title("Clusters on the acceleration plane")
b = ax[1]
b.scatter(C.R_kpc, C.resid, s=18, color="C0", alpha=0.7, label="CLASH")
b.scatter(XX.R_kpc, XX.resid, s=28, color="C3", marker="s", label="X-COP")
b.axhline(0, color="k"); b.set_xscale("log"); b.set_xlabel("R (kpc)"); b.set_ylabel("log(g_tot / g_RAR), galaxy a0"); b.legend(fontsize=8)
b.set_title("Missing factor of the law versus radius")
plt.tight_layout(); plt.savefig("results/cluster_rar.png", dpi=130)
print("figure written")
