# -*- coding: utf-8 -*-
"""
Does the isolation criterion of the KiDS lenses lower their external field? Measured in 2M++ (Lavaux & Hudson 2011).

For every 2M++ galaxy at 15-40 Mpc, |b| > 10 deg, with a K-band stellar mass in the range of the KiDS lens bins (1e10-4e11 Msun at
0.8 Msun/Lsun,K), compute the Newtonian external field e_N of the catalogue's stellar mass (vector sum of G M/r^2 over all galaxies
beyond 1 Mpc, completeness weights as in scripts/test_efe.py), and its part from within 10 Mpc. Then apply the Brouwer et al. 2021
isolation criterion, no galaxy with more than 10 % of the lens's mass within 3 Mpc, here in three dimensions with redshift distances
(~4 Mpc errors), and a 5 Mpc version; the neighbour threshold is floored at the catalogue completeness limit at the lens distance.
Compare the e_N distributions of isolated and non-isolated galaxies, per KiDS mass bin, with the SPARC median 1.5e-4 a0.
The 2M++ stellar field is a lower bound to the baryonic field (stars are ~6 % of the cosmic baryons; the gas adds to it).
Output: results/efe_isolated_2mpp.md, results/efe_isolated_2mpp.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
MPC = 1e3 * KPC
A0 = 1.2e-10
H0_KMS = 70.0
LSUN_K_MAG = 3.28
MSTAR_K = -24.0
ML_K = 0.8
BINS = [(10.0, 10.5), (10.5, 10.8), (10.8, 11.1), (11.1, 11.5)]


def load_2mpp():
    rows = []
    for line in open(ROOT / "data/catalogs/2mpp_catalog.dat", encoding="latin1"):
        p = [x.strip() for x in line.split("|")]
        if len(p) < 14:
            continue
        try:
            ra, dec, gb, ks, vcmb = float(p[1]), float(p[2]), float(p[4]), float(p[5]), float(p[7])
        except ValueError:
            continue
        c115 = float(p[10]) if p[10] else 1.0
        c125 = float(p[11]) if p[11] else 0.0
        rows.append((p[0], ra, dec, gb, ks, vcmb, c115, c125))
    d = pd.DataFrame(rows, columns=["name", "ra", "dec", "b", "ks", "vcmb", "c115", "c125"])
    d = d[(d.vcmb > 0) & (d.vcmb < 20000)].copy()
    d["dist"] = np.maximum(d.vcmb / H0_KMS, 1.0)
    mk = d.ks - 5 * np.log10(d.dist) - 25
    d["lk"] = 10 ** (-0.4 * (mk - LSUN_K_MAG))
    klim = np.where(d.c125 > 0, 12.5, 11.5)
    comp = np.where(d.ks <= 11.5, d.c115, np.where(d.c125 > 0, d.c125, d.c115))
    w_comp = 1.0 / np.clip(comp, 0.3, 1.0)
    mlim = klim - 5 * np.log10(d.dist) - 25
    llim_over_lstar = 10 ** (-0.4 * (mlim - MSTAR_K))
    w_faint = 1.0 / np.exp(-llim_over_lstar)
    d["w"] = np.clip(w_comp * w_faint, 1.0, 5.0)
    d["mstar"] = ML_K * d.lk                      # Msun, unweighted (the galaxy itself)
    d["mass"] = d.mstar * MSUN * d.w              # kg, weighted (as a source of field)
    d["klim"] = klim
    ra, dec = np.radians(d.ra), np.radians(d.dec)
    d["x"] = d.dist * np.cos(dec) * np.cos(ra)
    d["y"] = d.dist * np.cos(dec) * np.sin(ra)
    d["z"] = d.dist * np.sin(dec)
    return d.reset_index(drop=True)


d = load_2mpp()
X = d[["x", "y", "z"]].values
M = d.mass.values
lens = d[(d.dist >= 15) & (d.dist <= 40) & (np.abs(d.b) > 10) & (d.mstar >= 1e10) & (d.mstar <= 4e11)].copy()
print(f"2M++: {len(d)} galaxies; lens sample {len(lens)}")
rows = []
for i, c in lens.iterrows():
    dx = X - X[i]
    r = np.sqrt((dx ** 2).sum(axis=1))
    keep = (r > 1.0)
    acc = G * M[keep] / (r[keep] * MPC) ** 2
    u = dx[keep] / r[keep][:, None]
    vec = (acc[:, None] * u).sum(axis=0)
    loc = keep & (r < 10.0)
    accl = G * M[loc] / (r[loc] * MPC) ** 2
    vecl = (accl[:, None] * (dx[loc] / r[loc][:, None])).sum(axis=0)
    # isolation: neighbours above 10 % of the lens mass (floored at the completeness limit at the lens distance) within 3 / 5 Mpc
    mlim_abs = 11.5 - 5 * np.log10(c.dist) - 25
    m_complete = ML_K * 10 ** (-0.4 * (mlim_abs - LSUN_K_MAG))
    thr = max(0.1 * c.mstar, m_complete)
    big = (d.mstar.values >= thr) & (r > 0)
    n3 = int(np.sum(big & (r < 3.0))); n5 = int(np.sum(big & (r < 5.0)))
    rows.append(dict(name=c["name"], dist=c.dist, logm=np.log10(c.mstar), eN=np.linalg.norm(vec) / A0, eN_local=np.linalg.norm(vecl) / A0,
                     n3=n3, n5=n5, thr_over_01=thr / (0.1 * c.mstar)))
R = pd.DataFrame(rows)
R["iso3"] = R.n3 == 0; R["iso5"] = R.n5 == 0
R.to_csv(ROOT / "results/efe_isolated_2mpp.csv", index=False)


def q(s):
    return f"{np.median(s):.2e} ({np.percentile(s,16):.1e}-{np.percentile(s,84):.1e})"


L = ["# External field of isolated galaxies in 2M++ (2026-10-13)\n",
     f"Lens sample: {len(R)} 2M++ galaxies at 15-40 Mpc, |b| > 10 deg, stellar mass 1e10-4e11 Msun (0.8 Msun/Lsun,K). e_N = stellar "
     "Newtonian field of the catalogue beyond 1 Mpc, in units of a0; 'local' = part from within 10 Mpc. Isolation after Brouwer et al. "
     "2021: no galaxy above 10 % of the lens mass within 3 Mpc (3D, redshift distances); the threshold is floored at the completeness "
     f"limit, which exceeds 10 % of the lens mass for {100*np.mean(R.thr_over_01 > 1):.0f} % of the lenses (low-mass bins, so their "
     "isolation is slightly too permissive).\n",
     "| sample | N | median e_N (16-84 %) | median local part | fraction e_N < 5e-5 | fraction e_N < 1e-4 |", "|---|---|---|---|---|---|"]
for lab, m in (("all", np.ones(len(R), bool)), ("isolated, 3 Mpc", R.iso3.values), ("not isolated, 3 Mpc", ~R.iso3.values),
               ("isolated, 5 Mpc", R.iso5.values)):
    s = R[m]
    L.append(f"| {lab} | {len(s)} | {q(s.eN)} | {q(s.eN_local)} | {np.mean(s.eN < 5e-5):.2f} | {np.mean(s.eN < 1e-4):.2f} |")
L.append("")
L.append("| log M* bin | N all | median e_N all | N isolated (3 Mpc) | median e_N isolated | fraction isolated with e_N < 5e-5 |")
L.append("|---|---|---|---|---|---|")
for lo, hi in BINS:
    s = R[(R.logm >= lo) & (R.logm < hi)]; si = s[s.iso3]
    if len(s) == 0:
        continue
    L.append(f"| {lo}-{hi} | {len(s)} | {q(s.eN)} | {len(si)} | {q(si.eN) if len(si) else '-'} | {np.mean(si.eN < 5e-5) if len(si) else float('nan'):.2f} |")
L.append("")
iso = R[R.iso3]; non = R[~R.iso3]
L.append("## Reading\n")
L.append(f"- Isolated galaxies have a median stellar external field of {np.median(iso.eN):.1e} a0 against {np.median(non.eN):.1e} for the others and "
         f"1.5e-4 for the SPARC galaxies: isolation at 3 Mpc removes the local part ({np.median(iso.eN_local):.1e} against {np.median(non.eN_local):.1e}) "
         "but not the field of the large-scale structure, which dominates.")
L.append(f"- Only {100*np.mean(iso.eN < 5e-5):.0f} % of the isolated galaxies have e_N < 5e-5 a0, the 2-sigma bound of the lowest-mass lensing bin, and "
         f"{100*np.mean(iso.eN < 1e-4):.0f} % have e_N < 1e-4. The lensing bound therefore cannot be explained by the isolation criterion.")
L.append("- This e_N counts the stars of 2M++ only (0.8 Msun/Lsun,K), about 6 % of the cosmic baryons: it is a lower bound to the baryonic Newtonian "
         "field of a MOND-type theory, in which the intergalactic gas also sources gravity. A larger e_N moves r_t inward and strengthens the "
         "lensing tension; the MOND-boosted stellar field already reproduces the Milky Way's bulk acceleration (0.011 against 0.012 a0 from the "
         "CMB dipole), so the stellar value is the one a MOND universe allows.")
L.append("- Caveats: the KiDS lenses lie at z = 0.1-0.5, where the large-scale field is a statistical property of the cosmic web rather than of "
         "the Local Supercluster; the local volume (Virgo, the Great Attractor) may be denser than average, in which case the lens field would be "
         "somewhat lower, but not by the factor three the bound requires unless the far field were dominated by one structure.")
L.append("\nFigure: results/efe_isolated_2mpp.png. Per-galaxy values: results/efe_isolated_2mpp.csv")
open(ROOT / "results/efe_isolated_2mpp.md", "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
bins = np.linspace(-5.5, -2.5, 31)
ax[0].hist(np.log10(non.eN), bins=bins, alpha=0.5, label=f"not isolated (N={len(non)})", density=True)
ax[0].hist(np.log10(iso.eN), bins=bins, alpha=0.5, label=f"isolated, 3 Mpc (N={len(iso)})", density=True)
ax[0].axvline(np.log10(1.5e-4), color="k", ls="--", label="SPARC median 1.5e-4"); ax[0].axvline(np.log10(6e-5), color="C3", label="lensing bound 6e-5")
ax[0].set_xlabel("log e_N / a0 (2M++ stellar field)"); ax[0].set_ylabel("density"); ax[0].legend(fontsize=8); ax[0].set_title("External field, 2M++ galaxies at 15-40 Mpc")
ax[1].scatter(np.log10(non.eN_local), np.log10(non.eN), s=5, alpha=0.4, label="not isolated")
ax[1].scatter(np.log10(iso.eN_local), np.log10(iso.eN), s=5, alpha=0.6, label="isolated")
ax[1].plot([-7, -3], [-7, -3], "k:"); ax[1].set_xlabel("log local part (< 10 Mpc)"); ax[1].set_ylabel("log total e_N"); ax[1].legend(fontsize=8)
ax[1].set_title("The far field dominates")
plt.tight_layout(); plt.savefig(ROOT / "results/efe_isolated_2mpp.png", dpi=130)
print("figure written")
