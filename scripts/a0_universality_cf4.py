"""Universality of a0 with Cosmicflows-4 rotation-independent distances.

Same per-galaxy fit as a0_universality.py (a0 free, Y_disk / distance / inclination
marginalised on Gaussian priors), run twice: with the SPARC distances and with the
Cosmicflows-4 TRGB / Cepheid / maser / SBF / SN distances where they are tighter
(data/sparc_sample/cf4_distances.csv). Curves are rescaled with distance: R ~ D,
V_gas, V_disk, V_bul ~ sqrt(D). The comparison on the galaxies whose distance changed
tells whether the intrinsic scatter of a0 is a distance effect.
Writes results/a0_universality_cf4.md, results/a0_universality_cf4.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KPC = fb.KPC
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")
cf4 = pd.read_csv(ROOT / "data/sparc_sample/cf4_distances.csv").set_index("Galaxy")
A0S = np.logspace(-11, -9, 61)


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def fit_a0(df, ysed, D, eD, inc, einc, scale):
    r0 = df.R.values * KPC * scale
    sq = np.sqrt(scale)
    gas, disk, bul = df.Vgas.values * 1e3 * sq, df.Vdisk.values * 1e3 * sq, df.Vbul.values * 1e3 * sq
    vobs0, e0 = df.Vobs.values * 1e3, np.maximum(df.eVobs.values, 0.5) * 1e3
    prof = np.full(len(A0S), np.inf)
    for yd in ysed * 10 ** np.array([-0.2, -0.1, 0.0, 0.1, 0.2]):
        pen_y = (np.log10(yd / ysed) / 0.1) ** 2
        for kd in (-2, -1, 0, 1, 2):
            fd = max(1 + kd * eD / D, 0.3)
            r = r0 * fd
            vb2 = (gas * np.abs(gas) + yd * disk * np.abs(disk) + 1.4 * yd * bul * np.abs(bul)) * fd
            gN = np.maximum(vb2, 1.0) / r
            for ki in (-2, -1, 0, 1, 2):
                i1 = np.clip(inc + ki * einc, 10.0, 89.0)
                s = np.sin(np.radians(inc)) / np.sin(np.radians(i1))
                vobs, e = vobs0 * s, e0 * s
                pen = pen_y + kd**2 + ki**2
                for j, a0 in enumerate(A0S):
                    v = np.sqrt(gN * nu(gN / a0) * r)
                    chi2 = np.sum(((vobs - v) / e) ** 2) + pen
                    if chi2 < prof[j]:
                        prof[j] = chi2
    j = int(np.argmin(prof))
    sc = max(prof[j] / len(df), 1.0)
    ok = prof <= prof[j] + sc
    return A0S[j], A0S[ok].min(), A0S[ok].max(), (j == 0) or (j == len(A0S) - 1)


rows = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    n = f.stem.replace("_rotmod", "")
    c = cat.loc[n]
    ysed = float(c.Y36) if pd.notna(c.Y36) else 0.5
    D, eD, inc, einc = float(c.D_Mpc), float(c.e_D_Mpc), float(c.Inc_deg), float(c.e_Inc_deg)
    a_s, lo_s, hi_s, edge_s = fit_a0(df, ysed, D, eD, inc, einc, 1.0)
    replaced = False
    a_c, lo_c, hi_c, edge_c = a_s, lo_s, hi_s, edge_s
    if n in cf4.index and pd.notna(cf4.loc[n, "D_cf4_indep"]):
        Dc = float(cf4.loc[n, "D_cf4_indep"])
        eDc = Dc * (10 ** float(cf4.loc[n, "eD_cf4_indep_dex"]) - 1)
        if eDc / Dc < eD / D:
            replaced = True
            a_c, lo_c, hi_c, edge_c = fit_a0(df, ysed, Dc, eDc, inc, einc, Dc / D)
    rows.append(dict(Galaxy=n, Q=int(c.Q), Vmax=df.Vobs.max(), replaced=replaced, D_sparc=D, D_used=Dc if replaced else D,
                     a0_sparc=a_s, lo_s=lo_s, hi_s=hi_s, edge_s=edge_s, a0_cf4=a_c, lo_c=lo_c, hi_c=hi_c, edge_c=edge_c,
                     methods=cf4.loc[n, "indep_methods"] if replaced else ""))
T = pd.DataFrame(rows)
T.to_csv(ROOT / "results/a0_universality_cf4.csv", index=False)


def stats(a, lo, hi, edge):
    g = ~edge
    la = np.log10(a[g])
    err = 0.5 * (np.log10(hi[g]) - np.log10(lo[g]))
    intr = max(la.var() - np.mean(err**2), 0) ** 0.5
    w = 1 / np.maximum(err, 0.02) ** 2
    mean_w = np.sum(w * la) / np.sum(w)
    chi2 = np.sum(w * (la - mean_w) ** 2)
    return dict(N=int(g.sum()), median=10**np.median(la) * 1e10, wmean=10**mean_w * 1e10, scatter=la.std(), intr=intr, chi2=chi2,
                excl2=float(((hi[g] < 0.6e-10) | (lo[g] > 2.4e-10)).mean()))


L = [f"{len(T)} galaxies; {T.replaced.sum()} distances replaced by Cosmicflows-4 rotation-independent values (TRGB, Cepheids, SBF, SNe)."]
L.append("\n| sample | distances | N | median a0 | weighted mean | observed scatter | intrinsic scatter | chi2 of a universal a0 | exclude 1.2e-10 by x2 |")
L.append("|---|---|---|---|---|---|---|---|---|")
for lab, m in (("all", np.ones(len(T), bool)), ("replaced subset", T.replaced.values), ("quality 1, replaced", (T.replaced & (T.Q == 1)).values)):
    for dl, (a, lo, hi, ed) in (("SPARC", (T.a0_sparc, T.lo_s, T.hi_s, T.edge_s)), ("CF4", (T.a0_cf4, T.lo_c, T.hi_c, T.edge_c))):
        s = stats(a[m].values, lo[m].values, hi[m].values, ed[m].values)
        L.append(f"| {lab} | {dl} | {s['N']} | {s['median']:.2f} | {s['wmean']:.2f} | {s['scatter']:.2f} | {s['intr']:.2f} | {s['chi2']:.0f} | {s['excl2']:.0%} |")
rep = T[T.replaced & ~T.edge_s & ~T.edge_c]
d = np.log10(rep.a0_cf4 / rep.a0_sparc)
L.append(f"\nOn the replaced subset, change of log a0 per galaxy: median {d.median():+.3f}, rms {d.std():.3f} dex; "
         f"correlation with the distance change log(D_cf4/D_sparc): Spearman {spearmanr(np.log10(rep.D_used / rep.D_sparc), d)[0]:+.2f}.")
L.append("\nLargest moves:")
L.append("| Galaxy | methods | D SPARC | D CF4 | a0 SPARC (1e-10) | a0 CF4 |\n|---|---|---|---|---|---|")
for _, s in rep.assign(dd=np.abs(np.log10(rep.D_used / rep.D_sparc))).sort_values("dd", ascending=False).head(10).iterrows():
    L.append(f"| {s.Galaxy} | {s.methods} | {s.D_sparc:.2f} | {s.D_used:.2f} | {s.a0_sparc * 1e10:.2f} | {s.a0_cf4 * 1e10:.2f} |")
with open(ROOT / "results/a0_universality_cf4.md", "w", encoding="utf-8") as fh:
    fh.write("# Universality of a0 with Cosmicflows-4 distances (2026-10-10)\n\n" + "\n".join(L) + "\n")
print("\n".join(L))
