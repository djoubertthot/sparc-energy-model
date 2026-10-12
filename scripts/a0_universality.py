"""Is a0 universal? Per-galaxy a0 with the nuisance parameters marginalised.

For each galaxy, the acceleration law g = nu(g_bar/a0) g_bar is fitted with a0 free on a
log grid while Y_disk (0.1 dex prior on the SED value), distance (SPARC error) and
inclination (SPARC error) are varied on grids with Gaussian priors; the per-galaxy a0 is
the profile minimum, with a 1 sigma interval from Delta chi2 = 1 after rescaling the
errors so that chi2/N = 1 at the minimum. A universal law predicts a scatter consistent
with the intervals; correlations of a0 with galaxy properties would falsify it.
Writes results/a0_universality.md, results/a0_universality.csv
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
A0S = np.logspace(-11, -9, 61)


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


rows = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    n = f.stem.replace("_rotmod", "")
    c = cat.loc[n]
    ysed = float(c.Y36) if pd.notna(c.Y36) else 0.5
    D, eD, inc, einc = float(c.D_Mpc), float(c.e_D_Mpc), float(c.Inc_deg), float(c.e_Inc_deg)
    r0 = df.R.values * KPC
    gas, disk, bul = df.Vgas.values * 1e3, df.Vdisk.values * 1e3, df.Vbul.values * 1e3
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
    scale = max(prof[j] / len(df), 1.0)
    ok = prof <= prof[j] + scale
    rows.append(dict(Galaxy=n, T=int(c["T"]), Vmax=df.Vobs.max(), SB=float(c.SBdisk_Lsun_pc2), Rd=float(c.Rdisk_kpc), Q=int(c.Q),
                     fgas=float(1.33 * c.sparc_MHI_1e9_Msun / (1.33 * c.sparc_MHI_1e9_Msun + 0.5 * c.sparc_L36_1e9_Lsun)),
                     a0=A0S[j], a0_lo=A0S[ok].min(), a0_hi=A0S[ok].max(), chi2n=prof[j] / len(df), N=len(df),
                     edge=(j == 0) or (j == len(A0S) - 1)))
T = pd.DataFrame(rows)
T.to_csv(ROOT / "results/a0_universality.csv", index=False)
g = T[~T.edge]
la = np.log10(g.a0)
err = 0.5 * (np.log10(g.a0_hi) - np.log10(g.a0_lo))
w = 1 / np.maximum(err, 0.02) ** 2
mean_w = np.sum(w * la) / np.sum(w)
chi2_const = np.sum(w * (la - mean_w) ** 2)
intrinsic = max(la.var() - np.mean(err**2), 0) ** 0.5
L = [f"{len(T)} galaxies, {len(g)} with a0 inside the grid (1e-11 to 1e-9)."]
L.append(f"Per-galaxy a0: median {10**la.median() * 1e10:.2f} x 1e-10, 16-84 %: {10**la.quantile(.16) * 1e10:.2f}-{10**la.quantile(.84) * 1e10:.2f}, observed scatter {la.std():.2f} dex, median 1 sigma interval {err.median():.2f} dex.")
L.append(f"Intrinsic scatter after subtracting the per-galaxy intervals in quadrature: {intrinsic:.2f} dex.")
L.append(f"Weighted mean log a0: {10**mean_w * 1e10:.2f} x 1e-10; chi2 of a universal a0 = {chi2_const:.0f} for {len(g) - 1} degrees of freedom.")
L.append(f"Galaxies whose 1 sigma interval excludes 1.2e-10: {((g.a0_hi < 1.2e-10) | (g.a0_lo > 1.2e-10)).mean():.0%}; excludes it by a factor 2: {((g.a0_hi < 0.6e-10) | (g.a0_lo > 2.4e-10)).mean():.0%}.")
L.append("\nCorrelations of log a0 (Spearman, p):")
for col in ("Vmax", "SB", "Rd", "fgas", "T", "Q"):
    x = np.log10(g[col]) if col in ("Vmax", "SB", "Rd") else g[col]
    r, p = spearmanr(x, la)
    L.append(f"- {col}: {r:+.2f} (p = {p:.3f})")
hq = g[g.Q == 1]
L.append(f"\nQuality 1 only ({len(hq)}): median {10**np.log10(hq.a0).median() * 1e10:.2f}, scatter {np.log10(hq.a0).std():.2f} dex.")
with open(ROOT / "results/a0_universality.md", "w", encoding="utf-8") as fh:
    fh.write("# Universality of a0 on SPARC (2026-10-10)\n\n" + "\n".join(L) + "\n")
print("\n".join(L))
