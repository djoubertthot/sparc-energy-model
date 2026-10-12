"""Measure the interpolation function nu(g_bar/a0) on SPARC and test the universality of a0.

1. Binned acceleration relation (SED M/L, points with dV/V < 0.1): median log g_obs in
   20 bins of log g_bar, with bootstrap errors; compared with candidate functions, each
   with a0 fitted:
     RAR      nu = 1/(1 - exp(-sqrt y))                         McGaugh, Lelli & Schombert 2016
     simple   nu = 1/2 + sqrt(1/4 + 1/y)                        Famaey & Binney 2005
     standard nu = sqrt(1/2 + sqrt(1/4 + 1/y^2))                Milgrom 1983
     n-family nu_n = (1 - exp(-y^(n/2)))^(-1/n), n free         McGaugh 2008 (n = 1 is RAR)
     Verlinde g = g_bar + sqrt(g_bar a0_V / 6), a0_V = c H0     Verlinde 2016 (no free scale)
2. Per-galaxy a0: fitted with Y_disk free (0.1 dex prior); distribution and correlations
   with V_max, surface brightness, type, gas fraction. A universal law predicts a constant.
Writes results/measure_nu.md, results/measure_nu.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KPC = fb.KPC
C = 2.99792458e8
H0 = 70e3 / 3.0857e22
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")

gb, go, gal = [], [], []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    n = f.stem.replace("_rotmod", "")
    c = cat.loc[n]
    yd = float(c.Y36) if pd.notna(c.Y36) else 0.5
    r = df.R.values * KPC
    gN = fb.baryonic2(df, yd, 1.0, 1.4 * yd) / r
    gobs = (df.Vobs.values * 1e3) ** 2 / r
    ok = (gN > 0) & (df.eVobs.values / df.Vobs.values < 0.1)
    gb += list(gN[ok]); go += list(gobs[ok])
    gal.append(dict(name=n, df=df, ysed=yd, T=int(c["T"]), Vmax=df.Vobs.max(), SB=float(c.SBdisk_Lsun_pc2),
                    fgas=float(1.33 * c.sparc_MHI_1e9_Msun / (1.33 * c.sparc_MHI_1e9_Msun + 0.5 * c.sparc_L36_1e9_Lsun))))
gb, go = np.array(gb), np.array(go)
lx = np.log10(gb)
edges = np.linspace(-13.0, -8.5, 21)
rng = np.random.default_rng(0)
bins = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (lx >= lo) & (lx < hi)
    if m.sum() < 8:
        continue
    y = np.log10(go[m])
    boot = [np.median(rng.choice(y, len(y))) for _ in range(300)]
    bins.append(dict(lg=np.median(lx[m]), lo_=np.median(y), err=max(np.std(boot), 0.01), n=int(m.sum())))
B = pd.DataFrame(bins)

NU = {
    "RAR": lambda y: 1 / (1 - np.exp(-np.sqrt(y))),
    "simple": lambda y: 0.5 + np.sqrt(0.25 + 1 / y),
    "standard": lambda y: np.sqrt(0.5 + np.sqrt(0.25 + 1 / y**2)),
}


def chi2_nu(nu, a0):
    g = 10**B.lg
    pred = np.log10(g * nu(g / a0))
    return np.sum(((B.lo_ - pred) / B.err) ** 2)


L = [f"Binned relation: {len(B)} bins, {len(gb)} points, bootstrap errors on bin medians (floor 0.01 dex)."]
L.append("\n| function | best a0 (1e-10) | chi2 / bins | rms of bin residuals (dex) |\n|---|---|---|---|")
for k, nu in NU.items():
    res = minimize_scalar(lambda a: chi2_nu(nu, a), bounds=(0.5e-10, 3e-10), method="bounded")
    g = 10**B.lg
    rr = B.lo_ - np.log10(g * nu(g / res.x))
    L.append(f"| {k} | {res.x * 1e10:.2f} | {res.fun / len(B):.2f} | {rr.std():.3f} |")
best_n, best_v = None, np.inf
for n in np.linspace(0.5, 3.0, 26):
    nu_n = lambda y, n=n: (1 - np.exp(-y ** (n / 2))) ** (-1 / n)
    res = minimize_scalar(lambda a: chi2_nu(nu_n, a), bounds=(0.5e-10, 3e-10), method="bounded")
    if res.fun < best_v:
        best_n, best_v, best_a = n, res.fun, res.x
nu_n = lambda y: (1 - np.exp(-y ** (best_n / 2))) ** (-1 / best_n)
rr = B.lo_ - np.log10(10**B.lg * nu_n(10**B.lg / best_a))
L.append(f"| n-family, n = {best_n:.1f} | {best_a * 1e10:.2f} | {best_v / len(B):.2f} | {rr.std():.3f} |")
g = 10**B.lg
pred_v = np.log10(g + np.sqrt(g * C * H0 / 6))
rr = B.lo_ - pred_v
L.append(f"| Verlinde, a0 = cH0 fixed | (cH0/6 = {C * H0 / 6 * 1e10:.2f}) | {np.sum((rr / B.err) ** 2) / len(B):.2f} | {rr.std():.3f} |")
# deep-regime slope
deep = B[B.lg < np.log10(0.1 * 1.2e-10)]
sl = np.polyfit(deep.lg, deep.lo_, 1)[0]
L.append(f"\nDeep-regime slope d log g_obs / d log g_bar for g_bar < 0.1 a0 ({len(deep)} bins): {sl:.3f} (MOND: 0.500; Newton + baryon-tracing mass: 1.000).")

# per-galaxy a0
def chi2_g(df, yd, a0):
    r = df.R.values * KPC
    gN = np.maximum(fb.baryonic2(df, yd, 1.0, 1.4 * yd), 1.0) / r
    v = np.sqrt(gN * NU["RAR"](gN / a0) * r) / 1e3
    return np.sum(((df.Vobs.values - v) / np.maximum(df.eVobs.values, 0.5)) ** 2)


rows = []
for gl in gal:
    best = None
    for a0 in np.logspace(-11, -9, 41):
        for yd in gl["ysed"] * 10 ** np.linspace(-0.3, 0.3, 7):
            chi = chi2_g(gl["df"], yd, a0) + (np.log10(yd / gl["ysed"]) / 0.1) ** 2
            if best is None or chi < best[0]:
                best = (chi, a0, yd)
    rows.append(dict(Galaxy=gl["name"], a0=best[1], yd=best[2], T=gl["T"], Vmax=gl["Vmax"], SB=gl["SB"], fgas=gl["fgas"], chi2=best[0]))
P = pd.DataFrame(rows)
P.to_csv(ROOT / "results/a0_per_galaxy.csv", index=False)
inner = P[(P.a0 > 1.1e-11) & (P.a0 < 9e-10)]
la = np.log10(inner.a0)
L.append(f"\nPer-galaxy a0 (Y free with prior), {len(inner)} galaxies off the grid edges: median {10**la.median() * 1e10:.2f}, scatter {la.std():.2f} dex, 16-84 %: {10**la.quantile(.16) * 1e10:.2f}-{10**la.quantile(.84) * 1e10:.2f}.")
for col in ("Vmax", "SB", "T", "fgas"):
    r, p = spearmanr(np.log10(inner[col]) if col != "T" else inner[col], la)
    L.append(f"- Spearman(log a0, {col}) = {r:+.2f}, p = {p:.3f}")

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(lx, np.log10(go), ".", ms=2, alpha=0.25, color="gray")
ax[0].errorbar(B.lg, B.lo_, yerr=B.err, fmt="o", color="k", ms=5, label="binned medians")
xx = np.linspace(-13, -8.5, 200)
for k, nu in NU.items():
    res = minimize_scalar(lambda a: chi2_nu(nu, a), bounds=(0.5e-10, 3e-10), method="bounded")
    ax[0].plot(xx, np.log10(10**xx * nu(10**xx / res.x)), lw=1.2, label=f"{k}, a0={res.x * 1e10:.2f}")
ax[0].plot(xx, np.log10(10**xx + np.sqrt(10**xx * C * H0 / 6)), "--", lw=1.2, label="Verlinde, cH0/6")
ax[0].plot(xx, xx, ":", color="k", lw=0.8)
ax[0].set_xlabel("log g_bar (m/s^2)"); ax[0].set_ylabel("log g_obs"); ax[0].legend(fontsize=7); ax[0].set_title("Acceleration relation, SED M/L")
ax[1].hist(la, bins=np.arange(-11, -9, 0.1), color="tab:orange", alpha=0.7)
ax[1].axvline(np.log10(1.2e-10), color="k"); ax[1].set_xlabel("log a0 per galaxy"); ax[1].set_title(f"a0 per galaxy, scatter {la.std():.2f} dex")
fig.tight_layout(); fig.savefig(ROOT / "results/measure_nu.png", dpi=130)
with open(ROOT / "results/measure_nu.md", "w", encoding="utf-8") as fh:
    fh.write("# The interpolation function and the universality of a0 (2026-10-10)\n\n" + "\n".join(L) + "\n\nFigure: results/measure_nu.png\n")
print("\n".join(L))
