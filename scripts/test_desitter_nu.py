# -*- coding: utf-8 -*-
"""
The one 'quantum' derivation of the law: Milgrom (1999, Phys. Lett. A 253, 273) reads modified inertia as a vacuum effect in
de Sitter space, where an accelerated observer sees the Unruh temperature T(a) ~ sqrt(a^2 + a_L^2) with a_L = c sqrt(Lambda/3).
Inertia proportional to the excess T(a) - T(0) gives F = m a mu(a) with mu(a) = [sqrt(a^2 + a_L^2) - a_L] / a, i.e. deep-MOND
with a0 = 2 a_L and, in closed form, nu(y) = sqrt(1 + 1/y) with y = g_bar / a0. Two things are testable on SPARC:
  (1) the SHAPE of this nu against the measured relation (same binned relation as scripts/measure_nu.py), compared with the
      RAR, simple, standard and n-family functions and with Verlinde (2016);
  (2) the SCALE: 2 c sqrt(Lambda/3) = 2 c H0 sqrt(Omega_L) versus the measured a0, i.e. whether the prefactor is right.
Output: results/test_desitter_nu.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KPC = fb.KPC
C = 2.99792458e8
H0 = 70e3 / 3.0857e22
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")

gb, go = [], []
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
    "RAR (McGaugh 2016)": lambda y: 1 / (1 - np.exp(-np.sqrt(y))),
    "simple": lambda y: 0.5 + np.sqrt(0.25 + 1 / y),
    "standard": lambda y: np.sqrt(0.5 + np.sqrt(0.25 + 1 / y ** 2)),
    "n-family, n = 1.4": lambda y: (1 - np.exp(-y ** 0.7)) ** (-1 / 1.4),
    "de Sitter vacuum (Milgrom 1999): sqrt(1 + 1/y)": lambda y: np.sqrt(1 + 1 / y),
    "Verlinde 2016: 1/2 + sqrt(1/4 + 1/y)... same as simple with a0 = cH0/6, fixed": None,
}


def chi2_nu(nu, a0):
    g = 10 ** B.lg
    pred = np.log10(g * nu(g / a0))
    return np.sum(((B.lo_ - pred) / B.err) ** 2)


def pointwise_rms(nu, a0, lo=-13, hi=-8.5):
    m = (lx > lo) & (lx < hi)
    res = np.log10(go[m]) - np.log10(gb[m] * nu(gb[m] / a0))
    return np.std(res), np.median(res)


L = ["# The de Sitter-vacuum interpolation function on SPARC (2026-10-12)\n",
     f"Binned relation of scripts/measure_nu.py: {len(B)} bins, {len(gb)} points (Y36 per galaxy, errors < 10 %). The binned a0 is "
     "biased high by binning (see results/measure_nu.md); only the SHAPE comparison between functions is meaningful here. The "
     "point-wise columns use a0 = 1.2e-10 on all points.\n",
     "| function | best a0 binned (1e-10) | chi2 / bins | rms of bin residuals (dex) | point-wise scatter at a0 = 1.2e-10 (dex) | point-wise median residual (dex) |",
     "|---|---|---|---|---|---|"]
for name, nu in NU.items():
    if nu is None:
        continue
    res = minimize_scalar(lambda a: chi2_nu(nu, a), bounds=(0.5e-10, 3e-10), method="bounded")
    rr = B.lo_ - np.log10(10 ** B.lg * nu(10 ** B.lg / res.x))
    sc, md = pointwise_rms(nu, 1.2e-10)
    L.append(f"| {name} | {res.x * 1e10:.2f} | {res.fun / len(B):.2f} | {rr.std():.3f} | {sc:.3f} | {md:+.3f} |")
# high-acceleration tail: where the functions differ most (y > 3)
L.append("")
L.append("High-acceleration tail (g_bar > 3 a0), where the functions differ: nu - 1 predicted at g_bar = 3, 10 and 30 a0, and the "
         "measured median residual log(g_obs/g_bar) of the SPARC points in those ranges:")
for name, nu in NU.items():
    if nu is None:
        continue
    L.append(f"- {name}: nu - 1 = {nu(3.0)-1:.3f}, {nu(10.0)-1:.3f}, {nu(30.0)-1:.3f} -> log nu = {np.log10(nu(3.0)):.3f}, {np.log10(nu(10.0)):.3f}, {np.log10(nu(30.0)):.3f}")
for lo, hi in ((0.3, 0.5), (0.5, 1.0), (1.0, 1.5)):
    m = (lx > np.log10(1.2e-10) + lo) & (lx < np.log10(1.2e-10) + hi)
    L.append(f"- measured, g_bar/a0 in [{10**lo:.0f}, {10**hi:.0f}]: median log(g_obs/g_bar) = {np.median(np.log10(go[m]/gb[m])):+.3f} dex (N = {m.sum()})")
# the scale
aL = C * H0 * np.sqrt(0.7)
L.append("")
L.append("## The scale\n")
L.append(f"- a_Lambda = c sqrt(Lambda/3) = c H0 sqrt(Omega_Lambda) = {aL:.2e} m/s^2 (H0 = 70, Omega_Lambda = 0.7).")
L.append(f"- Milgrom 1999 vacuum inertia gives a0 = 2 a_Lambda = {2*aL:.2e}, i.e. {2*aL/1.2e-10:.1f} times the measured 1.2e-10: the derivation "
         "fixes the shape sqrt(1 + 1/y) and the order of magnitude, not the prefactor; a0 = 0.21 a_Lambda is measured, 2 a_Lambda predicted.")
L.append(f"- With the prefactor left free, the Unruh temperature at a0 is T_U(a0)/T_dS = a0/a_Lambda = {1.2e-10/aL:.2f}: the MOND regime starts "
         "where the Unruh temperature of the orbit falls to a fifth of the de Sitter temperature.")
L.append("")
L.append("## Reading\n")
L.append("- The vacuum-inertia function has no free parameter besides a0 and is a closed form; its fit quality relative to the RAR "
         "and n = 1.4 functions is given in the table. It is the only candidate derivation of the law that touches quantum physics "
         "(Unruh radiation in de Sitter space), and it is a heuristic: no field theory, no cosmology, no relativistic limit.")
open(ROOT / "results/test_desitter_nu.md", "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))
