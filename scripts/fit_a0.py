"""Measure a0 on SPARC and compare with the candidate cosmological expressions.

g_obs = nu(g_bar/a0) g_bar with nu(y) = 1/(1 - exp(-sqrt y)); a0 scanned globally.
Three treatments of the stellar M/L:
  (i)   Y_disk = SED value per galaxy (table enrichie), Y_bulge = 1.4 Y_disk, fixed
  (ii)  Y_disk = 0.5 fixed
  (iii) Y_disk free per galaxy with a 0.1 dex prior around the SED value (profile likelihood)
Statistic: total chi^2 over all points with eVobs >= 0.5 km/s; the 1 sigma interval is
Delta chi^2 = 1 on the profile; a systematic band is the spread between treatments.
Writes results/fit_a0.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, KPC = fb.G, fb.KPC
C = 2.99792458e8
H0 = 70e3 / 3.0857e22
OL = 0.7
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")

gal = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    c = cat.loc[f.stem.replace("_rotmod", "")]
    gal.append((df, float(c.Y36) if pd.notna(c.Y36) else 0.5))


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def chi2_galaxy(df, yd, a0):
    r = df.R.values * KPC
    gN = np.maximum(fb.baryonic2(df, yd, 1.0, 1.4 * yd), 1.0) / r
    v = np.sqrt(gN * nu(gN / a0) * r) / 1e3
    return np.sum(((df.Vobs.values - v) / np.maximum(df.eVobs.values, 0.5)) ** 2)


A0S = np.linspace(0.6e-10, 2.0e-10, 141)
npts = sum(len(df) for df, _ in gal)
res = {}
for label in ("Y = SED per galaxy", "Y = 0.5 fixed", "Y free, 0.1 dex prior on SED"):
    tot = np.zeros_like(A0S)
    for i, a0 in enumerate(A0S):
        s = 0.0
        for df, ysed in gal:
            if label.startswith("Y = SED"):
                s += chi2_galaxy(df, ysed, a0)
            elif label.startswith("Y = 0.5"):
                s += chi2_galaxy(df, 0.5, a0)
            else:
                best = min(chi2_galaxy(df, yd, a0) + (np.log10(yd / ysed) / 0.1) ** 2 for yd in ysed * 10 ** np.linspace(-0.3, 0.3, 13))
                s += best
        tot[i] = s
    j = int(np.argmin(tot))
    inside = A0S[tot <= tot[j] + tot[j] / npts]      # errors rescaled so that chi2/N = 1
    res[label] = (A0S[j], inside.min(), inside.max(), tot[j])

lines = ["| M/L treatment | best a0 (1e-10 m/s^2) | 1 sigma with errors rescaled to chi2/N = 1 | chi2_min / N_points |", "|---|---|---|---|"]
for k, (a, lo, hi, chi) in res.items():
    lines.append(f"| {k} | {a * 1e10:.3f} | {lo * 1e10:.3f} - {hi * 1e10:.3f} | {chi / npts:.2f} |")
vals = [v[0] for v in res.values()]
lines.append(f"\nSystematic band from the M/L treatment: {min(vals) * 1e10:.2f} - {max(vals) * 1e10:.2f} (x1e-10). Lelli et al. 2017: 1.20 +/- 0.02 (random) +/- 0.24 (systematic).")
lines.append("\n| candidate expression | value (1e-10 m/s^2), H0 = 70 | H0 = 73 |")
lines.append("|---|---|---|")
for name, f in (("c H0 / (2 pi)", lambda h: C * h / (2 * np.pi)), ("c H0 / 6  (Verlinde 2016)", lambda h: C * h / 6),
                ("c sqrt(Lambda/3) / (2 pi) = c H0 sqrt(Omega_L) / (2 pi)", lambda h: C * h * np.sqrt(OL) / (2 * np.pi)),
                ("c sqrt(Lambda/3) / 6", lambda h: C * h * np.sqrt(OL) / 6), ("c H0 / 5", lambda h: C * h / 5)):
    lines.append(f"| {name} | {f(H0) * 1e10:.2f} | {f(73e3 / 3.0857e22) * 1e10:.2f} |")
with open(ROOT / "results/fit_a0.md", "w", encoding="utf-8") as fh:
    fh.write("# a0 measured on SPARC (2026-10-10)\n\n" + "\n".join(lines) + "\n")
print("\n".join(lines))
