"""Scaling relations of the fitted haloes (from results/halo_fits_per_galaxy.csv).

- core radius vs disc scale length: log r_0 = a + b log R_d
- core surface density vs V_max:     log(rho_0 r_0) = a + b log V_max
- core/cusp preference vs central surface brightness and gas fraction
- RAR scatter with SED M/L vs fixed M/L (from the rotmod files)
Writes results/halo_relations.md
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
A0 = 1.2e-10
T = pd.read_csv(ROOT / "results/halo_fits_per_galaxy.csv").set_index("Galaxy")
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")
T["SB"] = cat.SBdisk_Lsun_pc2.reindex(T.index)
T["fgas"] = (1.33 * cat.sparc_MHI_1e9_Msun / (1.33 * cat.sparc_MHI_1e9_Msun + 0.5 * cat.sparc_L36_1e9_Lsun)).reindex(T.index)
T["dchi2"] = T.chi2_bur - T.chi2_nfw

L = []


def linfit(x, y, label):
    m = np.isfinite(x) & np.isfinite(y)
    b, a = np.polyfit(x[m], y[m], 1)
    sc = np.std(y[m] - (a + b * x[m]))
    rho = spearmanr(x[m], y[m])[0]
    L.append(f"- {label}: slope {b:+.2f}, intercept {a:+.2f}, scatter {sc:.2f} dex, Spearman {rho:+.2f}, N = {m.sum()}")


okb = (T.r_0 > 0.15) & (T.r_0 < 150)
linfit(np.log10(T.Rd[okb]), np.log10(T.r_0[okb]), "log r_0 (kpc) vs log R_d (kpc)")
linfit(np.log10(T.Vmax[okb]), np.log10(T.core_sd[okb]), "log rho_0 r_0 (Msun/pc^2) vs log V_max (km/s)")
linfit(np.log10(T.Vmax[okb]), np.log10(T.r_0[okb]), "log r_0 vs log V_max")
linfit(np.log10(T.Vmax[okb]), np.log10(T.rho_0[okb]), "log rho_0 (Msun/pc^3) vs log V_max")
for col, lab in (("SB", "central surface brightness"), ("fgas", "gas fraction"), ("Vmax", "V_max")):
    r, p = spearmanr(T[col], T.dchi2, nan_policy="omit")
    L.append(f"- core preference (chi2_Burkert - chi2_NFW, negative = core) vs {lab}: Spearman {r:+.2f}, p = {p:.3f}")
T["SBbin"] = pd.qcut(T.SB, 3, labels=["low SB", "mid", "high SB"])
L.append("\n| surface-brightness tercile | N | Burkert preferred | median Delta chi2 |\n|---|---|---|---|")
for t in ("low SB", "mid", "high SB"):
    s = T[T.SBbin == t]
    L.append(f"| {t} | {len(s)} | {(s.dchi2 < 0).mean():.0%} | {s.dchi2.median():+.1f} |")

# RAR scatter with SED vs fixed M/L
def rar_scatter(use_sed):
    res = []
    for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
        df = fb.load_rotmod(f)
        if len(df) < 5:
            continue
        c = cat.loc[f.stem.replace("_rotmod", "")]
        yd = float(c.Y36) if (use_sed and pd.notna(c.Y36)) else 0.5
        r = df.R.values * KPC
        gN = fb.baryonic2(df, yd, 1.0, 1.4 * yd if use_sed else 0.7) / r
        gobs = (df.Vobs.values * 1e3) ** 2 / r
        ok = (gN > 0) & (df.eVobs.values / df.Vobs.values < 0.1)
        res += list(np.log10(gobs[ok] / (gN[ok] / (1 - np.exp(-np.sqrt(gN[ok] / A0))))))
    res = np.array(res)
    return len(res), res.std(), res.mean()


for use_sed, lab in ((True, "SED M/L per galaxy"), (False, "fixed 0.5 / 0.7")):
    n, s, m = rar_scatter(use_sed)
    L.append(f"- RAR residual scatter, {lab}: {s:.3f} dex (mean {m:+.3f}), {n} points")

with open(ROOT / "results/halo_relations.md", "w", encoding="utf-8") as fh:
    fh.write("# Halo scaling relations on SPARC (2026-10-10)\n\n" + "\n".join(L) + "\n")
print("\n".join(L))
