"""Diversity of inner rotation-curve slopes (Oman et al. 2015) on SPARC.

For every galaxy with a measured point on both sides of 2 kpc, compares the
observed circular velocity at 2 kpc with two predictions:
- MOND / RAR from the baryons alone: g = g_N nu(g_N/a0), nu(y) = 1/(1-exp(-sqrt y)),
  a0 = 1.2e-10, Y*_disk = 0.5, Y*_bulge = 0.7  (zero free parameters)
- LambdaCDM NFW halo with the cosmological concentration-mass relation of
  Dutton & Maccio 2014 (Planck, z = 0, scatter 0.11 dex), halo V_max matched to
  the observed V_max after removing the baryons in quadrature, baryons added
  back at 2 kpc  (zero free parameters per galaxy, no adiabatic contraction)

Oman et al. found that, at fixed V_max, observed dwarfs show a much wider spread
of V(2 kpc) than LambdaCDM haloes predict. Here the question is which of the two
zero-parameter predictions tracks the observed V(2 kpc) with the smaller scatter.

Writes results/test_diversity.md, results/diversity_per_galaxy.csv, results/diversity.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
A0 = 1.2e-10
H0 = 70e3 / 3.0857e22
H = 0.7
R2 = 2.0 * KPC

cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def f_nfw(x):
    return np.log(1 + x) - x / (1 + x)


def c_dm14(m200):  # Dutton & Maccio 2014, Planck, z = 0
    return 10 ** (0.905 - 0.101 * np.log10(m200 / (1e12 / H * MSUN)))


def nfw_from_vmax(vmax_halo, dex=0.0):
    """M200, c, r_s for an NFW halo whose peak velocity is vmax_halo (m/s)."""
    c = 10.0
    for _ in range(40):
        v200 = vmax_halo / (0.465 * np.sqrt(c / f_nfw(c)))
        m200 = v200**3 / (10 * G * H0)
        c = c_dm14(m200) * 10**dex
    v200 = vmax_halo / (0.465 * np.sqrt(c / f_nfw(c)))
    m200 = v200**3 / (10 * G * H0)
    r200 = v200 / (10 * H0)
    return m200, c, r200 / c


rows = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    r = df.R.values
    if not (r.min() <= 2.0 <= r.max()):
        continue
    name = f.stem.replace("_rotmod", "")
    R = r * KPC
    vobs = df.Vobs.values * 1e3
    vbar2 = fb.baryonic2(df, 0.5, 1.0, 0.7)
    v2 = np.interp(2.0, r, vobs)
    e2 = np.interp(2.0, r, df.eVobs.values * 1e3)
    vb2_2 = np.interp(2.0, r, vbar2)
    gN = max(vb2_2, 1.0) / R2
    v_mond = np.sqrt(gN * nu(gN / A0) * R2)
    imax = int(np.argmax(vobs))
    vmax = vobs[imax]
    vh2 = max(vmax**2 - max(vbar2[imax], 0.0), (10e3) ** 2)
    preds = {}
    for tag, dex in (("nfw", 0.0), ("nfw_lo", -0.11), ("nfw_hi", +0.11)):
        m200, c, rs = nfw_from_vmax(np.sqrt(vh2), dex)
        m_in = m200 * f_nfw(R2 / rs) / f_nfw(c)
        preds[tag] = np.sqrt(G * m_in / R2 + max(vb2_2, 0.0))
    rows.append(dict(Galaxy=name, T=int(cat.loc[name, "T"]), Vmax=vmax / 1e3, V2_obs=v2 / 1e3, eV2=e2 / 1e3,
                     V2_bar=np.sqrt(max(vb2_2, 0)) / 1e3, V2_mond=v_mond / 1e3, V2_nfw=preds["nfw"] / 1e3,
                     V2_nfw_lo=preds["nfw_lo"] / 1e3, V2_nfw_hi=preds["nfw_hi"] / 1e3,
                     d_mond=np.log10(v2 / v_mond), d_nfw=np.log10(v2 / preds["nfw"])))
T = pd.DataFrame(rows)
T.to_csv(ROOT / "results/diversity_per_galaxy.csv", index=False)

bins = [(0, 60), (60, 100), (100, 150), (150, 400)]
lines = []
lines.append(f"N = {len(T)} galaxies with a measured point on both sides of 2 kpc.\n")
lines.append("| V_max bin (km/s) | N | obs V(2)/V_max: p16, median, p84 | scatter of log V_obs/V_MOND | median | scatter of log V_obs/V_NFW | median | galaxies > 0.15 dex below NFW |")
lines.append("|---|---|---|---|---|---|---|---|")
for lo, hi in bins:
    m = (T.Vmax >= lo) & (T.Vmax < hi)
    if m.sum() < 3:
        continue
    q = (T.V2_obs[m] / T.Vmax[m]).quantile([.16, .5, .84]).values
    lines.append(f"| {lo}-{hi} | {m.sum()} | {q[0]:.2f}, {q[1]:.2f}, {q[2]:.2f} | {T.d_mond[m].std():.3f} | {T.d_mond[m].median():+.3f} | "
                 f"{T.d_nfw[m].std():.3f} | {T.d_nfw[m].median():+.3f} | {(T.d_nfw[m] < -0.15).mean():.0%} |")
lines.append(f"| all | {len(T)} | | {T.d_mond.std():.3f} | {T.d_mond.median():+.3f} | {T.d_nfw.std():.3f} | {T.d_nfw.median():+.3f} | {(T.d_nfw < -0.15).mean():.0%} |")
lines.append("\nThe NFW band from the 0.11 dex scatter of the concentration-mass relation spans "
             f"{np.median(np.log10(T.V2_nfw_hi / T.V2_nfw_lo)):.2f} dex in V(2 kpc) (median), which is the LambdaCDM-predicted diversity.")

# figure
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(1, 2, figsize=(12, 5.2))
late = T["T"] >= 8
ax[0].errorbar(T.Vmax, T.V2_obs, yerr=T.eV2, fmt="o", ms=4, color="k", alpha=0.7, lw=0.6, label="observed V(2 kpc)")
o = np.argsort(T.Vmax.values)
ax[0].fill_between(T.Vmax.values[o], T.V2_nfw_lo.values[o], T.V2_nfw_hi.values[o], color="tab:blue", alpha=0.25, label="NFW + baryons, c(M) ± 0.11 dex")
ax[0].plot(T.Vmax, T.V2_mond, "s", ms=4, color="tab:red", alpha=0.8, label="MOND from baryons")
ax[0].plot([20, 350], [20, 350], ":", color="gray")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("V_max (km/s)"); ax[0].set_ylabel("V(2 kpc) (km/s)")
ax[0].legend(fontsize=8); ax[0].set_title("Inner velocity vs. maximum velocity")
ax[1].hist(T.d_nfw, bins=np.arange(-0.6, 0.61, 0.05), alpha=0.5, color="tab:blue", label=f"log V_obs/V_NFW  (σ={T.d_nfw.std():.2f})")
ax[1].hist(T.d_mond, bins=np.arange(-0.6, 0.61, 0.05), alpha=0.5, color="tab:red", label=f"log V_obs/V_MOND (σ={T.d_mond.std():.2f})")
ax[1].axvline(0, color="k", lw=0.8); ax[1].set_xlabel("residual at 2 kpc (dex)"); ax[1].legend(fontsize=8)
ax[1].set_title("Residuals of the two zero-parameter predictions")
fig.tight_layout(); fig.savefig(ROOT / "results/diversity.png", dpi=130)

with open(ROOT / "results/test_diversity.md", "w", encoding="utf-8") as fh:
    fh.write("# Inner-slope diversity test (Oman et al. 2015) on SPARC, 2026-10-10\n\n")
    fh.write("\n".join(lines) + "\n\nFigure: results/diversity.png. Per-galaxy values: results/diversity_per_galaxy.csv\n")
print("\n".join(lines))
print("wrote results/test_diversity.md, results/diversity.png")
