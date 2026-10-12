"""External-field-effect test on SPARC, after Chae et al. 2020 (ApJ 904, 51).

MOND predicts that a galaxy embedded in an external gravitational field g_ext has
weaker internal MOND boost; its outer rotation curve declines. Newtonian gravity
with dark matter predicts no such dependence (strong equivalence principle).

Two independent quantities per galaxy:
1. e_N,env = |g_ext,N| / a0, the Newtonian field of the baryonic matter around the
   galaxy, computed here from the 2M++ redshift catalogue (Lavaux & Hudson 2011):
   K_s luminosity -> baryonic mass 0.8 L_K (Msun/Lsun,K), weighted for redshift
   incompleteness and for galaxies below the magnitude limit (Schechter alpha=-1,
   M*_K = -24.0), vector sum of G M / r^2 over all catalogue galaxies farther than
   1 Mpc from the SPARC galaxy, the SPARC galaxy placed at its SPARC distance.
2. e_N,fit, the external field that best fits the rotation curve with the 1D
   QUMOND external-field formula (Chae et al. 2020, eq. 6):
       nu_e(z) = nu(z + e_N) + (e_N / z) [nu(z + e_N) - nu(e_N)],  z = g_N / a0
   with nu(y) = 1/(1 - exp(-sqrt y)), a0 = 1.2e-10, marginalising over Y*_disk
   (0.1 dex prior around 0.5), distance and inclination (SPARC errors as priors).

Tests: correlation of e_N,fit with e_N,env; dependence of the low-acceleration
RAR residual and of the outer slope on e_N,env; chi^2 improvement from the EFE.
Writes results/test_efe.md, results/efe_per_galaxy.csv, results/efe.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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

cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


# ------------------------------------------------------------------ 2M++ field
def load_2mpp():
    rows = []
    for line in open(ROOT / "data/catalogs/2mpp_catalog.dat", encoding="latin1"):
        p = [x.strip() for x in line.split("|")]
        if len(p) < 14:
            continue
        try:
            ra, dec, ks, vcmb = float(p[1]), float(p[2]), float(p[5]), float(p[7])
        except ValueError:
            continue
        c115 = float(p[10]) if p[10] else 1.0
        c125 = float(p[11]) if p[11] else 0.0
        rows.append((p[0], ra, dec, ks, vcmb, c115, c125, int(p[12] or 0)))
    d = pd.DataFrame(rows, columns=["name", "ra", "dec", "ks", "vcmb", "c115", "c125", "zoa"])
    d = d[(d.vcmb > 0) & (d.vcmb < 20000)].copy()
    d["dist"] = np.maximum(d.vcmb / H0_KMS, 1.0)
    mk = d.ks - 5 * np.log10(d.dist) - 25
    d["lk"] = 10 ** (-0.4 * (mk - LSUN_K_MAG))
    klim = np.where(d.c125 > 0, 12.5, 11.5)
    comp = np.where(d.ks <= 11.5, d.c115, np.where(d.c125 > 0, d.c125, d.c115))
    w_comp = 1.0 / np.clip(comp, 0.3, 1.0)
    mlim = klim - 5 * np.log10(d.dist) - 25
    llim_over_lstar = 10 ** (-0.4 * (mlim - MSTAR_K))
    w_faint = 1.0 / np.exp(-llim_over_lstar)           # Schechter alpha = -1: fraction above limit = exp(-L_lim/L*)
    d["w"] = np.clip(w_comp * w_faint, 1.0, 5.0)
    d["mass"] = ML_K * d.lk * MSUN * d.w
    ra, dec = np.radians(d.ra), np.radians(d.dec)
    d["x"] = d.dist * np.cos(dec) * np.cos(ra)
    d["y"] = d.dist * np.cos(dec) * np.sin(ra)
    d["z"] = d.dist * np.sin(dec)
    return d


def external_field(d, ra_deg, dec_deg, dist_mpc, rmin=1.0):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    p = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]) * dist_mpc
    dx, dy, dz = d.x.values - p[0], d.y.values - p[1], d.z.values - p[2]
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    sep = np.degrees(np.arccos(np.clip(np.sin(dec) * np.sin(np.radians(d.dec.values)) +
                                       np.cos(dec) * np.cos(np.radians(d.dec.values)) * np.cos(ra - np.radians(d.ra.values)), -1, 1)))
    keep = (r > rmin) & (sep > 0.05)
    acc = G * d.mass.values[keep] / (r[keep] * MPC) ** 2
    ux, uy, uz = dx[keep] / r[keep], dy[keep] / r[keep], dz[keep] / r[keep]
    vec = np.array([np.sum(acc * ux), np.sum(acc * uy), np.sum(acc * uz)])
    local = keep & (r < 10.0)
    acc_l = G * d.mass.values[local] / (r[local] * MPC) ** 2
    vec_l = np.array([np.sum(acc_l * dx[local] / r[local]), np.sum(acc_l * dy[local] / r[local]), np.sum(acc_l * dz[local] / r[local])])
    return np.linalg.norm(vec) / A0, np.linalg.norm(vec_l) / A0, vec


# ------------------------------------------------------------------ EFE rotation-curve model
def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def nu_e(z, e):
    z = np.maximum(z, 1e-12)
    if e <= 0:
        return nu(z)
    return nu(z + e) + (e / z) * (nu(z + e) - nu(e))


EN_GRID = np.concatenate([[0.0], np.logspace(-3, 0, 31)])


def fit_galaxy(df, D, eD, inc, einc):
    r0 = df.R.values * KPC
    gas, disk, bul = df.Vgas.values * 1e3, df.Vdisk.values * 1e3, df.Vbul.values * 1e3
    vobs0, e0 = df.Vobs.values * 1e3, np.maximum(df.eVobs.values, 0.5) * 1e3
    best = dict(chi2=np.inf)
    prof = np.full(len(EN_GRID), np.inf)      # chi2 profile in e_N, minimised over nuisances
    for yd in (0.3, 0.4, 0.5, 0.6, 0.7):
        pen_y = (np.log10(yd / 0.5) / 0.1) ** 2
        for kd in (-2, -1, 0, 1, 2):
            fd = max(1 + kd * eD / D, 0.3)
            pen_d = kd**2
            r = r0 * fd
            vb2 = (gas * np.abs(gas) + yd * disk * np.abs(disk) + 1.4 * yd * bul * np.abs(bul)) * fd
            gN = np.maximum(vb2, 1.0) / r
            z = gN / A0
            for ki in (-2, -1, 0, 1, 2):
                i1 = np.clip(inc + ki * einc, 10.0, 89.0)
                s = np.sin(np.radians(inc)) / np.sin(np.radians(i1))
                vobs, e = vobs0 * s, e0 * s
                pen = pen_y + pen_d + ki**2
                for j, en in enumerate(EN_GRID):
                    v = np.sqrt(gN * nu_e(z, en) * r)
                    chi2 = np.sum(((vobs - v) / e) ** 2) + pen
                    if chi2 < prof[j]:
                        prof[j] = chi2
                    if chi2 < best["chi2"]:
                        best = dict(chi2=chi2, en=en, yd=yd, fd=fd, inc=i1)
    ok = prof <= prof.min() + 1.0
    lo, hi = EN_GRID[ok].min(), EN_GRID[ok].max()
    return best, prof, lo, hi


# ------------------------------------------------------------------ main
print("loading 2M++ ...")
d2 = load_2mpp()
print(f"  {len(d2)} galaxies with distances, total baryonic mass {d2.mass.sum() / MSUN:.2e} Msun within 285 Mpc")
# calibration: field at the Milky Way (origin). The Local Group moves at 630 km/s w.r.t. the CMB,
# i.e. a mean acceleration of ~630 km/s / t_H = 1.4e-12 m/s^2 = 0.012 a0 toward l~270, b~30 (Virgo / Great Attractor)
mw_en, mw_loc, mw_vec = external_field(d2, 0.0, 0.0, 0.0, rmin=1.0)
mw_boost = nu(mw_en) * mw_en
vx, vy, vz = mw_vec / np.linalg.norm(mw_vec)
ra_v, dec_v = np.degrees(np.arctan2(vy, vx)) % 360, np.degrees(np.arcsin(vz))
print(f"  Milky Way calibration: baryonic Newtonian e_N = {mw_en:.2e}, MOND-boosted nu(e)e = {mw_boost:.3f} (expected ~0.012 from the 630 km/s CMB dipole), direction RA={ra_v:.0f} Dec={dec_v:.0f} (Virgo: RA 187, Dec +13; Great Attractor: RA ~200, Dec -45)")
MW_CAL = (mw_en, mw_boost, ra_v, dec_v)

rows = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    name = f.stem.replace("_rotmod", "")
    c = cat.loc[name]
    best, prof, lo, hi = fit_galaxy(df, float(c.D_Mpc), float(c.e_D_Mpc), float(c.Inc_deg), float(c.e_Inc_deg))
    dchi2 = prof[0] - prof.min()
    # RAR residual at low acceleration, no EFE, fixed Y*
    r = df.R.values * KPC
    vb2 = fb.baryonic2(df, 0.5, 1.0, 0.7)
    gN = np.maximum(vb2, 1.0) / r
    gobs = (df.Vobs.values * 1e3) ** 2 / r
    lowg = gN < 0.1 * A0
    res_low = float(np.mean(np.log10(gobs[lowg] / (gN[lowg] * nu(gN[lowg] / A0))))) if lowg.sum() >= 2 else np.nan
    x, y = np.log(df.R.values[-3:]), np.log(df.Vobs.values[-3:])
    slope = float(np.polyfit(x, y, 1)[0])
    en_env = en_loc = np.nan
    if pd.notna(c.RA_deg) and pd.notna(c.Dec_deg):
        en_env, en_loc, _ = external_field(d2, float(c.RA_deg), float(c.Dec_deg), float(c.D_Mpc))
    rows.append(dict(Galaxy=name, T=int(c["T"]), Q=int(c.Q), D=float(c.D_Mpc), en_env=en_env, en_local10Mpc=en_loc,
                     en_fit=best["en"], en_lo=lo, en_hi=hi, dchi2=dchi2, chi2=best["chi2"], N=len(df),
                     yd=best["yd"], res_lowg=res_low, outer_slope=slope))
T = pd.DataFrame(rows)
T["e_env_mond"] = nu(T.en_env.fillna(0).values) * T.en_env.fillna(0).values
T.loc[T.en_env.isna(), "e_env_mond"] = np.nan
T.to_csv(ROOT / "results/efe_per_galaxy.csv", index=False)
from scipy.stats import mannwhitneyu

# ------------------------------------------------------------------ statistics
L = []
L.append(f"{len(T)} galaxies fitted; {T.en_env.notna().sum()} with coordinates and an environmental field from 2M++.")
L.append(f"Calibration at the Milky Way: baryonic Newtonian e_N = {MW_CAL[0]:.2e}, MOND-boosted {MW_CAL[1]:.3f} (expected ~0.012), direction RA {MW_CAL[2]:.0f}, Dec {MW_CAL[3]:.0f}.")
L.append(f"e_N,env (Newtonian, baryons only): median {T.en_env.median():.2e}, 16-84 % {T.en_env.quantile(.16):.2e} to {T.en_env.quantile(.84):.2e} "
         f"({np.log10(T.en_env.quantile(.84) / T.en_env.quantile(.16)):.2f} dex of dynamic range); MOND-boosted nu(e)e: median {T.e_env_mond.median():.3f}; local part < 10 Mpc median {T.en_local10Mpc.median():.2e}.")
L.append(f"e_N,fit: median {T.en_fit.median():.3f}; galaxies preferring e_N > 0 at Delta chi2 > 1: {(T.dchi2 > 1).mean():.0%}, at > 4: {(T.dchi2 > 4).mean():.0%}; "
         f"galaxies with e_N = 0 excluded at 1 sigma (lower limit > 0): {(T.en_lo > 0).mean():.0%}.")
m = T.en_env.notna()
r1, p1 = spearmanr(T.en_fit[m], T.en_env[m])
L.append(f"Spearman(e_N,fit, e_N,env), all {m.sum()}: rho = {r1:+.2f}, p = {p1:.3f}")
m2 = m & (T.dchi2 > 1)
r2, p2 = spearmanr(T.en_fit[m2], T.en_env[m2]) if m2.sum() > 5 else (np.nan, np.nan)
L.append(f"Spearman(e_N,fit, e_N,env), galaxies where the curve constrains e_N (Delta chi2 > 1), N = {m2.sum()}: rho = {r2:+.2f}, p = {p2:.3f}")
m3 = m & T.res_lowg.notna()
r3, p3 = spearmanr(T.res_lowg[m3], T.en_env[m3])
L.append(f"Spearman(low-acceleration RAR residual, e_N,env), N = {m3.sum()}: rho = {r3:+.2f}, p = {p3:.3f}  (MOND predicts negative)")
r4, p4 = spearmanr(T.outer_slope[m], T.en_env[m])
L.append(f"Spearman(outer slope d ln V / d ln R, e_N,env), N = {m.sum()}: rho = {r4:+.2f}, p = {p4:.3f}  (MOND predicts negative)")
q1, q3 = T.en_env.quantile(.25), T.en_env.quantile(.75)
lo_q, hi_q = T[m & (T.en_env <= q1)], T[m & (T.en_env >= q3)]
u, pu = mannwhitneyu(hi_q.en_fit, lo_q.en_fit, alternative="greater")
L.append(f"Chae-style quartile test: median e_N,fit = {hi_q.en_fit.median():.3f} in the top quartile of e_N,env vs {lo_q.en_fit.median():.3f} in the bottom quartile; Mann-Whitney one-sided p = {pu:.3f}")
dec = T[m & (T.outer_slope < -0.1)]
L.append(f"Galaxies with declining outer curves (slope < -0.1): {len(dec)}; their e_N,env rank (0 = weakest field, 1 = strongest), median {T.en_env.rank(pct=True)[dec.index].median():.2f} (MOND predicts > 0.5)")
T["tercile"] = pd.qcut(T.en_env, 3, labels=["low", "mid", "high"])
L.append("\n| e_N,env tercile | N | median e_N,env | median e_N,fit | median Delta chi2 | median low-g residual (dex) | median outer slope |")
L.append("|---|---|---|---|---|---|---|")
for t in ("low", "mid", "high"):
    s = T[T.tercile == t]
    L.append(f"| {t} | {len(s)} | {s.en_env.median():.2e} | {s.en_fit.median():.3f} | {s.dchi2.median():.1f} | {s.res_lowg.median():+.3f} | {s.outer_slope.median():+.2f} |")
L.append("\nGalaxies with the strongest fitted external field (Delta chi2 > 4):")
L.append("| Galaxy | e_N,fit [1 sigma] | e_N,env | Delta chi2 | outer slope |")
L.append("|---|---|---|---|---|")
for _, s in T[T.dchi2 > 4].sort_values("en_fit", ascending=False).head(12).iterrows():
    L.append(f"| {s.Galaxy} | {s.en_fit:.3f} [{s.en_lo:.3f}, {s.en_hi:.3f}] | {s.en_env:.2e} | {s.dchi2:.1f} | {s.outer_slope:+.2f} |")

# ------------------------------------------------------------------ figure
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
s = T[m]
yerr = np.array([np.maximum(s.en_fit - s.en_lo, 1e-4), np.maximum(s.en_hi - s.en_fit, 1e-4)])
ax[0].errorbar(s.en_env, np.maximum(s.en_fit, 5e-4), yerr=yerr, fmt="o", ms=4, alpha=0.6, lw=0.5, color="gray", label="all")
hs = s[s.dchi2 > 4]
ax[0].plot(hs.en_env, np.maximum(hs.en_fit, 5e-4), "o", ms=6, color="tab:red", label="Delta chi2 > 4")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("e_N,env from 2M++"); ax[0].set_ylabel("e_N,fit from the rotation curve")
ax[0].set_title(f"Spearman rho = {r1:+.2f} (p = {p1:.3f})"); ax[0].legend(fontsize=8)
s3 = T[m3]
ax[1].plot(s3.en_env, s3.res_lowg, "o", ms=4, alpha=0.7, color="tab:blue")
ax[1].axhline(0, color="k", lw=0.8); ax[1].set_xscale("log")
ax[1].set_xlabel("e_N,env from 2M++"); ax[1].set_ylabel("mean log(g_obs / g_RAR) at g_N < 0.1 a0")
ax[1].set_title(f"low-acceleration residual, rho = {r3:+.2f} (p = {p3:.3f})")
fig.tight_layout(); fig.savefig(ROOT / "results/efe.png", dpi=130)

with open(ROOT / "results/test_efe.md", "w", encoding="utf-8") as fh:
    fh.write("# External-field-effect test on SPARC (Chae et al. 2020 method), 2026-10-10\n\n" + "\n".join(L) +
             "\n\nFigure: results/efe.png. Per-galaxy values: results/efe_per_galaxy.csv\n")
print("\n".join(L))
print("wrote results/test_efe.md, results/efe.png")
