"""Dark-halo fits on SPARC: NFW vs Burkert vs the acceleration law, at equal footing.

Per galaxy, fit V_obs(r) with  V^2 = V_bar^2(Y) + V_halo^2  where
  NFW     : rho = rho_s / [(r/r_s)(1 + r/r_s)^2],  M(r) = 4 pi rho_s r_s^3 [ln(1+x) - x/(1+x)]
  Burkert : rho = rho_0 r_0^3 / [(r + r_0)(r^2 + r_0^2)],
            M(r) = 2 pi rho_0 r_0^3 [ln(1+x) - arctan(x) + ln(1+x^2)/2]
and the stellar M/L Y_disk as a nuisance with a 0.1 dex log-normal prior centred on the
per-galaxy SED value Y36 (table enrichie) or 0.5, Y_bulge = 1.4 Y_disk.
The acceleration law g = nu(g_bar/a0) g_bar is fitted with the same Y prior (1 parameter).

Reported: chi^2, RMS, AIC = chi^2 + 2k (k = 3 for haloes, 1 for the law), fraction of
galaxies preferring each, NFW concentration-mass relation vs Dutton & Maccio 2014,
Burkert core surface density rho_0 r_0 (Donato et al. 2009 find ~140 Msun/pc^2),
core radius vs disc scale length.
Writes results/halo_fits.md, results/halo_fits_per_galaxy.csv, results/halo_fits.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
PC = KPC / 1e3
A0 = 1.2e-10
H0 = 70e3 / 3.0857e22
RHO_CRIT = 3 * H0**2 / (8 * np.pi * G)
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")


def m_nfw(r, rho_s, r_s):
    x = r / r_s
    return 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x / (1 + x))


def m_burkert(r, rho_0, r_0):
    x = r / r_0
    return 2 * np.pi * rho_0 * r_0**3 * (np.log(1 + x) - np.arctan(x) + 0.5 * np.log(1 + x**2))


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def vbar2(df, yd):
    return fb.baryonic2(df, yd, 1.0, 1.4 * yd)


def fit_halo(df, y0, mfun):
    r = df.R.values * KPC
    vobs, e = df.Vobs.values * 1e3, np.maximum(df.eVobs.values, 0.5) * 1e3

    def resid(p):
        lrho, lr, ly = p
        v2 = vbar2(df, 10**ly) + G * mfun(r, 10**lrho * MSUN / PC**3, 10**lr * KPC) / r
        v = np.sqrt(np.maximum(v2, 0))
        return np.concatenate([(vobs - v) / e, [np.log10(10**ly / y0) / 0.1]])

    best = None
    for lrho0 in (-2.5, -1.5, -0.5):
        for lr0 in (0.0, 0.7, 1.4):
            try:
                res = least_squares(resid, [lrho0, lr0, np.log10(y0)], bounds=([-5, -1, np.log10(0.1)], [1.5, 2.7, np.log10(1.5)]))
            except Exception:
                continue
            if best is None or res.cost < best.cost:
                best = res
    chi2 = np.sum(best.fun[:-1] ** 2)
    lrho, lr, ly = best.x
    v2 = vbar2(df, 10**ly) + G * mfun(r, 10**lrho * MSUN / PC**3, 10**lr * KPC) / r
    rms = np.sqrt(np.mean((vobs - np.sqrt(np.maximum(v2, 0))) ** 2)) / 1e3
    return chi2, rms, 10**lrho, 10**lr, 10**ly


def fit_law(df, y0):
    r = df.R.values * KPC
    vobs, e = df.Vobs.values * 1e3, np.maximum(df.eVobs.values, 0.5) * 1e3
    best = None
    for yd in np.linspace(0.15, 1.2, 43):
        gN = np.maximum(vbar2(df, yd), 1.0) / r
        v = np.sqrt(gN * nu(gN / A0) * r)
        chi2 = np.sum(((vobs - v) / e) ** 2)
        tot = chi2 + (np.log10(yd / y0) / 0.1) ** 2
        if best is None or tot < best[0]:
            best = (tot, chi2, np.sqrt(np.mean((vobs - v) ** 2)) / 1e3, yd)
    return best[1], best[2], best[3]


def nfw_c_m(rho_s, r_s):
    """c200 and M200 from (rho_s [Msun/pc^3], r_s [kpc])."""
    rs = r_s * KPC
    rhos = rho_s * MSUN / PC**3
    c = np.logspace(-0.5, 2.5, 3000)
    mean = 3 * rhos * (np.log(1 + c) - c / (1 + c)) / c**3
    i = int(np.argmin(np.abs(mean - 200 * RHO_CRIT)))
    c200 = c[i]
    m200 = 4 / 3 * np.pi * (c200 * rs) ** 3 * 200 * RHO_CRIT
    return c200, m200 / MSUN


rows = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    n = f.stem.replace("_rotmod", "")
    c = cat.loc[n]
    y0 = float(c.Y36) if pd.notna(c.Y36) else 0.5
    cn, rn, rhos, rs, yn = fit_halo(df, y0, m_nfw)
    cb, rb, rho0, r0, yb = fit_halo(df, y0, m_burkert)
    cl, rl, yl = fit_law(df, y0)
    c200, m200 = nfw_c_m(rhos, rs)
    rows.append(dict(Galaxy=n, T=int(c["T"]), N=len(df), Rd=float(c.Rdisk_kpc), Vmax=df.Vobs.max(), y0=y0,
                     chi2_nfw=cn, rms_nfw=rn, rho_s=rhos, r_s=rs, y_nfw=yn, c200=c200, M200=m200,
                     chi2_bur=cb, rms_bur=rb, rho_0=rho0, r_0=r0, y_bur=yb, core_sd=rho0 * r0 * 1e3,
                     chi2_law=cl, rms_law=rl, y_law=yl))
T = pd.DataFrame(rows)
T["aic_nfw"], T["aic_bur"], T["aic_law"] = T.chi2_nfw + 6, T.chi2_bur + 6, T.chi2_law + 2
T.to_csv(ROOT / "results/halo_fits_per_galaxy.csv", index=False)

L = [f"{len(T)} galaxies. Y_disk prior: 0.1 dex around the SED value (or 0.5). Haloes: 3 parameters each; law: 1 parameter."]
L.append("\n| model | median RMS (km/s) | median chi2/N | sum AIC | galaxies with lowest AIC | median Y_disk |")
L.append("|---|---|---|---|---|---|")
best = T[["aic_nfw", "aic_bur", "aic_law"]].idxmin(axis=1)
for k, lab in (("nfw", "NFW + baryons"), ("bur", "Burkert + baryons"), ("law", "g = nu(g_bar/a0) g_bar")):
    L.append(f"| {lab} | {T[f'rms_{k}'].median():.1f} | {(T[f'chi2_{k}'] / T.N).median():.2f} | {T[f'aic_{k}'].sum():.0f} | {(best == f'aic_{k}').mean():.0%} | {T[f'y_{k}'].median():.2f} |")
d = T.chi2_bur - T.chi2_nfw
L.append(f"\nBurkert vs NFW: Burkert better in {(d < 0).mean():.0%} of galaxies; median Delta chi2 (Burkert - NFW) = {d.median():+.1f}; "
         f"late types (T>=8): {(d[T['T'] >= 8] < 0).mean():.0%} prefer Burkert, early (T<=4): {(d[T['T'] <= 4] < 0).mean():.0%}.")
# concentration-mass relation
ok = (T.M200 > 1e9) & (T.M200 < 1e14) & (T.c200 > 1) & (T.c200 < 100)
logc_pred = 0.905 - 0.101 * np.log10(T.M200[ok] / (1e12 / 0.7))
dc = np.log10(T.c200[ok]) - logc_pred
L.append(f"NFW concentration-mass relation, {ok.sum()} galaxies with 1e9 < M200 < 1e14: median offset from Dutton & Maccio 2014 {dc.median():+.2f} dex, scatter {dc.std():.2f} dex (LambdaCDM: 0.11 dex); "
         f"fraction with c200 below the relation by > 0.3 dex: {(dc < -0.3).mean():.0%}.")
# core surface density and core radius
okb = (T.r_0 > 0.15) & (T.r_0 < 150)
L.append(f"Burkert core surface density rho_0 r_0: median {T.core_sd[okb].median():.0f} Msun/pc^2, 16-84 % {T.core_sd[okb].quantile(.16):.0f}-{T.core_sd[okb].quantile(.84):.0f} (Donato et al. 2009: ~140, scatter 0.2 dex); "
         f"Spearman(rho_0 r_0, V_max) = {T[okb][['core_sd', 'Vmax']].corr('spearman').iloc[0, 1]:+.2f}.")
L.append(f"Burkert core radius / disc scale length: median {(T.r_0 / T.Rd)[okb].median():.2f}, 16-84 % {(T.r_0 / T.Rd)[okb].quantile(.16):.2f}-{(T.r_0 / T.Rd)[okb].quantile(.84):.2f}; "
         f"Spearman(r_0, Rd) = {T[okb][['r_0', 'Rd']].corr('spearman').iloc[0, 1]:+.2f}.")
L.append(f"Stellar M/L pulled by the fits: NFW median Y/Y0 = {(T.y_nfw / T.y0).median():.2f}, Burkert {(T.y_bur / T.y0).median():.2f}, law {(T.y_law / T.y0).median():.2f}.")

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(2, 2, figsize=(12, 9))
ax[0, 0].hist(np.clip(d, -60, 60), bins=40, color="tab:purple", alpha=0.7)
ax[0, 0].axvline(0, color="k"); ax[0, 0].set_xlabel("chi2(Burkert) - chi2(NFW), clipped"); ax[0, 0].set_title("Core vs cusp, per galaxy")
mm = np.logspace(9, 14, 50)
ax[0, 1].plot(T.M200[ok], T.c200[ok], "o", ms=4, alpha=0.7)
ax[0, 1].plot(mm, 10 ** (0.905 - 0.101 * np.log10(mm / (1e12 / 0.7))), "k-", label="Dutton & Maccio 2014")
ax[0, 1].fill_between(mm, 10 ** (0.905 - 0.101 * np.log10(mm / (1e12 / 0.7)) - 0.11), 10 ** (0.905 - 0.101 * np.log10(mm / (1e12 / 0.7)) + 0.11), color="k", alpha=0.15)
ax[0, 1].set_xscale("log"); ax[0, 1].set_yscale("log"); ax[0, 1].set_xlabel("M200 (Msun)"); ax[0, 1].set_ylabel("c200"); ax[0, 1].legend(fontsize=8); ax[0, 1].set_title("NFW concentration-mass")
ax[1, 0].plot(T.Vmax[okb], T.core_sd[okb], "o", ms=4, alpha=0.7, color="tab:green")
ax[1, 0].axhline(140, color="k", ls="--", label="Donato et al. 2009"); ax[1, 0].set_xscale("log"); ax[1, 0].set_yscale("log")
ax[1, 0].set_xlabel("V_max (km/s)"); ax[1, 0].set_ylabel("rho_0 r_0 (Msun/pc^2)"); ax[1, 0].legend(fontsize=8); ax[1, 0].set_title("Burkert core surface density")
for k, lab, col in (("nfw", "NFW (3 par.)", "tab:blue"), ("bur", "Burkert (3 par.)", "tab:green"), ("law", "law (1 par.)", "tab:red")):
    ax[1, 1].hist(np.log10(T[f"chi2_{k}"] / T.N), bins=np.arange(-1.5, 2.6, 0.1), alpha=0.45, color=col, label=lab)
ax[1, 1].set_xlabel("log10 chi2/N"); ax[1, 1].legend(fontsize=8); ax[1, 1].set_title("Goodness of fit")
fig.tight_layout(); fig.savefig(ROOT / "results/halo_fits.png", dpi=130)

with open(ROOT / "results/halo_fits.md", "w", encoding="utf-8") as fh:
    fh.write("# Halo fits on SPARC: NFW, Burkert and the acceleration law (2026-10-10)\n\n" + "\n".join(L) + "\n\nFigure: results/halo_fits.png; per galaxy: results/halo_fits_per_galaxy.csv\n")
print("\n".join(L))
print("wrote results/halo_fits.md, results/halo_fits.png")
