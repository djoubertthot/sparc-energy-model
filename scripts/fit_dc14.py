"""DC14 haloes (Di Cintio et al. 2014) on SPARC, and the cosmological consistency of all halo fits.

DC14: rho = rho_s / [ x^gamma (1 + x^alpha)^((beta - gamma)/alpha) ],  x = r / r_s, with
alpha, beta, gamma functions of X = log10(M_star / M_halo) (feedback-flattened cusps):
  alpha = 2.94 - log10[(10^(X+2.33))^-1.08 + (10^(X+2.33))^2.29]
  beta  = 4.23 + 1.34 X + 0.26 X^2
  gamma = -0.06 + log10[(10^(X+2.56))^-0.68 + 10^(X+2.56)]          valid for -4.1 < X < -1.3
Fitted like NFW and Burkert (rho_s, r_s, Y_disk with a 0.1 dex prior), X updated self-consistently.
Cosmological tests on every halo model: concentration c200 = r200 / r_-2 against Dutton & Maccio 2014,
stellar-to-halo mass ratio against Moster et al. 2013, baryon fraction against 0.157.
Writes results/dc14.md, results/dc14_per_galaxy.csv, results/dc14.png
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
H0 = 70e3 / 3.0857e22
RHO_CRIT = 3 * H0**2 / (8 * np.pi * G)
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")
XG = np.logspace(-3, 3.5, 500)
LNX = np.log(XG)


def dc14_params(X):
    X = np.clip(X, -4.1, -1.3)
    a = 2.94 - np.log10((10 ** (X + 2.33)) ** -1.08 + (10 ** (X + 2.33)) ** 2.29)
    b = 4.23 + 1.34 * X + 0.26 * X**2
    g = -0.06 + np.log10((10 ** (X + 2.56)) ** -0.68 + 10 ** (X + 2.56))
    return a, b, g


def mass_profile(a, b, g):
    """F(x) = int_0^x x'^2 rho~(x') dx' on the grid, rho~ = x^-g (1+x^a)^(-(b-g)/a)."""
    rho = XG ** (-g) * (1 + XG**a) ** (-(b - g) / a)
    integrand = XG**3 * rho                       # x^2 rho dx = x^3 rho dlnx
    F = np.concatenate([[integrand[0] / (3 - g)], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(LNX))])
    F[1:] += integrand[0] / (3 - g)
    return F


def m_enclosed(r, rho_s, r_s, F):
    return 4 * np.pi * rho_s * r_s**3 * np.interp(np.log(r / r_s), LNX, F)


def r200_of(rho_s, r_s, F):
    mean = 3 * rho_s * F / XG**3
    i = np.where(mean < 200 * RHO_CRIT)[0]
    if len(i) == 0:
        return np.nan
    j = i[0]
    return XG[j] * r_s


def moster(mh):
    N, M1, b, g = 0.0351, 10**11.59, 1.376, 0.608
    return 2 * N / ((mh / M1) ** -b + (mh / M1) ** g)


def fit_dc14(df, y0, l36):
    r = df.R.values * KPC
    vobs, e = df.Vobs.values * 1e3, np.maximum(df.eVobs.values, 0.5) * 1e3
    X = -1.7
    best = None
    for _ in range(3):
        a, b, g = dc14_params(X)
        F = mass_profile(a, b, g)

        def resid(p):
            lrho, lr, ly = p
            v2 = fb.baryonic2(df, 10**ly, 1.0, 1.4 * 10**ly) + G * m_enclosed(r, 10**lrho * MSUN / PC**3, 10**lr * KPC, F) / r
            return np.concatenate([(vobs - np.sqrt(np.maximum(v2, 0))) / e, [np.log10(10**ly / y0) / 0.1]])

        best = None
        for lrho0 in (-2.5, -1.5, -0.5):
            for lr0 in (0.0, 0.7, 1.4):
                try:
                    res = least_squares(resid, [lrho0, lr0, np.log10(y0)], bounds=([-5, -1, np.log10(0.1)], [1.5, 2.7, np.log10(1.5)]))
                except Exception:
                    continue
                if best is None or res.cost < best.cost:
                    best = res
        lrho, lr, ly = best.x
        rho_s, r_s = 10**lrho * MSUN / PC**3, 10**lr * KPC
        r200 = r200_of(rho_s, r_s, F)
        m200 = 4 / 3 * np.pi * r200**3 * 200 * RHO_CRIT if np.isfinite(r200) else np.nan
        mstar = 10**ly * l36 * 1e9 * MSUN
        X_new = np.log10(mstar / m200) if np.isfinite(m200) else X
        if abs(X_new - X) < 0.05:
            X = X_new
            break
        X = X_new
    a, b, g = dc14_params(X)
    chi2 = np.sum(best.fun[:-1] ** 2)
    rm2 = r_s * ((2 - g) / (b - 2)) ** (1 / a)
    return dict(chi2_dc14=chi2, rho_s_dc14=rho_s / (MSUN / PC**3), r_s_dc14=r_s / KPC, y_dc14=10**ly, X=X,
                M200_dc14=m200 / MSUN if np.isfinite(m200) else np.nan, c200_dc14=r200 / rm2 if np.isfinite(r200) else np.nan,
                alpha=a, beta=b, gamma=g)


def burkert_m200(rho0, r0):
    rho0 = rho0 * MSUN / PC**3
    r0 = r0 * KPC
    rr = np.logspace(-1, 4, 2000) * KPC
    x = rr / r0
    M = 2 * np.pi * rho0 * r0**3 * (np.log(1 + x) - np.arctan(x) + 0.5 * np.log(1 + x**2))
    mean = M / (4 / 3 * np.pi * rr**3)
    i = np.where(mean < 200 * RHO_CRIT)[0]
    return (4 / 3 * np.pi * rr[i[0]] ** 3 * 200 * RHO_CRIT / MSUN) if len(i) else np.nan


def main():
    H = pd.read_csv(ROOT / "results/halo_fits_per_galaxy.csv").set_index("Galaxy")
    rows = []
    for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
        df = fb.load_rotmod(f)
        if len(df) < 5:
            continue
        n = f.stem.replace("_rotmod", "")
        c = cat.loc[n]
        y0 = float(c.Y36) if pd.notna(c.Y36) else 0.5
        l36 = float(c.sparc_L36_1e9_Lsun)
        d = fit_dc14(df, y0, l36)
        h = H.loc[n]
        d.update(Galaxy=n, T=int(c["T"]), N=len(df), L36=l36, Mgas=1.33 * float(c.sparc_MHI_1e9_Msun) * 1e9,
                 chi2_nfw=h.chi2_nfw, chi2_bur=h.chi2_bur, chi2_law=h.chi2_law, M200_nfw=h.M200, c200_nfw=h.c200, y_nfw=h.y_nfw,
                 M200_bur=burkert_m200(h.rho_0, h.r_0), y_bur=h.y_bur)
        rows.append(d)
    T = pd.DataFrame(rows)
    T.to_csv(ROOT / "results/dc14_per_galaxy.csv", index=False)

    L = [f"{len(T)} galaxies. DC14 converged X (log M*/M200) median {T.X.median():.2f}; inner slope gamma median {T.gamma.median():.2f} (NFW: 1)."]
    aic = pd.DataFrame({"NFW": T.chi2_nfw + 6, "Burkert": T.chi2_bur + 6, "DC14": T.chi2_dc14 + 6, "law": T.chi2_law + 2})
    best = aic.idxmin(axis=1)
    L.append("\n| model | median chi2/N | sum AIC | lowest AIC in |\n|---|---|---|---|")
    for k, col in (("NFW", "chi2_nfw"), ("Burkert", "chi2_bur"), ("DC14", "chi2_dc14"), ("law", "chi2_law")):
        L.append(f"| {k} | {(T[col] / T.N).median():.2f} | {aic[k].sum():.0f} | {(best == k).mean():.0%} |")
    L.append(f"DC14 vs NFW: DC14 better in {(T.chi2_dc14 < T.chi2_nfw).mean():.0%}; DC14 vs Burkert: DC14 better in {(T.chi2_dc14 < T.chi2_bur).mean():.0%}.")

    L.append("\n## Cosmological consistency (constrained fits only: 2.5 < c200 < 60, 1e9 < M200 < 3e12 Msun)\n")
    L.append("| halo | N | c-M offset (dex) | c-M scatter (dex), LCDM 0.11 | SHM offset log(M*/M200) vs Moster 2013 (dex) | SHM scatter | median (M*+Mgas)/M200, limit 0.157 |")
    L.append("|---|---|---|---|---|---|---|")
    for k, (mc, cc, yc) in {"NFW": ("M200_nfw", "c200_nfw", "y_nfw"), "DC14": ("M200_dc14", "c200_dc14", "y_dc14")}.items():
        ok = (T[cc] > 2.5) & (T[cc] < 60) & (T[mc] > 1e9) & (T[mc] < 3e12)
        dc = np.log10(T[cc][ok]) - (0.905 - 0.101 * np.log10(T[mc][ok] / (1e12 / 0.7)))
        mstar = T[yc][ok] * T.L36[ok] * 1e9
        ds = np.log10(mstar / T[mc][ok]) - np.log10(moster(T[mc][ok]))
        fb_ = (mstar + T.Mgas[ok]) / T[mc][ok]
        L.append(f"| {k} | {ok.sum()} | {dc.median():+.2f} | {dc.std():.2f} | {ds.median():+.2f} | {ds.std():.2f} | {fb_.median():.3f} |")
    okb = (T.M200_bur > 1e9) & (T.M200_bur < 3e12)
    mstar = T.y_bur[okb] * T.L36[okb] * 1e9
    ds = np.log10(mstar / T.M200_bur[okb]) - np.log10(moster(T.M200_bur[okb]))
    L.append(f"| Burkert (M200 by extrapolation) | {okb.sum()} | - | - | {ds.median():+.2f} | {ds.std():.2f} | {((mstar + T.Mgas[okb]) / T.M200_bur[okb]).median():.3f} |")

    import matplotlib  # noqa: E402
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    mm = np.logspace(9, 13, 60)
    for k, (mc, cc, yc), col in zip(("NFW", "DC14"), (("M200_nfw", "c200_nfw", "y_nfw"), ("M200_dc14", "c200_dc14", "y_dc14")), ("tab:blue", "tab:red")):
        ok = (T[cc] > 2.5) & (T[cc] < 60) & (T[mc] > 1e9) & (T[mc] < 3e12)
        ax[0].plot(T[mc][ok], T[cc][ok], "o", ms=4, alpha=0.6, color=col, label=k)
        ax[1].plot(T[mc][ok], T[yc][ok] * T.L36[ok] * 1e9 / T[mc][ok], "o", ms=4, alpha=0.6, color=col, label=k)
    ax[0].plot(mm, 10 ** (0.905 - 0.101 * np.log10(mm / (1e12 / 0.7))), "k-", label="Dutton & Maccio 2014")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("M200 (Msun)"); ax[0].set_ylabel("c200"); ax[0].legend(fontsize=8); ax[0].set_title("Concentration-mass")
    ax[1].plot(mm, moster(mm), "k-", label="Moster et al. 2013")
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("M200 (Msun)"); ax[1].set_ylabel("M* / M200"); ax[1].legend(fontsize=8); ax[1].set_title("Stellar-to-halo mass")
    fig.tight_layout(); fig.savefig(ROOT / "results/dc14.png", dpi=130)
    with open(ROOT / "results/dc14.md", "w", encoding="utf-8") as fh:
        fh.write("# DC14 haloes and cosmological consistency of the halo fits (2026-10-10)\n\n" + "\n".join(L) + "\n\nFigure: results/dc14.png\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
