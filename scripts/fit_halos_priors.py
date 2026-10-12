"""NFW and DC14 haloes with LambdaCDM priors (after Katz et al. 2017).

Parameters: log M200, log c200, log Y_disk. Priors: c200 from Dutton & Maccio 2014 with 0.11 dex
scatter; M* / M200 from Moster et al. 2013 with 0.15 dex scatter; Y_disk 0.1 dex around the SED
value. Chi^2 includes the prior terms. If the halo model is the LambdaCDM one, the fits should stay
good under these priors; a profile that needs cosmologically unlikely parameters degrades.
Writes results/halo_priors.md, results/halo_priors_per_galaxy.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402
import fit_dc14 as dc  # noqa: E402  (loads its own results on import; reuse the functions)

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
PC = KPC / 1e3
H0 = dc.H0
RHO_CRIT = dc.RHO_CRIT
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")


def logc_dm14(m200):
    return 0.905 - 0.101 * np.log10(m200 / (1e12 / 0.7))


def halo_mass_nfw(r, m200, c):
    r200 = (3 * m200 * MSUN / (4 * np.pi * 200 * RHO_CRIT)) ** (1 / 3)
    rs = r200 / c
    x = r / rs
    f = lambda y: np.log(1 + y) - y / (1 + y)
    return m200 * MSUN * f(x) / f(c)


def halo_mass_dc14(r, m200, c, X):
    a, b, g = dc.dc14_params(X)
    F = dc.mass_profile(a, b, g)
    r200 = (3 * m200 * MSUN / (4 * np.pi * 200 * RHO_CRIT)) ** (1 / 3)
    rm2_over_rs = ((2 - g) / (b - 2)) ** (1 / a)
    rs = r200 / (c * rm2_over_rs)
    F200 = np.interp(np.log(r200 / rs), dc.LNX, F)
    rho_s = m200 * MSUN / (4 * np.pi * rs**3 * F200)
    return dc.m_enclosed(r, rho_s, rs, F)


def fit(df, y0, l36, kind):
    r = df.R.values * KPC
    vobs, e = df.Vobs.values * 1e3, np.maximum(df.eVobs.values, 0.5) * 1e3

    def resid(p):
        lm, lc, ly = p
        m200, c, yd = 10**lm, 10**lc, 10**ly
        mstar = yd * l36 * 1e9
        X = np.log10(mstar / m200)
        mh = halo_mass_nfw(r, m200, c) if kind == "nfw" else halo_mass_dc14(r, m200, c, X)
        v2 = fb.baryonic2(df, yd, 1.0, 1.4 * yd) + G * mh / r
        pri = [(lc - logc_dm14(m200)) / 0.11, (np.log10(mstar / m200) - np.log10(dc.moster(m200))) / 0.15, np.log10(yd / y0) / 0.1]
        return np.concatenate([(vobs - np.sqrt(np.maximum(v2, 0))) / e, pri])

    best = None
    for lm0 in (10.0, 11.0, 12.0):
        for lc0 in (0.7, 1.0, 1.3):
            try:
                res = least_squares(resid, [lm0, lc0, np.log10(y0)], bounds=([8.5, 0.0, np.log10(0.1)], [13.5, 2.0, np.log10(1.5)]))
            except Exception:
                continue
            if best is None or res.cost < best.cost:
                best = res
    chi2_data = np.sum(best.fun[:-3] ** 2)
    chi2_tot = np.sum(best.fun**2)
    lm, lc, ly = best.x
    return chi2_data, chi2_tot, 10**lm, 10**lc, 10**ly


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
    cn, tn, mn, ccn, yn = fit(df, y0, l36, "nfw")
    cd, td, md, ccd, ydd = fit(df, y0, l36, "dc14")
    rows.append(dict(Galaxy=n, T=int(c["T"]), N=len(df), chi2_nfw_prior=cn, chi2tot_nfw_prior=tn, M200_nfw_p=mn, c_nfw_p=ccn,
                     chi2_dc14_prior=cd, chi2tot_dc14_prior=td, M200_dc14_p=md, c_dc14_p=ccd,
                     chi2_nfw_free=H.loc[n, "chi2_nfw"], chi2_bur_free=H.loc[n, "chi2_bur"], chi2_law=H.loc[n, "chi2_law"]))
T = pd.DataFrame(rows)
T.to_csv(ROOT / "results/halo_priors_per_galaxy.csv", index=False)
L = [f"{len(T)} galaxies; priors: c(M) +/- 0.11 dex, M*/M200 Moster +/- 0.15 dex, Y_disk +/- 0.1 dex. chi2 (data only) quoted."]
L.append("\n| model | median chi2/N | galaxies with chi2/N < 2 | sum chi2 |\n|---|---|---|---|")
for k, col in (("NFW, free (rho_s, r_s)", "chi2_nfw_free"), ("Burkert, free", "chi2_bur_free"), ("NFW, LCDM priors", "chi2_nfw_prior"),
               ("DC14, LCDM priors", "chi2_dc14_prior"), ("acceleration law", "chi2_law")):
    L.append(f"| {k} | {(T[col] / T.N).median():.2f} | {((T[col] / T.N) < 2).mean():.0%} | {T[col].sum():.0f} |")
L.append(f"\nUnder LCDM priors, DC14 beats NFW in {(T.chi2_dc14_prior < T.chi2_nfw_prior).mean():.0%} of galaxies (late types: "
         f"{(T.chi2_dc14_prior[T['T'] >= 8] < T.chi2_nfw_prior[T['T'] >= 8]).mean():.0%}).")
L.append(f"Degradation imposed by the priors, median chi2 ratio prior/free: NFW {(T.chi2_nfw_prior / T.chi2_nfw_free).median():.2f}, DC14 (vs free Burkert) {(T.chi2_dc14_prior / T.chi2_bur_free).median():.2f}.")
L.append(f"Law vs DC14 with priors: law better in {(T.chi2_law < T.chi2_dc14_prior).mean():.0%} of galaxies (1 parameter vs 3 with priors).")
with open(ROOT / "results/halo_priors.md", "w", encoding="utf-8") as fh:
    fh.write("# Haloes under LambdaCDM priors (2026-10-10)\n\n" + "\n".join(L) + "\n")
print("\n".join(L))
