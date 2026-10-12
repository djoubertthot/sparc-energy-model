"""Model v2: the formulation that survives the audit (see docs/THEORY.md).

Baryons   : V_bar^2 = Vgas|Vgas| + Yd Vdisk|Vdisk| + Yb Vbul|Vbul|,  g_bar = V_bar^2 / r
Law       : g_obs = nu(g_bar / a0) g_bar,     a0 = 1.2e-10 m/s^2
            nu_RAR(y)    = 1 / (1 - exp(-sqrt(y)))              (McGaugh, Lelli & Schombert 2016)
            nu_simple(y) = 1/2 + sqrt(1/4 + 1/y)                 (Famaey & Binney 2005)
Equivalent: M_eff(<r) = M_bar(<r) nu(g_bar/a0): the repository's "amplification as a function
            of a_bary", applied to the enclosed mass, with exponent 1/2 in the deep regime.
Derived   : dark-equivalent mass M_x(<r) = (nu - 1) M_bar(<r); BTFR V_flat^4 = G M_bar a0.

Runs the law on every SPARC curve with zero free parameters per galaxy, in three
baryonic settings (fixed Y, per-galaxy Y36 from the enriched table, Y free with a
0.1 dex prior), reports RMS / bias / outer residual / by type, the implied dark
component, and the BTFR prediction. Writes results/model_v2.md and results/model_v2_per_galaxy.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
PC = KPC / 1e3
A0 = 1.2e-10
C = 2.99792458e8
H0 = 70e3 / 3.0857e22

cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")


def nu_rar(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def nu_simple(y):
    y = np.maximum(y, 1e-12)
    return 0.5 + np.sqrt(0.25 + 1.0 / y)


NU = {"RAR": nu_rar, "simple": nu_simple}


def model(df, yd, yb, nu):
    r = df.R.values * KPC
    vb2 = fb.baryonic2(df, yd, 1.0, yb)
    gbar = np.maximum(vb2, 1.0) / r
    gobs = gbar * nu(gbar / A0)
    return np.sqrt(gobs * r) / 1e3, gbar, vb2


def metrics(df, v, rd):
    res = df.Vobs.values - v
    w = 1 / np.maximum(df.eVobs.values, 0.5) ** 2
    R = df.R.values
    return dict(rms=np.sqrt(np.mean(res**2)), wrms=np.sqrt(np.average(res**2, weights=w)), bias=res.mean(),
                outer=res[R > 20].mean() if (R > 20).any() else np.nan, inner=res[R < 2 * rd].mean() if (R < 2 * rd).any() else np.nan,
                chi2n=np.mean((res / np.maximum(df.eVobs.values, 0.5)) ** 2))


gal = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    n = f.stem.replace("_rotmod", "")
    c = cat.loc[n]
    gal.append(dict(name=n, df=df, T=int(c["T"]), rd=float(c.Rdisk_kpc), y36=float(c.Y36) if pd.notna(c.Y36) else np.nan,
                    vflat=float(c.Vflat_kms) if pd.notna(c.Vflat_kms) and c.Vflat_kms > 0 else np.nan,
                    l36=float(c.sparc_L36_1e9_Lsun), mhi=float(c.sparc_MHI_1e9_Msun)))

rows, per = [], []
for nu_name, nu in NU.items():
    for setting in ("Y fixed 0.5 / 0.7", "Y36 per galaxy", "Y free, 0.1 dex prior"):
        ms = []
        for g in gal:
            df = g["df"]
            if setting == "Y fixed 0.5 / 0.7":
                yd = 0.5
                v, gbar, vb2 = model(df, yd, 0.7, nu)
            elif setting == "Y36 per galaxy":
                yd = g["y36"] if np.isfinite(g["y36"]) else 0.5
                v, gbar, vb2 = model(df, yd, 1.4 * yd, nu)
            else:
                best = None
                for yd in np.linspace(0.2, 1.0, 17):
                    v, gbar, vb2 = model(df, yd, 1.4 * yd, nu)
                    chi2 = np.sum(((df.Vobs.values - v) / np.maximum(df.eVobs.values, 0.5)) ** 2) + (np.log10(yd / 0.5) / 0.1) ** 2
                    if best is None or chi2 < best[0]:
                        best = (chi2, yd, v, gbar, vb2)
                _, yd, v, gbar, vb2 = best
            m = metrics(df, v, g["rd"])
            m.update(T=g["T"])
            ms.append(m)
            if nu_name == "RAR" and setting == "Y fixed 0.5 / 0.7":
                R = df.R.values[-1] * KPC
                mbar = max(vb2[-1], 1.0) * R / G
                nu_last = nu(gbar[-1] / A0)
                mx = (nu_last - 1) * mbar
                per.append(dict(Galaxy=g["name"], T=g["T"], N=len(df), RMS=round(m["rms"], 2), bias=round(m["bias"], 2),
                                g_bar_last=gbar[-1], nu_last=nu_last, M_bar_Msun=mbar / MSUN, M_x_Msun=mx / MSUN,
                                Sigma_x_Msun_pc2=mx / (np.pi * R**2) / (MSUN / PC**2),
                                V_flat=g["vflat"], V_btfr=(G * (0.5 * g["l36"] + 1.33 * g["mhi"]) * 1e9 * MSUN * A0) ** 0.25 / 1e3))
        d = pd.DataFrame(ms)
        row = dict(law=nu_name, baryons=setting, N=len(d), rms=d.rms.median(), wrms=d.wrms.median(), bias=d.bias.median(),
                   outer=d.outer.median(), inner=d.inner.median(), chi2n=d.chi2n.median(),
                   early=d.rms[d["T"] <= 4].median(), inter=d.rms[(d["T"] >= 5) & (d["T"] <= 7)].median(), late=d.rms[d["T"] >= 8].median())
        rows.append(row)
        print(f"{nu_name:7s} {setting:24s} N={len(d)} RMS={row['rms']:5.1f} wRMS={row['wrms']:5.1f} bias={row['bias']:+5.1f} outer={row['outer']:+5.1f} "
              f"inner={row['inner']:+5.1f} chi2/N={row['chi2n']:5.1f} | early {row['early']:5.1f} inter {row['inter']:5.1f} late {row['late']:5.1f}")

T = pd.DataFrame(rows)
P = pd.DataFrame(per)
P.to_csv(ROOT / "results/model_v2_per_galaxy.csv", index=False)

# derived dark component and BTFR
b = P.dropna(subset=["V_flat"])
btfr_res = np.log10(b.V_flat / b.V_btfr)
print(f"\nderived dark component at R_last: M_x/M_bar median {(P.M_x_Msun / P.M_bar_Msun).median():.2f}; "
      f"Sigma_x median {P.Sigma_x_Msun_pc2.median():.0f} Msun/pc^2 (p16 {P.Sigma_x_Msun_pc2.quantile(.16):.0f}, p84 {P.Sigma_x_Msun_pc2.quantile(.84):.0f})")
print(f"BTFR prediction V_flat = (G M_bar a0)^(1/4), {len(b)} galaxies with V_flat: median log(V_obs/V_pred) = {btfr_res.median():+.3f}, scatter = {btfr_res.std():.3f} dex")
print(f"a0 = {A0:.2e} ; c H0 / 2 pi = {C * H0 / (2 * np.pi):.2e} (H0 = 70) ; ratio = {A0 / (C * H0 / (2 * np.pi)):.2f}")

with open(ROOT / "results/model_v2.md", "w", encoding="utf-8") as fh:
    fh.write("# Model v2 on SPARC (2026-10-10): g_obs = nu(g_bar/a0) g_bar, zero free parameters per galaxy\n\n")
    fh.write("| law | baryons | N | RMS | wRMS | bias | outer R>20 | inner <2Rd | chi2/N | early | inter | late |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        fh.write(f"| {r['law']} | {r['baryons']} | {r['N']} | {r['rms']:.1f} | {r['wrms']:.1f} | {r['bias']:+.1f} | {r['outer']:+.1f} | {r['inner']:+.1f} | {r['chi2n']:.1f} | {r['early']:.1f} | {r['inter']:.1f} | {r['late']:.1f} |\n")
    fh.write(f"\nDerived dark-equivalent component at R_last (RAR, Y fixed): M_x/M_bar median {(P.M_x_Msun / P.M_bar_Msun).median():.2f}; "
             f"Sigma_x median {P.Sigma_x_Msun_pc2.median():.0f} Msun/pc^2 (16-84 %: {P.Sigma_x_Msun_pc2.quantile(.16):.0f}-{P.Sigma_x_Msun_pc2.quantile(.84):.0f}).\n")
    fh.write(f"BTFR prediction V_flat = (G M_bar a0)^(1/4) on {len(b)} galaxies: median offset {btfr_res.median():+.3f} dex, scatter {btfr_res.std():.3f} dex.\n")
    fh.write(f"a0 = {A0:.2e} m/s^2 ; c H0 / 2 pi = {C * H0 / (2 * np.pi):.2e} m/s^2 for H0 = 70 km/s/Mpc.\n")
print("wrote results/model_v2.md, results/model_v2_per_galaxy.csv")
