"""RMS tests of the tag-1.0 baryonic model (scripts/final_baryonic_model.py) on the SPARC curves.

Reproduces every configuration the documents quote a number for (options table,
intermediate CGM versions, final parameters), measures the criteria the documents
use (RMS, weighted RMS, mean bias, outer residual beyond 20 kpc, inner residual
inside 2 R_disk, RMS by type), scans the CGM and baryonic parameters, and
quantifies the hot-gas density the CGM term implies. Two references for scale.

Outputs (results/): rms_tests_v1.md, rms_tests_v1.csv, rms_per_galaxy_v1.csv, grid_cgm_v1.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")

galaxies = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    name = f.stem.replace("_rotmod", "")
    c = cat.loc[name]
    galaxies.append(dict(name=name, df=df, T=int(c["T"]), Q=int(c.Q), inc=float(c.Inc_deg),
                         l36=float(c.sparc_L36_1e9_Lsun), mhi=float(c.sparc_MHI_1e9_Msun), rd=float(c.Rdisk_kpc)))
N_ALL = len(galaxies)


def metrics(out, g):
    res = out["residual"].values
    r = out["R"].values
    w = 1.0 / np.maximum(out["eVobs"].values, 0.1) ** 2
    outer = res[r > 20.0]
    inner = res[r < 2.0 * g["rd"]]
    return dict(rms=float(np.sqrt(np.mean(res**2))), wrms=float(np.sqrt(np.average(res**2, weights=w))),
                bias=float(np.mean(res)), outer=float(np.mean(outer)) if len(outer) else np.nan,
                inner=float(np.mean(inner)) if len(inner) else np.nan)


def model(g, **kw):
    return metrics(fb.compute(g["df"], g["l36"], g["mhi"], **kw), g)


def mond(g):
    df = g["df"]
    R = df.R.values * fb.KPC
    gn = np.maximum(fb.baryonic2(df, 0.5, 1.0, 0.7), 1.0) / R
    gm = gn / 2 + np.sqrt(gn**2 / 4 + gn * 1.2e-10)
    out = df.copy()
    out["residual"] = df.Vobs.values - np.sqrt(gm * R) / 1e3
    return metrics(out, g)


def run(fn, subset=None, **kw):
    ms = [fn(g, **kw) for g in galaxies if subset is None or subset(g)]
    d = pd.DataFrame(ms)
    return dict(N=len(d), median=d.rms.median(), wmedian=d.wrms.median(), mean=d.rms.mean(),
                p16=d.rms.quantile(0.16), p84=d.rms.quantile(0.84),
                bias=d.bias.median(), outer=d.outer.median(), inner=d.inner.median())


rows = []


def add(test, group, s, claimed=""):
    rows.append(dict(group=group, test=test, N=s["N"], median=round(s["median"], 1), wmedian=round(s["wmedian"], 1),
                     mean=round(s["mean"], 1), p16=round(s["p16"], 1), p84=round(s["p84"], 1),
                     bias=round(s["bias"], 1), outer_R20=round(s["outer"], 1), inner_2Rd=round(s["inner"], 1), claimed=claimed))
    print(f"{group:11s} {test:58s} N={s['N']:3d} RMS={s['median']:6.1f} w={s['wmedian']:6.1f} bias={s['bias']:+6.1f} "
          f"outer={s['outer']:+6.1f} inner={s['inner']:+6.1f}  claimed={claimed}")


FINAL = dict(ystar=0.88, f_gas=2.7, m_factor=3.8, r_c=6.0, beta=0.68)

# 1. final parameters under the three CGM normalisations ----------------------
for norm, label in [("last_std", "M_CGM(<R_last) = 3.8 x (0.5 L36 + 1.33 MHI)"),
                    ("last_model", "M_CGM(<R_last) = 3.8 x (0.88 L36 + 2.7 x 1.33 MHI)"),
                    ("r25_std", "M_CGM(<25 kpc) = 3.8 x (0.5 L36 + 1.33 MHI)")]:
    add(f"final 1.0, {label}", "final", run(model, norm=norm, **FINAL), "21")

# 2. options table, pure baryonic (docs/OPTIONS_COMPARISON.md, MODEL_OPTION_A.md) ----
for ys, fg, claimed in [(0.5, 1.0, "56"), (1.0, 1.0, "32 / 31.7"), (1.0, 2.0, "28"), (1.0, 2.5, "26"),
                        (1.2, 2.0, "25"), (1.5, 2.0, "26"), (0.88, 2.7, "")]:
    add(f"no CGM, Y*={ys}, gas x{fg}", "options", run(model, ystar=ys, f_gas=fg, m_factor=0.0), claimed)
add("no CGM, Y*=0.5 disk / 0.7 bulge, gas x1 (SPARC standard)", "options", run(model, ystar=0.5, f_gas=1.0, m_factor=0.0, ybul=0.7), "")

# 3. intermediate CGM versions quoted in the docs (norm last_std) ------------------
for label, kw, claimed in [
        ("Y*=1.0, gas x2, CGM 3.5, r_c=20, beta=0.6", dict(ystar=1.0, f_gas=2.0, m_factor=3.5, r_c=20.0, beta=0.6), "26"),
        ("Y*=1.0, gas x2, CGM 3.5, r_c=20, beta=0.7 (REVERT)", dict(ystar=1.0, f_gas=2.0, m_factor=3.5, r_c=20.0, beta=0.7), "~26, outer +30"),
        ("Y*=1.0, gas x2, M=4.0, r_c=10, beta=0.75 (FINAL_PHYSICAL)", dict(ystar=1.0, f_gas=2.0, m_factor=4.0, r_c=10.0, beta=0.75), "23"),
        ("Y*=1.0, gas x2, CGM 4.5, r_c=30, beta=0.5 (FINAL_ADJUSTMENTS)", dict(ystar=1.0, f_gas=2.0, m_factor=4.5, r_c=30.0, beta=0.5), "24, outer +7"),
        ("Y*=0.9, gas 2.5, M=4.0, r_c=6, beta=0.65 (EVEN_FINER)", dict(ystar=0.9, f_gas=2.5, m_factor=4.0, r_c=6.0, beta=0.65), "~22")]:
    add(label, "versions", run(model, norm="last_std", **kw), claimed)

# 4. grid over CGM parameters at Y*=0.88, gas 2.7 (norm last_std) -----------------
grid = []
for m in [1.0, 2.0, 3.0, 3.8, 5.0, 7.0]:
    for rc in [2.0, 6.0, 10.0, 20.0, 30.0]:
        for b in [0.5, 0.68, 0.8, 1.0]:
            s = run(model, norm="last_std", ystar=0.88, f_gas=2.7, m_factor=m, r_c=rc, beta=b)
            grid.append(dict(M=m, r_c=rc, beta=b, median_rms=round(s["median"], 2), bias=round(s["bias"], 2), outer=round(s["outer"], 2)))
grid = pd.DataFrame(grid)
grid.to_csv(ROOT / "results/grid_cgm_v1.csv", index=False)
best = grid.sort_values("median_rms").iloc[0]
print("\n--- CGM grid, 5 best by median RMS ---")
print(grid.sort_values("median_rms").head(5).to_string(index=False))
add(f"grid best CGM: M={best.M}, r_c={best.r_c}, beta={best.beta}", "grid",
    run(model, norm="last_std", ystar=0.88, f_gas=2.7, m_factor=best.M, r_c=best.r_c, beta=best.beta), "")

# 5. Y* and gas scan at the final CGM --------------------------------------------
print("\n--- Y* x gas scan at CGM 3.8 / 6 / 0.68 (median RMS) ---")
yg = []
for ys in [0.5, 0.7, 0.88, 1.0, 1.2]:
    for fg in [1.0, 1.5, 2.0, 2.7]:
        s = run(model, norm="last_std", ystar=ys, f_gas=fg, m_factor=3.8, r_c=6.0, beta=0.68)
        yg.append(dict(ystar=ys, gas=fg, median=round(s["median"], 1)))
yg = pd.DataFrame(yg)
print(yg.pivot(index="ystar", columns="gas", values="median").to_string())
by = yg.sort_values("median").iloc[0]
add(f"Y* x gas best at final CGM: Y*={by.ystar}, gas x{by.gas}", "grid",
    run(model, norm="last_std", ystar=by.ystar, f_gas=by.gas, m_factor=3.8, r_c=6.0, beta=0.68), "")

# 6. by type and subsamples (final, last_std) -------------------------------------
add("final, early types T<=4", "by type", run(model, subset=lambda g: g["T"] <= 4, norm="last_std", **FINAL), "47 (Y*=1, gas x2)")
add("final, intermediate 5<=T<=7", "by type", run(model, subset=lambda g: 5 <= g["T"] <= 7, norm="last_std", **FINAL), "31")
add("final, late types T>=8", "by type", run(model, subset=lambda g: g["T"] >= 8, norm="last_std", **FINAL), "22")
add("final, Q<=2 and 40<i<75", "subsample", run(model, subset=lambda g: g["Q"] <= 2 and 40 < g["inc"] < 75, norm="last_std", **FINAL), "")

# 7. references ---------------------------------------------------------------------
add("MOND simple mu, a0=1.2e-10, Y*=0.5/0.7, 0 free parameters", "reference", run(mond), "~15 (archive_v1/docs/MODEL_COMPARISON.md)")

# 8. implied hot-gas density of the CGM term (final, last_std) ----------------------
dens = []
for g in galaxies:
    out = fb.compute(g["df"], g["l36"], g["mhi"], norm="last_std", **FINAL)
    rl = g["df"].R.values[-1] * fb.KPC
    m = out["M_cgm_Msun"].values[-1] * fb.MSUN
    n_mean = m / (4 / 3 * np.pi * rl**3) / (1.4 * 1.6726e-27) * 1e-6   # mean H density, cm^-3
    dm = n_mean * 1.2 * (rl / 3.0857e16)                                 # electron DM through R_last, pc cm^-3 (n_e ~ 1.2 n_H)
    dens.append(dict(Galaxy=g["name"], R_last_kpc=g["df"].R.values[-1], M_cgm_Msun=out["M_cgm_Msun"].values[-1],
                     M_std_Msun=(0.5 * g["l36"] + 1.33 * g["mhi"]) * 1e9, n_mean_cm3=n_mean, DM_pc_cm3=dm))
dens = pd.DataFrame(dens)
print("\n--- CGM term implications (final, last_std) ---")
print(f"median M_cgm(<R_last) = {dens.M_cgm_Msun.median():.2e} Msun ; median mean n_H inside R_last = {dens.n_mean_cm3.median():.2e} cm^-3 "
      f"(p16 {dens.n_mean_cm3.quantile(.16):.1e}, p84 {dens.n_mean_cm3.quantile(.84):.1e})")
print(f"median dispersion measure through R_last = {dens.DM_pc_cm3.median():.0f} pc cm^-3 ; NGC2403: n={dens.set_index('Galaxy').loc['NGC2403','n_mean_cm3']:.2e}, "
      f"M_cgm={dens.set_index('Galaxy').loc['NGC2403','M_cgm_Msun']:.2e}")

# outputs ----------------------------------------------------------------------------
per = []
for g in galaxies:
    o = fb.compute(g["df"], g["l36"], g["mhi"], norm="last_std", **FINAL)
    m = metrics(o, g)
    per.append(dict(Galaxy=g["name"], N=len(g["df"]), T=g["T"], Q=g["Q"], Inc=g["inc"], R_last=g["df"].R.values[-1],
                    RMS_final=round(m["rms"], 2), bias=round(m["bias"], 2), outer=round(m["outer"], 2) if np.isfinite(m["outer"]) else "",
                    inner=round(m["inner"], 2) if np.isfinite(m["inner"]) else "", RMS_MOND=round(mond(g)["rms"], 2),
                    M_cgm_Msun=f"{o['M_cgm_Msun'].values[-1]:.3e}"))
pd.DataFrame(per).to_csv(ROOT / "results/rms_per_galaxy_v1.csv", index=False)
pd.DataFrame(rows).to_csv(ROOT / "results/rms_tests_v1.csv", index=False)
with open(ROOT / "results/rms_tests_v1.md", "w", encoding="utf-8") as fh:
    fh.write("# RMS tests of the tag-1.0 baryonic model (2026-10-10)\n\n")
    fh.write("Model: scripts/final_baryonic_model.py, Y*=0.88, gas x2.7, CGM beta profile M=3.8, r_c=6 kpc, beta=0.68, E0=0 (archive_v1/docs/CONSOLIDATED_NOTES.md). ")
    fh.write(f"Sample: {N_ALL} SPARC galaxies with at least 5 points. RMS in km/s; bias = median over galaxies of the mean residual Vobs-Vmodel; ")
    fh.write("outer_R20 = same for points beyond 20 kpc; inner_2Rd = same for points inside 2 R_disk; claimed = value quoted in the docs.\n\n")
    fh.write("| group | test | N | median | wmedian | mean | p16 | p84 | bias | outer_R20 | inner_2Rd | claimed |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        fh.write("| " + " | ".join(str(r[k]) for k in ["group", "test", "N", "median", "wmedian", "mean", "p16", "p84", "bias", "outer_R20", "inner_2Rd", "claimed"]) + " |\n")
    fh.write("\n## Y* x gas scan at the final CGM, median RMS\n\n" + yg.pivot(index="ystar", columns="gas", values="median").to_markdown() + "\n")
    fh.write("\n## CGM term implications (final, last_std)\n\n")
    fh.write(f"median M_cgm(<R_last) = {dens.M_cgm_Msun.median():.2e} Msun; median mean hydrogen density inside R_last = {dens.n_mean_cm3.median():.2e} cm^-3; ")
    fh.write(f"median electron dispersion measure through R_last = {dens.DM_pc_cm3.median():.0f} pc cm^-3.\n")
print("\nwrote results/rms_tests_v1.md, results/rms_tests_v1.csv, results/rms_per_galaxy_v1.csv, results/grid_cgm_v1.csv")
