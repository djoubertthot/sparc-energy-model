"""RMS tests of the final documented model (scripts/final_model.py) on the SPARC curves.

Runs every test the repository announces: budget vs calibrated normalisation,
beta x rho0 grid, accumulation factors (type, linear age, type + Rd), dependencies
of rho0 on galaxy properties, quality and inclination subsamples, weighted RMS,
100-galaxy subset, kinetic term, Monte-Carlo on inclination. Two reference rows
(baryons only, MOND simple mu) are appended for scale.

Outputs (results/): rms_tests.md, rms_tests.csv, grid_search_beta_rho.csv, rms_per_galaxy.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_model as fm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")

galaxies = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fm.load_rotmod(f)
    if len(df) < 5:
        continue
    name = f.stem.replace("_rotmod", "")
    c = cat.loc[name]
    vmax = c.Vflat_kms if pd.notna(c.Vflat_kms) and c.Vflat_kms > 0 else df.Vobs.max()
    mbar = fm.YSTAR * c.L36_1e9_Lsun + 1.33 * c.MHI_1e9_Msun
    galaxies.append(dict(name=name, df=df, T=int(c["T"]), Q=int(c.Q), inc=float(c.Inc_deg), e_inc=float(c.e_Inc_deg),
                         vmax=float(vmax), mbar=float(mbar), sb=float(c.SBdisk_Lsun_pc2), rd=float(c.Rdisk_kpc),
                         fgas=float(1.33 * c.MHI_1e9_Msun / mbar)))
N_ALL = len(galaxies)
MED_RD = float(np.median([g["rd"] for g in galaxies]))


def f_type(T, inverse=False):
    f = 1.2 if T <= 4 else (0.8 if T >= 8 else 1.0)
    return 1 / f if inverse else f


def f_linear(T):  # t_stellar / 10 Gyr: early ~11, intermediate ~9, late ~6.5 Gyr (docs/PHYSICAL_AGE_FACTOR.md)
    return 1.1 if T <= 4 else (0.65 if T >= 8 else 0.9)


def final(g, norm="calibrated", rho0=1e-24, beta=fm.BETA, f_acc=1.0, vobs=None, kinetic=False, **kw):
    out = fm.compute_final_model(g["df"], rdisk_kpc=g["rd"], f_acc=f_acc, beta=beta, norm=norm, rho0=rho0, **kw)
    if vobs is not None:
        out["residual"] = vobs - out["Vmodel"]
    if kinetic:  # m_k = (1/2) M_bar(<R) V_bar^2 / c^2 added as extra mass (EXTENSIONS.md)
        R = g["df"].R.values * fm.KPC
        vb2 = (out["Vbar"].values * 1e3) ** 2
        mk = 0.5 * (vb2 * R / fm.G) * vb2 / fm.C**2
        vm = np.sqrt((out["Vmodel"].values * 1e3) ** 2 + fm.G * mk / R) / 1e3
        out["residual"] = g["df"].Vobs.values - vm
    return fm.rms(out), fm.rms(out, weighted=True)


def baryons_only(g):
    df = g["df"]
    vb = np.sqrt(np.maximum(fm.baryonic(df), 0)) / 1e3
    res = df.Vobs.values - vb
    w = 1 / np.maximum(df.eVobs.values, 0.1) ** 2
    return float(np.sqrt(np.mean(res**2))), float(np.sqrt(np.average(res**2, weights=w)))


def mond(g):
    df = g["df"]
    R = df.R.values * fm.KPC
    gn = np.maximum(fm.baryonic(df), 1.0) / R
    gm = gn / 2 + np.sqrt(gn**2 / 4 + gn * fm.A0)
    res = df.Vobs.values - np.sqrt(gm * R) / 1e3
    w = 1 / np.maximum(df.eVobs.values, 0.1) ** 2
    return float(np.sqrt(np.mean(res**2))), float(np.sqrt(np.average(res**2, weights=w)))


def run(fn, subset=None, **kw):
    a = np.array([fn(g, **kw) for g in galaxies if subset is None or subset(g)])
    return dict(N=int(len(a)), median=float(np.median(a[:, 0])), mean=float(np.mean(a[:, 0])),
                p16=float(np.percentile(a[:, 0], 16)), p84=float(np.percentile(a[:, 0], 84)),
                wmedian=float(np.median(a[:, 1])))


rows = []


def add(test, group, s, note=""):
    rows.append(dict(group=group, test=test, N=s["N"], median=round(s["median"], 2), wmedian=round(s["wmedian"], 2),
                     mean=round(s["mean"], 1), p16=round(s["p16"], 1), p84=round(s["p84"], 1), note=note))
    print(f"{group:13s} {test:52s} N={s['N']:3d} median={s['median']:8.2f} wmedian={s['wmedian']:8.2f} "
          f"mean={s['mean']:9.1f} [{s['p16']:.1f}, {s['p84']:.1f}]")


# 1. normalisations ---------------------------------------------------------
add("budget: literal inventory, f_missing=80, no free parameter", "final", run(final, norm="budget"), "rho_uniform=1.6e-22")
add("budget, without internal stellar part", "final", run(final, norm="budget", include_internal=False))
add("calibrated: uniform part = rho0 = 1e-24 (archive_v1/THEORY_v1.md)", "final", run(final, rho0=1e-24), "hybrid, f_missing matched")
add("calibrated: uniform part = rho0 = 2e-24 (MEMORY.md)", "final", run(final, rho0=2e-24))
add("calibrated 1e-24, gamma = 1 exactly", "final", run(final, rho0=1e-24, gamma=False))
add("calibrated 1e-24 + kinetic term m_k = E_k/c^2", "final", run(final, rho0=1e-24, kinetic=True))

# 2. grid beta x rho0 (calibrated) -----------------------------------------
grid = []
for beta in [1.0, 1.2, 1.4, 1.6, 1.8]:
    for rho0 in [1e-25, 5e-25, 1e-24, 2e-24, 5e-24]:
        s = run(final, rho0=rho0, beta=beta)
        grid.append(dict(beta=beta, rho0=rho0, median_rms=round(s["median"], 2), wmedian=round(s["wmedian"], 2), mean=round(s["mean"], 1)))
grid = pd.DataFrame(grid)
grid.to_csv(ROOT / "results/grid_search_beta_rho.csv", index=False)
best = grid.sort_values("median_rms").iloc[0]
print("\n--- grid beta x rho0, median RMS (calibrated normalisation) ---")
print(grid.pivot(index="beta", columns="rho0", values="median_rms").to_string())
add(f"grid best: beta={best.beta}, rho0={best.rho0:.0e}", "grid", run(final, rho0=best.rho0, beta=best.beta), "25-point grid")
B_RHO, B_BETA = float(best.rho0), float(best.beta)

# 3. accumulation factor -----------------------------------------------------
add("f_acc = 1 (no accumulation proxy)", "accumulation", run(final, rho0=B_RHO, beta=B_BETA))
add("binary type proxy (early 1.2 / late 0.8)", "accumulation", run(lambda g: final(g, rho0=B_RHO, beta=B_BETA, f_acc=f_type(g["T"]))))
add("inverse binary (early 0.8 / late 1.2)", "accumulation", run(lambda g: final(g, rho0=B_RHO, beta=B_BETA, f_acc=f_type(g["T"], True))))
add("linear age t/10 Gyr (1.1 / 0.9 / 0.65)", "accumulation", run(lambda g: final(g, rho0=B_RHO, beta=B_BETA, f_acc=f_linear(g["T"]))))
best_a, best_s = None, None
for alpha in [-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0]:
    s = run(lambda g, a=alpha: final(g, rho0=B_RHO, beta=B_BETA, f_acc=f_type(g["T"]) * (g["rd"] / MED_RD) ** a))
    print(f"  type + Rd  alpha={alpha:+5.2f}  median={s['median']:6.2f}")
    if best_s is None or s["median"] < best_s["median"]:
        best_a, best_s = alpha, s
add(f"type + (Rd/median)^alpha, best alpha={best_a:+.2f}", "accumulation", best_s, "docs claim 43.8 km/s")
add("beta early/late = 1.6 / 1.2", "accumulation", run(lambda g: final(g, rho0=B_RHO, beta=1.6 if g["T"] <= 4 else (1.2 if g["T"] >= 8 else B_BETA))))
add("beta early/late = 1.2 / 1.6", "accumulation", run(lambda g: final(g, rho0=B_RHO, beta=1.2 if g["T"] <= 4 else (1.6 if g["T"] >= 8 else B_BETA))))

# 4. dependencies rho0_i = rho0 (X_i / median X)^alpha ------------------------
print("\n--- dependencies ---")
for key, label in [("mbar", "M_bary"), ("vmax", "Vmax"), ("sb", "surface density"), ("fgas", "gas fraction"), ("rd", "R_disk")]:
    med = float(np.median([g[key] for g in galaxies]))
    best_a, best_s = None, None
    for alpha in [-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0]:
        s = run(lambda g, a=alpha: final(g, rho0=B_RHO * (g[key] / med) ** a, beta=B_BETA))
        print(f"  {label:16s} alpha={alpha:+5.2f}  median={s['median']:6.2f}")
        if best_s is None or s["median"] < best_s["median"]:
            best_a, best_s = alpha, s
    add(f"rho0 ~ {label}^alpha, best alpha={best_a:+.2f}", "dependencies", best_s, "7-value scan")

# 5. subsamples ----------------------------------------------------------------
sub = dict(rho0=B_RHO, beta=B_BETA)
add("Q <= 2", "subsamples", run(final, subset=lambda g: g["Q"] <= 2, **sub))
add("40 < i < 75", "subsamples", run(final, subset=lambda g: 40 < g["inc"] < 75, **sub))
add("Q <= 2 and 40 < i < 75", "subsamples", run(final, subset=lambda g: g["Q"] <= 2 and 40 < g["inc"] < 75, **sub))
add("Q = 1 and 40 < i < 75", "subsamples", run(final, subset=lambda g: g["Q"] == 1 and 40 < g["inc"] < 75, **sub))
add("i > 75 only", "subsamples", run(final, subset=lambda g: g["inc"] > 75, **sub))
add("first 100 rotmod files", "subsamples", run(final, subset=lambda g: galaxies.index(g) < 100, **sub))

# 6. Monte-Carlo on inclination -----------------------------------------------
rng = np.random.default_rng(1)
meds = []
for _ in range(100):
    vals = []
    for g in galaxies:
        i0 = np.radians(g["inc"])
        i1 = np.radians(np.clip(g["inc"] + rng.normal(0, g["e_inc"]), 15, 89))
        vals.append(final(g, vobs=g["df"].Vobs.values * np.sin(i0) / np.sin(i1), **sub)[0])
    meds.append(np.median(vals))
add("MC inclination +/- e_Inc, 100 draws (p16/p84 over draws)", "uncertainty",
    dict(N=N_ALL, median=float(np.median(meds)), mean=float(np.mean(meds)), p16=float(np.percentile(meds, 16)),
         p84=float(np.percentile(meds, 84)), wmedian=float("nan")))

# 7. references ----------------------------------------------------------------
add("baryons only, Y* = 0.5, no dark component", "reference", run(baryons_only))
add("MOND simple mu, a0 = 1.2e-10, Y* = 0.5, 0 free parameters", "reference", run(mond))

# outputs --------------------------------------------------------------------
per = pd.DataFrame([dict(Galaxy=g["name"], N=len(g["df"]), T=g["T"], Q=g["Q"], Inc=g["inc"], Rd=g["rd"],
                         RMS_budget=round(final(g, norm="budget")[0], 2),
                         RMS_calibrated=round(final(g, rho0=B_RHO, beta=B_BETA)[0], 2),
                         WRMS_calibrated=round(final(g, rho0=B_RHO, beta=B_BETA)[1], 2),
                         RMS_baryons_only=round(baryons_only(g)[0], 2), RMS_MOND=round(mond(g)[0], 2))
                    for g in galaxies])
per.to_csv(ROOT / "results/rms_per_galaxy.csv", index=False)
pd.DataFrame(rows).to_csv(ROOT / "results/rms_tests.csv", index=False)
with open(ROOT / "results/rms_tests.md", "w", encoding="utf-8") as fh:
    fh.write("# RMS tests of the final documented model (2026-10-10)\n\n")
    fh.write("Model: scripts/final_model.py, the formula of archive_v1/THEORY_v1.md with the documented inventory. ")
    fh.write(f"Sample: {N_ALL} SPARC galaxies with at least 5 points (D512-2, NGC6789, UGC00634, UGC07232 dropped). ")
    fh.write("RMS in km/s, unweighted unless noted; wmedian = median of the 1/eVobs^2-weighted RMS; p16/p84 over galaxies.\n\n")
    fh.write("| group | test | N | median | wmedian | mean | p16 | p84 | note |\n|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        fh.write(f"| {r['group']} | {r['test']} | {r['N']} | {r['median']} | {r['wmedian']} | {r['mean']} | {r['p16']} | {r['p84']} | {r['note']} |\n")
    fh.write("\n## Grid beta x rho0, median RMS, calibrated normalisation\n\n")
    fh.write(grid.pivot(index="beta", columns="rho0", values="median_rms").to_markdown() + "\n")
print("\nwrote results/rms_tests.md, results/rms_tests.csv, results/grid_search_beta_rho.csv, results/rms_per_galaxy.csv")
