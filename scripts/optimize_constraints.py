"""Constrained optimisation of the baryonic + CGM model family of tag 1.0.

Finds, level by level, the parameters that minimise the median RMS while
satisfying an increasing set of physical criteria. Family:

    V^2 = g * Vgas|Vgas| + Y* (Vdisk|Vdisk| + Vbul|Vbul|) + G M_cgm(<R) / R
    M_cgm(<R) = M * M_std * f_beta(R / r_c) / f_beta(R_ref / r_c)
    f_beta(x) = int_0^x t^2 (1 + t^2)^(-3 beta / 2) dt   (beta profile)
    M_std = 0.5 L36 + 1.33 MHI   (standard SPARC baryons, 1e9 Msun)

Normalisation modes for R_ref:  "last" = R_last of each galaxy (hidden parameter),
"fixed" = R0 kpc for all galaxies, "rd" = k * R_disk (self-similar halo, r_c = k2 * R_disk).

Levels:
  L0  minimise RMS only, mode last
  L1  + |bias| <= 1 and |outer residual (R > 20 kpc)| <= 2 km/s
  L2  + no hidden parameter: mode fixed or rd
  L3  + Y* <= 0.7 and gas factor <= 1.3 (SPS and He + H2 budget; SPARC Vgas already includes He)
  L4  + hot-gas density: mean n_H inside R_last <= 3 x Milky-Way beta-model value (Miller & Bregman 2015)
  L5  reference: no halo (M = 0), Y* and gas free within L3 bounds

Writes results/optimize_constraints.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
MP = 1.6726e-27
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")

gal = []
for f in sorted((ROOT / "data/rotmod").glob("*_rotmod.dat")):
    df = fb.load_rotmod(f)
    if len(df) < 5:
        continue
    n = f.stem.replace("_rotmod", "")
    c = cat.loc[n]
    gas, d, b = df.Vgas.values * 1e3, df.Vdisk.values * 1e3, df.Vbul.values * 1e3
    r = df.R.values
    gal.append(dict(name=n, r=r, R=r * KPC, gas2=gas * np.abs(gas), star2=d * np.abs(d) + b * np.abs(b),
                    vobs=df.Vobs.values, T=int(c["T"]), rd=float(c.Rdisk_kpc),
                    mstd=(0.5 * float(c.sparc_L36_1e9_Lsun) + 1.33 * float(c.sparc_MHI_1e9_Msun)) * 1e9 * MSUN,
                    outer=r > 20.0, inner=r < 2.0 * float(c.Rdisk_kpc)))
EARLY = np.array([g["T"] <= 4 for g in gal])
RLAST = np.array([g["r"][-1] for g in gal])
# physical bound on the mean hydrogen density inside R (Milky-Way hot halo, n ~ 1.35e-2 r_kpc^-1.5 cm^-3), x3
N_MAX = 3.0 * 2.7e-2 * RLAST ** (-1.5)

_fcache = {}


def f_beta(beta):
    key = round(beta, 3)
    if key not in _fcache:
        x = np.concatenate([[0.0], np.logspace(-3, 4, 1400)])
        y = x**2 * (1 + x**2) ** (-1.5 * key)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])
        _fcache[key] = (x, cum)
    return _fcache[key]


def evaluate(ystar, gfac, M, beta, mode, rc=None, R0=None, k=None, k2=None):
    xg, cum = f_beta(beta)
    rms, bias, outer, inner, nh = [], [], [], [], []
    for g in gal:
        vbar2 = gfac * g["gas2"] + ystar * g["star2"]
        if M > 0:
            r_c = rc if mode != "rd" else k2 * g["rd"]
            r_ref = g["r"][-1] if mode == "last" else (R0 if mode == "fixed" else k * g["rd"])
            fr = np.interp(g["r"] / r_c, xg, cum)
            fref = np.interp(r_ref / r_c, xg, cum)
            mc = M * g["mstd"] * fr / max(fref, 1e-30)
        else:
            mc = np.zeros_like(g["R"])
        v = np.sqrt(np.maximum(vbar2 + G * mc / g["R"], 0.0)) / 1e3
        res = g["vobs"] - v
        rms.append(np.sqrt(np.mean(res**2)))
        bias.append(res.mean())
        outer.append(res[g["outer"]].mean() if g["outer"].any() else np.nan)
        inner.append(res[g["inner"]].mean() if g["inner"].any() else np.nan)
        rl = g["R"][-1]
        nh.append(mc[-1] / (4 / 3 * np.pi * rl**3) / (1.4 * MP) * 1e-6)
    rms, bias, outer, inner, nh = map(np.array, (rms, bias, outer, inner, nh))
    return dict(rms=np.median(rms), bias=np.median(bias), outer=np.nanmedian(outer), inner=np.nanmedian(inner),
                early_rms=np.median(rms[EARLY]), early_inner=np.nanmedian(inner[EARLY]), late_rms=np.median(rms[~EARLY]),
                nh=np.median(nh), nh_viol=np.mean(nh > N_MAX), nh_max_ratio=np.max(nh / N_MAX))


LEVELS = {
    "L0 RMS only, mode last": dict(mode="last", bias=False, phys=False, dens=False),
    "L1 + bias~0, outer~0, mode last": dict(mode="last", bias=True, phys=False, dens=False),
    "L2a + no hidden parameter, mode fixed R0": dict(mode="fixed", bias=True, phys=False, dens=False),
    "L2b + no hidden parameter, self-similar (R_ref, r_c prop. to Rd)": dict(mode="rd", bias=True, phys=False, dens=False),
    "L3a + Y*<=0.7, gas<=1.3, mode fixed": dict(mode="fixed", bias=True, phys=True, dens=False),
    "L3b + Y*<=0.7, gas<=1.3, self-similar": dict(mode="rd", bias=True, phys=True, dens=False),
    "L4a + hot-gas density bound, mode fixed": dict(mode="fixed", bias=True, phys=True, dens=True),
    "L4b + hot-gas density bound, self-similar": dict(mode="rd", bias=True, phys=True, dens=True),
    "L4c + hot-gas density bound, mode last (hidden parameter allowed)": dict(mode="last", bias=True, phys=True, dens=True),
    "L5 no halo (M=0), Y*<=0.7, gas<=1.3": dict(mode="none", bias=False, phys=True, dens=False),
}


def bounds_for(cfg):
    ymax, gmax = (0.7, 1.3) if cfg["phys"] else (1.5, 3.0)
    b = [(0.3, ymax), (1.0, gmax)]
    if cfg["mode"] == "none":
        return b
    b += [(0.0, 12.0), (0.4, 1.5)]
    if cfg["mode"] == "last":
        b += [(1.0, 60.0)]                      # rc
    elif cfg["mode"] == "fixed":
        b += [(1.0, 60.0), (5.0, 150.0)]        # rc, R0
    else:
        b += [(2.0, 40.0), (0.2, 15.0)]         # k, k2
    return b


def unpack(p, cfg):
    kw = dict(ystar=p[0], gfac=p[1])
    if cfg["mode"] == "none":
        return dict(kw, M=0.0, beta=0.7, mode="none")
    kw.update(M=p[2], beta=p[3], mode=cfg["mode"])
    if cfg["mode"] == "last":
        kw["rc"] = p[4]
    elif cfg["mode"] == "fixed":
        kw["rc"], kw["R0"] = p[4], p[5]
    else:
        kw["k"], kw["k2"] = p[4], p[5]
    return kw


def objective(p, cfg, bnds):
    for v, (lo, hi) in zip(p, bnds):
        if v < lo or v > hi:
            return 1e3 + 10 * sum(max(0, lo - v, v - hi) for v, (lo, hi) in zip(p, bnds))
    m = evaluate(**unpack(p, cfg))
    val = m["rms"]
    if cfg["bias"]:
        val += 3 * max(0, abs(m["bias"]) - 1.0) + 3 * max(0, abs(m["outer"]) - 2.0)
    if cfg["dens"]:
        val += 30 * m["nh_viol"] + 3 * max(0, np.log10(m["nh_max_ratio"]))
    return val


def search(cfg, n_random=600, seed=0):
    rng = np.random.default_rng(seed)
    bnds = bounds_for(cfg)
    samples = [np.array([rng.uniform(lo, hi) for lo, hi in bnds]) for _ in range(n_random)]
    scored = sorted(((objective(p, cfg, bnds), p) for p in samples), key=lambda t: t[0])
    best = None
    for _, p0 in scored[:3]:
        res = minimize(objective, p0, args=(cfg, bnds), method="Nelder-Mead",
                       options=dict(maxiter=600, xatol=1e-3, fatol=1e-3))
        if best is None or res.fun < best.fun:
            best = res
    return best, bnds


rows = []
for name, cfg in LEVELS.items():
    best, bnds = search(cfg)
    kw = unpack(best.x, cfg)
    m = evaluate(**kw)
    par = ", ".join(f"{k}={v:.2f}" for k, v in kw.items() if isinstance(v, float))
    row = dict(level=name, params=par, score=round(best.fun, 2), **{k: round(float(v), 2) for k, v in m.items()})
    rows.append(row)
    print(f"{name}\n   {par}\n   RMS={m['rms']:.1f} bias={m['bias']:+.1f} outer={m['outer']:+.1f} inner={m['inner']:+.1f} "
          f"early RMS={m['early_rms']:.1f} early inner={m['early_inner']:+.1f} late RMS={m['late_rms']:.1f} "
          f"n_H med={m['nh']:.1e} cm^-3  frac>bound={m['nh_viol']:.2f}  max n/bound={m['nh_max_ratio']:.1e}")

T = pd.DataFrame(rows)
T.to_csv(ROOT / "results/optimize_constraints.csv", index=False)
with open(ROOT / "results/optimize_constraints.md", "w", encoding="utf-8") as fh:
    fh.write("# Constrained optimisation of the tag-1.0 model family (2026-10-10)\n\n")
    fh.write("Median over 171 SPARC galaxies. RMS, bias, outer (R>20 kpc), inner (R<2 Rd) in km/s. n_H = median mean hydrogen density of the halo inside R_last; ")
    fh.write("bound = 3 x Milky-Way hot-halo beta model. frac>bound = fraction of galaxies violating it.\n\n")
    fh.write("| level | parameters | RMS | bias | outer | inner | early RMS | early inner | late RMS | n_H (cm^-3) | frac>bound | max n/bound |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        fh.write(f"| {r['level']} | {r['params']} | {r['rms']} | {r['bias']} | {r['outer']} | {r['inner']} | {r['early_rms']} | {r['early_inner']} | {r['late_rms']} | {r['nh']:.1e} | {r['nh_viol']} | {r['nh_max_ratio']:.1e} |\n")
print("\nwrote results/optimize_constraints.md")
