"""Circular velocities at 50-300 kpc extrapolated from the SPARC halo fits vs weak lensing.

Mistele et al. 2024 (ApJL 969, L3) deproject the KiDS-1000 excess surface density around isolated
galaxies into circular velocities and find them flat from 50 kpc to at least 300 kpc (no decline seen
to 1 Mpc), with V_flat per baryonic-mass bin on the kinematic baryonic Tully-Fisher relation. Here the
haloes fitted on the SPARC rotation curves (results/halo_fits_per_galaxy.csv: NFW and Burkert;
results/dc14_per_galaxy.csv: DC14) are extrapolated to 50, 100, 200 and 300 kpc, together with the
acceleration law (isolated galaxy, no external field), and compared bin by bin with the lensing V_flat.
Baryons are a point mass M_b = Y_fit L_36 (bulge x1.4) + 1.33 M_HI beyond 50 kpc.
Writes results/lensing_vflat.md, results/lensing_vflat.csv, results/lensing_vflat.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402
import fit_dc14 as dc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC, PC = fb.G, fb.MSUN, fb.KPC, fb.KPC / 1e3
A0 = 1.2e-10
RADII = np.array([50.0, 100.0, 200.0, 300.0])
# Mistele et al. 2024, Table 2, all isolated galaxies, 50 < R < 300 kpc
MIS = pd.DataFrame(dict(lo=[9.0, 10.5, 10.8, 11.1], hi=[10.5, 10.8, 11.1, 11.5], logMb=[10.10, 10.66, 10.96, 11.29],
                        vflat=[137.3, 182.1, 204.7, 260.2], e=[11.5, 10.8, 8.2, 6.2]))


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


E_N = 0.02  # external field in units of a0, typical large-scale-structure value for SPARC galaxies (Chae et al. 2020, 2021)


def nu_efe(z, e):
    """Chae et al. 2020 approximation of the external field effect on nu."""
    return nu(z + e) + (e / z) * (nu(z + e) - nu(e))


def m_nfw(r, rho_s, r_s):
    x = r / r_s
    return 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x / (1 + x))


def m_bur(r, rho_0, r_0):
    x = r / r_0
    return 2 * np.pi * rho_0 * r_0**3 * (np.log(1 + x) - np.arctan(x) + 0.5 * np.log(1 + x**2))


def main():
    cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")
    H = pd.read_csv(ROOT / "results/halo_fits_per_galaxy.csv").set_index("Galaxy")
    DC = pd.read_csv(ROOT / "results/dc14_per_galaxy.csv").set_index("Galaxy")
    HP = pd.read_csv(ROOT / "results/halo_priors_per_galaxy.csv").set_index("Galaxy")
    rows = []
    for g, h in H.iterrows():
        c = cat.loc[g]
        l36 = float(c.L36_1e9_Lsun) * 1e9
        lbul = float(c.Lbul_1e9_Lsun) * 1e9 if pd.notna(c.Lbul_1e9_Lsun) else 0.0
        mgas = 1.33 * float(c.MHI_1e9_Msun) * 1e9
        r = RADII * KPC
        out = dict(Galaxy=g, T=int(h["T"]), Vmax=h.Vmax, Vflat_kin=float(c.Vflat_kms) if pd.notna(c.Vflat_kms) else np.nan)
        for lab, y in (("nfw", h.y_nfw), ("nfwp", h.y_nfw), ("bur", h.y_bur), ("dc14", DC.loc[g, "y_dc14"] if g in DC.index else np.nan), ("law", h.y_law), ("lawefe", h.y_law)):
            mb = y * (l36 - lbul) + 1.4 * y * lbul + mgas
            out[f"logMb_{lab}"] = np.log10(mb)
            if lab == "nfw":
                mh = m_nfw(r, h.rho_s * MSUN / PC**3, h.r_s * KPC)
            elif lab == "nfwp":
                if g not in HP.index or not np.isfinite(HP.loc[g, "M200_nfw_p"]):
                    continue
                m200, cc = HP.loc[g, "M200_nfw_p"] * MSUN, HP.loc[g, "c_nfw_p"]
                r200 = (3 * m200 / (4 * np.pi * 200 * dc.RHO_CRIT)) ** (1 / 3)
                rs = r200 / cc
                mh = m_nfw(r, m200 / (4 * np.pi * rs**3 * (np.log(1 + cc) - cc / (1 + cc))), rs)
            elif lab == "bur":
                mh = m_bur(r, h.rho_0 * MSUN / PC**3, h.r_0 * KPC)
            elif lab == "dc14":
                if g not in DC.index or not np.isfinite(DC.loc[g, "rho_s_dc14"]):
                    continue
                d = DC.loc[g]
                F = dc.mass_profile(d.alpha, d.beta, d.gamma)
                mh = dc.m_enclosed(r, d.rho_s_dc14 * MSUN / PC**3, d.r_s_dc14 * KPC, F)
            if lab == "law":
                gb = G * mb * MSUN / r**2
                v = np.sqrt(gb * nu(gb / A0) * r)
            elif lab == "lawefe":
                gb = G * mb * MSUN / r**2
                v = np.sqrt(gb * nu_efe(gb / A0, E_N) * r)
            else:
                v = np.sqrt(G * (mh + mb * MSUN) / r)
            for R, vv in zip(RADII, v):
                out[f"V_{lab}_{R:.0f}"] = vv / 1e3
        # NFW: radius where V has dropped to 80 % of its maximum beyond 30 kpc, and R200
        rr = np.logspace(np.log10(30), 3.5, 400) * KPC
        vn = np.sqrt(G * (m_nfw(rr, h.rho_s * MSUN / PC**3, h.r_s * KPC) + (h.y_nfw * l36 + mgas) * MSUN) / rr)
        j = int(np.argmax(vn))
        k = np.where((np.arange(len(rr)) > j) & (vn < 0.8 * vn[j]))[0]
        out["R80_nfw_kpc"] = rr[k[0]] / KPC if len(k) else np.inf
        out["R200_nfw_kpc"] = (3 * h.M200 * MSUN / (4 * np.pi * 200 * dc.RHO_CRIT)) ** (1 / 3) / KPC if np.isfinite(h.M200) else np.nan
        vb = np.sqrt(G * (m_bur(rr, h.rho_0 * MSUN / PC**3, h.r_0 * KPC) + (h.y_bur * l36 + mgas) * MSUN) / rr)
        j = int(np.argmax(vb))
        k = np.where((np.arange(len(rr)) > j) & (vb < 0.8 * vb[j]))[0]
        out["R80_bur_kpc"] = rr[k[0]] / KPC if len(k) else np.inf
        rows.append(out)
    T = pd.DataFrame(rows)
    T.to_csv(ROOT / "results/lensing_vflat.csv", index=False)

    L = ["# Circular velocities at 50-300 kpc: SPARC halo fits extrapolated vs weak lensing (2026-10-11)\n",
         "Lensing: Mistele et al. 2024, KiDS-1000 isolated galaxies, V_flat over 50-300 kpc per baryonic-mass bin. "
         "SPARC galaxies are placed in the same bins with M_b from each fit (Y_fit L_36 + 1.33 M_HI); "
         "values are the median over the galaxies of the bin, with the 16-84 % range.\n"]
    L.append("| bin log M_b | lensing V_flat | N SPARC | kinematic V_flat (SPARC) | model | V(50 kpc) | V(100) | V(200) | V(300) | V(300)/V_lens | median V(300)/V_flat,kin |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    summary = {}
    for _, b in MIS.iterrows():
        for lab in ("nfw", "nfwp", "bur", "dc14", "law", "lawefe"):
            m = (T[f"logMb_{lab}"] >= b.lo) & (T[f"logMb_{lab}"] < b.hi) & np.isfinite(T[f"V_{lab}_300"])
            s = T[m]
            if len(s) < 3:
                continue
            vals = [s[f"V_{lab}_{R:.0f}"] for R in RADII]
            cell = " | ".join(f"{v.median():.0f} ({v.quantile(.16):.0f}-{v.quantile(.84):.0f})" for v in vals)
            vk = s.Vflat_kin.dropna()
            ratio = (s[f"V_{lab}_300"] / s.Vflat_kin).dropna()
            L.append(f"| {b.lo}-{b.hi} ({b.logMb}) | {b.vflat} +/- {b.e} | {len(s)} | {vk.median():.0f} ({vk.quantile(.16):.0f}-{vk.quantile(.84):.0f}) | {lab} | {cell} | "
                     f"{vals[-1].median() / b.vflat:.2f} | {ratio.median():.2f} ({ratio.quantile(.16):.2f}-{ratio.quantile(.84):.2f}) |")
            summary[(b.logMb, lab)] = (len(s), [v.median() for v in vals], [v.quantile(.16) for v in vals], [v.quantile(.84) for v in vals])
    L.append("")
    L.append("Shape test over all SPARC galaxies with a kinematic V_flat (lensing: V(50-300 kpc) = V_flat,kin within the errors, i.e. ratio 1):")
    L.append("| model | N | median V(50)/V_flat,kin | V(100)/V_flat,kin | V(200)/V_flat,kin | V(300)/V_flat,kin | 16-84 % at 300 kpc |")
    L.append("|---|---|---|---|---|---|---|")
    for lab, name in (("nfw", "NFW, free fit"), ("nfwp", "NFW, LambdaCDM priors c(M) + SHM"), ("bur", "Burkert"), ("dc14", "DC14, LambdaCDM priors"), ("law", "acceleration law, isolated"), ("lawefe", f"acceleration law + external field e_N = {E_N} a0")):
        ok = np.isfinite(T.Vflat_kin) & np.isfinite(T[f"V_{lab}_300"])
        rat = [(T[f"V_{lab}_{R:.0f}"] / T.Vflat_kin)[ok] for R in RADII]
        L.append(f"| {name} | {ok.sum()} | " + " | ".join(f"{x.median():.2f}" for x in rat) + f" | {rat[-1].quantile(.16):.2f}-{rat[-1].quantile(.84):.2f} |")
    ok = np.isfinite(T.Vflat_kin) & np.isfinite(T.logMb_law) & (T.logMb_law >= 10.0)   # the mass range of the lensing stacks
    mb = 10 ** T.logMb_law.values[ok] * MSUN
    r300 = 300 * KPC
    scan = []
    for e in np.logspace(-6, -1, 51):
        gb = G * mb / r300**2
        v = np.sqrt(gb * nu_efe(gb / A0, e) * r300) / 1e3
        scan.append((e, np.median(v / T.Vflat_kin.values[ok])))
    scan = np.array(scan)
    r0 = scan[0, 1]                                   # isolated limit
    e90 = np.interp(0.9 * r0, scan[::-1, 1], scan[::-1, 0])
    L.append(f"External field scan ({ok.sum()} galaxies with log M_b >= 10): median V(300 kpc)/V_flat,kin = " + ", ".join(f"{r:.2f} at e_N = {e:.0e}" for e, r in scan[::10]) +
             f"; the ratio stays within 10 % of its isolated value only for e_N < {e90:.1e} a0 (Chae et al. 2020 formula; the Milky Way sits in e_N = 0.011-0.012 a0, Chae et al. 2021).")
    L.append("")
    fin = np.isfinite(T.R80_nfw_kpc) & (T.R80_nfw_kpc < np.inf)
    L.append(f"NFW fits: median R200 = {T.R200_nfw_kpc.median():.0f} kpc; radius where V falls to 80 % of its peak: median "
             f"{T.R80_nfw_kpc[fin].median():.0f} kpc ({T.R80_nfw_kpc[fin].quantile(.16):.0f}-{T.R80_nfw_kpc[fin].quantile(.84):.0f}), "
             f"{(T.R80_nfw_kpc > 300).mean():.0%} of the fits stay above 80 % out to 300 kpc.")
    fin = np.isfinite(T.R80_bur_kpc) & (T.R80_bur_kpc < np.inf)
    L.append(f"Burkert fits: radius where V falls to 80 % of its peak: median {T.R80_bur_kpc[fin].median():.0f} kpc "
             f"({T.R80_bur_kpc[fin].quantile(.16):.0f}-{T.R80_bur_kpc[fin].quantile(.84):.0f}), {(T.R80_bur_kpc > 300).mean():.0%} stay above 80 % to 300 kpc.")
    L.append(f"Acceleration law: V is flat by construction once g_bar << a0 for an isolated galaxy; with an external field e_N = {E_N} a0 "
             "(Chae et al. 2020 formula) the law becomes Newtonian with a boosted G once g_bar << e_N a0 and V declines as R^-1/2 too.")
    with open(ROOT / "results/lensing_vflat.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n\nFigure: results/lensing_vflat.png\n")
    print("\n".join(L))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bins = [b for b in MIS.logMb if any(k[0] == b for k in summary)]
    fig, ax = plt.subplots(1, len(bins), figsize=(4.2 * len(bins), 4.2), squeeze=False)
    cols = dict(nfw="tab:blue", nfwp="tab:cyan", bur="tab:green", dc14="tab:purple", law="tab:orange", lawefe="tab:red")
    for a, b in zip(ax[0], bins):
        row = MIS[MIS.logMb == b].iloc[0]
        a.axhspan(row.vflat - row.e, row.vflat + row.e, color="gray", alpha=0.3, label="lensing V_flat, 50-300 kpc")
        for lab in ("nfw", "nfwp", "bur", "dc14", "law", "lawefe"):
            if (b, lab) in summary:
                n, med, lo, hi = summary[(b, lab)]
                a.plot(RADII, med, "o-", color=cols[lab], label=f"{lab} (N={n})")
                a.fill_between(RADII, lo, hi, color=cols[lab], alpha=0.12)
        a.set_title(f"log M_b = {row.lo}-{row.hi}")
        a.set_xlabel("R (kpc)")
        a.set_ylabel("V_circ (km/s)")
        a.set_xscale("log")
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(ROOT / "results/lensing_vflat.png", dpi=130)


if __name__ == "__main__":
    main()
