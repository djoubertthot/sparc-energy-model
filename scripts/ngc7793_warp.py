"""NGC 7793: tilted rings with the inclination free per ring on the MHONGOOSE velocity fields.

The 2D fit of scripts/mhongoose_rings.py held the inclination fixed at the value fitted inside R_last
(51 deg) and found the curve declining to V_bar beyond 8 kpc, residual -0.4 dex from the acceleration law.
Gentile et al. 2011 note that the THINGS inclination of this galaxy (de Blok et al. 2008: 49.6 deg on
average, PA 290) "is low and presents large variations in adjacent radii". Here each ring (2 beams wide,
1 beam step) is fitted with V_rot, V_rad, PA and inclination free (V_sys and centre fixed from the inner
geometry), on the 26" map inside 10 kpc and on the 60" map outside. The chi2 profile in inclination gives
its 1 sigma range with the number of independent beams in the ring as the number of degrees of freedom.
For each ring the inclination that would put the ring on the acceleration law is also given.
Writes results/ngc7793_warp.md, results/ngc7793_warp.csv, results/ngc7793_warp.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402
import mhongoose_rings as mh  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, KPC, A0 = fb.G, fb.KPC, 1.2e-10
SP, GAL = "NGC7793", "J2357-32"
c = mh.cat.loc[SP]
D = float(c.D_Mpc)
KAS = D * 1e3 * np.pi / 180 / 3600          # kpc per arcsec
df = fb.load_rotmod(ROOT / f"data/rotmod/{SP}_rotmod.dat")
r_last = df.R.values[-1]
prev = pd.read_csv(ROOT / f"results/mhongoose_{SP}.csv").sort_values("R_kpc")
vbar_of = lambda rk: np.interp(rk, prev.R_kpc.values, prev.Vbar.values)  # noqa: E731


def nu(y):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))


def geometry(v, ra, dec):
    """Same geometry stage as mhongoose_rings.main: PA scan, bounded refinement inside R_last."""
    ra0, dec0 = float(c.RA_deg), float(c.Dec_deg)
    best = None
    for pa_try in (0, 60, 120, 180, 240, 300):
        g = mh.fit_geometry(ra, dec, v, 0.5 * r_last / KAS, r_last / KAS, ra0, dec0, np.nanmedian(v), pa_try, float(c.Inc_deg))
        r, th = mh.ring_coords(ra, dec, ra0 + g[0] / 3600 / np.cos(np.radians(dec0)), dec0 + g[1] / 3600, g[3], g[4])
        m = np.isfinite(v) & (r < r_last / KAS) & (np.abs(np.cos(th)) > 0.3)
        chi = np.sum((v[m] - g[2] - g[5] * np.cos(th[m]) * np.sin(np.radians(g[4]))) ** 2)
        if best is None or chi < best[0]:
            best = (chi, g)
    g = best[1]
    if g[5] < 0:
        g[3] = (g[3] + 180) % 360
        g[5] = -g[5]
    return ra0 + g[0] / 3600 / np.cos(np.radians(dec0)), dec0 + g[1] / 3600, g[2], g[3], g[4]


def ring_scan(ra, dec, v, ra_c, dec_c, vsys, radii, width, beam, pix, inc_grid, pa_grid, min_pix=40):
    ok = np.isfinite(v)
    ra, dec, v = ra[ok], dec[ok], v[ok]
    beams_per_pix = pix**2 / (1.133 * beam**2)
    rows = {R: [] for R in radii}
    for pa in pa_grid:
        for inc in inc_grid:
            r, th = mh.ring_coords(ra, dec, ra_c, dec_c, pa, inc)
            si = np.sin(np.radians(inc))
            for R in radii:
                m = np.abs(r - R) < width / 2
                n = int(m.sum())
                if n < min_pix:
                    continue
                cover = len(np.unique((np.degrees(th[m]) % 360 // 15).astype(int))) / 24.0
                if cover < 0.5:
                    continue
                A = np.vstack([np.cos(th[m]) * si, np.sin(th[m]) * si]).T
                coef, *_ = np.linalg.lstsq(A, v[m] - vsys, rcond=None)
                if coef[0] <= 0:
                    continue
                res = v[m] - vsys - A @ coef
                rows[R].append((np.mean(res**2), pa, inc, coef[0], coef[1], n, cover, np.sum(A[:, 0] ** 2)))
    out = []
    for R in radii:
        if not rows[R]:
            continue
        t = np.array(rows[R])
        j = int(np.argmin(t[:, 0]))
        chi2n, pa, inc, vrot, vrad, n, cover, sc = t[j]
        nb = max(n * beams_per_pix, 2.0)
        prof = {}
        for row in t:
            prof[row[2]] = min(prof.get(row[2], np.inf), row[0])
        incs = np.array(sorted(prof))
        chis = np.array([prof[i] for i in incs])
        acc = chis <= chi2n * (1 + 1 / nb)            # delta chi2 = 1 with errors rescaled to chi2/N_beam = 1
        inc_lo, inc_hi = incs[acc].min(), incs[acc].max()
        e_stat = np.sqrt(chi2n / sc * (n / nb))
        v_lo = vrot * np.sin(np.radians(inc)) / np.sin(np.radians(inc_hi))
        v_hi = vrot * np.sin(np.radians(inc)) / np.sin(np.radians(inc_lo))
        out.append(dict(R_as=R, R_kpc=R * KAS, Vrot=vrot / 1e3, e_stat=e_stat / 1e3, V_lo=v_lo / 1e3, V_hi=v_hi / 1e3,
                        Vrad=vrad / 1e3, PA=pa, inc=inc, inc_lo=inc_lo, inc_hi=inc_hi, npix=int(n), nbeam=nb, cover=cover,
                        rms=np.sqrt(chi2n) / 1e3))
    return out


def main():
    h1, v1 = mh.load(GAL, "r10_t00", "mom1")
    ra1, dec1 = mh.sky_grid(h1, v1)
    ra_c, dec_c, vsys, pa_g, inc_g = geometry(v1, ra1, dec1)
    beam1, pix1 = h1["BMAJ"] * 3600, abs(h1["CDELT2"]) * 3600
    h6, v6 = mh.load(GAL, "r05_t60", "mom1")
    ra6, dec6 = mh.sky_grid(h6, v6)
    beam6, pix6 = h6["BMAJ"] * 3600, abs(h6["CDELT2"]) * 3600
    inc_grid = np.arange(20, 81, 2.0)
    pa_grid = np.arange(pa_g - 90, pa_g + 91, 4.0)       # the outer rings twist by more than 40 deg
    radii1 = np.arange(3.0 * beam1, 10.0 / KAS, beam1)   # rings closer than 3 beams are not resolved
    radii6 = np.arange(1.5 * beam6, 2.2 * r_last / KAS, beam6)
    r1 = ring_scan(ra1, dec1, v1, ra_c, dec_c, vsys, radii1, 2 * beam1, beam1, pix1, inc_grid, pa_grid, 60)
    r6 = ring_scan(ra6, dec6, v6, ra_c, dec_c, vsys, radii6, 2 * beam6, beam6, pix6, inc_grid, pa_grid, 30)
    T = pd.DataFrame([dict(src="26arcsec", **r) for r in r1] + [dict(src="60arcsec", **r) for r in r6])
    T["Vbar"] = vbar_of(T.R_kpc.values)
    g_bar = (T.Vbar.values * 1e3) ** 2 / (T.R_kpc.values * KPC)
    v_rar = np.sqrt(g_bar * nu(g_bar / A0) * T.R_kpc.values * KPC) / 1e3
    T["V_RAR"] = v_rar
    T["res_RAR"] = np.log10((T.Vrot / v_rar) ** 2)
    T["res_RAR_lo"] = np.log10((T.V_lo / v_rar) ** 2)
    T["res_RAR_hi"] = np.log10((T.V_hi / v_rar) ** 2)
    s = T.Vrot.values * np.sin(np.radians(T.inc.values)) / v_rar
    T["inc_RAR"] = np.where(s < 1, np.degrees(np.arcsin(np.clip(s, 0, 1))), np.nan)
    T["beyond_sparc"] = T.R_kpc > r_last
    T.to_csv(ROOT / "results/ngc7793_warp.csv", index=False)

    q = '"'
    L = ["# NGC 7793: inclination free per ring on the MHONGOOSE velocity fields (2026-10-11)\n",
         f"Geometry inside R_last = {r_last:.1f} kpc (26{q} map): Vsys = {vsys / 1e3:.1f} km/s, PA = {pa_g:.0f} deg, i = {inc_g:.0f} deg "
         f"(SPARC {c.Inc_deg} +/- {c.e_Inc_deg}; THINGS, de Blok et al. 2008: 49.6 deg, PA 290.1, D = 3.9 Mpc). SPARC D = {D} Mpc.",
         "Rings 2 beams wide, 1 beam step; V_rot, V_rad, PA and i free per ring (PA within +/- 90 deg of the inner value); i range: delta chi2 = 1 with chi2/N_beam = 1, statistical only; the rms column and V_rad/V_rot measure the departure from circular rotation.\n",
         "| map | R (kpc) | V_rot (km/s) | V range from i | V_rad | V_rad/V_rot | rms (km/s) | PA | i (deg) | i range | N_beam | log(g_obs/g_RAR) | range | i needed for the law |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, t in T.iterrows():
        need = f"{t.inc_RAR:.0f}" if np.isfinite(t.inc_RAR) else "none (V_proj > V_law)"
        L.append(f"| {t.src} | {t.R_kpc:.1f} | {t.Vrot:.1f} +/- {t.e_stat:.1f} | {t.V_lo:.0f}-{t.V_hi:.0f} | {t.Vrad:+.0f} | {abs(t.Vrad) / t.Vrot:.2f} | {t.rms:.1f} | {t.PA:.0f} | {t.inc:.0f} | "
                 f"{t.inc_lo:.0f}-{t.inc_hi:.0f} | {t.nbeam:.0f} | {t.res_RAR:+.2f} | {t.res_RAR_hi:+.2f} to {t.res_RAR_lo:+.2f} | {need} |")
    inn = T[~T.beyond_sparc & (T.src == "26arcsec")]
    out1 = T[T.beyond_sparc & (T.src == "26arcsec")]
    out6 = T[T.beyond_sparc & (T.src == "60arcsec")]
    L.append("")
    L.append(f"Inside R_last (26{q}): median i = {inn.inc.median():.0f} deg (16-84 %: {inn.inc.quantile(.16):.0f}-{inn.inc.quantile(.84):.0f}), "
             f"PA from {inn.PA.min():.0f} to {inn.PA.max():.0f}, median |V_rad|/V_rot = {(inn.Vrad.abs() / inn.Vrot).median():.2f}, median rms = {inn.rms.median():.1f} km/s, median residual {inn.res_RAR.median():+.3f} dex.")
    for lab, o in ((f"26{q}", out1), (f"60{q}", out6)):
        if len(o):
            L.append(f"Beyond R_last ({lab}, {len(o)} rings to {o.R_kpc.max():.1f} kpc): median i = {o.inc.median():.0f} deg "
                     f"({o.inc.quantile(.16):.0f}-{o.inc.quantile(.84):.0f}), PA from {o.PA.min():.0f} to {o.PA.max():.0f}, median |V_rad|/V_rot = {(o.Vrad.abs() / o.Vrot).median():.2f}, median rms = {o.rms.median():.1f} km/s; "
                     f"median residual {o.res_RAR.median():+.3f} dex (with i at its 1 sigma low: {o.res_RAR_hi.median():+.3f}); "
                     f"median i needed for the law: {o.inc_RAR.median():.0f} deg.")
    fixed = prev[prev.beyond_sparc & (prev.eV < 15)]
    L.append(f"For comparison, fixed i = {inc_g:.0f} deg (scripts/mhongoose_rings.py): median residual beyond R_last {fixed.res_RAR.median():+.3f} dex.")
    with open(ROOT / "results/ngc7793_warp.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n\nFigure: results/ngc7793_warp.png\n")
    print("\n".join(L))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for src, col in (("26arcsec", "tab:blue"), ("60arcsec", "tab:red")):
        t = T[T.src == src]
        ax[0].errorbar(t.R_kpc, t.inc, yerr=[np.maximum(t.inc - t.inc_lo, 0), np.maximum(t.inc_hi - t.inc, 0)], fmt="o", ms=3, color=col, label=f"i free, {src}")
        ax[0].plot(t.R_kpc, t.PA - pa_g + inc_g, "s", ms=2, color=col, alpha=0.5)
        ax[1].errorbar(t.R_kpc, t.Vrot, yerr=[np.maximum(t.Vrot - t.V_lo, 0), np.maximum(t.V_hi - t.Vrot, 0)], fmt="o", ms=3, color=col, label=f"i free, {src}")
    ax[0].axhline(inc_g, color="k", lw=0.8, ls="--", label=f"i fixed = {inc_g:.0f}")
    ax[0].plot(T.R_kpc, T.inc_RAR, "x", color="gray", label="i needed for the law")
    ax[0].set_xlabel("R (kpc)")
    ax[0].set_ylabel("inclination (deg); squares: PA - PA_in + i_in")
    ax[0].legend(fontsize=7)
    ax[1].errorbar(df.R, df.Vobs, yerr=df.eVobs, fmt="o", ms=3, color="k", label="SPARC")
    pf = prev[prev.eV < 15]
    ax[1].plot(pf.R_kpc, pf.Vrot, ".", color="gray", alpha=0.6, label=f"i fixed = {inc_g:.0f}")
    Ts = T.sort_values("R_kpc")
    ax[1].plot(Ts.R_kpc, Ts.Vbar, ":", color="green", label="V_bar")
    ax[1].plot(Ts.R_kpc, Ts.V_RAR, "-", color="orange", lw=1, label="acceleration law")
    ax[1].axvline(r_last, color="k", lw=0.6, ls="--")
    ax[1].set_xlabel("R (kpc)")
    ax[1].set_ylabel("V (km/s)")
    ax[1].legend(fontsize=7)
    ax[1].set_title("NGC 7793")
    fig.tight_layout()
    fig.savefig(ROOT / "results/ngc7793_warp.png", dpi=130)


if __name__ == "__main__":
    main()
