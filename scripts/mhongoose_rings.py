"""Tilted-ring rotation curves from the MHONGOOSE moment maps, extended beyond SPARC.

For each of the three deep galaxies shared with SPARC (DDO 161 = J1303-17b, ESO 444-G084 = J1337-28,
NGC 7793 = J2357-32):
1. On the 26" moment-1 map, fit centre, systemic velocity, position angle and inclination with
   rings inside the SPARC radial range (geometry stage), inclination bounded to +/- 15 deg of SPARC.
2. With the geometry fixed, fit V_rot ring by ring on the 26" map, then on the 60" map
   (10^18 cm^-2 sensitivity) out to the last ring with enough pixels, using only |cos theta| > 0.4.
3. Gas surface density from the primary-beam-corrected moment-0 (N_HI = 1.106e24 S / (theta_maj theta_min)
   with S in Jy/beam km/s, theta in arcsec), azimuthally averaged in the same rings, M_HI(<R) x 1.33,
   V_gas extended beyond SPARC's last point with V_gas^2 = G M_gas(<R)/R scaled to match SPARC at R_last.
4. Compare with the SPARC curve (same distance as SPARC) and compute the RAR residual in the
   extension, with the SPARC stellar components held at their last value (negligible there).
Writes results/mhongoose_rings.md, results/mhongoose_<gal>.csv, results/mhongoose_rings.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import final_baryonic_model as fb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G, MSUN, KPC = fb.G, fb.MSUN, fb.KPC
PC = KPC / 1e3
A0 = 1.2e-10
cat = pd.read_csv(ROOT / "data/sparc_sample/sparc_enriched.csv").set_index("Galaxy")
GAL = {"J1303-17b": "DDO161", "J1337-28": "ESO444-G084", "J2357-32": "NGC7793"}


def load(gal, res, kind):
    p = ROOT / f"data/external/mhongoose/{gal}/{res}_{res}_{kind}.fits"
    h = fits.getheader(p)
    d = np.squeeze(fits.getdata(p)).astype(float)
    return h, d


def sky_grid(h, d):
    w = WCS(h).celestial
    ny, nx = d.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    ra, dec = w.all_pix2world(xx, yy, 0)
    return np.mod(ra, 360.0), dec


def ring_coords(ra, dec, ra0, dec0, pa_deg, inc_deg):
    """Deprojected radius (arcsec) and azimuth (rad) in the disc plane."""
    dx = (ra - ra0) * np.cos(np.radians(dec0)) * 3600.0    # +x east (RA grows eastward)
    dy = (dec - dec0) * 3600.0
    pa = np.radians(pa_deg)
    xr = dx * np.sin(pa) + dy * np.cos(pa)       # along the major axis (PA measured N through E of the receding side)
    yr = -dx * np.cos(pa) + dy * np.sin(pa)
    r = np.sqrt(xr**2 + (yr / np.cos(np.radians(inc_deg))) ** 2)
    theta = np.arctan2(yr / np.cos(np.radians(inc_deg)), xr)
    return r, theta


def fit_geometry(ra, dec, v, r_in, r_out, ra0, dec0, vsys0, pa0, inc0):
    """Coarse PA scan with the SPARC centre and inclination, then a bounded refinement."""
    r0_, _ = ring_coords(ra, dec, ra0, dec0, 0.0, inc0)
    m = np.isfinite(v) & (r0_ > r_in) & (r0_ < r_out)
    ra_m, dec_m, v_m = ra[m], dec[m], v[m]
    best = None
    for pa in np.arange(0, 360, 5.0):
        _, th = ring_coords(ra_m, dec_m, ra0, dec0, pa, inc0)
        x = np.cos(th) * np.sin(np.radians(inc0))
        A = np.vstack([np.ones_like(x), x]).T
        coef, *_ = np.linalg.lstsq(A, v_m, rcond=None)
        chi = np.sum((v_m - A @ coef) ** 2)
        if coef[1] > 0 and (best is None or chi < best[0]):
            best = (chi, pa, coef[0], coef[1])
    _, pa_best, vsys_best, vrot_best = best

    def resid(p):
        dra, ddec, vsys, pa, inc, vrot = p
        _, th = ring_coords(ra_m, dec_m, ra0 + dra / 3600 / np.cos(np.radians(dec0)), dec0 + ddec / 3600, pa, inc)
        w = np.clip(np.abs(np.cos(th)), 0.0, 1.0)
        return w * (v_m - vsys - vrot * np.cos(th) * np.sin(np.radians(inc))) / 1e3

    p0 = [0.0, 0.0, vsys_best, pa_best, inc0, vrot_best]
    lo = [-20, -20, vsys_best - 30e3, pa_best - 25, max(inc0 - 10, 20), 1e3]
    hi = [20, 20, vsys_best + 30e3, pa_best + 25, min(inc0 + 10, 85), 400e3]
    return least_squares(resid, p0, bounds=(lo, hi)).x


def fit_rings(ra, dec, v, geom, radii, width, min_pix, pa_free=False, nhi_ring=None, nhi_min=1e18):
    """V_rot per ring; with pa_free the position angle is scanned within +/- 45 deg per ring (warps)."""
    dra, ddec, vsys, pa_g, inc, _ = geom
    ra_c, dec_c = geom_ra0 + dra / 3600 / np.cos(np.radians(geom_dec0)), geom_dec0 + ddec / 3600
    out = []
    for R in radii:
        best = None
        for pa in (np.arange(pa_g - 60, pa_g + 61, 5.0) if pa_free else [pa_g]):
            r, th = ring_coords(ra, dec, ra_c, dec_c, pa, inc)
            m = np.isfinite(v) & (np.abs(r - R) < width / 2) & (np.abs(np.cos(th)) > 0.4)
            if m.sum() < min_pix:
                continue
            cover = len(np.unique((np.degrees(th[m]) % 360 // 15).astype(int))) / 24.0
            if cover < 0.5:
                continue
            if nhi_ring is not None and nhi_ring(R) < nhi_min:
                continue
            c = np.cos(th[m]) * np.sin(np.radians(inc))
            vr = np.sum(c * (v[m] - vsys)) / np.sum(c**2)
            res = v[m] - vsys - vr * c
            chi = np.sum(res**2) / m.sum()
            if best is None or chi < best[0]:
                err = np.sqrt(np.sum(res**2) / max(m.sum() - 1, 1) / np.sum(c**2))
                ma, mr = (np.cos(th[m]) < 0), (np.cos(th[m]) > 0)
                va = np.sum(c[ma] * (v[m][ma] - vsys)) / np.sum(c[ma] ** 2) if ma.sum() > 3 else vr
                vb = np.sum(c[mr] * (v[m][mr] - vsys)) / np.sum(c[mr] ** 2) if mr.sum() > 3 else vr
                best = (chi, R, vr / 1e3, max(err, abs(va - vb) / 2) / 1e3, int(m.sum()), pa)
        if best is not None:
            out.append(best[1:])
    return out



def main():
    global geom_ra0, geom_dec0
    lines, figs = [], []
    for gal, sp in GAL.items():
        c = cat.loc[sp]
        D = float(c.D_Mpc)
        kpc_per_arcsec = D * 1e3 * np.pi / 180 / 3600
        df = fb.load_rotmod(ROOT / f"data/rotmod/{sp}_rotmod.dat")
        r_last = df.R.values[-1]
        # ---- geometry on the 26" map
        h1, v1 = load(gal, "r10_t00", "mom1")
        ra1, dec1 = sky_grid(h1, v1)
        geom_ra0, geom_dec0 = float(c.RA_deg), float(c.Dec_deg)
        vsys0 = np.nanmedian(v1)
        pa0 = 90.0
        best = None
        for pa_try in (0, 60, 120, 180, 240, 300):
            g = fit_geometry(ra1, dec1, v1, 0.5 * r_last / kpc_per_arcsec, r_last / kpc_per_arcsec, geom_ra0, geom_dec0, vsys0, pa_try, float(c.Inc_deg))
            r, th = ring_coords(ra1, dec1, geom_ra0 + g[0] / 3600 / np.cos(np.radians(geom_dec0)), geom_dec0 + g[1] / 3600, g[3], g[4])
            m = np.isfinite(v1) & (r < r_last / kpc_per_arcsec) & (np.abs(np.cos(th)) > 0.3)
            chi = np.sum((v1[m] - g[2] - g[5] * np.cos(th[m]) * np.sin(np.radians(g[4]))) ** 2)
            if best is None or chi < best[0]:
                best = (chi, g)
        geom = best[1]
        if geom[5] < 0:
            geom[3] = (geom[3] + 180) % 360
            geom[5] = -geom[5]
        # ---- gas from the 60" moment-0 (pb corrected)
        h0, m0 = load(gal, "r05_t60", "mom0_pb")
        bm = h0["BMAJ"] * 3600 * h0["BMIN"] * 3600
        nhi = 1.106e24 * (m0 / 1e3) / bm                        # cm^-2, S in Jy/beam km/s
        sig = nhi * 1.6726e-24 / 1.989e33 * (3.0857e18) ** 2     # Msun/pc^2
        ra0m, dec0m = sky_grid(h0, m0)
        rm, _ = ring_coords(ra0m, dec0m, geom_ra0 + geom[0] / 3600 / np.cos(np.radians(geom_dec0)), geom_dec0 + geom[1] / 3600, geom[3], geom[4])
        nhi_prof = lambda r_as: np.nanmean(nhi[(np.abs(rm - r_as) < 0.5 * h0['BMAJ'] * 3600) & np.isfinite(nhi)]) if np.any(np.abs(rm - r_as) < 0.5 * h0['BMAJ'] * 3600) else 0.0
        pix_area = (abs(h0["CDELT2"]) * 3600 * kpc_per_arcsec * 1e3) ** 2 / np.cos(np.radians(geom[4]))  # pc^2 deprojected
        edges = np.arange(0, rm.max(), h0['BMAJ'] * 3600 / 2)
        prof = []
        mcum = 0.0
        for lo_, hi_ in zip(edges[:-1], edges[1:]):
            mm = (rm >= lo_) & (rm < hi_) & np.isfinite(sig)
            if mm.sum() == 0:
                prof.append((hi_, 0.0, mcum)); continue
            mcum += np.sum(np.clip(sig[mm], 0, None)) * pix_area * 1.33
            prof.append((hi_, np.nanmean(np.clip(sig[mm], 0, None)), mcum))
        prof = pd.DataFrame(prof, columns=["r_arcsec", "sig", "mgas_cum"])
        vgas_sph = lambda r_as: np.sqrt(G * np.interp(r_as, prof.r_arcsec, prof.mgas_cum) * MSUN / max(r_as * kpc_per_arcsec * KPC, 1e15)) / 1e3
        scale_gas = abs(df.Vgas.values[-1]) / max(vgas_sph(r_last / kpc_per_arcsec), 1e-3)
        # ---- rings on the 26" map, then the 60" map
        beam1 = h1["BMAJ"] * 3600
        radii1 = np.arange(beam1, 2.5 * r_last / kpc_per_arcsec, beam1 / 2)
        rings = fit_rings(ra1, dec1, v1, geom, radii1, beam1, 25, pa_free=True, nhi_ring=nhi_prof, nhi_min=5e18)
        h6, v6 = load(gal, "r05_t60", "mom1")
        ra6, dec6 = sky_grid(h6, v6)
        beam6 = h6["BMAJ"] * 3600
        radii6 = np.arange(beam6, 6 * r_last / kpc_per_arcsec, beam6 / 2)
        rings6 = fit_rings(ra6, dec6, v6, geom, radii6, beam6, 12, pa_free=True, nhi_ring=nhi_prof, nhi_min=1e18)
        # ---- assemble
        rows = []
        for src, rr in (("26arcsec", rings), ("60arcsec", rings6)):
            for R, vr, er, n, pa_ring in rr:
                rk = R * kpc_per_arcsec
                vgas = df.Vgas.values[-1] * 0 + (np.interp(rk, df.R.values, df.Vgas.values) if rk <= r_last else vgas_sph(R) * scale_gas)
                vdisk = np.interp(rk, df.R.values, df.Vdisk.values) if rk <= r_last else df.Vdisk.values[-1] * np.sqrt(r_last / rk)
                vbul = np.interp(rk, df.R.values, df.Vbul.values) if rk <= r_last else df.Vbul.values[-1] * np.sqrt(r_last / rk)
                yd = float(c.Y36) if pd.notna(c.Y36) else 0.5
                vb2 = (vgas * abs(vgas) + yd * vdisk * abs(vdisk) + 1.4 * yd * vbul * abs(vbul)) * 1e6
                gbar = max(vb2, 1.0) / (rk * KPC)
                gobs = (vr * 1e3) ** 2 / (rk * KPC)
                grar = gbar / (1 - np.exp(-np.sqrt(gbar / A0)))
                rows.append(dict(src=src, R_kpc=rk, Vrot=vr, eV=er, npix=n, PA=pa_ring, dPA=pa_ring - geom[3], Vbar=np.sqrt(max(vb2, 0)) / 1e3, log_gbar=np.log10(gbar),
                                 res_RAR=np.log10(gobs / grar), beyond_sparc=rk > r_last))
        T = pd.DataFrame(rows)
        T.to_csv(ROOT / f"results/mhongoose_{sp}.csv", index=False)
        ext = T[T.beyond_sparc & (T.eV < 15)]
        inside = T[~T.beyond_sparc & (T.src == "26arcsec")]
        vs = np.interp(inside.R_kpc, df.R.values, df.Vobs.values) if len(inside) else np.array([])
        lines.append(f"\n## {sp} ({gal}), SPARC D = {D} Mpc, i = {c.Inc_deg} +/- {c.e_Inc_deg}, R_last = {r_last:.1f} kpc, N_sparc = {len(df)}")
        lines.append(f"geometry fitted on the 26\" map inside R_last: Vsys = {geom[2] / 1e3:.1f} km/s, PA = {geom[3]:.0f} deg, i = {geom[4]:.0f} deg (SPARC {c.Inc_deg}), centre offset ({geom[0]:.0f}\", {geom[1]:.0f}\")")
        if len(inside):
            lines.append(f"inside R_last: MHONGOOSE 26\" rings vs SPARC, median V_MH / V_SPARC = {np.median(inside.Vrot.values / vs):.3f} ({len(inside)} rings)")
        lines.append(f"extension: {len(ext)} rings beyond R_last up to {ext.R_kpc.max() if len(ext) else 0:.1f} kpc = {ext.R_kpc.max() / r_last if len(ext) else 0:.2f} R_last; "
                     f"outer V = {ext.Vrot.values[-1] if len(ext) else float('nan'):.1f} +/- {ext.eV.values[-1] if len(ext) else float('nan'):.1f} km/s; "
                     f"log g_bar there = {ext.log_gbar.values[-1] if len(ext) else float('nan'):.2f} ({10 ** ext.log_gbar.values[-1] / A0 if len(ext) else float('nan'):.3f} a0)")
        if len(ext):
            lines.append(f"RAR residual log(g_obs/g_RAR) in the extension: median {ext.res_RAR.median():+.3f}, range {ext.res_RAR.min():+.2f} to {ext.res_RAR.max():+.2f}; inside R_last (26\"): median {inside.res_RAR.median():+.3f}")
        figs.append((sp, df, T, r_last))

    import matplotlib  # noqa: E402
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    for a, (sp, df, T, r_last) in zip(ax, figs):
        a.errorbar(df.R, df.Vobs, yerr=df.eVobs, fmt="o", color="k", ms=4, label="SPARC")
        for src, col in (("26arcsec", "tab:blue"), ("60arcsec", "tab:red")):
            s = T[(T.src == src) & (T.eV < 15)]
            a.errorbar(s.R_kpc, s.Vrot, yerr=s.eV, fmt="s", ms=3, color=col, alpha=0.8, label=f"MHONGOOSE {src}")
        s = T[T.src == "60arcsec"]
        rr = np.linspace(0.2, max(s.R_kpc.max(), r_last), 100)
        a.plot(T.R_kpc[T.src == "60arcsec"], T.Vbar[T.src == "60arcsec"], ":", color="gray", label="V_bar")
        a.axvline(r_last, color="k", lw=0.6, ls="--")
        a.set_title(sp); a.set_xlabel("R (kpc)"); a.set_ylabel("V (km/s)"); a.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(ROOT / "results/mhongoose_rings.png", dpi=130)
    with open(ROOT / "results/mhongoose_rings.md", "w", encoding="utf-8") as fh:
        fh.write("# MHONGOOSE tilted-ring curves vs SPARC (2026-10-10)\n" + "\n".join(lines) + "\n\nFigure: results/mhongoose_rings.png\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
