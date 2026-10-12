# -*- coding: utf-8 -*-
"""
NGC 7793: 3D tilted-ring model of the MHONGOOSE 26" cube (DR2, r15_t00), to decide whether the decline of the rotation
curve beyond 8 kpc seen in the 2D velocity fields (results/ngc7793_warp.md) is a genuine fall of the circular velocity or
the signature of a warped / disturbed outer disc.

Method (a compact 3DBarolo-like fit):
  - rings of width 0.9 kpc from 0.9 to 13.5 kpc; each ring is a cloud of points with surface density set by the observed
    flux in the ring (azimuthal normalisation), line-of-sight velocity Vsys + sin i (Vrot cos theta + Vrad sin theta),
    Gaussian velocity dispersion sigma, scale height 0.15 kpc; the model cube is convolved with the beam (34.8" x 26.3",
    PA 139) and compared with the data cube, channel by channel, in the ring's projected annulus (+/- 1 beam).
  - stage 1: inner disc (R <= 7.5 kpc): Vrot and sigma per ring, then the global geometry (centre, Vsys, PA, i), then
    Vrot and sigma again.
  - stage 2: outer disc (R > 7.5 kpc), three models fitted ring by ring outward with the inner disc fixed:
      M1 'decline': i and PA fixed to the inner values, Vrot, sigma, Vrad free;
      M2 'warp'   : i, PA, Vrot, sigma, Vrad free;
      M3 'law'    : Vrot fixed to the prediction of the law (V_RAR from results/ngc7793_warp.csv), i, PA, sigma, Vrad free.
    The residual of each model in the outer region tells whether a warp can reconcile the cube with the law.
Outputs: results/ngc7793_cube3d.md, results/ngc7793_cube3d.csv, results/ngc7793_cube3d.png
"""
import os
import sys
import time

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import minimize
from scipy.signal import fftconvolve

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = "data/external/mhongoose/J2357-32/cube"
D_MPC = 3.61
KPC_AS = D_MPC * 1e3 * np.pi / 180 / 3600          # kpc per arcsec
PIX_AS = 7.0
PIX_KPC = PIX_AS * KPC_AS
RMS = 1.0e-4                                       # Jy/beam per 5.5 km/s channel (edge channels 0.8-1.4e-4)
Z0_KPC = 0.15
RING_W = 0.9
R_INNER_MAX = 7.5
R_MAX = 12.0            # rings beyond carry < 1.5 Jy km/s at 26" and are unconstrained
SEED = 12345
CKPT = "results/ngc7793_cube3d_stage1.json"
RESUME = "--resume" in sys.argv

# ---------------------------------------------------------------- data
cube = fits.getdata(f"{P}/ngc7793_r15_sub.fits").astype(np.float64)
hdr = fits.getheader(f"{P}/ngc7793_r15_sub.fits")
mask = fits.getdata(f"{P}/ngc7793_r15_sub_mask.fits") > 0
NCH, NY, NX = cube.shape
VAX = hdr["CRVAL3"] + np.arange(NCH) * hdr["CDELT3"]           # km/s, descending
DV = abs(hdr["CDELT3"])
BMAJ, BMIN, BPA = hdr["BMAJ"] * 3600, hdr["BMIN"] * 3600, hdr["BPA"]
BEAM_PIX2 = np.pi * (BMAJ / PIX_AS) * (BMIN / PIX_AS) / (4 * np.log(2))
XC0, YC0 = 150.0, 150.0                                        # catalogue centre pixel
yy, xx = np.mgrid[0:NY, 0:NX]


def beam_kernel():
    n = 25
    gy, gx = np.mgrid[-n:n + 1, -n:n + 1].astype(float)
    dx, dy = -gx * PIX_AS, gy * PIX_AS                         # east (arcsec), north
    t = np.radians(BPA)                                        # BPA from north through east
    u = dx * np.sin(t) + dy * np.cos(t)                        # along major axis
    w = -dx * np.cos(t) + dy * np.sin(t)
    sa, sb = BMAJ / 2.3548, BMIN / 2.3548
    return np.exp(-0.5 * (u / sa) ** 2 - 0.5 * (w / sb) ** 2)   # peak 1 -> Jy/pixel in, Jy/beam out


KERNEL = beam_kernel()


def sky_to_disc(pa, inc, xc, yc):
    """deprojected radius (kpc) and azimuth for every pixel"""
    dx = -(xx - xc) * PIX_AS; dy = (yy - yc) * PIX_AS
    p = np.radians(pa)
    xr = dx * np.sin(p) + dy * np.cos(p); yr = -dx * np.cos(p) + dy * np.sin(p)
    r = np.sqrt(xr ** 2 + (yr / np.cos(np.radians(inc))) ** 2) * KPC_AS
    th = np.arctan2(yr / np.cos(np.radians(inc)), xr)
    return r, th


# ---------------------------------------------------------------- ring model
class Ring:
    def __init__(self, R, flux, npts):
        self.R, self.flux, self.npts = R, flux, npts
        rng = np.random.default_rng(SEED + int(R * 100))
        self.u_r = np.sqrt(rng.uniform((R - RING_W / 2) ** 2, (R + RING_W / 2) ** 2, npts))
        self.u_t = rng.uniform(0, 2 * np.pi, npts)
        self.u_z = rng.normal(0, Z0_KPC, npts)
        self.u_g = rng.normal(0, 1, npts)                       # velocity jitter draw (sub-channel smoothing)
        self.par = dict(vrot=50.0, vrad=0.0, sigma=10.0, inc=50.0, pa=290.0)
        self.cube = None

    def render(self, geo):
        vrot, vrad, sig, inc, pa = (self.par[k] for k in ("vrot", "vrad", "sigma", "inc", "pa"))
        xc, yc, vsys = geo
        ci, si = np.cos(np.radians(inc)), np.sin(np.radians(inc))
        xr = self.u_r * np.cos(self.u_t)
        yr = self.u_r * np.sin(self.u_t) * ci - self.u_z * si
        p = np.radians(pa)
        dx = (xr * np.sin(p) - yr * np.cos(p)) / KPC_AS         # arcsec east
        dy = (xr * np.cos(p) + yr * np.sin(p)) / KPC_AS
        ix = np.rint(xc - dx / PIX_AS).astype(int); iy = np.rint(yc + dy / PIX_AS).astype(int)
        ok = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        vlos = vsys + si * (vrot * np.cos(self.u_t) + vrad * np.sin(self.u_t))
        sig_t = np.sqrt(sig ** 2 + (DV / np.sqrt(12)) ** 2)
        f = self.flux / self.npts                                 # Jy km/s per point
        grid = np.zeros((NCH, NY, NX))
        jc = (vlos - VAX[0]) / (VAX[1] - VAX[0])                  # fractional channel
        K = int(np.ceil(3 * sig_t / DV)) + 1
        for k in range(-K, K + 1):
            j = np.rint(jc).astype(int) + k
            vj = VAX[0] + j * (VAX[1] - VAX[0])
            w = np.exp(-0.5 * ((vj - vlos) / sig_t) ** 2) / (np.sqrt(2 * np.pi) * sig_t)   # 1/(km/s)
            sel = ok & (j >= 0) & (j < NCH)
            np.add.at(grid, (j[sel], iy[sel], ix[sel]), f * w[sel])
        # convolve the channels that carry flux
        out = np.zeros_like(grid)
        chs = np.where(grid.reshape(NCH, -1).any(axis=1))[0]
        for j in chs:
            out[j] = fftconvolve(grid[j], KERNEL, mode="same")
        self.cube = out
        return out


def ring_flux(R, geo, pa, inc):
    """Jy km/s inside the ring annulus, from the masked data"""
    r, _ = sky_to_disc(pa, inc, geo[0], geo[1])
    sel = (np.abs(r - R) < RING_W / 2)
    m = mask & sel[None, :, :]
    return float((cube * m).sum() * DV / BEAM_PIX2), int(sel.sum())


def region(R, geo, pa, inc, margin_kpc):
    r, _ = sky_to_disc(pa, inc, geo[0], geo[1])
    return np.abs(r - R) < RING_W / 2 + margin_kpc


def residual(model, reg):
    """mean |data - model| / rms over the region, all channels (3DBarolo-like absolute residual)"""
    d = cube[:, reg]; m = model[:, reg]
    return float(np.abs(d - m).sum() / (RMS * d.size)), float(((d - m) ** 2).sum() / (RMS ** 2 * d.size))


def total_model(rings, geo, skip=None):
    tot = np.zeros_like(cube)
    for rg in rings:
        if rg is skip or rg.cube is None:
            continue
        tot += rg.cube
    return tot


# ---------------------------------------------------------------- fitting helpers
def fit_ring(rg, rings, geo, free, bounds, reg):
    """Nelder-Mead on the free parameters of one ring, others fixed; returns best residual."""
    others = total_model(rings, geo, skip=rg)
    x0 = np.array([rg.par[k] for k in free])

    def cost(x):
        for k, v, (lo, hi) in zip(free, x, bounds):
            rg.par[k] = float(np.clip(v, lo, hi))
        m = others + rg.render(geo)
        return residual(m, reg)[0]
    res = minimize(cost, x0, method="Nelder-Mead", options=dict(xatol=0.5, fatol=1e-4, maxiter=250, initial_simplex=None))
    cost(res.x)
    return res.fun


def main():
    t0 = time.time()
    warp = pd.read_csv("results/ngc7793_warp.csv")
    w26 = warp[warp.src == "26arcsec"].sort_values("R_kpc")
    w60 = warp[warp.src == "60arcsec"].sort_values("R_kpc")
    vrar_tab = pd.concat([w26[["R_kpc", "V_RAR", "Vbar"]], w60[w60.R_kpc > w26.R_kpc.max()][["R_kpc", "V_RAR", "Vbar"]]]).sort_values("R_kpc")
    v_rar = lambda R: np.interp(R, vrar_tab.R_kpc, vrar_tab.V_RAR)   # noqa: E731
    v_bar = lambda R: np.interp(R, vrar_tab.R_kpc, vrar_tab.Vbar)    # noqa: E731
    rings26 = pd.read_csv("results/mhongoose_NGC7793.csv"); rings26 = rings26[rings26.src == "26arcsec"]
    v2d = lambda R: np.interp(R, rings26.R_kpc, rings26.Vrot)        # noqa: E731

    # initial geometry from the 2D analysis: centre offset (+16" east, -4" north), Vsys 225.6, PA 291, i 51
    geo = [XC0 - 16.0 / PIX_AS, YC0 - 4.0 / PIX_AS, 225.6]
    PA0, INC0 = 291.0, 51.0
    radii = np.arange(RING_W, R_MAX + 0.01, RING_W)
    rings = []
    for R in radii:
        fl, npix = ring_flux(R, geo, PA0, INC0)
        rg = Ring(R, max(fl, 1e-4), int(min(40000, max(6000, 2500 * R))))
        rg.par.update(vrot=float(v2d(R)), sigma=9.0, inc=INC0, pa=PA0, vrad=0.0)
        rg.render(geo)
        rings.append(rg)
        print(f"ring R={R:.1f} flux {fl:.3f} Jy km/s npix {npix} vrot0 {v2d(R):.1f}", flush=True)
    inner = [rg for rg in rings if rg.R <= R_INNER_MAX]
    outer = [rg for rg in rings if rg.R > R_INNER_MAX]

    # ---------------- stage 1: inner disc
    import json
    if RESUME and os.path.exists(CKPT):
        ck = json.load(open(CKPT))
        geo, PA0, INC0 = ck["geo"], ck["PA0"], ck["INC0"]
        for rg in rings:
            rg.par["pa"], rg.par["inc"] = PA0, INC0
            if str(round(rg.R, 2)) in ck["inner"]:
                rg.par.update(ck["inner"][str(round(rg.R, 2))])
            fl, _ = ring_flux(rg.R, geo, PA0, INC0); rg.flux = max(fl, 1e-4)
            rg.render(geo)
        print("stage1 restored from", CKPT, flush=True)
    for it in range(0 if RESUME and os.path.exists(CKPT) else 2):
        for rg in inner:
            reg = region(rg.R, geo, PA0, INC0, 0.6)
            fit_ring(rg, rings, geo, ["vrot", "sigma"], [(0, 250), (3, 40)], reg)
        print(f"stage1 pass {it}: vrot inner {[round(r.par['vrot'],1) for r in inner]} ({time.time()-t0:.0f} s)", flush=True)
        if it == 0:
            # global geometry: centre, vsys, PA, i on the inner region
            regin = region(0.5 * R_INNER_MAX, geo, PA0, INC0, 0.5 * R_INNER_MAX)

            def gcost(x):
                g = [x[0], x[1], x[2]]
                for rg in inner:
                    rg.par["pa"], rg.par["inc"] = float(np.clip(x[3], 240, 340)), float(np.clip(x[4], 20, 80))
                    rg.render(g)
                return residual(total_model(inner, g), regin)[0]
            res = minimize(gcost, np.array([geo[0], geo[1], geo[2], PA0, INC0]), method="Nelder-Mead",
                           options=dict(xatol=0.3, fatol=1e-4, maxiter=120))
            geo = [float(res.x[0]), float(res.x[1]), float(res.x[2])]
            PA0, INC0 = float(np.clip(res.x[3], 240, 340)), float(np.clip(res.x[4], 20, 80))
            for rg in rings:
                rg.par["pa"], rg.par["inc"] = PA0, INC0
                fl, _ = ring_flux(rg.R, geo, PA0, INC0); rg.flux = max(fl, 1e-4)
                rg.render(geo)
            print(f"stage1 geometry: centre ({(XC0-geo[0])*PIX_AS:+.1f}\" E, {(geo[1]-YC0)*PIX_AS:+.1f}\" N), Vsys {geo[2]:.1f}, PA {PA0:.1f}, i {INC0:.1f} ({time.time()-t0:.0f} s)", flush=True)
    json.dump(dict(geo=geo, PA0=PA0, INC0=INC0, inner={str(round(rg.R, 2)): rg.par for rg in inner}), open(CKPT, "w"), indent=1)
    inner_res = residual(total_model(rings, geo), region(0.5 * R_INNER_MAX, geo, PA0, INC0, 0.5 * R_INNER_MAX))
    rows = [dict(R=rg.R, model="inner", vrot=rg.par["vrot"], sigma=rg.par["sigma"], inc=rg.par["inc"], pa=rg.par["pa"], vrad=0.0,
                 flux=rg.flux, resid=np.nan) for rg in inner]

    # ---------------- stage 2: outer disc, three models
    results = {}
    saved_inner = {rg.R: dict(rg.par) for rg in inner}
    for name, free, bounds in (("M1 decline", ["vrot", "sigma", "vrad"], [(0, 200), (3, 40), (-80, 80)]),
                               ("M2 warp", ["vrot", "sigma", "vrad", "inc", "pa"], [(0, 200), (3, 40), (-80, 80), (20, 80), (PA0 - 70, PA0 + 70)]),
                               ("M3 law", ["sigma", "vrad", "inc", "pa"], [(3, 40), (-80, 80), (20, 80), (PA0 - 70, PA0 + 70)]),
                               ("M4 law i33", ["sigma", "vrad", "pa"], [(3, 40), (-80, 80), (PA0 - 70, PA0 + 70)])):
        for rg in outer:                                  # reset
            rg.par.update(vrot=float(v2d(rg.R)) if not name.startswith("M3") and not name.startswith("M4") else float(v_rar(rg.R)),
                          sigma=9.0, inc=INC0 if not name.startswith("M4") else 33.0, pa=PA0, vrad=-15.0)   # vrad starts at the 2D value so that the simplex explores it
            rg.render(geo)
        tot_res = []
        for rg in outer:
            if name.startswith("M3") or name.startswith("M4"):
                rg.par["vrot"] = float(v_rar(rg.R))
            if name.startswith("M4"):
                rg.par["inc"] = 33.0
            reg = region(rg.R, geo, PA0, INC0, 0.7)
            r1 = fit_ring(rg, rings, geo, free, bounds, reg)
            if name in ("M2 warp", "M3 law"):             # second pass from the inner geometry as alternative start
                keep = dict(rg.par); rg.par.update(inc=INC0, pa=PA0); r2 = fit_ring(rg, rings, geo, free, bounds, reg)
                if r2 > r1:
                    rg.par.update(keep); rg.render(geo)
            rr = residual(total_model(rings, geo), reg)
            tot_res.append(rr[0])
            rows.append(dict(R=rg.R, model=name, vrot=rg.par["vrot"], sigma=rg.par["sigma"], inc=rg.par["inc"], pa=rg.par["pa"],
                             vrad=rg.par["vrad"], flux=rg.flux, resid=rr[0], chi2n=rr[1], v_rar=float(v_rar(rg.R)), v_bar=float(v_bar(rg.R))))
            print(f"{name} R={rg.R:.1f}: vrot {rg.par['vrot']:.1f} sigma {rg.par['sigma']:.1f} i {rg.par['inc']:.0f} pa {rg.par['pa']:.0f} vrad {rg.par['vrad']:+.0f} resid {rr[0]:.4f} ({time.time()-t0:.0f} s)", flush=True)
        regout = region(0.5 * (R_INNER_MAX + R_MAX), geo, PA0, INC0, 0.5 * (R_MAX - R_INNER_MAX) + 0.7)
        results[name] = dict(res=residual(total_model(rings, geo), regout), pars={rg.R: dict(rg.par) for rg in outer},
                             model=total_model(rings, geo).copy())
        print(f"== {name}: outer-region residual {results[name]['res'][0]:.4f}, chi2/N {results[name]['res'][1]:.3f} ({time.time()-t0:.0f} s)", flush=True)
    T = pd.DataFrame(rows); T.to_csv("results/ngc7793_cube3d.csv", index=False)

    # ---------------- report
    L = ["# NGC 7793: 3D tilted-ring model of the MHONGOOSE 26\" cube (2026-10-12)\n",
         f"Cube: DR2 r15_t00, beam {BMAJ:.1f}\" x {BMIN:.1f}\", channels rebinned to {DV:.1f} km/s, rms {RMS*1e3:.2f} mJy/beam, D = {D_MPC} Mpc. "
         f"Rings {RING_W} kpc wide, cloud model with {Z0_KPC} kpc scale height, azimuthal flux normalisation, beam-convolved, residual = mean |data - model| / rms "
         "over the ring annulus (+/- 1 beam) and all channels.\n",
         f"Inner disc (R <= {R_INNER_MAX} kpc): centre offset ({(XC0-geo[0])*PIX_AS:+.1f}\" E, {(geo[1]-YC0)*PIX_AS:+.1f}\" N) from the catalogue position, "
         f"Vsys = {geo[2]:.1f} km/s, PA = {PA0:.1f}, i = {INC0:.1f} (2D velocity fields: 225.6, 291, 51; THINGS: 290.1, 49.6). Residual {inner_res[0]:.4f}, chi2/N {inner_res[1]:.3f}.\n",
         "| R (kpc) | Vrot inner (km/s) | sigma | 2D Vrot | V_RAR |", "|---|---|---|---|---|"]
    for rg in inner:
        L.append(f"| {rg.R:.1f} | {rg.par['vrot']:.1f} | {rg.par['sigma']:.1f} | {v2d(rg.R):.1f} | {v_rar(rg.R):.1f} |")
    L.append("\n## Outer disc (R > 7.5 kpc): three models\n")
    L.append("| R (kpc) | flux (Jy km/s) | V_RAR | M1 decline: Vrot, sigma, Vrad, resid | M2 warp: Vrot, sigma, i, PA, Vrad, resid | M3 law: Vrot (fixed), sigma, i, PA, Vrad, resid | M4 law, i = 33: sigma, PA, Vrad, resid |")
    L.append("|---|---|---|---|---|---|---|")
    for rg in outer:
        a = results["M1 decline"]["pars"][rg.R]; b = results["M2 warp"]["pars"][rg.R]; c = results["M3 law"]["pars"][rg.R]; d4 = results["M4 law i33"]["pars"][rg.R]
        rr = {k: T[(T.model == k) & (T.R == rg.R)].resid.iloc[0] for k in results}
        L.append(f"| {rg.R:.1f} | {rg.flux:.3f} | {v_rar(rg.R):.0f} | {a['vrot']:.0f}, {a['sigma']:.0f}, {a['vrad']:+.0f}, {rr['M1 decline']:.4f} | "
                 f"{b['vrot']:.0f}, {b['sigma']:.0f}, {b['inc']:.0f}, {b['pa']:.0f}, {b['vrad']:+.0f}, {rr['M2 warp']:.4f} | "
                 f"{c['vrot']:.0f}, {c['sigma']:.0f}, {c['inc']:.0f}, {c['pa']:.0f}, {c['vrad']:+.0f}, {rr['M3 law']:.4f} | "
                 f"{d4['sigma']:.0f}, {d4['pa']:.0f}, {d4['vrad']:+.0f}, {rr['M4 law i33']:.4f} |")
    L.append("")
    for name in results:
        L.append(f"- {name}: residual over the whole outer region {results[name]['res'][0]:.4f} (chi2/N {results[name]['res'][1]:.3f})")
    r1, r2, r3, r4 = (results[k]["res"][0] for k in ("M1 decline", "M2 warp", "M3 law", "M4 law i33"))
    L.append(f"- M3 (law, geometry free) relative to M2 (free warp): {100*(r3/r2-1):+.1f} % in residual; M4 (law, i = 33 deg) relative to M2: {100*(r4/r2-1):+.1f} %; "
             f"M1 (fixed geometry) relative to M2: {100*(r1/r2-1):+.1f} %.")
    # pressure-support corrected circular velocity of the free-warp solution
    L.append("\nAsymmetric-drift correction of the M2 solution (V_c^2 = Vrot^2 - sigma^2 dln(Sigma sigma^2)/dlnR, Sigma from the ring fluxes):")
    Rs = np.array([rg.R for rg in outer]); Sig = np.array([rg.flux / rg.R for rg in outer]); sg = np.array([results["M2 warp"]["pars"][rg.R]["sigma"] for rg in outer])
    vr = np.array([results["M2 warp"]["pars"][rg.R]["vrot"] for rg in outer])
    dln = np.gradient(np.log(Sig * sg ** 2), np.log(Rs))
    vc = np.sqrt(np.maximum(vr ** 2 - sg ** 2 * dln, 0))
    for R, v, s, c_, vrr in zip(Rs, vr, sg, vc, [v_rar(R) for R in Rs]):
        L.append(f"- R = {R:.1f} kpc: Vrot {v:.0f}, sigma {s:.0f}, V_c {c_:.0f} km/s, law {vrr:.0f}, log(g_obs/g_RAR) = {2*np.log10(max(c_,1)/vrr):+.2f} dex")
    L.append("\nReading: if M3 reaches the residual of M2, a warp reconciles the cube with the law and the 2D decline was geometric; "
             "if M3 stays well above M2 while M1 is close to M2, the decline is genuine and no inclination or position-angle change rescues the law.")
    L.append("\nFigure: results/ngc7793_cube3d.png (major-axis position-velocity diagram, data versus M2 and M3; rotation curves).")
    open("results/ngc7793_cube3d.md", "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))

    # ---------------- figure: PV along the inner major axis + rotation curves
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p = np.radians(PA0)
    L_as = np.arange(-900, 901, PIX_AS)                        # offset along the major axis (arcsec)
    px = geo[0] - (L_as * np.sin(p)) / PIX_AS; py = geo[1] + (L_as * np.cos(p)) / PIX_AS
    ix = np.clip(np.rint(px).astype(int), 0, NX - 1); iy = np.clip(np.rint(py).astype(int), 0, NY - 1)
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    pv_d = cube[:, iy, ix]
    for a, (name, lab) in zip(ax[:2], (("M2 warp", "M2: free warp"), ("M3 law", "M3: law's curve, geometry free"))):
        pv_m = results[name]["model"][:, iy, ix]
        a.imshow(pv_d, aspect="auto", origin="lower", extent=[L_as[0] * KPC_AS, L_as[-1] * KPC_AS, VAX[-1], VAX[0]], cmap="Greys", vmin=0, vmax=0.003)
        a.contour(L_as * KPC_AS, VAX, pv_m, levels=[0.0003, 0.0007, 0.0015, 0.003], colors="r", linewidths=0.8)
        a.contour(L_as * KPC_AS, VAX, pv_d, levels=[0.0003, 0.0007, 0.0015, 0.003], colors="k", linewidths=0.5)
        a.set_xlabel("offset along the major axis (kpc)"); a.set_ylabel("V (km/s)"); a.set_title(f"PV: data (grey, black) vs {lab} (red)")
        a.axvline(R_INNER_MAX, ls=":", color="b"); a.axvline(-R_INNER_MAX, ls=":", color="b")
    a = ax[2]
    a.plot([rg.R for rg in inner], [rg.par["vrot"] for rg in inner], "ko-", label="3D inner disc")
    for name, col, mk in (("M1 decline", "C0", "s"), ("M2 warp", "C1", "o"), ("M3 law", "C3", "^")):
        a.plot([rg.R for rg in outer], [results[name]["pars"][rg.R]["vrot"] for rg in outer], marker=mk, color=col, ls="-", label=name)
    a.plot(vrar_tab.R_kpc, vrar_tab.V_RAR, "r--", label="law (V_RAR)"); a.plot(vrar_tab.R_kpc, vrar_tab.Vbar, "g:", label="baryons")
    a.plot(w26.R_kpc, w26.Vrot, "x", color="grey", label="2D rings, i free")
    a.set_xlabel("R (kpc)"); a.set_ylabel("V (km/s)"); a.legend(fontsize=7); a.set_ylim(0, 160); a.set_title("NGC 7793 rotation curves")
    plt.tight_layout(); plt.savefig("results/ngc7793_cube3d.png", dpi=130)
    print("figure written", f"({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
