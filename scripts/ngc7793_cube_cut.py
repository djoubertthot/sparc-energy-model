# -*- coding: utf-8 -*-
"""Cut the MHONGOOSE 26" cube of NGC 7793 (J2357-32, DR2, r15_t00) to the galaxy: +/- 150 pixels (17.5') around the
centre, velocity window 90-360 km/s, channels rebinned by 4 (5.5 km/s). Same for the SoFiA mask (max over the bin).
Output: data/external/mhongoose/J2357-32/cube/ngc7793_r15_sub.fits (+ _mask.fits), small enough to model in memory."""
import gzip
import os
import shutil

import numpy as np
from astropy.io import fits

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = "data/external/mhongoose/J2357-32"
h = fits.getheader(f"{P}/cube/J2357-32_r15_t00_cube.fits")
d = fits.open(f"{P}/cube/J2357-32_r15_t00_cube.fits", memmap=True)[0].data
n3 = h["NAXIS3"]
v = (h["CRVAL3"] + (np.arange(n3) + 1 - h["CRPIX3"]) * h["CDELT3"]) / 1e3   # km/s, radio, descending
k0 = int(np.argmin(np.abs(v - 360.0))); k1 = int(np.argmin(np.abs(v - 90.0)))
k0 -= (k1 - k0 + 1) % 4 and 0
nch = ((k1 - k0 + 1) // 4) * 4
k1 = k0 + nch
x0, y0 = 386, 386; half = 150
sub = np.array(d[k0:k1, y0 - half:y0 + half + 1, x0 - half:x0 + half + 1], dtype=np.float32)
vsub = v[k0:k1]
sub4 = sub.reshape(nch // 4, 4, sub.shape[1], sub.shape[2]).mean(axis=1)
v4 = vsub.reshape(nch // 4, 4).mean(axis=1)
hdr = fits.Header()
hdr["BUNIT"] = "JY/BEAM"; hdr["BMAJ"] = h["BMAJ"]; hdr["BMIN"] = h["BMIN"]; hdr["BPA"] = h["BPA"]
hdr["CDELT1"] = h["CDELT1"]; hdr["CDELT2"] = h["CDELT2"]
hdr["CRPIX1"] = half + 1; hdr["CRPIX2"] = half + 1
hdr["CRVAL1"] = h["CRVAL1"] + (x0 - (h["CRPIX1"] - 1)) * h["CDELT1"] / np.cos(np.radians(h["CRVAL2"]))
hdr["CRVAL2"] = h["CRVAL2"] + (y0 - (h["CRPIX2"] - 1)) * h["CDELT2"]
hdr["CTYPE1"] = "RA---SIN"; hdr["CTYPE2"] = "DEC--SIN"; hdr["CTYPE3"] = "VRAD"; hdr["CUNIT3"] = "km/s"
hdr["CRPIX3"] = 1; hdr["CRVAL3"] = float(v4[0]); hdr["CDELT3"] = float(v4[1] - v4[0])
hdr["RESTFREQ"] = h["RESTFREQ"]; hdr["COMMENT"] = "cut of the MHONGOOSE DR2 r15_t00 cube, channels rebinned x4"
fits.writeto(f"{P}/cube/ngc7793_r15_sub.fits", sub4, hdr, overwrite=True)
print("sub-cube", sub4.shape, "v from", v4[0], "to", v4[-1], "dv", v4[1] - v4[0], "rms edge channels", float(np.nanstd(sub4[:3])), float(np.nanstd(sub4[-3:])))
# mask
gz = f"{P}/r15_t00_mask.fits.gz"
with gzip.open(gz, "rb") as fin, open(f"{P}/cube/mask_full.fits", "wb") as fout:
    shutil.copyfileobj(fin, fout)
hm = fits.getheader(f"{P}/cube/mask_full.fits")
m = fits.open(f"{P}/cube/mask_full.fits", memmap=True)[0].data
print("mask shape", m.shape, "naxis", hm["NAXIS"], "ctype3", hm.get("CTYPE3"), "crval3", hm.get("CRVAL3"), "cdelt3", hm.get("CDELT3"), "crpix3", hm.get("CRPIX3"))
if m.ndim == 3 and m.shape[0] == n3:
    msub = np.array(m[k0:k1, y0 - half:y0 + half + 1, x0 - half:x0 + half + 1]) > 0
    m4 = msub.reshape(nch // 4, 4, msub.shape[1], msub.shape[2]).max(axis=1).astype(np.int16)
    fits.writeto(f"{P}/cube/ngc7793_r15_sub_mask.fits", m4, hdr, overwrite=True)
    print("mask voxels", int(m4.sum()), "of", m4.size)
else:
    print("mask grid differs from the cube; not cut")
os.remove(f"{P}/cube/mask_full.fits")
