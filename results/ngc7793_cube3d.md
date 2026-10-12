# NGC 7793: 3D tilted-ring model of the MHONGOOSE 26" cube (2026-10-12)

Cube: DR2 r15_t00, beam 34.8" x 26.3", channels rebinned to 5.5 km/s, rms 0.10 mJy/beam, D = 3.61 Mpc. Rings 0.9 kpc wide, cloud model with 0.15 kpc scale height, azimuthal flux normalisation, beam-convolved, residual = mean |data - model| / rms over the ring annulus (+/- 1 beam) and all channels.

Inner disc (R <= 7.5 kpc): centre offset (+0.3" E, +4.8" N) from the catalogue position, Vsys = 227.4 km/s, PA = 291.6, i = 50.5 (2D velocity fields: 225.6, 291, 51; THINGS: 290.1, 49.6). Residual 9.9199, chi2/N 715.089.

| R (kpc) | Vrot inner (km/s) | sigma | 2D Vrot | V_RAR |
|---|---|---|---|---|
| 0.9 | 57.0 | 12.3 | 50.8 | 68.0 |
| 1.8 | 79.7 | 11.7 | 74.1 | 77.4 |
| 2.7 | 98.3 | 11.6 | 95.7 | 92.2 |
| 3.6 | 111.0 | 11.4 | 104.9 | 99.4 |
| 4.5 | 113.6 | 10.9 | 106.4 | 101.2 |
| 5.4 | 104.5 | 10.9 | 101.6 | 102.8 |
| 6.3 | 97.4 | 13.1 | 93.0 | 103.1 |
| 7.2 | 86.1 | 19.0 | 85.3 | 102.6 |

## Outer disc (R > 7.5 kpc): three models

| R (kpc) | flux (Jy km/s) | V_RAR | M1 decline: Vrot, sigma, Vrad, resid | M2 warp: Vrot, sigma, i, PA, Vrad, resid | M3 law: Vrot (fixed), sigma, i, PA, Vrad, resid | M4 law, i = 33: sigma, PA, Vrad, resid |
|---|---|---|---|---|---|---|
| 8.1 | 11.189 | 101 | 89, 10, -19, 3.3735 | 89, 9, 51, 296, -13, 3.2717 | 101, 13, 50, 295, -13, 3.7852 | 15, 309, -8, 3.8346 |
| 9.0 | 8.627 | 101 | 68, 16, -23, 2.6262 | 69, 16, 51, 295, -17, 2.4537 | 101, 12, 27, 362, +66, 3.0031 | 18, 282, -9, 3.0508 |
| 9.9 | 6.499 | 101 | 60, 20, -18, 2.2202 | 70, 13, 50, 293, -23, 2.1368 | 101, 26, 53, 295, -8, 2.5379 | 26, 362, +15, 2.3670 |
| 10.8 | 3.642 | 100 | 27, 27, -18, 1.6757 | 39, 13, 58, 362, +33, 1.5955 | 100, 40, 49, 347, -26, 1.8907 | 40, 362, -28, 1.8091 |
| 11.7 | 2.071 | 100 | 34, 22, -17, 1.1682 | 28, 11, 68, 342, -8, 1.1120 | 100, 40, 60, 334, +30, 1.2745 | 40, 362, +44, 1.3689 |

- M1 decline: residual over the whole outer region 2.2792 (chi2/N 42.292)
- M2 warp: residual over the whole outer region 2.2295 (chi2/N 41.510)
- M3 law: residual over the whole outer region 2.4984 (chi2/N 48.369)
- M4 law i33: residual over the whole outer region 2.5772 (chi2/N 53.527)
- M3 (law, geometry free) relative to M2 (free warp): +12.1 % in residual; M4 (law, i = 33 deg) relative to M2: +15.6 %; M1 (fixed geometry) relative to M2: +2.2 %.

Asymmetric-drift correction of the M2 solution (V_c^2 = Vrot^2 - sigma^2 dln(Sigma sigma^2)/dlnR, Sigma from the ring fluxes):
- R = 8.1 kpc: Vrot 89, sigma 9, V_c 86 km/s, law 101, log(g_obs/g_RAR) = -0.15 dex
- R = 9.0 kpc: Vrot 69, sigma 16, V_c 71 km/s, law 101, log(g_obs/g_RAR) = -0.31 dex
- R = 9.9 kpc: Vrot 70, sigma 13, V_c 79 km/s, law 101, log(g_obs/g_RAR) = -0.21 dex
- R = 10.8 kpc: Vrot 39, sigma 13, V_c 59 km/s, law 100, log(g_obs/g_RAR) = -0.46 dex
- R = 11.7 kpc: Vrot 28, sigma 11, V_c 49 km/s, law 100, log(g_obs/g_RAR) = -0.62 dex

Reading: if M3 reaches the residual of M2, a warp reconciles the cube with the law and the 2D decline was geometric; if M3 stays well above M2 while M1 is close to M2, the decline is genuine and no inclination or position-angle change rescues the law.

Figure: results/ngc7793_cube3d.png (major-axis position-velocity diagram, data versus M2 and M3; rotation curves).

## Interpretation (2026-10-12)

- The decline is genuine. With the inner disc fixed (PA 291.6, i 50.5, Vsys 227.4, all within 2 deg and 2 km/s of THINGS), the
  circular velocity falls from 86 km/s at 8.1 kpc to 71-79 km/s at 9-10 kpc (asymmetric drift included), where the law predicts
  101 km/s: log(g_obs/g_RAR) = -0.15 to -0.31 dex at g_bar = 0.03-0.04 a0. The law's curve fitted with i, PA, sigma and V_rad free
  (M3) or with i forced to the 33 deg the 2D analysis required (M4) leaves residuals 12 and 16 % above the free-warp model, ring by
  ring, and needs sigma of 26-40 km/s beyond 9.9 kpc to spread the model velocities over the observed ones. The 2D result (-0.35 dex)
  is confirmed and slightly reduced.
- The outer disc is not a clean warp. Between 8 and 10 kpc the free fit keeps PA within 5 deg and i within 1 deg of the inner disc
  (M1 and M2 differ by 2 %), but a radial term of -13 to -23 km/s (20-30 % of the rotation) is required by every model, and beyond
  10.8 kpc (fluxes < 4 Jy km/s) the position angle swings by 50-70 deg: the outermost gas is disturbed, as the 2D fields showed.
- Not an external-field effect. The deficit has the size the Chae et al. (2020) function gives for the galaxy's own 2M++ field
  (e_N = 1.18e-4 a0, boosted e = 0.011 a0: -0.22 dex at 0.035 a0), and the SPARC curve alone preferred e_N = 0.02 (Delta chi2 6.3).
  But NGC 7793's environment field is at the 17th percentile of SPARC (median 1.47e-4), while DDO 161 and ESO 444-G084, with
  stronger fields (1.5e-4), show no deficit at 0.007-0.010 a0 (ESO 444-G084: +0.08 dex) where the same function would give -0.4 dex;
  and the lensing flatness of isolated galaxies (results/lensing_vflat.md) excludes e_N > 4e-5 a0. The deficit belongs to NGC 7793,
  not to its environment.
- What remains: a disc whose outer HI rotates 20-30 % slower than the law predicts at 8-10 kpc, with 20-30 % radial motions and a
  disturbed envelope beyond 10 kpc. Either the law fails here by 0.2-0.3 dex (two to three times the RAR scatter) or the outer gas
  is not in equilibrium; the radial term and the position-angle swing favour the second without proving it. Among the three
  MHONGOOSE discs with deep extensions this is the only such case; the two others stay on the law to 1.3-1.4 R_last.
- Method limits: cloud model with azimuthal normalisation and a thin disc (0.15 kpc), Nelder-Mead fits ring by ring from inside
  out (local minima possible; the 9.9 kpc ring moved between runs), residual = mean absolute deviation dominated by the bright
  pixels; the ring fluxes beyond 10.8 kpc are 1-4 Jy km/s, so those rings are indicative only.