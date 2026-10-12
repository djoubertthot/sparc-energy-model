# a0 at high redshift, second pass: measured cold gas and curve shapes (2026-10-12)

Same test as results/a0_highz.md (r = a0(z)/a0(0); H_A: r = H(z)/H0; H_B: r = 1). For the f_DM samples the observable is now nu_obs = 1/(1 - f_DM) against nu(g_bar/(r a0)) with g_bar from the authors' own baryonic masses and sizes (exponential disc + Hernquist bulge); velocities enter only through the authors' fits. New samples: NOEMA3D (CO gas measured), CRISTAL ([CII] discs at z = 4.4-5.7), ALPAKA+JWST (three discs with curves to 4 Re and fitted gas normalisation), the Rizzo et al. 2021 lensed discs with their own baryonic masses. Systematics as before (sys_M: all masses +/- 0.2 dex; 'mass shift for H_A': uniform change making a0 ~ H(z) the best fit).

## r = a0(z)/a0(0) per sample

| sample | z | N | g_bar/a0 med | nu_obs med | r (68 %) | r, errors rescaled | r_A | chi2 min | H_B: chi2 - min | H_A: chi2 - min | (A-B) rescaled | sys_M +0.2/-0.2 | mass shift for H_A |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RC100 | 0.6-1.0 | 26 | 1.08 | 1.61 | 1.12 (0.89-1.41) | 0.95-1.35 | 1.58 | 15.6 | +0.3 | +2.6 | +4 | -0.68/+0.43 | -0.06 dex |
| RC100 | 1.0-1.7 | 33 | 1.72 | 1.39 | 1.12 (0.91-1.35) | 0.95-1.32 | 2.26 | 20.5 | +0.3 | +13.7 | +22 | -1.02/+0.47 | -0.12 dex |
| RC100 | 1.7-2.6 | 41 | 1.66 | 1.35 | 0.89 (0.71-1.07) | 0.76-1.02 | 3.27 | 23.3 | +0.3 | +51.7 | +91 | -0.95/+0.51 | -0.24 dex |
| MSA-3D | 0.5-1.7 | 30 | 0.83 | 2.54 | 2.29 (2.00-2.63) | 1.95-2.69 | 1.90 | 39.2 | +25.7 | +0.6 | -19 | -0.39/+0.31 | +0.04 dex |
| KMOS3D | 0.6-1.1 | 65 | 0.94 | 2.32 | 2.75 (2.34-3.24) | 2.45-3.16 | 1.66 | 45.5 | +26.1 | +8.7 | -25 | -0.47/+0.32 | +0.10 dex |
| KMOS3D | 1.1-1.9 | 24 | 2.36 | 1.45 | 0.98 (0.47-1.55) | 0.66-1.32 | 2.36 | 9.1 | +0.0 | +4.0 | +10 | -0.99/+0.52 | -0.14 dex |
| KMOS3D | 1.9-2.7 | 46 | 4.14 | 0.94 | 0.10 (0.10-0.25) | 0.10-0.18 | 3.28 | 23.2 | +9.5 | +45.2 | +71 | +0.00/+1.20 | -0.32 dex |
| KROSS | 0.6-1.1 | 161 | 0.60 | 2.42 | 2.14 (1.86-2.40) | 1.95-2.34 | 1.62 | 84.1 | +23.0 | +3.5 | -37 | -0.37/+0.29 | +0.08 dex |
| NOEMA3D (CO gas) | 1.0-1.7 | 10 | 3.54 | 1.27 | 2.04 (1.23-3.02) | 1.45-2.69 | 2.02 | 5.2 | +1.8 | -0.6 | -5 | -1.31/+0.49 | +0.00 dex |
| CRISTAL ([CII], z~5) | 4.0-6.0 | 14 | 5.52 | 1.21 | 2.57 (1.74-3.72) | 2.04-3.24 | 8.48 | 5.6 | +4.2 | +12.8 | +22 | -0.84/+0.46 | -0.24 dex |
| z4.5 curves | 4.0-5.0 | 4 | 11.35 | 1.00 | 0.49 (0.10-2.19) | 0.35-0.65 | 6.82 | 0.1 | +0.2 | +4.8 | +229 | -0.69/+0.69 | -0.36 dex |
| z4.5 assumed R_out | 4.0-5.0 | 8 | 9.36 | 1.02 | 0.10 (0.10-2.45) | 0.10-0.72 | 7.22 | 1.3 | +0.3 | +3.8 | +21 | +0.00/+1.60 | -0.30 dex |

## z >= 1.7, all samples (N = 60)

- RC100 z>1.7: N = 41, chi2(H_B) = 23.6, chi2(H_A) = 75.0
- CRISTAL: N = 14, chi2(H_B) = 9.8, chi2(H_A) = 18.5
- ALPAKA ID13: N = 1, chi2(H_B) = 0.8, chi2(H_A) = 0.0
- z4.5 curves: N = 4, chi2(H_B) = 0.2, chi2(H_A) = 4.9
- combined: r = 1.05 (0.87-1.23); chi2(H_B) = 34.5, chi2(H_A) = 98.4, chi2(MUSE) = 108.7; Delta(A - B) = +64.0 with quoted errors
- all baryonic masses -0.1 dex: r = 2.04, Delta(A - B) = -2.9
- all baryonic masses -0.2 dex: r = 3.47, Delta(A - B) = -80.3
- all baryonic masses -0.3 dex: r = 5.25, Delta(A - B) = -165.8

## Shape test on the four z ~ 4.4 curves (Roman-Oliveira et al. 2024), baryonic mass free with a 0.3 dex prior

| galaxy | z | N rings | R_max (kpc) | r (68 %) | best dlog M | chi2 min | chi2 H_B | chi2 H_A | r_A |
|---|---|---|---|---|---|---|---|---|---|
| BRI1335-0417 | 4.41 | 5 | 4.6 | 0.10 (0.10-0.37) | +0.20 | 3.5 | 6.9 | 16.5 | 6.9 |
| J081740 | 4.54 | 4 | 3.1 | 0.10 (0.10-1.95) | +0.18 | 3.1 | 3.5 | 5.9 | 7.2 |
| SGP38326-1 | 4.28 | 5 | 4.0 | 0.10 (0.10-5.62) | +0.10 | 3.4 | 3.4 | 4.6 | 6.7 |
| SGP38326-2 | 4.28 | 3 | 2.0 | 0.10 (0.10-19.95) | +0.04 | 1.2 | 1.2 | 1.6 | 6.7 |

Combined: r = 0.10 (0.10-0.35); chi2(H_B) = 15.1, chi2(H_A at z = 4.4) = 28.5. The mass is marginalised, so this is the curve shape (rise of nu with radius) speaking, not the amplitude.

## Checks

- NOEMA3D: g_obs from the authors' velocity at Re (V_c, or V_rot^2 + 3.36 sigma0^2) minus nu_obs g_bar(mass model): median -0.07 dex, 16-84 %: -0.10 to -0.03.
- CRISTAL: g_obs from the authors' velocity at Re (V_c, or V_rot^2 + 3.36 sigma0^2) minus nu_obs g_bar(mass model): median -0.08 dex, 16-84 %: -0.09 to -0.03.
- NOEMA3D inclinations: G4_38065 18 deg, G4_38232 20 deg, G4_20371 32 deg, G4_23011 51 deg, GN4_24517 41 deg, GN4_18574 12 deg, G4_24078 76 deg, GN4_32842 49 deg, G4_17555 25 deg, G4_37375 50 deg; the three discs below 21 deg carry 0.3-0.4 dex inclination errors on nu_obs.
- CRISTAL: sigma0 = 48-108 km/s against V_rot(Re) = 70-184: the pressure-support term is comparable to the rotation term, and the f_DM values have 68 % intervals of 0.3-0.5.
- ALPAKA+JWST galaxy by galaxy (outer ring, 3.8-4.2 Re):
  - GOODS-S 15503 z = 0.56, R = 19.0 kpc, i = 47 deg, g_bar/a0 = 0.15, nu_obs = 3.86 (+/- 1.40), law H_B nu = 3.08, H_A nu = 3.50
  - COSMOS 1648673 z = 1.45, R = 7.3 kpc, i = 28 deg, g_bar/a0 = 2.30, nu_obs = 2.69 (+/- 1.42), law H_B nu = 1.28, H_A nu = 1.57
  - COSMOS 3182 z = 2.10, R = 21.4 kpc, i = 24 deg, g_bar/a0 = 0.16, nu_obs = 5.60 (+/- 4.01), law H_B nu = 3.03, H_A nu = 4.93

Figure: results/a0_highz_coldgas.png. Per-galaxy values: results/a0_highz_coldgas.csv