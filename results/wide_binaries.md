# Wide binaries below a0 (2026-10-13)

El-Badry et al. 2021 pairs: 132829 main-sequence pairs within 250 pc with R < 0.1 downloaded; 33286 after the quality cuts (both parallaxes > 4 mas at > 20 sigma, RUWE < 1.4, chance alignment < 1 %, parallaxes consistent at 3 sigma, 2 < M_G < 14, both stars within 0.5 mag of the main sequence in M_G(BP-RP), separation < 200 arcsec, no multi-peak detections, error on v~ < 0.2). Radial velocity available for the perspective correction in 46 % of the pairs.

Calibration (g_N > 10 a0, 24204 pairs): median v~ = 0.563. Newtonian Monte Carlo of the median v~: uniform eccentricities 0.608, thermal 0.554, superthermal (p(e) ~ e^1.5) 0.531.

| g_N / a0 | N | median s (kAU) | median M_tot | median v~ | ratio to calibration | Newton | law, no EFE | QUMOND, e_N = 1.2 | fraction v~ > 1 | fraction v~ > 1.5 | median error on v~ | pairs with two agreeing RVs: N, ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30.0-1000.0 | 17155 | 0.69 | 0.97 | 0.560 ± 0.003 | 0.995 ± 0.006 | 1.000 | 1.000 | 1.000 | 0.113 | 0.019 | 0.030 | 1189, 0.927 ± 0.020 |
| 10.0-30.0 | 6455 | 1.70 | 1.05 | 0.572 ± 0.006 | 1.015 ± 0.010 | 1.000 | 1.007 | 1.007 | 0.137 | 0.030 | 0.059 | 374, 0.996 ± 0.035 |
| 3.0-10.0 | 4652 | 3.06 | 1.10 | 0.569 ± 0.006 | 1.011 ± 0.011 | 1.000 | 1.048 | 1.047 | 0.135 | 0.039 | 0.069 | 299, 0.958 ± 0.050 |
| 1.0-3.0 | 2378 | 5.56 | 1.16 | 0.589 ± 0.010 | 1.046 ± 0.019 | 1.000 | 1.158 | 1.137 | 0.124 | 0.026 | 0.079 | 164, 1.006 ± 0.086 |
| 0.3-1.0 | 1380 | 9.98 | 1.21 | 0.595 ± 0.012 | 1.057 ± 0.022 | 1.000 | 1.356 | 1.167 | 0.135 | 0.034 | 0.089 | 118, 1.055 ± 0.064 |
| 0.1-0.3 | 554 | 17.55 | 1.23 | 0.594 ± 0.023 | 1.055 ± 0.041 | 1.000 | 1.663 | 1.168 | 0.148 | 0.029 | 0.112 | 56, 1.013 ± 0.112 |
| 0.03-0.1 | 118 | 29.18 | 1.38 | 0.730 ± 0.034 | 1.297 ± 0.061 | 1.000 | 2.018 | 1.168 | 0.220 | 0.042 | 0.129 | 11, 1.360 ± 0.065 |

chi2 of the ratio on the bins below 1 a0 (3 bins): Newton 32.4, law without EFE 549.7, QUMOND with the Galactic field 37.5.
Lowest bin (0.03-0.1 a0, 118.0 pairs, median s = 29.2 kAU): ratio 1.297 ± 0.061; the law without EFE predicts 2.02, i.e. 12 sigma away; QUMOND with EFE 1.17 (-2.1 sigma); Newton 1 (-4.9 sigma).

## Reading

- From 1000 to 0.1 a0 the ratio stays within 0.99-1.06 of the calibration. The law applied locally without external-field effect, the only form in which a modified-dynamics reading survived the lensing test, predicts 1.36 at 0.3-1 a0 and 1.66 at 0.1-0.3 a0 against 1.057 ± 0.022 and 1.055 ± 0.041 measured: excluded at 14 and 15 sigma, far beyond the 10-20 % that unresolved companions and the eccentricity distribution can move the median (Monte Carlo above; Chae 2023 versus Banik et al. 2024).
- QUMOND with the Galactic external field predicts 1.17 in the same two bins, 5.1 and 2.7 sigma above the data; below 0.1 a0 the 118 widest pairs give 1.30 ± 0.06, above every prediction but the law's, with a fraction of v~ > 1 that rises from 0.13 to 0.22: the regime of unbound and contaminated pairs. The known systematics are of the size of the QUMOND effect, so this quick analysis sides with Banik et al. 2024 (Newton) rather than Chae 2023 (a boost of 1.3-1.5 in g) without settling their dispute.
- The fraction of pairs with v~ > 1.5 stays at 2-4 % in every bin, so chance alignments do not drive the result; the subsample with two agreeing radial velocities gives the same ratios where it is large enough. The dark-component reading predicts Newton at every separation and is contradicted by nothing here.
- Caveats: hidden companions add mass and photocentre motion, which raise v~ at wide separations more than at close ones (an s-dependent contamination that mimics a MOND-like rise, estimated at up to 20 % by Banik et al.); the eccentricity distribution of wide binaries steepens with s (Hwang et al. 2022), which changes the Newtonian median by the amount the Monte Carlo gives; reddening is neglected within 250 pc; masses from a single main-sequence relation (0.1-0.2 Msun errors, 5-10 % in v_c).

Figure: results/wide_binaries.png. Per-pair values: results/wide_binaries.csv