# RMS tests of the tag-1.0 baryonic model (2026-10-10)

Model: scripts/final_baryonic_model.py, Y*=0.88, gas x2.7, CGM beta profile M=3.8, r_c=6 kpc, beta=0.68, E0=0 (archive_v1/docs/CONSOLIDATED_NOTES.md). Sample: 171 SPARC galaxies with at least 5 points. RMS in km/s; bias = median over galaxies of the mean residual Vobs-Vmodel; outer_R20 = same for points beyond 20 kpc; inner_2Rd = same for points inside 2 R_disk; claimed = value quoted in the docs.

| group | test | N | median | wmedian | mean | p16 | p84 | bias | outer_R20 | inner_2Rd | claimed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| final | final 1.0, M_CGM(<R_last) = 3.8 x (0.5 L36 + 1.33 MHI) | 171 | 18.9 | 17.6 | 29.0 | 8.9 | 48.3 | -5.0 | -12.5 | -2.7 | 21 |
| final | final 1.0, M_CGM(<R_last) = 3.8 x (0.88 L36 + 2.7 x 1.33 MHI) | 171 | 28.6 | 28.4 | 46.6 | 11.6 | 85.2 | -23.3 | -61.1 | -10.0 | 21 |
| final | final 1.0, M_CGM(<25 kpc) = 3.8 x (0.5 L36 + 1.33 MHI) | 171 | 20.1 | 19.7 | 34.9 | 8.2 | 67.7 | 2.8 | -56.8 | 0.1 | 21 |
| options | no CGM, Y*=0.5, gas x1.0 | 171 | 43.0 | 41.0 | 44.6 | 21.3 | 65.5 | 37.9 | 82.4 | 20.8 | 56 |
| options | no CGM, Y*=1.0, gas x1.0 | 171 | 33.1 | 31.1 | 34.0 | 15.0 | 51.9 | 18.6 | 55.9 | 3.7 | 32 / 31.7 |
| options | no CGM, Y*=1.0, gas x2.0 | 171 | 29.2 | 26.3 | 30.6 | 11.2 | 49.4 | 15.5 | 49.2 | 2.1 | 28 |
| options | no CGM, Y*=1.0, gas x2.5 | 171 | 27.0 | 25.5 | 29.2 | 9.7 | 47.3 | 12.7 | 45.1 | 1.4 | 26 |
| options | no CGM, Y*=1.2, gas x2.0 | 171 | 27.9 | 26.1 | 31.9 | 10.6 | 52.6 | 7.9 | 39.2 | -1.0 | 25 |
| options | no CGM, Y*=1.5, gas x2.0 | 171 | 25.3 | 24.6 | 36.0 | 10.4 | 67.3 | 1.9 | 27.0 | -5.1 | 26 |
| options | no CGM, Y*=0.88, gas x2.7 | 171 | 26.9 | 25.1 | 28.5 | 9.7 | 46.5 | 15.7 | 46.8 | 3.7 |  |
| options | no CGM, Y*=0.5 disk / 0.7 bulge, gas x1 (SPARC standard) | 171 | 41.1 | 39.7 | 43.5 | 21.3 | 62.1 | 35.0 | 81.8 | 18.8 |  |
| versions | Y*=1.0, gas x2, CGM 3.5, r_c=20, beta=0.6 | 171 | 20.0 | 19.2 | 28.3 | 9.9 | 47.2 | -1.6 | 3.6 | -1.3 | 26 |
| versions | Y*=1.0, gas x2, CGM 3.5, r_c=20, beta=0.7 (REVERT) | 171 | 20.0 | 19.1 | 28.4 | 9.9 | 45.9 | -1.7 | 2.6 | -1.4 | ~26, outer +30 |
| versions | Y*=1.0, gas x2, M=4.0, r_c=10, beta=0.75 (FINAL_PHYSICAL) | 171 | 20.5 | 19.8 | 31.0 | 9.0 | 52.3 | -4.2 | -17.3 | -2.8 | 23 |
| versions | Y*=1.0, gas x2, CGM 4.5, r_c=30, beta=0.5 (FINAL_ADJUSTMENTS) | 171 | 21.1 | 19.6 | 30.2 | 10.8 | 48.4 | -3.6 | -4.8 | -1.9 | 24, outer +7 |
| versions | Y*=0.9, gas 2.5, M=4.0, r_c=6, beta=0.65 (EVEN_FINER) | 171 | 19.2 | 17.6 | 29.7 | 9.0 | 49.5 | -5.8 | -14.7 | -3.0 | ~22 |
| grid | grid best CGM: M=3.0, r_c=2.0, beta=0.68 | 171 | 17.3 | 15.9 | 29.2 | 7.6 | 54.8 | -7.1 | -1.8 | -8.1 |  |
| grid | Y* x gas best at final CGM: Y*=0.7, gas x2.0 | 171 | 17.8 | 17.3 | 25.8 | 9.3 | 39.5 | 0.1 | 1.2 | 2.1 |  |
| by type | final, early types T<=4 | 46 | 35.9 | 32.7 | 54.6 | 21.4 | 101.3 | -30.6 | -17.2 | -32.3 | 47 (Y*=1, gas x2) |
| by type | final, intermediate 5<=T<=7 | 48 | 17.6 | 16.3 | 25.6 | 9.6 | 41.3 | -0.6 | -25.4 | 3.5 | 31 |
| by type | final, late types T>=8 | 77 | 12.8 | 12.0 | 15.8 | 7.1 | 26.8 | 0.5 | 10.9 | 0.7 | 22 |
| subsample | final, Q<=2 and 40<i<75 | 93 | 15.3 | 14.9 | 26.7 | 8.2 | 41.2 | -3.1 | -8.4 | -1.7 |  |
| reference | MOND simple mu, a0=1.2e-10, Y*=0.5/0.7, 0 free parameters | 171 | 14.0 | 12.4 | 17.3 | 6.8 | 29.5 | -4.2 | -10.8 | -3.5 | ~15 (archive_v1/docs/MODEL_COMPARISON.md) |

## Y* x gas scan at the final CGM, median RMS

|   ystar |   1.0 |   1.5 |   2.0 |   2.7 |
|--------:|------:|------:|------:|------:|
|    0.5  |  19.8 |  18.8 |  18.6 |  18.4 |
|    0.7  |  19.6 |  18.5 |  17.8 |  17.8 |
|    0.88 |  19.7 |  19.3 |  18.6 |  18.9 |
|    1    |  21   |  20.9 |  20.5 |  19.3 |
|    1.2  |  21   |  20.6 |  19.8 |  19.2 |

## CGM term implications (final, last_std)

median M_cgm(<R_last) = 1.98e+10 Msun; median mean hydrogen density inside R_last = 7.50e-02 cm^-3; median electron dispersion measure through R_last = 1195 pc cm^-3.
