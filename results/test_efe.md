# External-field-effect test on SPARC (Chae et al. 2020 method), 2026-10-10

171 galaxies fitted; 152 with coordinates and an environmental field from 2M++.
Calibration at the Milky Way: baryonic Newtonian e_N = 1.26e-04, MOND-boosted 0.011 (expected ~0.012), direction RA 171, Dec -32.
e_N,env (Newtonian, baryons only): median 1.47e-04, 16-84 % 1.16e-04 to 2.04e-04 (0.24 dex of dynamic range); MOND-boosted nu(e)e: median 0.012; local part < 10 Mpc median 3.65e-05.
e_N,fit: median 0.002; galaxies preferring e_N > 0 at Delta chi2 > 1: 56%, at > 4: 42%; galaxies with e_N = 0 excluded at 1 sigma (lower limit > 0): 56%.
Spearman(e_N,fit, e_N,env), all 152: rho = +0.14, p = 0.084
Spearman(e_N,fit, e_N,env), galaxies where the curve constrains e_N (Delta chi2 > 1), N = 83: rho = +0.02, p = 0.823
Spearman(low-acceleration RAR residual, e_N,env), N = 113: rho = +0.00, p = 0.959  (MOND predicts negative)
Spearman(outer slope d ln V / d ln R, e_N,env), N = 152: rho = -0.07, p = 0.381  (MOND predicts negative)
Chae-style quartile test: median e_N,fit = 0.003 in the top quartile of e_N,env vs 0.000 in the bottom quartile; Mann-Whitney one-sided p = 0.033
Galaxies with declining outer curves (slope < -0.1): 38; their e_N,env rank (0 = weakest field, 1 = strongest), median 0.49 (MOND predicts > 0.5)

| e_N,env tercile | N | median e_N,env | median e_N,fit | median Delta chi2 | median low-g residual (dex) | median outer slope |
|---|---|---|---|---|---|---|
| low | 51 | 1.16e-04 | 0.000 | 0.0 | -0.041 | +0.13 |
| mid | 50 | 1.47e-04 | 0.002 | 3.9 | -0.043 | +0.01 |
| high | 51 | 2.03e-04 | 0.002 | 2.5 | -0.035 | +0.08 |

Galaxies with the strongest fitted external field (Delta chi2 > 4):
| Galaxy | e_N,fit [1 sigma] | e_N,env | Delta chi2 | outer slope |
|---|---|---|---|---|
| CamB | 1.000 [0.501, 1.000] | 1.34e-04 | 318.3 | +0.82 |
| F574-2 | 1.000 [0.050, 1.000] | nan | 6.7 | +0.57 |
| PGC51017 | 1.000 [0.631, 1.000] | 1.92e-04 | 327.7 | -0.06 |
| F579-V1 | 0.501 [0.063, 0.794] | nan | 11.4 | +0.00 |
| F563-V1 | 0.316 [0.200, 0.794] | nan | 16.1 | -0.19 |
| UGC07577 | 0.316 [0.200, 0.631] | 1.32e-04 | 80.9 | +0.94 |
| NGC3953 | 0.251 [0.200, 1.000] | 2.00e-04 | 13.1 | -0.14 |
| NGC5371 | 0.158 [0.158, 0.158] | 2.71e-04 | 333.3 | +0.14 |
| UGC04305 | 0.158 [0.158, 1.000] | 1.33e-04 | 15.8 | -0.68 |
| NGC2976 | 0.158 [0.100, 0.398] | 1.36e-04 | 15.5 | +0.11 |
| NGC4051 | 0.126 [0.100, 0.501] | 2.04e-04 | 5.5 | -0.16 |
| UGC01230 | 0.079 [0.063, 0.158] | 6.86e-05 | 19.0 | -0.06 |

Figure: results/efe.png. Per-galaxy values: results/efe_per_galaxy.csv
