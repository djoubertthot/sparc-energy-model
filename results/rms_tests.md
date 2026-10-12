# RMS tests of the final documented model (2026-10-10)

Model: scripts/final_model.py, the formula of docs/THEORY.md with the documented inventory. Sample: 171 SPARC galaxies with at least 5 points (D512-2, NGC6789, UGC00634, UGC07232 dropped). RMS in km/s, unweighted unless noted; wmedian = median of the 1/eVobs^2-weighted RMS; p16/p84 over galaxies.

| group | test | N | median | wmedian | mean | p16 | p84 | note |
|---|---|---|---|---|---|---|---|---|
| final | budget: literal inventory, f_missing=80, no free parameter | 171 | 191.83 | 167.01 | 331.7 | 62.9 | 544.1 | rho_uniform=1.6e-22 |
| final | budget, without internal stellar part | 171 | 189.96 | 166.47 | 330.9 | 62.5 | 543.9 |  |
| final | calibrated: uniform part = rho0 = 1e-24 (THEORY.md) | 171 | 34.98 | 33.86 | 39.6 | 15.4 | 60.0 | hybrid, f_missing matched |
| final | calibrated: uniform part = rho0 = 2e-24 (MEMORY.md) | 171 | 32.07 | 32.55 | 40.6 | 14.3 | 59.0 |  |
| final | calibrated 1e-24, gamma = 1 exactly | 171 | 34.98 | 33.86 | 39.6 | 15.4 | 60.0 |  |
| final | calibrated 1e-24 + kinetic term m_k = E_k/c^2 | 171 | 34.98 | 33.86 | 39.6 | 15.4 | 60.0 |  |
| grid | grid best: beta=1.2, rho0=5e-24 | 171 | 31.88 | 31.42 | 40.4 | 13.4 | 57.8 | 25-point grid |
| accumulation | f_acc = 1 (no accumulation proxy) | 171 | 31.88 | 31.42 | 40.4 | 13.4 | 57.8 |  |
| accumulation | binary type proxy (early 1.2 / late 0.8) | 171 | 30.56 | 30.52 | 40.2 | 13.5 | 58.3 |  |
| accumulation | inverse binary (early 0.8 / late 1.2) | 171 | 33.72 | 32.53 | 41.1 | 13.6 | 59.8 |  |
| accumulation | linear age t/10 Gyr (1.1 / 0.9 / 0.65) | 171 | 31.8 | 30.56 | 39.6 | 14.1 | 58.1 |  |
| accumulation | type + (Rd/median)^alpha, best alpha=-0.25 | 171 | 30.32 | 31.88 | 39.5 | 13.5 | 57.8 | docs claim 43.8 km/s |
| accumulation | beta early/late = 1.6 / 1.2 | 171 | 31.88 | 31.15 | 47.3 | 13.4 | 63.8 |  |
| accumulation | beta early/late = 1.2 / 1.6 | 171 | 39.2 | 39.19 | 54.6 | 18.1 | 85.8 |  |
| dependencies | rho0 ~ M_bary^alpha, best alpha=+0.25 | 171 | 30.64 | 29.79 | 44.3 | 14.3 | 66.7 | 7-value scan |
| dependencies | rho0 ~ Vmax^alpha, best alpha=+0.50 | 171 | 29.64 | 30.38 | 41.1 | 13.5 | 59.2 | 7-value scan |
| dependencies | rho0 ~ surface density^alpha, best alpha=+1.00 | 171 | 29.49 | 30.54 | 57.6 | 14.2 | 73.3 | 7-value scan |
| dependencies | rho0 ~ gas fraction^alpha, best alpha=-0.50 | 171 | 30.51 | 30.51 | 40.3 | 13.5 | 60.3 | 7-value scan |
| dependencies | rho0 ~ R_disk^alpha, best alpha=-0.25 | 171 | 31.35 | 31.79 | 39.9 | 14.1 | 57.8 | 7-value scan |
| subsamples | Q <= 2 | 159 | 33.41 | 32.4 | 41.8 | 14.2 | 58.7 |  |
| subsamples | 40 < i < 75 | 100 | 28.74 | 30.23 | 41.4 | 13.3 | 58.6 |  |
| subsamples | Q <= 2 and 40 < i < 75 | 93 | 29.94 | 31.08 | 43.0 | 13.9 | 62.5 |  |
| subsamples | Q = 1 and 40 < i < 75 | 57 | 38.28 | 36.95 | 49.5 | 19.1 | 69.9 |  |
| subsamples | i > 75 only | 36 | 38.4 | 39.32 | 41.1 | 17.4 | 55.7 |  |
| subsamples | first 100 rotmod files | 100 | 35.3 | 34.92 | 43.9 | 16.1 | 56.3 |  |
| uncertainty | MC inclination +/- e_Inc, 100 draws (p16/p84 over draws) | 171 | 32.52 | nan | 32.6 | 31.9 | 33.6 |  |
| reference | baryons only, Y* = 0.5, no dark component | 171 | 43.04 | 40.97 | 44.6 | 21.3 | 65.5 |  |
| reference | MOND simple mu, a0 = 1.2e-10, Y* = 0.5, 0 free parameters | 171 | 13.66 | 13.06 | 17.3 | 6.7 | 28.2 |  |

## Grid beta x rho0, median RMS, calibrated normalisation

|   beta |   1e-25 |   5e-25 |   1e-24 |   2e-24 |   5e-24 |
|-------:|--------:|--------:|--------:|--------:|--------:|
|    1   |   42.97 |   40.08 |   36.76 |   36    |   32.83 |
|    1.2 |   42.63 |   36.95 |   36.32 |   34.29 |   31.88 |
|    1.4 |   41.32 |   36.73 |   34.98 |   32.07 |   35.85 |
|    1.6 |   38.53 |   35.5  |   33.74 |   36.22 |   40.71 |
|    1.8 |   37.61 |   34.53 |   37.16 |   41.74 |   48.67 |
