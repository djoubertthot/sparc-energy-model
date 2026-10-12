# The interpolation function and the universality of a0 (2026-10-10)

Binned relation: 16 bins, 2795 points, bootstrap errors on bin medians (floor 0.01 dex).

| function | best a0 (1e-10) | chi2 / bins | rms of bin residuals (dex) |
|---|---|---|---|
| RAR | 1.45 | 17.70 | 0.029 |
| simple | 1.45 | 24.18 | 0.034 |
| standard | 1.45 | 9.58 | 0.048 |
| n-family, n = 1.4 | 1.45 | 5.54 | 0.038 |
| Verlinde, a0 = cH0 fixed | (cH0/6 = 1.13) | 73.51 | 0.066 |

Deep-regime slope d log g_obs / d log g_bar for g_bar < 0.1 a0 (5 bins): 0.498 (MOND: 0.500; Newton + baryon-tracing mass: 1.000).

Per-galaxy a0 (Y free with prior), 163 galaxies off the grid edges: median 0.89, scatter 0.33 dex, 16-84 %: 0.50-1.78.
- Spearman(log a0, Vmax) = +0.22, p = 0.006
- Spearman(log a0, SB) = +0.12, p = 0.126
- Spearman(log a0, T) = -0.15, p = 0.061
- Spearman(log a0, fgas) = -0.05, p = 0.541

Figure: results/measure_nu.png

## Caveat on the binned a0

The a0 of 1.45 returned by every function in the binned fit is an artefact of binning: medians of log g_obs in bins
of noisy log g_bar are biased upward where the relation has slope 1/2 (Eddington bias), and the error floor of 0.01 dex
lets the many deep-regime bins dominate. The point-wise deep regime (g_bar < 0.1 a0) gives a median residual of
-0.008 dex with SED M/L (implied a0 = 1.16e-10) and -0.021 dex with fixed M/L (1.09e-10), consistent with
results/fit_a0.md. The binned fit is retained only for the comparison of the SHAPE of nu between functions.
