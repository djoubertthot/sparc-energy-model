# The de Sitter-vacuum interpolation function on SPARC (2026-10-12)

Binned relation of scripts/measure_nu.py: 16 bins, 2795 points (Y36 per galaxy, errors < 10 %). The binned a0 is biased high by binning (see results/measure_nu.md); only the SHAPE comparison between functions is meaningful here. The point-wise columns use a0 = 1.2e-10 on all points.

| function | best a0 binned (1e-10) | chi2 / bins | rms of bin residuals (dex) | point-wise scatter at a0 = 1.2e-10 (dex) | point-wise median residual (dex) |
|---|---|---|---|---|---|
| RAR (McGaugh 2016) | 1.45 | 17.70 | 0.029 | 0.148 | -0.012 |
| simple | 1.45 | 24.18 | 0.034 | 0.149 | -0.019 |
| standard | 1.45 | 9.58 | 0.048 | 0.151 | +0.050 |
| n-family, n = 1.4 | 1.45 | 5.54 | 0.038 | 0.149 | +0.029 |
| de Sitter vacuum (Milgrom 1999): sqrt(1 + 1/y) | 1.45 | 8.65 | 0.044 | 0.151 | +0.027 |

High-acceleration tail (g_bar > 3 a0), where the functions differ: nu - 1 predicted at g_bar = 3, 10 and 30 a0, and the measured median residual log(g_obs/g_bar) of the SPARC points in those ranges:
- RAR (McGaugh 2016): nu - 1 = 0.215, 0.044, 0.004 -> log nu = 0.085, 0.019, 0.002
- simple: nu - 1 = 0.264, 0.092, 0.032 -> log nu = 0.102, 0.038, 0.014
- standard: nu - 1 = 0.049, 0.005, 0.001 -> log nu = 0.021, 0.002, 0.000
- n-family, n = 1.4: nu - 1 = 0.092, 0.005, 0.000 -> log nu = 0.038, 0.002, 0.000
- de Sitter vacuum (Milgrom 1999): sqrt(1 + 1/y): nu - 1 = 0.155, 0.049, 0.017 -> log nu = 0.062, 0.021, 0.007
- measured, g_bar/a0 in [2, 3]: median log(g_obs/g_bar) = +0.086 dex (N = 124)
- measured, g_bar/a0 in [3, 10]: median log(g_obs/g_bar) = +0.020 dex (N = 221)
- measured, g_bar/a0 in [10, 32]: median log(g_obs/g_bar) = -0.127 dex (N = 75)

## The scale

- a_Lambda = c sqrt(Lambda/3) = c H0 sqrt(Omega_Lambda) = 5.69e-10 m/s^2 (H0 = 70, Omega_Lambda = 0.7).
- Milgrom 1999 vacuum inertia gives a0 = 2 a_Lambda = 1.14e-09, i.e. 9.5 times the measured 1.2e-10: the derivation fixes the shape sqrt(1 + 1/y) and the order of magnitude, not the prefactor; a0 = 0.21 a_Lambda is measured, 2 a_Lambda predicted.
- With the prefactor left free, the Unruh temperature at a0 is T_U(a0)/T_dS = a0/a_Lambda = 0.21: the MOND regime starts where the Unruh temperature of the orbit falls to a fifth of the de Sitter temperature.

## Reading

- The vacuum-inertia function has no free parameter besides a0 and is a closed form; its fit quality relative to the RAR and n = 1.4 functions is given in the table. It is the only candidate derivation of the law that touches quantum physics (Unruh radiation in de Sitter space), and it is a heuristic: no field theory, no cosmology, no relativistic limit.