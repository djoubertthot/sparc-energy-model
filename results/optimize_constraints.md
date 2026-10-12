# Constrained optimisation of the tag-1.0 model family (2026-10-10)

Median over 171 SPARC galaxies. RMS, bias, outer (R>20 kpc), inner (R<2 Rd) in km/s. n_H = median mean hydrogen density of the halo inside R_last; bound = 3 x Milky-Way hot-halo beta model. frac>bound = fraction of galaxies violating it.

| level | parameters | RMS | bias | outer | inner | early RMS | early inner | late RMS | n_H (cm^-3) | frac>bound | max n/bound |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L0 RMS only, mode last | ystar=0.52, gfac=2.19, M=3.04, beta=0.74, rc=1.42 | 14.97 | -2.67 | 8.5 | -5.23 | 36.26 | -26.85 | 12.7 | 6.0e-02 | 1.0 | 6.4e+02 |
| L1 + bias~0, outer~0, mode last | ystar=0.95, gfac=3.00, M=3.55, beta=0.55, rc=56.12 | 18.64 | -1.0 | -0.81 | -0.83 | 45.82 | -31.51 | 15.35 | 7.0e-02 | 1.0 | 7.5e+02 |
| L2a + no hidden parameter, mode fixed R0 | ystar=1.42, gfac=3.00, M=10.67, beta=0.50, rc=5.28, R0=142.79 | 21.27 | -1.0 | -2.0 | -5.8 | 67.46 | -76.64 | 16.89 | 0.0e+00 | 0.6 | 1.4e+02 |
| L2b + no hidden parameter, self-similar (R_ref, r_c prop. to Rd) | ystar=1.16, gfac=3.00, M=6.21, beta=1.35, k=19.32, k2=10.73 | 18.28 | 0.83 | -1.86 | -2.01 | 47.38 | -46.81 | 13.51 | 1.0e-02 | 0.92 | 1.5e+02 |
| L3a + Y*<=0.7, gas<=1.3, mode fixed | ystar=0.40, gfac=1.30, M=3.26, beta=1.42, rc=3.69, R0=98.09 | 22.48 | 1.0 | -1.89 | -4.59 | 62.28 | -76.96 | 17.37 | 4.0e-02 | 0.96 | 4.7e+02 |
| L3b + Y*<=0.7, gas<=1.3, self-similar | ystar=0.48, gfac=1.04, M=3.93, beta=1.34, k=25.11, k2=1.70 | 18.79 | 0.56 | 1.69 | -3.63 | 36.24 | -13.76 | 15.59 | 5.0e-02 | 1.0 | 4.7e+02 |
| L4a + hot-gas density bound, mode fixed | ystar=0.42, gfac=1.30, M=3.23, beta=1.40, rc=3.62, R0=103.30 | 22.36 | 0.97 | -1.72 | -4.41 | 61.71 | -78.07 | 17.71 | 4.0e-02 | 0.96 | 4.6e+02 |
| L4b + hot-gas density bound, self-similar | ystar=0.56, gfac=1.30, M=3.80, beta=1.33, k=26.05, k2=1.73 | 17.85 | -0.1 | 1.43 | -3.72 | 33.91 | -19.98 | 14.95 | 5.0e-02 | 0.99 | 4.4e+02 |
| L4c + hot-gas density bound, mode last (hidden parameter allowed) | ystar=0.36, gfac=1.25, M=3.88, beta=1.01, rc=4.41 | 18.76 | 0.32 | -1.99 | 1.75 | 38.55 | -11.94 | 13.33 | 8.0e-02 | 1.0 | 8.2e+02 |
| L5 no halo (M=0), Y*<=0.7, gas<=1.3 | ystar=0.70, gfac=1.30, M=0.00, beta=0.70 | 34.26 | 25.45 | 67.12 | 10.41 | 46.45 | 5.66 | 27.83 | 0.0e+00 | 0.0 | 0.0e+00 |
