# External-field effect from the QUMOND field equation versus the 1D formula (2026-10-12)

Point mass M in a uniform Newtonian external field e_N (units of a0); effective enclosed mass from Gauss's law on the QUMOND field; V_eff(r)/V_flat with V_flat = (G M a0)^(1/4), the quantity the weak-lensing stacks measure (Mistele et al. 2024: flat to 5-8 % from 50 to 300 kpc, no decline seen to ~1 Mpc). The 1D column is the Chae et al. (2020) formula g = nu_e(g_N/a0, e) g_N with the MOND-boosted e = nu(e_N) e_N, as used in results/lensing_vflat.md.

## M_b = 1e+10 Msun, V_flat = 112 km/s

| e_N / a0 | boosted e / a0 | r where g_N = e_N (kpc) | V(50)/V_flat | V(100)/V_flat | V(200)/V_flat | V(300)/V_flat | V(500)/V_flat | V(1000)/V_flat | V(2000)/V_flat | 1D formula at 300 kpc | 1D at 1000 kpc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0e+00 | 0.0000 | inf | 1.02 | 1.01 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4.0e-05 | 0.0063 | 539 | 1.02 | 1.01 | 1.00 | 0.99 | 0.93 | 0.67 | 0.47 | 0.28 | 0.15 |
| 1.5e-04 | 0.0123 | 278 | 1.02 | 1.01 | 0.98 | 0.87 | 0.68 | 0.48 | 0.34 | 0.24 | 0.13 |
| 5.0e-04 | 0.0226 | 152 | 1.02 | 0.99 | 0.80 | 0.65 | 0.51 | 0.36 | 0.25 | 0.21 | 0.11 |
| 2.0e-03 | 0.0457 | 76 | 1.00 | 0.80 | 0.57 | 0.47 | 0.36 | 0.26 | 0.18 | 0.18 | 0.10 |
| 1.0e-02 | 0.1051 | 34 | 0.77 | 0.55 | 0.39 | 0.32 | 0.25 | 0.17 | 0.12 | 0.15 | 0.08 |

## M_b = 5e+10 Msun, V_flat = 168 km/s

| e_N / a0 | boosted e / a0 | r where g_N = e_N (kpc) | V(50)/V_flat | V(100)/V_flat | V(200)/V_flat | V(300)/V_flat | V(500)/V_flat | V(1000)/V_flat | V(2000)/V_flat | 1D formula at 300 kpc | 1D at 1000 kpc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0e+00 | 0.0000 | inf | 1.04 | 1.02 | 1.01 | 1.01 | 1.00 | 1.00 | 1.00 | 1.01 | 1.00 |
| 4.0e-05 | 0.0063 | 1205 | 1.04 | 1.02 | 1.01 | 1.01 | 1.00 | 0.95 | 0.71 | 0.41 | 0.23 |
| 1.5e-04 | 0.0123 | 622 | 1.04 | 1.02 | 1.01 | 1.00 | 0.96 | 0.72 | 0.51 | 0.35 | 0.20 |
| 5.0e-04 | 0.0226 | 341 | 1.04 | 1.02 | 1.00 | 0.94 | 0.76 | 0.54 | 0.38 | 0.31 | 0.17 |
| 2.0e-03 | 0.0457 | 170 | 1.04 | 1.01 | 0.85 | 0.70 | 0.54 | 0.38 | 0.27 | 0.27 | 0.15 |
| 1.0e-02 | 0.1051 | 76 | 1.02 | 0.82 | 0.58 | 0.47 | 0.37 | 0.26 | 0.18 | 0.23 | 0.13 |

## M_b = 2e+11 Msun, V_flat = 238 km/s

| e_N / a0 | boosted e / a0 | r where g_N = e_N (kpc) | V(50)/V_flat | V(100)/V_flat | V(200)/V_flat | V(300)/V_flat | V(500)/V_flat | V(1000)/V_flat | V(2000)/V_flat | 1D formula at 300 kpc | 1D at 1000 kpc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0e+00 | 0.0000 | inf | 1.08 | 1.04 | 1.02 | 1.01 | 1.01 | 1.00 | 1.00 | 1.01 | 1.00 |
| 4.0e-05 | 0.0063 | 2410 | 1.08 | 1.04 | 1.02 | 1.01 | 1.01 | 1.00 | 0.95 | 0.56 | 0.32 |
| 1.5e-04 | 0.0123 | 1244 | 1.08 | 1.04 | 1.02 | 1.01 | 1.00 | 0.96 | 0.72 | 0.49 | 0.28 |
| 5.0e-04 | 0.0226 | 682 | 1.08 | 1.04 | 1.02 | 1.01 | 0.98 | 0.76 | 0.54 | 0.44 | 0.24 |
| 2.0e-03 | 0.0457 | 341 | 1.08 | 1.04 | 1.01 | 0.95 | 0.76 | 0.54 | 0.38 | 0.38 | 0.21 |
| 1.0e-02 | 0.1051 | 152 | 1.08 | 1.02 | 0.82 | 0.67 | 0.52 | 0.37 | 0.26 | 0.32 | 0.18 |

## Reading

- In the field equation the external field only matters where the galaxy's own Newtonian field falls below e_N, i.e. beyond r_t = sqrt(G M / e_N): 300-700 kpc for e_N = 1.5e-4 a0 (the 2M++ field of SPARC galaxies) and M_b = 1e10 to 2e11. Inside r_t the effective circular velocity stays within a few per cent of the isolated value.
- The 1D formula applied in lensing_vflat.py puts the decline two to ten times closer in, because it adds the boosted external field linearly to the internal one at every point; the bound e_N < 4e-5 a0 derived from it is a bound on that formula, not on QUMOND.
- What the lensing flatness tests in a MOND field theory is therefore the regime beyond r_t: the 500 kpc to 1 Mpc points of Mistele et al. for the massive bins, where QUMOND with e_N = 1.5e-4 predicts a decline of the size given in the tables.
- Consequence for the theory landscape: a QUMOND-type theory, and its relativistic version RMOND, is not excluded by the flatness at 50-300 kpc; it is tested by the flatness at 0.5-1 Mpc and by the external-field dependence of the discs, which SPARC cannot resolve (0.24 dex of dynamic range in e_N).