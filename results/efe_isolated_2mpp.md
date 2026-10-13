# External field of isolated galaxies in 2M++ (2026-10-13)

Lens sample: 1214 2M++ galaxies at 15-40 Mpc, |b| > 10 deg, stellar mass 1e10-4e11 Msun (0.8 Msun/Lsun,K). e_N = stellar Newtonian field of the catalogue beyond 1 Mpc, in units of a0; 'local' = part from within 10 Mpc. Isolation after Brouwer et al. 2021: no galaxy above 10 % of the lens mass within 3 Mpc (3D, redshift distances); the threshold is floored at the completeness limit, which exceeds 10 % of the lens mass for 54 % of the lenses (low-mass bins, so their isolation is slightly too permissive).

| sample | N | median e_N (16-84 %) | median local part | fraction e_N < 5e-5 | fraction e_N < 1e-4 |
|---|---|---|---|---|---|
| all | 1214 | 1.33e-04 (9.5e-05-1.9e-04) | 4.61e-05 (1.4e-05-1.3e-04) | 0.02 | 0.20 |
| isolated, 3 Mpc | 112 | 1.14e-04 (9.1e-05-1.4e-04) | 1.21e-05 (4.7e-06-3.1e-05) | 0.00 | 0.37 |
| not isolated, 3 Mpc | 1102 | 1.37e-04 (9.6e-05-1.9e-04) | 4.99e-05 (1.8e-05-1.4e-04) | 0.02 | 0.18 |
| isolated, 5 Mpc | 21 | 9.90e-05 (8.5e-05-1.3e-04) | 5.16e-06 (2.1e-06-1.4e-05) | 0.00 | 0.52 |

| log M* bin | N all | median e_N all | N isolated (3 Mpc) | median e_N isolated | fraction isolated with e_N < 5e-5 |
|---|---|---|---|---|---|
| 10.0-10.5 | 657 | 1.33e-04 (9.5e-05-1.9e-04) | 62 | 1.16e-04 (8.9e-05-1.4e-04) | 0.00 |
| 10.5-10.8 | 324 | 1.34e-04 (9.6e-05-1.9e-04) | 25 | 1.09e-04 (8.9e-05-1.4e-04) | 0.00 |
| 10.8-11.1 | 170 | 1.29e-04 (9.6e-05-2.1e-04) | 17 | 1.14e-04 (9.5e-05-1.4e-04) | 0.00 |
| 11.1-11.5 | 61 | 1.30e-04 (8.8e-05-2.2e-04) | 8 | 1.01e-04 (9.1e-05-1.2e-04) | 0.00 |

## Reading

- Isolated galaxies have a median stellar external field of 1.1e-04 a0 against 1.4e-04 for the others and 1.5e-4 for the SPARC galaxies: isolation at 3 Mpc removes the local part (1.2e-05 against 5.0e-05) but not the field of the large-scale structure, which dominates.
- Only 0 % of the isolated galaxies have e_N < 5e-5 a0, the 2-sigma bound of the lowest-mass lensing bin, and 37 % have e_N < 1e-4. The lensing bound therefore cannot be explained by the isolation criterion.
- This e_N counts the stars of 2M++ only (0.8 Msun/Lsun,K), about 6 % of the cosmic baryons: it is a lower bound to the baryonic Newtonian field of a MOND-type theory, in which the intergalactic gas also sources gravity. A larger e_N moves r_t inward and strengthens the lensing tension; the MOND-boosted stellar field already reproduces the Milky Way's bulk acceleration (0.011 against 0.012 a0 from the CMB dipole), so the stellar value is the one a MOND universe allows.
- Caveats: the KiDS lenses lie at z = 0.1-0.5, where the large-scale field is a statistical property of the cosmic web rather than of the Local Supercluster; the local volume (Virgo, the Great Attractor) may be denser than average, in which case the lens field would be somewhat lower, but not by the factor three the bound requires unless the far field were dominated by one structure.

Figure: results/efe_isolated_2mpp.png. Per-galaxy values: results/efe_isolated_2mpp.csv