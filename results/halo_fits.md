# Halo fits on SPARC: NFW, Burkert and the acceleration law (2026-10-10)

171 galaxies. Y_disk prior: 0.1 dex around the SED value (or 0.5). Haloes: 3 parameters each; law: 1 parameter.

| model | median RMS (km/s) | median chi2/N | sum AIC | galaxies with lowest AIC | median Y_disk |
|---|---|---|---|---|---|
| NFW + baryons | 4.5 | 0.88 | 11340 | 30% | 0.46 |
| Burkert + baryons | 3.7 | 0.59 | 10439 | 63% | 0.48 |
| g = nu(g_bar/a0) g_bar | 8.9 | 3.55 | 42602 | 8% | 0.43 |

Burkert vs NFW: Burkert better in 69% of galaxies; median Delta chi2 (Burkert - NFW) = -2.1; late types (T>=8): 75% prefer Burkert, early (T<=4): 59%.
NFW concentration-mass relation, 144 galaxies with 1e9 < M200 < 1e14: median offset from Dutton & Maccio 2014 -0.05 dex, scatter 0.36 dex (LambdaCDM: 0.11 dex); fraction with c200 below the relation by > 0.3 dex: 27%.
Burkert core surface density rho_0 r_0: median 169 Msun/pc^2, 16-84 % 79-377 (Donato et al. 2009: ~140, scatter 0.2 dex); Spearman(rho_0 r_0, V_max) = +0.62.
Burkert core radius / disc scale length: median 1.98, 16-84 % 0.97-3.76; Spearman(r_0, Rd) = +0.62.
Stellar M/L pulled by the fits: NFW median Y/Y0 = 0.97, Burkert 0.99, law 0.86.

Figure: results/halo_fits.png; per galaxy: results/halo_fits_per_galaxy.csv
