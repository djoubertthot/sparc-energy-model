# The law in galaxy clusters: the missing factor (2026-10-12)

g_RAR = nu(g_bar/a0) g_bar with the galaxy a0 = 1.2e-10 m/s^2; residual = log(g_tot/g_RAR). Two ways to express the missing mass: M_extra/M_bar = (g_tot - g_RAR)/g_bar is the Newtonian mass to add to the baryons so that the law, applied to the baryons only, plus that mass in Newtonian gravity, gives g_tot (a dark-component reading); M_MOND/M_bar is the total mass the law itself would need, in units of the observed baryons (the usual 'MOND missing mass' of the cluster literature, the quantity behind the 'factor 2').

## CLASH (Tian et al. 2020), 20 clusters, 14-600 kpc

- 10-30 kpc: N = 19, median residual +0.43 dex (16-84 %: +0.26 to +0.53), median g_bar/a0 = 0.82, M_extra/M_bar = 3.2, M_MOND/M_bar = 3.9 (2.2-5.1)
- 30-100 kpc: N = 1, median residual +0.53 dex (16-84 %: +0.53 to +0.53), median g_bar/a0 = 0.52, M_extra/M_bar = 4.7, M_MOND/M_bar = 5.4 (5.4-5.4)
- 100-300 kpc: N = 38, median residual +0.57 dex (16-84 %: +0.46 to +0.64), median g_bar/a0 = 0.30, M_extra/M_bar = 6.7, M_MOND/M_bar = 7.0 (4.5-8.4)
- 300-700 kpc: N = 26, median residual +0.47 dex (16-84 %: +0.35 to +0.55), median g_bar/a0 = 0.18, M_extra/M_bar = 5.4, M_MOND/M_bar = 5.3 (3.5-7.0)
- all points: N = 84, median residual +0.50 dex (16-84 %: +0.38 to +0.61), median g_bar/a0 = 0.28, M_extra/M_bar = 5.3, M_MOND/M_bar = 5.4 (3.4-7.9)

## X-COP (Eckert et al. 2019, 2022), 11 clusters, R2500 to R200

- R2500 (median 546 kpc): N = 11, median residual +0.35 dex (16-84 %: +0.31 to +0.40), median g_bar/a0 = 0.10, M_extra/M_bar = 4.4, M_MOND/M_bar = 3.7 (3.3-4.6)
- R500 (median 1214 kpc): N = 11, median residual +0.15 dex (16-84 %: +0.14 to +0.18), median g_bar/a0 = 0.07, M_extra/M_bar = 1.9, M_MOND/M_bar = 1.8 (1.7-2.0)
- R200 (median 1938 kpc): N = 11, median residual +0.06 dex (16-84 %: +0.03 to +0.08), median g_bar/a0 = 0.04, M_extra/M_bar = 0.7, M_MOND/M_bar = 1.3 (1.1-1.4)

## Scale of the relation in clusters

- CLASH: best a0 with the RAR shape = 1.76e-09 m/s^2 = 14.7 x the galaxy value, scatter 0.11 dex (Tian et al. 2020: 2.0e-9 with a power law of slope 0.5)
- X-COP: best a0 with the RAR shape = 3.15e-10 m/s^2 = 2.6 x the galaxy value, scatter 0.13 dex (Tian et al. 2020: 2.0e-9 with a power law of slope 0.5)
- both: best a0 with the RAR shape = 1.05e-09 m/s^2 = 8.8 x the galaxy value, scatter 0.20 dex (Tian et al. 2020: 2.0e-9 with a power law of slope 0.5)

## Reading

- With the galaxy a0 the law under-predicts the cluster accelerations everywhere: median residual +0.57 dex at 100-300 kpc (CLASH), +0.35 at R2500, +0.15 at R500 and +0.06 at R200 (X-COP). In MOND mass the law needs 7.0 times the observed baryons at 100-300 kpc, 3.7 at R2500, 1.8 at R500 and 1.3 at R200.
- The 'factor 2' of the cluster lore is the R500 value. It is not a constant: the missing factor is 4-7 between 15 and 600 kpc (where g_bar is 0.1-0.8 a0 and the law gives little boost), 2 at R500 and about 1.3 at R200. The mass the law lacks is therefore concentrated inside ~0.5 Mpc, the known result of Sanders 2003 and Angus et al. 2008: it must be a component more centrally concentrated than the hot gas, with a mass of order 10^13-10^14 Msun per cluster.
- The cluster points form their own sequence with the RAR shape and a larger scale (CLASH: 15 x a0, Tian et al. 2020: 17 x; X-COP outskirts: 2.6 x). The two samples do not even agree on a cluster a0, because the residual depends on radius: there is no single a0 for galaxies and clusters, nor one for clusters.
- What this means for a MOND-type theory such as RMOND: an additional component of 1-5 times the observed baryons, concentrated in the cluster cores. Candidates in the literature are undetected baryons (cool gas, excluded at that level by X-ray and SZ data), sterile neutrinos of a few eV (the 'nuHDM' model, which also helps the CMB in RMOND-like theories; but a few-eV neutrino cannot cluster enough to reach a factor 5 in the cores), or a failure of the law on cluster scales. The data here cannot distinguish these; they fix the size and the location of the required component.
- Caveats: X-COP hydrostatic masses may be biased low by 5-15 % by non-thermal pressure (Eckert et al. 2019: median 6 % at R500), which would make the missing factor slightly larger; CLASH accelerations come from lensing and X-ray and are free of that bias; the baryon budgets include BCG stars and ICM gas, not the intracluster light, which adds at most a few per cent.

Figure: results/cluster_rar.png. Per-point values: results/cluster_rar.csv