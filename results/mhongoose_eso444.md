# ESO 444-G084: SPARC vs the MHONGOOSE 3D rotation curve (Namumba et al. 2025, arXiv:2506.04101)

| | SPARC (Lelli+2016) | MHONGOOSE (PyFAT / TiRiFiC 3D) |
|---|---|---|
| distance | 4.83 Mpc (TRGB) | 4.6 +/- 0.4 Mpc |
| inclination | 32 +/- 6 deg | 49.0 +/- 1.2 deg |
| last radius | 4.44 kpc | 3.87 kpc |
| velocity at 3.4-3.9 kpc | 63 km/s | 51-52 km/s (asymmetric-drift corrected) |
| N points | 7 | 16 |

Velocity ratio MHONGOOSE / SPARC at common radii: 0.48 at 0.4 kpc rising to 0.83 at 3.9 kpc (sin 32 / sin 49 = 0.70).
The inner discrepancy exceeds the inclination factor: the SPARC curve is also beam-smeared / less resolved inside 2 kpc.

a0 fitted on the same radial range with the SPARC baryons (Y* = 0.5 and 0.3):
- SPARC curve: a0 = 3.5-3.8e-10 (chi2/N 1.5-2.0)
- MHONGOOSE curve: a0 = 0.32-0.38e-10 (chi2/N 46-56: the RAR over-predicts the inner curve by ~4 km/s and under-predicts the
  outer point by ~8 km/s with these small quoted errors)

Lesson: one inclination revision of 17 deg moves this galaxy's a0 by a factor 10. The per-galaxy a0 scatter of 0.33 dex measured
on SPARC includes such cases; the 2D-fitted MHONGOOSE inclination (42 deg, results/mhongoose_rings.md) agrees with the 3D value.
Data: data/external/eso444g084_mhongoose.csv (Table 5 of the paper).
