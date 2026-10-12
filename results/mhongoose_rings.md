# MHONGOOSE tilted-ring curves vs SPARC (2026-10-10)

## DDO161 (J1303-17b), SPARC D = 7.5 Mpc, i = 70 +/- 10, R_last = 13.4 kpc, N_sparc = 31
geometry fitted on the 26" map inside R_last: Vsys = 738.3 km/s, PA = 293 deg, i = 62 deg (SPARC 70), centre offset (-20", 17")
inside R_last: MHONGOOSE 26" rings vs SPARC, median V_MH / V_SPARC = 0.922 (26 rings)
extension: 9 rings beyond R_last up to 17.8 kpc = 1.33 R_last; outer V = 73.7 +/- 2.7 km/s; log g_bar there = -11.94 (0.010 a0)
RAR residual log(g_obs/g_RAR) in the extension: median -0.220, range -0.28 to -0.09; inside R_last (26"): median -0.289

## ESO444-G084 (J1337-28), SPARC D = 4.83 Mpc, i = 32 +/- 6, R_last = 4.4 kpc, N_sparc = 7
geometry fitted on the 26" map inside R_last: Vsys = 585.3 km/s, PA = 84 deg, i = 42 deg (SPARC 32), centre offset (19", 2")
inside R_last: MHONGOOSE 26" rings vs SPARC, median V_MH / V_SPARC = 0.686 (13 rings)
extension: 6 rings beyond R_last up to 6.2 kpc = 1.39 R_last; outer V = 49.0 +/- 4.8 km/s; log g_bar there = -12.07 (0.007 a0)
RAR residual log(g_obs/g_RAR) in the extension: median +0.081, range -0.02 to +0.28; inside R_last (26"): median -0.129

## NGC7793 (J2357-32), SPARC D = 3.61 Mpc, i = 47 +/- 9, R_last = 7.9 kpc, N_sparc = 46
geometry fitted on the 26" map inside R_last: Vsys = 225.6 km/s, PA = 291 deg, i = 51 deg (SPARC 47), centre offset (16", -4")
inside R_last: MHONGOOSE 26" rings vs SPARC, median V_MH / V_SPARC = 0.925 (33 rings)
extension: 24 rings beyond R_last up to 13.2 kpc = 1.68 R_last; outer V = 29.7 +/- 13.7 km/s; log g_bar there = -11.38 (0.034 a0)
RAR residual log(g_obs/g_RAR) in the extension: median -0.406, range -1.05 to -0.17; inside R_last (26"): median -0.034

Figure: results/mhongoose_rings.png

## Reading of the three results (first pass, 2D tilted rings on moment maps)

- Inside the SPARC radial range the 26" rings reproduce the SPARC velocities to 8 % for DDO 161 and NGC 7793 with the
  centre, PA and inclination recovered independently (DDO 161: i = 62 vs 70 +/- 10; NGC 7793: 51 vs 47 +/- 9).
- ESO 444-G084: the fitted inclination is 42 deg, the MHONGOOSE 3D model of Namumba et al. 2025 gives 49 +/- 1.2, SPARC adopted
  32 +/- 6. The SPARC velocities are 1.3 to 1.5 times too high for this galaxy; with the corrected inclination the galaxy moves from
  a0 = 3.5e-10 (SPARC curve) to ~0.3-0.4e-10 (published curve, see results/mhongoose_eso444.md). This single galaxy illustrates the
  inclination term of the a0 scatter.
- DDO 161: 9 rings beyond R_last to 17.8 kpc (1.33 R_last), g_bar down to 0.010 a0; the outer velocities rise from 62 to 74 km/s. The RAR
  residual in the extension (-0.22) equals the residual inside (-0.29) within 0.07 dex: the curve stays on the relation. The rise itself
  needs a 3D model (edge-of-map rings, 60" beam).
- NGC 7793: 24 rings to 13.2 kpc (1.68 R_last); the known decline continues, V falls to ~40 km/s at 12 kpc against V_bar ~ 45, i.e. the
  curve approaches the baryonic curve and the RAR residual reaches -0.4 dex at 0.034 a0. A warp is known in this galaxy and the ring
  fit holds the inclination fixed at 51 deg; a more face-on outer disc would raise the true velocities. This is the galaxy where a full
  3D model (cube) is worth the effort: a rotation curve that falls to the baryonic value at 1.7 R_last is incompatible with the
  acceleration law and requires a truncated halo in the dark-matter reading.

Method limits: velocity fields (moment 1) rather than cubes, so no beam-smearing or asymmetric-drift correction; PA free per ring
within +/- 60 deg but inclination fixed; rings accepted only with > 50 % azimuthal coverage and mean N_HI > 1e18 (60") or 5e18 (26").

Note (2026-10-11): the east-west sign of the sky offsets was corrected in ring_coords (PA was mirrored, 360 - PA); velocities, inclinations and residuals are unchanged. NGC 7793: PA = 291 vs 290.1 in THINGS (de Blok et al. 2008). NGC 7793 with inclination free per ring: results/ngc7793_warp.md.
