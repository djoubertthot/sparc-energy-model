# High-redshift rotation-curve data used by scripts/a0_highz.py (downloaded 2026-10-11)

- `rc100_nestorshachar2023_table3.csv`: Nestor Shachar et al. 2023, ApJ 944, 78, Table 3 (arXiv:2209.12199). The table is an image in the PDF;
  transcribed by hand from 200 dpi renders and checked column by column: log Sigma_DM(<Re) recomputed from f_DM, V_c(Re), Re agrees with the
  printed column to 0.005 dex (median), the 7 larger deviations (< 0.1 dex) are rounding of f_DM <= 0.06.
- `price2021_rc41_*.csv`: Price et al. 2021 (arXiv:2109.02659) Tables 1 and 3 for the RC41 subset, parsed from the arXiv HTML; used as a cross-check.
- `msa3d_2026.csv`: MSA-3D JWST/NIRSpec rotation curves (arXiv:2606.27853), galaxy_parameters.tbl + fitted_parameters.tbl from the arXiv source.
- `ubler2017_kmos3d_tfr.*`: Ubler et al. 2017, ApJ 842, 121, VizieR J/ApJ/842/121 (mirror vizier.cfa.harvard.edu).
- `harrison2017_kross.*`: Harrison et al. 2017, MNRAS 467, 1965, VizieR J/MNRAS/467/1965 (KROSS v2 catalogue).
- `lelli2016_btfr.*`: Lelli et al. 2016, ApJL 816, L14, VizieR J/ApJ/816/L14 (local bTFR calibration).
- `z45/circular_speed_romanoliveira2024.csv`: Roman-Oliveira et al. 2024 (A&A 684, A118), asymmetric-drift-corrected circular speeds of
  BRI1335-0417, J081740, SGP38326-1 and SGP38326-2 from Zenodo record 10707348 (files circular_speed_*.dat).
  Masses and sizes of these and of the Rizzo et al. 2021 / Fraternali et al. 2021 discs are hard-coded in scripts/a0_highz.py with their sources.
- `noema3d_jolly2026_dysmalpy.csv`: Jolly et al. 2026 (arXiv:2604.18503) Tables 1 and 3 and Chen et al. 2026 (arXiv:2604.18504), transcribed
  from the LaTeX sources: CO gas masses, DysmalPy log M_bary, sigma0, V_c(Re), f_DM(Re), Re_disk, Re_bulge, B/T, inclination.
- `cristal_lee2025_dysmalpy.csv`: Lee et al. 2025 (arXiv:2507.11600) Tables 1, 2 and 4 (DysmalPy fits of the 14 CRISTAL discs at z = 4.4-5.7).
- `alpaka_romanoliveira2026.csv`: Roman-Oliveira, Rizzo & Fraternali 2026 (arXiv:2601.03338) Tables 1, 4, 5 with V_ext, sigma_ext, V_max and
  inclinations from Rizzo et al. 2023 (arXiv:2303.16227) Tables 3 and 4. Rizzo et al. 2021 (arXiv:2102.05671) baryonic masses, half-mass
  radii and asymptotic velocities are hard-coded in scripts/a0_highz_coldgas.py.
- `kurvs_puglisi2023.csv`: Puglisi et al. 2023 (MNRAS 524, 2814; arXiv:2305.04382) Tables 1-4 transcribed from the LaTeX source: 22 KURVS-CDFS
  galaxies at z ~ 1.5, M*, SFR, Reff, inclinations, sigma0, v/sigma0, outer slope t = v(6 Rd)/v(3 Rd), v at 3 and 6 Rd, f_DM(<Re) for ten discs.
- ALPAKA kinematic profiles: the survey site (alpaka-survey.github.io) only serves a catalogue of global properties (ALPAKA_II_v1.fits,
  57 galaxies, sigma, Vmax, V/sigma, M*, SFR); the per-ring profiles announced as public were not found online.
- `alpaka_rizzo2023_rings_byeye.csv`: 3DBarolo ring rotation velocities of the 19 ALPAKA discs READ BY EYE from Fig. 6 of Rizzo et al. 2023
  (figures/vrot.png in the arXiv source, 937 x 1020 pixels): 2-5 rings per galaxy, radii in arcsec, velocities to ~5 %, errors from the
  shaded bands. Not the authors' numbers.
- `lelli2023_rings_byeye.csv`: ring velocities of zC-400569 and zC-488879 READ BY EYE from Fig. 7 of Lelli et al. 2023 (Rotcur.jpg in the
  arXiv source, 1984 x 982 pixels), CO(3-2), CO(4-3), CO(2-1) and H-alpha points with their error bars. Not the authors' numbers.
- Rizzo et al. 2020 (SPT0418-47) and 2021 (five SPT discs): rotation-curve parameters (V_t, R_t, beta, xi), dispersion profiles (sigma0, R_sigma),
  stellar Sersic parameters and gas scale lengths are taken from their tables and hard-coded in scripts/a0_highz_shape2.py; the circular
  velocities are rebuilt analytically with their asymmetric-drift formula.
