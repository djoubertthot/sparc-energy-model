# Galaxy-cluster tables used by scripts/cluster_rar.py (downloaded 2026-10-12)

- `tian2020_clash_rar_fig2.tsv`, `tian2020_clash_table1.tsv`: Tian, Umetsu, Ko, Donahue & Chiu 2020, ApJ 896, 70, VizieR J/ApJ/896/70
  (asu-tsv from the CfA mirror): 20 CLASH clusters, log g_bar (stars of the central galaxy + X-ray gas) and log g_tot (strong and weak
  lensing masses of the CLASH programme) at 14-600 kpc with errors, 84 points. Parsed into `tian2020_clash_rar.csv` by
  `scripts/parse_cluster_tables.py`.
- `eckert2019.html`, `eckert2022.html`: ar5iv HTML of Eckert et al. 2019, A&A 621, A40 (non-thermal pressure in X-COP: hydrostatic
  M500 and M200, gas fractions at R500 and R200) and Eckert et al. 2022, A&A 662, A123 (the gravitational field of X-COP clusters: NFW
  M200 and c200 from X-ray and SZ profiles; gas and stellar fractions at R2500, R500 and R200 for five clusters). Parsed into
  `xcop_profiles.csv` (12 clusters; Zw1215 has no gas fractions and is not used) by the same script. Where the stellar fraction is not
  published, 1.5 % is assumed; where the R2500 gas fraction is not published, 0.8 times the R500 value; R2500 is taken as 0.45 R500.
