# Wide-binary data used by scripts/wide_binaries.py (downloaded 2026-10-13)

- `elbadry2021_d250.tsv`: El-Badry, Rix & Heintz 2021, MNRAS 506, 2269, VizieR J/MNRAS/506/2269/catalog (asu-tsv from the CfA mirror),
  the pairs with parallax of the primary > 4 mas (within 250 pc), binary type MSMS (two main-sequence stars) and chance-alignment
  probability R < 0.1: 132 829 pairs. Columns kept: Gaia EDR3 source ids, positions, parallaxes, proper motions with errors, RUWE,
  multi-peak fractions, duplicate flags, G magnitudes, BP-RP colours, Gaia DR2 radial velocities with errors, Galactic latitude, angular
  separation, projected separation in AU, binary type and R. Query:
  `-source=J/MNRAS/506/2269/catalog&-out.max=unlimited&Plx1=4..1000&BinType=MSMS&R=<0.1` with the column list in the file header.
- `EEM_dwarf_UBVIJHK_colors_Teff.txt`: the mean dwarf sequence of Pecaut & Mamajek 2013 (ApJS 208, 9), version 2022.04.16 from
  E. Mamajek's page, used for the mass-M_G relation and the main-sequence locus in M_G(BP-RP).
