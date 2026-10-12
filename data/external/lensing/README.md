# Weak-lensing circular velocities used by scripts/efe_lensing_mpc.py (downloaded 2026-10-12)

- `mistele2024_vc_mrt.txt`: Mistele, McGaugh, Lelli, Schombert & Li 2024, ApJL 969, L3, machine-readable Table 1 from the AAS/IOP
  article page (fixed width, byte layout in the header): sample (All, LTG, ETG, GAMA), R in kpc, circular velocity with statistical and
  systematic errors for the four baryonic-mass bins log M_b = 9.0-10.5, 10.5-10.8, 10.8-11.1 and 11.1-11.5 (bin means 10.10, 10.66,
  10.96, 11.29), from 34 kpc to 2.5 Mpc. The velocities are the spherical deprojection of the stacked KiDS-1000 excess surface density
  around isolated lenses (isolation criterion of Brouwer et al. 2021, which the authors trust to about 1 Mpc). The "All" sample is used.
- `mistele2024_vc.csv`: the part of the same table present in the arXiv source (plots/Vc-table.tex of arXiv:2406.09685), parsed by
  `scripts/parse_mistele_vc.py`; kept as a cross-check, the MRT above is the reference.
