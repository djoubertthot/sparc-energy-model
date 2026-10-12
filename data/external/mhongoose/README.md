# MHONGOOSE data used here (de Blok et al. 2024, A&A 688, A109; survey site mhongoose.astron.nl)

All FITS files are ignored by git (`*.fits` in .gitignore); fetch them with `gdown <id>` from the public Google Drive
releases listed on the survey site (Deep Data Release 2 folder 1znmFHW71juiWP6QBY6-aE6CyAAVvZjyK, single-track release
folder 1KMvBhAUflzMTrPsrWgMll6tFD7Ggaakw).

- `J1303-17b/` (DDO 161), `J1337-28/` (ESO 444-G084), `J2357-32/` (NGC 7793): moment maps (mom0, mom0_pb, mom1, mom2) at the
  r05_t60, r10_t00 and r15_t00 resolutions of the deep releases, used by `scripts/mhongoose_rings.py` and `scripts/ngc7793_warp.py`.
- `J2357-32/r15_t00_mask.fits.gz`, `r15_t00_mask-2d.fits`, `r15_t00_primbeam.fits`: SoFiA mask and primary beam of the 26" data
  (DR2 ids 1dZtTPmnABWCMpPA1u0xJKXkJS2DUgfid, 1snSfBU7WQR0VUua3pxI14Xi6e2j4Rjaw, 1PdLALBkFWdvjx410kqbuGgrFwCNJHf_5).
- `J2357-32/cube/J2357-32_r15_t00_cube.fits`: the 26" data cube of NGC 7793, DR2 file
  `cubes/J2357-32/comb_10tracks_J2357_32_HI_r15_t0_mg095_2.image-ics.fits` (id 1P7cV2ZR5hZ4HF8L1DJwOx2us3lkBjvbB, 1.63 GB,
  772 x 772 pixels of 7", 684 channels of 1.38 km/s, beam 34.8" x 26.3" at PA 139, Jy/beam, rms 0.15 mJy/beam per channel).
  `scripts/ngc7793_cube_cut.py` cuts it to `ngc7793_r15_sub.fits` (301 x 301 pixels around the galaxy, 90-360 km/s, channels
  rebinned to 5.5 km/s) and the matching mask; `scripts/ngc7793_cube3d.py` fits the 3D tilted-ring models on the cut cube.
- `DeepDataRelease2_README.pdf`: the release notes (id 1O34h8T16OsyUq-Q-xQ9GbO3HkU41GnEK).
