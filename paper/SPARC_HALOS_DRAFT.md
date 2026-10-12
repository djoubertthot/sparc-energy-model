# Mass models of 171 SPARC galaxies with per-galaxy SED stellar masses: the acceleration relation, core-cusp preference and halo scaling relations

Draft v0, 2026-10-10. Every number below is produced by a script in `scripts/` from the data in `data/`; the results files in `results/` are the primary record.

## Abstract

We fit the rotation curves of the 171 SPARC galaxies with at least five measured points, using the baryonic mass models of Lelli, McGaugh & Schombert (2016) and a stellar mass-to-light ratio per galaxy from the SED fits of Marasco et al. (2025) where available (139 galaxies; median Υ★ = 0.47 at 3.6 µm). We measure the radial acceleration relation (2 795 points, scatter 0.147 dex, no residual correlation with luminosity, type, surface brightness, disc size, star-formation rate or gas fraction at |ρ| ≤ 0.2), the baryonic Tully-Fisher relation predicted from it with no free parameter (135 galaxies, scatter 0.088 dex), and the acceleration scale a₀ = 1.06 to 1.19 × 10⁻¹⁰ m s⁻², the range being set entirely by the stellar mass-to-light ratio. We then fit NFW and Burkert haloes with three parameters per galaxy. Burkert haloes are preferred in 69 % of the galaxies (75 % of late types), and the NFW concentrations scatter by 0.26 dex around the cosmological concentration-mass relation, against 0.11 dex expected. The halo core radius scales with the disc scale length (r₀ ∝ R_d^0.74) and with the maximum velocity (r₀ ∝ V_max^1.0), while the central core density is independent of V_max (slope −0.11, ρ₀ ≈ 0.06 M☉ pc⁻³); the core surface density ρ₀r₀ rises as V_max^0.89 and is therefore not universal. A test of the external field effect with an environmental field computed from the 2M++ catalogue finds no signal, with the caveat that the baryonic environmental field varies by only 0.24 dex across the sample. Extrapolated to 50-300 kpc, every halo fitted on the curves predicts a 30-50 % decline of the circular velocity by 300 kpc (median V(300 kpc)/V_flat = 0.68 for NFW, 0.57 for DC14, 0.46 for Burkert), where the weak-lensing profiles of isolated galaxies (Mistele et al. 2024) are flat to 5 %; the law alone gives 0.94, and its external field effect is excluded in its 1D form for e_N > 4e-5 a0, below the environmental field of the sample. Deep MeerKAT and 3D-modelled curves for six galaxies in common with SPARC show that one inclination revision (ESO 444-G084, 32 to 42-49 deg) moves a0 by a factor ten, that two well-measured dwarfs stay a factor two below the sample a0, and that the declining outer curve of NGC 7793 is confined to one annulus beyond which the gas is not in circular rotation. We discuss what these relations require from dark-matter haloes and from a modified acceleration law.

## 1. Data

- SPARC rotation curves and mass models (Lelli et al. 2016): 175 galaxies, 4 excluded with fewer than five points (D512-2, NGC6789, UGC00634, UGC07232).
- Stellar masses: Marasco et al. (2025) SED fits for 139 galaxies, Duey et al. (2025) WISE colour M/L for 17, z0MGS for 1, missing for 18 (Υ★ = 0.5 adopted). Dust masses (37), star-formation rates (67) and gas metallicities (30) from DustPedia and De Vis et al. (2019). Cross-matching and the adjustment log: `scripts/enrich_sparc.py`, `data/sparc_sample/`.
- Distances and inclinations: SPARC values; one distance uncertainty tightened (NGC 4559). Rotation curves byte-identical to the SPARC release.
- Environment: 2M++ (Lavaux & Hudson 2011), 67 402 galaxies with redshift distances, K-band luminosities converted to baryonic mass (0.8 M☉/L☉,K), completeness and faint-end corrections.

## 2. Methods

Baryonic velocity: V_bar² = V_gas|V_gas| + Υ_d V_disk|V_disk| + Υ_b V_bul|V_bul|, Υ_b = 1.4 Υ_d; g_bar = V_bar²/r.

Acceleration law: g_obs = ν(g_bar/a₀) g_bar, ν(y) = 1/(1 − e^{−√y}) (McGaugh et al. 2016); the simple function ν = ½ + √(¼ + 1/y) gives indistinguishable results.

Halo fits: V² = V_bar²(Υ_d) + G M_halo(<r)/r with NFW (ρ_s, r_s) or Burkert (ρ₀, r₀), Υ_d a nuisance parameter with a 0.1 dex log-normal prior on the SED value; least squares from nine starting points. Model comparison with AIC = χ² + 2k (k = 3 haloes, k = 1 law). NFW c₂₀₀ and M₂₀₀ from (ρ_s, r_s); fits with c₂₀₀ at the bound or M₂₀₀ > 3 × 10¹² M☉ (58 galaxies whose curves do not reach the turnover) are excluded from the concentration-mass statistics.

External field: e_N = |Σ G M_i/r_i²|/a₀ over 2M++ galaxies farther than 1 Mpc, galaxy placed at its SPARC distance; rotation curves refitted with the 1D QUMOND external-field function (Chae et al. 2020), marginalising Υ_d, distance and inclination. Calibration at the Milky Way: baryonic Newtonian e_N = 1.3 × 10⁻⁴, MOND-boosted 0.011 a₀, against 0.012 expected from the CMB dipole.

Scripts: `model_v2.py`, `fit_a0.py`, `fit_halos.py`, `halo_relations.py`, `test_efe.py`, `test_diversity.py`.

## 3. Results

### 3.1 The acceleration relation and its scale

| Quantity | Value |
|---|---|
| Points, galaxies | 2 795, 154 (points with δV/V < 0.1) |
| Scatter of log(g_obs/g_RAR), SED Υ★ | 0.147 dex (mean −0.015) |
| Scatter with fixed Υ★ = 0.5 / 0.7 | 0.144 dex |
| Residual vs log L[3.6], type, log SB, log R_d, log sSFR, f_gas (Spearman) | −0.08, +0.02, +0.02, −0.04, +0.19 (N = 64), +0.13 |
| a₀, Υ★ SED fixed / 0.5 fixed / free with prior | 1.06 / 1.19 / 1.09 × 10⁻¹⁰ m s⁻² |
| BTFR V_flat⁴ = G M_b a₀, zero free parameters, 135 galaxies | offset +0.026 dex, scatter 0.088 dex |

Per-galaxy SED masses do not reduce the scatter relative to a constant Υ★, which bounds the contribution of M/L variations to the relation's width.

### 3.2 Rotation curves with the law alone

| Baryons | N | median RMS | χ²/N | bias | early | intermediate | late |
|---|---|---|---|---|---|---|---|
| Υ★ = 0.5 / 0.7 | 171 | 13.7 km s⁻¹ | 11.3 | −4.1 | 24.9 | 10.7 | 10.8 |
| Υ★ SED per galaxy | 171 | 14.3 | 10.4 | −4.7 | 27.6 | 11.5 | 10.8 |
| Υ_d free, 0.1 dex prior | 171 | 9.2 | 3.7 | −1.7 | 14.3 | 7.5 | 8.1 |

The early types carry most of the residual; their bulge M/L is the least constrained input.

### 3.3 Halo fits

| Model | median RMS | median χ²/N | Σ AIC | lowest AIC in | median Υ_d/Υ_SED |
|---|---|---|---|---|---|
| Burkert + baryons | 3.7 km s⁻¹ | 0.59 | 10 439 | 63 % | 0.99 |
| NFW + baryons | 4.5 | 0.88 | 11 340 | 30 % | 0.97 |
| Acceleration law | 8.9 | 3.55 | 42 602 | 8 % | 0.86 |

Burkert is preferred over NFW in 69 % of galaxies (median Δχ² = −2.1), 75 % of late types, 59 % of early types; the preference does not correlate with central surface brightness (ρ = +0.02), gas fraction (−0.07) or V_max (+0.13, p = 0.09). By surface-brightness tercile: 72 %, 75 %, 60 % prefer a core.

### 3.4 Scaling relations

| Relation | Slope | Scatter | Spearman | N |
|---|---|---|---|---|
| log r₀ vs log R_d | +0.74 | 0.32 dex | +0.62 | 163 |
| log r₀ vs log V_max | +1.00 | 0.30 | +0.66 | 163 |
| log ρ₀ vs log V_max | −0.11 | 0.58 | −0.02 | 163 |
| log ρ₀r₀ vs log V_max | +0.89 | 0.33 | +0.62 | 163 |
| NFW log c₂₀₀, offset from Dutton & Macciò (2014) | +0.04 dex | 0.26 dex (0.11 expected) | | 113 |
| idem, late types / early types | −0.11 / +0.26 | 0.25 / 0.22 | | 47 / 31 |

The Burkert central density is constant across the sample, ρ₀ ≈ 0.06 M☉ pc⁻³, and the core radius grows with the galaxy; the product ρ₀r₀ (median 169 M☉ pc⁻², 16-84 %: 79-377) is therefore not the constant of Donato et al. (2009) but rises with V_max.

### 3.5 External field effect

152 galaxies with an environmental field: e_N,env median 1.5 × 10⁻⁴ (16-84 %: 1.2-2.0 × 10⁻⁴), local part within 10 Mpc 3.7 × 10⁻⁵. Spearman(e_N,fit, e_N,env) = +0.14 (p = 0.08); restricted to the 83 galaxies whose curves constrain e_N, +0.02 (p = 0.82). Low-acceleration residual vs e_N,env: 0.00. Outer slope vs e_N,env: −0.07 (p = 0.38). Declining curves (38) have a median environmental rank of 0.49. No detection; the dynamic range of the baryonic environmental field across SPARC is 0.24 dex, which limits the test. The large fitted e_N (0.1-1) belong to rising-curve dwarfs, where the external field absorbs the law's inner over-prediction rather than an outer decline.

### 3.6 Inner diversity

For 150 galaxies with data on both sides of 2 kpc, V(2 kpc)/V_max spans 0.46-0.86 (16-84 %) at V_max = 60-100 km s⁻¹, against 0.06 dex predicted by the cosmological concentration-mass scatter. The residual scatter at 2 kpc is 0.12 dex for the law and 0.10 dex for NFW matched to V_max; the law over-predicts the smallest dwarfs by 0.12 dex.

### 3.7 Feedback-shaped haloes and ΛCDM priors

The DC14 profile (Di Cintio et al. 2014), whose inner slope depends on M★/M₂₀₀, fitted with the same freedom as NFW and Burkert, converges to a median X = log(M★/M₂₀₀) = −1.81 and a median inner slope γ = 0.84. Its summed AIC is the lowest of the four models (DC14 10 289, Burkert 10 439, NFW 11 340, law 42 602) although Burkert wins more individual galaxies (51 % against 22 %). With ΛCDM priors imposed on the fits, concentration from Dutton & Macciò (2014) ± 0.11 dex and stellar-to-halo mass from Moster et al. (2013) ± 0.15 dex:

| Model | median χ²/N | galaxies with χ²/N < 2 |
|---|---|---|
| NFW, free | 0.88 | 70 % |
| NFW, ΛCDM priors | 1.50 | 61 % |
| DC14, ΛCDM priors | 0.89 | 77 % |
| Burkert, free | 0.59 | 82 % |

Under the priors DC14 beats NFW in 71 % of galaxies (77 % of late types), reproducing Katz et al. (2017): the cusp is incompatible with the cosmological parameter distributions, the feedback-flattened profile is not. Without priors, halo masses extrapolated from the curves are too uncertain to test the stellar-to-halo mass relation (scatter 0.9-1.0 dex); the fitted baryon fractions (M★+M_gas)/M₂₀₀ are 0.03-0.06, below the cosmic 0.157.

### 3.8 The shape of ν and the universality of a₀

Binned in g_bar, the relation is reproduced by the McGaugh function to 0.029 dex rms per bin, by Milgrom's standard function to 0.048 dex, by the n-family at n = 1.4 to 0.038 dex, and by Verlinde's relation, which has no free scale, to 0.066 dex; the deep-regime slope over five bins is 0.498. The bin medians are biased by 0.04 dex in a₀ (Eddington bias on the noisy abscissa), so a₀ is taken from the point-wise fits of § 3.1.

Fitted galaxy by galaxy with Υ_d, distance and inclination marginalised over Gaussian priors, a₀ has a median of 0.93 × 10⁻¹⁰ m s⁻², a weighted mean of 1.23, an observed scatter of 0.43 dex and an intrinsic scatter of 0.33 dex after subtracting the per-galaxy intervals (median 0.13 dex); a universal a₀ gives χ² = 6 842 for 165 degrees of freedom, and 28 % of the galaxies exclude 1.2 × 10⁻¹⁰ by a factor of two. The per-galaxy a₀ correlates weakly with V_max (ρ = +0.20, p = 0.01) and with the quality flag (ρ = −0.22, p = 0.005). This reproduces the result of Rodrigues et al. (2018); its interpretation as a genuine non-universality was disputed by Kroupa et al. (2018) and Li et al. (2018), who attribute the excess to the error model of SPARC (non-circular motions, distance and inclination systematics beyond Gaussian). The correlation with the quality flag supports the latter reading; the correlation with V_max does not settle it.

A direct test uses Cosmicflows-4 (Tully et al. 2023): 110 SPARC galaxies are matched, 53 have a rotation-independent distance (TRGB for 36, Cepheids, SBF, supernovae) and 44 of these are tighter than the SPARC value; the independent distances agree with SPARC's to −0.010 dex in the median with 0.048 dex scatter, but six galaxies move by more than 20 %. Refitting a₀ with the new distances:

| Sample | Distances | N | observed scatter | intrinsic scatter | χ² of a universal a₀ | exclude 1.2 × 10⁻¹⁰ by ×2 |
|---|---|---|---|---|---|---|
| replaced galaxies | SPARC | 43 | 0.35 dex | 0.28 | 1 308 | 28 % |
| replaced galaxies | CF4 | 43 | 0.32 | 0.29 | 948 | 28 % |
| replaced, quality 1 | SPARC | 22 | 0.30 | 0.25 | 944 | 23 % |
| replaced, quality 1 | CF4 | 22 | 0.21 | 0.19 | 496 | 14 % |

The per-galaxy change of log a₀ anticorrelates with the change of distance (ρ = −0.64), as expected from a₀ ∝ 1/D in the deep regime, and the two most discrepant galaxies with SPARC distances, NGC 2903 (a₀ = 3.7 × 10⁻¹⁰) and UGC 7603 (3.2), fall to 1.3 and 1.6 with their TRGB distances. Better distances reduce the scatter of a₀ by a third in the best-quality subset without removing it; the remaining 0.19 dex is to be compared with the inclination and M/L systematics that the Gaussian priors do not capture.

Two independent kinematic datasets test the inclination term. (i) MHONGOOSE (de Blok et al. 2024) shares five galaxies with SPARC, three with public deep data; tilted rings fitted on the MeerKAT velocity fields (`mhongoose_rings.py`) recover the SPARC geometry for DDO 161 (i = 62° vs 70 ± 10) and NGC 7793 (51° vs 47 ± 9) and reproduce the SPARC velocities to 8 %, but give i = 42° for ESO 444-G084 against 32 ± 6 in SPARC, in agreement with the 3D model of Namumba et al. (2025), 49 ± 1.2°. For this galaxy the SPARC curve implies a₀ = 3.5 × 10⁻¹⁰ and the 3D curve 0.3-0.4 × 10⁻¹⁰: one inclination revision moves a₀ by a factor ten. The deep maps extend the curves to 1.3-1.7 R_last and g_bar = 0.007-0.03 a₀: DDO 161 stays on the relation (extension residual −0.22 dex against −0.29 inside), ESO 444-G084 as well within its errors, while NGC 7793's known decline (Carignan & Puche 1990; de Blok et al. 2008) continues beyond R_last. Refitting its rings with V_rot, V_rad, position angle and inclination all free (`ngc7793_warp.py`) separates two zones. Inside R_last the inclination wanders between 40° and 52° from ring to ring, as Gentile et al. (2011) noted for the THINGS data, the non-circular term is 10 % of the rotation, the rms per ring 3-7 km s⁻¹, and the law holds with a median residual of +0.02 dex. Between 8 and 9.6 kpc the decline survives the free geometry (residual −0.35 dex; the law would need i ≈ 33° where the rings give 48-60°), but the non-circular term grows to 20-40 % of the rotation and the rms doubles; beyond 10 kpc the radial term equals or exceeds the rotation (|V_rad|/V_rot = 0.8-2.4, position angle swinging over 80°): that gas is not in circular rotation and its velocities do not measure the potential. NGC 7793 therefore remains the one unfavourable case of the sample, confined to one annulus at 8-10 kpc and 0.03-0.04 a₀, where only a 3D model of the cube with a warp can tell a genuine decline from a disturbed outer disc. (ii) The 3D-Barolo curves of Iorio et al. (2017) for DDO 154, DDO 168 and NGC 2366, with inclinations to 1-4° and asymmetric-drift corrections of 9-20 km s⁻¹ at the last point, agree with the SPARC velocities to 1-3 % and give a₀ = 1.06 (1.06-1.19), 0.56 (0.42-0.84) and 0.63 (0.53-0.79) × 10⁻¹⁰ with χ²/N of 1.0-3.3, against 0.94, 0.89 and 0.63 from the SPARC curves. For these three galaxies the geometry is not the issue: two of them sit a factor two below the sample value with the best kinematics available. In the deep regime a 0.3 dex offset in a₀ is a 0.15 dex offset in g_obs, the same size as the scatter of the relation itself; the per-galaxy a₀ scatter and the RAR scatter are one quantity seen from two sides.

### 3.9 Beyond the disc: weak lensing at 50-300 kpc

Mistele et al. (2024) deproject the KiDS-1000 excess surface density around the isolated galaxies of Brouwer et al. (2021) into circular velocities and find them flat from 50 kpc to at least 300 kpc, with no sign of a decline out to 1 Mpc: V_flat = 137 ± 12, 182 ± 11, 205 ± 8 and 260 ± 6 km s⁻¹ for log M_b = 9.0-10.5, 10.5-10.8, 10.8-11.1 and 11.1-11.5, on the kinematic baryonic Tully-Fisher relation. The haloes fitted on the SPARC curves (§ 3.3, 3.7) extrapolate to these radii with no new freedom (`lensing_vflat.py`), the baryons being a point mass M_b = Υ_fit L₃₆ + 1.33 M_HI. Median over the 135 galaxies with a kinematic V_flat:

| Model | V(50)/V_flat | V(100)/V_flat | V(200)/V_flat | V(300)/V_flat | 16-84 % at 300 kpc |
|---|---|---|---|---|---|
| NFW, free | 0.97 | 0.91 | 0.75 | 0.68 | 0.47-1.81 |
| NFW, ΛCDM priors | 0.96 | 0.87 | 0.74 | 0.66 | 0.55-0.84 |
| DC14, ΛCDM priors | 0.83 | 0.73 | 0.63 | 0.57 | 0.41-0.91 |
| Burkert | 0.81 | 0.67 | 0.53 | 0.46 | 0.32-0.65 |
| law, isolated | 0.96 | 0.96 | 0.95 | 0.94 | 0.86-1.05 |
| law + external field, e_N = 0.02 a₀ | 0.49 | 0.35 | 0.25 | 0.20 | 0.14-0.33 |
| lensing (Mistele et al. 2024) | 1 | 1 | 1 | 1 | ± 0.04-0.08 per bin |

Bin by bin, for log M_b = 10.5-10.8 (21 SPARC galaxies, lensing 182 ± 11 km s⁻¹) the models give at 300 kpc: NFW 114 (97-305), NFW with priors 123 (109-166), DC14 106 (70-192), Burkert 95 (77-143), law 158 (154-170); for 10.8-11.1 (15 galaxies, lensing 205 ± 8): 152, 158, 135, 135 and 197. Every halo fitted on the rotation curves predicts a decline of 30-50 % by 300 kpc that the lensing profiles do not show. The NFW fits have a median R₂₀₀ of 165 kpc: the lensing data straddle the virial radius, where the √(ln r / r) decline of a ΛCDM halo must set in, and see nothing. The isolated law reproduces the flatness to 5-10 %.

The external field effect removes this agreement. Computed with the Chae et al. (2020) function, for the 63 galaxies with log M_b ≥ 10 the ratio at 300 kpc stays within 10 % of its isolated value only for e_N < 4 × 10⁻⁵ a₀, below the baryonic environmental field of the SPARC galaxies themselves (1.5 × 10⁻⁴, § 3.5) and two orders of magnitude below the e_N ≈ 0.01-0.03 of the Milky Way or of the detection claimed by Chae et al. (2021). The 1D function is the strongest form of the effect and the bound applies to that form, but the sign is general: an external field makes the law Newtonian with a boosted G beyond the radius where g_bar falls under e_N a₀, so circular velocities must turn over, and the lensing curves do not. Two independent null results therefore point the same way, the absence of an environmental signature inside the discs (§ 3.5) and the flat lensing velocities outside them: the law holds as an isolated-galaxy relation, without the external field effect that a modified-inertia reading requires. Caveats: the lensing M_b excludes hot gas, which could raise g_bar by a factor 2-3 at these radii (Brouwer et al. 2021) and would move the law prediction up by 20-30 % and the halo predictions by less than 10 %, in neither case towards a decline; the two-halo term is negligible below 300 kpc for isolated lenses; the SPARC bins are not mass-matched to the lensing stacks, which is why the shape ratio V(R)/V_flat is the quantity compared.

### 3.10 The scale a₀ across cosmic time

Two readings of a₀ ≈ cH₀/2π make opposite predictions: an expansion-tied scale, a₀ ∝ H(z), is 1.8 times larger at z = 1, 3.3 at z = 2.2 and 7 at z = 4.5; a vacuum-tied scale, a₀ ∝ c√(Λ/3), is constant. We test r = a₀(z)/a₀(0) on every published rotation curve with a total and a baryonic acceleration at the same radius (`a0_highz.py`, `a0_highz_coldgas.py`, `a0_highz_shape2.py`): RC100 (Nestor Shachar et al. 2023), KMOS³ᴰ (Übler et al. 2017), KROSS (Harrison et al. 2017), MSA-3D (JWST/NIRSpec), KURVS (Puglisi et al. 2023), NOEMA3D (Jolly et al. 2026), CRISTAL (Lee et al. 2025), and the cold-gas discs of Rizzo et al. (2020, 2021), Roman-Oliveira et al. (2024), Lelli et al. (2023) and ALPAKA (Rizzo et al. 2023). Two tests. (i) Amplitude: the mass discrepancy at the measured radius against ν(g_bar/(r a₀)). At z ≥ 1.7, 60 galaxies give r = 1.05 (0.87-1.23) where the expansion-tied scale requires 3.3-8.5; but a uniform 0.1 dex overestimate of the baryonic masses makes the two hypotheses equivalent and 0.2 dex reverses the verdict, a precision that SED stellar masses and α_CO or α_[CII] do not reach. At z ≈ 1 the samples sit in the regime g_bar ≈ 0.6-0.9 a₀, where a₀ ∝ 1/M_bar, and the evolution claimed by Ciocan et al. (2026) from MUSE at z < 1.44 (a₀ rising from 2.0 to 2.7 × 10⁻¹⁰ m s⁻²) disappears once the molecular gas their mass budget omits is added with the Tacconi et al. (2018) scaling (a₀ = 1.3-1.6 → 0.9-1.3 × 10⁻¹⁰). (ii) Shape: the rise of ν with radius along each cold-gas curve, with the total baryonic mass free (0.3 dex prior), which does not depend on α_CO, M/L or the gas fraction. On 12 curves at z = 1.5-4.8 (six lensed [CII] discs rebuilt from the published kinematic and structural parameters of Rizzo et al., four [CII] discs of Roman-Oliveira et al., two CO discs of Lelli et al.), the expansion-tied scale lies Δχ² = 21 above the best fit and the constant scale 0.9 above it; with the 19 ALPAKA curves added, a₀ ∝ H(z)ⁿ gives n = 0.05 (−0.25 to 0.30). Ten of the twelve are at z > 4, where a disc rotating at 200-500 km s⁻¹ within 2-4 kpc would have to see ν double between R_e and 3 R_e; none does. Ionised-gas slopes (KURVS, v/σ₀ ≈ 2) cannot enter this test: the pressure correction, 6-12 σ₀², exceeds the 5 % maximum outer rise the law allows for an exponential disc. The scale of the law does not follow the expansion.

## 4. Discussion

What a dark-matter halo must reproduce: a relation between g_obs and g_bar with 0.15 dex scatter and no memory of luminosity, star formation or gas content; a core radius locked to the disc scale length (ρ = +0.62) and to V_max with unit slope; a constant central density; a concentration scatter 2.4 times the cosmological one, with late types below and early types above the relation; and circular velocities that stay flat to 300 kpc across the virial radius, where every profile fitted on the discs declines by 30-50 % (§ 3.9). These are the known challenges of galaxy formation in ΛCDM (cusp-core, diversity, the baryon-halo conspiracy); here they are measured on one homogeneous sample with individual stellar masses.

What a modified acceleration law must reproduce: nothing further on these data, it is the law that generates the relations above; its open problems lie elsewhere, in clusters, in the CMB, and in the external field effect, which the disc test does not detect and the lensing flatness excludes at the level the environment implies (§ 3.9); a modified-inertia reading predicts that effect, a dark-component reading does not, and the data side with the second.

The value of a₀ sits in the band 1.06-1.19 × 10⁻¹⁰ m s⁻², which contains cH₀/2π (1.08-1.13) and cH₀/6 (1.13-1.18) and lies 10-15 % above c√(Λ/3)/2π; the prefactors being heuristic, the value alone does not decide between an expansion-tied and a vacuum-tied scale. The shape of the cold-gas rotation curves at z = 1.5-4.8 does (§ 3.10): a₀ is constant to n = 0.05 (−0.25 to 0.30) in a₀ ∝ H(z)ⁿ, which excludes the expansion-tied reading at 4.6-5.2σ and leaves the vacuum-tied one, a₀ ≈ 0.2 c√(Λ/3), with a prefactor to explain. The baryonic Tully-Fisher zero points of Übler et al. (2017) at z ≈ 0.9 and 2.3 (−0.44 and −0.27 dex) are reproduced to 0.02 dex by the law with a constant a₀ once v_circ,max is evaluated at 2.2 R_d in the high-acceleration regime of those discs.

Limitations: SPARC distances and inclinations are taken at face value in the law fits; bulge M/L is tied to the disc value; 58 NFW fits are unconstrained; the environmental field uses redshift distances with ~4 Mpc errors on neighbours; the 2 kpc test is affected by non-circular motions and beam smearing.

## 5. Reproducibility

`python scripts/enrich_sparc.py` builds the table; `model_v2.py`, `fit_a0.py`, `fit_halos.py`, `halo_relations.py`, `test_efe.py` (needs `data/catalogs/2mpp_catalog.dat`), `test_diversity.py`, `fit_dc14.py`, `fit_halos_priors.py`, `measure_nu.py`, `a0_universality.py`, `a0_universality_cf4.py`, `mhongoose_rings.py` and `ngc7793_warp.py` (need the MHONGOOSE moment maps in `data/external/mhongoose/`), `littlethings3d_a0.py`, `lensing_vflat.py`, `a0_highz.py`, `a0_highz_coldgas.py`, `a0_highz_shape.py` and `a0_highz_shape2.py` (need the tables in `data/external/highz/`, documented in its README) produce every table in `results/`.

## References

Brouwer M. M. et al. 2021, A&A 650, A113 · Carignan C., Puche D. 1990, AJ 100, 394 · Chae K.-H. et al. 2020, ApJ 904, 51 · Chae K.-H. et al. 2021, ApJ 921, 104 · Ciocan B. I. et al. 2026, A&A (arXiv:2604.22613) · de Blok W. J. G. et al. 2008, AJ 136, 2648 · de Blok W. J. G. et al. 2024, A&A 688, A109 (MHONGOOSE) · De Vis P. et al. 2019, A&A 623, A5 · Di Cintio A. et al. 2014, MNRAS 441, 2986 · Gentile G., Famaey B., de Blok W. J. G. 2011, A&A 527, A76 · Harrison C. M. et al. 2017, MNRAS 467, 1965 · Iorio G. et al. 2017, MNRAS 466, 4159 · Jolly J.-B. et al. 2026, A&A (arXiv:2604.18503) · Namumba B. et al. 2025, A&A (arXiv:2506.04101) · Donato F. et al. 2009, MNRAS 397, 1169 · Dutton A., Macciò A. 2014, MNRAS 441, 3359 · Katz H. et al. 2017, MNRAS 466, 1648 · Kroupa P. et al. 2018, Nat. Astron. 2, 925 · Lavaux G., Hudson M. 2011, MNRAS 416, 2840 · Lee L. L. et al. 2025, A&A (arXiv:2507.11600) · Lelli F., McGaugh S., Schombert J. 2016, AJ 152, 157 · Lelli F. et al. 2017, ApJ 836, 152 · Lelli F. et al. 2023, A&A 672, A106 · Li P. et al. 2018, A&A 615, A3 · Marasco A. et al. 2025, A&A 695, L23 · Milgrom M. 2017, arXiv:1703.06110 · Mistele T., McGaugh S., Lelli F., Schombert J., Li P. 2024, ApJL 969, L3 · McGaugh S., Lelli F., Schombert J. 2016, PRL 117, 201101 · Moster B. et al. 2013, MNRAS 428, 3121 · Nersesian A. et al. 2019, A&A 624, A80 · Nestor Shachar A. et al. 2023, ApJ 944, 78 · Oman K. et al. 2015, MNRAS 452, 3650 · Price S. H. et al. 2021, ApJ 922, 143 · Puech M. et al. 2010, A&A 510, A68 · Puglisi A. et al. 2023, MNRAS 524, 2814 · Rizzo F. et al. 2020, Nature 584, 201 · Rizzo F. et al. 2021, MNRAS 507, 3952 · Rizzo F. et al. 2023, A&A 679, A129 · Rodrigues D. et al. 2018, Nat. Astron. 2, 668 · Roman-Oliveira F. et al. 2024, A&A (arXiv:2403.00904) · Tacconi L. J. et al. 2018, ApJ 853, 179 · Tully R. B. et al. 2023, ApJ 944, 94 · Übler H. et al. 2017, ApJ 842, 121 · van der Wel A. et al. 2014, ApJ 788, 28 · Verlinde E. 2016, SciPost Phys. 2, 016.
