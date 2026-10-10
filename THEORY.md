# SPARC Energy-Mass Model

## Fundamental Motivation
The model originates from the idea that energy trapped or traversing baryonic matter over long timescales contributes an apparent mass via the relativistic relation:
\[
E = \gamma m c^2 \quad \Rightarrow \quad \Delta m = \frac{E}{\gamma c^2}
\]

Photons (and other particles) emitted or absorbed inside stars and galaxies spend significant time interacting with matter. During that time they contribute an effective mass. External radiation (cosmic, radio, neutrino background) can similarly accumulate an apparent mass density.

## Density Amplification Formula
\[
\rho_{\rm eq}(r) = \rho_0 \times \left(\frac{a_0}{a_{\rm bary}(r)}\right)^\beta
\]

### Origin of each term
- **ρ₀**: Characteristic energy density scale. Related to the average energy density of trapped radiation or cosmic backgrounds that can be converted into equivalent mass. Empirically tuned; order of magnitude consistent with baryonic energy densities.
- **a₀**: Acceleration scale (≈ 1.2×10^{-10} m/s²). Anchored to the cosmological acceleration scale (c × H₀ / 2π or similar). Below this acceleration the trapping/amplification becomes significant.
- **(a₀ / a_bary)^β**: Amplification factor. In low-acceleration regions (outer disks), the relative contribution of trapped energy grows. β ≈ 1.4 is an empirical exponent that best matches SPARC residuals; physically it may relate to optical depth, scattering, or time-of-flight effects.
- **Radial integration**: The equivalent density is integrated as dM = 4π R² ρ_eq dR to obtain the cumulative dynamical mass contribution.

## Age / Time Factor
A binary morphological proxy (early-type 1.2, late-type 0.8) is used as a crude stand-in for stellar age / accumulation time. Real ages would allow a more physical ∫ E(t) dt term.

## Methodology
1. Compute baryonic circular velocity and acceleration from SPARC gas + disk + bulge components.
2. Apply the amplification to obtain ρ_eq(r).
3. Integrate to obtain M_eq(r).
4. Add the resulting V_eq to the baryonic velocity and compute RMS residual against observed rotation curves.

## Current Performance
Median RMS ≈ 44.9 km/s on ~120 galaxies (binary age proxy).

## Limitations
- The mapping from microscopic trapping time to the macroscopic (a₀/a_bary)^β form is phenomenological.
- Limited real stellar age data for SPARC prevents a fully time-dependent accumulation model.
- Not a complete relativistic derivation; exploratory link between energy trapping and apparent mass.
