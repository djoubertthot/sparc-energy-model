# SPARC Energy-Mass Model

## Theory
The model explores the idea that apparent dynamical mass discrepancies in galaxies can be related to energy accumulation from radiation (photons, neutrinos, etc.) trapped or traversing baryonic matter over long timescales.

Key idea: 
\[
\rho_{\rm eq}(r) = \rho_0 \times \left(\frac{a_0}{a_{\rm bary}(r)}\right)^\beta
\]

where the equivalent density amplifies the baryonic contribution in a way that mimics dark matter effects, anchored to known physical scales (a0 related to cosmological acceleration).

## Parameters
- ρ0 ≈ 2×10^{-24} kg/m³ (single free global parameter)
- β ≈ 1.4
- a0 = 1.2×10^{-10} m/s²
- Y_star ≈ 0.5–0.6
- Age proxy: binary morphological factor (early-type 1.2, late-type 0.8)

## Methodology
1. Load SPARC rotation curves (*_rotmod.dat)
2. Compute baryonic acceleration a_bary(r)
3. Apply radial density amplification
4. Integrate cumulative equivalent mass
5. Derive model velocity and compute RMS residual

## Results
Median RMS ≈ 44.9 km/s on ~120 SPARC galaxies with the binary age proxy.
Comparable to simple phenomenological models; further improvement possible with better age estimates.

## Usage
See README.md for installation and running instructions.

## Limitations
- Limited real stellar ages available for SPARC
- Large external catalogs have poor automatic overlap
- Model is exploratory and not a complete physical theory
