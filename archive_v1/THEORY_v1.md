# SPARC Energy-Mass Model

## Current Equation
\[
\rho_{\rm eq}(r) = \rho_0 \times \left( \frac{a_0}{a_{\rm bary}(r)} \right)^\beta
\]

The equivalent density is then radially integrated to obtain an additional dynamical mass contribution.

## Parameter Explanations

| Parameter | Optimal value | Meaning |
|-----------|---------------|---------|
| **ρ₀** | 1.0 × 10^{-24} kg/m³ | Reference energy/mass density scale. Sets the overall amplitude of the equivalent mass. Empirically calibrated on SPARC. |
| **β** | 1.4 | Exponent that controls how strongly the equivalent density rises when the baryonic acceleration drops. Higher β = stronger amplification in the outer disk. |
| **a₀** | 1.2 × 10^{-10} m/s² | Characteristic acceleration scale. Anchored to the cosmological acceleration (order of c × H₀). Below this scale the amplification becomes significant. |
| **Y★** | 0.5 | Stellar mass-to-light ratio. Converts the observed stellar light (disk + bulge) into baryonic mass. Standard value for SPARC. |

## Physical Motivation
- ρ₀ and β together determine how much apparent mass is added by trapped or traversing energy (via E = γ m c²).
- a₀ marks the transition between baryon-dominated and energy-amplification-dominated regimes.
- Y★ is a standard astrophysical parameter, not unique to this model.

## Performance
With the optimal values above, the median RMS on 145 SPARC galaxies is approximately 44.8 km/s.

## Notes
These parameters remain phenomenological. A first-principles derivation linking ρ₀ and β directly to photon/neutrino residence times has not been completed.
