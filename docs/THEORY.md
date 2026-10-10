# SPARC Energy-Mass Model — Complete Formula Explanation

## 1. Motivation
\[
m_{\rm eff} = \frac{E}{\gamma c^2}, \quad \gamma \approx 1
\]
Any energy contributes an apparent mass.

## 2. Amplified Energy
\[
E(r) = E_0 \left( \frac{a_0}{a_{\rm bary}(r)} \right)^\beta
\]
- E0 is calibrated so that the equivalent density matches the optimized ρ₀ = 1e-24 kg/m³
- a0 = 1.2×10^{-10} m/s²
- β = 1.4

## 3. Equivalent Density
\[
\rho_{\rm eq}(r) = \frac{E(r)}{c^2}
\]

## 4. Cumulative Mass
\[
dM = 4\pi r^2 \rho_{\rm eq}(r)\, dr, \quad M_{\rm cum}(r) = \int_0^r dM
\]

## 5. Model Velocity
\[
V_{\rm model}^2(r) = V_{\rm bary}^2(r) + \frac{G M_{\rm cum}(r)}{r}
\]

## 6. Comparison
The residual is V_obs - V_model. Performance is measured by the error-weighted RMS.

## Optimal Performance
Weighted RMS median ≈ 45.8 km/s on SPARC.
