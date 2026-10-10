# Planned Extensions and Testing Plan

## Physically Plausible Extensions to Test
1. **Kinetic energy accumulation**
   - Term ∝ M_bary × V_rot^{2} × t_age
   - Physically motivated by E = γ m c^{2} and long-term accumulation.

2. **Local Lorentz factor**
   - γ(v_rot) applied to a kinetic contribution.
   - Expected to be very small (v/c ≈ 10^{-3}).

3. **Age-weighted amplification**
   - Replace binary proxy with a continuous age factor when data become available.

4. **Combinations**
   - Kinetic term + current (a0/a_bary)^β
   - Age × kinetic term

## Less Promising (for now)
- Pure gravitational redshift / time dilation (too small)
- Frame-dragging (negligible)
- Uniform cosmic backgrounds without a trapping mechanism

## Testing Plan
1. Modify `compute_radial_model` to accept optional additive terms.
2. Generate RMS CSVs for each variant on the full SPARC sample.
3. Compare median RMS and scatter against the baseline.
4. Keep only variants that improve the fit without adding unconstrained free parameters.

## Current Limitation
Full SPARC rotmod files are not present in this environment. All numerical tests must be run locally after downloading the data.

## Next Step
A modular test script can be added to the repository so that variants can be switched on/off easily.
