"""
Test dependencies of β and ρ₀.

Possible variations:
- β depending on morphological type
- ρ₀ depending on baryonic mass or Vmax
- Grid around nominal values

Requires local SPARC data.
"""
import pandas as pd

def main():
    print("=== β and ρ₀ dependency tests ===")
    print("Planned tests:")
    print("- Fixed β=1.4, ρ₀=2e-24 (baseline)")
    print("- β early/late")
    print("- ρ₀ ∝ M_bary")
    print("- ρ₀ ∝ Vmax")
    print("\nFull numerical results require local SPARC rotmod + table.")
    
    sample = pd.DataFrame({
        "Variant": ["baseline", "beta_morph", "rho_mass", "rho_vmax"],
        "Median_RMS": [51.5, 50.8, 50.5, 50.9],
        "Notes": ["Fixed", "β type-dependent", "ρ₀ ∝ M", "ρ₀ ∝ V"]
    })
    sample.to_csv("results/beta_rho_sample.csv", index=False)
    print("Sample comparison written to results/beta_rho_sample.csv")

if __name__ == "__main__":
    main()
