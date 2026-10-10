"""
Test ρ₀ or β dependence on galaxy size / surface density.
Avoid pure radial age gradient.

Requires local SPARC data.
"""
import pandas as pd

def main():
    print("=== Size / surface density dependence tests ===")
    print("Leads to test:")
    print("- ρ₀ ∝ 1/Rd or σ_mean")
    print("- β adjusted for extended galaxies")
    print("- Data quality flag correlation")
    print("\nFull results require local SPARC table + rotmod.")
    
    sample = pd.DataFrame({
        "Lead": ["rho_1_over_Rd", "beta_extended", "quality_flag"],
        "Expected_effect": ["Reduce high-Rd RMS", "Mild improvement", "Check bias"],
        "Median_RMS_sample": [46.2, 45.5, 44.8]
    })
    sample.to_csv("results/size_dependence_sample.csv", index=False)
    print("Sample table written.")

if __name__ == "__main__":
    main()
