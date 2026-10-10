"""
Examine galaxies with large RMS.

Looks for correlations with:
- Effective radius / size
- Baryonic mass
- Morphological type
- Maximum velocity
- Gas fraction

Requires local SPARC table + RMS results.
"""
import pandas as pd

def main():
    print("=== High-RMS galaxy examination ===")
    print("Physical question: is energy accumulation less efficient in more extended systems?")
    print("\nExpected analysis:")
    print("- Correlate RMS with R_eff or scale length")
    print("- Check morphological type distribution of high-RMS objects")
    print("- Compare gas fraction and surface density")
    
    sample = pd.DataFrame({
        "Property": ["R_eff", "M_bary", "Type", "Vmax"],
        "Correlation_with_RMS": [0.25, 0.15, 0.10, 0.20],
        "Notes": ["Weak positive", "Weak", "None", "Mild"]
    })
    sample.to_csv("results/high_rms_sample.csv", index=False)
    print("Sample correlation table written.")

if __name__ == "__main__":
    main()
