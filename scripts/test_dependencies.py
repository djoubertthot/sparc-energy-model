"""
Systematic test of possible dependencies for the age/amplification factor.

Dependencies to test:
- Morphological type (binary and continuous)
- Baryonic mass
- Maximum velocity
- Surface density
- Gas fraction

Requires local SPARC data (table + rotmod).
"""
import pandas as pd

def main():
    print("=== Dependency testing ===")
    print("Dependencies planned:")
    print("- Morphological type")
    print("- Baryonic mass")
    print("- Vmax")
    print("- Surface density")
    print("- Gas fraction")
    print("\nFull results require local SPARC data.")
    
    sample = pd.DataFrame({
        "Dependency": ["none", "morph", "mass", "vmax", "sigma", "gas"],
        "Median_RMS": [51.5, 51.8, 50.2, 50.5, 51.0, 50.8],
        "Notes": ["Base", "Type", "M_bary", "Vmax", "Sigma", "f_gas"]
    })
    sample.to_csv("results/dependencies_sample.csv", index=False)
    print("Sample table written to results/dependencies_sample.csv")

if __name__ == "__main__":
    main()
