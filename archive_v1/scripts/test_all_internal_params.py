"""
Test dependencies of source factors on internal galaxy parameters.

Possible parameters:
- Surface density
- Morphological type
- Luminosity
- Gas fraction
- etc.

Requires local SPARC table + additional catalogs.
"""
import pandas as pd

def main():
    print("=== Test all internal parameters ===")
    print("Parameters to test: surface density, type, luminosity, gas fraction")
    print("Full results require local data.")
    
    sample = pd.DataFrame({
        "Parameter": ["surface_density", "morph_type", "luminosity"],
        "Expected_impact": ["mild", "small", "small"]
    })
    sample.to_csv("results/internal_params_sample.csv", index=False)
    print("Sample written.")

if __name__ == "__main__":
    main()
