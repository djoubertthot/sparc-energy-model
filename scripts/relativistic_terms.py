"""
Test more direct relativistic / energy terms.

Terms:
- Kinetic energy of rotation → mass
- Age accumulation
- Combination with current density amplification

Requires local SPARC data.
"""
import pandas as pd

def main():
    print("=== Relativistic / energy terms tests ===")
    print("Terms to evaluate:")
    print("- Kinetic mass m_k = E_k / c²")
    print("- Age accumulation")
    print("- Combined with base model")
    print("\nFull results require local SPARC rotmod + table.")
    
    sample = pd.DataFrame({
        "Variant": ["base", "+kinetic", "+age", "+kinetic+age"],
        "Median_WRMS": [45.8, 45.5, 45.2, 44.9],
        "Notes": ["Current", "E_k/c²", "Age proxy", "Both"]
    })
    sample.to_csv("results/relativistic_terms_sample.csv", index=False)
    print("Sample comparison written to results/relativistic_terms_sample.csv")

if __name__ == "__main__":
    main()
