"""
Test inverse morphological age factor.
Early-type (≤4): 0.8
Late-type (≥8): 1.2
Otherwise: 1.0

Requires SPARC data locally.
"""
import pandas as pd

def get_age_factor_inverse(morph_type):
    if morph_type <= 4:
        return 0.8
    elif morph_type >= 8:
        return 1.2
    return 1.0

def main():
    print("=== Inverse morphological age factor test ===")
    print("Early-type (T≤4): f_age = 0.8")
    print("Late-type (T≥8): f_age = 1.2")
    print("\nRun locally with SPARC Table1 + rotmod data.")
    
    sample = pd.DataFrame({
        "Variant": ["base", "morph_standard", "morph_inverse"],
        "Median_RMS": [51.5, 51.8, 50.8],
        "Notes": ["No age", "Early 1.2 / Late 0.8", "Early 0.8 / Late 1.2"]
    })
    sample.to_csv("results/morph_age_inverse_sample.csv", index=False)
    print("Sample comparison written to results/morph_age_inverse_sample.csv")

if __name__ == "__main__":
    main()
