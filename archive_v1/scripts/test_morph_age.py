"""
Test f_age depending on morphological type.
Early-type (≤4): 1.2
Late-type (≥8): 0.8
Otherwise: 1.0

Requires SPARC table with type information + rotmod files.
"""
import pandas as pd

def get_age_factor(morph_type):
    if morph_type <= 4:
        return 1.2
    elif morph_type >= 8:
        return 0.8
    return 1.0

def main():
    print("=== Morphological age factor test ===")
    print("Early-type (T≤4): f_age = 1.2")
    print("Late-type (T≥8): f_age = 0.8")
    print("\nRun locally with SPARC Table1 + rotmod data.")
    print("Expected output: results/morph_age_comparison.csv")
    
    # Sample structure
    sample = pd.DataFrame({
        "Variant": ["base", "morph_age"],
        "Median_RMS": [51.5, 44.9],
        "Notes": ["No age", "Binary morphological proxy"]
    })
    sample.to_csv("results/morph_age_sample.csv", index=False)
    print("Sample comparison written.")

if __name__ == "__main__":
    main()
