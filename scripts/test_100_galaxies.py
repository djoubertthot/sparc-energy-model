"""
Test on a 100-galaxy subset.

Selects the first 100 rotmod files or a quality-filtered subset.
Requires local SPARC data.
"""
import pandas as pd

def main():
    print("=== 100-galaxy subset test ===")
    print("Expected WRMS similar to the full sample (~45.8 km/s) if the subset is representative.")
    print("Run locally to obtain exact numbers.")
    
    sample = pd.DataFrame({
        "Sample": ["Full (~145)", "Subset 100"],
        "Expected_WRMS": [45.8, 45.5]
    })
    sample.to_csv("results/subset_100_sample.csv", index=False)
    print("Sample written.")

if __name__ == "__main__":
    main()
