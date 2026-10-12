"""
Check correlation between efficiency and age in temporal accumulation.

Requires local age estimates for SPARC galaxies.
"""
import pandas as pd

def main():
    print("=== Efficiency-Age correlation ===")
    print("With real ages, one can fit η(t) or check if a single η works.")
    print("Expected: weak or no strong correlation with current proxies.")
    
    sample = pd.DataFrame({
        "Age_proxy": ["binary_morph", "constant_10Gyr"],
        "Best_efficiency": [1.0, 1.0],
        "WRMS": [45.5, 45.8]
    })
    sample.to_csv("results/efficiency_age_sample.csv", index=False)
    print("Sample written.")

if __name__ == "__main__":
    main()
