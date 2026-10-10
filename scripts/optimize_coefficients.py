"""
Simple grid search for age and kinetic coefficients.
Requires SPARC data in data/rotmod/.
"""
import pandas as pd
import numpy as np
from radial_model import process_directory

def main():
    print("=== Coefficient optimization (run locally) ===\n")
    
    # Example grid (user can expand)
    age_factors = [0.8, 1.0, 1.2]
    kin_alphas = [0.0, 1e-6, 1e-5, 1e-4]
    
    results = []
    for age in age_factors:
        for alpha in kin_alphas:
            # Placeholder: replace with actual model call that accepts parameters
            median_rms = 45.0 + np.random.normal(0, 1)  # demo only
            results.append({
                "age_factor": age,
                "kin_alpha": alpha,
                "median_rms": median_rms
            })
    
    df = pd.DataFrame(results)
    df.to_csv("results/coefficient_search_sample.csv", index=False)
    print(df.sort_values("median_rms").head())
    print("\nFull search requires local SPARC data and parameter-aware model.")
    print("Best result will be written to results/coefficient_search.csv")

if __name__ == "__main__":
    main()
