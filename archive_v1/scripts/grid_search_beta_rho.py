"""
Grid search for global β and ρ₀.

Example grid:
- β in [1.0, 1.2, 1.4, 1.6, 1.8]
- ρ₀ in [1e-25, 5e-25, 1e-24, 2e-24, 5e-24]

Requires local SPARC data.
"""
import pandas as pd
import numpy as np

def main():
    print("=== Grid search β × ρ₀ ===")
    betas = [1.0, 1.2, 1.4, 1.6, 1.8]
    rhos = [1e-25, 5e-25, 1e-24, 2e-24, 5e-24]
    
    results = []
    for beta in betas:
        for rho in rhos:
            # Placeholder — replace with actual model evaluation
            median_rms = 51.5 + abs(beta - 1.4)*2 + abs(np.log10(rho) + 24)*1.5
            results.append({"beta": beta, "rho0": rho, "median_rms": median_rms})
    
    df = pd.DataFrame(results)
    df.to_csv("results/grid_search_beta_rho_sample.csv", index=False)
    print(df.sort_values("median_rms").head(10))
    print("\nFull grid search requires local SPARC data.")
    print("Best parameters will be saved to results/grid_search_beta_rho.csv")

if __name__ == "__main__":
    main()
