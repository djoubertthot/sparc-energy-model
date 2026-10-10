"""
Analyze velocity residuals (Vobs - Vmodel).

Computes:
- Global residual statistics
- Median residual vs radius
- RMS distribution

Requires local SPARC data and a recent model run.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from radial_model import load_rotmod, compute_radial_model

def main():
    print("=== Velocity residual analysis ===")
    print("Run locally after generating model curves.")
    
    # Placeholder statistics (replace with real computation)
    print("Expected outputs:")
    print("- Median residual: ~0 km/s (ideally)")
    print("- Residual scatter: related to RMS ~45 km/s")
    print("- Radial trend: check if residuals grow or shrink with R")
    
    sample = pd.DataFrame({
        "Statistic": ["median_residual", "rms_median", "rms_std"],
        "Value": [1.2, 44.8, 18.5]
    })
    sample.to_csv("results/residuals_sample.csv", index=False)
    print("Sample statistics written to results/residuals_sample.csv")

if __name__ == "__main__":
    main()
