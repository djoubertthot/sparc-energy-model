"""
Test different physical extensions of the energy-mass model.
Requires SPARC rotmod files in data/rotmod/.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from radial_model import load_rotmod, process_directory

def main():
    print("=== Variant Testing (requires local SPARC data) ===\n")
    print("Variants to be implemented:")
    print("- base")
    print("- base + age")
    print("- base + kinetic")
    print("- base + age + kinetic")
    print("- combinations with surface density")
    print("\nRun this script after placing *_rotmod.dat files in data/rotmod/")
    print("Results will be saved to results/variants_comparison.csv")

    # Placeholder: actual computation goes here once data is available
    results = pd.DataFrame({
        "Variant": ["base", "base+age", "base+kin", "base+age+kin"],
        "Median_RMS": [45.0, 44.5, 44.8, 44.2],
        "Notes": ["Reference", "Age proxy", "Kinetic term", "Combined"]
    })
    results.to_csv("results/variants_comparison_sample.csv", index=False)
    print("\nSample comparison table written to results/variants_comparison_sample.csv")

if __name__ == "__main__":
    main()
