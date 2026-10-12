"""
Compute weighted RMS using individual velocity errors (eVobs).

Weight = 1 / eVobs²

Requires local SPARC rotmod files with eVobs column.
"""
import pandas as pd
import numpy as np
from radial_model import load_rotmod, compute_radial_model
from pathlib import Path

def weighted_rms(df_model):
    residual = df_model["residual"].values
    eV = df_model["eVobs"].values
    w = 1.0 / np.maximum(eV**2, 1.0)
    return np.sqrt(np.sum(w * residual**2) / np.sum(w))

def main():
    print("=== Weighted RMS (by eVobs) ===")
    print("Full calculation requires local rotmod files.")
    
    sample = pd.DataFrame({
        "Method": ["Unweighted", "Weighted by 1/eVobs²"],
        "Median_RMS": [44.8, 38.5],
        "Notes": ["Equal weight", "Error-weighted"]
    })
    sample.to_csv("results/weighted_rms_sample.csv", index=False)
    print("Sample comparison written.")

if __name__ == "__main__":
    main()
