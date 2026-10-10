import pandas as pd
import numpy as np

def main():
    df = pd.read_csv("results/rms_all_galaxies.csv")
    print("=== RMS Statistics ===")
    print(df["RMS"].describe())
    print("\nMedian RMS:", df["RMS"].median())
    print("\nTop 5 best:")
    print(df.nsmallest(5, "RMS")[["Galaxy", "RMS"]])
    print("\nTop 5 worst:")
    print(df.nlargest(5, "RMS")[["Galaxy", "RMS"]])
    # Simple correlation example if more columns exist
    if "Vobs_last" in df.columns and "Vmodel_last" in df.columns:
        corr = df[["Vobs_last", "Vmodel_last"]].corr().iloc[0,1]
        print(f"\nCorrelation Vobs_last vs Vmodel_last: {corr:.3f}")

if __name__ == "__main__":
    main()
