"""
Analyse des résidus pour le modèle avec type + Rd.
"""
import numpy as np
import pandas as pd
from pathlib import Path

def analyze_residuals(results_csv="results/rms_with_type_rd.csv"):
    try:
        df = pd.read_csv(results_csv)
        print("=== Analyse des résidus ===")
        print(f"WRMS médian : {df['WRMS'].median():.1f} km/s")
        print(f"WRMS moyen : {df['WRMS'].mean():.1f} km/s")
        print(f"Écart-type : {df['WRMS'].std():.1f} km/s")
        print("\nGalaxies avec plus grand RMS :")
        print(df.nlargest(5, 'WRMS')[['Galaxy', 'WRMS']].to_string(index=False))
    except FileNotFoundError:
        print("Fichier de résultats introuvable.")
        print("Lance d'abord le modèle avec les données SPARC.")

if __name__ == "__main__":
    analyze_residuals()
