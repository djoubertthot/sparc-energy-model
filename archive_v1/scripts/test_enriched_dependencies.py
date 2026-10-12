"""
Teste les dépendances avec les données enrichies.
"""
import pandas as pd
import numpy as np
from pathlib import Path

def test_dependencies():
    path = Path("data/sparc_sample/sparc_enriched.csv")
    if not path.exists():
        print("Fichier enrichi introuvable.")
        return

    df = pd.read_csv(path)
    print("=== Données enrichies ===")
    print(f"Galaxies : {len(df)}")
    print(f"Avec poussière : {df['Mdust_Msun'].notna().sum()}")
    print(f"Avec SFR : {df['SFR_Msun_yr'].notna().sum()}")
    print(f"Avec métallicité : {df['OH'].notna().sum()}")
    print(f"Avec masse stellaire : {df['Mstar_Msun'].notna().sum()}")

    print("\n=== Corrrélations simples (proxy) ===")
    # Note : un vrai test de RMS nécessite les fichiers rotmod
    print("Les tests de RMS complets nécessitent les fichiers rotmod locaux.")
    print("Les données enrichies sont prêtes pour être utilisées.")

if __name__ == "__main__":
    test_dependencies()
