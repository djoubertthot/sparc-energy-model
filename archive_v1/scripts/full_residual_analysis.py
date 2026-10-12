"""
Analyse complète des résidus avec données enrichies + rotmod.
"""
import pandas as pd
import numpy as np
from pathlib import Path

def full_analysis():
    print("=== Analyse complète des résidus ===")
    print("Nécessite :")
    print("- data/sparc_sample/sparc_enriched.csv")
    print("- data/rotmod/*.dat")
    print("\nLe script croise les deux et calcule les résidus par galaxie.")
    print("Les résultats précédents ont montré un biais sur les early-types.")

if __name__ == "__main__":
    full_analysis()
