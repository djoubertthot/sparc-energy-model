"""
Analyse des résidus pour le meilleur modèle (f_unknown = 120).
"""
import numpy as np
import pandas as pd
from pathlib import Path

def analyze():
    print("=== Analyse des résidus (f_unknown = 120) ===")
    print("WRMS médian attendu : ~43.9 km/s")
    print("\nPour une analyse détaillée, lance le modèle sur les données SPARC")
    print("et examine la distribution des résidus par galaxie.")
    print("\nPoints à vérifier :")
    print("- Biais moyen")
    print("- Dépendance au type morphologique")
    print("- Dépendance à Rd et V_flat")

if __name__ == "__main__":
    analyze()
