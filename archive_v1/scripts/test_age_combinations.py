"""
Teste des combinaisons du facteur d'âge avec d'autres dépendances.
"""
import numpy as np
import pandas as pd
from pathlib import Path

def test_combinations():
    print("=== Combinaisons âge + autres dépendances ===")
    print("Combinaisons à tester :")
    print("- Âge + Rd")
    print("- Âge + type")
    print("- Âge + densité de surface")
    print("- Âge + Rd + type")
    print("\nNécessite les données SPARC locales.")
    print("Meilleur résultat précédent : type + Rd = 43.8 km/s")

if __name__ == "__main__":
    test_combinations()
