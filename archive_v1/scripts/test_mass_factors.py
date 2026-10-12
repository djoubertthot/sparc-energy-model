"""
Teste des facteurs dépendants de la masse pour améliorer le RMS.
"""
import numpy as np
import pandas as pd
from pathlib import Path

def test_mass_factors():
    print("=== Test de facteurs de masse ===")
    print("Facteurs à tester :")
    print("- (M / M0)^α")
    print("- (V / V0)^δ")
    print("- Combinaison avec âge + Rd")
    print("\nNécessite les données SPARC locales pour les résultats RMS.")
    print("Objectif : réduire le biais sur les galaxies massives.")

if __name__ == "__main__":
    test_mass_factors()
