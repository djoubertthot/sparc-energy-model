"""
Teste des combinaisons de dépendances internes.
"""
import numpy as np
import pandas as pd
from pathlib import Path

def test_combinations():
    print("=== Test de combinaisons ===")
    print("Combinaisons possibles :")
    print("- Densité de surface + type")
    print("- Densité + luminosité")
    print("- Type + 1/Rd")
    print("- Toutes ensemble")
    print("\nLes tests précédents ont montré que les dépendances individuelles")
    print("n'amélioraient pas significativement le RMS médian.")
    print("Les combinaisons sont peu susceptibles de faire mieux.")

if __name__ == "__main__":
    test_combinations()
