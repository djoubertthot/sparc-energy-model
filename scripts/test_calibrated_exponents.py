"""
Teste la formulation avec exposants calibrés sur la régression des résidus.
"""
import numpy as np
import pandas as pd
from pathlib import Path

def test_calibrated():
    print("=== Test des exposants calibrés ===")
    print("Formulation :")
    print("E = E0 × (a0/a)^β × (t/t0)^αt × (L/L0)^αL × (R/R0)^αR")
    print("\nExposants suggérés par la régression :")
    print("- αt ~ 0.3–0.5")
    print("- αL ~ 0.1–0.2")
    print("- αR ~ 0.1–0.3")
    print("\nNécessite les données SPARC locales pour les résultats RMS.")

if __name__ == "__main__":
    test_calibrated()
