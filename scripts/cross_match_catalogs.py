"""
Script de départ pour croiser SPARC avec des catalogues externes.
"""
import pandas as pd
from pathlib import Path

def cross_match():
    print("=== Croisement de catalogues ===")
    print("Catalogues cibles :")
    print("- GSWLC (âges stellaires)")
    print("- z0MGS Dust (masse de poussière)")
    print("- SDSS MPA-JHU (métallicité)")
    print("\nProcédure :")
    print("1. Télécharger les catalogues")
    print("2. Matcher par coordonnées ou noms")
    print("3. Fusionner avec la table SPARC")
    print("\nNécessite les fichiers locaux.")

if __name__ == "__main__":
    cross_match()
