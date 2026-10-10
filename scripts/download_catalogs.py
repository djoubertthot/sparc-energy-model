"""
Script pour télécharger les catalogues utiles.
Note : les catalogues complets (MaNGA, CALIFA) sont très volumineux.
"""
import urllib.request
import zipfile
import os
from pathlib import Path

def download_sparc():
    """Télécharge les fichiers SPARC de base"""
    url = "http://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt"
    print("Téléchargement de la table SPARC...")
    # Exemple simplifié — adapter selon les besoins
    print("Visite http://astroweb.case.edu/SPARC/ pour télécharger les fichiers rotmod manuellement.")
    print("Ou utilise le CDS : https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/152/157")

def note_other_catalogs():
    print("\n=== Autres catalogues ===")
    print("CALIFA / PyCASSO : http://pycasso.ufsc.br/")
    print("MaNGA FIREFLY : https://www.sdss.org/dr17/manga/manga-data/manga-firefly-value-added-catalog/")
    print("Ces catalogues sont volumineux (plusieurs Go) et ont peu d'overlap avec SPARC.")
    print("Téléchargement manuel recommandé.")

if __name__ == "__main__":
    Path("data").mkdir(exist_ok=True)
    download_sparc()
    note_other_catalogs()
