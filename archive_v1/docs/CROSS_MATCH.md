# Croisement de catalogues

## Scripts
1. `scripts/download_catalogs.py` télécharge SPARC, Marasco+2025, z0MGS, DustPedia et De Vis+2019.
2. `scripts/enrich_sparc.py` croise par nom, puis par coordonnées à 20 arcsec, et écrit `data/sparc_sample/sparc_enriched.csv`.

`scripts/cross_match_catalogs.py` appelle cette seconde étape.

## Règle
Une valeur SPARC n'est remplacée que si le catalogue externe donne une mesure plus précise qui exclut l'ancienne, ou si la colonne SPARC est vide. Les courbes de rotation restent celles du fichier `Rotmod_LTG.zip` tant que la distance adoptée ne change pas.
