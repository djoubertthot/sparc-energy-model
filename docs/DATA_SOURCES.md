# Sources de données supplémentaires

La table retenue est `data/sparc_sample/sparc_enriched.csv`, produite par
`scripts/download_catalogs.py` puis `scripts/enrich_sparc.py`.
Les valeurs SPARC d'origine restent dans les colonnes `sparc_*`.
Le journal des remplacements est `data/sparc_sample/adjustment_log.csv`.

## Catalogues utilisés
- SPARC (Lelli, McGaugh & Schombert 2016) : table, modèles de masse, courbes
- Marasco et al. 2025, A&A 695, L23 : masse stellaire SED pour 154 galaxies
- Duey et al. 2025, table WISE II du site SPARC : rapport M/L selon g−W1
- Leroy et al. 2019, z0MGS : distance, masse stellaire, SFR
- Nersesian et al. 2019, DustPedia THEMIS : masse de poussière et SFR
- De Vis et al. 2019 : métallicité 12+log(O/H), calibration PG16 S
- SIMBAD : coordonnées, absentes de la table SPARC

## Ce qui est ajusté
- M/L constant 0,5 : remplacé par la masse SED, sinon par le M/L WISE
- Distance Hubble-flow : l'incertitude est resserrée seulement si un catalogue
  est plus précis et couvre encore la valeur SPARC. La seule mise à jour
  est NGC 4559, 2,7 Mpc vers 0,83 Mpc, distance 9,0 Mpc inchangée.
- Vflat = 0 : laissé vide. Ni Lelli+2019 ni Marasco+2025 ne donnent de
  vitesse plate pour ces galaxies.

## Ce qui reste absent
- Âge stellaire en Gyr : aucun des catalogues publiés ne le donne.
  Le rapport M/L à 3,6 µm et le SFR en tiennent lieu.
- Coordonnées de 20 galaxies LSB (F, D et KK98) que SIMBAD ne résout pas.
- Masse stellaire de 18 galaxies sans SED exploitable, sans WISE et sans z0MGS.
- Poussière et métallicité hors de l'empreinte Herschel / DustPedia.
