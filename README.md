# sparc-energy-model

Analyse des courbes de rotation SPARC et de leurs prolongements à haut redshift. État au
12 octobre 2026.

Le dépôt a commencé comme un modèle d'« accumulation d'énergie » censé expliquer la masse
manquante des galaxies. L'audit du 10 octobre 2026 a montré que ce modèle n'avait ni le
budget d'énergie (sept ordres de grandeur manquants) ni la forme fonctionnelle requise, et
que ce qui marchait dans ses ajustements était la relation d'accélération
g_obs = ν(g_bar/a₀) g_bar, sans paramètre par galaxie. Le dépôt teste maintenant cette loi,
ce qu'elle implique pour les halos, et l'origine de son échelle a₀. L'ancien matériel est
dans `archive_v1/` ; la release 1.0 en reste l'état de référence et n'est pas touchée.

## Résultats en bref

- **La loi tient** sur 171 galaxies SPARC avec un RMS de 13,7 km/s et zéro paramètre par
  galaxie, 9,2 km/s avec le rapport M/L libre dans sa marge SED (`results/model_v2.md`).
- **Les halos qu'elle impose** ont un cœur lié à la taille du disque, une densité centrale
  constante, une dispersion de concentration 2,4 fois la cosmologique, et des vitesses
  circulaires plates jusqu'à 300 kpc là où tout halo ajusté sur les disques décline de 30 à
  50 % (`paper/SPARC_HALOS_DRAFT.md`, `results/halo_*.md`, `results/lensing_vflat.md`).
- **L'effet de champ externe** n'est pas détecté dans les disques, où l'équation de champ ne
  prédit rien de mesurable ; la lentille à 300 kpc ne le teste pas non plus (moins de 13 % de
  déclin prédit par QUMOND) ; la lentille à 0,3-1 Mpc le teste et ne le voit pas : les lentilles
  de 10¹⁰ M☉ restent plates là où QUMOND avec le champ 2M++ prédit 30 à 50 % de déclin,
  e_N < 6·10⁻⁵ a₀, champ 2M++ exclu à ~3σ ; dans 2M++, les galaxies isolées au sens de Brouwer
  et al. gardent un champ de 1,1·10⁻⁴ a₀ (`results/test_efe.md`, `results/efe_qumond_lensing.md`,
  `results/efe_lensing_mpc.md`, `results/efe_isolated_2mpp.md`).
- **Les binaires larges de Gaia** : de 1000 à 0,1 a₀ le rapport vitesse relative sur vitesse
  circulaire newtonienne reste à 1,00-1,06 de sa calibration, là où une loi locale sans effet de
  champ externe prédit ×1,36 à ×1,66 (exclue à 15σ) et QUMOND avec le champ galactique ×1,17
  (3 à 5σ au-dessus des données, systématiques de même taille) (`results/wide_binaries.md`).
- **Les amas** : avec le a₀ des galaxies, la loi manque un facteur 4-7 en masse entre 15 et
  600 kpc (CLASH), 2 à R500 et 1,3 à R200 (X-COP) ; le « facteur 2 » est la valeur à R500, et
  la masse manquante est concentrée dans les cœurs (`results/cluster_rar.md`).
- **a₀ ne suit pas l'expansion.** Sur 380 galaxies à z = 0,6-5,7 en amplitude et 31 courbes
  en gaz froid à z = 1,5-4,8 en forme, a₀ ∝ H(z)ⁿ donne n = 0,05 (−0,25 à 0,30) ; a₀ est
  constant, l'échelle est celle de Λ (`docs/A0_ORIGIN.md`, `results/a0_highz*.md`).

## Installation

```bash
git clone https://github.com/djoubertthot/sparc-energy-model.git
cd sparc-energy-model
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

## Données

```bash
python scripts/download_catalogs.py   # SPARC, modèles de masse, catalogues auxiliaires
python scripts/enrich_sparc.py        # table croisée data/sparc_sample/sparc_enriched.csv
```

Les sources de la table enrichie sont décrites dans `docs/DATA_SOURCES.md`. Les données
externes utilisées par les tests (cartes de moments MHONGOOSE, LITTLE THINGS 3D, catalogue
2M++, distances CF4, tables à haut redshift, vitesses de lentille KiDS, profils d'amas, binaires
larges Gaia) sont dans `data/external/` et `data/catalogs/` avec, pour le haut redshift, la
lentille, les amas et les binaires, un README qui donne la provenance de chaque table et signale
celles qui sont des lectures de figures.

## Reproduire

Chaque script écrit sa table dans `results/` (fichier `.md` du même nom) et, le cas échéant,
une figure et un CSV par galaxie.

| Script | Produit |
|---|---|
| `model_v2.py` | la loi sur SPARC, RMS, biais, composante sombre équivalente |
| `fit_a0.py`, `measure_nu.py` | valeur de a₀, forme de la fonction ν |
| `fit_halos.py`, `fit_halos_priors.py`, `fit_dc14.py`, `halo_relations.py` | halos NFW, Burkert, DC14, priors ΛCDM, relations d'échelle |
| `test_efe.py`, `test_diversity.py` | effet de champ externe (catalogue 2M++), diversité à 2 kpc |
| `a0_universality.py`, `a0_universality_cf4.py` | a₀ par galaxie, distances CF4 |
| `mhongoose_rings.py`, `ngc7793_warp.py`, `littlethings3d_a0.py` | cinématique HI profonde et 3D |
| `ngc7793_cube_cut.py`, `ngc7793_cube3d.py` | modèle 3D du cube MHONGOOSE de NGC 7793 (anneaux, faisceau, modèles de déclin, gauchissement, loi) |
| `lensing_vflat.py` | halos extrapolés à 50-300 kpc contre la lentille faible |
| `efe_qumond_lensing.py`, `efe_lensing_mpc.py` | effet de champ externe QUMOND : prédiction, puis test bin par bin sur les vitesses de lentille à 0,3-2,5 Mpc |
| `efe_isolated_2mpp.py` | champ externe des galaxies isolées dans 2M++ |
| `parse_cluster_tables.py`, `cluster_rar.py` | la loi dans les amas CLASH et X-COP, facteur manquant par rayon |
| `wide_binaries.py` | binaires larges Gaia EDR3 sous a₀ : loi sans effet de champ externe, QUMOND, Newton |
| `test_desitter_nu.py` | fonction ν du vide de de Sitter (Milgrom 1999) sur SPARC |
| `a0_highz.py`, `a0_highz_coldgas.py`, `a0_highz_shape.py`, `a0_highz_shape2.py` | a₀ à haut redshift, amplitude puis forme des courbes |
| `final_model.py`, `final_baryonic_model.py`, `run_rms_tests*.py`, `optimize_constraints.py`, `temporary_mass_budget.py` | l'audit du modèle 1.0 |

## Documents

- `docs/THEORY.md` : la loi, ses deux lectures, l'échelle a₀ et ce qui tuerait la théorie.
- `docs/A0_ORIGIN.md`, `docs/A0_CANDIDATES.md` : chantier sur l'origine de a₀, verdict sur
  sa dépendance temporelle.
- `paper/SPARC_HALOS_DRAFT.md` : brouillon d'article sur les halos imposés par SPARC.
- `archive_v1/` : le modèle abandonné, avec son propre README.
