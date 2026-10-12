# a₀ mesuré sur SPARC et les expressions candidates

Mesure : `scripts/fit_a0.py`, résultats `results/fit_a0.md`. Prédictions en redshift :
`docs/A0_ORIGIN.md`.

## 1. La mesure

Loi g_obs = ν(g_bar/a₀) g_bar, ν(y) = 1/(1 − e^{−√y}), a₀ balayé globalement sur les
171 courbes, 2 795 points.

| Traitement de Υ★ | a₀ optimal, 10⁻¹⁰ m/s² | Intervalle statistique, erreurs renormalisées | χ²/N |
|---|---|---|---|
| Υ★ SED par galaxie, fixé | 1,06 | 1,05 à 1,06 | 47 |
| Υ★ = 0,5 fixé | 1,19 | 1,18 à 1,20 | 47 |
| Υ★ libre, prior 0,1 dex autour du SED | 1,09 | 1,09 à 1,10 | 14 |

L'erreur statistique est négligeable ; l'incertitude est entièrement systématique et
vient du rapport masse-luminosité : a₀ et Υ★ sont dégénérés, a₀ ∝ Υ★^{−1} environ dans le
régime profond. Résultat : **a₀ = 1,06 à 1,19 × 10⁻¹⁰ m/s²**, en accord avec Lelli et al.
2017, 1,20 ± 0,02 (statistique) ± 0,24 (systématique). Le χ²/N de 47 à Υ★ fixé signale
que la dispersion de la relation dépasse les erreurs de mesure déclarées, 0,146 dex contre
0,1 dex attendu ; avec Υ★ libre il tombe à 14.

## 2. Les expressions candidates

| Expression | H₀ = 70 | H₀ = 73 | Dans la bande 1,06 à 1,19 ? |
|---|---|---|---|
| c H₀ / 2π | 1,08 | 1,13 | oui |
| c H₀ / 6, Verlinde 2016 | 1,13 | 1,18 | oui |
| c √(Λ/3) / 2π | 0,91 | 0,94 | non, 15 % trop bas |
| c √(Λ/3) / 6 | 0,95 | 0,99 | non, 10 % trop bas |
| c √(Λ/3) / 5 | 1,14 | 1,19 | oui |

Avec le préfacteur 2π ou 6, les expressions liées à l'expansion H₀ tombent dans la bande
et celles liées au vide Λ tombent en dessous. Mais aucune théorie ne fixe le préfacteur à
mieux que l'ordre de grandeur, et un préfacteur 1/5 remet la version Λ dans la bande. La
valeur seule ne départage donc pas H_A et H_B ; c'est la dépendance temporelle qui le fait.
Le test point par point sur 380 galaxies à z = 0,6-4,6 (`scripts/a0_highz.py`,
`docs/A0_ORIGIN.md` § 4-5) donne a₀(z)/a₀(0) = 1,0 (0,7-1,5) à z ≈ 2,2 et 0,5-0,9 à z ≈ 4,5,
là où H_A exige 3,3 et 7 ; mais ce test d'amplitude bascule pour 0,1 à 0,2 dex de
surestimation des masses baryoniques, et l'évolution revendiquée par MUSE-DARK III à z < 1,4
est le gaz moléculaire absent de son bilan. Le test de forme des courbes, indépendant des
masses (`scripts/a0_highz_shape2.py`), exclut H_A à 4,6σ sur douze courbes en gaz froid à
z = 1,5-4,8 et à 5,2σ avec les dix-neuf courbes ALPAKA, et donne n = 0,05 (−0,25 à 0,30) dans
a₀ ∝ H(z)ⁿ. Le tableau cohérent est donc : a₀ ≈ 0,2 c √(Λ/3), constant, avec un préfacteur à
expliquer ; a₀ ne suit pas l'expansion.

## 3. Les trois dérivations existantes

**Milgrom 1999, températures.** Un corps d'accélération a voit un bain d'Unruh à
T_U = ħa/(2πck) ; l'horizon de de Sitter rayonne à T_dS = ħ√(Λ/3)/(2πk). Pour a = a₀ :
T_U = 4,9·10⁻³¹ K, T_dS = 2,3·10⁻³⁰ K, rapport 0,21. L'argument dit que l'inertie change
quand T_U devient inférieure à T_dS, c'est-à-dire pour a < c√(Λ/3) = 5,7·10⁻¹⁰ ; il
place a₀ un facteur 5 en dessous sans dériver ce facteur.

**Verlinde 2016, gravité émergente.** L'entropie d'intrication de l'espace de de Sitter,
déplacée par la matière, produit une gravité supplémentaire g_D² = (cH₀/6) g_B. Dans le
régime profond cela équivaut à a₀ = cH₀/6 = 1,13 à 1,18, dans la bande. Mais la relation
complète de Verlinde n'est pas la relation d'accélération : Lelli et al. 2017 trouvent une
tension dans le régime de transition, et la théorie ne prédit pas de fonction ν.

**Smolin 2017, gravité quantique.** Une inertie modifiée dérivée d'un couplage entre la
longueur de Planck et Λ ; a₀ ∝ c√Λ avec un coefficient d'ordre unité non fixé.

## 4. Deux mesures supplémentaires (`scripts/measure_nu.py`, `scripts/a0_universality.py`)

**La forme de ν.** Relation binnée en g_bar, 16 bins, 2 795 points. Pente mesurée dans le
régime profond, g_bar < 0,1 a₀ : 0,498, pour 0,500 attendu. Écart quadratique des bins aux
fonctions candidates : McGaugh 0,029 dex, famille n à n = 1,4 0,038, Milgrom standard 0,048,
Famaey-Binney simple 0,034, Verlinde 0,066 avec son échelle fixée à cH₀. La relation de
Verlinde, sans paramètre, est la seule nettement rejetée par la forme ; les autres se
distinguent par moins de 0,02 dex, ce qui est la précision actuelle sur ν. Une théorie
doit prédire ν à mieux que cela pour être testée par la forme.

**L'universalité de a₀.** Ajusté galaxie par galaxie avec Υ★, distance et inclinaison
marginalisés sur des priors gaussiens : médiane 0,93·10⁻¹⁰, moyenne pondérée 1,23,
dispersion observée 0,43 dex, intervalles individuels 0,13 dex, dispersion intrinsèque
0,33 dex ; χ² = 6 842 pour 165 degrés de liberté si a₀ est universel ; 28 % des galaxies
excluent 1,2·10⁻¹⁰ d'un facteur 2. Corrélations : V_max +0,20, drapeau de qualité −0,22.
C'est le résultat de Rodrigues et al. 2018, contesté par Kroupa et al. 2018 et Li et al.
2018 au motif que les erreurs de SPARC ne sont pas gaussiennes. La corrélation avec la
qualité va dans leur sens ; celle avec V_max reste à expliquer. Le point est ouvert, et il
est décisif : une loi exige a₀ constant, un halo ne l'exige pas.

**Avec les distances de Cosmicflows-4** (`scripts/a0_universality_cf4.py`) : 44 galaxies
reçoivent une distance indépendante de la rotation, TRGB pour la plupart, plus précise que
celle de SPARC. Sur les 22 de qualité 1, la dispersion intrinsèque de a₀ passe de 0,25 à
0,19 dex, le χ² d'un a₀ universel de 944 à 496, la fraction de galaxies excluant 1,2·10⁻¹⁰
d'un facteur 2 de 23 % à 14 %. Les deux galaxies les plus discordantes, NGC 2903 et
UGC 7603, à 3,7 et 3,2·10⁻¹⁰ avec les distances de flot de Hubble, reviennent à 1,3 et 1,6
avec leurs distances TRGB. Le décalage de a₀ par galaxie anticorrèle avec le décalage de
distance, ρ = −0,64, comme a₀ ∝ 1/D l'impose en régime profond. Une partie de la
non-universalité apparente est donc un effet de distance, le tiers environ ; le reste,
0,19 dex, attend un modèle d'erreur sur les inclinaisons et les mouvements non circulaires
que les priors gaussiens ne capturent pas.

**Avec des cinématiques 3D** (`scripts/mhongoose_rings.py`, `scripts/littlethings3d_a0.py`) :
- MHONGOOSE, trois galaxies profondes communes avec SPARC. La géométrie retrouvée sur les
  champs de vitesse MeerKAT confirme SPARC pour DDO 161 et NGC 7793, et donne 42° pour
  ESO 444-G084 contre 32° dans SPARC, en accord avec le modèle 3D de Namumba et al. 2025,
  49°. Cette seule révision fait passer a₀ de 3,5 à 0,3-0,4·10⁻¹⁰ pour cette galaxie.
  Les cartes profondes prolongent les courbes à 1,3-1,7 R_dernier, jusqu'à 0,007 a₀ : DDO 161
  et ESO 444-G084 restent sur la relation ; NGC 7793 continue de décliner. Avec la géométrie
  libre par anneau (`scripts/ngc7793_warp.py` : V_rot, V_rad, PA, i), la loi tient à l'intérieur
  de R_dernier (résidu +0,02 dex, inclinaison oscillant de 40 à 52° entre anneaux voisins, terme
  non circulaire 10 %) ; entre 8 et 9,6 kpc le déclin subsiste (−0,35 dex, la loi exigerait
  i ≈ 33° contre 48-60° ajustés) mais le terme non circulaire monte à 20-40 % et le rms double ;
  au-delà de 10 kpc le terme radial égale ou dépasse la rotation : ce gaz ne mesure pas le
  potentiel. Le modèle 3D du cube MHONGOOSE à 26″ (`scripts/ngc7793_cube3d.py`,
  `results/ngc7793_cube3d.md`) tranche : entre 8 et 10 kpc la vitesse circulaire vaut 71-86 km/s
  contre 101 pour la loi (−0,15 à −0,31 dex), et la courbe de la loi, ajustée avec inclinaison,
  angle de position, dispersion et terme radial libres, ou avec i = 33° imposé, reste 12 à 16 %
  au-dessus du meilleur modèle. Ce n'est pas un effet de champ externe (champ 2M++ de NGC 7793 au
  17ᵉ centile de SPARC, et ESO 444-G084, dans un champ plus fort, sans déficit à 0,007 a₀). Le
  déficit est réel, propre à cette galaxie, sur un anneau dont le gaz porte 20-30 % de mouvements
  radiaux : défaut de la loi ou gaz hors d'équilibre, le second favorisé sans être prouvé.
- LITTLE THINGS en 3D, Iorio et al. 2017, trois naines communes. Inclinaisons à 1-4°,
  dérive asymétrique de 9 à 20 km/s au dernier point. Les vitesses 3D reproduisent SPARC à
  1-3 %, et a₀ vaut 1,06, 0,56 et 0,63·10⁻¹⁰ pour DDO 154, DDO 168 et NGC 2366, avec χ²/N
  de 1 à 3. Deux de ces galaxies restent un facteur 2 sous la valeur de l'échantillon avec
  la meilleure cinématique disponible : là, la géométrie n'est pas en cause.

**Au-delà des disques, lentille faible à 50-300 kpc** (`scripts/lensing_vflat.py`) : Mistele et
al. 2024 déprojettent le signal de lentille KiDS-1000 autour des galaxies isolées en vitesses
circulaires, plates de 50 à 300 kpc au moins et sans déclin visible jusqu'à 1 Mpc, sur la
relation de Tully-Fisher baryonique cinématique. Les halos ajustés sur SPARC, extrapolés sans
liberté nouvelle, donnent en médiane V(300 kpc)/V_plat = 0,68 (NFW libre), 0,66 (NFW avec
priors ΛCDM), 0,57 (DC14), 0,46 (Burkert), contre 1,00 ± 0,05 mesuré ; la loi isolée donne 0,94.
La loi avec effet de champ externe, formule 1D de Chae et al. 2020, ne garde la platitude à 10 %
près que pour e_N < 4·10⁻⁵ a₀ ; mais cette formule surestime l'effet. Résolue avec l'équation
de champ QUMOND (`scripts/efe_qumond_lensing.py`, `results/efe_qumond_lensing.md`), l'effet
n'apparaît qu'au-delà de r_t = √(G M_b/e_N), soit 300 à 1200 kpc pour le champ 2M++ des
galaxies SPARC (1,5·10⁻⁴ a₀) : V(300 kpc)/V_plat = 0,87, 1,00 et 1,01 pour M_b = 10¹⁰, 5·10¹⁰
et 2·10¹¹, et le déclin se produit entre 0,5 et 1 Mpc. La platitude à 300 kpc n'exclut donc pas
une théorie de champ de type MOND ; elle la teste à 0,5-1 Mpc, où les données de Mistele et al.
ne voient pas de déclin mais perdent leur fiabilité. Les deux lectures de la loi, dynamique
modifiée avec son effet de champ externe ou composante sombre sans lui, ne sont pas séparées par
les données actuelles (correction du 2026-10-12 d'une conclusion antérieure). Le rayon viriel médian
des ajustements NFW est 165 kpc : les lentilles traversent le rayon viriel sans voir le déclin
en √(ln r / r) qu'un halo ΛCDM y impose.

Bilan : la dispersion de a₀ par galaxie se décompose en distances, un tiers ; inclinaisons,
dominantes dans quelques cas comme ESO 444-G084 ; et un résidu de l'ordre de 0,2-0,3 dex
dans des naines bien mesurées, qui équivaut à 0,1-0,15 dex en g_obs, soit la dispersion de
la relation elle-même. La non-universalité de a₀ et la dispersion de la relation
d'accélération sont la même quantité vue de deux côtés ; savoir si elle est intrinsèque
demande les bilans de masse de gaz à 10 % et des échantillons 3D de la taille de BIG-SPARC.

## 5. Ce qu'une dérivation doit livrer

1. Le préfacteur numérique devant c√(Λ/3) ou cH₀, à 10 % près, car la mesure l'exige.
2. La fonction d'interpolation ν entière, pas seulement le régime profond : la relation
   d'accélération est mesurée de 0,003 a₀ à 30 a₀ et sa forme est une donnée.
3. Le signe et l'amplitude de l'effet de champ externe, absent dans une lecture
   « composante sombre », présent dans une lecture « inertie ». Non détecté dans les disques,
   où l'équation de champ ne prédit rien de mesurable à e_N ≈ 10⁻⁴ a₀ ; non testé par la lentille
   à 300 kpc, où QUMOND prédit moins de 13 % de déclin ; décidé par la lentille à 0,5-1 Mpc des
   lentilles de faible masse, où QUMOND prédit 30 à 50 % de déclin.
4. La dépendance temporelle, H(z) ou Λ, testable à z = 1 à 2.

Aucune des trois dérivations ne livre les quatre. C'est l'état exact du problème ouvert.
