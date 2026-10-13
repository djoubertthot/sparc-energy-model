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
place a₀ un facteur 5 en dessous sans dériver ce facteur. Testé le 2026-10-12
(`scripts/test_desitter_nu.py`, `results/test_desitter_nu.md`) : l'inertie proportionnelle à
l'excès de température T(a) − T(0) avec T ∝ √(a² + a_Λ²) donne en forme fermée
ν(y) = √(1 + 1/y). Sur la relation binnée de SPARC, χ²/bins 8,65 et 0,044 dex de résidu,
entre la famille n = 1,4 (5,5 ; 0,038) et la RAR (17,7 ; 0,029), loin devant Verlinde
(73 ; 0,066) : la forme passe. L'échelle prédite, a₀ = 2a_Λ = 1,14·10⁻⁹, est 9,5 fois trop
grande : le vide fixe la forme et l'ordre de grandeur, pas le préfacteur.

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
une théorie de champ de type MOND ; elle la teste à 0,3-1 Mpc (correction du 2026-10-12 d'une
conclusion antérieure tirée de la formule 1D). Test fait le 2026-10-13 bin par bin sur la table
de Mistele et al. (`scripts/efe_lensing_mpc.py`, `results/efe_lensing_mpc.md`), amplitude libre
par bin pour que seule la forme au-delà de r_t compte : le bin de plus faible masse
(log M_b = 10,1, r_t = 300 kpc pour le champ 2M++) reste plat jusqu'à 1 Mpc et donne
e_N < 5·10⁻⁵ a₀ à 2σ, le champ 2M++ exclu à Δχ² = 14 dans ce bin seul ; les trois bins plus
massifs ont leur r_t à 0,7-1,5 Mpc, là où les données s'arrêtent, et acceptent aussi bien 0 que
1,5·10⁻⁴. Les quatre bins combinés à moins de 1 Mpc donnent e_N < 6·10⁻⁵ a₀ et excluent
1,5·10⁻⁴ à 3,9σ formels, environ 3σ si les systématiques sont corrélées entre rayons. Comme la
part du champ 2M++ venant de moins de 10 Mpc n'est que 3,7·10⁻⁵, l'isolement des lentilles ne
peut pas abaisser e_N sous 10⁻⁴ (vérifié dans 2M++ le 2026-10-13, `scripts/efe_isolated_2mpp.py` :
les 112 galaxies de 15-40 Mpc qui satisfont le critère d'isolement de Brouwer et al. ont un champ
médian de 1,14·10⁻⁴ a₀, 16-84 % 0,9-1,4·10⁻⁴, aucune sous 5·10⁻⁵ ; l'isolement retire la part
locale, 1,2·10⁻⁵ contre 5·10⁻⁵, pas le champ de la structure à grande échelle). Les deux
lectures de la loi sont donc séparées par ces données :
la lecture « composante sombre », ou une loi sans effet de champ externe, tient ; QUMOND et RMOND
avec leur effet de champ externe sont en tension à 3σ. Le rayon viriel médian
des ajustements NFW est 165 kpc : les lentilles traversent le rayon viriel sans voir le déclin
en √(ln r / r) qu'un halo ΛCDM y impose.

**Au-delà des galaxies, les amas** (`scripts/cluster_rar.py`, `results/cluster_rar.md`) : sur
les 20 amas CLASH de Tian et al. 2020 (14-600 kpc, lentille forte et faible contre étoiles de la
galaxie centrale et gaz X) et les 11 amas X-COP d'Eckert et al. 2019, 2022 (masses
hydrostatiques et SZ, fractions de gaz à R2500, R500, R200), la loi avec le a₀ des galaxies
manque +0,57 dex en accélération à 100-300 kpc, +0,35 à R2500, +0,15 à R500 et +0,06 à R200,
soit en masse 4 à 7 fois les baryons observés entre 15 et 600 kpc, 1,8 fois à R500 et 1,3 fois
à R200. Le « facteur 2 » de la littérature est la valeur à R500 ; la masse manquante est
concentrée dans les cœurs, où g_bar vaut 0,1-0,8 a₀ et où la loi n'amplifie presque rien.
Les amas forment leur propre séquence de forme RAR mais d'échelle plus grande (CLASH 15 a₀,
X-COP 2,6 a₀ : pas de a₀ unique même entre amas). Une théorie de type MOND exige donc une
composante supplémentaire, 10¹³-10¹⁴ M☉ par amas, plus concentrée que le gaz chaud ; une lecture
« composante sombre » n'exige rien de plus.

**Sous a₀ à l'échelle stellaire, les binaires larges de Gaia** (`scripts/wide_binaries.py`,
`results/wide_binaries.md`, 2026-10-13) : 33 286 paires de séquence principale à moins de 250 pc
(El-Badry et al. 2021 ; parallaxes à 20σ, RUWE < 1,4, alignement fortuit < 1 %, les deux étoiles
sur la séquence principale), masses par la relation de Pecaut & Mamajek, statistique
ṽ = v_p/v_c calibrée sur les paires à g_N > 10 a₀ (médiane 0,563, celle d'une distribution
thermique d'excentricités en gravité newtonienne, 0,556). De 1000 à 0,1 a₀ le rapport à la
calibration reste entre 1,00 et 1,06. Une loi locale sans effet de champ externe prédit 1,36 à
0,3-1 a₀ et 1,66 à 0,1-0,3 a₀, contre 1,057 ± 0,022 et 1,055 ± 0,041 mesurés : exclue à 14-15σ,
bien au-delà des 10-20 % que compagnons cachés et excentricités peuvent déplacer. QUMOND avec le
champ galactique (e_N = 1,2 a₀, soit 1,8 a₀ après amplification) prédit 1,17 : 3 à 5σ au-dessus
des données, mais les systématiques connues sont de cette taille, et la littérature se partage
(Chae 2023 : amplification 1,3-1,5 en g ; Banik et al. 2024 : Newton à 16σ). Sous 0,1 a₀,
118 paires à 30 kAU donnent 1,30 ± 0,06, entre Newton et la loi, avec une fraction de ṽ > 1 qui
passe de 0,13 à 0,22 : le régime des paires contaminées ou non liées, qui ne tranche rien. Bilan
pour les deux lectures : « composante sombre » prédit Newton partout et n'est contredite par
rien ; une dynamique modifiée sans effet de champ externe est morte ici, avec effet de champ
externe elle est en tension à 3σ à 1 Mpc et à 3-5σ ici.

**L'inertie modifiée de Milgrom sur ces deux tests** (`scripts/mi_milgrom_predictions.py`,
`results/mi_milgrom_predictions.md`, 2026-10-13) : dans les modèles de Milgrom 2022 (PRD 106,
064060), l'équation du mouvement est â(ω) μ[𝒜(ω)/a₀] = â_N(ω) avec 𝒜(ω_in) = a_in + θ(ω_ex/ω_in) a_ex,
θ(1) = 1, θ décroissante, θ(0) ≈ 2-3 : le champ externe entre avec son accélération observée et
un poids fixé par le rapport des fréquences, et l'effet est isotrope. Binaires : a_ex = 1,8 a₀,
ω_ex/ω_in ≈ 0,01, poids θ(0). Le modèle prédit 1,09 et 1,10 (θ(0) = 2) ou 1,055 et 1,058
(θ(0) = 3) dans les bins 0,3-1 et 0,1-0,3 a₀, contre 1,057 ± 0,022 et 1,055 ± 0,041 mesurés ;
χ² sur les trois bins sous 1 a₀ : 13,4 au meilleur θ(0) = 2,3 (1,9-2,7), contre 32 pour Newton et
38 pour QUMOND, le reste venant du dernier bin que rien n'ajuste. Lentille : avec a_ex = 0,0105 a₀
(champ des galaxies isolées) et un poids θ_eff libre, les quatre bins donnent θ_eff < 0,2 à 2σ et
excluent θ_eff = 1, la force d'une théorie de champ, à 7σ, parce que la forme algébrique de
l'inertie modifiée fait sentir le champ externe dès que a_ex est une fraction de a, sans l'écran
de la masse fantôme de QUMOND à l'intérieur de r_t. Le modèle fournit ce faible poids si
ω_ex/ω_in > 3, c'est-à-dire si le champ cosmique varie en moins de 15 Gyr à 1 Mpc (période
orbitale 46 Gyr) et de 5 Gyr à 300 kpc (14 Gyr) : plausible pour la première, marginal pour la
seconde. Réserves : pas de version relativiste, donc pas de prédiction de lentille propre (on
suppose que la lumière suit les orbites test circulaires) ; construction de Fourier qui
présuppose des trajectoires stationnaires, ce que des orbites de 50 Gyr ne sont pas ; θ inconnue
hors θ(1) = 1 ; l'excès de 5 % des binaires est aussi à la portée des compagnons cachés. Les amas
ne l'atteignent pas : ni la lentille ni le gaz hydrostatique n'ont de version en inertie
modifiée. Bilan : l'inertie modifiée survit aux deux tests, au premier en le prédisant, au second
en n'y prédisant rien de ferme.

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
   à 300 kpc, où QUMOND prédit moins de 13 % de déclin ; tranché par la lentille à 0,3-1 Mpc des
   lentilles de faible masse, où QUMOND prédit 30 à 50 % de déclin et où la mesure est plate :
   e_N < 6·10⁻⁵ a₀, champ 2M++ exclu à ~3σ (`results/efe_lensing_mpc.md`) ; et les binaires
   larges de Gaia excluent à 15σ une loi locale sans effet de champ externe
   (`results/wide_binaries.md`). Une dérivation doit donc livrer un effet de champ externe qui
   n'agit pas à 1 Mpc sur une galaxie de 10¹⁰ M☉ dans un champ de 10⁻⁴ a₀ mais efface la loi
   dans une binaire plongée dans le champ galactique de 1,8 a₀ : aucune théorie de champ locale
   ne fait les deux ; il reste l'inertie modifiée non locale de Milgrom, qui fait les deux
   (`results/mi_milgrom_predictions.md` : +5 à +9 % prédit dans les binaires pour θ(0) = 2-3 contre
   +6 ± 2 % mesuré ; effet affaibli à 1 Mpc parce que le poids du champ externe décroît quand
   l'orbite interne est plus lente que la variation du champ), sans version relativiste ni
   prédiction pour les amas, et la lecture « composante sombre ». Et dire ce qui remplit les
   cœurs d'amas (`results/cluster_rar.md`).
4. La dépendance temporelle, H(z) ou Λ, testable à z = 1 à 2.

Aucune des trois dérivations ne livre les quatre. C'est l'état exact du problème ouvert.
