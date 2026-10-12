# L'origine de a₀ : deux hypothèses, leurs prédictions, et ce que les données disent

Calculs : `scripts/a0_evolution.py`, résultats `results/a0_evolution.md`.

## 1. Le fait

La loi d'accélération de SPARC a une seule échelle, a₀ = 1,2·10⁻¹⁰ m/s². Elle est
numériquement proche de deux grandeurs cosmologiques :

| Expression | Valeur (H₀ = 70) | Rapport à a₀ |
|---|---|---|
| c H₀ / 2π | 1,08·10⁻¹⁰ | 0,90 |
| c √(Λ/3) / 2π = c H₀ √Ω_Λ / 2π | 9,1·10⁻¹¹ | 0,75 |

Avec H₀ = 73, la première devient 1,13·10⁻¹⁰. Le facteur 2π est conventionnel ; ce qui
est significatif est l'ordre de grandeur, a₀ ≈ cH₀/6, alors que rien dans la physique
galactique ne contient c ni H₀. L'argument de Milgrom 1999 l'exprime en températures :
un corps accéléré à a₀ a une température d'Unruh de 4,9·10⁻³¹ K, l'horizon de de Sitter
une température de 2,3·10⁻³⁰ K ; le régime MOND commence quand la première passe sous
la seconde.

## 2. Les deux hypothèses

**H_A, a₀ suit l'expansion** : a₀(z) = c H(z)/2π. C'est la version « le fond cosmologique
agit depuis la naissance » : l'échelle était plus grande quand l'Univers s'expandait plus vite.

**H_B, a₀ suit le vide** : a₀ = c √(Λ/3)/2π, constante dans le temps. Milgrom 1999,
Verlinde 2016, Smolin 2017 sont de cette famille : c'est la constante cosmologique, pas
le taux d'expansion, qui fixe l'échelle.

Elles coïncident aujourd'hui à 20 % près et divergent dans le passé.

## 3. Prédictions, à distribution baryonique fixée

| z | H(z)/H₀ = a₀_A(z)/a₀ | Tully-Fisher baryonique, Δlog M_b à V fixe | Relation d'accélération, Δlog g_obs à g_bar fixe | H_B |
|---|---|---|---|---|
| 0,5 | 1,31 | −0,12 | +0,06 | 0 |
| 1,0 | 1,76 | −0,25 | +0,12 | 0 |
| 2,0 | 2,97 | −0,47 | +0,24 | 0 |
| 2,3 | 3,39 | −0,53 | +0,27 | 0 |

Sous H_A, une galaxie de z = 2,3 a, à vitesse plate égale, deux fois moins de baryons
qu'aujourd'hui, et deux fois plus d'accélération manquante dans le régime profond. Sous
H_B, aucune évolution. Les deux prédictions supposent le régime profond, g_bar ≪ a₀, donc
les parties externes des disques, au-delà de 2 à 3 rayons effectifs.

## 4. Ce que les données disent

**Tully-Fisher baryonique à haut redshift**, Übler et al. 2017, KMOS³ᴰ, 135 disques massifs,
vitesse circulaire maximale corrigée de la pression, référence locale Lelli et al. 2016 :

| | Observé, Δlog M_b à V fixe | H_A | H_B |
|---|---|---|---|
| z ≈ 0,9 | −0,44 | −0,25 | 0 |
| z ≈ 2,3 | −0,27 | −0,53 | 0 |

L'évolution observée n'est pas monotone : plus forte à z ≈ 0,9 qu'à z ≈ 2,3. Ces décalages
sont reproduits à 0,02 dex près par `scripts/a0_highz.py` (−0,46 et −0,30) et s'expliquent
sans évolution de a₀ : ces disques compacts ont g_bar ≈ 1 à 4 a₀ au rayon de mesure, leur
vitesse circulaire maximale à 2,2 R_d dépasse la vitesse plate asymptotique, et la relation
de Tully-Fisher les place mécaniquement sous la relation locale. Le test point par point
ci-dessous remplace cette lecture globale.

**Relation d'accélération point par point, z = 0,6 à 4,6** (`scripts/a0_highz.py`,
`results/a0_highz.md`, figure `results/a0_highz.png`) : pour chaque galaxie ayant une
accélération totale g_obs = V_c²/R et une accélération baryonique g_bar au même rayon, la
loi g_obs = ν(g_bar/a₀) g_bar est ajustée pour r = a₀(z)/a₀(0). Échantillons : RC100
(Nestor Shachar et al. 2023, 100 disques massifs à z = 0,6-2,5, f_DM(R_e), V_c(R_e) et R_e
transcrits de leur table 3, `data/external/highz/`), MSA-3D (JWST/NIRSpec 2026, 30 galaxies
à z = 0,5-1,7, log M★ 9-11), KMOS³ᴰ (Übler 2017, 135, tailles de van der Wel 2014), KROSS
(Harrison 2017, 161 galaxies à z ≈ 0,9, gaz moléculaire de Tacconi 2018), quatre disques
ALMA à z ≈ 4,4 avec courbes publiées (Roman-Oliveira et al. 2024) et sept autres à rayon
externe supposé (Rizzo 2021, Fraternali 2021).

| échantillon | z médian | N | g_bar/a₀ médian | r ajusté (68 %) | H_A exige | masses à retirer pour que H_A soit le meilleur ajustement |
|---|---|---|---|---|---|---|
| RC100 | 0,8 | 26 | 1,1 | 1,35 (0,9-2,0) | 1,6 | −0,04 dex |
| RC100 | 1,45 | 33 | 1,7 | 1,15 (0,7-1,7) | 2,3 | −0,12 dex |
| RC100 | 2,2 | 41 | 1,7 | 1,05 (0,7-1,5) | 3,3 | −0,22 dex |
| MSA-3D | 0,9 et 1,3 | 30 | 0,8 | 2,6 (2,0-3,3) | 1,7 et 2,1 | +0,1 dex |
| KMOS³ᴰ | 0,9 | 65 | 0,9 | 2,75 (2,3-3,2) | 1,7 | +0,1 dex |
| KROSS | 0,9 | 161 | 0,6 | 2,1 (1,9-2,4) | 1,6 | +0,1 dex |
| z ≈ 4,4, quatre courbes | 4,35 | 4 | 11 | 0,5 (0,1-2,2) | 6,8 | −0,36 dex |

Trois lectures :

1. À z ≥ 1,7, les données telles quelles donnent r = 1,00 (0,66-1,45), soit H_B, là où H_A
   exige 3,3 à 6,8 ; H_A y est rejetée à Δχ² = 17 avec les erreurs publiées. Mais une baisse
   uniforme de 0,2 dex de toutes les masses baryoniques ramène r à 2,95 et fait de H_A le
   meilleur ajustement. H_A n'est donc pas exclue : elle exige que les masses baryoniques
   soient surestimées de 0,22 dex à z ≈ 2,2 et de 0,36 dex à z ≈ 4,5, partout du même signe.
   Les masses à ces redshifts, étoiles par SED et gaz par relations d'échelle ou par α_[CII],
   ne sont pas connues à ce niveau.
2. À z ≈ 0,9-1,3, les quatre échantillons donnent r ≈ 2 à 2,7, au-dessus même de H_A ; +0,1 dex
   de baryons les ramène sur H_A, +0,2 dex sur H_B. Ils sont en régime g_bar ≈ 0,6-0,9 a₀,
   où a₀ ∝ 1/M_bar : le gaz non mesuré (H₂, HI) et le rapport M/L décident, pas a₀.
   Non discriminant.
3. MUSE-DARK III (Ciocan et al. 2026) revendique a₀ croissant de 2,0 à 2,7·10⁻¹⁰ entre
   z = 0,45 et 1,3, soit a₀(z) = (1,0 + 1,59 z)·10⁻¹⁰, sur 79 galaxies de log M★ > 8,8. Leur
   masse baryonique est étoiles + un disque HI symbolique, sans gaz moléculaire ; ils écrivent
   eux-mêmes qu'il faudrait +0,2 à +0,45 dex de masse pour retrouver 1,2·10⁻¹⁰. Le H₂ de
   Tacconi 2018 pour log M★ = 9,3-10,1 vaut justement +0,1 à +0,47 dex à ces redshifts :
   corrigé, leur a₀ passe de 1,3-1,6 à z = 0,45 à 0,9-1,3 à z = 1,3, aucune évolution,
   compatible avec 1,2 ± 0,26. L'évolution annoncée est le gaz qu'ils omettent.

**Second passage, gaz mesuré et forme des courbes, 2026-10-12** (`scripts/a0_highz_coldgas.py`,
`results/a0_highz_coldgas.md`, figure `results/a0_highz_coldgas.png`) : pour les échantillons à
f_DM, l'observable devient ν_obs = 1/(1 − f_DM) contre ν(g_bar/(r a₀)), avec g_bar construit sur
les masses et tailles des auteurs ; les vitesses n'entrent que par leurs ajustements. Échantillons
ajoutés : NOEMA3D (Jolly et al. 2026, 10 disques à z = 1,1-1,6, gaz CO mesuré), CRISTAL (Lee et
al. 2025, 14 disques [CII] à z = 4,4-5,7), ALPAKA + JWST (Roman-Oliveira et al. 2026, trois
disques suivis jusqu'à 4 R_e, normalisation du gaz ajustée), Rizzo et al. 2021 avec leurs masses
baryoniques. Et un test nouveau : la forme des quatre courbes publiées à z ≈ 4,4, masse
baryonique marginalisée avec un prior de 0,3 dex, qui ne teste que la montée de ν avec le rayon.

| échantillon | z | N | r (68 %) | H_A exige | H_B : Δχ² | H_A : Δχ² |
|---|---|---|---|---|---|---|
| RC100 | 2,2 | 41 | 0,89 (0,71-1,07) | 3,3 | +0,3 | +52 |
| NOEMA3D, gaz CO | 1,3 | 10 | 2,0 (1,2-3,0) | 2,0 | +1,8 | −0,6 |
| CRISTAL, [CII] | 5,2 | 14 | 2,6 (1,7-3,7) | 8,5 | +4,2 | +13 |
| ALPAKA ID13 à 21 kpc, g_bar = 0,16 a₀ | 2,1 | 1 | ν_obs 5,6 ± 4,0 | ν = 4,9 | +0,8 | 0 |
| quatre courbes z ≈ 4,4, amplitude au dernier anneau | 4,35 | 4 | 0,5 (0,1-2,2) | 6,8 | +0,2 | +4,8 |
| mêmes courbes, forme seule, masse libre | 4,4 | 17 anneaux | ≤ 0,35 | 6,9 | +3,9 | +17 |

Lecture : à z ≥ 1,7, 60 galaxies donnent r = 1,05 (0,87-1,23) et rejettent H_A à Δχ² = 64 au
pied de la lettre ; mais −0,1 dex sur toutes les masses baryoniques rend H_A et H_B
équivalentes, et −0,2 dex renverse le verdict. Le test d'amplitude repose donc sur l'échelle
des masses baryoniques à 0,1 dex près, précision que ni les masses stellaires par SED ni α_CO
ou α_[CII] n'atteignent. NOEMA3D, seul échantillon à gaz CO mesuré à z ≈ 1,3, est en régime
newtonien (g_bar ≈ 3,5 a₀) et son r bascule de 0,1 à 3 pour ±0,2 dex de masse. Le test de
forme, lui, ne dépend pas de cette échelle : sur les quatre courbes à z ≈ 4,4 il exclut H_A à
4σ et reste compatible avec H_B à 2σ, les courbes déclinant même un peu plus vite que la loi
locale ne le prévoit. Il ne porte que sur 17 anneaux.

**Troisième passage, pentes externes de KURVS, z ≈ 1,5** (`scripts/a0_highz_shape.py`,
`results/a0_highz_shape.md`) : Puglisi et al. 2023 tabulent, pour 22 galaxies suivies en Hα
jusqu'à 4 R_e, la vitesse de rotation à 3 et 6 longueurs d'échelle, donc la pente externe
t = v(6 R_d)/v(3 R_d), et f_DM(R_e) pour dix disques dominés par la rotation. En amplitude,
f_DM donne r = 2,75 (2,0-3,6), qui bascule de 0,8 à 5,9 pour ±0,2 dex de masse : même
dégénérescence qu'à z ≈ 1. En forme, la loi borne la pente d'un disque exponentiel seul à
t ≤ 1,05 quel que soit a₀ ; H_B prédit 0,90-0,95 et H_A 0,93-0,98 pour ces galaxies, un
écart de 0,03 inférieur aux erreurs. Les pentes publiées valent 0,92-1,33 (médiane 1,03), et
1,06-1,40 une fois la correction de pression de Burkert et al. 2010 appliquée (6 à 12 σ₀²
ajoutés à V_c², avec σ₀ = 40-75 km/s contre v = 70-280) : ces valeurs mesurent la correction
de pression et l'extrapolation de la courbe ajustée au-delà de R_Hα,max, pas a₀. Avec la
moitié des baryons dans un disque de gaz deux fois plus étendu que les étoiles et sans
correction de pression, la loi donne t = 1,02 (H_B) et 1,04 (H_A) contre 1,03 observé :
compatible avec les deux, r = 1,7 (0,9-3,4). Leçon : les pentes en gaz ionisé à v/σ₀ < 5 ne
tranchent pas ; il faut le gaz froid.

**Quatrième passage, forme de toutes les courbes en gaz froid récupérables, 2026-10-12**
(`scripts/a0_highz_shape2.py`, `results/a0_highz_shape2.md`, figure `results/a0_highz_shape2.png`) :
les courbes n'existant qu'en figures ou en paramètres, elles ont été reconstruites : Rizzo et
al. 2020 et 2021, six disques [CII] lentillés à z = 4,2-4,8, courbes analytiques depuis leurs
paramètres publiés (V_t, R_t, β, ξ, σ₀, R_σ) avec leur correction de dérive asymétrique et
leur structure baryonique (Sérsic stellaire, disque de gaz exponentiel) ; Roman-Oliveira et
al. 2024, quatre disques [CII] à z = 4,3-4,5 (Zenodo) ; Lelli et al. 2023, deux disques CO à
z = 1,47 et 2,24, anneaux lus sur leur figure 7, structure bulbe + disque de leurs ajustements
HST ; ALPAKA, 19 disques CO/[CI] à z = 0,6-3,6, anneaux lus sur une figure de 937 pixels de
large, taille et masse libres. Pour chaque galaxie, la courbe est confrontée à la loi avec la
masse baryonique totale libre (prior 0,3 dex) : seule la montée de ν avec le rayon compte.

| jeu | z | N | H_B : Δχ² au-dessus du meilleur r | H_A : Δχ² | A − B |
|---|---|---|---|---|---|
| Rizzo 2020/2021, [CII] lentillé | 4,2-4,8 | 6 | +3,5 | +8,5 | +5,0 |
| Roman-Oliveira 2024, [CII] | 4,3-4,5 | 4 | +3,9 | +17,4 | +13,5 |
| Lelli 2023, CO | 1,5 et 2,2 | 2 | 0,0 | +2,0 | +2,0 |
| jeu propre, total | 1,5-4,8 | 12 | +7,4 | +27,9 | +20,5 |
| ALPAKA, CO/[CI], tailles libres | 0,6-3,6 | 19 | +18,4 | +24,8 | +6,4 |

Ajustement global a₀ ∝ H(z)ⁿ sur la forme des courbes : n = −0,65 (−1,0 à 0,0) sur le jeu
propre, n = 0,05 (−0,25 à 0,30) avec ALPAKA ; H_A (n = 1) est à Δχ² = +21 (jeu propre) et +27
(tout), soit 4,6 et 5,2σ ; H_B (n = 0) est à +0,9 et 0,0. Les erreurs des courbes de Rizzo
sont surestimées par construction (χ²_min/N = 0,02), ce qui rend ces Δχ² conservateurs. Dix
des douze galaxies propres sont à z > 4, où H_A exige r = 6,6 à 7,6 et où un disque qui tourne
à 200-500 km/s dans 2-4 kpc devrait voir ν doubler entre R_e et 3 R_e ; aucune ne le montre.

La courbe r(z) prise au pied de la lettre, 2,5 à z ≈ 1, 1,0 à z ≈ 2,2, 0,5-0,9 à z ≈ 4,5,
n'est monotone pour aucune hypothèse : c'est la signature de systématiques de masse qui
diffèrent d'un échantillon à l'autre, pas d'un a₀ cosmologique. Les contrôles passent :
g_bar des auteurs de RC100 contre un modèle disque + bulbe construit sur leurs masses,
−0,07 dex ; Milgrom 2017 sur les six galaxies de Genzel 2017 (4 a₀ à z ≈ 2 « all but
excluded ») et Lelli et al. 2023 (deux disques CO à z = 1,5 et 2,2 ajustés avec le a₀ local)
vont dans le même sens que RC100.

**Fractions sombres dans le rayon effectif**, Nestor Shachar et al. 2023, RC100, 100 disques
à z = 0,6 à 2,5 : f_DM(R_e) médiane 0,38 ± 0,23 à z ≈ 1, 0,27 ± 0,18 à z ≈ 2, anticorrélée
à la densité de surface baryonique ; la moitié des galaxies à z ≈ 2 sont des disques
maximaux. Ces galaxies ont des densités de surface de l'ordre de 10⁹ M☉/kpc², soit
g_bar ≈ 5 à 10 a₀ dans R_e : elles sont en régime newtonien à l'intérieur, et une loi qui
ne dépend que de g_bar y prédit une fraction sombre faible et décroissante avec Σ, ce qui
est observé. Ce test ne discrimine ni H_A ni H_B de ΛCDM, qui l'explique par la
concentration des baryons ; il confirme seulement que f_DM est une fonction de g_bar à
toutes les époques.

## 5. Verdict (2026-10-12, après le quatrième passage)

Deux tests indépendants, deux réponses cohérentes. L'amplitude : à z ≥ 2, 60 galaxies donnent
r = 1,05 (0,87-1,23) là où H_A exige 3,3 à 8,5, mais ce test bascule pour 0,1 à 0,2 dex de
surestimation des masses baryoniques, précision hors de portée des masses SED et des facteurs
α_CO ou α_[CII]. La forme : sur 12 courbes en gaz froid à z = 1,5-4,8, masse baryonique
libre, H_A est à Δχ² = +21 au-dessus du meilleur ajustement et H_B à +0,9 ; avec les 19
courbes ALPAKA, n = 0,05 (−0,25 à 0,30) dans a₀ ∝ H(z)ⁿ, H_A à 5,2σ. Le test de forme ne
dépend ni de α_CO, ni de M/L, ni de la fraction de gaz ; il dépend de la structure baryonique
adoptée (Sérsic, disque exponentiel), de la correction de dérive asymétrique de Rizzo et al.
et de la lecture des figures pour ALPAKA et Lelli. Verdict : a₀ ne suit pas l'expansion. Dans
les termes du chantier, c'est H_B, la version « échelle fixée par Λ », qui survit ; l'intuition
« le fond cosmologique agit depuis la naissance » exigerait des courbes à z ≈ 4,5 deux fois
plus montantes que celles que l'on observe. À z ≤ 1,5 l'évolution publiée par MUSE-DARK III
reste l'effet du gaz moléculaire absent de son bilan.

## 6. L'observation qui tranche

Le test de forme est fait, à 4,6-5,2σ, sur 12 à 31 courbes reconstruites. Ce qui le
consoliderait : les mêmes courbes sous forme tabulée par leurs auteurs (ALPAKA, CRISTAL,
NOEMA3D, RC100) au lieu de lectures de figures et de paramètres ; des courbes en gaz froid
(CO, [CII], v/σ₀ > 5) résolues sur 4 à 6 anneaux au-delà de R_e à z ≥ 2, où H_A prédit une
montée de ν deux à trois fois plus forte que la loi locale ; et une structure baryonique
mesurée (JWST) plutôt que supposée. KURVS montre qu'en gaz ionisé à v/σ₀ ≈ 2 la correction de
pression domine la pente et interdit la mesure ; le gaz froid est la seule voie. Les données existent en partie
: les 19 disques ALPAKA dont les profils cinématiques sont annoncés publics, les 14 courbes
CRISTAL, les courbes externes de RC100 et de Genzel et al. 2020, les deux disques CO de Lelli
et al. 2023, aucune n'étant encore tabulée. Une quarantaine de courbes de ce type, à 17
anneaux pour quatre galaxies on est déjà à 4σ, trancherait sans dépendre des masses. Le test
d'amplitude, lui, ne tranchera qu'avec des masses baryoniques à 0,1 dex, ce qui demande α_CO
calibré galaxie par galaxie et des rapports M/L dynamiques.
