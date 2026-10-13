# Théorie — la formulation qui tient

Ce document remplace la théorie d'origine du dépôt (modèle d'accumulation d'énergie, puis
modèle baryonique à halo CGM, archivés dans `archive_v1/`). Il ne garde que ce qui a survécu
à l'audit du 10 octobre 2026 (voir `results/rms_tests*.md`, `results/optimize_constraints.md`,
`results/temporary_mass_budget.md`, `results/test_efe.md`, `results/model_v2.md`).
Implémentation : `scripts/model_v2.py`. Aucun paramètre libre par galaxie.

## 1. Ce que les données établissent

Sur 171 galaxies SPARC, avec les rapports masse-luminosité individuels de la table
enrichie, l'accélération observée ne dépend que de l'accélération baryonique locale :

- relation d'accélération radiale sur 2 784 points, dispersion 0,146 dex ;
- aucune corrélation du résidu avec la luminosité, le type, la brillance de surface,
  le rayon de disque, le taux de formation stellaire ou la fraction de gaz (|ρ| ≤ 0,2) ;
- la masse manquante vaut r·√(M_bary a₀/G) − M_bary à 0,005 dex près en médiane,
  pente 0,56 en log M_bary à rayon fixé (1,0 pour une masse proportionnelle aux baryons).

Toute énergie absorbée, stockée ou accumulée est proportionnelle à ses sources et
bornée par ce que la galaxie a reçu (≤ 4·10⁻⁴ M_bary c²) : elle ne peut produire ni
la racine carrée ni l'amplitude. Le mécanisme énergétique est abandonné.

## 2. Formules

### 2.1 Baryons

```
V_bar²(r) = V_gas|V_gas| + Υ_d V_disk|V_disk| + Υ_b V_bul|V_bul|
g_bar(r)  = V_bar²(r) / r
Υ_d = 0,5 (ou Y36 par galaxie, table enrichie),  Υ_b = 1,4 Υ_d
```

V_gas, V_disk, V_bul sont les vitesses des modèles de masse SPARC (Υ = 1, hélium
inclus dans V_gas). Υ multiplie une masse, donc V², jamais V.

### 2.2 Loi d'accélération

```
g_obs(r) = ν( g_bar(r) / a₀ ) · g_bar(r)          a₀ = 1,2·10⁻¹⁰ m/s²
V_model(r) = √( g_obs(r) · r )
```

Deux fonctions d'interpolation, équivalentes sur SPARC :

```
ν_RAR(y)    = 1 / (1 − e^{−√y})          (McGaugh, Lelli & Schombert 2016)
ν_simple(y) = 1/2 + √(1/4 + 1/y)         (Famaey & Binney 2005)
```

Limites : y ≫ 1, ν → 1, Newton ; y ≪ 1, g_obs → √(g_bar a₀).

### 2.3 Forme équivalente : amplification de la masse

```
M_eff(<r) = M_bar(<r) · ν( g_bar / a₀ )
```

C'est l'idée d'origine du dépôt, « amplification en fonction de a_bary », appliquée
à la masse incluse et non à une densité locale. Dans le régime profond,
M_eff = M_bar · (a₀/a_bar)^{1/2} : l'exposant est ½ sur a_bar, ou 1 sur
l'accélération totale avec la boucle de rétroaction (a = g_bar · a₀/a ⇒ a = √(g_bar a₀)).
Appliquée à une densité locale, la même formule fait croître la densité vers l'extérieur
et ne peut pas donner de courbe plate ; c'est l'erreur qui a forcé β = 1,4.

### 2.4 Composante sombre équivalente

```
M_x(<r)  = [ν(g_bar/a₀) − 1] · M_bar(<r)
ρ_x(r)   = (1/4πG) ∇·[ (ν − 1) g_N ]
```

Au dernier point mesuré : M_x/M_bar médian 4,1 ; densité surfacique
Σ_x = M_x/(πr²) médiane 54 M☉/pc², 16ᵉ à 84ᵉ centile 39 à 80. C'est le halo trouvé
par les recherches sur grille du modèle 1.0 (cœur ≈ 1,7 R_d, profil externe entre r⁻²
et r⁻⁴), ici sans aucun paramètre.

## 3. Prédictions sans paramètre, vérifiées

| Prédiction | Résultat sur SPARC |
|---|---|
| Courbes de rotation, Υ fixés | RMS médian 13,7 km/s, χ²/N 11 ; précoces 25, intermédiaires 11, tardives 11 |
| Idem, Υ_d libre avec prior SPS de 0,1 dex | 9,2 km/s, χ²/N 3,7, biais −1,7, résidu externe −4 |
| Tully-Fisher baryonique, V_flat⁴ = G M_bary a₀ | 135 galaxies, décalage +0,026 dex, dispersion 0,088 dex |
| Masse manquante en √M_bary · r | pente 0,56, décalage −0,005 dex |

Le modèle 1.0 atteignait 19 à 21 km/s avec cinq paramètres et une normalisation par
galaxie ; le modèle énergétique 35 km/s ; les baryons seuls 43.

## 4. Deux lectures physiques, mêmes courbes de rotation

**A. Inertie modifiée.** L'inertie d'un corps dépend de son accélération :
m_i = m · μ(a/a₀), avec μ(x) = x/(1+x) pour ν_simple. C'est la seule façon de rester
fidèle au point de départ E = γmc² : la relation entre énergie et inertie n'est pas
celle de la relativité restreinte quand a ≪ a₀. La correction est une fonction de
l'accélération, pas de la vitesse ; γ(v) vaut 1 + 5·10⁻⁷ dans une galaxie.

**B. Composante sombre.** Une masse non baryonique de distribution M_x(<r) ci-dessus.
Elle ne peut pas être du gaz (densité 100 fois supérieure aux mesures X et de dispersion),
ni de l'énergie (pression 10⁷ fois celle du milieu interstellaire) : elle doit être
froide et sans pression, ce qui est la définition de la matière noire froide.

Les deux lectures ne diffèrent que hors des courbes de rotation : effet de champ
externe, amas, CMB, binaires larges.

## 5. L'échelle a₀ et le fond cosmologique

```
a₀ = 1,2·10⁻¹⁰ m/s²        c H₀ / 2π = 1,08·10⁻¹⁰ (H₀ = 70)       rapport 1,11
```

C'est l'élément du programme initial qui survit : le fond cosmologique fixe la
dynamique locale. Deux hypothèses testables, incompatibles entre elles :

- **a₀ ∝ H(z)** : a₀(z=1) = 1,76 a₀, a₀(z=2) = 3,0 a₀. Prédit plus de masse manquante
  à haut redshift à baryons fixés, ×1,3 en accélération à z = 1, ×1,7 à z = 2.
- **a₀ ∝ c√(Λ/3)** : a₀ constant dans le temps. Prédit la même relation d'accélération
  à toute époque.

Tranché le 2026-10-12 (`docs/A0_ORIGIN.md`, `scripts/a0_highz_shape2.py`) : la forme des
courbes de rotation en gaz froid à z = 1,5-4,8, à masse baryonique libre, donne
n = 0,05 (−0,25 à 0,30) dans a₀ ∝ H(z)ⁿ et place a₀ ∝ H(z) à 4,6-5,2σ au-dessus du meilleur
ajustement. L'amplitude de la relation à z ≥ 2 dit la même chose (r = 1,05, 0,87-1,23) mais
dépend des masses baryoniques à 0,1 dex près. a₀ est constant : c'est la version
« échelle fixée par Λ » qui survit, avec un préfacteur ≈ 0,2 à expliquer.

## 6. Ce qui tue la théorie

| Observable | Lecture A, inertie modifiée | Lecture B, composante sombre |
|---|---|---|
| Effet de champ externe dans SPARC | prédit, mais sous le seuil de mesure à e_N ≈ 10⁻⁴ a₀ d'après l'équation de champ ; non détecté dans notre reproduction (plage dynamique 0,24 dex seulement) | absent |
| Vitesses de lentille à 0,3-1 Mpc, lentilles de faible masse | déclin de 30 à 50 % prédit par QUMOND avec le champ 2M++ ; mesurées plates : e_N < 6·10⁻⁵ a₀, champ 2M++ exclu à 3,9σ formels, ~3σ avec les systématiques corrélées (`results/efe_lensing_mpc.md`) | plates, comme observé |
| Binaires larges Gaia sous a₀ | loi locale sans effet de champ externe : ×1,36 à ×1,66 prédit entre 1 et 0,1 a₀, mesuré 1,05 ± 0,02, exclue à 15σ ; QUMOND avec le champ galactique : ×1,17, mesuré 1,06 ± 0,02, 3 à 5σ sous la prédiction mais systématiques de même taille (`results/wide_binaries.md`) | newtoniennes, comme observé jusqu'à 0,1 a₀ |
| Amas de galaxies | manque, avec le a₀ des galaxies, un facteur 4-7 en masse entre 15 et 600 kpc, 2 à R500, 1,3 à R200 (`results/cluster_rar.md`) | ajustés |
| Troisième pic du CMB, BAO | exige une extension relativiste (Skordis & Złośnik 2021) | prédits |
| Amas du Boulet | difficile | naturel |
| Courbes de rotation à z > 1 | selon a₀(z) | évolution par concentration des baryons |

La théorie v2 est morte si l'une de ces lignes est tranchée contre les deux lectures
à la fois : par exemple des binaires larges newtoniennes **et** une relation d'accélération
dont la dispersion intrinsèque reste nulle à toutes les masses, ce qu'aucun halo ne produit
sans réglage.

État au 13 octobre 2026 : deux lignes sont tranchées contre la lecture A seule. La lentille à
0,3-1 Mpc (`scripts/efe_lensing_mpc.py`) : le bin de plus faible masse de Mistele et al. 2024
(log M_b = 10,1, r_t = 300 kpc pour le champ 2M++) reste plat jusqu'à 1 Mpc là où QUMOND
prédit 30 à 50 % de déclin, et les quatre bins combinés donnent e_N < 6·10⁻⁵ a₀ ; comme la part
du champ 2M++ venant de moins de 10 Mpc n'est que 3,7·10⁻⁵, l'isolement des lentilles ne peut
pas expliquer un champ si faible (vérifié dans 2M++ : les galaxies qui satisfont le critère
d'isolement de Brouwer et al. gardent un champ médian de 1,1·10⁻⁴ a₀, aucune sous 5·10⁻⁵,
`results/efe_isolated_2mpp.md`). Les amas (`scripts/cluster_rar.py`) : avec le a₀ des galaxies,
la loi manque un facteur 4-7 en masse entre 15 et 600 kpc (CLASH), 2 à R500 et 1,3 à R200
(X-COP), la masse manquante étant concentrée dans les cœurs. La lecture B absorbe les deux. La
lecture A ne survivait que sous une forme sans effet de champ externe ; les binaires larges de
Gaia (`scripts/wide_binaries.py`) ferment cette porte : une loi locale sans effet de champ
externe prédit ×1,4 à ×1,7 sur la vitesse relative entre 1 et 0,1 a₀, et la mesure donne
1,05 ± 0,02. Aucune théorie de champ locale ne peut à la fois ignorer un champ de 10⁻⁴ a₀ à
1 Mpc et effacer la loi dans le champ galactique de 1,8 a₀. Reste à la lecture A l'inertie
modifiée non locale de Milgrom, dont les prédictions pour ces deux tests ne sont pas calculées.
La lecture B tient partout, au prix de devoir produire la loi.

## 7. Ce que devient E = γmc²

Une identité de conversion. Appliquée à une galaxie, elle rend M_bary à un millionième
près : énergie cinétique 2·10⁻⁷, liaison −4·10⁻⁷, rayonnement présent 10⁻⁸. Elle ne
contraint aucune des formules ci-dessus. Le contenu physique du programme est dans a₀
et dans ν, pas dans c².
