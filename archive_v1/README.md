# Archive du modèle 1.0

Ce répertoire contient, tels qu'ils étaient, les documents, scripts et résultats du
modèle d'« accumulation d'énergie » développé jusqu'au 10 octobre 2026 et abandonné après
l'audit décrit dans `results/rms_tests.md`, `results/rms_tests_v1.md`,
`results/optimize_constraints.md` et `results/temporary_mass_budget.md` à la racine du dépôt.

Rien ici n'est maintenu ni cité par le travail actuel. L'état de référence de ce modèle est
la release 1.0 (tag `1.0`), laissée intacte.

- `THEORY_v1.md` : la théorie d'origine (densité équivalente ρ_eq ∝ (a₀/a_bary)^β intégrée
  radialement, puis modèle baryonique à halo CGM).
- `root/` : les notes de développement qui étaient à la racine (MEMORY, COMPARISON,
  DEPROJECTION, EXTENSIONS, INCLINATION, SPARC_QUALITY).
- `docs/` : les quelque 230 notes courtes de la session de conception (hypothèses, plans,
  relectures, variantes de paramètres). La plupart sont des coquilles de quelques lignes avec
  des nombres non reproduits par le code.
- `scripts/` : les scripts de cette période. L'audit a montré que la majorité sont des
  coquilles aux valeurs codées en dur ; `radial_model.py` est le seul qui calcule
  réellement, et son résultat (RMS 38 km/s) ne correspond pas aux valeurs annoncées.
- `results/` : les tables de résultats commitées avant que les données ne le soient.

Pourquoi le modèle a été abandonné, en deux chiffres : la masse temporaire que l'énergie en
transit peut représenter vaut 10⁻⁷ de la masse baryonique quand il en faudrait 3 fois plus
que les baryons, et sa formule, une fois E₀ supprimé et l'amplification portée sur la masse
incluse, se réduit à la relation d'accélération de Milgrom, qui est ce que le dépôt teste
désormais (`docs/THEORY.md`).
