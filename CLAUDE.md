# CLAUDE.md — Projet Hema (Stage Pierre Fabre)

## Contexte du projet

Stage de fin d'études M2 Science des Données de Santé — **Omar Hassak**, Pierre Fabre (2026).
**Intitulé :** Prédiction de l'hématotoxicité d'un médicament.

**Objectif :** Reproduire le modèle PK/PD semi-mécaniste de Fornari et al. 2019, validé sur le rat et transposé à l'humain pour le carboplatine, puis adapter ce pipeline à de nouvelles molécules développées en interne chez Pierre Fabre (avec estimation des paramètres en R).

**Livrables :**
- Mémoire de recherche format IMRaD (Master M2)
- Scripts/modèle documentés (Pierre Fabre)
- Slides hebdomadaires d'avancement (réunions de suivi)

---

## Exécution

Tous les scripts R doivent être lancés **depuis le répertoire `scripts/`** (chemins relatifs pour `source()` et `read.csv()`).

```bash
# Simulation rat — Figure 3 + Figure S1 + VPC + grades
cd scripts && Rscript run_CORRECT.R

# Simulation humaine — Figure 4 + VPC + Figure 4c (grades)
cd scripts && Rscript run_human.R

# Calibration neutrophiles (humain) — à lancer EN PREMIER
cd scripts && Rscript calibrate_neut.R

# Calibration plaquettes (humain) — à lancer APRÈS calibrate_neut.R
cd scripts && Rscript calibrate_plt.R

# Générer les présentations PowerPoint
python3 make_pptx.py
python3 reunion_point/make_pptx.py
```

**Packages R requis :** `deSolve`, `ggplot2`, `gridExtra`, `grid`, `scales`
**Package Python requis :** `python-pptx`
**Script d'installation :** `scripts/install_deps.R`

---

## Architecture

### Chaîne de `source()` (ordre critique)

```
run_CORRECT.R (rat)                run_human.R (humain)
  ├── parameters_rat.R               ├── parameters_human.R
  ├── parameters_FORNARI_CORRECT.R   ├── parameters_FORNARI_CORRECT.R
  ├── pkpd_model_FORNARI.R           ├── pkpd_model_FORNARI.R
  ├── data_fornari.R                 ├── plots_human.R
  ├── plots.R                        ├── plots_grades.R
  ├── plots_grades.R                 └── vpc.R
  └── vpc.R
```

**CRITIQUE :** `parameters_FORNARI_CORRECT.R` doit être sourcé APRÈS `parameters_rat.R` ou `parameters_human.R`. Il lit `init_pars` et `init_state` depuis l'environnement global et dérive tous les taux de transit (`k_tr_*`), taux de prolifération (`k_prol_*`), `k_stem`, et les conditions initiales des compartiments transit (Eq. S4).

### Fichiers clés

| Fichier | Rôle |
|---------|------|
| `pkpd_model_FORNARI.R` | Système ODE (`pkpd_fornari`) + wrapper `simulate_all()` |
| `parameters_rat.R` | Paramètres rat — baselines, PK, effets du médicament (Tables 1-2) |
| `parameters_human.R` | Paramètres humain — allométrie, PK Zandvliet 2008, Slopes IC50-scalés |
| `parameters_FORNARI_CORRECT.R` | Paramètres dérivés via Eq. S4 (lambda, k_tr, k_prol, k_stem, états initiaux transit) |
| `data_fornari.R` | Chargement CSVs digitalisés depuis `../data/` → liste `obs_fornari` |
| `vpc.R` | VPC Monte Carlo (1000 simulations, erreur résiduelle log-additive, sigma Table S4) |
| `plots_grades.R` | Grading NCI-CTCAE v5.0 (neutropénie + thrombocytopénie) |
| `calibrate_neut.R` | Calibration 1D de `Slope_CMP` — nadir Neut cycles 1/2 vs Figure 4 — **modifie `parameters_human.R`** |
| `calibrate_plt.R` | Calibration 2D de `delta_Plt` et `gamma_prolTrans` — **modifie `parameters_human.R`** |
| `compare_options.R` | Comparaison de 3 stratégies IIV pour la Figure 4c |

### Répertoires de sortie

- `scripts/results_CORRECT/` — PDFs rat (Figure 3, Figure S1, VPC, grades)
- `scripts/results_HUMAN/` — PDFs humain (Figure 4, VPC, Figure 4c)
- `reunion_point/` — Présentations PowerPoint

---

## Structure du modèle ODE

**23 variables d'état :**
`C1, C2, Damage, MPP, CMP, MEP, T1_Neut–T3_Neut, Neut, T1_Mono–T3_Mono, Mono, T1_Ret–T3_Ret, Ret, RBC, T1_Plt–T3_Plt, Plt`

**Mécanismes clés :**

- **PK (2 compartiments) :** `C1` (plasma) et `C2` (tissus), infusion via `pars$rate_fun(time)` (closure issue de `make_repeated_infusion()`).
- **Damage :** `dDamage = k_dam × C_free_µM − k_rep × Damage`. Fraction libre : `fu_t = fu_inf + (fu0 − fu_inf) × exp(−k_bind × t)`. Pour le carboplatine, `fu = 1` (liaison protéique négligeable).
- **Effet du médicament sur les progéniteurs :** Facteur multiplicatif `(1 − Slope_X × Damage)` sur la prolifération (MPP, CMP, MEP).
- **Effet sur les transits Ret/Plt :** `delta_X × Slope_MEP × k_prol_X × Damage` appliqué aux compartiments transit prolifératifs T1 et T2.
- **Feedbacks (Eq. 11-13) :** Ratios normalisés (clampés 0.2–50 pour stem/mat, 0.2–10 pour prol) élevés aux puissances gamma.
- **Amplification lambda :** Les transits Ret et Plt ont 2 compartiments prolifératifs (lambda1 = lambda2 = 2). Les conditions initiales de T1 diffèrent de T2/T3 d'un facteur `1/lambda` (Eq. S3).

---

## Conventions

- **Langue :** Commentaires en **français**, noms de variables/fonctions/labels en **anglais**.
- **Unités internes :** Temps en **heures**, concentrations cellulaires en **10⁹ cellules/L**, concentrations médicament en **mg/L** (conversion vers µM : `mgL_to_uM = 1000/371.25`).
- **Unités pour les plots :** Temps en **jours** (`out$days <- out$time / 24`).
- **Nommage des fichiers :** `parameters_*.R` (jeux de paramètres), `plots_*.R` (visualisation), `run_*.R` (points d'entrée), `results_*/` (sorties).
- **Plots :** Axe Y en `log10`, `theme_bw`, palette de couleurs cohérente entre `plots.R`, `vpc.R`, `plots_grades.R`.

---

## Pièges courants

1. **Ordre des `source()` :** `parameters_FORNARI_CORRECT.R` doit être sourcé APRÈS le fichier de paramètres espèce (rat ou humain). Si sourcé en premier, les variables `init_pars`/`init_state` n'existent pas encore.

2. **Working directory :** Les scripts utilisent des `source()` et `read.csv()` avec chemins relatifs. **Toujours lancer depuis `scripts/`**, sinon les fichiers ne sont pas trouvés.

3. **Stiffness du solveur :** Le système peut devenir raide. `lsoda` est configuré avec `rtol=1e-7, atol=1e-9` pour les runs déterministes, et `rtol=1e-4, atol=1e-6` pour les boucles VPC. `maxsteps` peut aller jusqu'à 200 000.

4. **Kill dominance :** Si `Slope_X × Damage_max > 1`, l'effet cytotoxique dépasse la prolifération → nadirs irréalistement bas. `simulate_all()` affiche un diagnostic pour ce cas.

5. **Lambda amplification :** Les compartiments transit Ret et Plt sont prolifératifs (lambda=2). Les conditions initiales de T1_Ret et T1_Plt ne sont PAS identiques à T2/T3 (facteur lambda, Eq. S3).

6. **Calvert dosing (humain) :** Dose (mg) = AUC_cible × (GFR + 25). Patient de référence : GFR = 125 mL/min (Supplementary S11).

7. **Ordre de calibration :** Lancer `calibrate_neut.R` (Slope_CMP) **avant** `calibrate_plt.R` (delta_Plt, gamma_prolTrans). Les deux scripts réécrivent `parameters_human.R` à la convergence — garder une copie de sauvegarde. Si `calibrate_neut.R` change Slope_CMP, relancer ensuite `calibrate_plt.R` car delta_Plt peut légèrement interagir via MPP.

8. **VPC vs IIV :** Le VPC standard utilise uniquement l'erreur résiduelle log-additive (sigma de la Table S4), pas d'IIV. L'IIV (sur Slopes et baselines) n'est ajoutée que dans `save_grade_figure4c()` dans `plots_grades.R`.

---

## Données

Fichiers CSV dans `data/` : digitalisés depuis les figures de Fornari 2019 (WebPlotDigitizer).
- **Format :** 2 colonnes — `time` (jours), `value` (10⁹/L)
- **Données rat :** `MPP.csv`, `CMP.csv`, `MEP.csv`, `MEPF.csv`, `Neut.csv`, `MONO.csv`, `Plt.csv`, `Ret.csv`, `RBC.csv`
- **Données humaines :** `Neut_H.csv`, `Plt_H.csv`

---

## Roadmap du stage

| Phase | Statut | Description |
|-------|--------|-------------|
| 1. Reproduction Fornari | ✅ En cours | Validation rat (Fig. 3/S1) + transposition humain (Fig. 4) |
| 2. Application à de nouvelles molécules | 🔜 À venir | Adaptation du pipeline PK/PD + estimation des paramètres en R (`nlmixr2` ou `saemix`) |

---

## Références clés

- **Fornari et al. 2019** — CPT: Pharmacometrics & Systems Pharmacology (modèle principal)
- **Zandvliet et al. 2008** — Br J Clin Pharmacol 66:485-497 (PK humain carboplatine)
- **De Carlo et al. 2025** — Br J Clin Pharmacol (dosage de précision guidé par modèle)
- **NCI-CTCAE v5.0** — Seuils de grading de toxicité
