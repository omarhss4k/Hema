# ADC-HematoTox — outil de prédiction de la myélotoxicité des ADC

Outil semi-mécanistique (PK/PD, modèle de Fornari) pour **prédire la
myélosuppression** (neutrophiles + monocytes) d'un anticorps-médicament
conjugué (ADC) à partir de sa **PK** et de son **IC50**, sans nécessiter
d'étude in vivo pour le nouveau composé.

---

## Principe en 2 temps

```
   ADC de REFERENCE                          NOUVEAU composé
   (PK + IC50 + profils)                     (PK + IC50 seulement)
          │                                          │
          ▼                                          ▼
   scripts/calibrate.R   ──calibration.rds──▶  scripts/predict.R
   (cale le modèle)                            (prédit sa toxicité)
```

1. **Calibrer** le modèle sur un ADC de référence dont on a toutes les
   données (PK, IC50, profils hémato observés).
2. **Prédire** un nouveau composé : on branche sa PK + son IC50, la
   *potency* est translatée depuis la référence via le rapport des IC50.

> La prédiction suppose un mécanisme **cytotoxique** : le payload tue les
> cellules en division (progéniteurs myéloïdes).

---

## Installation

- **R** (≥ 4.0) avec le package **deSolve** :
  ```r
  install.packages("deSolve")
  ```

---

## Utilisation

Toujours lancer **depuis la racine de l'outil** (`adc_hematox_tool/`).

### Sous RStudio / Posit Workbench
Ouvre le dossier `adc_hematox_tool/` **comme projet RStudio** (ou
`Session ▸ Set Working Directory ▸ To Source File Location`), puis dans la
console :
```r
source("scripts/predict.R")     # ou calibrate.R / run_all.R
```
Le répertoire de travail doit être la racine de l'outil. Les scripts ne
contiennent **aucun `quit()`** : ils ne fermeront jamais ta session.

### 1. Prédire un composé (cas courant)
```bash
Rscript scripts/predict.R
```
→ `results/prediction_<nom>.pdf` + grades CTCAE dans la console.
Utilise la calibration existante (ou les valeurs par défaut si aucune).

### Importer des données FGFR2 comme référence
Si ta référence est l'étude FGFR2 (fichier `nhp_hema_data.csv` avec colonnes
`Neut_1e3_uL`, `Mono_1e3_uL`, …), place-le dans `data/` puis :
```bash
Rscript scripts/import_fgfr2.R      # -> data/reference_profiles.csv
```
Le script mappe les colonnes (Neut/Mono en 10³/µL = 10⁹/L, copie directe) et
affiche les baselines moyennes à reporter dans `config/compounds.R`.

### 2. (Re)calibrer sur l'ADC de référence
```bash
Rscript scripts/calibrate.R
```
→ `results/calibration.rds` (lue automatiquement par `predict.R`)
+ `results/calibration_fit.pdf` (validation vs profils observés).
Nécessite les profils de référence dans `data/` (voir plus bas).

### 2 bis. Calibration AUTOMATIQUE (optionnel, recommandé)
```bash
Rscript scripts/autocalibrate.R
```
Ajuste automatiquement les 4 paramètres PD clés (`IC50ADC_scale`, `Emax_CMP`,
`k_depl_direct`, `k_rep`) pour **minimiser l'écart moindres-carrés** modèle↔profils
observés, au lieu de les régler à la main. Sauvegarde la calibration optimisée
dans `results/calibration.rds` (lue par `predict.R`), affiche les valeurs
avant→après (à recopier dans `config/compounds.R`) et trace `calibration_fit.pdf`.

> ⚠️ **Identifiabilité** : avec peu de points, plusieurs jeux de paramètres
> peuvent coller aux données. Les paramètres sont **bornés** à des plages
> physiques pour éviter les solutions aberrantes, mais **regarde toujours la
> figure** : si un paramètre bute sur sa borne, le mécanisme est mal contraint.

### 3. Tout enchaîner (calage → prédiction)
```bash
Rscript scripts/run_all.R
```
Lance `calibrate.R` puis `predict.R` d'affilée.

### Sorties
- `results/prediction_<nom>.pdf` — courbes neutro + monocytes
- `results/prediction_<nom>_grades.csv` — tableau (dose × IC50 → nadir,
  jour, valeur absolue, grade CTCAE) réutilisable dans Excel
- `results/calibration_fit.pdf` — validation du calage
- `results/calibration.rds` — calibration persistée

---

## Le seul fichier à éditer : `config/compounds.R`

Tout se règle là, sans toucher au moteur :

- **`REFERENCE`** — PK, IC50 et fichier de profils de l'ADC de référence.
- **`CALIB`** — paramètres PD calibrés (profondeur/timing du nadir). Pour
  re-caler, ajuste-les et relance `calibrate.R`.
- **`NEW_COMPOUND`** — PK + IC50 du composé à prédire. `IC50_myelo` accepte
  **plusieurs valeurs** (fourchette d'incertitude → plusieurs courbes).

### Allométrie (PK d'une autre espèce)
Si la PK du composé a été mesurée dans une **autre espèce** que la cible
(ex. PK **souris**, prédiction **NHP**), mets dans `NEW_COMPOUND` :
```r
allometry = TRUE,
pk_bw     = 0.02,   # BW (kg) de l'espèce où la PK a été mesurée
```
Le modèle scale automatiquement : **CL, Q ∝ BW^0.75** et **V1, V2 ∝ BW^1.0**
(Boxenbaum), de `pk_bw` vers `BW_kg`.

> ⚠️ Vérifie d'abord ta PK source : un V1 de souris (~20 g) doit être
> **~0.002–0.004 L**. Si le V1 source est bien plus grand, le fit amont est
> faux (souvent unités de concentration) et l'allométrie donnera des volumes
> aberrants. Un **garde-fou** alerte si le V scalé est implausible.

### Unités (important)
- PK : **CL, Q en L/h** ; **V1, V2 en L** ; doses en **mg/kg**.
- Concentrations (pour l'ajustement PK amont) : **µg/mL** (= mg/L).
  Si ton dosage est en ng/mL → **diviser par 1000**.
- IC50 : **nM**.

---

## Données confidentielles

Les profils hémato observés (référence) vont dans **`data/`**, au format :

| Animal_Id | Dose_mgkg | time_day | Neut | Mono | ... |
|-----------|-----------|----------|------|------|-----|

Fichier attendu : `data/reference_profiles.csv` (chemin réglable dans
`config/compounds.R`). **`data/` et `results/` sont gitignorés** — aucune
donnée confidentielle n'est versionnée.

---

## Structure

```
adc_hematox_tool/
├── README.md
├── config/compounds.R      ← À ÉDITER : composés + calibration
├── model/                  ← moteur (ne pas modifier)
│   ├── parameters_baseline.R   biologie de base (baselines PD)
│   ├── parameters_fornari.R    feedbacks Fornari
│   ├── parameters_adc.R        constantes PK/PD ADC
│   ├── parameters_infusion.R   fonction de perfusion
│   ├── model_ode.R             système d'ODE (Fornari + déplétion directe)
│   └── engine.R                fonctions haut niveau (load/build/simulate)
├── scripts/
│   ├── calibrate.R         1) cale sur la référence
│   └── predict.R           2) prédit un nouveau composé
├── data/                   profils confidentiels (gitignored)
└── results/                sorties (gitignored)
```

---

## Ce que le modèle capture / ne capture pas

**Capture** : myélosuppression dose-dépendante (neutro + monocytes),
timing du nadir, sélectivité myéloïde (érythroïde épargné selon l'IC50),
grades CTCAE.

**Ne capture pas** : la neutrophilie de stress précoce (un cytotoxique ne
fait que déprimer) et les points terminaux/sacrifice.

---

## Limites & bonnes pratiques

- La **valeur absolue** de l'IC50 pilote la prédiction : en cas
  d'incertitude (ex. 2 mesures), lancer avec les deux → fourchette.
- La calibration porte une **constante d'intensité** (`CALIB`) calée sur la
  référence ; extrapoler à un composé de linker/DAR très différent réclame
  une re-calibration.
- Modèle **descriptif** sur la référence, **prédictif** sur un nouveau
  composé via translation d'IC50. Pour une prédiction *validée*, confronter
  a posteriori aux données in vivo quand elles arrivent.
