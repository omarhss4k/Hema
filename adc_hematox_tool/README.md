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

> La prédiction suppose un mécanisme **cytotoxique** (le payload tue les
> cellules en division). Pour un ADC **immuno/STING**, la toxicité est
> immuno-médiée : la prédiction cytotoxique devient un **plancher**
> (la toxicité réelle peut être supérieure).

---

## Installation

- **R** (≥ 4.0) avec le package **deSolve** :
  ```r
  install.packages("deSolve")
  ```

---

## Utilisation

Toujours lancer **depuis la racine de l'outil** (`adc_hematox_tool/`).

### 1. Prédire un composé (cas courant)
```bash
Rscript scripts/predict.R
```
→ `results/prediction_<nom>.pdf` + grades CTCAE dans la console.
Utilise la calibration existante (ou les valeurs par défaut si aucune).

### 2. (Re)calibrer sur l'ADC de référence
```bash
Rscript scripts/calibrate.R
```
→ `results/calibration.rds` (lue automatiquement par `predict.R`)
+ `results/calibration_fit.pdf` (validation vs profils observés).
Nécessite les profils de référence dans `data/` (voir plus bas).

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
fait que déprimer), les points terminaux/sacrifice, et — pour un composé
STING — l'amplification immunitaire (→ prédiction = plancher).

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
