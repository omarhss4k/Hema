# Hema — PK/PD Modeling of Carboplatin-Induced Hematological Toxicity

> Stage de fin d'études M2 Science des Données de Santé — Pierre Fabre, 2026
> **Omar Hassak**

## Description

Ce projet reproduit et étend le modèle PK/PD semi-mécaniste de **Fornari et al. 2019** (CPT: Pharmacometrics & Systems Pharmacology) pour prédire la toxicité hématologique induite par le carboplatine. Le modèle est d'abord validé sur données de rat, puis transposé à l'humain via une mise à l'échelle allométrique.

**Pipeline translationnel :** Rat (validation) → Humain (prédiction)

Le modèle décrit la dynamique de **8 lignées cellulaires** (MPP, CMP, MEP, Neutrophiles, Monocytes, Réticulocytes, GR, Plaquettes) via un système de **23 équations différentielles ordinaires (ODE) couplées**, incluant :
- Un modèle PK à 2 compartiments (plasma + tissus)
- Un modèle de dommage ADN (liaison du platine libre)
- Des compartiments de transit hématopoïétiques (moelle osseuse → sang périphérique)
- Des boucles de rétroaction physiologiques (Eq. 11-13, Fornari 2019)

---

## Installation

### Prérequis R (>= 4.0)

```r
# Lancer depuis R ou via le script dédié
source("scripts/install_deps.R")
```

Packages requis : `deSolve`, `ggplot2`, `gridExtra`, `grid`, `scales`

### Prérequis Python (>= 3.8)

```bash
pip install -r requirements.txt
```

---

## Utilisation rapide

Tous les scripts R doivent être exécutés **depuis le répertoire `scripts/`**.

```bash
# 1. Validation rat (Figure 3 + Figure S1 + VPC + grading NCI-CTCAE)
cd scripts && Rscript run_CORRECT.R

# 2. Simulation humaine (Figure 4 + VPC + distribution de grades)
cd scripts && Rscript run_human.R

# 3. Générer les slides de réunion
python3 reunion_point/make_pptx.py
```

Les résultats (PDFs) sont sauvegardés dans :
- `scripts/results_CORRECT/` — Rat
- `scripts/results_HUMAN/` — Humain

---

## Structure du projet

```
Hema/
├── data/                          # Données expérimentales digitalisées (CSV)
│   ├── MPP.csv, CMP.csv, MEP.csv  # Progéniteurs (rat)
│   ├── Neut.csv, MONO.csv         # Lignée myéloïde (rat)
│   ├── Ret.csv, RBC.csv, Plt.csv  # Lignées érythroïde/plaquettes (rat)
│   └── Neut_H.csv, Plt_H.csv      # Données humaines
├── documents/                     # Publications de référence (PDFs)
├── scripts/                       # Scripts R de modélisation
│   ├── pkpd_model_FORNARI.R       # Système ODE (Eq. 1-13)
│   ├── parameters_rat.R           # Paramètres rat (Tables 1-2)
│   ├── parameters_human.R         # Paramètres humain (allométrie)
│   ├── parameters_FORNARI_CORRECT.R # Paramètres dérivés (Eq. S4)
│   ├── data_fornari.R             # Chargement des données CSV
│   ├── run_CORRECT.R              # Exécution principale (rat)
│   ├── run_human.R                # Exécution principale (humain)
│   ├── plots.R, plots_human.R     # Visualisation des simulations
│   ├── plots_grades.R             # Grading NCI-CTCAE v5.0
│   ├── vpc.R                      # Visual Predictive Check (1000 sims)
│   ├── calibrate_plt.R            # Calibration paramètres plaquettaires
│   ├── results_CORRECT/           # Sorties PDF (rat)
│   └── results_HUMAN/             # Sorties PDF (humain)
├── reunion_point/                 # Matériel pour les réunions
│   ├── make_pptx.py               # Génération des slides
│   └── img/                       # Figures pour les présentations
├── make_pptx.py                   # Générateur de présentation principal
├── requirements.txt               # Dépendances Python
├── CLAUDE.md                      # Instructions pour Claude Code
└── README.md                      # Ce fichier
```

---

## Modèle PK/PD

### Pharmacocinétique

Modèle à 2 compartiments avec infusion IV :
- **C1** (plasma central) — **C2** (compartiment périphérique)
- Paramètres : CL, V1, V2, Q (spécifiques à l'espèce)
- Fraction libre de platine : `fu(t) = fu_inf + (fu0 - fu_inf) × exp(-k_bind × t)`

### Effet du médicament

Dommage ADN : `dDamage/dt = k_dam × C_free_µM - k_rep × Damage`

Kill sur les progéniteurs (MPP, CMP, MEP) : `Kill = Slope × Damage`

### Dosage humain (formule de Calvert)

```
Dose (mg) = AUC_cible × (DFG + 25)
```

DFG de référence : 125 mL/min (Supplementary S11, Fornari 2019)

---

## Grading NCI-CTCAE v5.0

| Grade | Neutropénie (10⁹/L) | Thrombocytopénie (10⁹/L) |
|-------|---------------------|--------------------------|
| G1 | < 2.0 | < 150 |
| G2 | < 1.5 | < 75 |
| G3 | < 1.0 | < 50 |
| G4 | < 0.5 | < 25 |

---

## Références

1. **Fornari et al. 2019** — A Mathematical Model of Hematological Toxicities Induced by Repeated Cycles of Chemotherapy. *CPT: Pharmacometrics & Systems Pharmacology*, 8(9):653-663.
2. **Zandvliet et al. 2008** — Semiphysiological and empirical modelling of the population pharmacokinetics of the prodrug eniluracil and its metabolite uracil. *Br J Clin Pharmacol*, 66(4):485-497.
3. **De Carlo et al. 2025** — Model-informed precision dosing of carboplatin. *Br J Clin Pharmacol*.
4. **NCI-CTCAE v5.0** — Common Terminology Criteria for Adverse Events.
