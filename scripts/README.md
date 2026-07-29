# Modèle PK/PD TGI — Simeoni (2004)
Script R générique pour l'ajustement séquentiel du modèle de croissance tumorale inhibée (TGI) de Simeoni (2004) sur données de xénogreffe souris, dose unique IV bolus.

---

## Utilisation rapide

```r
Rscript pkpd_simeoni_tgi.R
```

Sorties générées dans `scripts/` :
- `plot_PKPD_simeoni_<COHORTE>.png` — graphe volume tumoral observé vs prédit
- `resultats_PKPD_simeoni_<COHORTE>.RData` — tous les paramètres estimés

---

## Fichiers requis

| Fichier | Format | Contenu |
|---------|--------|---------|
| CSV PK | CSV | Colonnes : `Animal`, `CL`, `V1`, `V2`, `Q`, `cmax` (unités mL/kg et /h) |
| Excel tumeur | `.xlsx` | Volume tumoral moyen + SEM, format TALL ou WIDE (voir ci-dessous) |

### Format Excel TALL (temps en colonne 1)
```
Time (day) | Ctrl | Dose 1 | Dose 2 | ...
    0       |  182 |   182  |   182  |
    7       |  398 |   350  |   280  |
   ...
--- SEM ---
    0       |   12 |    12  |    12  |
   ...
```

### Format Excel WIDE (temps en ligne d'en-tête)
```
Group   |  0  |  7  | 14  | ...
Ctrl    | 182 | 398 | 812 |
Dose 1  | 182 | 350 | 600 |
--- SEM ---
Ctrl    |  12 |  28 |  45 |
...
```
Les deux formats sont détectés automatiquement. Dans les deux cas, 2 blocs sont attendus (moyennes puis SEM).

---

## Configuration

**Seule la section `CONFIG` (lignes 34–83) doit être modifiée** entre produits.

```r
# Fichiers
FICHIER_PK    <- "pk_resultats_20260720_150418.csv"
FICHIER_TUMOR <- "Tumor_volum_SNU.xlsx"

# Identité
NOM_PRODUIT <- "Fc-silent FGFR2-huBPA-LP1"
NOM_COHORTE <- "SNU"

# Groupes traités — autant d'éléments que de doses
NOM_DOSES    <- c("groupe 4", "grp 5", "grp 7")   # noms dans le CSV PK (colonne Animal)
DOSES_UGKG   <- c(1200, 5600, 11800)               # doses en µg/kg
W_DOSES      <- c(0, 1, 0.1)                       # poids dans le fit (0 = exclu)

# Position dans l'Excel (offset depuis la ligne d'en-tête)
LIGNE_CTRL   <- 1
LIGNES_DOSES <- c(2, 3, 4)

# Bornes biologiques pour l'Etape B
MTT_MIN_J <- 5
MTT_MAX_J <- 40
```

`NOM_DOSES`, `DOSES_UGKG`, `W_DOSES` et `LIGNES_DOSES` doivent tous avoir la même longueur — un `stopifnot()` le vérifie au démarrage.

---

## Le modèle

### PK — 2 compartiments IV bolus
```
dA1/dt = -(CL/V1 + Q/V1)·A1 + (Q/V2)·A2
dA2/dt =  (Q/V1)·A1 - (Q/V2)·A2
C1 = A1 / V1
```
Paramètres `CL`, `V1`, `V2`, `Q` fixés depuis le CSV PK (convertis mL/kg·/h → L/kg·/j).

### PD — Transit (Simeoni 2004)
```
dx1/dt = [g(w)/w − k2·C1]·x1
dx2/dt =  k2·C1·x1 − k1·x2
dx3/dt =  k1·x2 − k1·x3
dx4/dt =  k1·x3 − k1·x4

g(w) = L0·w / (1 + (L0·w/L1)^20)^(1/20)   [ψ = 20]
```

| Paramètre | Unité | Estimé / Dérivé |
|-----------|-------|-----------------|
| `L0` | /j | Estimé — taux de croissance exponentielle |
| `L1` | mm³/j | Estimé — plateau de croissance linéaire |
| `k1` | /j | Estimé — taux de transit (MTT = 4/k1) |
| `k2` | L/µg/j | Estimé — taux de cytotoxicité |
| `TSC` | µg/L | **Dérivé** = L0/k2 — concentration seuil |
| `MTT` | j | **Dérivé** = 4/k1 — temps de transit moyen |

---

## Ajustement séquentiel

**Étape A** — Véhicule uniquement : estime `L0` et `L1`

**Étape B** — Groupes traités : estime `k1` et `k2` avec `L0`/`L1` fixés

Fonction de coût (inverse-variance, échelle log) :
```
cost = Σᵢ Wᵢ × Σₜ [1/SEMᵢₜ² × (log(obs) − log(pred))²]
```

Optimiseur : DEoptim (global) suivi de nlminb (local) — le meilleur des deux est retenu.

Les timepoints sans SEM valide sont **exclus** du fit (évite un biais de poids artificiel).

### Poids W_DOSES
- `W = 1.0` → poids normal dans le fit
- `W = 0.1` → poids réduit (groupe moins prioritaire)
- `W = 0`   → **exclu du fit** (courbe prédite affichée quand même)

Utile quand un seul `k2` ne peut pas fitter toutes les doses simultaneously (structural limitation).

---

## Avertissements console

| Message | Cause | Action |
|---------|-------|--------|
| `[AVERT] k2 sur borne inférieure` | TSC très élevé — drug trop faible | Augmenter `K2_TSC_MAX_X` |
| `[AVERT] k2 sur borne supérieure` | TSC très bas — drug très puissant | Augmenter `K2_TSC_MIN_DIV` |
| `[AVERT] k1 sur borne sup/inf` | MTT sur limite | Élargir `MTT_MIN_J` / `MTT_MAX_J` |
| `Groupe 'X' introuvable dans le CSV PK` | Nom dans `NOM_DOSES` incorrect | Vérifier colonne `Animal` du CSV |

---

## Dépendances R

```r
install.packages(c("deSolve", "DEoptim", "ggplot2", "readxl"))
```

Sur Debian/Ubuntu sans accès CRAN :
```bash
apt-get install r-cran-desolve r-cran-deoptim r-cran-ggplot2 r-cran-readxl
```

---

## Ne jamais modifier

Les fonctions suivantes implémentent les équations de Simeoni — les modifier casserait le modèle :
- `simeoni_rhs` — ODEs groupes traités
- `simeoni_ctrl_rhs` — ODEs véhicule
- `sim_ctrl_fn` — simulation véhicule
- `sim_treated` — simulation groupe traité

---

## Résultats SNU — Fc-silent FGFR2-huBPA-LP1

| Paramètre | Valeur | Source |
|-----------|--------|--------|
| L0 | 0.125 /j | Étape A |
| L1 | 1 372 mm³/j | Étape A |
| TSC | 7 927 µg/L | Étape B |
| MTT | 17.6 j | Étape B |
| k2 | 1.58×10⁻⁵ L/µg/j | Étape B |

Configuration utilisée : `W_DOSES = c(0, 1, 0.1)` — 1.2 mg/kg exclu du fit, 5.6 mg/kg prioritaire.

---

*Référence : Simeoni M. et al., Cancer Research, 2004.*
