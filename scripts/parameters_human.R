############################################################
# parameters_human.R
#
# Paramètres physiologiques de base — HOMME
# Baselines, MTTs, circulating rates, fu, damage
# Source principale : Fornari et al., 2019 (Tables 1–2)
############################################################

init_pars <- list()

# ── PK (valeurs humaines de référence) ──
# Les valeurs ci-dessous sont celles utilisées comme références
# dans la modélisation et l’allométrie (CL, V1, Q, V2) pour l’adulte ~70 kg. [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
BW_h <- 70
init_pars$CL <- 6.15   # L/h   (humain) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$V1 <- 16    # L     (humain) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$Q  <- 6     # L/h   (humain) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$V2 <- 20    # L     (humain) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)

# ── Fraction libre
# Carboplatine : fraction libre identique entre espèces → pas
# d’ajustement interespèces requis pour fu ; ici fu = 1. [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$fu0    <- 1.0
init_pars$fu_inf <- 1.0
init_pars$k_bind <- 0.0

# ── Damage (identique rat ↔ homme, Eq. 3) ──
# Cinétiques de création/réparation du dommage (adduits Pt–ADN). [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$k_dam <- 0.017   # 1/h
init_pars$k_rep <- 0.017   # 1/h

# ── Baselines (Table 1, unités 10^9 cellules/L) ──
# MPP, CMP, MEP humains issus d’allométrie/proportions ; autres baselines de la littérature. [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$MPP0  <- 1.3     # MPP0 (allo. depuis rat) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$CMP0  <- 20.8    # Proportions humaines vs MPP [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$MEP0  <- 14.3    # Proportions humaines vs MPP [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)

# Pour les lignées circulantes, Fornari Table 1 fournit des plages typiques.
# Ci-dessous, on fixe des valeurs centrales/pratiques cohérentes.
init_pars$Neut0 <- 4.3     # 2–7 ref. (valeur médiane pratique) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$Mono0 <- 0.4     # ~0.2–1.0 typique; Table 1 montre 0.04–? (OCR), on choisit 0.4 [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$Ret0  <- 80.0    # 40–115 (valeur médiane) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$RBC0  <- 5000.0  # 4100–5900 (valeur médiane) [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$Plt0  <- 360.0   # valeur reportée pour l’homme [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)

# ── MTTs (Table 1, heures) ──
# Définis/obtenus par mise à l’échelle humaine des durées de maturation. [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$MTT_Neut <- 210.0
init_pars$MTT_Mono <- 121.5
init_pars$MTT_Ret  <- 66.0
init_pars$MTT_Plt  <- 168.0

# ── Circulating rates (Table 1, 1/h) ──
# Dérivés des demi‑vies de circulation chez l’homme. [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$k_circ_Neut <- 0.10
init_pars$k_circ_Mono <- 0.04
init_pars$k_circ_Plt  <- 0.0052
init_pars$k_circ_RBC  <- 0.00037
# k_circ_Ret est dérivé dans votre script de paramètres détaillé (Eq S4). [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)

# ── Drug effects (Table 2, ajustés pour la sensibilité humaine) ──
# Les pentes (Slope) humaines proviennent du scaling par IC50 (Eq. 10).
# δ_Ret et δ_Plt conservés identiques au rat. [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$Slope_MPP <- 0.79  # 1/µM
init_pars$Slope_CMP <- 0.57  # 1/µM
init_pars$Slope_MEP <- 0.66  # 1/µM
init_pars$delta_Ret <- 2.8
init_pars$delta_Plt <- 2
init_pars$delta_Neut <- 0.003   # à ajuster si nadir trop bas

# ── Feedback powers (identiques rat ↔ homme ; Table 1) ──
init_pars$gamma_stem      <- 0.07  # [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$gamma_mat_CMP   <- 0.60  # [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$gamma_mat_MEP   <- 0.30  # [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)
init_pars$gamma_prolTrans <- 0.60  # [1](https://pfabre-my.sharepoint.com/personal/pe155390_pierre-fabre_com/Documents/Fichiers%20Microsoft%20Copilot%20Chat/Fornari%202019%20-%20Quantifying_Drug-Induced_Bone_Marrow_Toxicity_Usin.pdf)

# ── Masse molaire carboplatine ──
init_pars$MW_carboplatin <- 371.25
init_pars$mgL_to_uM      <- 1000 / init_pars$MW_carboplatin

# ── État initial (placeholder ; peut être réajusté via steady‑state solver) ──
init_state <- c(
  C1 = 0, C2 = 0, Damage = 0,
  MPP = init_pars$MPP0,
  CMP = init_pars$CMP0,
  MEP = init_pars$MEP0,
  T1_Neut = 0, T2_Neut = 0, T3_Neut = 0, Neut = init_pars$Neut0,
  T1_Mono = 0, T2_Mono = 0, T3_Mono = 0, Mono = init_pars$Mono0,
  T1_Ret  = 0, T2_Ret  = 0, T3_Ret  = 0, Ret  = init_pars$Ret0,
  RBC = init_pars$RBC0,
  T1_Plt  = 0, T2_Plt  = 0, T3_Plt  = 0, Plt  = init_pars$Plt0
)

# ── Fonction perfusion répétée (identique) ──
make_repeated_infusion <- function(dose_mg, Tinfu_h = 1, interval_h, n_cycles) {
  rate     <- dose_mg / Tinfu_h
  t_starts <- seq(0, by = interval_h, length.out = n_cycles)
  function(t) {
    if (any(t >= t_starts & t < (t_starts + Tinfu_h))) rate else 0
  }
}