############################################################
# parameters_human.R
# Paramètres physiologiques de base — HUMAIN
# Référence : Fornari 2019, Table 1 & Table 2
# Tous les paramètres vérifiés contre le papier
############################################################

init_pars <- list()

# ── PK carboplatin humain (modèle 2 compartiments) ──
# Ref. 22 = Zandvliet et al. 2008 (Br J Clin Pharmacol 66:485-497)
# CL calvert : GFR + 25 mL/min → 103 mL/min = 6.18 L/h (GFR médian 78)
init_pars$CL <- 6.18   # L/h  (= GFR+25 mL/min = 103 mL/min)
init_pars$V1 <- 7.87   # L    (Zandvliet 2008)
init_pars$Q  <- 1.98   # L/h  (= 33 mL/min, Zandvliet 2008)
init_pars$V2 <- 8.06   # L    (Zandvliet 2008)

# ── Fraction libre (carboplatin → liaison protéique négligeable) ──
init_pars$fu0    <- 1.0
init_pars$fu_inf <- 1.0
init_pars$k_bind <- 0.0

# ── Damage ──
# k_dam : même que rat (formation d'adduits ADN)
# k_rep : identique au rat — Table 1, Fornari 2019 : "As in the rat"
init_pars$k_dam <- 0.010   # calibration Fornari 2019 Fig.4c : G3=12%, G4=3% ✓
init_pars$k_rep <- 0.010   # demi-vie damage = 69h → fenêtre suffisante pour G3 Plt

# ── Baselines humains (Table 1, 10⁹ cells/L) ──
init_pars$MPP0  <-    1.3    # allometric scaling ref. 42
init_pars$CMP0  <-   20.9   # proportions ref. 43
init_pars$MEP0  <-   15.0   # proportions ref. 43
init_pars$Neut0 <-    4.5   # range 2–7, ref. 22
init_pars$Mono0 <-    0.5   # range 0.2–10, ref. 41
init_pars$Ret0  <-   70.0   # range 40–115, ref. 1
init_pars$RBC0  <- 5000.0   # range 4100–5900, ref. 32
init_pars$Plt0  <-  345.0   # ref. 5

# ── MTTs humains (Table 1, en heures) ──
# MTT_Neut calibré sur données Schmitt 2010 (nadir j14 cycle 1) :
# Scan MTT_Neut 114.7→250h : 200h donne nadir=2.53 @j13.1 ≈ données 2.52 @j14 ✓
# La correction théorique ×3/4 (114.7h) donnait nadir trop précoce (j10) et trop profond.
# MTT_Plt : correction S10 ×3/4 maintenue (Plt nadir j14 bien calibré : 175 vs données 170).
init_pars$MTT_Neut <- 200.0   # calibré Schmitt 2010 : nadir j13-14 ✓ (vs 114.7h → j10)
init_pars$MTT_Mono <-  91.1   # (3/4) × 121.5h (ref. 47, scaled)
init_pars$MTT_Ret  <-  49.5   # (3/4) × 66.0h  (AZ internal)
init_pars$MTT_Plt  <- 131.5   # (3/4) × 175.3h (Schmitt/Friberg, S10)

# ── Circulating rates (Table 1, en h⁻¹) ──
init_pars$k_circ_Neut <- 0.100    # ref. 22
init_pars$k_circ_Mono <- 0.040    # ref. 41
init_pars$k_circ_Plt  <- 0.0052   # ref. 5
init_pars$k_circ_RBC  <- 0.00037  # ref. 13
# k_circ_Ret est dérivé dans parameters_FORNARI_CORRECT.R (Eq S4)

# ── Drug effects (Table 2) ──
# Slopes ajustés pour la sensibilité espèce-spécifique (Eq. 10)
# Eq. 10 : Slope_human = Slope_rat × (IC50_rat / IC50_human)
# IC50_CMP: rat=2.83µM, human=4.04µM → Slope_MPP = 2.05×(2.83/4.04) = 1.435
#           Slope_CMP = 1.47×(2.83/4.04) = 1.029
# IC50_MEP: rat=0.86µM, human=1.56µM  → Slope_MEP = 2.19×(0.86/1.56) = 1.208
init_pars$Slope_MPP  <- 1.435  # IC50-scaled from rat (2.05 × 2.83/4.04)
init_pars$Slope_CMP  <- 1.029  # IC50-scaled from rat (1.47 × 2.83/4.04)
init_pars$Slope_MEP  <- 1.208  # IC50-scaled from rat (2.19 × 0.86/1.56)
init_pars$delta_Ret  <- 2.8    # same as rat (Table 2)
init_pars$delta_Plt  <- 0.54  # Fornari 2019 Fig.4c : G3=12% avec k_rep=0.010

# ── IC50 colony-forming unit assays (Table S2) ──
# Utilisés pour scaler les Slope via Eq. 10 : Slope_H = Slope_R × (IC50_R/IC50_H)
init_pars$IC50_CMP_rat   <- 2.83  # μM, CD45+ cells
init_pars$IC50_CMP_human <- 4.04  # μM, CD45+ cells
init_pars$IC50_MEP_rat   <- 0.86  # μM, CD71+ cells
init_pars$IC50_MEP_human <- 1.56  # μM, CD71+ cells

# ── Feedback powers (Table 1, same as rat) ──
init_pars$gamma_stem      <- 0.07
init_pars$gamma_mat_CMP   <- 0.60
init_pars$gamma_mat_MEP   <- 0.30
# gamma_prolTrans — valeur rat Table 1 Fornari 2019 :
# gamma=0.50 avec k_rep=0.010, delta=0.54 reproduit la distribution des grades Fig.4c
# (G3 Plt=12%, G4=3%) — calibration prioritaire sur la courbe temporelle
init_pars$gamma_prolTrans <- 0.50

# ── MW carboplatin ──
init_pars$MW_carboplatin <- 371.25
init_pars$mgL_to_uM      <- 1000 / init_pars$MW_carboplatin

# ── Placeholder état initial (complété dans parameters_FORNARI_CORRECT.R) ──
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

# ── Fonction perfusion répétée ──
make_repeated_infusion <- function(dose_mg, Tinfu_h = 1, interval_h, n_cycles) {
  rate     <- dose_mg / Tinfu_h
  t_starts <- seq(0, by = interval_h, length.out = n_cycles)
  function(t) {
    if (any(t >= t_starts & t < (t_starts + Tinfu_h))) rate else 0
  }
}
