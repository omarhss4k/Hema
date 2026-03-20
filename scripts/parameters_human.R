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
init_pars$k_dam <- 0.017
init_pars$k_rep <- 0.017

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
init_pars$MTT_Neut <- 210.0   # scaled from ref. 5
init_pars$MTT_Mono <- 121.5   # scaled from ref. 47
init_pars$MTT_Ret  <-  66.0   # scaled from internal AZ study
init_pars$MTT_Plt  <- 168.0   # scaled from ref. 5

# ── Circulating rates (Table 1, en h⁻¹) ──
init_pars$k_circ_Neut <- 0.100    # ref. 22
init_pars$k_circ_Mono <- 0.040    # ref. 41
init_pars$k_circ_Plt  <- 0.0052   # ref. 5
init_pars$k_circ_RBC  <- 0.00037  # ref. 13
# k_circ_Ret est dérivé dans parameters_FORNARI_CORRECT.R (Eq S4)

# ── Drug effects (Table 2) ──
# Slopes ajustés pour la sensibilité espèce-spécifique (Eq. 10)
init_pars$Slope_MPP  <- 0.79   # IC50-scaled from rat (2.05)
init_pars$Slope_CMP  <- 0.57   # IC50-scaled from rat (1.47)
init_pars$Slope_MEP  <- 0.66   # IC50-scaled from rat (2.19)
init_pars$delta_Ret  <- 2.8    # same as rat
init_pars$delta_Plt  <- 0.80  # calibration auto : nadir1~160 nadir2~146
# delta_Neut absent du modèle Fornari (pas de drug effect sur transit neutrophiles)

# ── Feedback powers (Table 1, same as rat) ──
init_pars$gamma_stem      <- 0.07
init_pars$gamma_mat_CMP   <- 0.60
init_pars$gamma_mat_MEP   <- 0.30
init_pars$gamma_prolTrans <- 0.40

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
