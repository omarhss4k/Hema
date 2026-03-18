############################################################
# parameters_human.R
# Paramètres physiologiques de base — HUMAIN
# Baselines, MTTs, circulating rates, PK carboplatin
# Référence : Fornari 2019, Table 1 & Calvert formula
############################################################

init_pars <- list()

# ── PK carboplatin humain (modèle 2 compartiments) ──
# CL = GFR + 25 mL/min (formule Calvert, GFR médian = 78 mL/min)
# CL = 103 mL/min = 6.18 L/h
init_pars$CL <- 6.18   # L/h
init_pars$V1 <- 17.0   # L
init_pars$Q  <- 7.5    # L/h
init_pars$V2 <- 25.0   # L

# ── Fraction libre (carboplatin → liaison protéique négligeable) ──
init_pars$fu0    <- 1.0
init_pars$fu_inf <- 1.0
init_pars$k_bind <- 0.0

# ── Damage (Table 1, mêmes que rat) ──
init_pars$k_dam <- 0.017
init_pars$k_rep <- 0.017

# ── Baselines humains (Table 1, 10⁹ cells/L) ──
init_pars$MPP0  <-    6.2
init_pars$CMP0  <-   14.3
init_pars$MEP0  <-    6.3
init_pars$Neut0 <-    4.5
init_pars$Mono0 <-    0.5
init_pars$Ret0  <-   70.0
init_pars$RBC0  <- 5000.0
init_pars$Plt0  <-  250.0

# ── MTTs humains (Table 1, en heures) ──
init_pars$MTT_Neut <- 130.0   # granulopoïèse : ~5-6 jours
init_pars$MTT_Mono <-  80.0
init_pars$MTT_Ret  <-  72.0   # réticulocytes : ~3 jours
init_pars$MTT_Plt  <- 200.0   # thrombopoïèse : ~8-10 jours

# ── Circulating rates (Table 1, en h⁻¹) ──
init_pars$k_circ_Neut <- 0.130   # demi-vie ~7.6 h
init_pars$k_circ_Mono <- 0.020   # demi-vie ~2 jours
init_pars$k_circ_Plt  <- 0.0042  # durée de vie ~10 jours
init_pars$k_circ_RBC  <- 0.00035 # durée de vie ~120 jours
# k_circ_Ret est dérivé dans parameters_FORNARI_CORRECT.R (Eq S4)

# ── Drug effects (Table 2) ──
init_pars$Slope_MPP  <- 2.05
init_pars$Slope_CMP  <- 1.47
init_pars$Slope_MEP  <- 2.19
init_pars$delta_Ret  <- 2.8
init_pars$delta_Plt  <- 1.00  # calibration auto : nadir1~150 nadir2~130
init_pars$delta_Neut <- 0.003

# ── Feedback powers (Table 1) ──
init_pars$gamma_stem      <- 0.07
init_pars$gamma_mat_CMP   <- 0.60
init_pars$gamma_mat_MEP   <- 0.30
init_pars$gamma_prolTrans <- 0.55

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
