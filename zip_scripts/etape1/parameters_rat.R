############################################################
# parameters_rat.R
# Paramètres physiologiques de base -- RAT
# Baselines, MTTs, circulating rates, fu, damage
############################################################

init_pars <- list()

# -- PK (allométrie rat 250g depuis humain) --
BW_h <- 70; BW_r <- 0.25
init_pars$CL <- 0.42 
init_pars$V1 <- 0.06
init_pars$Q  <- 0.26
init_pars$V2 <- 0.30

# -- Fraction libre
# Note : platine libre déjà modélisé dans C1 → fu = 1
init_pars$fu0    <- 1.0
init_pars$fu_inf <- 1.0
init_pars$k_bind <- 0.0

# -- Damage (Table 1) --
init_pars$k_dam <- 0.017
init_pars$k_rep <- 0.017

# -- Baselines (Table 1) --
init_pars$MPP0  <- 69.0
init_pars$CMP0  <- 214.0
init_pars$MEP0  <- 156.0
init_pars$Neut0 <- 1.32
init_pars$Mono0 <- 0.16
init_pars$Ret0  <- 235.0
init_pars$RBC0  <- 8161.0
init_pars$Plt0  <- 967.0

# -- MTTs (Table 1) --
init_pars$MTT_Neut <- 61.5
init_pars$MTT_Mono <- 62.2
init_pars$MTT_Ret  <- 75.0
init_pars$MTT_Plt  <- 74.0

# -- Circulating rates (Table 1) --
init_pars$k_circ_Neut <- 0.17
init_pars$k_circ_Mono <- 0.06
init_pars$k_circ_Plt  <- 0.01
init_pars$k_circ_RBC  <- 0.0007
# k_circ_Ret est dérivé dans parameters_FORNARI_CORRECT.R (Eq S4)

# -- Drug effects (Table 2) --
init_pars$Slope_MPP <- 2.05
init_pars$Slope_CMP <- 1.47
init_pars$Slope_MEP <- 2.19
init_pars$delta_Ret <- 2.8
init_pars$delta_Plt <- 0.54

# -- IC50 colony-forming unit assays (Table S2) --
init_pars$IC50_CMP_rat   <- 2.83  # μM, CD45+ cells
init_pars$IC50_CMP_human <- 4.04  # μM, CD45+ cells
init_pars$IC50_MEP_rat   <- 0.86  # μM, CD71+ cells
init_pars$IC50_MEP_human <- 1.56  # μM, CD71+ cells

# -- Feedback powers (Table 1) --
init_pars$gamma_stem      <- 0.07
init_pars$gamma_mat_CMP   <- 0.60
init_pars$gamma_mat_MEP   <- 0.30
init_pars$gamma_prolTrans <- 0.70

# -- MW carboplatin --
init_pars$MW_carboplatin <- 371.25
init_pars$mgL_to_uM      <- 1000 / init_pars$MW_carboplatin

# -- Placeholder état initial (complété dans parameters_FORNARI_CORRECT.R) --
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

# -- Fonction perfusion répétée --
make_repeated_infusion <- function(dose_mg, Tinfu_h = 1, interval_h, n_cycles) {
  rate     <- dose_mg / Tinfu_h
  t_starts <- seq(0, by = interval_h, length.out = n_cycles)
  function(t) {
    if (any(t >= t_starts & t < (t_starts + Tinfu_h))) rate else 0
  }
}
