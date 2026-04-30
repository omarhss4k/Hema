############################################################
# parameters_tdxd_nhp.R
# Paramètres PK — T-DXd (DS-8201a) — SINGE CYNOMOLGUS (NHP)
#
# Méthode :
#   1. Allométrie depuis paramètres rat calibrés FDA (Table 6)
#      scale_CL = (BW_nhp/BW_rat)^0.75   [=  8.0  pour 4/0.25 kg]
#      scale_V  = (BW_nhp/BW_rat)^1.0    [= 16.0  pour 4/0.25 kg]
#   2. Calibration directe sur FDA BLA 761139 Table 7
#      "3-Month Intermittent IV Dose Toxicity Study in Cynomolgus Monkeys"
#      Doses Q3W : 3, 10, 30 mg/kg
#   3. Krel calibré sur [DXd]plasma / [ADC] ratio Table 7
#
# Différence clé vs rat :
#   → Singe exprime HER2 (Kd ≈ 7.46 ng/mL, étude ELISA FDA p.44)
#   → k_int_nhp estimé depuis la non-linéarité CL(dose) Table 7
#
# Sources :
#   FDA BLA 761139 Multi-Discipline Review (2019), Table 6 (rat) & Table 7 (NHP)
#   Yin et al. 2020 (PK humain, allométrie point de départ)
#   Vasalou et al. 2024 (k_int HER2, référence humain)
#
# Unités : temps [h] | concentrations [mg/L = µg/mL]
############################################################

tdxd_nhp <- list()

# ── Poids corporels ──────────────────────────────────────
BW_rat   <- 0.25   # kg  (étude rat FDA Table 6)
BW_nhp   <- 4.0    # kg  (singe cynomolgus adulte, protocole FDA : 3–10 ans)

# Facteurs allométriques NHP/rat
allo_CL  <- (BW_nhp / BW_rat)^0.75   # 16^0.75 = 8.0
allo_V   <- (BW_nhp / BW_rat)^1.0    # 16.0

# ── Propriétés ADC (identiques rat/humain) ───────────────
tdxd_nhp$DAR          <- 8
tdxd_nhp$MW_ADC       <- 148000  # g/mol
tdxd_nhp$MW_DXd       <- 718.8   # g/mol
tdxd_nhp$mass_frac_DXd <- 8 * 718.8 / 148000  # ≈ 0.03885

# ── Paramètres rat calibrés FDA Table 6 (point de départ) ─
CL_rat  <- 1.01e-4    # L/h
V1_rat  <- 0.01105    # L
V2_rat  <- 0.01803    # L
Q_rat   <- (0.174 / 24) * ((BW_rat / 70)^0.75)  # L/h

CL_DXd_rat <- 19.2 * ((BW_rat / 70)^0.75)   # L/h
V_DXd_rat  <- 29.41 * (BW_rat / 70)          # L
krel_rat   <- 1.19e-3   # h⁻¹

# ── Allométrie rat → NHP ─────────────────────────────────
CL_ADC_allom  <- CL_rat   * allo_CL
V1_ADC_allom  <- V1_rat   * allo_V
V2_ADC_allom  <- V2_rat   * allo_V
Q_ADC_allom   <- Q_rat    * allo_CL
CL_DXd_allom  <- CL_DXd_rat * allo_CL
V_DXd_allom   <- V_DXd_rat  * allo_V

# ── FDA BLA 761139 Table 7 — Données TK singe (Day 1) ────
fda_tk_nhp <- list(
  list(dose_mgkg  = 3,
       C0_ADC     = (101   + 95.3) / 2,
       AUC21d_ADC = (317   + 268)  / 2,
       t_half_d   = (3.95  + 3.85) / 2,
       C0_DXd_ng  = (0.242 + 0.248)/ 2),

  list(dose_mgkg  = 10,
       C0_ADC     = (295   + 339)  / 2,
       AUC21d_ADC = (1220  + 1080) / 2,
       t_half_d   = (5.56  + 5.13) / 2,
       C0_DXd_ng  = (0.656 + 1.02) / 2),

  list(dose_mgkg  = 30,
       C0_ADC     = (877   + 899)  / 2,
       AUC21d_ADC = (4090  + 3770) / 2,
       t_half_d   = (7.71  + 6.53) / 2,
       C0_DXd_ng  = (2.71  + 3.9)  / 2)
)

interval_h <- 21 * 24   # 504 h = Q3W

# ── Calibration V1, CL, V2 depuis Table 7 ────────────────
V1_v <- numeric(3); CL_v <- numeric(3); V2_v <- numeric(3)

for (i in seq_along(fda_tk_nhp)) {
  d      <- fda_tk_nhp[[i]]
  dose_mg <- d$dose_mgkg * BW_nhp

  V1 <- dose_mg / d$C0_ADC

  beta    <- log(2) / (d$t_half_d * 24)
  frac    <- 1 - exp(-beta * interval_h)
  AUC_inf <- d$AUC21d_ADC / frac

  CL <- dose_mg / (AUC_inf * 24)

  Vss <- CL / beta
  V2  <- max(Vss - V1, 0)

  V1_v[i] <- V1;  CL_v[i] <- CL;  V2_v[i] <- V2
}

V1_fda_nhp  <- sum(V1_v^2) / sum(V1_v)
CL_fda_nhp  <- sum(CL_v^2) / sum(CL_v)

t_half_tgt_d <- max(sapply(fda_tk_nhp, function(x) x$t_half_d)) * 1.25
beta_tgt     <- log(2) / (t_half_tgt_d * 24)
k10_nhp      <- CL_fda_nhp  / V1_fda_nhp
k12_nhp      <- Q_ADC_allom / V1_fda_nhp
k21_tgt      <- beta_tgt * (beta_tgt - k10_nhp - k12_nhp) / (beta_tgt - k10_nhp)
V2_fda_nhp   <- Q_ADC_allom / k21_tgt

k_int_nhp     <- 0

CL_DXd_nhp <- CL_DXd_allom
V_DXd_nhp  <- V_DXd_allom
V_ic_nhp   <- 0.003 * (BW_nhp / BW_rat)

krel_v <- numeric(3)
for (i in seq_along(fda_tk_nhp)) {
  d <- fda_tk_nhp[[i]]
  C_DXd_mgL  <- d$C0_DXd_ng * 1e-3
  C_ADC_mgL  <- d$C0_ADC
  krel_v[i]  <- C_DXd_mgL * CL_DXd_nhp /
                (tdxd_nhp$mass_frac_DXd * C_ADC_mgL * V1_v[i])
}
krel_fda_nhp <- mean(krel_v)

# ── Consolidation des paramètres T-DXd ───────────────────
tdxd_nhp$CL_ADC        <- CL_fda_nhp
tdxd_nhp$CL_lin        <- CL_fda_nhp
tdxd_nhp$Vmax_MM       <- 0.060
tdxd_nhp$Km_MM         <- 4.0
tdxd_nhp$V1_ADC        <- V1_fda_nhp
tdxd_nhp$V2_ADC        <- V2_fda_nhp
tdxd_nhp$Q_ADC         <- Q_ADC_allom
tdxd_nhp$k_int         <- k_int_nhp
tdxd_nhp$CL_DXd        <- CL_DXd_nhp
tdxd_nhp$V_DXd         <- V_DXd_nhp
tdxd_nhp$V_ic          <- V_ic_nhp
tdxd_nhp$k_rel_c1      <- krel_fda_nhp
tdxd_nhp$krel_power    <- -0.137
tdxd_nhp$krel_factor   <- 0.830
tdxd_nhp$interval_h    <- interval_h
tdxd_nhp$k_inD         <- 0.7
tdxd_nhp$k_effD        <- 0.7
tdxd_nhp$IC50_DXd_uM   <- 0.31
tdxd_nhp$k_dam         <- 0.017
tdxd_nhp$k_rep         <- 0.017
tdxd_nhp$mgL_to_uM_DXd <- 1000 / tdxd_nhp$MW_DXd

# ── Fonction d'administration T-DXd ─────────────────────
make_nhp_infusion <- function(dose_mgkg, BW_kg = 4.0,
                               Tinfu_h = 0.5, interval_h = NULL,
                               n_cycles = 1) {
  dose_mg  <- dose_mgkg * BW_kg
  rate_mgh <- dose_mg / Tinfu_h

  if (is.null(interval_h) || n_cycles == 1) {
    function(t) if (t >= 0 & t < Tinfu_h) rate_mgh else 0
  } else {
    t_starts <- seq(0, by = interval_h, length.out = n_cycles)
    function(t) {
      if (any(t >= t_starts & t < (t_starts + Tinfu_h))) rate_mgh else 0
    }
  }
}

# ════════════════════════════════════════════════════════
# CARBOPLATIN NHP PK — fit rxode2 (nca_analysis.R)
# Animaux : Animal_01 / Animal_02  |  Dose : 3 mg/kg IV bolus
# ════════════════════════════════════════════════════════
# Paramètres lus depuis pk2cmt_params.csv (généré par nca_analysis.R)
# Si le fichier n'existe pas, des valeurs de secours sont utilisées.

pk2cmt_file <- file.path(dirname(sys.frame(1)$ofile %||% "."), "pk2cmt_params.csv")
if (!file.exists(pk2cmt_file)) pk2cmt_file <- "pk2cmt_params.csv"

if (file.exists(pk2cmt_file)) {
  pk2cmt_res    <- read.csv(pk2cmt_file, stringsAsFactors = FALSE)
  mean_pars     <- pk2cmt_res[pk2cmt_res$Animal_Id == "Mean", ]
  CL_carbo_mLhkg <- mean_pars$CL_mL_h_kg
  V1_carbo_mLkg  <- mean_pars$V1_mL_kg
  Q_carbo_mLhkg  <- mean_pars$Q_mL_h_kg
  V2_carbo_mLkg  <- mean_pars$V2_mL_kg
  cat(sprintf("Carboplatin PK chargés depuis %s\n", pk2cmt_file))
  cat(sprintf("  CL=%.4f  V1=%.4f  Q=%.4f  V2=%.4f  mL/h/kg ou mL/kg\n",
              CL_carbo_mLhkg, V1_carbo_mLkg, Q_carbo_mLhkg, V2_carbo_mLkg))
} else {
  warning("pk2cmt_params.csv introuvable — valeurs de secours utilisées. Lancez nca_analysis.R d'abord.")
  CL_carbo_mLhkg <- 1.73
  V1_carbo_mLkg  <- 43.0
  Q_carbo_mLhkg  <- NA
  V2_carbo_mLkg  <- 132.0
}

carbo_nhp <- list(
  CL   = CL_carbo_mLhkg * BW_nhp / 1000,   # L/h   = 0.00692
  V1   = V1_carbo_mLkg  * BW_nhp / 1000,   # L     = 0.172
  Q    = if (is.na(Q_carbo_mLhkg)) NA_real_
         else Q_carbo_mLhkg * BW_nhp / 1000, # L/h  (placeholder)
  V2   = V2_carbo_mLkg  * BW_nhp / 1000,   # L     = 0.528
  t_half_h       = 70,                       # h     (régression 96-168h)
  fu0            = 1.0,   # fraction libre initiale (1 = platine libre mesuré)
  fu_inf         = 1.0,   # fraction libre terminale
  k_bind         = 0.0,   # liaison protéique (h⁻¹)
  k_dam          = 0.017, # formation adduits ADN (h⁻¹ µM⁻¹)
  k_rep          = 0.017, # réparation (h⁻¹)
  MW             = 371.25,
  mgL_to_uM      = 1000 / 371.25
)

# Fonction perfusion carboplatin NHP
make_nhp_carbo_infusion <- function(dose_mgkg, BW_kg = BW_nhp,
                                     Tinfu_h = 1, interval_h = NULL,
                                     n_cycles = 1) {
  dose_mg <- dose_mgkg * BW_kg
  rate    <- dose_mg / Tinfu_h
  if (is.null(interval_h) || n_cycles == 1) {
    function(t) if (t >= 0 & t < Tinfu_h) rate else 0
  } else {
    t_starts <- seq(0, by = interval_h, length.out = n_cycles)
    function(t) {
      if (any(t >= t_starts & t < (t_starts + Tinfu_h))) rate else 0
    }
  }
}

# État initial — inclut les compartiments carboplatin
tdxd_nhp_state0 <- c(
  C_ADC1   = 0,
  C_ADC2   = 0,
  C_DXd    = 0,
  C_DXd_ic = 0,
  C_carbo1 = 0,
  C_carbo2 = 0,
  Damage   = 0
)

cat("PK T-DXd NHP  : V1=", round(tdxd_nhp$V1_ADC,4), "L  CL=",
    formatC(tdxd_nhp$CL_ADC, format="e", digits=3), "L/h\n")
cat("PK Carboplatin: V1=", carbo_nhp$V1, "L  CL=", carbo_nhp$CL, "L/h",
    " t½=", carbo_nhp$t_half_h, "h\n")
