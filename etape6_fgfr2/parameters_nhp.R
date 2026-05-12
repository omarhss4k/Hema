############################################################
# parameters_nhp.R
# Paramètres physiologiques — NHP (macaque cynomolgus/rhésus)
# PK : issus du fit rxode2 2-compartiments sur données NHP
#       (nca_analysis.R — Animal_01 / Animal_02)
############################################################

init_pars <- list()

# ══════════════════════════════════════════════════════════
# POIDS CORPOREL NHP
# ══════════════════════════════════════════════════════════
# Ajustez BW_nhp avec le poids réel de chaque animal.
# Les paramètres PK absolus (L, L/h) en dépendent directement.

BW_nhp <- 5.0   # kg  ← REMPLACEZ par le poids réel du NHP

# ══════════════════════════════════════════════════════════
# PK CARBOPLATIN — modèle 2 compartiments IV bolus
# ══════════════════════════════════════════════════════════
# Source : fit rxode2 (nca_analysis.R) sur Animal_01 / Animal_02
# Paramètres lus depuis pk2cmt_params.csv (généré par nca_analysis.R).
# Si le fichier est absent, des valeurs de secours sont utilisées.
#
# IMPORTANT : ces concentrations mesurent le PLATINE TOTAL.
# Si le modèle PD utilise le platine libre, ajustez fu0/fu_inf/k_bind.

pk2cmt_file_nhp <- "pk2cmt_params.csv"

if (file.exists(pk2cmt_file_nhp)) {
  pk2cmt_res  <- read.csv(pk2cmt_file_nhp, stringsAsFactors = FALSE)
  anim_pars   <- pk2cmt_res[pk2cmt_res$Animal_Id == "Animal_01", ]
  CL_mL_h_kg <- anim_pars$CL_mL_h_kg
  V1_mL_kg   <- anim_pars$V1_mL_kg
  Q_mL_h_kg  <- anim_pars$Q_mL_h_kg
  V2_mL_kg   <- anim_pars$V2_mL_kg
} else {
  warning("pk2cmt_params.csv introuvable — valeurs de secours. Lancez nca_analysis.R d'abord.")
  CL_mL_h_kg <- 1.73
  V1_mL_kg   <- 43.0
  Q_mL_h_kg  <- NA
  V2_mL_kg   <- 132.0
}

init_pars$CL <- CL_mL_h_kg * BW_nhp / 1000   # L/h
init_pars$V1 <- V1_mL_kg   * BW_nhp / 1000   # L
init_pars$Q  <- Q_mL_h_kg  * BW_nhp / 1000   # L/h
init_pars$V2 <- V2_mL_kg   * BW_nhp / 1000   # L

cat(sprintf("PK NHP (BW = %.1f kg) — source : %s\n", BW_nhp,
            ifelse(file.exists(pk2cmt_file_nhp), "pk2cmt_params.csv (rxode2 fit)", "valeurs de secours")))
cat(sprintf("  CL = %.4f L/h  |  V1 = %.4f L\n", init_pars$CL, init_pars$V1))
cat(sprintf("  Q  = %.4f L/h  |  V2 = %.4f L\n", init_pars$Q,  init_pars$V2))
cat(sprintf("  t1/2 terminal ≈ %.1f h\n\n", log(2) / (init_pars$CL / init_pars$V2)))

# ── Liaison protéique ──────────────────────────────────────
# Si les concentrations mesurées = platine total :
#   fu0    = fraction libre initiale  (ex. 0.3–0.5 pour NHP)
#   fu_inf = fraction libre terminale (ex. 0.05–0.1)
#   k_bind = constante de liaison     (ex. 0.01–0.05 /h)
# Si concentrations = platine libre → laisser fu = 1, k_bind = 0
init_pars$fu0    <- 1.0    # ← à ajuster si platine total mesuré
init_pars$fu_inf <- 1.0
init_pars$k_bind <- 0.0

# ══════════════════════════════════════════════════════════
# DAMAGE (ADN) — identique rat / humain (Fornari 2019 Table 1)
# ══════════════════════════════════════════════════════════
init_pars$k_dam <- 0.075   # formation d'adduits (h⁻¹ µM⁻¹) — augmenté ×4.4 pour compenser AUC réel vs FDA
init_pars$k_rep <- 0.017   # réparation ADN (h⁻¹)

# ══════════════════════════════════════════════════════════
# BASELINES NHP (10⁹ cellules/L sauf MPP/CMP/MEP)
# Références : Davies et al. 2017, Fabian et al. 2018 (macaque)
# ══════════════════════════════════════════════════════════
init_pars$MPP0  <-    1.5    # ← à ajuster (interpolation rat/humain)
init_pars$CMP0  <-   18.0   # ← à ajuster
init_pars$MEP0  <-   13.0   # ← à ajuster
init_pars$Neut0 <-    1.98   # 10⁹/L  — moyenne pré-dose (t=-3) 8 animaux
init_pars$Mono0 <-    0.4    # 10⁹/L
init_pars$Ret0  <-   82.5   # 10⁹/L  — moyenne pré-dose (t=-3) 8 animaux
init_pars$RBC0  <- 5806.0   # 10⁹/L  — moyenne pré-dose (t=-3) 8 animaux ×1000
init_pars$Plt0  <-  436.0   # 10⁹/L  — moyenne pré-dose (t=-3) 8 animaux

# ══════════════════════════════════════════════════════════
# MTTs NHP (en heures) — proches humain, légèrement réduits
# ══════════════════════════════════════════════════════════
init_pars$MTT_Neut <- 168.0   # h  (7 jours, humain = 210h)
init_pars$MTT_Mono <- 110.0   # h
init_pars$MTT_Ret  <-  66.0   # h  (identique humain)
init_pars$MTT_Plt  <- 168.0   # h

# ══════════════════════════════════════════════════════════
# TAUX DE CIRCULATION (h⁻¹)
# ══════════════════════════════════════════════════════════
init_pars$k_circ_Neut <- 0.100
init_pars$k_circ_Mono <- 0.040
init_pars$k_circ_Plt  <- 0.0052
init_pars$k_circ_RBC  <- 0.00037

# ══════════════════════════════════════════════════════════
# EFFETS DU MÉDICAMENT (slopes) — à calibrer sur données NHP
# Valeur initiale : humain (Fornari 2019 Table 2)
# ══════════════════════════════════════════════════════════
init_pars$Slope_MPP <- 0.79   # ← à recalibrer sur données NHP
init_pars$Slope_CMP <- 0.57
init_pars$Slope_MEP <- 0.66
init_pars$delta_Ret <- 2.8
init_pars$delta_Plt <- 0.54

# ── IC50 NHP (à renseigner si données disponibles) ────────
init_pars$IC50_CMP_nhp <- NA   # µM
init_pars$IC50_MEP_nhp <- NA   # µM

# ── Feedback powers (Table 1, identiques toutes espèces) ──
init_pars$gamma_stem      <- 0.07
init_pars$gamma_mat_CMP   <- 0.60
init_pars$gamma_mat_MEP   <- 0.30
init_pars$gamma_prolTrans <- 0.70

# ══════════════════════════════════════════════════════════
# MW carboplatin + conversion mg/L → µM
# ══════════════════════════════════════════════════════════
init_pars$MW_carboplatin <- 371.25
init_pars$mgL_to_uM      <- 1000 / init_pars$MW_carboplatin

# ══════════════════════════════════════════════════════════
# ÉTAT INITIAL (complété par parameters_FORNARI_CORRECT.R)
# ══════════════════════════════════════════════════════════
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

# ══════════════════════════════════════════════════════════
# FONCTION PERFUSION
# ══════════════════════════════════════════════════════════
# dose_mg = dose en mg ABSOLUE (= dose_mg_kg × BW_nhp)
make_repeated_infusion <- function(dose_mg, Tinfu_h = 1, interval_h, n_cycles) {
  rate     <- dose_mg / Tinfu_h
  t_starts <- seq(0, by = interval_h, length.out = n_cycles)
  function(t) {
    if (any(t >= t_starts & t < (t_starts + Tinfu_h))) rate else 0
  }
}
