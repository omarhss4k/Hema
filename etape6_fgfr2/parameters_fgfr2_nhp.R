############################################################
# parameters_fgfr2_nhp.R
# Paramètres PK — Inhibiteur FGFR2 — SINGE CYNOMOLGUS (NHP)
#
# Méthode :
#   1. Allométrie depuis paramètres rat calibrés (Table 6 données précliniques)
#      scale_CL = (BW_nhp/BW_rat)^0.75   [=  8.0  pour 4/0.25 kg]
#      scale_V  = (BW_nhp/BW_rat)^1.0    [= 16.0  pour 4/0.25 kg]
#   2. Calibration directe sur données précliniques NHP
#      "Étude toxicité IV — singe cynomolgus"
#      Doses Q3W : 3, 10, 30 mg/kg
#
# Différence clé vs rat :
#   → Singe exprime FGFR2 (non-linéarité CL dose-dépendante)
#   → k_int estimé depuis la non-linéarité CL(dose) observée
#
# Sources :
#   Données précliniques internes — singe cynomolgus (pk2cmt_params.csv)
#   Allométrie de référence : Yin et al. 2020
#
# Unités : temps [h] | concentrations [mg/L = µg/mL]
############################################################

fgfr2_nhp <- list()

# ── Poids corporels ──────────────────────────────────────
BW_rat   <- 0.25   # kg  (étude rat FDA Table 6)
BW_nhp   <- 4.0    # kg  (singe cynomolgus adulte, protocole FDA : 3–10 ans)

# Facteurs allométriques NHP/rat
allo_CL  <- (BW_nhp / BW_rat)^0.75   # 16^0.75 = 8.0
allo_V   <- (BW_nhp / BW_rat)^1.0    # 16.0

# ── Propriétés moléculaires (placeholder — à adapter pour FGFR2) ────
fgfr2_nhp$DAR          <- 8
fgfr2_nhp$MW_ADC       <- 148000  # g/mol
fgfr2_nhp$MW_DXd       <- 718.8   # g/mol
fgfr2_nhp$mass_frac_DXd <- 8 * 718.8 / 148000  # ≈ 0.03885

# ── Paramètres rat calibrés FDA Table 6 (point de départ) ─
CL_rat  <- 1.01e-4    # L/h
V1_rat  <- 0.01105    # L
V2_rat  <- 0.01803    # L
Q_rat   <- (0.174 / 24) * ((BW_rat / 70)^0.75)  # L/h

# ── Allométrie rat → NHP ─────────────────────────────────
CL_ADC_allom  <- CL_rat   * allo_CL
V1_ADC_allom  <- V1_rat   * allo_V
V2_ADC_allom  <- V2_rat   * allo_V
Q_ADC_allom   <- Q_rat    * allo_CL
# Valeurs DXd rat (Yin 2020, allométrie depuis humain)
CL_DXd_rat   <- 19.2 * (BW_rat / 70)^0.75   # L/h
V_DXd_rat    <- (17 * 1.73) * (BW_rat / 70)^1.0  # L
CL_DXd_allom  <- CL_DXd_rat * allo_CL
V_DXd_allom   <- V_DXd_rat  * allo_V

# ── Cibles de calibration TK — données réelles NHP (pk2cmt_params.csv) ───
# Mapping : Animal_01 = 3 mg/kg | Animal_03 = 10 mg/kg | Animal_04 = 30 mg/kg
# Fallback : FDA BLA 761139 Table 7 si le fichier est absent

pk2cmt_tk_file <- "pk2cmt_params.csv"

if (file.exists(pk2cmt_tk_file)) {
  pk2cmt_tk <- read.csv(pk2cmt_tk_file, stringsAsFactors = FALSE)
  pk2cmt_tk <- pk2cmt_tk[!is.na(pk2cmt_tk$dose_mg_kg), ]

  # Facteur de correction dose : 1.297 (correction dose réelle vs dose nominale)
  # Doses corrigées : 3→4 | 10→13 | 20→26 | 30→39 mg/kg
  # 2 animaux par dose — moyennés
  dose_factor <- 1.297
  doses_nhp   <- c(4, 13, 26, 39)
  dose_animal_map <- list(
    "4"  = c("Animal_01", "Animal_02"),
    "13" = c("Animal_03", "Animal_04"),
    "26" = c("Animal_05", "Animal_06"),
    "39" = c("Animal_07", "Animal_08")
  )

  fda_tk_nhp <- lapply(doses_nhp, function(d) {
    animals <- dose_animal_map[[as.character(d)]]
    rows    <- pk2cmt_tk[pk2cmt_tk$Animal_Id %in% animals, ]
    if (nrow(rows) == 0) return(NULL)   # dose sans données — ignorée
    beta_h  <- mean(log(2) / rows$t12_beta_h)
    C0      <- mean(d * 1000 / rows$V1_mL_kg)
    AUCinf  <- mean(d * 1000 / (rows$CL_mL_h_kg * 24))
    list(dose_mgkg  = d,
         C0_ADC     = C0,
         AUC21d_ADC = AUCinf * (1 - exp(-beta_h * 21 * 24)),
         t_half_d   = mean(rows$t12_beta_h) / 24)
  })
  fda_tk_nhp <- Filter(Negate(is.null), fda_tk_nhp)   # retirer les NULL
  doses_avec_data <- sapply(fda_tk_nhp, function(x) x$dose_mgkg)
  cat(sprintf("Cibles TK : données réelles NHP — doses disponibles : %s mg/kg\n",
              paste(doses_avec_data, collapse = ", ")))
} else {
  warning("pk2cmt_params.csv introuvable — fallback valeurs corrigées ×1.297")
  # Valeurs originales × 1.297 (dose et concentrations)
  fda_tk_nhp <- list(
    list(dose_mgkg=4,  C0_ADC=(101+95.3)/2  *1.297, AUC21d_ADC=(317+268)/2   *1.297, t_half_d=(3.95+3.85)/2),
    list(dose_mgkg=13, C0_ADC=(295+339)/2   *1.297, AUC21d_ADC=(1220+1080)/2 *1.297, t_half_d=(5.56+5.13)/2),
    list(dose_mgkg=26, C0_ADC=(586+638)/2   *1.297, AUC21d_ADC=(2655+2425)/2 *1.297, t_half_d=(6.64+5.83)/2),
    list(dose_mgkg=39, C0_ADC=(877+899)/2   *1.297, AUC21d_ADC=(4090+3770)/2 *1.297, t_half_d=(7.71+6.53)/2)
  )
}

interval_h <- 21 * 24   # 504 h = Q3W

# ── Calibration V1, CL, V2 ───────────────────────────────
n_doses <- length(fda_tk_nhp)
V1_v <- numeric(n_doses); CL_v <- numeric(n_doses); V2_v <- numeric(n_doses)

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
denom_k21    <- beta_tgt - k10_nhp
# Si beta_tgt <= k10_nhp la formule analytique diverge : fallback allométrique
if (is.finite(denom_k21) && abs(denom_k21) > 1e-9) {
  k21_tgt    <- beta_tgt * (beta_tgt - k10_nhp - k12_nhp) / denom_k21
  V2_fda_nhp <- if (is.finite(k21_tgt) && k21_tgt > 0)
                  Q_ADC_allom / k21_tgt
                else V2_ADC_allom
} else {
  V2_fda_nhp <- V2_ADC_allom   # fallback allométrique
}

# ── Consolidation des paramètres FGFR2 ────────────────────
# PK 2-cmt linéaire — pas de TMDD :
#   FGFR2 n'est pas exprimé sur les HSC/progéniteurs hématopoïétiques.
#   L'hématotoxicité est due au payload (effet bystander), pas à
#   l'internalisation ADC via FGFR2.  CL calibré sur données NCA NHP.
fgfr2_nhp$CL_ADC        <- CL_fda_nhp
fgfr2_nhp$V1_ADC        <- V1_fda_nhp
fgfr2_nhp$V2_ADC        <- V2_fda_nhp
fgfr2_nhp$Q_ADC         <- Q_ADC_allom
fgfr2_nhp$interval_h    <- interval_h
fgfr2_nhp$k_dam          <- 0.075
fgfr2_nhp$k_rep          <- 0.017
fgfr2_nhp$mgL_to_uM_ADC  <- 1000 / fgfr2_nhp$MW_ADC   # mg/L → µM  (MW=148000 g/mol)

# ── Fonction d'administration FGFR2 ─────────────────────
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

pk2cmt_file <- "pk2cmt_params.csv"

if (file.exists(pk2cmt_file)) {
  pk2cmt_res    <- read.csv(pk2cmt_file, stringsAsFactors = FALSE)
  anim_pars      <- pk2cmt_res[pk2cmt_res$Animal_Id == "Animal_01", ]
  CL_carbo_mLhkg <- anim_pars$CL_mL_h_kg
  V1_carbo_mLkg  <- anim_pars$V1_mL_kg
  Q_carbo_mLhkg  <- anim_pars$Q_mL_h_kg
  V2_carbo_mLkg  <- anim_pars$V2_mL_kg
  cat(sprintf("Carboplatin PK chargés depuis %s (Animal_01)\n", pk2cmt_file))
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
fgfr2_nhp_state0 <- c(
  C_ADC1   = 0,
  C_ADC2   = 0,
  C_DXd    = 0,
  C_DXd_ic = 0,
  C_carbo1 = 0,
  C_carbo2 = 0,
  Damage   = 0
)

cat("PK FGFR2 NHP  : V1=", round(fgfr2_nhp$V1_ADC,4), "L  CL=",
    formatC(fgfr2_nhp$CL_ADC, format="e", digits=3), "L/h\n")
cat("PK Carboplatin: V1=", carbo_nhp$V1, "L  CL=", carbo_nhp$CL, "L/h",
    " t½=", carbo_nhp$t_half_h, "h\n")
