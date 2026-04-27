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
# Valeurs retranscrites depuis parameters_tdxd_rat.R
CL_rat  <- 1.01e-4    # L/h  (FDA rat calibré)
V1_rat  <- 0.01105    # L
V2_rat  <- 0.01803    # L
Q_rat   <- (0.174 / 24) * ((BW_rat / 70)^0.75)  # L/h  (allométrie depuis humain)

CL_DXd_rat <- 19.2 * ((BW_rat / 70)^0.75)   # L/h
V_DXd_rat  <- 29.41 * (BW_rat / 70)          # L
krel_rat   <- 1.19e-3   # h⁻¹ (calibré FDA rat Table 6)

# ── Allométrie rat → NHP ─────────────────────────────────
CL_ADC_allom  <- CL_rat   * allo_CL   # L/h
V1_ADC_allom  <- V1_rat   * allo_V    # L
V2_ADC_allom  <- V2_rat   * allo_V    # L
Q_ADC_allom   <- Q_rat    * allo_CL   # L/h
CL_DXd_allom  <- CL_DXd_rat * allo_CL  # L/h
V_DXd_allom   <- V_DXd_rat  * allo_V   # L

# ── FDA BLA 761139 Table 7 — Données TK singe (Day 1) ────
# Doses Q3W × 5 cycles (3 mois) : 3, 10, 30 mg/kg
# Valeurs = moyennes mâles/femelles, Day 1
# DS-8201a : C0 [µg/mL], AUC0-21d [µg.d/mL], T½ [jours]
# MAAA-1181a (DXd) : C0 [ng/mL]
fda_tk_nhp <- list(
  list(dose_mgkg  = 3,
       C0_ADC     = (101   + 95.3) / 2,   # µg/mL   M+F mean (M=101,   F=95.3)
       AUC21d_ADC = (317   + 268)  / 2,   # µg.d/mL M+F mean (M=317,   F=268)
       t_half_d   = (3.95  + 3.85) / 2,   # jours   M+F mean (M=3.95,  F=3.85)
       C0_DXd_ng  = (0.242 + 0.248)/ 2),  # ng/mL   M+F mean (M=0.242, F=0.248)

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
# Même méthode que parameters_tdxd_rat.R (FDA Table 6 → rat)
V1_v <- numeric(3); CL_v <- numeric(3); V2_v <- numeric(3)

for (i in seq_along(fda_tk_nhp)) {
  d      <- fda_tk_nhp[[i]]
  dose_mg <- d$dose_mgkg * BW_nhp    # mg total (BW = 4 kg)

  # V1 [L] = dose [mg] / C0 [mg/L]  (C0 µg/mL = mg/L)
  V1 <- dose_mg / d$C0_ADC

  # Correction AUC0-21d → AUC0-inf via T½
  beta    <- log(2) / (d$t_half_d * 24)        # h⁻¹
  frac    <- 1 - exp(-beta * interval_h)         # fraction éliminée en 21 j
  AUC_inf <- d$AUC21d_ADC / frac               # µg.d/mL

  # CL [L/h] = dose [mg] / (AUC_inf [µg.d/mL × 24 h/d])
  #          = dose [mg] / (AUC_inf [mg.d/L] × 24)
  CL <- dose_mg / (AUC_inf * 24)

  # V_ss [L] = CL / β ;   V2 = V_ss − V1
  Vss <- CL / beta
  V2  <- max(Vss - V1, 0)   # plancher à 0

  V1_v[i] <- V1;  CL_v[i] <- CL;  V2_v[i] <- V2
}

# WLS-optimal: minimise Σ(pred/obs − 1)²  →  w* = Σwi² / Σwi
V1_fda_nhp  <- sum(V1_v^2) / sum(V1_v)   # L    WLS-optimal C0
CL_fda_nhp  <- sum(CL_v^2) / sum(CL_v)   # L/h  WLS-optimal AUC

# V2 analytically → T½ linéaire cible = max(T½_obs) × 1.20
# Le TMDD ne peut qu'accélérer l'élimination (T½ < T½_linéaire).
# Donc T½_linéaire doit être > max T½ observé (7.12 j) pour que le TMDD
# puisse réduire T½ de 8.5 j → 7.12 j (30 mg/kg) et → 3.90 j (3 mg/kg).
t_half_tgt_d <- max(sapply(fda_tk_nhp, function(x) x$t_half_d)) * 1.20
beta_tgt     <- log(2) / (t_half_tgt_d * 24)    # h⁻¹
k10_nhp      <- CL_fda_nhp  / V1_fda_nhp
k12_nhp      <- Q_ADC_allom / V1_fda_nhp
# Solve 2-cpt beta quadratic for k21 given target beta:
#   k21 = beta*(beta − k10 − k12) / (beta − k10)
k21_tgt      <- beta_tgt * (beta_tgt - k10_nhp - k12_nhp) / (beta_tgt - k10_nhp)
V2_fda_nhp   <- Q_ADC_allom / k21_tgt            # L

# ── k_int NHP : k_int = 0 ────────────────────────────────
# Le T½ augmente avec la dose (3.9 → 5.3 → 7.1 j) = signature TMDD
# Un k_int linéaire ne capte pas la saturation HER2 → Option retenue :
#   k_int = 0,  V1/CL WLS-optimaux (moindres-carrés pondérés),
#   V2 analytique pour T½ = moyenne géométrique des T½ observés
k_int_nhp     <- 0                    # h⁻¹  (TMDD via MM, pas d'internalisation linéaire)

# ── DXd — CL_DXd et V_DXd allométriques depuis rat ──────
CL_DXd_nhp <- CL_DXd_allom   # L/h  (allom. rat → NHP, pas de calibration DXd)
V_DXd_nhp  <- V_DXd_allom    # L

# Volume intracellulaire (moelle osseuse) — proportionnel BW
# V_ic_rat = 0.003 L (0.25 kg)  →  V_ic_nhp = 0.003 × (4/0.25) = 0.048 L
V_ic_nhp   <- 0.003 * (BW_nhp / BW_rat)   # L

# ── Krel NHP calibré depuis Table 7 ──────────────────────
# Pseudo-SS : C_DXd = krel × mass_frac × C_ADC × V1 / CL_DXd
# → krel = C_DXd [mg/L] × CL_DXd / (mass_frac × C_ADC [mg/L] × V1 [L])
# NOTE : C0_DXd en ng/mL → × 1e-3 pour mg/L
krel_v <- numeric(3)
for (i in seq_along(fda_tk_nhp)) {
  d <- fda_tk_nhp[[i]]
  C_DXd_mgL  <- d$C0_DXd_ng * 1e-3        # ng/mL → mg/L
  C_ADC_mgL  <- d$C0_ADC                  # µg/mL = mg/L
  krel_v[i]  <- C_DXd_mgL * CL_DXd_nhp /
                (tdxd_nhp$mass_frac_DXd * C_ADC_mgL * V1_v[i])
}
krel_fda_nhp <- mean(krel_v)   # h⁻¹

# ── Consolidation des paramètres finaux ───────────────────
tdxd_nhp$CL_ADC        <- CL_fda_nhp      # WLS-optimal (référence linéaire)
tdxd_nhp$CL_lin        <- CL_fda_nhp      # L/h  (départ optim TMDD = CL_wls)
tdxd_nhp$Vmax_MM       <- 0.060           # mg/L/h (départ grille analytique 1-cpt)
tdxd_nhp$Km_MM         <- 4.0             # mg/L   (départ grille analytique 1-cpt)
tdxd_nhp$V1_ADC        <- V1_fda_nhp
tdxd_nhp$V2_ADC        <- V2_fda_nhp
tdxd_nhp$Q_ADC         <- Q_ADC_allom     # allométrie (pas de calibration Q)
tdxd_nhp$k_int         <- k_int_nhp
tdxd_nhp$CL_DXd        <- CL_DXd_nhp
tdxd_nhp$V_DXd         <- V_DXd_nhp
tdxd_nhp$V_ic          <- V_ic_nhp
tdxd_nhp$k_rel_c1      <- krel_fda_nhp
tdxd_nhp$krel_power    <- -0.137    # Yin 2020 (identique rat/humain)
tdxd_nhp$krel_factor   <- 0.830
tdxd_nhp$interval_h    <- interval_h
tdxd_nhp$k_inD         <- 0.7       # h⁻¹  (Vasalou 2024)
tdxd_nhp$k_effD        <- 0.7
tdxd_nhp$IC50_DXd_uM   <- 0.31      # µM
tdxd_nhp$k_dam         <- 0.017
tdxd_nhp$k_rep         <- 0.017
tdxd_nhp$mgL_to_uM_DXd <- 1000 / tdxd_nhp$MW_DXd

# ── Cibles de validation FDA Table 7 ────────────────────
fda_nhp_targets <- fda_tk_nhp

# ── Fonction d'administration IV (identique rat) ─────────
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

# ── État initial PK ──────────────────────────────────────
tdxd_nhp_state0 <- c(
  C_ADC1   = 0,
  C_ADC2   = 0,
  C_DXd    = 0,
  C_DXd_ic = 0,
  Damage   = 0
)

# ════════════════════════════════════════════════════════
# RÉSUMÉ
# ════════════════════════════════════════════════════════
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  PK T-DXd — SINGE CYNOMOLGUS (calibré FDA BLA Table 7)  ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  BW_nhp = %.1f kg   |  allo_CL = (4/0.25)^0.75 = %.1f   allo_V = %.0f\n\n",
            BW_nhp, allo_CL, allo_V))

cat("── ADC (2-cpt) : allométrie rat → NHP ──────────────────────\n")
cat(sprintf("  Allométrique  : V1=%.4f L  V2=%.4f L  CL=%.4e L/h\n",
            V1_ADC_allom, V2_ADC_allom, CL_ADC_allom))
cat(sprintf("  Calibré FDA   : V1=%.5f L  V2=%.5f L  CL_wls=%.4e L/h\n",
            tdxd_nhp$V1_ADC, tdxd_nhp$V2_ADC, tdxd_nhp$CL_ADC))
cat(sprintf("  T½ linéaire cible : %.2f j  (max obs × 1.20)\n",
            t_half_tgt_d))
cat(sprintf("  Ratio FDA/allo: V1×%.2f            CL×%.2f\n\n",
            tdxd_nhp$V1_ADC / V1_ADC_allom,
            tdxd_nhp$CL_ADC / CL_ADC_allom))

cat("── TMDD (Michaelis-Menten) : calibration via optim() ──────────\n")
cat(sprintf("  k_int = 0  |  CL_lin départ=CL_wls, Vmax_MM=0.060, Km_MM=4.0 (grille 1-cpt)\n"))
cat(sprintf("  Calibration finale (C0+AUC+T½ × 3 doses) → run_pkpd_tdxd_nhp.R\n"))
cat(sprintf("  V2 analytique → T½ linéaire cible = %.2f j\n\n", t_half_tgt_d))

cat("── DXd (1-cpt) : allométrie rat → NHP ─────────────────────\n")
cat(sprintf("  CL_DXd = %.4f L/h   V_DXd = %.4f L   V_ic = %.4f L\n",
            tdxd_nhp$CL_DXd, tdxd_nhp$V_DXd, tdxd_nhp$V_ic))
cat(sprintf("  Krel   = %.4e h⁻¹  (calibré C0_DXd/C0_ADC Table 7 M+F)\n\n",
            tdxd_nhp$k_rel_c1))

cat("── Validation NCA Table 7 (2-cpt linéaire, avant TMDD) ────────\n")
cat(sprintf("  %-8s  %-7s %-7s %-7s │ %-8s %-8s %-7s │ %-6s\n",
            "Dose","C0_obs","C0_sim","ratio","AUC_obs","AUC_sim","ratio","T½_sim"))
cat(sprintf("  %-8s  %-7s %-7s %-7s │ %-8s %-8s %-7s │ %-6s\n",
            "mg/kg","µg/mL","µg/mL","","µg.d/mL","µg.d/mL","","jours"))
cat(sprintf("  %s\n", paste(rep("─", 80), collapse="")))

for (i in seq_along(fda_tk_nhp)) {
  d       <- fda_tk_nhp[[i]]
  dose_mg <- d$dose_mgkg * BW_nhp

  # 2-cpt analytical prediction
  V1 <- tdxd_nhp$V1_ADC; CL <- tdxd_nhp$CL_ADC
  V2 <- tdxd_nhp$V2_ADC; Q  <- tdxd_nhp$Q_ADC
  k10 <- CL/V1; k12 <- Q/V1; k21 <- Q/V2
  s   <- k10+k12+k21
  b   <- (s - sqrt(s^2 - 4*k10*k21))/2
  a   <- s - b
  A   <- (a-k21)/(a-b) * dose_mg/V1
  B   <- (k21-b)/(a-b) * dose_mg/V1
  T   <- 21*24
  C0_sim  <- dose_mg / V1
  AUC_sim <- (A/a*(1-exp(-a*T)) + B/b*(1-exp(-b*T))) / 24
  t12_sim <- log(2)/b/24

  cat(sprintf("  %-8s  %-7.1f %-7.1f %-7.3f │ %-8.1f %-8.1f %-7.3f │ %-6.2f\n",
              paste0(d$dose_mgkg, " mg/kg"),
              d$C0_ADC, C0_sim, C0_sim/d$C0_ADC,
              d$AUC21d_ADC, AUC_sim, AUC_sim/d$AUC21d_ADC,
              t12_sim))
}
cat(sprintf("  (T½ linéaire — modèle TMDD calibré dans run_pkpd_tdxd_nhp.R)\n\n"))

cat("── Krel par dose ───────────────────────────────────────────\n")
for (i in seq_along(fda_tk_nhp)) {
  d <- fda_tk_nhp[[i]]
  cat(sprintf("  %2d mg/kg : C0_DXd=%.3f ng/mL → Krel=%.4e h⁻¹\n",
              d$dose_mgkg, d$C0_DXd_ng, krel_v[i]))
}
cat(sprintf("  Moyenne  : Krel = %.4e h⁻¹  (rat=1.19e-3, Krel_nhp/rat=%.2f)\n\n",
            tdxd_nhp$k_rel_c1, tdxd_nhp$k_rel_c1 / 1.19e-3))
