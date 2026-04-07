############################################################
# run_tdxd_rat.R
# Simulations PK/Damage T-DXd — RAT
#
# Chaîne modélisée :
#   ADC sérum (2 cpt) → DXd plasma → DXd intracell. → Damage ADN
#
# Scénarios :
#   1. Dose unique 5 mg/kg IV (perfusion 30 min)
#   2. Dose unique 10 mg/kg IV (perfusion 30 min)
#   3. Q3W × 3 cycles, 5 mg/kg
############################################################
source("pk_model_tdxd.R")
source("parameters_tdxd_rat.R")

if (!dir.exists("results_TDXD")) dir.create("results_TDXD")

# ──────────────────────────────────────────────────────────
# 1. Dose unique — 5 mg/kg
# ──────────────────────────────────────────────────────────
cat("=== Scénario 1 : 5 mg/kg dose unique (IV 30 min) ===\n")

pars_s1          <- tdxd_pars
pars_s1$rate_fun <- make_tdxd_infusion(dose_mgkg = 5, BW_kg = 0.25,
                                       Tinfu_h = 0.5, n_cycles = 1)
times_s1 <- seq(0, 21 * 24, by = 0.5)   # 21 jours, pas de 30 min

sim_s1 <- simulate_pk_tdxd(times_s1, pars_s1, tdxd_state0)

plot_pk_tdxd(
  sim       = sim_s1,
  pars      = tdxd_pars,
  titre     = "T-DXd 5 mg/kg dose unique — Rat",
  dose_times_h = 0,
  file      = "results_TDXD/PK_5mgkg_single.pdf"
)

# ──────────────────────────────────────────────────────────
# 2. Dose unique — 10 mg/kg
# ──────────────────────────────────────────────────────────
cat("\n=== Scénario 2 : 10 mg/kg dose unique (IV 30 min) ===\n")

pars_s2          <- tdxd_pars
pars_s2$rate_fun <- make_tdxd_infusion(dose_mgkg = 10, BW_kg = 0.25,
                                       Tinfu_h = 0.5, n_cycles = 1)
times_s2 <- seq(0, 21 * 24, by = 0.5)

sim_s2 <- simulate_pk_tdxd(times_s2, pars_s2, tdxd_state0)

plot_pk_tdxd(
  sim          = sim_s2,
  pars         = tdxd_pars,
  titre        = "T-DXd 10 mg/kg dose unique — Rat",
  dose_times_h = 0,
  file         = "results_TDXD/PK_10mgkg_single.pdf"
)

# ──────────────────────────────────────────────────────────
# 3. Q3W × 3 cycles — 5 mg/kg
# ──────────────────────────────────────────────────────────
cat("\n=== Scénario 3 : 5 mg/kg Q3W × 3 cycles ===\n")

interval_Q3W  <- 21 * 24   # 504 h
pars_s3          <- tdxd_pars
pars_s3$rate_fun <- make_tdxd_infusion(dose_mgkg = 5, BW_kg = 0.25,
                                       Tinfu_h = 0.5,
                                       interval_h = interval_Q3W,
                                       n_cycles = 3)
times_s3     <- seq(0, 63 * 24, by = 1)   # 63 jours
dose_days_s3 <- c(0, 21, 42)              # jours des doses

sim_s3 <- simulate_pk_tdxd(times_s3, pars_s3, tdxd_state0)

plot_pk_tdxd(
  sim          = sim_s3,
  pars         = tdxd_pars,
  titre        = "T-DXd 5 mg/kg Q3W × 3 cycles — Rat",
  dose_times_h = dose_days_s3 * 24,
  file         = "results_TDXD/PK_5mgkg_Q3Wx3.pdf"
)

# ──────────────────────────────────────────────────────────
# 4. Validation PK humaine — FDA BLA 761139
# Simule 5.4 mg/kg Q3W chez humain (70 kg, params non scalés)
# Vérifie Cmax ADC et Cmax DXd vs cibles FDA
# ──────────────────────────────────────────────────────────
cat("\n=== Validation PK humaine (FDA BLA 761139) ===\n")
cat("    Dose: 5.4 mg/kg × 70 kg = 378 mg, perfusion 1.5h, Q3W × 1 cycle\n")

pars_hu          <- tdxd_pars_human
pars_hu$rate_fun <- make_tdxd_infusion(
  dose_mgkg = tdxd_fda_targets$dose_mgkg,
  BW_kg     = tdxd_fda_targets$BW_kg,
  Tinfu_h   = tdxd_fda_targets$Tinfu_h,
  n_cycles  = 1
)

times_hu  <- seq(0, 21 * 24, by = 0.25)   # 21 jours, pas 15 min

sim_hu <- simulate_pk_tdxd(times_hu, pars_hu, tdxd_state0)

# Extraction des métriques
Cmax_ADC_sim   <- max(sim_hu$C_ADC1,   na.rm = TRUE)   # mg/L
Cmax_DXd_sim   <- max(sim_hu$C_DXd,    na.rm = TRUE)   # mg/L
Cmax_DXd_sim_ng <- Cmax_DXd_sim * 1000                 # ng/mL
AUC_ADC_trap   <- sum(diff(sim_hu$time_h) *
                       (head(sim_hu$C_ADC1, -1) + tail(sim_hu$C_ADC1, -1)) / 2)

cat(sprintf("\n  Résultats simulation (humain, 5.4 mg/kg Q3W) :\n"))
cat(sprintf("  ADC  Cmax  simulé  = %6.1f mg/L    |  cible FDA = %6.1f  |  écart = %+.1f%%\n",
            Cmax_ADC_sim,
            tdxd_fda_targets$Cmax_ADC_mgL,
            100 * (Cmax_ADC_sim - tdxd_fda_targets$Cmax_ADC_mgL) /
              tdxd_fda_targets$Cmax_ADC_mgL))
cat(sprintf("  DXd  Cmax  simulé  = %6.4f ng/mL   |  cible FDA = %6.1f  |  écart = %+.1f%%\n",
            Cmax_DXd_sim_ng,
            tdxd_fda_targets$Cmax_DXd_ngmL,
            100 * (Cmax_DXd_sim_ng - tdxd_fda_targets$Cmax_DXd_ngmL) /
              tdxd_fda_targets$Cmax_DXd_ngmL))
cat(sprintf("  ADC  AUC0-21d sim. = %6.0f mg·h/L  |  cible ~20000 (Dose/CL)\n",
            AUC_ADC_trap))

# Verdict
pass_Cmax_ADC <- abs(Cmax_ADC_sim - tdxd_fda_targets$Cmax_ADC_mgL) /
                 tdxd_fda_targets$Cmax_ADC_mgL <= tdxd_fda_targets$tol_Cmax
pass_Cmax_DXd <- abs(Cmax_DXd_sim_ng - tdxd_fda_targets$Cmax_DXd_ngmL) /
                 tdxd_fda_targets$Cmax_DXd_ngmL <= tdxd_fda_targets$tol_Cmax

cat(sprintf("\n  [%s] Cmax ADC  dans ±%.0f%% de la cible FDA\n",
            ifelse(pass_Cmax_ADC, "OK", "!!"), tdxd_fda_targets$tol_Cmax * 100))
cat(sprintf("  [%s] Cmax DXd  dans ±%.0f%% de la cible FDA\n",
            ifelse(pass_Cmax_DXd, "OK", "!!"), tdxd_fda_targets$tol_Cmax * 100))

if (!pass_Cmax_ADC || !pass_Cmax_DXd)
  cat("  >>> Vérifier paramètres PK humains (CL_ADC, V1_ADC, Krel)\n")

# ──────────────────────────────────────────────────────────
# Tableau récapitulatif
# ──────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════\n")
cat("  PARAMÈTRES PK UTILISÉS (RAT)\n")
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("  CL_ADC  = %.4e L/h\n", tdxd_pars$CL_ADC))
cat(sprintf("  V1_ADC  = %.5f L\n",   tdxd_pars$V1_ADC))
cat(sprintf("  V2_ADC  = %.5f L\n",   tdxd_pars$V2_ADC))
cat(sprintf("  Q_ADC   = %.4e L/h\n", tdxd_pars$Q_ADC))
cat(sprintf("  k_int   = %.5f h⁻¹ (t½=%.1fh)\n",
            tdxd_pars$k_int, log(2)/tdxd_pars$k_int))
cat(sprintf("  CL_DXd  = %.5f L/h\n", tdxd_pars$CL_DXd))
cat(sprintf("  V_DXd   = %.5f L\n",   tdxd_pars$V_DXd))
cat(sprintf("  IC50    = %.2f µM\n",  tdxd_pars$IC50_DXd_uM))
cat(sprintf("  V_ic    = %.4f L   (ratio ic/plasma = %.0fx)\n",
            tdxd_pars$V_ic, tdxd_pars$V_DXd / tdxd_pars$V_ic))
cat(sprintf("  k_dam   = %.4f h⁻¹  k_rep = %.4f h⁻¹\n",
            tdxd_pars$k_dam, tdxd_pars$k_rep))
cat("═══════════════════════════════════════════════════════════\n")
cat("  Fichiers : results_TDXD/\n")
cat("    -> PK_5mgkg_single.pdf   (6 panels : ADC + DXd + DXd_ic + Ratio + Damage + Emax)\n")
cat("    -> PK_10mgkg_single.pdf\n")
cat("    -> PK_5mgkg_Q3Wx3.pdf\n")
cat("  Validation FDA BLA 761139 :\n")
cat("    -> Cmax ADC et DXd simulés vs cibles GeomMean FDA (5.4 mg/kg Q3W humain)\n")
cat("    -> Cibles : ADC Cmax=122 mg/L | DXd Cmax=4.4 ng/mL\n")
cat("═══════════════════════════════════════════════════════════\n")
