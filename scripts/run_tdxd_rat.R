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
cat("═══════════════════════════════════════════════════════════\n")
