############################################################
# run_human.R
# Simulation humaine — Carboplatine AUC=5, Q21D x 2 cycles
# Figure 4 + VPC + Figure 4c grades (Fornari 2019)
############################################################
source("pkpd_model_FORNARI.R")
source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")
source("plots_human.R")
source("plots_grades.R")
source("vpc.R")

if (!dir.exists("results_HUMAN")) dir.create("results_HUMAN")

# ── Dose Calvert ─────────────────────────────────────────
AUC_target   <- 5
GFR_mLmin    <- 125   # normal renal function (Fornari 2019 S11)
dose_calvert <- AUC_target * (GFR_mLmin + 25)
cat(sprintf("=== Dose Calvert : AUC=%g, GFR=%g => Dose = %.0f mg ===\n",
            AUC_target, GFR_mLmin, dose_calvert))

# ── Simulation déterministe Q21D x 2 ────────────────────
cat("\n=== Figure 4 : simulation Q21D x 2 ===\n")
pars_hu          <- init_pars
pars_hu$rate_fun <- make_repeated_infusion(
  dose_mg = dose_calvert, Tinfu_h = 1,
  interval_h = 21*24, n_cycles = 2
)
times_hu     <- seq(0, 63*24, by = 1)
dose_days_hu <- c(0, 21)

sim_hu <- simulate_all(times_hu, pars_hu, init_state)

save_human_neut_plt(
  sim = sim_hu, pars = init_pars,
  file = "results_HUMAN/Figure4_Q21D_x2.pdf",
  titre = "Hematopoiesis – Carboplatin (Human) Q21D × 2",
  dose_days = dose_days_hu, obs_list = obs_list_human
)

# ── Grades courbe déterministe ───────────────────────────
cat("\n--- Grades NCI-CTCAE : simulation déterministe ---\n")
grades_det <- save_grade_plots(
  sim = sim_hu, pars = init_pars,
  file = "results_HUMAN/Figure4_grades_deterministe.pdf",
  titre_base = "Carboplatine (Human) AUC=5 Q21D x 2",
  dose_days = dose_days_hu
)
cat(sprintf("  Neut : nadir=%.3f → Grade %d\n",
            grades_det$nadir_neut, grades_det$grade_nadir_neut))
cat(sprintf("  Plt  : nadir=%.1f  → Grade %d\n",
            grades_det$nadir_plt, grades_det$grade_nadir_plt))

# ── VPC population ───────────────────────────────────────
cat("\n--- VPC humain : AUC=5, Q21D x 2 ---\n")
save_vpc_human(
  sim = sim_hu, pars = init_pars,
  file = "results_HUMAN/Figure4_VPC.pdf",
  titre = "VPC - Carboplatin (Human) AUC=5 Q21D x 2 (1000 patients, GFR~N(78,20))",
  n_sim = 1000, dose_days = dose_days_hu,
  obs_list = obs_list_human,
  auc_target = AUC_target, times = times_hu,
  interval_h = 21*24, n_cycles = 2
)

# ── Figure 4c — % patients par grade ────────────────────
cat("\n=== Figure 4c : % patients par grade (AUC=5) ===\n")
save_grade_figure4c(
  base_pars  = init_pars,
  init_state = init_state,
  auc_target = 5,
  n_cycles   = 2,
  interval_h = 21*24,
  n_patients = 1000,
  file       = "results_HUMAN/Figure4c_grades_AUC5.pdf",
  titre      = "Carboplatine AUC=5, Q21D x 2 — % patients par grade",
  seed       = 42
)

# ── Sensibilité AUC 4 / 5 / 6 ───────────────────────────
cat("\n=== Sensibilité AUC 4/5/6 ===\n")
for (auc in c(4, 5, 6)) {
  cat(sprintf("\n  AUC=%d (dose=%.0f mg)\n", auc, auc*(GFR_mLmin+25)))

  pars_auc          <- init_pars
  pars_auc$rate_fun <- make_repeated_infusion(
    dose_mg = auc*(GFR_mLmin+25), Tinfu_h = 1,
    interval_h = 21*24, n_cycles = 2
  )
  sim_auc <- simulate_all(times_hu, pars_auc, init_state)

  save_vpc_human(
    sim = sim_auc, pars = init_pars,
    file = sprintf("results_HUMAN/Figure4_VPC_AUC%d.pdf", auc),
    titre = sprintf("VPC - Carboplatin (Human) AUC=%d Q21D x 2", auc),
    n_sim = 1000, dose_days = dose_days_hu,
    auc_target = auc, times = times_hu,
    interval_h = 21*24, n_cycles = 2
  )

  save_grade_figure4c(
    base_pars  = init_pars,
    init_state = init_state,
    auc_target = auc,
    n_cycles   = 2,
    interval_h = 21*24,
    n_patients = 1000,
    file       = sprintf("results_HUMAN/Figure4c_AUC%d.pdf", auc),
    titre      = sprintf("Carboplatine AUC=%d, Q21D x 2 — %% patients", auc),
    seed       = 42
  )
}

cat("\n=======================================================\n")
cat("Fichiers dans results_HUMAN/ :\n")
cat("  -> Figure4_Q21D_x2.pdf              (courbes Neut+Plt)\n")
cat("  -> Figure4_grades_deterministe.pdf  (grades courbe det.)\n")
cat("  -> Figure4_VPC.pdf                  (VPC population)\n")
cat("  -> Figure4c_grades_AUC5.pdf         (Figure 4c, 500 patients)\n")
cat("  -> Figure4_VPC_AUC4/5/6.pdf         (VPC sensibilite)\n")
cat("  -> Figure4c_AUC4/5/6.pdf            (grades sensibilite)\n")
cat("=======================================================\n")
