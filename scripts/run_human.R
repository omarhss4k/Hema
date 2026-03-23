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
GFR_mLmin    <- 78    # GFR médian population (Fornari 2019)
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
  titre = "VPC - Carboplatin (Human) AUC=5 Q21D x 2 (1000 patients, GFR=78 fixe — médian population)",
  n_sim = 1000, dose_days = dose_days_hu,
  obs_list = obs_list_human,
  auc_target = AUC_target, times = times_hu,
  interval_h = 21*24, n_cycles = 2
)

# ── VPC élargi — variabilité baselines ──────────────────
cat("\n--- VPC élargi : Neut0~LogU[2,7], Plt0~LogU[150,400] ---\n")
save_vpc_human_extended(
  sim = sim_hu, pars = init_pars,
  file = "results_HUMAN/Figure4_VPC_extended.pdf",
  titre = "VPC élargi — AUC=5 Q21D×2 | Neut0~LogU[2,7] | Plt0~LogU[150,400] ×10⁹/L",
  n_sim = 1000, dose_days = dose_days_hu,
  obs_list = obs_list_human,
  auc_target = AUC_target, times = times_hu,
  interval_h = 21*24, n_cycles = 2,
  gfr_fixed = GFR_mLmin
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

cat("\n=======================================================\n")
cat("Fichiers dans results_HUMAN/ :\n")
cat("  -> Figure4_Q21D_x2.pdf              (courbes Neut+Plt)\n")
cat("  -> Figure4_grades_deterministe.pdf  (grades courbe det.)\n")
cat("  -> Figure4_VPC.pdf                  (VPC population)\n")
cat("  -> Figure4_VPC_extended.pdf         (VPC élargi, baselines variables)\n")
cat("  -> Figure4c_grades_AUC5.pdf         (Figure 4c, 1000 patients)\n")
cat("=======================================================\n")
