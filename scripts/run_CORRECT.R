############################################################
# run_CORRECT.R
# Simulations rat + VPC + grades NCI-CTCAE
############################################################
source("pkpd_model_FORNARI.R")
source("parameters_FORNARI_CORRECT.R")
source("plots.R")
source("plots_grades.R")
source("vpc.R")
source("data_fornari.R")
if (!dir.exists("results_CORRECT")) dir.create("results_CORRECT")

# ============================================================
# Figure 3 — 40 mg/kg Q14D x 8 cycles (120 jours)
# ============================================================
cat("=== Figure 3 : 40 mg/kg Q14D x 8 cycles ===\n")
pars_f3          <- init_pars
pars_f3$rate_fun <- make_repeated_infusion(
  dose_mg = 40, Tinfu_h = 1, interval_h = 14*24, n_cycles = 8
)
times_f3     <- seq(0, 120*24, by = 1)
dose_days_f3 <- seq(0, 7*14,   by = 14)

sim_f3 <- simulate_all(times_f3, pars_f3, init_state)
print_summary(sim_f3, init_pars)

save_all_cells(
  sim = sim_f3, pars = init_pars,
  file = "results_CORRECT/Figure3_simulation.pdf",
  titre = "Figure 3 - Carboplatin 40 mg/kg Q14D x 8 (simulation)",
  dose_days = dose_days_f3
)

save_all_cells(
  sim = sim_f3, pars = init_pars,
  file = "results_CORRECT/Figure3_overlay.pdf",
  titre = "Figure 3 - Simulation vs donnees Fornari 2019",
  dose_days = dose_days_f3, obs_list = obs_fornari
)

cat("\n--- VPC : 40 mg/kg Q14D x 8 cycles ---\n")
save_vpc(
  sim = sim_f3, pars = init_pars,
  file = "results_CORRECT/Figure3_VPC.pdf",
  titre = "VPC - Carboplatin 40 mg/kg Q14D x 8 (1000 simulations)",
  n_sim = 1000, dose_days = dose_days_f3, obs_list = obs_fornari
)

cat("\n--- Grades NCI-CTCAE : 40 mg/kg Q14D x 8 ---\n")
grades_f3 <- save_grade_plots(
  sim = sim_f3, pars = init_pars,
  file = "results_CORRECT/Figure3_grades.pdf",
  titre_base = "Carboplatine 40 mg/kg Q14D x 8 cycles",
  dose_days = dose_days_f3
)
cat(sprintf("  Neut : nadir=%.3f → Grade %d | G3: %.1f%%  G4: %.1f%%\n",
            grades_f3$nadir_neut, grades_f3$grade_nadir_neut,
            grades_f3$pct_neut[3], grades_f3$pct_neut[4]))
cat(sprintf("  Plt  : nadir=%.1f  → Grade %d | G3: %.1f%%  G4: %.1f%%\n",
            grades_f3$nadir_plt, grades_f3$grade_nadir_plt,
            grades_f3$pct_plt[3], grades_f3$pct_plt[4]))

# ============================================================
# Figure S1 — 30 mg/kg dose unique
# ============================================================
cat("\n=== Figure S1 : 30 mg/kg dose unique ===\n")
pars_s1          <- init_pars
pars_s1$rate_fun <- make_repeated_infusion(
  dose_mg = 30, Tinfu_h = 1, interval_h = 9999, n_cycles = 1
)
times_s1 <- seq(0, 30*24, by = 0.5)

sim_s1 <- simulate_all(times_s1, pars_s1, init_state)
print_summary(sim_s1, init_pars)

save_all_cells(
  sim = sim_s1, pars = init_pars,
  file = "results_CORRECT/FigureS1_simulation.pdf",
  titre = "Figure S1 - Carboplatin 30 mg/kg single dose",
  dose_days = 0
)

save_vpc(
  sim = sim_s1, pars = init_pars,
  file = "results_CORRECT/FigureS1_VPC.pdf",
  titre = "VPC - Carboplatin 30 mg/kg single dose (1000 simulations)",
  n_sim = 1000, dose_days = 0
)

cat("\n--- Grades NCI-CTCAE : 30 mg/kg dose unique ---\n")
grades_s1 <- save_grade_plots(
  sim = sim_s1, pars = init_pars,
  file = "results_CORRECT/FigureS1_grades.pdf",
  titre_base = "Carboplatine 30 mg/kg dose unique",
  dose_days = 0
)
cat(sprintf("  Neut : nadir=%.3f → Grade %d\n",
            grades_s1$nadir_neut, grades_s1$grade_nadir_neut))
cat(sprintf("  Plt  : nadir=%.1f  → Grade %d\n",
            grades_s1$nadir_plt, grades_s1$grade_nadir_plt))

cat("\n=======================================================\n")
cat("Fichiers dans results_CORRECT/ :\n")
cat("  -> Figure3_simulation.pdf\n")
cat("  -> Figure3_overlay.pdf\n")
cat("  -> Figure3_VPC.pdf\n")
cat("  -> Figure3_grades.pdf\n")
cat("  -> FigureS1_simulation.pdf\n")
cat("  -> FigureS1_VPC.pdf\n")
cat("  -> FigureS1_grades.pdf\n")
cat("=======================================================\n")
