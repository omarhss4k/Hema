# =============================================================================
# run_simulation.R
# Script principal d'exécution — Modèle PKPD Carboplatine (Fornari 2019)
#
# Usage :
#   Rscript run_simulation.R
#   ou sourcer dans RStudio : source("pkpd_model/run_simulation.R")
#
# Paramètres modifiables ci-dessous (SECTION CONFIGURATION)
# =============================================================================

# --- Chemin vers les modules ---
.dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "pkpd_model")

source(file.path(.dir, "01_model_ode.R"))
source(file.path(.dir, "02_parameters_library.R"))
source(file.path(.dir, "03_simulation_core.R"))
source(file.path(.dir, "04_scoring_toxicity.R"))
source(file.path(.dir, "05_visualization.R"))

# =============================================================================
# SECTION CONFIGURATION — Modifiez ici
# =============================================================================

N_PATIENTS   <- 1000         # Nombre de patients dans la population
SPECIES      <- "human"      # "human" ou "rat"
AUC_TARGET   <- 7.5          # AUC cible carboplatine (mg.min/mL)
GFR          <- 125          # GFR moyen de la population (mL/min)
INFUSION_H   <- 0.5          # Durée de perfusion (h)

# Cycles : doses à J0, J21, J42 (3 cycles q21j)
SCHEDULE_H   <- c(0, 21*24, 42*24)

T_END_DAYS   <- 90           # Durée de simulation (jours)
T_RES_H      <- 2            # Résolution (h)

SEED         <- 42           # Reproductibilité
OUT_DIR      <- "results_pkpd"

# =============================================================================
# EXÉCUTION
# =============================================================================

message("=== Modèle PKPD Carboplatine / Toxicité Hématologique ===")
message(sprintf("Espèce : %s | AUC cible : %g | %d cycles | %d patients",
                SPECIES, AUC_TARGET, length(SCHEDULE_H), N_PATIENTS))

# 1. Schéma posologique
dosing <- make_dosing(
  AUC_target = AUC_TARGET,
  GFR        = GFR,
  schedule   = SCHEDULE_H,
  infusion_h = INFUSION_H,
  species    = SPECIES
)
message(sprintf("Dose calculée : %.1f mg (%.0f µmol)", dosing$dose_mg, dosing$dose_umol))

# 2. Simulation populationnelle
pop_result <- simulate_population(
  n_patients  = N_PATIENTS,
  species     = SPECIES,
  dosing      = dosing,
  t_end       = T_END_DAYS * 24,
  t_res       = T_RES_H,
  seed        = SEED,
  patient_cov = data.frame(age=60, weight=70, creat=80, sex="M", GFR=GFR),
  verbose     = TRUE
)

# 3. Résumé populationnel (percentiles)
message("Calcul des percentiles populationnels...")
all_vars    <- c("Neut", "Plt", "Ret", "RBC", "Mono", "MPP", "CMP", "MEP", "Damage")
pop_summary <- summarize_population(pop_result, variables = all_vars)

# 4. Scores de toxicité
message("Calcul des scores de toxicité NCI-CTCAE...")
scores_df  <- score_population(pop_result, lineages = c("Neut", "Plt", "Ret", "RBC"))
grade_dist <- grade_distribution(scores_df)
nadir_summ <- nadir_summary(scores_df)

# 5. Affichage des résultats
message("\n--- Résumé des nadirs ---")
print(nadir_summ[, c("lineage","p50","p05","p95","pct_grade3_4","pct_grade4")])

message("\n--- Distribution des grades ---")
print(grade_dist)

# 6. Génération des graphiques
message("\nGénération des graphiques...")
save_all_plots(
  pop_result  = pop_result,
  scores_df   = scores_df,
  grade_dist  = grade_dist,
  pop_summary = pop_summary,
  out_dir     = OUT_DIR,
  prefix      = paste0(SPECIES, "_AUC", AUC_TARGET)
)

message("\nTerminé. Résultats dans : ", OUT_DIR)
