############################################################
# table6_metrics.R  --  Etape 1 : rat carboplatin
# Metriques de validation : RMSE log par lignee
# (coherent avec la fonction objectif SSR log du modele)
#
# Prerequis : run_CORRECT.R (genere sim_f3 en memoire)
#             OU lancer ce script depuis run_CORRECT.R
############################################################

source("../shared/pkpd_model_FORNARI.R")
source("parameters_rat.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("data_fornari.R")

# -- Simulation 40 mg/kg Q14D x 8 -------------------------
pars_f3 <- init_pars
pars_f3$rate_fun <- make_repeated_infusion(
  dose_mg = 40, Tinfu_h = 1, interval_h = 14*24, n_cycles = 8
)
times_f3 <- seq(0, 120*24, by = 1)
sim_f3   <- simulate_all(times_f3, pars_f3, init_state)

# -- Interpolation simulation aux temps observes -----------
interp_sim <- function(sim, var, t_obs) {
  approx(sim$days * 24, sim[[var]], xout = t_obs * 24)$y
}

# -- RMSE log-scale par lignee -----------------------------
rmse_log <- function(sim, var, obs_df) {
  # Exclut RBC t=0 artefact (valeur < 1000 x10^9/L)
  if (var == "RBC") obs_df <- obs_df[obs_df$value > 1000, ]
  pred <- interp_sim(sim, var, obs_df$time)
  obs  <- obs_df$value
  keep <- !is.na(pred) & pred > 0 & obs > 0
  sqrt(mean((log(pred[keep]) - log(obs[keep]))^2))
}

lineages <- c("Neut", "Plt", "Ret", "RBC", "Mono", "MPP", "CMP")
obs_map  <- list(
  Neut = obs_fornari$Neut,
  Plt  = obs_fornari$Plt,
  Ret  = obs_fornari$Ret,
  RBC  = obs_fornari$RBC,
  Mono = obs_fornari$Mono,
  MPP  = obs_fornari$MPP,
  CMP  = obs_fornari$CMP
)

rmse_vals <- sapply(lineages, function(v) rmse_log(sim_f3, v, obs_map[[v]]))
n_pts     <- sapply(lineages, function(v) {
  obs <- obs_map[[v]]
  if (v == "RBC") obs <- obs[obs$value > 1000, ]
  nrow(obs)
})

# -- Affichage ---------------------------------------------
cat("\n=== TABLE 6 -- Qualite d'ajustement (rat carboplatin 40 mg/kg Q14D x8) ===\n\n")
cat(sprintf("%-10s  %6s  %8s\n", "Lignee", "n obs", "RMSE (log)"))
cat(strrep("-", 28), "\n")
for (v in lineages) {
  cat(sprintf("%-10s  %6d  %8.4f\n", v, n_pts[v], rmse_vals[v]))
}
cat(strrep("-", 28), "\n")
cat(sprintf("%-10s  %6d  %8.4f\n", "GLOBAL",
            sum(n_pts), sqrt(mean(rmse_vals^2))))

# -- Export CSV --------------------------------------------
dir.create("results", showWarnings = FALSE)
tab <- data.frame(
  Lignee   = lineages,
  n_obs    = n_pts,
  RMSE_log = round(rmse_vals, 4),
  row.names = NULL
)
write.csv(tab, "results/table6_metrics.csv", row.names = FALSE)
cat("\n-> results/table6_metrics.csv\n\n")
