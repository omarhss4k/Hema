############################################################
# table6_metrics.R  --  Etape 1 : rat carboplatin
# Table 6 : RMSE rel., Residu median, Biais max par lignee
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
interp_sim <- function(sim, var, t_obs_days) {
  approx(sim$days, sim[[var]], xout = t_obs_days)$y
}

# -- Calcul des metriques relatives ------------------------
compute_metrics <- function(sim, var, obs_df, exclude_rbc_artefact = FALSE) {
  if (exclude_rbc_artefact) obs_df <- obs_df[obs_df$value > 1000, ]
  pred <- interp_sim(sim, var, obs_df$time)
  obs  <- obs_df$value
  keep <- !is.na(pred) & pred > 0 & obs > 0
  pred <- pred[keep]; obs <- obs[keep]
  rel  <- abs(pred - obs) / obs
  list(
    rmse_rel   = 100 * sqrt(mean(rel^2)),
    median_res = 100 * median(rel),
    biais_max  = 100 * max(rel),
    n          = length(obs)
  )
}

lineages <- c("MPP", "CMP", "MEP", "Neut", "Mono", "Plt", "Ret", "RBC")
obs_map  <- list(
  MPP  = obs_fornari$MPP,
  CMP  = obs_fornari$CMP,
  MEP  = obs_fornari$MEP,
  Neut = obs_fornari$Neut,
  Mono = obs_fornari$Mono,
  Plt  = obs_fornari$Plt,
  Ret  = obs_fornari$Ret,
  RBC  = obs_fornari$RBC
)

results_list <- lapply(lineages, function(v) {
  excl <- (v == "RBC")
  compute_metrics(sim_f3, v, obs_map[[v]], exclude_rbc_artefact = excl)
})
names(results_list) <- lineages

# -- Statut ------------------------------------------------
get_statut <- function(v, rmse) {
  if (v == "RBC")        return("Satisfaisant*")
  if (rmse < 10)         return("Satisfaisant")
  return("Partiel")
}

# -- Affichage ---------------------------------------------
cat("\n=== TABLE 6 -- Metriques de validation (rat carboplatin 40 mg/kg Q14D x8) ===\n\n")
cat(sprintf("%-6s  %14s  %17s  %14s  %s\n",
    "Lignee", "RMSE rel. (%)", "Residu median (%)", "Biais max (%)", "Statut"))
cat(strrep("-", 72), "\n")

tab_rows <- list()
for (v in lineages) {
  m <- results_list[[v]]
  statut <- get_statut(v, m$rmse_rel)
  symb   <- ifelse(grepl("Satisfaisant", statut), "OK", "!!")
  cat(sprintf("%-6s  %14.2f  %17.2f  %14.2f  %s %s\n",
      v, m$rmse_rel, m$median_res, m$biais_max, symb, statut))
  tab_rows[[v]] <- data.frame(
    Lignee       = v,
    RMSE_rel_pct = round(m$rmse_rel, 2),
    Residu_med   = round(m$median_res, 2),
    Biais_max    = round(m$biais_max, 2),
    Statut       = statut,
    stringsAsFactors = FALSE
  )
}
cat(strrep("-", 72), "\n")
cat("* RBC : biais max lie a l'artefact de digitalisation t=0 (exclu du calcul)\n\n")

# -- Export CSV --------------------------------------------
dir.create("results", showWarnings = FALSE)
tab <- do.call(rbind, tab_rows)
write.csv(tab, "results/table6_metrics.csv", row.names = FALSE)
cat("-> results/table6_metrics.csv\n\n")
print(tab)
