############################################################
# table6_metrics.R
# Génère le Tableau 6 : métriques de validation
# (proportions simulées vs FDA + RMSE/MAE)
#
# Prérequis : lancer d'abord run_pkpd_tdxd_human_population.R
#             -> results/population_results.rds
############################################################

# -- Chargement ----------------------------------------
rds_path <- "results/population_results.rds"
if (!file.exists(rds_path))
  stop("Lancez d'abord run_pkpd_tdxd_human_population.R")

results <- readRDS(rds_path)
n <- nrow(results)

# -- Référence FDA (BLA 761139, DESTINY-Breast01, n=184) ----
fda_neut   <- c(G0=71, G1=7,  G2=7,  G3=13, G4=3)
fda_anemia <- c(G0=30, G1=37, G2=24, G3=8,  G4=1)
fda_plt    <- c(G0=63, G1=30, G2=4,  G3=2,  G4=1)

fda_tg  <- c(Neut=100-fda_neut["G0"],   Anemia=100-fda_anemia["G0"], Plt=100-fda_plt["G0"])
fda_g34 <- c(Neut=fda_neut["G3"]+fda_neut["G4"],
             Anemia=fda_anemia["G3"]+fda_anemia["G4"],
             Plt=fda_plt["G3"]+fda_plt["G4"])

# -- Proportions modèle ------------------------------------
grade_order <- c("G0","G1","G2","G3","G4")
pct <- function(col) {
  tab <- table(factor(results[[col]], levels=grade_order))
  setNames(as.numeric(100 * tab / n), grade_order)
}
pn <- pct("Grade_Neut")
pa <- pct("Grade_Anemia")
pp <- pct("Grade_Plt")

mod_tg  <- c(Neut=sum(pn[c("G1","G2","G3","G4")]),
             Anemia=sum(pa[c("G1","G2","G3","G4")]),
             Plt=sum(pp[c("G1","G2","G3","G4")]))
mod_g34 <- c(Neut=sum(pn[c("G3","G4")]),
             Anemia=sum(pa[c("G3","G4")]),
             Plt=sum(pp[c("G3","G4")]))

# -- IC95% Wilson ------------------------------------------
wilson <- function(pct_val, N) {
  p <- pct_val / 100
  m <- 1.96 * sqrt(p*(1-p)/N) * 100
  sprintf("[%.1f - %.1f]", max(0, pct_val-m), min(100, pct_val+m))
}

# -- RMSE / MAE --------------------------------------------
pred_vec <- c(mod_tg["Neut"],   mod_g34["Neut"],
              mod_tg["Anemia"], mod_g34["Anemia"],
              mod_tg["Plt"],    mod_g34["Plt"])
obs_vec  <- c(fda_tg["Neut"],   fda_g34["Neut"],
              fda_tg["Anemia"], fda_g34["Anemia"],
              fda_tg["Plt"],    fda_g34["Plt"])

rmse_val <- sqrt(mean((pred_vec - obs_vec)^2))
mae_val  <- mean(abs(pred_vec - obs_vec))

# -- Tableau 6 (console) -----------------------------------
tox_names <- c("Neutropenie", "Neutropenie", "Anemie", "Anemie", "Thrombocytopenie", "Thrombocytopenie")
metriques  <- c("Tout grade", "G3-4", "Tout grade", "G3-4", "Tout grade", "G3-4")
mod_vals   <- c(mod_tg["Neut"],   mod_g34["Neut"],
                mod_tg["Anemia"], mod_g34["Anemia"],
                mod_tg["Plt"],    mod_g34["Plt"])
fda_vals   <- c(fda_tg["Neut"],   fda_g34["Neut"],
                fda_tg["Anemia"], fda_g34["Anemia"],
                fda_tg["Plt"],    fda_g34["Plt"])

cat("\n=== TABLEAU 6 -- Metriques de validation (N =", n, ") ===\n\n")
cat(sprintf("%-18s %-12s  %8s  %-16s  %8s\n",
            "Toxicite", "Metrique", "Modele%", "IC95% [Wilson]", "FDA%"))
cat(strrep("-", 68), "\n")
for (i in seq_along(tox_names)) {
  cat(sprintf("%-18s %-12s  %7.1f%%  %-16s  %7.1f%%\n",
              tox_names[i], metriques[i],
              mod_vals[i], wilson(mod_vals[i], n),
              fda_vals[i]))
}
cat(strrep("-", 68), "\n")
cat(sprintf("%-32s  RMSE = %.2f pp\n", "Performance globale", rmse_val))
cat(sprintf("%-32s  MAE  = %.2f pp\n", "", mae_val))

# -- Export CSV --------------------------------------------
dir.create("results", showWarnings=FALSE)
tab <- data.frame(
  Toxicite  = tox_names,
  Metrique  = metriques,
  Modele_pct = round(mod_vals, 1),
  IC95_Wilson = sapply(mod_vals, wilson, N=n),
  FDA_pct    = fda_vals,
  row.names  = NULL
)
write.csv(tab, "results/table6_metrics.csv", row.names=FALSE)
cat("\n-> results/table6_metrics.csv\n")
cat(sprintf("-> RMSE = %.2f pp  |  MAE = %.2f pp\n\n", rmse_val, mae_val))
