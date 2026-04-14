############################################################
# calibrate_slope_mep_human.R
# Calibration de Slope_MEP pour T-DXd humain
#
# Objectif : trouver Slope_MEP tel que la simulation de
#   population (N=100) reproduise G3-4 anémie ~9%
#   (FDA BLA 761139, DESTINY-Breast01, n=184)
#
# Méthode : scan de Slope_MEP, même protocole que la
#   calibration de Slope_CMP (Slope=12 → G3-4 Neut=15.7%)
#
# Usage : Rscript calibrate_slope_mep_human.R
############################################################
library(deSolve)

source("parameters_tdxd_rat.R")
source("parameters_tdxd_human.R")
source("pkpd_tdxd_rat.R")

set.seed(123)
N_scan      <- 100                        # patients par valeur (rapide)
slope_grid  <- seq(0.5, 8.0, by = 0.5)   # valeurs à tester

# ── Paramètres typiques ──────────────────────────────────
pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL
pars_typ   <- c(pars_pd_hu, tdxd_pars_hu)
pars_typ$Slope_CMP <- Slope_CMP_tdxd_human

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)

times_hu <- seq(0, 126 * 24, by = 6)   # 6 cycles Q3W, pas 6h

# ── IIV (mêmes ω que la simulation population) ──────────
omega_CL   <- 0.35
omega_V1   <- 0.20
omega_SMEP <- 0.33

# ── Fonction proxy anémie (identique à la population) ───
ctcae_anemia <- function(rbc, rbc0) {
  frac <- rbc / rbc0
  if      (frac < 0.54) "G4"
  else if (frac < 0.67) "G3"
  else if (frac < 0.83) "G2"
  else if (frac < 0.90) "G1"
  else                   "G0"
}

# ── Scan ─────────────────────────────────────────────────
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  Calibration Slope_MEP — T-DXd humain                  ║\n")
cat("║  Cible : G3-4 anémie ~9% (FDA DESTINY-Breast01)        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("  %-10s  %-12s  %-12s  %s\n",
            "Slope_MEP", "G3-4 anémie", "G1-4 anémie", "Statut"))
cat(sprintf("  %-10s  %-12s  %-12s  %s\n",
            "──────────","────────────","────────────","──────"))

results_scan <- data.frame(
  Slope_MEP   = slope_grid,
  pct_G34     = NA_real_,
  pct_G14     = NA_real_
)

for (s_idx in seq_along(slope_grid)) {
  sl <- slope_grid[s_idx]
  pars_sl <- pars_typ
  pars_sl$Slope_MEP <- sl

  grades <- character(N_scan)

  for (i in 1:N_scan) {
    eta_CL   <- rnorm(1, 0, omega_CL)
    eta_V1   <- rnorm(1, 0, omega_V1)
    eta_SMEP <- rnorm(1, 0, omega_SMEP)

    pars_i <- pars_sl
    pars_i$CL_ADC    <- pars_sl$CL_ADC  * exp(eta_CL)
    pars_i$V1_ADC    <- pars_sl$V1_ADC  * exp(eta_V1)
    pars_i$Slope_MEP <- pars_sl$Slope_MEP * exp(eta_SMEP)
    pars_i$rate_fun  <- make_tdxd_infusion(
      dose_mgkg = 5.4, BW_kg = 70, Tinfu_h = 1.5,
      interval_h = 21 * 24, n_cycles = 6
    )

    out <- tryCatch(
      suppressMessages(suppressWarnings(
        as.data.frame(lsoda(
          y = state0_hu, times = times_hu,
          func = pkpd_tdxd_fornari, parms = pars_i,
          rtol = 1e-5, atol = 1e-7, maxsteps = 300000
        ))
      )),
      error = function(e) NULL
    )

    if (!is.null(out) && nrow(out) > 10) {
      grades[i] <- ctcae_anemia(min(out$RBC, na.rm=TRUE), pars_typ$RBC0)
    } else {
      grades[i] <- "FAIL"
    }
  }

  n_ok  <- sum(grades != "FAIL")
  p_G34 <- 100 * sum(grades %in% c("G3","G4")) / n_ok
  p_G14 <- 100 * sum(grades %in% c("G1","G2","G3","G4")) / n_ok

  results_scan$pct_G34[s_idx] <- p_G34
  results_scan$pct_G14[s_idx] <- p_G14

  # Statut par rapport à la cible FDA
  statut <- if (p_G34 >= 7 && p_G34 <= 11) "✓ dans la cible" else
            if (p_G34 < 7)  "↓ sous la cible" else
                             "↑ au-dessus"

  cat(sprintf("  %-10.1f  %9.1f%%    %9.1f%%    %s\n",
              sl, p_G34, p_G14, statut))
}

# ── Résultat optimal ─────────────────────────────────────
idx_best <- which.min(abs(results_scan$pct_G34 - 9.0))
best_sl  <- results_scan$Slope_MEP[idx_best]
best_p   <- results_scan$pct_G34[idx_best]

cat(sprintf("\n  → Slope_MEP optimal : %.1f  (G3-4 anémie = %.1f%%)\n",
            best_sl, best_p))
cat("  → Mettre à jour parameters_tdxd_human.R :\n")
cat(sprintf("      Slope_MEP_tdxd_human <- %.1f\n\n", best_sl))

# Graphique
pdf("results_PKPD_human/calibrate_slope_mep_human.pdf", width = 8, height = 5)
par(mar = c(4, 4.2, 3, 1))
plot(results_scan$Slope_MEP, results_scan$pct_G34,
     type = "b", pch = 19, lwd = 2, col = "#d6604d",
     xlab = "Slope_MEP", ylab = "G3-4 anémie (%)",
     main = "Calibration Slope_MEP — T-DXd humain",
     ylim = c(0, max(results_scan$pct_G34) * 1.2))
abline(h = 9,  lty = 2, col = "red",    lwd = 1.5)
abline(h = 7,  lty = 3, col = "orange", lwd = 1.0)
abline(h = 11, lty = 3, col = "orange", lwd = 1.0)
text(max(slope_grid) * 0.8, 9.5, "FDA cible = 9%", col = "red", cex = 0.85)
abline(v = best_sl, lty = 2, col = "#2166ac", lwd = 1.5)
text(best_sl + 0.1, max(results_scan$pct_G34) * 0.9,
     sprintf("Optimal\n%.1f", best_sl), col = "#2166ac", cex = 0.85)
dev.off()
cat("  -> results_PKPD_human/calibrate_slope_mep_human.pdf\n")
