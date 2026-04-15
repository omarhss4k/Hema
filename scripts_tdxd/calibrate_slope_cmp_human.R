############################################################
# calibrate_slope_cmp_human.R
# Calibration Slope_sensitive — T-DXd HUMAIN (MIXTURE MODEL)
#
# Modèle de mélange bimodal (FDA DESTINY-Breast01, n=184) :
#   71% patients "résistants" : Slope_resist ≈ 0.10 → G0 garanti
#   29% patients "sensibles"  : Slope_sensitive >> 1 → G3-4 possible
#
# Cible :
#   G3-4 TOTAL    ≈ 20%  (FDA DESTINY-Breast01)
#   G3-4 SENSIBLE ≈ 69%  (= 20% / 29%)
#   Tout grade    ≈ 29%  (seuls les sensibles peuvent tomber < 2.0)
#
# Méthode : scan N=50 patients SENSIBLES uniquement
#   IIV : CL_ADC (ω=0.35), V1_ADC (ω=0.20), Slope_CMP (ω=0.33)
#   Target : G3-4_sensitive ≈ 69%
############################################################
library(deSolve)

setwd("/home/user/Hema/scripts_tdxd")
source("parameters_tdxd_rat.R")
source("parameters_tdxd_human.R")
source("pkpd_tdxd_rat.R")

set.seed(42)
N  <- 50   # patients sensibles par valeur

# ── Paramètres base ──────────────────────────────────────
pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL
pars_base <- c(pars_pd_hu, tdxd_pars_hu)

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)
times_hu    <- seq(0, 126 * 24, by = 12)  # pas 12h pour rapidité

# IIV
omega_CL   <- 0.35
omega_V1   <- 0.20
omega_SCMP <- 0.33

# Pre-tirer les etas (mêmes pour tous les Slope scannés → comparaison propre)
eta_CL   <- rnorm(N, 0, omega_CL)
eta_V1   <- rnorm(N, 0, omega_V1)
eta_SCMP <- rnorm(N, 0, omega_SCMP)

# ── Fonction : G3-4% pour N patients SENSIBLES avec Slope_sensitive_typ ──
pct_g34_sensitive <- function(slope_sensitive_typ) {
  neut_min_vec <- numeric(N)
  for (i in 1:N) {
    pars_i <- pars_base
    pars_i$CL_ADC    <- pars_base$CL_ADC * exp(eta_CL[i])
    pars_i$V1_ADC    <- pars_base$V1_ADC * exp(eta_V1[i])
    pars_i$Slope_CMP <- slope_sensitive_typ * exp(eta_SCMP[i])
    pars_i$rate_fun  <- make_tdxd_infusion(
      dose_mgkg = 5.4, BW_kg = 70, Tinfu_h = 1.5,
      interval_h = 21 * 24, n_cycles = 6
    )
    out <- tryCatch(
      suppressMessages(suppressWarnings(
        as.data.frame(lsoda(
          y = state0_hu, times = times_hu,
          func = pkpd_tdxd_fornari, parms = pars_i,
          rtol = 1e-4, atol = 1e-6, maxsteps = 500000
        ))
      )),
      error = function(e) NULL
    )
    neut_min_vec[i] <- if (!is.null(out)) min(out$Neut, na.rm=TRUE) else NA_real_
  }
  n_ok        <- sum(!is.na(neut_min_vec))
  g34_sens    <- sum(neut_min_vec < 1.0, na.rm=TRUE) / n_ok * 100
  g_any_sens  <- sum(neut_min_vec < 2.0, na.rm=TRUE) / n_ok * 100
  g34_total   <- g34_sens  * p_sensitive_tdxd  # projection tout-venant
  g_any_total <- g_any_sens * p_sensitive_tdxd
  list(g34_sens   = g34_sens,
       g_any_sens = g_any_sens,
       g34_total  = g34_total,
       g_any_total= g_any_total,
       med_nadir  = median(neut_min_vec, na.rm=TRUE))
}

# ── Scan Slope_sensitive — cible G3-4_sensible ≈ 69% (= 20%/29%) ──────────
slope_vals <- c(25, 35, 50, 70, 100, 150)
TARGET_G34_SENS <- 100 * 0.20 / p_sensitive_tdxd   # ≈ 69%

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  SCAN Slope_sensitive — Calibration mixture T-DXd (N=50 sensibles)\n")
cat("  kill = Slope × D_kill (linéaire, D0=0.05)\n")
cat(sprintf("  Cible G3-4 sensibles : %.1f%%  (= 20%% / %.0f%% sensibles)\n",
            TARGET_G34_SENS, 100 * p_sensitive_tdxd))
cat("  Cible G3-4 total     : ~20%%  (DESTINY-Breast01)\n")
cat("───────────────────────────────────────────────────────────────────\n")
cat(sprintf("  %-12s  %-12s  %-12s  %-12s  %s\n",
            "Slope_sens", "G3-4 sens%", "G3-4 tot%", "Neut méd.", "Statut"))
cat(sprintf("  %-12s  %-12s  %-12s  %-12s  %s\n",
            "───────────","──────────","──────────","──────────","──────"))

res_list <- list()
for (s in slope_vals) {
  cat(sprintf("  Slope=%5.0f ...", s)); flush.console()
  r <- pct_g34_sensitive(s)
  res_list[[as.character(s)]] <- r
  diff <- abs(r$g34_sens - TARGET_G34_SENS)
  status <- if (diff < 7) "CIBLE ✓" else if (r$g34_sens < TARGET_G34_SENS) "trop bas" else "trop haut"
  cat(sprintf("\r  %-12.0f  %-12.1f  %-12.1f  %-12.2f  %s\n",
              s, r$g34_sens, r$g34_total, r$med_nadir, status))
}

cat("═══════════════════════════════════════════════════════════════════\n")

g34_sens_vec <- sapply(res_list, function(r) r$g34_sens)
names(g34_sens_vec) <- slope_vals

below <- slope_vals[g34_sens_vec <= TARGET_G34_SENS]
above <- slope_vals[g34_sens_vec >= TARGET_G34_SENS]

if (length(below) > 0 && length(above) > 0) {
  s_low  <- max(below);  g_low  <- g34_sens_vec[as.character(s_low)]
  s_high <- min(above);  g_high <- g34_sens_vec[as.character(s_high)]
  if (s_high != s_low) {
    slope_opt <- s_low + (TARGET_G34_SENS - g_low) / (g_high - g_low) * (s_high - s_low)
  } else {
    slope_opt <- s_low
  }
  cat(sprintf("\n  → Slope_sensitive_tdxd = %.1f\n", slope_opt))
  cat(sprintf("    (entre %.0f [G3-4_sens=%.1f%%] et %.0f [G3-4_sens=%.1f%%])\n",
              s_low, g_low, s_high, g_high))
  cat(sprintf("    → G3-4 total projeté ≈ %.1f%%  (cible 20%%)\n",
              slope_opt * 0 + TARGET_G34_SENS * p_sensitive_tdxd))
} else if (length(below) == 0) {
  slope_opt <- slope_vals[which.min(abs(g34_sens_vec - TARGET_G34_SENS))]
  cat(sprintf("\n  ATTENTION : G3-4_sens > %.0f%% même au Slope le plus faible\n", TARGET_G34_SENS))
  cat(sprintf("  → Slope_sensitive_tdxd ≈ %.1f (minimum)\n", slope_opt))
} else {
  slope_opt <- max(slope_vals)
  cat(sprintf("\n  ATTENTION : G3-4_sens < %.0f%% pour toutes valeurs → augmenter plage\n", TARGET_G34_SENS))
  cat(sprintf("  → Slope_sensitive_tdxd = %.1f (max testé)\n", slope_opt))
}
cat("═══════════════════════════════════════════════════════════════════\n")
