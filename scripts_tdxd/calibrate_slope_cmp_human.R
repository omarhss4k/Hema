############################################################
# calibrate_slope_cmp_human.R
# Calibration Slope_CMP — T-DXd HUMAIN
#
# Kill linéaire (original Fornari) :
#   kill_CMP = Slope_CMP × D_kill
#   D_kill   = max(0, Damage - D0)  avec D0 = Damage_threshold = 0.05
#
# Slope_CMP >> 1 est empiriquement nécessaire pour reproduire G3-4 ~20% FDA :
#   CMP_ss = MPP_input / (Slope × k_prol_CMP + k_out)
#   À Slope=25 : CMP_ss ≈ 7% baseline → Neut nadir ≈ 0.3 × 10⁹/L (G4)
#   Avec IIV (ω_SCMP=0.33) → distribution G0/G1-2/G3-4 réaliste
#
# Cible : Grade 3-4 neutropénie ~20% (DESTINY-Breast01, n=184)
#   Grade 3 : Neut < 1.0 × 10⁹/L
#   Grade 4 : Neut < 0.5 × 10⁹/L
#
# Méthode : scan N=50 patients par valeur de Slope_CMP
#   IIV : CL_ADC (ω=0.35), V1_ADC (ω=0.20), Slope_CMP (ω=0.33)
#   Pas de temps = 12h (rapidité)
############################################################
library(deSolve)

setwd("/home/user/Hema/scripts_tdxd")
source("parameters_tdxd_rat.R")
source("parameters_tdxd_human.R")
source("pkpd_tdxd_rat.R")

set.seed(42)
N  <- 50   # patients par valeur

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

# ── Fonction : G3-4% pour un Slope_CMP_typ donné ─────────
pct_g34 <- function(slope_cmp_typ) {
  neut_min_vec <- numeric(N)
  for (i in 1:N) {
    pars_i <- pars_base
    pars_i$CL_ADC    <- pars_base$CL_ADC * exp(eta_CL[i])
    pars_i$V1_ADC    <- pars_base$V1_ADC * exp(eta_V1[i])
    pars_i$Slope_CMP <- slope_cmp_typ    * exp(eta_SCMP[i])
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
  n_ok      <- sum(!is.na(neut_min_vec))
  g34       <- sum(neut_min_vec < 1.0, na.rm=TRUE) / n_ok * 100
  g_any     <- sum(neut_min_vec < 2.0, na.rm=TRUE) / n_ok * 100
  g0        <- sum(neut_min_vec >= 2.0, na.rm=TRUE) / n_ok * 100
  list(g34 = g34, g_any = g_any, g0 = g0,
       med_nadir = median(neut_min_vec, na.rm=TRUE))
}

# ── Scan — D0=0.05, ω_SCMP=0.33 ──────────────────────────
slope_vals <- c(16, 20, 25, 28, 30, 35)

cat("═══════════════════════════════════════════════════════════\n")
cat("  SCAN Slope_CMP — Calibration T-DXd Humain (N=50/valeur)\n")
cat("  kill = Slope × D_kill (linéaire, D0=0.05)\n")
cat("  Cible : G3-4 neutropénie ~20%  (DESTINY-Breast01)\n")
cat("─────────────────────────────────────────────────────────\n")
cat(sprintf("  %-10s  %-10s  %-12s  %-10s  %-12s  %s\n",
            "Slope_CMP", "G3-4 %", "Tout grade %", "G0 %", "Neut médiane", "Statut"))
cat(sprintf("  %-10s  %-10s  %-12s  %-10s  %-12s  %s\n",
            "─────────","──────────","────────────","──────────","────────────","──────"))

res_list <- list()
for (s in slope_vals) {
  cat(sprintf("  Slope=%4.0f ...", s)); flush.console()
  r <- pct_g34(s)
  res_list[[as.character(s)]] <- r
  status <- if (abs(r$g34 - 20) < 5) "CIBLE ✓" else if (r$g34 < 20) "trop bas" else "trop haut"
  cat(sprintf("\r  %-10.0f  %-10.1f  %-12.1f  %-10.1f  %-12.2f  %s\n",
              s, r$g34, r$g_any, r$g0, r$med_nadir, status))
}

cat("═══════════════════════════════════════════════════════════\n")

g34_vec <- sapply(res_list, function(r) r$g34)
names(g34_vec) <- slope_vals

below <- slope_vals[g34_vec <= 20]
above <- slope_vals[g34_vec >= 20]

if (length(below) > 0 && length(above) > 0) {
  s_low  <- max(below);  g_low  <- g34_vec[as.character(s_low)]
  s_high <- min(above);  g_high <- g34_vec[as.character(s_high)]
  if (s_high != s_low) {
    slope_opt <- s_low + (20 - g_low) / (g_high - g_low) * (s_high - s_low)
  } else {
    slope_opt <- s_low
  }
  cat(sprintf("\n  → Slope_CMP_tdxd_human = %.1f\n", slope_opt))
  cat(sprintf("    (entre %.0f [G3-4=%.1f%%] et %.0f [G3-4=%.1f%%])\n",
              s_low, g_low, s_high, g_high))
} else if (length(below) == 0) {
  slope_opt <- slope_vals[which.min(abs(g34_vec - 20))]
  cat(sprintf("\n  ATTENTION : G3-4 > 20%% même au Slope le plus faible\n"))
  cat(sprintf("  → Slope_CMP_tdxd_human ≈ %.1f\n", slope_opt))
} else {
  slope_opt <- max(slope_vals)
  cat(sprintf("\n  ATTENTION : G3-4 < 20%% pour toutes valeurs → augmenter plage\n"))
  cat(sprintf("  → Slope_CMP_tdxd_human = %.1f (max)\n", slope_opt))
}
cat("═══════════════════════════════════════════════════════════\n")
