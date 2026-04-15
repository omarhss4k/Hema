############################################################
# calibrate_slope_mep_human.R
# Calibration Slope_MEP_sensitive — T-DXd HUMAIN (MIXTURE ANÉMIE)
#
# Modèle de mélange bimodal anémie (FDA DESTINY-Breast01, n=184) :
#   30% patients "résistants MEP" : Slope_MEP_resist ≈ 0.05 → G0 garanti
#   70% patients "sensibles MEP"  : Slope_MEP_sensitive → G1/G2/G3-4
#
# Cibles (dans le sous-groupe sensibles, 70%) :
#   G3-4_total    ≈ 9%    FDA DESTINY-Breast01
#   G3-4_sensible ≈ 13%   (= 9% / 70%)
#   G1_sensible   ≈ 53%   (= 37% / 70%)
#
# IIV sensibles : ω_SMEP = 0.50 (élargi vs 0.33 pour spread G1/G2/G3)
# Méthode : scan N=50 patients SENSIBLES uniquement
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
times_hu    <- seq(0, 126 * 24, by = 12)  # pas 12h

# IIV
omega_CL   <- 0.35
omega_V1   <- 0.20
omega_SMEP <- omega_Slope_MEP_sensitive   # 0.50 (paramètre chargé depuis parameters)

# Pre-tirer les etas (fixes pour tous les Slope scannés)
eta_CL   <- rnorm(N, 0, omega_CL)
eta_V1   <- rnorm(N, 0, omega_V1)
eta_SMEP <- rnorm(N, 0, omega_SMEP)

TARGET_G34_SENS <- 100 * 0.09 / p_sensitive_mep_tdxd   # 9% / 70% ≈ 12.9%
TARGET_G1_SENS  <- 100 * 0.37 / p_sensitive_mep_tdxd   # 37% / 70% ≈ 52.9%

# ── Seuils anémie (proxy RBC) ────────────────────────────
# G0 : RBC/RBC0 ≥ 0.90  (< 10% drop)
# G1 : 0.83 ≤ frac < 0.90
# G2 : 0.67 ≤ frac < 0.83
# G3 : 0.54 ≤ frac < 0.67
# G4 : frac < 0.54
ctcae_anemia_frac <- function(frac) {
  if      (frac < 0.54) "G4"
  else if (frac < 0.67) "G3"
  else if (frac < 0.83) "G2"
  else if (frac < 0.90) "G1"
  else                   "G0"
}

# ── Fonction : distribution anémie pour N patients SENSIBLES ──
scan_mep <- function(slope_mep_sens) {
  rbc_min_frac <- numeric(N)
  for (i in 1:N) {
    pars_i <- pars_base
    pars_i$CL_ADC    <- pars_base$CL_ADC * exp(eta_CL[i])
    pars_i$V1_ADC    <- pars_base$V1_ADC * exp(eta_V1[i])
    pars_i$Slope_MEP <- slope_mep_sens * exp(eta_SMEP[i])
    # Inclure la mixture CMP pour reproduire l'interaction MPP→MEP dans la pop. complète
    pars_i$Slope_CMP <- if (runif(1) < p_sensitive_tdxd) {
      Slope_sensitive_tdxd * exp(rnorm(1, 0, 0.33))
    } else {
      Slope_resist_tdxd * exp(rnorm(1, 0, 0.33))
    }
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
    if (!is.null(out)) {
      rbc_min <- min(out$RBC, na.rm = TRUE)
      rbc_min_frac[i] <- rbc_min / pars_base$RBC0
    } else {
      rbc_min_frac[i] <- NA_real_
    }
  }
  n_ok    <- sum(!is.na(rbc_min_frac))
  grades  <- sapply(rbc_min_frac[!is.na(rbc_min_frac)], ctcae_anemia_frac)
  g0_s    <- sum(grades == "G0") / n_ok * 100
  g1_s    <- sum(grades == "G1") / n_ok * 100
  g2_s    <- sum(grades == "G2") / n_ok * 100
  g34_s   <- sum(grades %in% c("G3","G4")) / n_ok * 100
  g34_tot <- g34_s   * p_sensitive_mep_tdxd
  g1_tot  <- g1_s    * p_sensitive_mep_tdxd
  list(g0_s = g0_s, g1_s = g1_s, g2_s = g2_s, g34_s = g34_s,
       g34_tot = g34_tot, g1_tot = g1_tot,
       med_frac = median(rbc_min_frac, na.rm = TRUE))
}

# ── Scan ────────────────────────────────────────────────
slope_vals <- c(0.2, 0.3, 0.4, 0.55, 0.66, 0.80)

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  SCAN Slope_MEP_sensitive — Calibration mixture anémie (N=50 sens.)\n")
cat(sprintf("  ω_SMEP = %.2f  |  Slope_MEP_resist = %.2f\n",
            omega_SMEP, Slope_MEP_resist_tdxd))
cat(sprintf("  Cible G3-4 sensibles : %.1f%%  (= 9%% / 70%%)\n", TARGET_G34_SENS))
cat(sprintf("  Cible G1  sensibles  : %.1f%%  (= 37%% / 70%%)\n", TARGET_G1_SENS))
cat("──────────────────────────────────────────────────────────────────────\n")
cat(sprintf("  %-10s  %-8s  %-8s  %-8s  %-10s  %-10s  %-10s\n",
            "Slope_MEP", "G0_s%", "G1_s%", "G2_s%", "G3-4_s%", "G3-4_tot%", "RBC_frac"))
cat(sprintf("  %-10s  %-8s  %-8s  %-8s  %-10s  %-10s  %-10s\n",
            "─────────","──────","──────","──────","────────","──────────","────────"))

res_list <- list()
for (s in slope_vals) {
  cat(sprintf("  Slope=%.2f ...", s)); flush.console()
  r <- scan_mep(s)
  res_list[[as.character(s)]] <- r
  diff34 <- abs(r$g34_s - TARGET_G34_SENS)
  status <- if (diff34 < 4) "CIBLE ✓" else if (r$g34_s < TARGET_G34_SENS) "trop bas" else "trop haut"
  cat(sprintf("\r  %-10.2f  %-8.1f  %-8.1f  %-8.1f  %-10.1f  %-10.1f  %-10.3f  %s\n",
              s, r$g0_s, r$g1_s, r$g2_s, r$g34_s, r$g34_tot, r$med_frac, status))
}
cat("══════════════════════════════════════════════════════════════════════\n")

g34_s_vec <- sapply(res_list, function(r) r$g34_s)
names(g34_s_vec) <- slope_vals

below <- slope_vals[g34_s_vec <= TARGET_G34_SENS]
above <- slope_vals[g34_s_vec >= TARGET_G34_SENS]

if (length(below) > 0 && length(above) > 0) {
  s_low  <- max(below);  g_low  <- g34_s_vec[as.character(s_low)]
  s_high <- min(above);  g_high <- g34_s_vec[as.character(s_high)]
  slope_opt <- if (s_high != s_low) {
    s_low + (TARGET_G34_SENS - g_low) / (g_high - g_low) * (s_high - s_low)
  } else s_low
  cat(sprintf("\n  → Slope_MEP_sensitive_tdxd = %.2f\n", slope_opt))
  cat(sprintf("    (entre %.2f [G3-4_s=%.1f%%] et %.2f [G3-4_s=%.1f%%])\n",
              s_low, g_low, s_high, g_high))
  cat(sprintf("    → G3-4 anémie total projeté ≈ %.1f%%  (cible 9%%)\n",
              TARGET_G34_SENS * p_sensitive_mep_tdxd))
} else if (length(below) == 0) {
  slope_opt <- slope_vals[which.min(abs(g34_s_vec - TARGET_G34_SENS))]
  cat(sprintf("\n  ATTENTION : G3-4_s > %.0f%% même au Slope minimum\n", TARGET_G34_SENS))
  cat(sprintf("  → Slope_MEP_sensitive_tdxd ≈ %.2f (minimum)\n", slope_opt))
} else {
  slope_opt <- max(slope_vals)
  cat(sprintf("\n  ATTENTION : G3-4_s < %.0f%% pour toutes valeurs\n", TARGET_G34_SENS))
  cat(sprintf("  → Slope_MEP_sensitive_tdxd = %.2f (max testé)\n", slope_opt))
}
cat("══════════════════════════════════════════════════════════════════════\n")
