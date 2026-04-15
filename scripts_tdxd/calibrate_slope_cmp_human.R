############################################################
# calibrate_slope_cmp_human.R
# Calibration Emax_CMP_kill — T-DXd HUMAIN
#
# Cible : Grade 3-4 neutropénie ~20% (DESTINY-Breast01, n=184)
#   Grade 3 : Neut < 1.0 × 10⁹/L
#   Grade 4 : Neut < 0.5 × 10⁹/L
#
# Kill Emax (nouveau) :
#   kill_CMP = Emax_CMP × D_kill / (ED50_kill + D_kill) ∈ [0, 1]
#   D_kill = max(0, Damage - D0)  avec D0 = 0.05
#   → kill borné à 1 → pas de prolifération négative
#   → IIV sur Emax_CMP crée distribution de sensibilité
#
# Méthode : scan N=50 patients par valeur de Emax_CMP
#   IIV : CL_ADC (ω=0.35), V1_ADC (ω=0.20), Emax_CMP (ω=0.40)
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
# Override ED50_kill : 0.05 (plus potent que 0.15)
# kill(D_kill=0.41) = Emax × 0.41/(0.05+0.41) = Emax × 0.891
# vs ED50=0.15 : Emax × 0.732 → insuffisant pour G3-4
# ED50_kill = 0.01 : courbe quasi-seuil (kill ≈ 0.98 au pic D_kill=0.41)
# kill(D_kill=0.41) = Emax × 0.41/0.42 = Emax × 0.976
# Damage_trough ≈ 0.0001 → kill_trough ≈ 0.001 (négligeable) → pas besoin de D0
pars_base$ED50_kill <- 0.01
pars_base$Damage_threshold <- 0   # D0=0 : Damage_trough≈0 → protection naturelle

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)
times_hu    <- seq(0, 126 * 24, by = 12)  # pas 12h pour rapidité

# IIV
omega_CL   <- 0.35
omega_V1   <- 0.20
# omega_Emax = 0.80 : IIV large → distribution bimodale "résistants" vs "sensibles"
# FDA DESTINY-Breast01 : ~45% G0, ~35% G1-2, ~20% G3-4
# P(Emax_eff>1) = P(eta > ln(1/emax_typ)/0.80) → fraction de patients à kill max
omega_Emax <- 0.80

# Pre-tirer les etas (mêmes pour tous les Emax scannés → comparaison propre)
eta_CL   <- rnorm(N, 0, omega_CL)
eta_V1   <- rnorm(N, 0, omega_V1)
eta_Emax <- rnorm(N, 0, omega_Emax)

# ── Fonction : G3-4% pour un Emax_CMP_typ donné ──────────
pct_g34 <- function(emax_cmp_typ) {
  neut_min_vec <- numeric(N)
  for (i in 1:N) {
    pars_i <- pars_base
    pars_i$CL_ADC      <- pars_base$CL_ADC * exp(eta_CL[i])
    pars_i$V1_ADC      <- pars_base$V1_ADC * exp(eta_V1[i])
    # Emax capé à 1 (kill max = arrêt complet de prolifération)
    emax_i <- min(1.0, emax_cmp_typ * exp(eta_Emax[i]))
    pars_i$Emax_CMP_kill <- emax_i
    # Emax_MPP et MEP scalés selon ratios Fornari (Table 2)
    pars_i$Emax_MPP_kill <- min(1.0, 0.55 * emax_i)
    pars_i$Emax_MEP_kill <- min(1.0, 0.85 * emax_i)
    pars_i$rate_fun <- make_tdxd_infusion(
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

# ── Scan Emax_CMP_kill ────────────────────────────────────
# D_kill_max = Damage_max ≈ 0.46 (D0=0), ED50=0.01 → E(0.46) = 0.979
# kill_max = Emax × 0.979 → quasi-seuil : kill ON dès Damage > 0.01
# IIV large (ω=0.80) crée distribution bimodale résistants/sensibles
# Plage : P(Emax_eff>1) = P(eta > -ln(emax)/0.80)
emax_vals <- c(0.25, 0.35, 0.45, 0.55, 0.65, 0.75)

cat("═══════════════════════════════════════════════════════════\n")
cat("  SCAN Emax_CMP_kill — Calibration T-DXd Humain (N=50/val)\n")
cat("  Kill Emax : kill = Emax × D_kill / (ED50 + D_kill) ≤ 1\n")
cat(sprintf("  ED50_kill = %.2f  D_kill_max ≈ 0.41  E_max_factor = %.3f\n",
            pars_base$ED50_kill, 0.41 / (pars_base$ED50_kill + 0.41)))
cat("  Cible : G3-4 neutropénie ~20%  (DESTINY-Breast01)\n")
cat("─────────────────────────────────────────────────────────\n")
cat(sprintf("  %-12s  %-10s  %-12s  %-10s  %-12s  %s\n",
            "Emax_CMP", "G3-4 %", "Tout grade %", "G0 %", "Neut médiane", "Statut"))
cat(sprintf("  %-12s  %-10s  %-12s  %-10s  %-12s  %s\n",
            "────────────","──────────","────────────","──────────","────────────","──────"))

res_list <- list()
for (e in emax_vals) {
  cat(sprintf("  Emax=%5.2f ...", e)); flush.console()
  r <- pct_g34(e)
  res_list[[as.character(e)]] <- r
  status <- if (abs(r$g34 - 20) < 5) "CIBLE ✓" else if (r$g34 < 20) "trop bas" else "trop haut"
  cat(sprintf("\r  %-12.2f  %-10.1f  %-12.1f  %-10.1f  %-12.2f  %s\n",
              e, r$g34, r$g_any, r$g0, r$med_nadir, status))
}

cat("═══════════════════════════════════════════════════════════\n")

# ── Trouver la valeur optimale par interpolation ──────────
g34_vec <- sapply(res_list, function(r) r$g34)
names(g34_vec) <- emax_vals

below <- emax_vals[g34_vec <= 20]
above <- emax_vals[g34_vec >= 20]

if (length(below) > 0 && length(above) > 0) {
  e_low  <- max(below);  g_low  <- g34_vec[as.character(e_low)]
  e_high <- min(above);  g_high <- g34_vec[as.character(e_high)]
  if (e_high != e_low) {
    emax_opt <- e_low + (20 - g_low) / (g_high - g_low) * (e_high - e_low)
  } else {
    emax_opt <- e_low
  }
  cat(sprintf("\n  Interpolation : Emax_CMP_tdxd_human ≈ %.2f\n", emax_opt))
  cat(sprintf("  (entre %.2f [G3-4=%.1f%%] et %.2f [G3-4=%.1f%%])\n",
              e_low, g_low, e_high, g_high))
} else if (length(below) == 0) {
  cat("\n  ATTENTION : même Emax=0.50 donne G3-4 > 20% — réduire la plage\n")
  emax_opt <- emax_vals[which.min(abs(g34_vec - 20))]
} else {
  cat("\n  ATTENTION : G3-4 < 20% pour toutes les valeurs — augmenter la plage\n")
  emax_opt <- max(emax_vals)
}

g0_at_opt <- res_list[[as.character(emax_vals[which.min(abs(emax_vals - emax_opt))])]]$g0

cat(sprintf("\n  → Emax_CMP_kill_tdxd_human = %.2f\n", emax_opt))
cat(sprintf("    G0 (Neut ≥ 2.0)    ≈ %.1f%%  (FDA : ~45%%)\n", g0_at_opt))
cat(sprintf("    G3-4 (Neut < 1.0)  ≈ 20%%   (FDA : ~16-20%%)\n"))
cat(sprintf("    ED50_kill = %.2f  ω_Emax = %.2f\n", pars_base$ED50_kill, omega_Emax))
cat("═══════════════════════════════════════════════════════════\n")
