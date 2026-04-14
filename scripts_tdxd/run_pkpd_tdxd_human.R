############################################################
# run_pkpd_tdxd_human.R
# Simulation PK/PD T-DXd — HUMAIN
#
# Scénarios :
#   1. 5.4 mg/kg Q3W × 6 cycles  (dose approuvée FDA)
#   2. 6.4 mg/kg Q3W × 6 cycles  (dose supérieure, exploration)
#
# Validation PK : FDA BLA 761139 (Cmax ADC, Cmax DXd)
# Validation PD : clinique (Grade 1-2 Neut en médiane)
############################################################
library(deSolve)

# setwd() supprimé — utiliser Rscript depuis scripts_tdxd/
# ou lancer via : source("scripts_tdxd/run_pkpd_tdxd_human.R")

source("parameters_tdxd_rat.R")     # pour make_tdxd_infusion
source("parameters_tdxd_human.R")   # PK humain + baselines Fornari humain
source("pkpd_tdxd_rat.R")           # ODE fusionné (même modèle, params différents)

dir.create("results_PKPD_human", showWarnings = FALSE)

# ── Fusion paramètres PK + PD ────────────────────────────
pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL

pars_full_hu <- c(pars_pd_hu, tdxd_pars_hu)

# ── Override Slopes calibrés T-DXd (DESTINY-Breast01) ────
# Slope_CMP : calibré à 12.0 pour G3-4 neutropénie ~20%
pars_full_hu$Slope_CMP <- Slope_CMP_tdxd_human
# Slope_MEP : à calibrer pour G3-4 anémie ~9%
#   Tant que Slope_MEP_tdxd_human n'est pas calibré, on conserve
#   la valeur carboplatin humain de init_pars (comportement précédent)
if (!is.na(Slope_MEP_tdxd_human)) {
  pars_full_hu$Slope_MEP <- Slope_MEP_tdxd_human
} else {
  warning("Slope_MEP_tdxd_human non calibré — valeur carboplatin utilisée")
}

# ── État initial complet ─────────────────────────────────
state_pd_hu  <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu    <- c(tdxd_hu_state0, state_pd_hu)

# ── Grille temporelle — 6 cycles Q3W (126 jours) ─────────
times_hu <- seq(0, 126 * 24, by = 4)   # pas de 4h

# ══════════════════════════════════════════════════════════
# Scénario 1 — 5.4 mg/kg Q3W × 6 cycles (dose approuvée)
# ══════════════════════════════════════════════════════════
cat("=== Scénario 1 : T-DXd 5.4 mg/kg Q3W × 6 cycles (dose approuvée) ===\n")

pars_s1_hu          <- pars_full_hu
pars_s1_hu$rate_fun <- make_tdxd_infusion(
  dose_mgkg  = 5.4, BW_kg = 70,
  Tinfu_h    = 1.5,
  interval_h = 21 * 24,
  n_cycles   = 6
)

sim_hu1 <- simulate_pkpd_tdxd(times_hu, pars_s1_hu, state0_hu)

# ── Validation PK ────────────────────────────────────────
Cmax_ADC_sim  <- max(sim_hu1$C_ADC1,  na.rm = TRUE)
Cmax_DXd_sim  <- max(sim_hu1$C_DXd,   na.rm = TRUE) * 1000  # mg/L → ng/mL

cat(sprintf("\n  Validation PK vs FDA BLA 761139 :\n"))
cat(sprintf("  ADC Cmax  : sim=%.1f  FDA=%.1f µg/mL  (%+.1f%%)\n",
            Cmax_ADC_sim, tdxd_hu_targets$Cmax_ADC_mgL,
            100 * (Cmax_ADC_sim - tdxd_hu_targets$Cmax_ADC_mgL) / tdxd_hu_targets$Cmax_ADC_mgL))
cat(sprintf("  DXd Cmax  : sim=%.2f  FDA=%.2f ng/mL  (%+.1f%%)\n",
            Cmax_DXd_sim, tdxd_hu_targets$Cmax_DXd_ngmL,
            100 * (Cmax_DXd_sim - tdxd_hu_targets$Cmax_DXd_ngmL) / tdxd_hu_targets$Cmax_DXd_ngmL))

# ══════════════════════════════════════════════════════════
# Scénario 2 — 6.4 mg/kg Q3W × 6 cycles
# ══════════════════════════════════════════════════════════
cat("\n=== Scénario 2 : T-DXd 6.4 mg/kg Q3W × 6 cycles ===\n")

pars_s2_hu          <- pars_full_hu
pars_s2_hu$rate_fun <- make_tdxd_infusion(
  dose_mgkg  = 6.4, BW_kg = 70,
  Tinfu_h    = 1.5,
  interval_h = 21 * 24,
  n_cycles   = 6
)

sim_hu2 <- simulate_pkpd_tdxd(times_hu, pars_s2_hu, state0_hu)

# ══════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════
pdf("results_PKPD_human/PKPD_human_5mg4_6mg4_Q3Wx6.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1))

times_d <- sim_hu1$time_d
col_54  <- "#2166ac"
col_64  <- "#b2182b"

# ── 1. ADC sérum ──
plot(times_d, sim_hu1$C_ADC1, type = "l", col = col_54, lwd = 2,
     xlab = "Temps (jours)", ylab = "ADC [µg/mL]",
     main = "ADC sérum — 2 compartiments")
lines(times_d, sim_hu2$C_ADC1, col = col_64, lwd = 2)
abline(h = tdxd_hu_targets$Cmax_ADC_mgL, lty = 2, col = "grey40")
text(5, tdxd_hu_targets$Cmax_ADC_mgL * 1.05, "FDA Cmax 122 µg/mL", cex = 0.8, col = "grey40")
abline(v = seq(0, 125*24, by = 21*24)/24, lty = 3, col = "grey80")
legend("topright", c("5.4 mg/kg", "6.4 mg/kg"), col = c(col_54, col_64),
       lwd = 2, bty = "n", cex = 0.9)

# ── 2. DXd plasma ──
plot(times_d, sim_hu1$C_DXd * 1000, type = "l", col = col_54, lwd = 2,
     xlab = "Temps (jours)", ylab = "DXd plasma [ng/mL]",
     main = "DXd plasma libre")
lines(times_d, sim_hu2$C_DXd * 1000, col = col_64, lwd = 2)
abline(h = tdxd_hu_targets$Cmax_DXd_ngmL, lty = 2, col = "grey40")
text(5, tdxd_hu_targets$Cmax_DXd_ngmL * 1.1, "FDA Cmax 4.4 ng/mL", cex = 0.8, col = "grey40")
abline(v = seq(0, 125*24, by = 21*24)/24, lty = 3, col = "grey80")
legend("topright", c("5.4 mg/kg", "6.4 mg/kg"), col = c(col_54, col_64),
       lwd = 2, bty = "n", cex = 0.9)

# ── 3. Dommages ADN ──
plot(times_d, sim_hu1$Damage, type = "l", col = col_54, lwd = 2,
     xlab = "Temps (jours)", ylab = "Damage [normalisé]",
     main = "Dommages ADN (γH2AX)")
lines(times_d, sim_hu2$Damage, col = col_64, lwd = 2)
abline(v = seq(0, 125*24, by = 21*24)/24, lty = 3, col = "grey80")
legend("topright", c("5.4 mg/kg", "6.4 mg/kg"), col = c(col_54, col_64),
       lwd = 2, bty = "n", cex = 0.9)

# ── 4. Neutrophiles ──
Neut0_hu <- pars_full_hu$Neut0
plot(times_d, sim_hu1$Neut, type = "l", col = col_54, lwd = 2,
     ylim = c(min(sim_hu2$Neut, na.rm=TRUE) * 0.9, Neut0_hu * 1.1),
     xlab = "Temps (jours)", ylab = "Neut [×10⁹/L]",
     main = "Neutrophiles")
lines(times_d, sim_hu2$Neut, col = col_64, lwd = 2)
abline(h = Neut0_hu, lty = 3, col = "grey40")
abline(h = 2.0, lty = 2, col = "orange", lwd = 1.5)   # Grade 1
abline(h = 1.0, lty = 2, col = "red",    lwd = 1.5)   # Grade 3
text(130, 2.1, "G1 (<2.0)", col = "orange", cex = 0.8)
text(130, 1.1, "G3 (<1.0)", col = "red",    cex = 0.8)
abline(v = seq(0, 125*24, by = 21*24)/24, lty = 3, col = "grey80")
legend("bottomright", c("5.4 mg/kg", "6.4 mg/kg"), col = c(col_54, col_64),
       lwd = 2, bty = "n", cex = 0.9)

# ── 5. Réticulocytes ──
Ret0_hu <- pars_full_hu$Ret0
plot(times_d, sim_hu1$Ret, type = "l", col = col_54, lwd = 2,
     ylim = c(min(sim_hu2$Ret, na.rm=TRUE) * 0.9, Ret0_hu * 1.1),
     xlab = "Temps (jours)", ylab = "Ret [×10⁹/L]",
     main = "Réticulocytes")
lines(times_d, sim_hu2$Ret, col = col_64, lwd = 2)
abline(h = Ret0_hu, lty = 3, col = "grey40")
abline(v = seq(0, 125*24, by = 21*24)/24, lty = 3, col = "grey80")
legend("bottomright", c("5.4 mg/kg", "6.4 mg/kg"), col = c(col_54, col_64),
       lwd = 2, bty = "n", cex = 0.9)

# ── 6. Plaquettes ──
Plt0_hu <- pars_full_hu$Plt0
plot(times_d, sim_hu1$Plt, type = "l", col = col_54, lwd = 2,
     ylim = c(min(sim_hu2$Plt, na.rm=TRUE) * 0.9, Plt0_hu * 1.15),
     xlab = "Temps (jours)", ylab = "Plt [×10⁹/L]",
     main = "Plaquettes")
lines(times_d, sim_hu2$Plt, col = col_64, lwd = 2)
abline(h = Plt0_hu, lty = 3, col = "grey40")
abline(h = 75, lty = 2, col = "red", lwd = 1.5)   # Grade 3 thrombopénie
text(130, 80, "G3 (<75)", col = "red", cex = 0.8)
abline(v = seq(0, 125*24, by = 21*24)/24, lty = 3, col = "grey80")
legend("bottomright", c("5.4 mg/kg", "6.4 mg/kg"), col = c(col_54, col_64),
       lwd = 2, bty = "n", cex = 0.9)

dev.off()
cat("  -> results_PKPD_human/PKPD_human_5mg4_6mg4_Q3Wx6.pdf\n")

# ══════════════════════════════════════════════════════════
# Résumé
# ══════════════════════════════════════════════════════════
cat("\n═══════════════════════════════════════════════════════════\n")
cat("  RÉSUMÉ PKPD — T-DXd HUMAIN (6 cycles Q3W)\n")
cat("─────────────────────────────────────────────────────────\n")

for (dose_str in c("5.4", "6.4")) {
  sim_i <- if (dose_str == "5.4") sim_hu1 else sim_hu2
  cat(sprintf("\n  Dose : %s mg/kg\n", dose_str))
  cat(sprintf("    ADC Cmax    = %.1f µg/mL      (FDA = 122)\n",
              max(sim_i$C_ADC1, na.rm=TRUE)))
  cat(sprintf("    DXd Cmax    = %.2f ng/mL      (FDA = 4.4)\n",
              max(sim_i$C_DXd, na.rm=TRUE) * 1000))
  cat(sprintf("    Damage max  = %.4f\n",
              max(sim_i$Damage, na.rm=TRUE)))
  cat(sprintf("    Neut nadir  = %.2f × 10⁹/L  (%+.1f%%)\n",
              min(sim_i$Neut, na.rm=TRUE),
              100*(min(sim_i$Neut,na.rm=TRUE) - Neut0_hu)/Neut0_hu))
  cat(sprintf("    Ret  nadir  = %.1f × 10⁹/L   (%+.1f%%)\n",
              min(sim_i$Ret,  na.rm=TRUE),
              100*(min(sim_i$Ret, na.rm=TRUE) - Ret0_hu)/Ret0_hu))
  cat(sprintf("    Plt  nadir  = %.1f × 10⁹/L   (%+.1f%%)\n",
              min(sim_i$Plt,  na.rm=TRUE),
              100*(min(sim_i$Plt, na.rm=TRUE) - Plt0_hu)/Plt0_hu))

  neut_min <- min(sim_i$Neut, na.rm=TRUE)
  grade <- if (neut_min < 0.5) "Grade 4" else if (neut_min < 1.0) "Grade 3" else
           if (neut_min < 2.0) "Grade 2/1" else "pas de neutropénie"
  cat(sprintf("    → Neutropénie : %s\n", grade))
}

cat("\n─────────────────────────────────────────────────────────\n")
cat("  Référence FDA (DESTINY-Breast01, n=184) :\n")
cat("  Grade 3-4 neutropénie : ~20% patients à 5.4 mg/kg\n")
cat("  Grade 1-2 neutropénie : ~35% patients\n")
cat("═══════════════════════════════════════════════════════════\n")
