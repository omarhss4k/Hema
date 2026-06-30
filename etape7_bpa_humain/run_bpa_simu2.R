############################################################
# run_bpa_simu2.R
# Simulation de population -- anti-FLT3_chBPA-STINGa20 ADC
# Scenario 2 (interne) : IC50_MEP=155 nM, IC50_CMP=0.015 nM
#
# PK : placeholder T-DXd (Yin 2020) -- à remplacer par PK BPA
# PD : Fornari 2019 + Slopes calibrés par ratio IC50 vs T-DXd
#
# ATTENTION :
#   IC50_CMP=0.015 nM → Slope_CMP ≈ 748 000 (EXTREME)
#   Toxicité myéloïde quasi-complète attendue
#   Résultat = borne supérieure théorique
#
# IIV (log-normal, sans mixture) :
#   CL_ADC  : ω = 0.35 (Yin 2020)
#   V1_ADC  : ω = 0.20 (Yin 2020)
#   Slope   : ω = 0.33 (Fornari Table S4)
############################################################
library(deSolve)

source("../etape2_carboplatin_humain/parameters_human.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("../etape5_tdxd_humain/parameters_tdxd_human.R")
source("../etape3_tdxd_rat/parameters_tdxd_rat.R")
source("../etape3_tdxd_rat/pkpd_tdxd_rat.R")
source("parameters_bpa_human.R")

dir.create("results", showWarnings = FALSE)

cat(sprintf("\n!!! AVERTISSEMENT : Slope_CMP=%.0f (IC50=0.015 nM interne)\n",
            Slope_CMP_BPA_s2))
cat("    Toxicite myeloide quasi-complete attendue.\n")
cat("    Ce scenario est une borne superieure theorique.\n\n")

set.seed(42)
N_patients <- 300

# -- Paramètres typiques BPA ---------------------------------
pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL
pars_typ <- c(pars_pd_hu, bpa_pars)

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)

times_hu <- seq(0, 126 * 24, by = 6)
n_times  <- length(times_hu)
times_j  <- times_hu / 24

# -- Grilles CTCAE -----------------------------------------
ctcae_neut <- function(neut) {
  if      (neut < 0.5) "G4"
  else if (neut < 1.0) "G3"
  else if (neut < 1.5) "G2"
  else if (neut < 2.0) "G1"
  else                  "G0"
}
ctcae_anemia <- function(rbc, rbc0) {
  frac <- rbc / rbc0
  if      (frac < 0.54) "G4"
  else if (frac < 0.67) "G3"
  else if (frac < 0.80) "G2"
  else if (frac < 0.90) "G1"
  else                   "G0"
}
ctcae_plt <- function(plt) {
  if      (plt < 25)  "G4"
  else if (plt < 50)  "G3"
  else if (plt < 75)  "G2"
  else if (plt < 150) "G1"
  else                "G0"
}

# ══════════════════════════════════════════════════════════
# Simulation population
# ══════════════════════════════════════════════════════════
cat(sprintf("[Simu 2 -- interne] IC50_MEP=155 nM  IC50_CMP=0.015 nM\n"))
cat(sprintf("Slope_MEP=%.3f  Slope_CMP=%.0f\n", Slope_MEP_BPA_s2, Slope_CMP_BPA_s2))
cat(sprintf("Simulation : N=%d patients, %.1f mg/kg Q3W x %d cycles\n\n",
            N_patients, DOSE_MGKG, N_CYCLES))
cat("(patience ~3-5 min)\n\n")

results <- data.frame(
  id           = 1:N_patients,
  CL_ADC       = NA_real_,
  V1_ADC       = NA_real_,
  Slope_CMP    = NA_real_,
  Slope_MEP    = NA_real_,
  Cmax_ADC     = NA_real_,
  Damage_max   = NA_real_,
  Neut_nadir   = NA_real_,
  Ret_nadir    = NA_real_,
  RBC_nadir    = NA_real_,
  Plt_nadir    = NA_real_,
  Grade_Neut   = NA_character_,
  Grade_Anemia = NA_character_,
  Grade_Plt    = NA_character_
)

mat_neut <- matrix(NA_real_, nrow = n_times, ncol = N_patients)
mat_rbc  <- matrix(NA_real_, nrow = n_times, ncol = N_patients)
mat_plt  <- matrix(NA_real_, nrow = n_times, ncol = N_patients)
mat_adc  <- matrix(NA_real_, nrow = n_times, ncol = N_patients)

pb_step <- floor(N_patients / 10)

for (i in 1:N_patients) {

  eta_CL   <- rnorm(1, 0, omega_CL)
  eta_V1   <- rnorm(1, 0, omega_V1)
  eta_SCMP <- rnorm(1, 0, omega_Slope_CMP)
  eta_SMEP <- rnorm(1, 0, omega_Slope_MEP)

  pars_i            <- pars_typ
  pars_i$CL_ADC    <- pars_typ$CL_ADC   * exp(eta_CL)
  pars_i$V1_ADC    <- pars_typ$V1_ADC   * exp(eta_V1)
  pars_i$Slope_CMP <- max(0, Slope_CMP_BPA_s2 * exp(eta_SCMP))
  pars_i$Slope_MEP <- max(0, Slope_MEP_BPA_s2 * exp(eta_SMEP))
  pars_i$rate_fun  <- make_tdxd_infusion(
    dose_mgkg  = DOSE_MGKG, BW_kg = BW_KG,
    Tinfu_h    = TINFU_H,   interval_h = INTERVAL_H, n_cycles = N_CYCLES
  )

  out <- tryCatch(
    suppressMessages(suppressWarnings(
      as.data.frame(lsoda(
        y = state0_hu, times = times_hu,
        func = pkpd_tdxd_fornari, parms = pars_i,
        rtol = 1e-5, atol = 1e-7, maxsteps = 500000
      ))
    )),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) == n_times) {
    results$CL_ADC[i]       <- pars_i$CL_ADC
    results$V1_ADC[i]       <- pars_i$V1_ADC
    results$Slope_CMP[i]    <- pars_i$Slope_CMP
    results$Slope_MEP[i]    <- pars_i$Slope_MEP
    results$Cmax_ADC[i]     <- max(out$C_ADC1, na.rm = TRUE)
    results$Damage_max[i]   <- max(out$Damage, na.rm = TRUE)
    results$Neut_nadir[i]   <- min(out$Neut,   na.rm = TRUE)
    results$Ret_nadir[i]    <- min(out$Ret,    na.rm = TRUE)
    results$RBC_nadir[i]    <- min(out$RBC,    na.rm = TRUE)
    results$Plt_nadir[i]    <- min(out$Plt,    na.rm = TRUE)
    results$Grade_Neut[i]   <- ctcae_neut(results$Neut_nadir[i])
    results$Grade_Anemia[i] <- ctcae_anemia(results$RBC_nadir[i], pars_typ$RBC0)
    results$Grade_Plt[i]    <- ctcae_plt(results$Plt_nadir[i])

    mat_neut[, i] <- out$Neut
    mat_rbc[, i]  <- out$RBC
    mat_plt[, i]  <- out$Plt
    mat_adc[, i]  <- out$C_ADC1
  }

  if (i %% pb_step == 0)
    cat(sprintf("  %3d/%d patients simules...\n", i, N_patients))
}

results <- results[!is.na(results$Neut_nadir), ]
n_ok    <- nrow(results)
ok_cols <- !is.na(mat_neut[1, ])
mat_neut <- mat_neut[, ok_cols, drop = FALSE]
mat_rbc  <- mat_rbc[,  ok_cols, drop = FALSE]
mat_plt  <- mat_plt[,  ok_cols, drop = FALSE]
mat_adc  <- mat_adc[,  ok_cols, drop = FALSE]

saveRDS(results, "results/bpa_simu2_results.rds")

# ══════════════════════════════════════════════════════════
# Résultats console
# ══════════════════════════════════════════════════════════
grade_order <- c("G0", "G1", "G2", "G3", "G4")
tab_neut   <- table(factor(results$Grade_Neut,   levels = grade_order))
tab_anemia <- table(factor(results$Grade_Anemia, levels = grade_order))
tab_plt    <- table(factor(results$Grade_Plt,    levels = grade_order))
pct_n <- round(100 * tab_neut  / n_ok, 1)
pct_a <- round(100 * tab_anemia/ n_ok, 1)
pct_p <- round(100 * tab_plt   / n_ok, 1)

cat("\n=============================================================\n")
cat(sprintf("  BPA SIMU 2 (interne) -- %.1f mg/kg Q3W x %d (N=%d)\n",
            DOSE_MGKG, N_CYCLES, n_ok))
cat(sprintf("  IC50_MEP=155 nM | IC50_CMP=0.015 nM  [EXTREME]\n"))
cat(sprintf("  Slope_MEP=%.3f | Slope_CMP=%.0f\n", Slope_MEP_BPA_s2, Slope_CMP_BPA_s2))
cat("-------------------------------------------------------------\n")
cat("  NEUTROPENIE :\n")
for (g in grade_order) cat(sprintf("    %-5s  %6d  %6.1f%%\n", g, tab_neut[g], pct_n[g]))
cat(sprintf("    TOTAL G3-4 : %.1f%%\n", pct_n["G3"] + pct_n["G4"]))

cat("\n  ANEMIE (proxy RBC) :\n")
for (g in grade_order) cat(sprintf("    %-5s  %6d  %6.1f%%\n", g, tab_anemia[g], pct_a[g]))
cat(sprintf("    TOTAL G3-4 : %.1f%%\n", pct_a["G3"] + pct_a["G4"]))

cat("\n  THROMBOCYTOPENIE :\n")
for (g in grade_order) cat(sprintf("    %-5s  %6d  %6.1f%%\n", g, tab_plt[g], pct_p[g]))
cat(sprintf("    TOTAL G3-4 : %.1f%%\n", pct_p["G3"] + pct_p["G4"]))

quants <- function(x) quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE)
cat("\n  Exposition mediane [P10-P90] :\n")
q_adc  <- quants(results$Cmax_ADC)
q_neut <- quants(results$Neut_nadir)
q_rbc  <- quants(results$RBC_nadir / pars_typ$RBC0)
cat(sprintf("    ADC Cmax   : %.0f [%.0f-%.0f] ug/mL\n", q_adc[2], q_adc[1], q_adc[3]))
cat(sprintf("    Neut nadir : %.3f [%.3f-%.3f] x 10^9/L\n", q_neut[2], q_neut[1], q_neut[3]))
cat(sprintf("    RBC/RBC0   : %.2f [%.2f-%.2f]\n", q_rbc[2], q_rbc[1], q_rbc[3]))
cat("\n  NOTE : Toxicite myeloide quasi-complete.\n")
cat("         Confirmer IC50_CMP=0.015 nM avant interpretation.\n")
cat("=============================================================\n")

# ══════════════════════════════════════════════════════════
# Helper : bande VPC
# ══════════════════════════════════════════════════════════
.vpc_band <- function(mat, times_j, col_med, col_band,
                      ylab, main, ylim = NULL,
                      hlines = NULL, hcols = NULL, hlabs = NULL) {
  p10 <- apply(mat, 1, quantile, 0.10, na.rm = TRUE)
  p50 <- apply(mat, 1, quantile, 0.50, na.rm = TRUE)
  p90 <- apply(mat, 1, quantile, 0.90, na.rm = TRUE)
  if (is.null(ylim)) ylim <- c(0, max(p90, na.rm = TRUE) * 1.05)

  plot(times_j, p50, type = "n", ylim = ylim,
       xlab = "Temps (jours)", ylab = ylab, main = main)
  polygon(c(times_j, rev(times_j)), c(p10, rev(p90)),
          col = col_band, border = NA)
  lines(times_j, p50, col = col_med, lwd = 2)

  if (!is.null(hlines)) {
    abline(h = hlines, lty = 2, col = hcols, lwd = 1.2)
    if (!is.null(hlabs))
      text(rep(max(times_j) * 0.98, length(hlines)), hlines + ylim[2] * 0.02,
           hlabs, col = hcols, cex = 0.75, adj = 1)
  }
  cycle_days <- (0:(N_CYCLES - 1)) * (INTERVAL_H / 24)
  abline(v = cycle_days, lty = 3, col = "grey70")
  legend("topright", legend = c("Mediane", "P10-P90"),
         col = c(col_med, col_band), lwd = c(2, 8), bty = "n", cex = 0.8)
}

# ══════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════
grade_cols <- c(G0 = "#2166ac", G1 = "#74add1", G2 = "#f4a582",
                G3 = "#d6604d", G4 = "#b2182b")

pdf("results/BPA_simu2_VPC.pdf", width = 14, height = 10)

# ── Page 1 : profils temporels ───────────────────────────
par(mfrow = c(2, 2), mar = c(4, 4.5, 3.5, 1.5))

.vpc_band(mat_adc, times_j,
          col_med  = "#1a1a2e",
          col_band = adjustcolor("#1a1a2e", alpha.f = 0.20),
          ylab = "ADC [ug/mL]",
          main = sprintf("BPA Simu2 : PK ADC (%.1f mg/kg Q3W x%d)", DOSE_MGKG, N_CYCLES))

.vpc_band(mat_neut, times_j,
          col_med  = "#b2182b",
          col_band = adjustcolor("#b2182b", alpha.f = 0.20),
          ylab = "Neutrophiles [x 10^9/L]",
          main = "BPA Simu2 : Neutrophiles [EXTREME]",
          hlines = c(2.0, 1.0, 0.5),
          hcols  = c("grey50", "red", "darkred"),
          hlabs  = c("G1", "G3", "G4"))
abline(h = pars_typ$Neut0, lty = 1, col = "grey30", lwd = 1)

rbc0 <- pars_typ$RBC0
.vpc_band(mat_rbc, times_j,
          col_med  = "#d6604d",
          col_band = adjustcolor("#d6604d", alpha.f = 0.20),
          ylab = "RBC [x 10^9/L]",
          main = "BPA Simu2 : Erythrocytes",
          hlines = c(rbc0 * 0.67, rbc0 * 0.80),
          hcols  = c("darkred", "orange"),
          hlabs  = c("G3", "G2"))
abline(h = rbc0, lty = 1, col = "grey30", lwd = 1)

.vpc_band(mat_plt, times_j,
          col_med  = "#4dac26",
          col_band = adjustcolor("#4dac26", alpha.f = 0.20),
          ylab = "Plaquettes [x 10^9/L]",
          main = "BPA Simu2 : Plaquettes",
          hlines = c(150, 75, 50, 25),
          hcols  = c("grey50", "orange", "red", "darkred"),
          hlabs  = c("G1", "G2", "G3", "G4"))
abline(h = pars_typ$Plt0, lty = 1, col = "grey30", lwd = 1)

# ── Page 2 : grades + comparaison Simu1 vs Simu2 ─────────
par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1))

# Grades Neut
bp <- barplot(pct_n[grade_order], col = grade_cols[grade_order],
              names.arg = grade_order, ylim = c(0, 100),
              main = "Grades Neutropenie BPA Simu2", ylab = "Patients (%)")
text(bp, pct_n[grade_order] + 2, sprintf("%.0f%%", pct_n[grade_order]), cex = 0.85)

# Grades Anemie
bp <- barplot(pct_a[grade_order], col = grade_cols[grade_order],
              names.arg = grade_order, ylim = c(0, 100),
              main = "Grades Anemie BPA Simu2", ylab = "Patients (%)")
text(bp, pct_a[grade_order] + 2, sprintf("%.0f%%", pct_a[grade_order]), cex = 0.85)

# Grades Thrombocytopenie
bp <- barplot(pct_p[grade_order], col = grade_cols[grade_order],
              names.arg = grade_order, ylim = c(0, 100),
              main = "Grades Thrombocytopenie BPA Simu2", ylab = "Patients (%)")
text(bp, pct_p[grade_order] + 2, sprintf("%.0f%%", pct_p[grade_order]), cex = 0.85)

# Comparaison Simu1 vs Simu2 (Neut)
simu1_path <- "results/bpa_simu1_results.rds"
if (file.exists(simu1_path)) {
  res1  <- readRDS(simu1_path)
  pct1  <- round(100 * table(factor(res1$Grade_Neut, levels=grade_order)) / nrow(res1), 1)
  cmp   <- rbind(pct1[grade_order], pct_n[grade_order])
  bp    <- barplot(cmp, beside = TRUE, col = c("#2166ac","#b2182b"),
                   names.arg = grade_order, ylim = c(0, 105),
                   main = "Neutropenie : Simu1 vs Simu2", ylab = "Patients (%)")
  legend("topright", legend = c("Simu1 (sous-traitant)","Simu2 (interne)"),
         fill = c("#2166ac","#b2182b"), bty = "n", cex = 0.85)
} else {
  plot.new()
  text(0.5, 0.5, "Lancer run_bpa_simu1.R\npour la comparaison", cex = 1.2)
}

# Nadir Neut distribution
hist(results$Neut_nadir, breaks = 30, col = "#d6604d", border = "white",
     main = "Distribution nadir Neutrophiles (Simu2)",
     xlab = "Neut nadir [x 10^9/L]", ylab = "Frequence")
abline(v = c(2.0, 1.0, 0.5), lty = 2, col = c("grey60","red","darkred"), lwd = 1.5)

# Nadir RBC distribution
pct_drop_rbc <- 100 * (results$RBC_nadir - pars_typ$RBC0) / pars_typ$RBC0
hist(pct_drop_rbc, breaks = 30, col = "#f4a582", border = "white",
     main = "Chute RBC (%Delta baseline, Simu2)",
     xlab = "%Delta RBC nadir", ylab = "Frequence")
abline(v = median(pct_drop_rbc, na.rm = TRUE), col = "#d6604d", lwd = 2)

dev.off()
cat("  -> results/BPA_simu2_VPC.pdf\n")
