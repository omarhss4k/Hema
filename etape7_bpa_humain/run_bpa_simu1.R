############################################################
# run_bpa_simu1.R
# Simulation de population -- anti-FLT3_chBPA-STINGa20 ADC
# Scenario 1 (sous-traitant) : IC50_MEP=155 nM, IC50_CMP=7.8 nM
#
# PK : placeholder T-DXd (Yin 2020) -- à remplacer par PK BPA
# PD : Fornari 2019 + Slopes calibrés par ratio IC50 vs T-DXd
#
# Résultat attendu :
#   Slope_CMP ≈ 1440  (IC50 BPA 7.8 nM vs IC50 T-DXd ~190 nM)
#   Slope_MEP ≈ 1.19  (IC50 BPA 155 nM vs IC50 T-DXd ~184 nM)
#   → Toxicité myéloïde marquée attendue
#   → Toxicité érythroïde proche T-DXd (MEP similaire)
#
# IIV (log-normal, sans mixture -- pas de données cliniques) :
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
cat(sprintf("\n[Simu 1 -- sous-traitant] IC50_MEP=155 nM  IC50_CMP=7.8 nM\n"))
cat(sprintf("Slope_MEP=%.3f  Slope_CMP=%.1f\n", Slope_MEP_BPA_s1, Slope_CMP_BPA_s1))
cat(sprintf("Simulation : N=%d patients, %.1f mg/kg Q3W × %d cycles\n\n",
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

pb_step <- floor(N_patients / 10)

for (i in 1:N_patients) {

  # Tirage IIV log-normal (sans mixture)
  eta_CL   <- rnorm(1, 0, omega_CL)
  eta_V1   <- rnorm(1, 0, omega_V1)
  eta_SCMP <- rnorm(1, 0, omega_Slope_CMP)
  eta_SMEP <- rnorm(1, 0, omega_Slope_MEP)

  pars_i           <- pars_typ
  pars_i$CL_ADC   <- pars_typ$CL_ADC   * exp(eta_CL)
  pars_i$V1_ADC   <- pars_typ$V1_ADC   * exp(eta_V1)
  pars_i$Slope_CMP <- max(0, Slope_CMP_BPA_s1 * exp(eta_SCMP))
  pars_i$Slope_MEP <- max(0, Slope_MEP_BPA_s1 * exp(eta_SMEP))

  pars_i$rate_fun <- make_tdxd_infusion(
    dose_mgkg  = DOSE_MGKG,
    BW_kg      = BW_KG,
    Tinfu_h    = TINFU_H,
    interval_h = INTERVAL_H,
    n_cycles   = N_CYCLES
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

  if (!is.null(out) && nrow(out) > 10) {
    results$CL_ADC[i]     <- pars_i$CL_ADC
    results$V1_ADC[i]     <- pars_i$V1_ADC
    results$Slope_CMP[i]  <- pars_i$Slope_CMP
    results$Slope_MEP[i]  <- pars_i$Slope_MEP
    results$Cmax_ADC[i]   <- max(out$C_ADC1, na.rm = TRUE)
    results$Damage_max[i] <- max(out$Damage, na.rm = TRUE)
    results$Neut_nadir[i] <- min(out$Neut, na.rm = TRUE)
    results$Ret_nadir[i]  <- min(out$Ret,  na.rm = TRUE)
    results$RBC_nadir[i]  <- min(out$RBC,  na.rm = TRUE)
    results$Plt_nadir[i]  <- min(out$Plt,  na.rm = TRUE)
    results$Grade_Neut[i]  <- ctcae_neut(results$Neut_nadir[i])
    results$Grade_Anemia[i] <- ctcae_anemia(results$RBC_nadir[i], pars_typ$RBC0)
    results$Grade_Plt[i]    <- ctcae_plt(results$Plt_nadir[i])
  }

  if (i %% pb_step == 0)
    cat(sprintf("  %3d/%d patients simulés...\n", i, N_patients))
}

results <- results[!is.na(results$Neut_nadir), ]
n_ok <- nrow(results)
saveRDS(results, "results/bpa_simu1_results.rds")

# ══════════════════════════════════════════════════════════
# Résultats
# ══════════════════════════════════════════════════════════
grade_order <- c("G0", "G1", "G2", "G3", "G4")

tab_neut  <- table(factor(results$Grade_Neut,   levels = grade_order))
tab_anemia<- table(factor(results$Grade_Anemia, levels = grade_order))
tab_plt   <- table(factor(results$Grade_Plt,    levels = grade_order))

pct_n <- round(100 * tab_neut  / n_ok, 1)
pct_a <- round(100 * tab_anemia/ n_ok, 1)
pct_p <- round(100 * tab_plt   / n_ok, 1)

cat("\n═══════════════════════════════════════════════════════════\n")
cat(sprintf("  BPA SIMU 1 (sous-traitant) -- %.1f mg/kg Q3W × %d (N=%d)\n",
            DOSE_MGKG, N_CYCLES, n_ok))
cat(sprintf("  IC50_MEP=155 nM | IC50_CMP=7.8 nM\n"))
cat(sprintf("  Slope_MEP=%.3f | Slope_CMP=%.1f\n", Slope_MEP_BPA_s1, Slope_CMP_BPA_s1))
cat("---------------------------------------------------------\n")
cat("  NEUTROPENIE :\n")
for (g in grade_order)
  cat(sprintf("    %-5s  %6d  %6.1f%%\n", g, tab_neut[g], pct_n[g]))
cat(sprintf("    TOTAL G3-4 : %.1f%%\n", pct_n["G3"] + pct_n["G4"]))
cat(sprintf("    TOUT GRADE : %.1f%%\n",
            pct_n["G1"] + pct_n["G2"] + pct_n["G3"] + pct_n["G4"]))

cat("\n  ANEMIE (proxy RBC) :\n")
for (g in grade_order)
  cat(sprintf("    %-5s  %6d  %6.1f%%\n", g, tab_anemia[g], pct_a[g]))
cat(sprintf("    TOTAL G3-4 : %.1f%%\n", pct_a["G3"] + pct_a["G4"]))

cat("\n  THROMBOCYTOPENIE :\n")
for (g in grade_order)
  cat(sprintf("    %-5s  %6d  %6.1f%%\n", g, tab_plt[g], pct_p[g]))
cat(sprintf("    TOTAL G3-4 : %.1f%%\n", pct_p["G3"] + pct_p["G4"]))

cat("\n---------------------------------------------------------\n")
quants <- function(x) quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE)
q_adc  <- quants(results$Cmax_ADC)
q_neut <- quants(results$Neut_nadir)
q_rbc  <- quants(results$RBC_nadir / pars_typ$RBC0)
cat("  Exposition médiane [P10-P90] :\n")
cat(sprintf("    ADC Cmax   : %.0f [%.0f-%.0f] µg/mL\n", q_adc[2], q_adc[1], q_adc[3]))
cat(sprintf("    Neut nadir : %.2f [%.2f-%.2f] × 10⁹/L\n", q_neut[2], q_neut[1], q_neut[3]))
cat(sprintf("    RBC/RBC0   : %.2f [%.2f-%.2f]\n", q_rbc[2], q_rbc[1], q_rbc[3]))
cat("═══════════════════════════════════════════════════════════\n")

# ══════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════
pdf("results/BPA_simu1_VPC.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1))

grade_cols <- c(G0 = "#2166ac", G1 = "#74add1", G2 = "#f4a582",
                G3 = "#d6604d", G4 = "#b2182b")

# 1. Distribution nadir Neut
hist(results$Neut_nadir, breaks = 30, col = "#74add1", border = "white",
     main = "BPA Simu1 : Nadir Neutrophiles",
     xlab = "Neut nadir [× 10⁹/L]", ylab = "Frequence")
abline(v = c(2.0, 1.5, 1.0, 0.5), lty = 2,
       col = c("grey60","orange","red","darkred"), lwd = 1.5)
text(c(2.0, 1.5, 1.0, 0.5) + 0.05, par("usr")[4] * 0.95,
     c("G1","G2","G3","G4"),
     col = c("grey60","orange","red","darkred"), cex = 0.85, adj = 0)
abline(v = median(results$Neut_nadir, na.rm = TRUE), col = "#2166ac", lwd = 2)

# 2. Barplot grades Neutropenie
bp <- barplot(pct_n[grade_order],
              col = grade_cols[grade_order],
              names.arg = grade_order,
              ylim = c(0, 100),
              main = "Grades Neutropenie BPA Simu1",
              ylab = "Patients (%)")
text(bp, pct_n[grade_order] + 2,
     sprintf("%.0f%%", pct_n[grade_order]), cex = 0.85)

# 3. Distribution RBC nadir (anemie)
pct_drop_rbc <- 100 * (results$RBC_nadir - pars_typ$RBC0) / pars_typ$RBC0
hist(pct_drop_rbc, breaks = 30, col = "#f4a582", border = "white",
     main = "BPA Simu1 : Chute RBC (%Δ baseline)",
     xlab = "%Δ RBC nadir", ylab = "Frequence")
abline(v = median(pct_drop_rbc, na.rm = TRUE), col = "#d6604d", lwd = 2)

# 4. Barplot grades Anemie
bp <- barplot(pct_a[grade_order],
              col = grade_cols[grade_order],
              names.arg = grade_order,
              ylim = c(0, 100),
              main = "Grades Anemie BPA Simu1",
              ylab = "Patients (%)")
text(bp, pct_a[grade_order] + 2,
     sprintf("%.0f%%", pct_a[grade_order]), cex = 0.85)

# 5. Neut nadir vs Slope_CMP
plot(results$Slope_CMP, results$Neut_nadir,
     pch = 16, cex = 0.6,
     col = grade_cols[results$Grade_Neut],
     xlab = "Slope_CMP (sensibilite PD)",
     ylab = "Neut nadir [× 10⁹/L]",
     main = "Neut nadir vs Slope_CMP")
abline(h = c(2.0, 1.5, 1.0, 0.5), lty = 2,
       col = c("grey60","orange","red","darkred"))
legend("topright", names(grade_cols), pch = 16,
       col = grade_cols, bty = "n", cex = 0.8)

# 6. Neut nadir vs CL_ADC
plot(results$CL_ADC * 24, results$Neut_nadir,
     pch = 16, cex = 0.6,
     col = grade_cols[results$Grade_Neut],
     xlab = "CL_ADC [L/jour]",
     ylab = "Neut nadir [× 10⁹/L]",
     main = "Neut nadir vs CL_ADC")
abline(h = c(2.0, 1.5, 1.0, 0.5), lty = 2,
       col = c("grey60","orange","red","darkred"))

dev.off()
cat("  -> results/BPA_simu1_VPC.pdf\n")
