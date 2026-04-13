############################################################
# run_pkpd_tdxd_human_population.R
# Simulation de population — T-DXd 5.4 mg/kg Q3W × 6 cycles
#
# Objectif : reproduire la distribution de toxicité CTCAE
#   DESTINY-Breast01 (FDA BLA 761139, n=184) :
#   Neutropénie G3-4 : ~20%    Anémie G3-4 : ~9%
#   Neutropénie tout grade : ~29%  Anémie tout grade : ~70%
#
# IIV (log-normal) :
#   CL_ADC     : ω = 0.35 (Yin 2020 PopPK)
#   V1_ADC     : ω = 0.20 (Yin 2020 PopPK)
#   Slope_CMP  : ω = 0.33 (Fornari 2019 Table S4)
#   Slope_MEP  : ω = 0.33 (Fornari 2019 Table S4)
############################################################
library(deSolve)

setwd("/home/user/Hema/scripts_tdxd")
source("parameters_tdxd_rat.R")
source("parameters_tdxd_human.R")
source("pkpd_tdxd_rat.R")

dir.create("results_PKPD_human", showWarnings = FALSE)

set.seed(42)
N_patients <- 300

# ── IIV (ω sur log-normal) ───────────────────────────────
omega_CL    <- 0.35   # Yin 2020
omega_V1    <- 0.20   # Yin 2020
omega_Slope_CMP <- 0.33   # Fornari Table S4
omega_Slope_MEP <- 0.33   # Fornari Table S4

# ── Paramètres typiques ──────────────────────────────────
pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL
pars_typ   <- c(pars_pd_hu, tdxd_pars_hu)

# ── Override Slope_CMP calibré T-DXd (DESTINY-Breast01) ──
pars_typ$Slope_CMP <- Slope_CMP_tdxd_human

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)

times_hu <- seq(0, 126 * 24, by = 6)   # pas 6h (rapidité)

# ── CTCAE Grade Neutropénie ──────────────────────────────
ctcae_neut <- function(neut) {
  if      (neut < 0.5) "G4"
  else if (neut < 1.0) "G3"
  else if (neut < 1.5) "G2"
  else if (neut < 2.0) "G1"
  else                  "G0"
}

# ── CTCAE Grade Anémie (proxy Hgb via RBC) ───────────────
# RBC0 = 5000 × 10⁹/L → Hgb ≈ 12 g/dL
# Grade basé sur % drop de RBC (proxy Hgb) :
#   G1 : Hgb < LNL (≈10 g/dL → RBC -17%)
#   G2 : Hgb < 10.0 → RBC < 83% baseline
#   G3 : Hgb < 8.0  → RBC < 67% baseline
#   G4 : Hgb < 6.5  → RBC < 54% baseline
ctcae_anemia <- function(rbc, rbc0) {
  frac <- rbc / rbc0
  if      (frac < 0.54) "G4"
  else if (frac < 0.67) "G3"
  else if (frac < 0.83) "G2"
  else if (frac < 0.90) "G1"
  else                   "G0"
}

# ══════════════════════════════════════════════════════════
# Simulation population
# ══════════════════════════════════════════════════════════
cat(sprintf("Simulation population : N=%d patients, 5.4 mg/kg Q3W × 6 cycles\n", N_patients))
cat("(patience ~3-5 min)\n\n")

results <- data.frame(
  id          = 1:N_patients,
  CL_ADC      = NA_real_,
  V1_ADC      = NA_real_,
  Slope_CMP   = NA_real_,
  Slope_MEP   = NA_real_,
  Cmax_ADC    = NA_real_,
  Cmax_DXd_ng = NA_real_,
  Damage_max  = NA_real_,
  Neut_nadir  = NA_real_,
  Ret_nadir   = NA_real_,
  RBC_nadir   = NA_real_,
  Plt_nadir   = NA_real_,
  Grade_Neut  = NA_character_,
  Grade_Anemia= NA_character_
)

pb_step <- floor(N_patients / 10)

for (i in 1:N_patients) {

  # Tirage IIV
  eta_CL   <- rnorm(1, 0, omega_CL)
  eta_V1   <- rnorm(1, 0, omega_V1)
  eta_SCMP <- rnorm(1, 0, omega_Slope_CMP)
  eta_SMEP <- rnorm(1, 0, omega_Slope_MEP)

  pars_i <- pars_typ
  pars_i$CL_ADC    <- pars_typ$CL_ADC  * exp(eta_CL)
  pars_i$V1_ADC    <- pars_typ$V1_ADC  * exp(eta_V1)
  pars_i$Slope_CMP <- pars_typ$Slope_CMP * exp(eta_SCMP)
  pars_i$Slope_MEP <- pars_typ$Slope_MEP * exp(eta_SMEP)

  pars_i$rate_fun  <- make_tdxd_infusion(
    dose_mgkg  = 5.4, BW_kg = 70,
    Tinfu_h    = 1.5,
    interval_h = 21 * 24,
    n_cycles   = 6
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
    results$CL_ADC[i]      <- pars_i$CL_ADC
    results$V1_ADC[i]      <- pars_i$V1_ADC
    results$Slope_CMP[i]   <- pars_i$Slope_CMP
    results$Slope_MEP[i]   <- pars_i$Slope_MEP
    results$Cmax_ADC[i]    <- max(out$C_ADC1,  na.rm=TRUE)
    results$Cmax_DXd_ng[i] <- max(out$C_DXd,   na.rm=TRUE) * 1000
    results$Damage_max[i]  <- max(out$Damage,  na.rm=TRUE)
    results$Neut_nadir[i]  <- min(out$Neut,    na.rm=TRUE)
    results$Ret_nadir[i]   <- min(out$Ret,     na.rm=TRUE)
    results$RBC_nadir[i]   <- min(out$RBC,     na.rm=TRUE)
    results$Plt_nadir[i]   <- min(out$Plt,     na.rm=TRUE)
    results$Grade_Neut[i]  <- ctcae_neut(results$Neut_nadir[i])
    results$Grade_Anemia[i]<- ctcae_anemia(results$RBC_nadir[i], pars_typ$RBC0)
  }

  if (i %% pb_step == 0)
    cat(sprintf("  %3d/%d patients simulés...\n", i, N_patients))
}

results <- results[!is.na(results$Neut_nadir), ]
n_ok <- nrow(results)

# ══════════════════════════════════════════════════════════
# Résultats
# ══════════════════════════════════════════════════════════
grade_order <- c("G0", "G1", "G2", "G3", "G4")

tab_neut  <- table(factor(results$Grade_Neut,   levels = grade_order))
tab_anemia<- table(factor(results$Grade_Anemia, levels = grade_order))

pct_n <- round(100 * tab_neut  / n_ok, 1)
pct_a <- round(100 * tab_anemia/ n_ok, 1)

cat("\n═══════════════════════════════════════════════════════════\n")
cat(sprintf("  SIMULATION POPULATION — T-DXd 5.4 mg/kg Q3W × 6 (N=%d)\n", n_ok))
cat("─────────────────────────────────────────────────────────\n")
cat("  NEUTROPÉNIE :\n")
cat(sprintf("    %-5s  %6s  %7s   %s\n", "Grade", "n", "%", "FDA obs."))
for (g in grade_order) {
  fda_ref <- switch(g, G0="~71%", G1="~7%", G2="~7%", G3="~13%", G4="~3%")
  cat(sprintf("    %-5s  %6d  %6.1f%%   %s\n", g, tab_neut[g], pct_n[g], fda_ref))
}
cat(sprintf("    TOTAL G3-4 : %.1f%%  (FDA : ~16-20%%)\n",
            pct_n["G3"] + pct_n["G4"]))
cat(sprintf("    TOUT GRADE : %.1f%%  (FDA : ~29%%)\n",
            pct_n["G1"] + pct_n["G2"] + pct_n["G3"] + pct_n["G4"]))

cat("\n  ANÉMIE (proxy RBC) :\n")
cat(sprintf("    %-5s  %6s  %7s   %s\n", "Grade", "n", "%", "FDA obs."))
for (g in grade_order) {
  fda_ref <- switch(g, G0="~30%", G1="~37%", G2="~24%", G3="~8%", G4="~1%")
  cat(sprintf("    %-5s  %6d  %6.1f%%   %s\n", g, tab_anemia[g], pct_a[g], fda_ref))
}
cat(sprintf("    TOTAL G3-4 : %.1f%%  (FDA : ~9%%)\n",
            pct_a["G3"] + pct_a["G4"]))

cat("\n─────────────────────────────────────────────────────────\n")
cat("  Statistiques exposition (médiane [P10-P90]) :\n")
quants <- function(x) quantile(x, c(0.1, 0.5, 0.9), na.rm=TRUE)
q_adc  <- quants(results$Cmax_ADC)
q_dxd  <- quants(results$Cmax_DXd_ng)
q_neut <- quants(results$Neut_nadir)
cat(sprintf("    ADC Cmax  : %.0f [%.0f–%.0f] µg/mL  (FDA=122)\n",
            q_adc[2], q_adc[1], q_adc[3]))
cat(sprintf("    DXd Cmax  : %.1f [%.1f–%.1f] ng/mL  (FDA=4.4)\n",
            q_dxd[2], q_dxd[1], q_dxd[3]))
cat(sprintf("    Neut nadir: %.2f [%.2f–%.2f] × 10⁹/L\n",
            q_neut[2], q_neut[1], q_neut[3]))
cat("═══════════════════════════════════════════════════════════\n")

# ══════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════
pdf("results_PKPD_human/PKPD_human_population_VPC.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1))

grade_cols <- c(G0="#2166ac", G1="#74add1", G2="#f4a582", G3="#d6604d", G4="#b2182b")

# ── 1. Distribution Neut nadir ──
hist(results$Neut_nadir, breaks = 30,
     col = "#74add1", border = "white",
     main = "Distribution nadir Neutrophiles",
     xlab = "Neut nadir [× 10⁹/L]",
     ylab = "Fréquence")
abline(v = c(2.0, 1.5, 1.0, 0.5), lty = 2,
       col = c("grey60", "orange", "red", "darkred"), lwd = 1.5)
text(c(2.0, 1.5, 1.0, 0.5) + 0.05, par("usr")[4] * 0.95,
     c("G1", "G2", "G3", "G4"),
     col = c("grey60", "orange", "red", "darkred"), cex = 0.85, adj = 0)
abline(v = median(results$Neut_nadir, na.rm=TRUE), col = "#2166ac", lwd = 2)

# ── 2. Distribution % chute Neut ──
pct_drop_neut <- 100 * (results$Neut_nadir - pars_typ$Neut0) / pars_typ$Neut0
hist(pct_drop_neut, breaks = 30,
     col = "#74add1", border = "white",
     main = "Chute Neutrophiles (%Δ baseline)",
     xlab = "%Δ Neut nadir", ylab = "Fréquence")
abline(v = median(pct_drop_neut, na.rm=TRUE), col = "#2166ac", lwd = 2)
text(median(pct_drop_neut)-1, par("usr")[4]*0.9,
     sprintf("Médiane\n%.1f%%", median(pct_drop_neut)), col="#2166ac", adj=1, cex=0.85)

# ── 3. Barplot grades neutropénie : modèle vs FDA ──
grades_mod <- pct_n[grade_order]
grades_fda <- c(71, 7, 7, 13, 3)   # approximatifs FDA

bar_data <- rbind(grades_mod, grades_fda)
bp <- barplot(bar_data,
              beside = TRUE,
              col    = c("#2166ac", "#d6604d"),
              names.arg = grade_order,
              ylim   = c(0, 100),
              main   = "Grades Neutropénie : modèle vs FDA",
              ylab   = "Patients (%)",
              legend.text = c("Modèle", "FDA (DESTINY-B01)"),
              args.legend = list(bty = "n", cex = 0.85))
text(bp[1,], grades_mod + 2, sprintf("%.0f%%", grades_mod), cex = 0.8, col = "#2166ac")

# ── 4. Distribution % chute RBC (anémie) ──
pct_drop_rbc <- 100 * (results$RBC_nadir - pars_typ$RBC0) / pars_typ$RBC0
hist(pct_drop_rbc, breaks = 30,
     col = "#f4a582", border = "white",
     main = "Chute RBC (%Δ baseline, proxy anémie)",
     xlab = "%Δ RBC nadir", ylab = "Fréquence")
abline(v = median(pct_drop_rbc, na.rm=TRUE), col = "#d6604d", lwd = 2)

# ── 5. Neut nadir vs CL_ADC ──
plot(results$CL_ADC * 24, results$Neut_nadir,
     pch = 16, cex = 0.6,
     col = grade_cols[results$Grade_Neut],
     xlab = "CL_ADC [L/jour]", ylab = "Neut nadir [× 10⁹/L]",
     main = "Neut nadir vs CL_ADC")
abline(h = c(2.0, 1.5, 1.0, 0.5), lty = 2, col = c("grey60","orange","red","darkred"))
legend("topright", names(grade_cols), pch = 16,
       col = grade_cols, bty = "n", cex = 0.8)

# ── 6. Neut nadir vs Slope_CMP ──
plot(results$Slope_CMP, results$Neut_nadir,
     pch = 16, cex = 0.6,
     col = grade_cols[results$Grade_Neut],
     xlab = "Slope_CMP (sensibilité PD)", ylab = "Neut nadir [× 10⁹/L]",
     main = "Neut nadir vs Slope_CMP")
abline(h = c(2.0, 1.5, 1.0, 0.5), lty = 2, col = c("grey60","orange","red","darkred"))
legend("topright", names(grade_cols), pch = 16,
       col = grade_cols, bty = "n", cex = 0.8)

dev.off()
cat("  -> results_PKPD_human/PKPD_human_population_VPC.pdf\n")
