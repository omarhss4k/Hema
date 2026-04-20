############################################################
# run_pkpd_tdxd_human_population.R
# Simulation de population — T-DXd 5.4 mg/kg Q3W × 6 cycles
#
# Objectif : reproduire la distribution de toxicité CTCAE
#   FDA BLA 761139 (cibles primaires DESTINY-Breast01, n=184) :
#   Neutropénie G3-4 : ~20%    Anémie G3-4 : ~9%
#   Neutropénie tout grade : ~29%  Anémie tout grade : ~70%
#   FDA pooled 5.4 mg/kg (N=234, Table 44) :
#   Neutropénie G3-4 : 16.2%  Anémie G3-4 : 7.3%  Thrombocy. G3-4 : 3.4%
#
# MODÈLE DE MÉLANGE BIMODAL (mixture) :
#   71% résistants  : Slope_CMP = Slope_resist_tdxd   × exp(η)  → G0 garanti
#   29% sensibles   : Slope_CMP = Slope_sensitive_tdxd × exp(η)  → G3-4 possible
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

# Note : Slope_CMP est géré par la mixture (Slope_resist / Slope_sensitive)
# Slope_CMP_tdxd_human = Slope_sensitive_tdxd (compat. backward, non utilisé ici)

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

# ── CTCAE Grade Thrombocytopénie (Plt en ×10⁹/L) ─────────
# Plt0 = 345 ×10⁹/L (Fornari 2019 humain)
# CTCAE v5 : G1 : <LLN-75  G2 : 50-<75  G3 : 25-<50  G4 : <25
ctcae_plt <- function(plt) {
  if      (plt < 25)  "G4"
  else if (plt < 50)  "G3"
  else if (plt < 75)  "G2"
  else if (plt < 150) "G1"   # LLN ≈ 150 ×10⁹/L
  else                "G0"
}

# ══════════════════════════════════════════════════════════
# Simulation population
# ══════════════════════════════════════════════════════════
cat(sprintf("Simulation population : N=%d patients, 5.4 mg/kg Q3W × 6 cycles\n", N_patients))
cat("(patience ~3-5 min)\n\n")

results <- data.frame(
  id               = 1:N_patients,
  group_cmp        = NA_character_,   # "resistant" / "moderate" / "sensitive"
  group_mep        = NA_character_,   # "resistant" / "light" / "sensitive"
  CL_ADC       = NA_real_,
  V1_ADC       = NA_real_,
  Slope_CMP    = NA_real_,
  Slope_MEP    = NA_real_,
  Cmax_ADC     = NA_real_,
  Cmax_DXd_ng  = NA_real_,
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

  # Tirage IIV
  eta_CL   <- rnorm(1, 0, omega_CL)
  eta_V1   <- rnorm(1, 0, omega_V1)
  eta_SCMP <- rnorm(1, 0, omega_Slope_CMP)
  eta_SMEP <- rnorm(1, 0, omega_Slope_MEP)

  # Mixture TRIMODALE CMP (neutropénie) : tirage indépendant
  p_resist <- 1 - p_sensitive_tdxd - p_moderate_tdxd
  rand_cmp <- runif(1)
  group_cmp_i <- if (rand_cmp < p_resist) {
    "resistant"
  } else if (rand_cmp < p_resist + p_moderate_tdxd) {
    "moderate"
  } else {
    "sensitive"
  }

  # Mixture TRIMODALE MEP (anémie) : tirage indépendant
  p_mep_resist <- 1 - p_MEP_light_tdxd - p_sensitive_mep_tdxd
  rand_mep <- runif(1)
  group_mep_i <- if (rand_mep < p_mep_resist) {
    "resistant"
  } else if (rand_mep < p_mep_resist + p_MEP_light_tdxd) {
    "light"
  } else {
    "sensitive"
  }

  eta_SMEP_mep <- rnorm(1, 0, omega_Slope_MEP_sensitive)

  pars_i <- pars_typ
  pars_i$CL_ADC  <- pars_typ$CL_ADC * exp(eta_CL)
  pars_i$V1_ADC  <- pars_typ$V1_ADC * exp(eta_V1)
  pars_i$Slope_CMP <- switch(group_cmp_i,
    resistant = Slope_resist_tdxd    * exp(eta_SCMP),
    moderate  = Slope_moderate_tdxd  * exp(eta_SCMP),
    sensitive = Slope_sensitive_tdxd * exp(eta_SCMP)
  )
  pars_i$Slope_MEP <- switch(group_mep_i,
    resistant = Slope_MEP_resist_tdxd    * exp(eta_SMEP_mep),
    light     = Slope_MEP_light_tdxd     * exp(eta_SMEP_mep),
    sensitive = Slope_MEP_sensitive_tdxd * exp(eta_SMEP_mep)
  )

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
    results$group_cmp[i] <- group_cmp_i
    results$group_mep[i] <- group_mep_i
    results$CL_ADC[i]           <- pars_i$CL_ADC
    results$V1_ADC[i]       <- pars_i$V1_ADC
    results$Slope_CMP[i]    <- pars_i$Slope_CMP
    results$Slope_MEP[i]    <- pars_i$Slope_MEP
    results$Cmax_ADC[i]    <- max(out$C_ADC1,  na.rm=TRUE)
    results$Cmax_DXd_ng[i] <- max(out$C_DXd,   na.rm=TRUE) * 1000
    results$Damage_max[i]  <- max(out$Damage,  na.rm=TRUE)
    results$Neut_nadir[i]  <- min(out$Neut,    na.rm=TRUE)
    results$Ret_nadir[i]   <- min(out$Ret,     na.rm=TRUE)
    results$RBC_nadir[i]   <- min(out$RBC,     na.rm=TRUE)
    results$Plt_nadir[i]   <- min(out$Plt,     na.rm=TRUE)
    results$Grade_Neut[i]  <- ctcae_neut(results$Neut_nadir[i])
    results$Grade_Anemia[i]<- ctcae_anemia(results$RBC_nadir[i], pars_typ$RBC0)
    results$Grade_Plt[i]   <- ctcae_plt(results$Plt_nadir[i])
  }

  if (i %% pb_step == 0)
    cat(sprintf("  %3d/%d patients simulés...\n", i, N_patients))
}

results <- results[!is.na(results$Neut_nadir), ]
n_ok <- nrow(results)
saveRDS(results, "results_PKPD_human/population_results.rds")

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
cat(sprintf("  SIMULATION POPULATION — T-DXd 5.4 mg/kg Q3W × 6 (N=%d)\n", n_ok))
cat("─────────────────────────────────────────────────────────\n")
cat("  NEUTROPÉNIE :\n")
cat(sprintf("    %-5s  %6s  %7s   %s\n", "Grade", "n", "%", "FDA obs."))
for (g in grade_order) {
  fda_ref <- switch(g, G0="~71%", G1="~7%", G2="~7%", G3="~13%", G4="~3%")
  cat(sprintf("    %-5s  %6d  %6.1f%%   %s\n", g, tab_neut[g], pct_n[g], fda_ref))
}
cat(sprintf("    TOTAL G3-4 : %.1f%%  (FDA : ~20%% U201 / 16.2%% pooled N=234)\n",
            pct_n["G3"] + pct_n["G4"]))
cat(sprintf("    TOUT GRADE : %.1f%%  (FDA : ~29%% U201 / 29.5%% pooled)\n",
            pct_n["G1"] + pct_n["G2"] + pct_n["G3"] + pct_n["G4"]))

cat("\n  ANÉMIE (proxy RBC) :\n")
cat(sprintf("    %-5s  %6s  %7s   %s\n", "Grade", "n", "%", "FDA obs."))
for (g in grade_order) {
  fda_ref <- switch(g, G0="~30%", G1="~37%", G2="~24%", G3="~8%", G4="~1%")
  cat(sprintf("    %-5s  %6d  %6.1f%%   %s\n", g, tab_anemia[g], pct_a[g], fda_ref))
}
cat(sprintf("    TOTAL G3-4 : %.1f%%  (FDA : ~9%%)\n",
            pct_a["G3"] + pct_a["G4"]))

cat("\n  THROMBOCYTOPÉNIE (Plt) :\n")
cat(sprintf("    %-5s  %6s  %7s   %s\n", "Grade", "n", "%", "FDA obs. (N=234)"))
for (g in grade_order) {
  fda_ref <- switch(g, G0="~63%", G1="~30%", G2="~4%", G3="~2%", G4="~1%")
  cat(sprintf("    %-5s  %6d  %6.1f%%   %s\n", g, tab_plt[g], pct_p[g], fda_ref))
}
cat(sprintf("    TOTAL G3-4 : %.1f%%  (FDA : ~3.4%%)\n",
            pct_p["G3"] + pct_p["G4"]))

cat("\n─────────────────────────────────────────────────────────\n")
# ── Breakdown par sous-groupe trimodal ──
res_r <- results[results$group_cmp == "resistant", ]
res_m <- results[results$group_cmp == "moderate",  ]
res_s <- results[results$group_cmp == "sensitive",  ]
n_r <- nrow(res_r); n_m <- nrow(res_m); n_s <- nrow(res_s)
cat(sprintf("  Sous-groupes (trimodal) : %d résistants (%.0f%%) | %d modérés (%.0f%%) | %d sensibles (%.0f%%)\n",
            n_r, 100*n_r/n_ok, n_m, 100*n_m/n_ok, n_s, 100*n_s/n_ok))
pct_g34 <- function(x) if (nrow(x)>0) 100*sum(x$Neut_nadir<1.0)/nrow(x) else 0
pct_g12 <- function(x) if (nrow(x)>0) 100*sum(x$Neut_nadir>=1.0 & x$Neut_nadir<2.0)/nrow(x) else 0
nadir_q  <- function(x) {
  sprintf("%.2f [%.2f-%.2f]",
          median(x$Neut_nadir,na.rm=T),
          quantile(x$Neut_nadir,0.1,na.rm=T),
          quantile(x$Neut_nadir,0.9,na.rm=T))
}
cat(sprintf("  Résistants  : G3-4=%.1f%%  G1-2=%.1f%%  nadir=%s\n",
            pct_g34(res_r), pct_g12(res_r), nadir_q(res_r)))
cat(sprintf("  Modérés     : G3-4=%.1f%%  G1-2=%.1f%%  nadir=%s\n",
            pct_g34(res_m), pct_g12(res_m), nadir_q(res_m)))
cat(sprintf("  Sensibles   : G3-4=%.1f%%  G1-2=%.1f%%  nadir=%s\n",
            pct_g34(res_s), pct_g12(res_s), nadir_q(res_s)))

# ── Breakdown anémie par sous-groupe MEP trimodal ──
res_mep_r <- results[results$group_mep == "resistant", ]
res_mep_l <- results[results$group_mep == "light",     ]
res_mep_s <- results[results$group_mep == "sensitive",  ]
n_mr <- nrow(res_mep_r); n_ml <- nrow(res_mep_l); n_ms <- nrow(res_mep_s)
rbc0 <- pars_typ$RBC0
pct_anemia_grade <- function(x, lo, hi) {
  if (nrow(x)==0) return(0)
  100 * sum(x$RBC_nadir/rbc0 >= lo & x$RBC_nadir/rbc0 < hi) / nrow(x)
}
cat(sprintf("\n  Anémie — trimodal MEP : %d résistants (%.0f%%) | %d légers (%.0f%%) | %d sensibles (%.0f%%)\n",
            n_mr, 100*n_mr/n_ok, n_ml, 100*n_ml/n_ok, n_ms, 100*n_ms/n_ok))
for (grp in list(list(res_mep_r,"Résistants"), list(res_mep_l,"Légers"), list(res_mep_s,"Sensibles"))) {
  x <- grp[[1]]; nm <- grp[[2]]
  cat(sprintf("  %-11s: G0=%.0f%% G1=%.0f%% G2=%.0f%% G3-4=%.0f%%\n", nm,
    pct_anemia_grade(x,0.90,Inf), pct_anemia_grade(x,0.83,0.90),
    pct_anemia_grade(x,0.67,0.83), pct_anemia_grade(x,0,0.67)))
}

cat("\n─────────────────────────────────────────────────────────\n")
cat("  Statistiques exposition (médiane [P10-P90]) :\n")
quants <- function(x) quantile(x, c(0.1, 0.5, 0.9), na.rm=TRUE)
q_adc  <- quants(results$Cmax_ADC)
q_dxd  <- quants(results$Cmax_DXd_ng)
q_neut <- quants(results$Neut_nadir)
q_plt  <- quants(results$Plt_nadir)
cat(sprintf("    ADC Cmax  : %.0f [%.0f–%.0f] µg/mL  (FDA SS=122)\n",
            q_adc[2], q_adc[1], q_adc[3]))
cat(sprintf("    DXd Cmax  : %.1f [%.1f–%.1f] ng/mL  (FDA SS=4.4)\n",
            q_dxd[2], q_dxd[1], q_dxd[3]))
cat(sprintf("    Neut nadir: %.2f [%.2f–%.2f] × 10⁹/L\n",
            q_neut[2], q_neut[1], q_neut[3]))
cat(sprintf("    Plt nadir : %.0f [%.0f–%.0f] × 10⁹/L  (Plt0=%.0f)\n",
            q_plt[2], q_plt[1], q_plt[3], pars_typ$Plt0))
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

# ══════════════════════════════════════════════════════════
# Graphique focalisé : grades Neutropénie + Anémie
# ══════════════════════════════════════════════════════════
source("../scripts_tdxd/plot_grades_human.R")
