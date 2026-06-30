############################################################
# run_bpa_simu2.R
# Simulation de population -- anti-FLT3_chBPA-STINGa20 ADC
# Scenario 2 (interne) : IC50_MEP=155 nM, IC50_CMP=0.015 nM
#
# ATTENTION : Slope_CMP ~ 747 000 000 (EXTREME -- Slope_tdxd_hu x 190/0.015)
#             kill_CMP=1 garanti des la 1ere dose pour 100% des patients
# PK : placeholder T-DXd (Yin 2020)
# IIV log-normal (sans mixture) : ω_CL=0.35 ω_V1=0.20 ω_Slope=0.33
############################################################
library(deSolve)

source("../etape2_carboplatin_humain/parameters_human.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("../etape5_tdxd_humain/parameters_tdxd_human.R")
source("../etape3_tdxd_rat/parameters_tdxd_rat.R")
source("../etape3_tdxd_rat/pkpd_tdxd_rat.R")
source("parameters_bpa_human.R")

dir.create("results", showWarnings = FALSE)

cat(sprintf("\n!!! Slope_CMP=%.0f (IC50=0.015 nM) -- borne superieure theorique\n\n",
            Slope_CMP_BPA_s2))

set.seed(42)
N_patients <- 300

pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL
pars_typ <- c(pars_pd_hu, bpa_pars)

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)

times_hu <- seq(0, 126 * 24, by = 6)
n_times  <- length(times_hu)
times_j  <- times_hu / 24

ctcae_neut <- function(x) {
  if (x < 0.5) "G4" else if (x < 1.0) "G3" else if (x < 1.5) "G2" else if (x < 2.0) "G1" else "G0"
}
ctcae_anemia <- function(x, x0) {
  f <- x / x0
  if (f < 0.54) "G4" else if (f < 0.67) "G3" else if (f < 0.80) "G2" else if (f < 0.90) "G1" else "G0"
}
ctcae_plt <- function(x) {
  if (x < 25) "G4" else if (x < 50) "G3" else if (x < 75) "G2" else if (x < 150) "G1" else "G0"
}

cat(sprintf("[Simu 2 -- interne] IC50_MEP=155 nM  IC50_CMP=0.015 nM\n"))
cat(sprintf("Slope_MEP=%.3f  Slope_CMP=%.0f\n", Slope_MEP_BPA_s2, Slope_CMP_BPA_s2))
cat(sprintf("N=%d patients, %.1f mg/kg Q3W x %d cycles  (~3-5 min)\n\n",
            N_patients, DOSE_MGKG, N_CYCLES))

results <- data.frame(
  id = 1:N_patients,
  CL_ADC = NA_real_, V1_ADC = NA_real_,
  Slope_CMP = NA_real_, Slope_MEP = NA_real_,
  Cmax_ADC = NA_real_, Damage_max = NA_real_,
  Neut_nadir = NA_real_, Ret_nadir = NA_real_,
  RBC_nadir  = NA_real_, Plt_nadir  = NA_real_,
  Grade_Neut = NA_character_, Grade_Anemia = NA_character_, Grade_Plt = NA_character_
)

mk <- function() matrix(NA_real_, nrow = n_times, ncol = N_patients)
mat_MPP  <- mk(); mat_CMP  <- mk(); mat_MEP  <- mk()
mat_Neut <- mk(); mat_Mono <- mk()
mat_Ret  <- mk(); mat_RBC  <- mk(); mat_Plt  <- mk()

pb_step <- floor(N_patients / 10)

for (i in 1:N_patients) {
  pars_i <- pars_typ
  pars_i$CL_ADC    <- pars_typ$CL_ADC   * exp(rnorm(1, 0, omega_CL))
  pars_i$V1_ADC    <- pars_typ$V1_ADC   * exp(rnorm(1, 0, omega_V1))
  pars_i$Slope_CMP <- max(0, Slope_CMP_BPA_s2 * exp(rnorm(1, 0, omega_Slope_CMP)))
  pars_i$Slope_MEP <- max(0, Slope_MEP_BPA_s2 * exp(rnorm(1, 0, omega_Slope_MEP)))
  pars_i$rate_fun  <- make_tdxd_infusion(
    dose_mgkg = DOSE_MGKG, BW_kg = BW_KG,
    Tinfu_h = TINFU_H, interval_h = INTERVAL_H, n_cycles = N_CYCLES)

  out <- tryCatch(suppressMessages(suppressWarnings(as.data.frame(lsoda(
    y = state0_hu, times = times_hu, func = pkpd_tdxd_fornari, parms = pars_i,
    rtol = 1e-5, atol = 1e-7, maxsteps = 500000)))), error = function(e) NULL)

  if (!is.null(out) && nrow(out) > 10) {
    results$CL_ADC[i]      <- pars_i$CL_ADC
    results$V1_ADC[i]      <- pars_i$V1_ADC
    results$Slope_CMP[i]   <- pars_i$Slope_CMP
    results$Slope_MEP[i]   <- pars_i$Slope_MEP
    results$Cmax_ADC[i]    <- max(out$C_ADC1, na.rm = TRUE)
    results$Damage_max[i]  <- max(out$Damage,  na.rm = TRUE)
    results$Neut_nadir[i]  <- min(out$Neut,    na.rm = TRUE)
    results$Ret_nadir[i]   <- min(out$Ret,     na.rm = TRUE)
    results$RBC_nadir[i]   <- min(out$RBC,     na.rm = TRUE)
    results$Plt_nadir[i]   <- min(out$Plt,     na.rm = TRUE)
    results$Grade_Neut[i]  <- ctcae_neut(results$Neut_nadir[i])
    results$Grade_Anemia[i]<- ctcae_anemia(results$RBC_nadir[i], pars_typ$RBC0)
    results$Grade_Plt[i]   <- ctcae_plt(results$Plt_nadir[i])

    if (nrow(out) >= n_times) {
      n <- n_times
      mat_MPP[, i]  <- out$MPP[1:n];   mat_CMP[, i]  <- out$CMP[1:n]
      mat_MEP[, i]  <- out$MEP[1:n];   mat_Neut[, i] <- out$Neut[1:n]
      mat_Mono[, i] <- out$Mono[1:n];  mat_Ret[, i]  <- out$Ret[1:n]
      mat_RBC[, i]  <- out$RBC[1:n];   mat_Plt[, i]  <- out$Plt[1:n]
    }
  }
  if (i %% pb_step == 0) cat(sprintf("  %d/%d...\n", i, N_patients))
}

results <- results[!is.na(results$Neut_nadir), ]
n_ok    <- nrow(results)
ok_cols <- colSums(!is.na(mat_Neut)) > (n_times / 2)
for (nm in c("mat_MPP","mat_CMP","mat_MEP","mat_Neut",
             "mat_Mono","mat_Ret","mat_RBC","mat_Plt"))
  assign(nm, get(nm)[, ok_cols, drop = FALSE])
has_profiles <- sum(ok_cols) > 0
cat(sprintf("  Profils complets : %d/%d patients\n", sum(ok_cols), n_ok))

saveRDS(results, "results/bpa_simu2_results.rds")

# ══════════════════════════════════════════════════════════
# Console
# ══════════════════════════════════════════════════════════
grade_order <- c("G0","G1","G2","G3","G4")
tab_n <- table(factor(results$Grade_Neut,   levels = grade_order))
tab_a <- table(factor(results$Grade_Anemia, levels = grade_order))
tab_p <- table(factor(results$Grade_Plt,    levels = grade_order))
pct_n <- round(100 * tab_n / n_ok, 1)
pct_a <- round(100 * tab_a / n_ok, 1)
pct_p <- round(100 * tab_p / n_ok, 1)

cat("\n=============================================================\n")
cat(sprintf("  BPA SIMU 2 (interne) -- %.1f mg/kg Q3W x %d (N=%d)  [EXTREME]\n", DOSE_MGKG, N_CYCLES, n_ok))
cat(sprintf("  Slope_MEP=%.3f | Slope_CMP=%.0f\n", Slope_MEP_BPA_s2, Slope_CMP_BPA_s2))
cat("  NEUTROPENIE :")
for (g in grade_order) cat(sprintf("  %s=%.0f%%", g, pct_n[g]))
cat(sprintf("\n  -> G3-4=%.1f%%\n", pct_n["G3"]+pct_n["G4"]))
cat("  ANEMIE      :")
for (g in grade_order) cat(sprintf("  %s=%.0f%%", g, pct_a[g]))
cat(sprintf("\n  -> G3-4=%.1f%%\n", pct_a["G3"]+pct_a["G4"]))
cat("  THROMBO     :")
for (g in grade_order) cat(sprintf("  %s=%.0f%%", g, pct_p[g]))
cat(sprintf("\n  -> G3-4=%.1f%%\n", pct_p["G3"]+pct_p["G4"]))
cat("=============================================================\n")

# ══════════════════════════════════════════════════════════
# Helper profil log-scale
# ══════════════════════════════════════════════════════════
.plot_profile <- function(mat, times_j, baseline, col, main, ylab,
                          hlines = NULL, hcols = NULL) {
  p10 <- apply(mat, 1, quantile, 0.10, na.rm = TRUE)
  p50 <- apply(mat, 1, quantile, 0.50, na.rm = TRUE)
  p90 <- apply(mat, 1, quantile, 0.90, na.rm = TRUE)

  pos <- p50[p50 > 0 & is.finite(p50)]
  ymin <- if (length(pos) > 0) max(1e-6, min(pos)*0.3) else 1e-3
  ymax <- max(p90[is.finite(p90)], na.rm = TRUE) * 3
  if (!is.finite(ymax) || ymax <= ymin) ymax <- ymin * 1000

  plot(times_j, pmax(p50, ymin), type = "n", log = "y",
       ylim = c(ymin, ymax), xlim = range(times_j),
       xlab = "Time (d)", ylab = ylab, main = main,
       panel.first = {
         abline(v = seq(0, max(times_j), by = 25), col = "grey90")
         abline(h = 10^(-6:6), col = "grey90")
       })

  px <- c(times_j, rev(times_j))
  py <- c(pmax(p10, ymin), rev(pmax(p90, ymin)))
  valid <- is.finite(py) & py > 0
  if (sum(valid) > 4)
    polygon(px[valid | c(tail(valid, length(times_j)), head(valid, length(times_j)))],
            py, col = adjustcolor(col, alpha.f = 0.25), border = NA)

  lines(times_j, pmax(p50, ymin), col = col, lwd = 2)

  if (!is.null(baseline) && is.finite(baseline) && baseline > 0)
    abline(h = baseline, lty = 2, col = "grey40", lwd = 1.2)
  if (!is.null(hlines))
    abline(h = hlines[hlines > ymin], lty = 2, col = hcols[hlines > ymin], lwd = 1.2)

  cycle_days <- (0:(N_CYCLES - 1)) * (INTERVAL_H / 24)
  abline(v = cycle_days, lty = 3, col = "grey60", lwd = 0.8)
}

# ══════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════
grade_cols <- c(G0="#2166ac",G1="#74add1",G2="#f4a582",G3="#d6604d",G4="#b2182b")

pdf("results/BPA_simu2_VPC.pdf", width = 16, height = 9)

# ── Page 1 : profils 2×4 style Fornari ───────────────────
if (has_profiles) {
  par(mfrow = c(2, 4), mar = c(3.5, 3.8, 2.8, 0.5), oma = c(0, 0, 2.5, 0))

  .plot_profile(mat_MPP,  times_j, pars_typ$MPP0,  "#4575b4", "Multi-potent progenitors", "10^9 cells L-1")
  .plot_profile(mat_Neut, times_j, pars_typ$Neut0, "#d73027", "Neutrophils",              "10^9 cells L-1",
                hlines=c(2,1,0.5), hcols=c("grey50","red","darkred"))
  .plot_profile(mat_CMP,  times_j, pars_typ$CMP0,  "#1a9641", "Common myeloid progenitors","10^9 cells L-1")
  .plot_profile(mat_Mono, times_j, pars_typ$Mono0, "#8c510a", "Monocytes",               "10^9 cells L-1")
  .plot_profile(mat_MEP,  times_j, pars_typ$MEP0,  "#762a83", "MEP",                     "10^9 cells L-1")
  .plot_profile(mat_Plt,  times_j, pars_typ$Plt0,  "#35978f", "Platelets",               "10^9 cells L-1",
                hlines=c(150,75,50,25), hcols=c("grey50","orange","red","darkred"))
  .plot_profile(mat_Ret,  times_j, pars_typ$Ret0,  "#f1b6da", "Reticulocytes",           "10^9 cells L-1")
  .plot_profile(mat_RBC,  times_j, pars_typ$RBC0,  "#c51b7d", "Red blood cells",         "10^9 cells L-1",
                hlines=c(pars_typ$RBC0*0.67, pars_typ$RBC0*0.80), hcols=c("darkred","orange"))

  mtext(sprintf("BPA Simu2 [EXTREME] -- %.1f mg/kg Q3W x%d  |  Slope_CMP=%.0f  Slope_MEP=%.3f  (N=%d)",
                DOSE_MGKG, N_CYCLES, Slope_CMP_BPA_s2, Slope_MEP_BPA_s2, n_ok),
        outer = TRUE, cex = 1.1, font = 2)
}

# ── Page 2 : grades + comparaison S1 vs S2 ───────────────
par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1), oma = c(0,0,0,0))

bp <- barplot(pct_n[grade_order], col=grade_cols[grade_order], names.arg=grade_order,
              ylim=c(0,100), main="Grades Neutropenie (Simu2)", ylab="Patients (%)")
text(bp, pct_n[grade_order]+2, sprintf("%.0f%%",pct_n[grade_order]), cex=0.85)

bp <- barplot(pct_a[grade_order], col=grade_cols[grade_order], names.arg=grade_order,
              ylim=c(0,100), main="Grades Anemie (Simu2)", ylab="Patients (%)")
text(bp, pct_a[grade_order]+2, sprintf("%.0f%%",pct_a[grade_order]), cex=0.85)

bp <- barplot(pct_p[grade_order], col=grade_cols[grade_order], names.arg=grade_order,
              ylim=c(0,100), main="Grades Thrombocytopenie (Simu2)", ylab="Patients (%)")
text(bp, pct_p[grade_order]+2, sprintf("%.0f%%",pct_p[grade_order]), cex=0.85)

# Comparaison Simu1 vs Simu2
simu1_path <- "results/bpa_simu1_results.rds"
if (file.exists(simu1_path)) {
  res1 <- readRDS(simu1_path)
  pct1 <- round(100*table(factor(res1$Grade_Neut,levels=grade_order))/nrow(res1),1)
  cmp  <- rbind(pct1[grade_order], pct_n[grade_order])
  bp   <- barplot(cmp, beside=TRUE, col=c("#2166ac","#b2182b"),
                  names.arg=grade_order, ylim=c(0,105),
                  main="Neutropenie : Simu1 vs Simu2", ylab="Patients (%)")
  legend("topright", c("Simu1 (sous-traitant)","Simu2 (interne)"),
         fill=c("#2166ac","#b2182b"), bty="n", cex=0.85)
} else {
  plot.new(); text(0.5,0.5,"Lancer run_bpa_simu1.R\npour la comparaison", cex=1.2)
}

hist(results$Neut_nadir, breaks=30, col="#d6604d", border="white",
     main="Nadir Neutrophiles (Simu2)", xlab="Neut nadir [x 10^9/L]", ylab="Frequence")
abline(v=c(2,1,0.5), lty=2, col=c("grey60","red","darkred"))

pct_drop <- 100*(results$RBC_nadir - pars_typ$RBC0)/pars_typ$RBC0
hist(pct_drop, breaks=30, col="#f4a582", border="white",
     main="Chute RBC (Simu2)", xlab="%Delta RBC", ylab="Frequence")
abline(v=median(pct_drop,na.rm=TRUE), col="#d6604d", lwd=2)

dev.off()
cat("  -> results/BPA_simu2_VPC.pdf\n")
