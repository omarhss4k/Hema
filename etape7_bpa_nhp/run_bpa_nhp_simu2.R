############################################################
# run_bpa_nhp_simu2.R
# Simulation -- anti-FLT3_chBPA-STINGa20 ADC -- SINGE CYNOMOLGUS
# Scenario 1 (interne) : IC50_MEP=155 nM, IC50_CMP=7.8 nM
#
# PK : allometrie souris->singe (BW=4 kg, 1-cmt)
#      CL_nhp ~ 0.000475 L/h, V1_nhp ~ 0.315 L, T½ ~ 19 jours
# PD : Fornari 2019 (parametres humain -- NHP ~ humain)
# IIV log-normal : omega_CL=0.35 omega_V1=0.20 omega_Slope=0.33
############################################################
library(deSolve)

source("../etape2_carboplatin_humain/parameters_human.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("../etape5_tdxd_humain/parameters_tdxd_human.R")
source("../etape3_tdxd_rat/parameters_tdxd_rat.R")
source("pkpd_bpa_1cmt.R")
source("parameters_bpa_nhp.R")

dir.create("results", showWarnings = FALSE)

set.seed(42)
N_animals <- 100   # animaux virtuels (VPC preclinique)

pars_pd <- init_pars
pars_pd[["k_dam"]] <- NULL
pars_pd[["k_rep"]] <- NULL
pars_typ <- c(pars_pd, bpa_pars_nhp)

bpa_state0_1cmt <- c(C_ADC1 = 0, C_DXd = 0, C_DXd_ic = 0, Damage = 0)
state_pd <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0   <- c(bpa_state0_1cmt, state_pd)

times <- seq(0, (N_CYCLES * INTERVAL_H + 7 * 24), by = 6)
n_times <- length(times)
times_j <- times / 24

ctcae_neut   <- function(x) {
  if (x < 0.5) "G4" else if (x < 1.0) "G3" else if (x < 1.5) "G2" else if (x < 2.0) "G1" else "G0"
}
ctcae_anemia <- function(x, x0) {
  f <- x / x0
  if (f < 0.54) "G4" else if (f < 0.67) "G3" else if (f < 0.80) "G2" else if (f < 0.90) "G1" else "G0"
}
ctcae_plt <- function(x) {
  if (x < 25) "G4" else if (x < 50) "G3" else if (x < 75) "G2" else if (x < 150) "G1" else "G0"
}

cat(sprintf("\n[NHP Simu 2] IC50_CMP=%.1f nM  Slope_CMP=%.0f\n",
            IC50_CMP_BPA_s2, Slope_CMP_BPA_s2))
cat(sprintf("BW=%.0f kg  %.1f mg/kg Q3W x%d  N=%d animaux\n\n",
            BW_KG, DOSE_MGKG, N_CYCLES, N_animals))

results <- data.frame(
  id          = 1:N_animals,
  CL_ADC      = NA_real_, V1_ADC    = NA_real_,
  Slope_CMP   = NA_real_, Slope_MEP = NA_real_,
  Cmax_ADC    = NA_real_,
  Neut_nadir  = NA_real_, Ret_nadir = NA_real_,
  RBC_nadir   = NA_real_, Plt_nadir  = NA_real_,
  Grade_Neut  = NA_character_, Grade_Anemia = NA_character_,
  Grade_Plt   = NA_character_
)

mk <- function() matrix(NA_real_, nrow = n_times, ncol = N_animals)
mat_MPP  <- mk(); mat_CMP  <- mk(); mat_MEP  <- mk()
mat_Neut <- mk(); mat_Mono <- mk()
mat_Ret  <- mk(); mat_RBC  <- mk(); mat_Plt  <- mk()

pb_step <- max(1, floor(N_animals / 10))

for (i in 1:N_animals) {
  pars_i <- pars_typ
  pars_i$CL_ADC    <- pars_typ$CL_ADC   * exp(rnorm(1, 0, omega_CL))
  pars_i$V1_ADC    <- pars_typ$V1_ADC   * exp(rnorm(1, 0, omega_V1))
  pars_i$Slope_CMP <- max(0, Slope_CMP_BPA_s2 * exp(rnorm(1, 0, omega_Slope_CMP)))
  pars_i$Slope_MEP <- max(0, Slope_MEP_BPA_s2 * exp(rnorm(1, 0, omega_Slope_MEP)))
  pars_i$rate_fun  <- make_tdxd_infusion(
    dose_mgkg = DOSE_MGKG, BW_kg = BW_KG,
    Tinfu_h = TINFU_H, interval_h = INTERVAL_H, n_cycles = N_CYCLES)

  out <- tryCatch(suppressMessages(suppressWarnings(as.data.frame(lsoda(
    y = state0, times = times, func = pkpd_bpa_fornari_1cmt, parms = pars_i,
    rtol = 1e-5, atol = 1e-7, maxsteps = 500000)))), error = function(e) NULL)

  if (!is.null(out) && nrow(out) > 10) {
    results$CL_ADC[i]      <- pars_i$CL_ADC
    results$V1_ADC[i]      <- pars_i$V1_ADC
    results$Slope_CMP[i]   <- pars_i$Slope_CMP
    results$Slope_MEP[i]   <- pars_i$Slope_MEP
    results$Cmax_ADC[i]    <- max(out$C_ADC1, na.rm = TRUE)
    results$Neut_nadir[i]  <- min(out$Neut,   na.rm = TRUE)
    results$Ret_nadir[i]   <- min(out$Ret,    na.rm = TRUE)
    results$RBC_nadir[i]   <- min(out$RBC,    na.rm = TRUE)
    results$Plt_nadir[i]   <- min(out$Plt,    na.rm = TRUE)
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
  if (i %% pb_step == 0) cat(sprintf("  %d/%d...\n", i, N_animals))
}

results  <- results[!is.na(results$Neut_nadir), ]
n_ok     <- nrow(results)
ok_cols  <- colSums(!is.na(mat_Neut)) > (n_times / 2)
for (nm in c("mat_MPP","mat_CMP","mat_MEP","mat_Neut",
             "mat_Mono","mat_Ret","mat_RBC","mat_Plt"))
  assign(nm, get(nm)[, ok_cols, drop = FALSE])
has_profiles <- sum(ok_cols) > 0
cat(sprintf("  Profils complets : %d/%d animaux\n", sum(ok_cols), n_ok))

saveRDS(results, "results/bpa_nhp_simu2_results.rds")

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
cat(sprintf("  BPA NHP SIMU 1 -- %.1f mg/kg Q3W x%d (N=%d animaux)\n",
            DOSE_MGKG, N_CYCLES, n_ok))
cat(sprintf("  Slope_MEP=%.3f | Slope_CMP=%.1f\n", Slope_MEP_BPA_s2, Slope_CMP_BPA_s2))
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
# Figures
# ══════════════════════════════════════════════════════════
.plot_profile <- function(mat, times_j, baseline, col, main, ylab,
                          hlines = NULL, hcols = NULL) {
  p10 <- apply(mat, 1, quantile, 0.10, na.rm = TRUE)
  p50 <- apply(mat, 1, quantile, 0.50, na.rm = TRUE)
  p90 <- apply(mat, 1, quantile, 0.90, na.rm = TRUE)

  ymin <- max(1e-3, min(p10[p10 > 0], na.rm = TRUE) * 0.5)
  ymax <- max(p90, na.rm = TRUE) * 2

  plot(times_j, p50, type = "n", log = "y",
       ylim = c(ymin, ymax), xlim = range(times_j),
       xlab = "Time (d)", ylab = ylab, main = main,
       panel.first = {
         abline(v = seq(0, max(times_j), by = 25), col = "grey90", lty = 1)
         abline(h = c(1e-3,1e-2,1e-1,1,10,1e2,1e3,1e4,1e5),
                col = "grey90", lty = 1)
       })

  px <- c(times_j, rev(times_j))
  py <- c(pmax(p10, ymin*0.9), rev(pmax(p90, ymin*0.9)))
  polygon(px, py, col = adjustcolor(col, alpha.f = 0.25), border = NA)
  lines(times_j, pmax(p50, ymin*0.9), col = col, lwd = 2)

  if (!is.null(baseline) && baseline > 0)
    abline(h = baseline, lty = 2, col = "grey40", lwd = 1.2)
  if (!is.null(hlines))
    abline(h = hlines, lty = 2, col = hcols, lwd = 1.2)

  cycle_days <- (0:(N_CYCLES - 1)) * (INTERVAL_H / 24)
  abline(v = cycle_days, lty = 3, col = "grey60", lwd = 0.8)
}

grade_cols <- c(G0="#2166ac",G1="#74add1",G2="#f4a582",G3="#d6604d",G4="#b2182b")

pdf("results/BPA_NHP_simu2_VPC.pdf", width = 16, height = 9)

if (has_profiles) {
  par(mfrow = c(2, 4), mar = c(3.5, 3.8, 2.8, 0.5), oma = c(0, 0, 2.5, 0))

  .plot_profile(mat_MPP,  times_j, pars_typ$MPP0,  "#4575b4",
                "Multi-potent progenitors", "10^9 cells L-1")
  .plot_profile(mat_Neut, times_j, pars_typ$Neut0, "#d73027",
                "Neutrophils", "10^9 cells L-1",
                hlines = c(2.0, 1.0, 0.5), hcols = c("grey50","red","darkred"))
  .plot_profile(mat_CMP,  times_j, pars_typ$CMP0,  "#1a9641",
                "Common myeloid progenitors", "10^9 cells L-1")
  .plot_profile(mat_Mono, times_j, pars_typ$Mono0, "#8c510a",
                "Monocytes", "10^9 cells L-1")

  .plot_profile(mat_MEP,  times_j, pars_typ$MEP0,  "#762a83",
                "MEP", "10^9 cells L-1")
  .plot_profile(mat_Plt,  times_j, pars_typ$Plt0,  "#35978f",
                "Platelets", "10^9 cells L-1",
                hlines = c(150, 75, 50, 25), hcols = c("grey50","orange","red","darkred"))
  .plot_profile(mat_Ret,  times_j, pars_typ$Ret0,  "#f1b6da",
                "Reticulocytes", "10^9 cells L-1")
  .plot_profile(mat_RBC,  times_j, pars_typ$RBC0,  "#c51b7d",
                "Red blood cells", "10^9 cells L-1",
                hlines = c(pars_typ$RBC0*0.67, pars_typ$RBC0*0.80),
                hcols  = c("darkred","orange"))

  mtext(sprintf("BPA NHP Simu2 -- %.1f mg/kg Q3W x%d  |  Slope_CMP=%.0f  Slope_MEP=%.3f  (N=%d)",
                DOSE_MGKG, N_CYCLES, Slope_CMP_BPA_s2, Slope_MEP_BPA_s2, n_ok),
        outer = TRUE, cex = 1.1, font = 2)
}

par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1), oma = c(0,0,0,0))

bp <- barplot(pct_n[grade_order], col = grade_cols[grade_order],
              names.arg = grade_order, ylim = c(0,100),
              main = "Grades Neutropenie (NHP)", ylab = "Animaux (%)")
text(bp, pct_n[grade_order]+2, sprintf("%.0f%%", pct_n[grade_order]), cex=0.85)

bp <- barplot(pct_a[grade_order], col = grade_cols[grade_order],
              names.arg = grade_order, ylim = c(0,100),
              main = "Grades Anemie (NHP)", ylab = "Animaux (%)")
text(bp, pct_a[grade_order]+2, sprintf("%.0f%%", pct_a[grade_order]), cex=0.85)

bp <- barplot(pct_p[grade_order], col = grade_cols[grade_order],
              names.arg = grade_order, ylim = c(0,100),
              main = "Grades Thrombocytopenie (NHP)", ylab = "Animaux (%)")
text(bp, pct_p[grade_order]+2, sprintf("%.0f%%", pct_p[grade_order]), cex=0.85)

hist(results$Neut_nadir, breaks=20, col="#74add1", border="white",
     main="Nadir Neutrophiles (NHP)", xlab="Neut nadir [x 10^9/L]", ylab="Frequence")
abline(v=c(2,1.5,1,0.5), lty=2, col=c("grey60","orange","red","darkred"))
abline(v=median(results$Neut_nadir,na.rm=TRUE), col="#2166ac", lwd=2)

plot(results$Slope_CMP, results$Neut_nadir, pch=16, cex=0.5,
     col=grade_cols[results$Grade_Neut],
     xlab="Slope_CMP", ylab="Neut nadir [x 10^9/L]",
     main="Neut nadir vs Slope_CMP (NHP)")
abline(h=c(2,1.5,1,0.5), lty=2, col=c("grey60","orange","red","darkred"))
legend("topright", names(grade_cols), pch=16, col=grade_cols, bty="n", cex=0.8)

pct_drop <- 100*(results$RBC_nadir - pars_typ$RBC0)/pars_typ$RBC0
hist(pct_drop, breaks=20, col="#f4a582", border="white",
     main="Chute RBC % (NHP)", xlab="%Delta RBC", ylab="Frequence")
abline(v=median(pct_drop,na.rm=TRUE), col="#d6604d", lwd=2)

dev.off()
cat("  -> results/BPA_NHP_simu2_VPC.pdf\n")
