############################################################
# run_pkpd_tdxd_nhp.R
# Simulation PK T-DXd — SINGE CYNOMOLGUS (NHP)
#
# Validation vs FDA BLA 761139 Table 7
#   "3-Month Intermittent IV Dose Toxicity Study in Cynomolgus Monkeys"
#   Doses Q3W : 3, 10, 30 mg/kg
#
# ODE : ADC 2-cpt + DXd 1-cpt + Carboplatin 2-cpt (NHP fit)
############################################################
library(deSolve)

source("parameters_tdxd_nhp.R")   # tdxd_nhp, carbo_nhp, fda_tk_nhp,
                                   # make_nhp_infusion, make_nhp_carbo_infusion,
                                   # tdxd_nhp_state0

if (!dir.exists("results_TDXD")) dir.create("results_TDXD")

# ════════════════════════════════════════════════════════
# ODE PK — ADC 2-cpt + DXd 1-cpt + Carboplatin 2-cpt
# ════════════════════════════════════════════════════════
pk_nhp_ode <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    rate_in       <- if (!is.null(pars$rate_fun))       pars$rate_fun(time)       else 0
    rate_in_carbo <- if (!is.null(pars$rate_fun_carbo)) pars$rate_fun_carbo(time) else 0

    # Krel cycle-dépendant (Yin 2020)
    cycle_num <- max(1L, floor(time / interval_h) + 1L)
    krel_t    <- k_rel_c1 * cycle_num^krel_power *
                 ifelse(cycle_num > 1L, krel_factor, 1.0)

    # ADC — 2 compartiments + élimination TMDD (Michaelis-Menten)
    dC_ADC1 <- rate_in / V1_ADC +
               (Q_ADC / V2_ADC) * C_ADC2 -
               (CL_lin / V1_ADC + Q_ADC / V1_ADC + k_int) * C_ADC1 -
               Vmax_MM * C_ADC1 / (Km_MM + C_ADC1)

    dC_ADC2 <- (Q_ADC / V1_ADC) * C_ADC1 -
               (Q_ADC / V2_ADC) * C_ADC2

    # DXd plasma
    release   <- mass_frac_DXd * krel_t * C_ADC1 * V1_ADC / V_DXd
    flux_in   <- k_inD  * C_DXd
    flux_out  <- k_effD * C_DXd_ic * (V_ic / V_DXd)
    dC_DXd    <- release - (CL_DXd / V_DXd) * C_DXd - flux_in + flux_out

    # DXd intracellulaire
    dC_DXd_ic <- k_inD * C_DXd * (V_DXd / V_ic) - k_effD * C_DXd_ic

    # Carboplatin — 2 compartiments (fit rxode2 NHP, nca_analysis.R)
    Q_c <- if (is.na(Q_carbo) || is.null(Q_carbo)) 0 else Q_carbo
    dC_carbo1 <- rate_in_carbo / V1_carbo -
                 (CL_carbo / V1_carbo + Q_c / V1_carbo) * C_carbo1 +
                 (Q_c / V2_carbo) * C_carbo2
    dC_carbo2 <- (Q_c / V1_carbo) * C_carbo1 - (Q_c / V2_carbo) * C_carbo2

    # Damage combiné : DXd (Emax, saturable) + carboplatin (linéaire, Fornari)
    fu_t         <- fu0_carbo + (fu_inf_carbo - fu0_carbo) * exp(-k_bind_carbo * time)
    C_carbo_free <- fu_t * C_carbo1 * mgL_to_uM_carbo   # µM platine libre
    E_DXd        <- (C_DXd_ic * mgL_to_uM_DXd) /
                    (IC50_DXd_uM + C_DXd_ic * mgL_to_uM_DXd)
    dDamage      <- k_dam_DXd   * E_DXd +
                    k_dam_carbo * C_carbo_free -
                    k_rep       * Damage

    list(c(dC_ADC1, dC_ADC2, dC_DXd, dC_DXd_ic, dC_carbo1, dC_carbo2, dDamage))
  })
}

# ── Paramètres carboplatin dans une liste pars ───────────
add_carbo_pars <- function(pars) {
  pars$CL_carbo        <- carbo_nhp$CL
  pars$V1_carbo        <- carbo_nhp$V1
  pars$Q_carbo         <- carbo_nhp$Q
  pars$V2_carbo        <- carbo_nhp$V2
  pars$fu0_carbo       <- carbo_nhp$fu0
  pars$fu_inf_carbo    <- carbo_nhp$fu_inf
  pars$k_bind_carbo    <- carbo_nhp$k_bind
  pars$k_dam_DXd       <- pars$k_dam
  pars$k_dam_carbo     <- carbo_nhp$k_dam
  pars$mgL_to_uM_carbo <- carbo_nhp$mgL_to_uM
  pars$rate_fun_carbo  <- NULL   # pas de carboplatin par défaut
  pars
}

# ── Simulation helper ────────────────────────────────────
simulate_pk_nhp <- function(dose_mgkg, n_cycles = 3,
                             Tinfu_h = 0.5, BW_kg = 4.0) {
  pars <- add_carbo_pars(tdxd_nhp)
  pars$rate_fun <- make_nhp_infusion(
    dose_mgkg  = dose_mgkg,
    BW_kg      = BW_kg,
    Tinfu_h    = Tinfu_h,
    interval_h = tdxd_nhp$interval_h,
    n_cycles   = n_cycles
  )

  times <- seq(0, n_cycles * 21 * 24, by = 1)

  sol <- as.data.frame(ode(
    y      = tdxd_nhp_state0,
    times  = times,
    func   = pk_nhp_ode,
    parms  = pars,
    method = "lsoda"
  ))
  sol$time_d       <- sol$time / 24
  sol$C_DXd_ngmL   <- sol$C_DXd * 1e3
  sol$C_DXd_ic_uM  <- sol$C_DXd_ic * tdxd_nhp$mgL_to_uM_DXd
  sol
}

# ════════════════════════════════════════════════════════
# Calibration TMDD via optim() — Michaelis-Menten
# ════════════════════════════════════════════════════════
cat("\nGrid search TMDD (2-cpt + Michaelis-Menten) ...\n")

doses_cal <- c(3, 10, 30)

sim_one_tmdd <- function(dose_mgkg, CL_lin, Vmax_MM, Km_MM) {
  p <- add_carbo_pars(tdxd_nhp)
  p$CL_lin  <- CL_lin
  p$Vmax_MM <- Vmax_MM
  p$Km_MM   <- Km_MM
  p$rate_fun <- make_nhp_infusion(dose_mgkg = dose_mgkg, BW_kg = 4.0,
                                   Tinfu_h = 0.5, interval_h = NULL, n_cycles = 1)
  times <- c(seq(0, 2, by = 0.1), seq(3, 504, by = 1))
  sol <- tryCatch(
    suppressWarnings(as.data.frame(ode(
      y = tdxd_nhp_state0, times = times,
      func = pk_nhp_ode, parms = p, method = "lsoda"))),
    error = function(e) NULL)
  if (is.null(sol) || any(is.nan(sol$C_ADC1)) || min(sol$C_ADC1) < -1e-6)
    return(NULL)
  sol$C_ADC1 <- pmax(sol$C_ADC1, 1e-12)
  sol
}

nca_tmdd <- function(sol, fda) {
  if (is.null(sol)) return(list(C0=NA, AUC=NA, t12=NA, ok=FALSE))
  C0  <- max(sol$C_ADC1[sol$time <= 1])
  idx <- sol$time <= 504
  AUC <- sum(diff(sol$time[idx]) *
             (sol$C_ADC1[idx][-sum(idx)] + sol$C_ADC1[idx][-1]) / 2) / 24
  idt <- sol$time >= 100 & sol$time <= 480 & sol$C_ADC1 > 0
  if (sum(idt) < 5) return(list(C0=C0, AUC=AUC, t12=NA, ok=FALSE))
  lm_f <- tryCatch(lm(log(C_ADC1) ~ time, data = sol[idt, ]), error=function(e) NULL)
  if (is.null(lm_f) || coef(lm_f)[2] >= 0) return(list(C0=C0, AUC=AUC, t12=NA, ok=FALSE))
  t12 <- log(2) / (-coef(lm_f)[2]) / 24
  list(C0=C0, AUC=AUC, t12=t12, ok=TRUE)
}

w_AUC <- c(2, 2, 1)
w_t12 <- c(1, 1, 2)

wrss_tmdd <- function(CL_lin, Vmax_MM, Km_MM) {
  total <- 0
  for (i in seq_along(doses_cal)) {
    fda <- fda_tk_nhp[[i]]
    sol <- sim_one_tmdd(doses_cal[i], CL_lin, Vmax_MM, Km_MM)
    nca <- nca_tmdd(sol, fda)
    if (!nca$ok || nca$C0 <= 0 || nca$AUC <= 0 || nca$t12 <= 0) return(1e8)
    total <- total +
      (log(nca$C0  / fda$C0_ADC))^2 +
      w_AUC[i] * (log(nca$AUC / fda$AUC21d_ADC))^2 +
      w_t12[i] * (log(nca$t12 / fda$t_half_d))^2
  }
  total
}

CL_grid <- exp(seq(log(5e-5),  log(2e-3),  length.out = 9))
VM_grid <- exp(seq(log(0.05),  log(6.0),   length.out = 9))
Km_grid <- exp(seq(log(50.0),  log(2000.0), length.out = 11))

best_val <- 1e8; best_par <- c(1.115e-3, 0.200, 400.0)
for (cl in CL_grid) for (vm in VM_grid) for (km in Km_grid) {
  v <- wrss_tmdd(cl, vm, km)
  if (v < best_val) { best_val <- v; best_par <- c(cl, vm, km) }
}
cat(sprintf("Coarse: CL=%.3e Vmax=%.3f Km=%.1f  RMSE=%.1f%%\n",
            best_par[1], best_par[2], best_par[3], 100*sqrt(best_val/9)))

CL_f <- exp(seq(log(best_par[1]*0.7), log(best_par[1]*1.4), length.out = 7))
VM_f <- exp(seq(log(best_par[2]*0.5), log(best_par[2]*2.0), length.out = 7))
Km_f <- exp(seq(log(best_par[3]*0.4), log(best_par[3]*2.5), length.out = 7))

best_val2 <- best_val; best_par2 <- best_par
for (cl in CL_f) for (vm in VM_f) for (km in Km_f) {
  v <- wrss_tmdd(cl, vm, km)
  if (v < best_val2) { best_val2 <- v; best_par2 <- c(cl, vm, km) }
}
cat(sprintf("Fine:   CL=%.4e Vmax=%.4f Km=%.2f  RMSE=%.1f%%\n",
            best_par2[1], best_par2[2], best_par2[3], 100*sqrt(best_val2/9)))

obj_tmdd <- function(theta) {
  CL_lin  <- exp(theta[1])
  Vmax_MM <- exp(theta[2])
  Km_MM   <- exp(theta[3])
  if (any(c(CL_lin, Vmax_MM, Km_MM) <= 0)) return(1e8)
  wrss_tmdd(CL_lin, Vmax_MM, Km_MM)
}

opt <- optim(log(best_par2), obj_tmdd, method = "Nelder-Mead",
             control = list(maxit = 5000, reltol = 1e-10))

CL_lin_cal  <- exp(opt$par[1])
Vmax_MM_cal <- exp(opt$par[2])
Km_MM_cal   <- exp(opt$par[3])
rmse_cal    <- 100 * sqrt(opt$value / 9)

tdxd_nhp$CL_lin  <- CL_lin_cal
tdxd_nhp$Vmax_MM <- Vmax_MM_cal
tdxd_nhp$Km_MM   <- Km_MM_cal

cat("\n═══ VALIDATION FINALE — 2-cpt + TMDD ═══\n")
for (i in seq_along(doses_cal)) {
  fda <- fda_tk_nhp[[i]]
  sol <- sim_one_tmdd(doses_cal[i], CL_lin_cal, Vmax_MM_cal, Km_MM_cal)
  nca <- nca_tmdd(sol, fda)
  cat(sprintf("  %2d mg/kg: C0 %.1f/%.1f (%+.1f%%)  AUC %d/%d (%+.1f%%)  T½ %.2f/%.2f (%+.1f%%)\n",
              doses_cal[i],
              nca$C0,  fda$C0_ADC,      100*(nca$C0  / fda$C0_ADC      - 1),
              round(nca$AUC), round(fda$AUC21d_ADC), 100*(nca$AUC / fda$AUC21d_ADC - 1),
              nca$t12, fda$t_half_d,    100*(nca$t12 / fda$t_half_d    - 1)))
}
cat(sprintf("\n  RMSE(C0,AUC,T½) = %.1f%%\n", rmse_cal))
cat(sprintf("    CL_lin  <- %.5e  # L/h\n",   CL_lin_cal))
cat(sprintf("    Vmax_MM <- %.5f  # mg/L/h\n", Vmax_MM_cal))
cat(sprintf("    Km_MM   <- %.2f      # mg/L\n\n", Km_MM_cal))

# ════════════════════════════════════════════════════════
# Simulations doses FDA Table 7 : 3, 10, 30 mg/kg Q3W × 3
# ════════════════════════════════════════════════════════
cat("\n── Simulation Q3W × 3 cycles ──────────────────────────────\n")
doses <- c(3, 10, 30)
sims  <- lapply(doses, function(d) simulate_pk_nhp(d, n_cycles = 3))

# ════════════════════════════════════════════════════════
# Graphiques PK — ADC + DXd par dose
# ════════════════════════════════════════════════════════
dose_cols <- c("#2166ac", "#4dac26", "#d6604d")
dose_days <- c(0, 21, 42)

pdf("results_TDXD/NHP_PK_validation_Table7.pdf", width = 14, height = 9)
par(mfrow = c(2, 3), mar = c(4, 4.5, 3, 1.5))

for (i in seq_along(doses)) {
  s   <- sims[[i]]
  col <- dose_cols[i]
  d   <- doses[i]

  plot(s$time_d, s$C_ADC1,
       type = "l", lwd = 2.5, col = col, log = "y",
       xlab = "Temps (jours)", ylab = "ADC [µg/mL]",
       main = sprintf("ADC sérum — %d mg/kg Q3W", d),
       ylim = c(1, max(s$C_ADC1) * 2))
  abline(v = dose_days, lty = 2, col = "grey60", lwd = 0.8)
  fda  <- fda_tk_nhp[[i]]
  points(0.02, fda$C0_ADC, pch = 19, cex = 1.5, col = "black")
  legend("topright", c("Simulation", "C0 Table 7"),
         col = c(col, "black"), lwd = c(2.5, NA), pch = c(NA, 19),
         bty = "n", cex = 0.85)

  plot(s$time_d, s$C_DXd_ngmL,
       type = "l", lwd = 2.5, col = col,
       xlab = "Temps (jours)", ylab = "DXd [ng/mL]",
       main = sprintf("DXd plasma — %d mg/kg Q3W", d))
  abline(v = dose_days, lty = 2, col = "grey60", lwd = 0.8)
  points(0.02, fda$C0_DXd_ng, pch = 17, cex = 1.5, col = "black")
  legend("topright", c("Simulation", "C0 Table 7"),
         col = c(col, "black"), lwd = c(2.5, NA), pch = c(NA, 17),
         bty = "n", cex = 0.85)
}
dev.off()
cat("  -> results_TDXD/NHP_PK_validation_Table7.pdf\n")

pdf("results_TDXD/NHP_PK_3doses_comparison.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1.5))

y_min <- min(sapply(sims, function(s) min(s$C_ADC1[s$C_ADC1 > 0])))
y_max <- max(sapply(sims, function(s) max(s$C_ADC1)))
plot(sims[[1]]$time_d, sims[[1]]$C_ADC1, type = "n", log = "y",
     ylim = c(y_min, y_max * 3), xlab = "Temps (jours)", ylab = "ADC [µg/mL]",
     main = "ADC sérum — NHP Q3W × 3")
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$C_ADC1, col = dose_cols[i], lwd = 2)
abline(v = dose_days, lty = 2, col = "grey70")
for (i in seq_along(doses))
  points(0.02, fda_tk_nhp[[i]]$C0_ADC, pch = 19, col = dose_cols[i], cex = 1.5)
legend("topright", paste0(doses, " mg/kg"), col = dose_cols, lwd = 2, bty = "n", cex = 0.9)

y_max_dxd <- max(sapply(sims, function(s) max(s$C_DXd_ngmL)))
plot(sims[[1]]$time_d, sims[[1]]$C_DXd_ngmL, type = "n",
     ylim = c(0, y_max_dxd * 1.15), xlab = "Temps (jours)", ylab = "DXd [ng/mL]",
     main = "DXd plasma — NHP Q3W × 3")
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$C_DXd_ngmL, col = dose_cols[i], lwd = 2)
abline(v = dose_days, lty = 2, col = "grey70")
for (i in seq_along(doses))
  points(0.02, fda_tk_nhp[[i]]$C0_DXd_ng, pch = 17, col = dose_cols[i], cex = 1.5)
legend("topright", paste0(doses, " mg/kg"), col = dose_cols, lwd = 2, bty = "n", cex = 0.9)
dev.off()
cat("  -> results_TDXD/NHP_PK_3doses_comparison.pdf\n")

# ════════════════════════════════════════════════════════
# Validation NCA
# ════════════════════════════════════════════════════════
cat("\n═══ VALIDATION NCA — NHP Q3W × 3 (Day 1) ═══\n")
cat(sprintf("  %-10s │ %-7s %-7s %-6s │ %-9s %-9s %-6s │ %-8s %-8s\n",
            "Dose","C0_obs","C0_sim","ratio","AUC_obs","AUC_sim","ratio","T½_obs","T½_sim"))
cat(sprintf("  %s\n", paste(rep("─", 80), collapse="")))

for (i in seq_along(doses)) {
  s   <- sims[[i]]
  fda <- fda_tk_nhp[[i]]
  C0_sim  <- max(s$C_ADC1[s$time <= 24])
  idx     <- s$time <= 504
  AUC_sim <- sum(diff(s$time[idx]) *
                 (s$C_ADC1[idx][-sum(idx)] + s$C_ADC1[idx][-1]) / 2) / 24
  idx_t   <- s$time >= 100 & s$time <= 480 & s$C_ADC1 > 0
  lm_fit  <- lm(log(C_ADC1) ~ time, data = s[idx_t, ])
  t12_sim <- log(2) / abs(coef(lm_fit)[2]) / 24
  cat(sprintf("  %-10s │ %-7.1f %-7.1f %-6.3f │ %-9.0f %-9.0f %-6.3f │ %-8.2f %-8.2f\n",
              paste0(doses[i], " mg/kg"),
              fda$C0_ADC, C0_sim, C0_sim/fda$C0_ADC,
              fda$AUC21d_ADC, AUC_sim, AUC_sim/fda$AUC21d_ADC,
              fda$t_half_d, t12_sim))
}
cat("═══════════════════════════════════════════════════════════\n\n")
cat("Fichiers générés dans results_TDXD/\n")
