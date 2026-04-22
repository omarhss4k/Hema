############################################################
# run_pkpd_tdxd_nhp.R
# Simulation PK T-DXd — SINGE CYNOMOLGUS (NHP)
#
# Validation vs FDA BLA 761139 Table 7
#   "3-Month Intermittent IV Dose Toxicity Study in Cynomolgus Monkeys"
#   Doses Q3W : 3, 10, 30 mg/kg
#
# PK uniquement (ADC + DXd) — sans PD hématologique
# ODE : réutilise la structure de pkpd_tdxd_rat.R
############################################################
library(deSolve)

# ── Chargement ───────────────────────────────────────────
source("parameters_tdxd_nhp.R")      # tdxd_nhp, fda_tk_nhp, make_nhp_infusion

if (!dir.exists("results_TDXD")) dir.create("results_TDXD")

# ════════════════════════════════════════════════════════
# ODE PK — ADC 2-cpt + DXd 1-cpt  (pas de PD Fornari)
# ════════════════════════════════════════════════════════
pk_nhp_ode <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0

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

    # Damage (non utilisé en PK-only mais conservé pour compatibilité)
    dDamage   <- k_dam * (C_DXd_ic * mgL_to_uM_DXd) /
                 (IC50_DXd_uM + C_DXd_ic * mgL_to_uM_DXd) - k_rep * Damage

    list(c(dC_ADC1, dC_ADC2, dC_DXd, dC_DXd_ic, dDamage))
  })
}

# ── Simulation helper ────────────────────────────────────
simulate_pk_nhp <- function(dose_mgkg, n_cycles = 3,
                             Tinfu_h = 0.5, BW_kg = 4.0) {
  pars <- tdxd_nhp
  pars$rate_fun <- make_nhp_infusion(
    dose_mgkg  = dose_mgkg,
    BW_kg      = BW_kg,
    Tinfu_h    = Tinfu_h,
    interval_h = tdxd_nhp$interval_h,
    n_cycles   = n_cycles
  )

  times <- seq(0, n_cycles * 21 * 24, by = 1)   # pas = 1 h

  sol <- as.data.frame(ode(
    y     = tdxd_nhp_state0,
    times = times,
    func  = pk_nhp_ode,
    parms = pars,
    method = "lsoda"
  ))
  sol$time_d     <- sol$time / 24
  sol$C_DXd_ngmL <- sol$C_DXd * 1e3      # mg/L → ng/mL
  sol$C_DXd_ic_uM <- sol$C_DXd_ic * tdxd_nhp$mgL_to_uM_DXd
  sol
}

# ════════════════════════════════════════════════════════
# Calibration TMDD via optim() — Michaelis-Menten non-linéaire
# Paramètres : CL_lin [L/h], Vmax_MM [mg/L/h], Km_MM [mg/L]
# Objectif   : WRSS log-relatif sur C0, AUC0-21d, T½ × 3 doses
# ════════════════════════════════════════════════════════
cat("\nGrid search TMDD (2-cpt + Michaelis-Menten) ...\n")

doses_cal <- c(3, 10, 30)

# ── Grille grossière : CL × Vmax × Km ────────────────────────────────────────
sim_one_tmdd <- function(dose_mgkg, CL_lin, Vmax_MM, Km_MM,
                          V2 = tdxd_nhp$V2_ADC) {
  p <- tdxd_nhp
  p$CL_lin  <- CL_lin
  p$Vmax_MM <- Vmax_MM
  p$Km_MM   <- Km_MM
  p$V2_ADC  <- V2
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
  # t>=200h : élimine la contamination phase-alpha (T½_α≈31h → 1.2% restant à 200h)
  idt <- sol$time >= 200 & sol$time <= 480 & sol$C_ADC1 > 0
  if (sum(idt) < 5) return(list(C0=C0, AUC=AUC, t12=NA, ok=FALSE))
  lm_f <- tryCatch(lm(log(C_ADC1) ~ time, data = sol[idt, ]), error=function(e) NULL)
  if (is.null(lm_f) || coef(lm_f)[2] >= 0) return(list(C0=C0, AUC=AUC, t12=NA, ok=FALSE))
  t12 <- log(2) / (-coef(lm_f)[2]) / 24
  list(C0=C0, AUC=AUC, t12=t12, ok=TRUE)
}

wrss_tmdd <- function(CL_lin, Vmax_MM, Km_MM, V2 = tdxd_nhp$V2_ADC) {
  total <- 0
  for (i in seq_along(doses_cal)) {
    fda <- fda_tk_nhp[[i]]
    sol <- sim_one_tmdd(doses_cal[i], CL_lin, Vmax_MM, Km_MM, V2)
    nca <- nca_tmdd(sol, fda)
    if (!nca$ok || nca$C0 <= 0 || nca$AUC <= 0 || nca$t12 <= 0) return(1e8)
    total <- total +
      1.0 * (log(nca$C0  / fda$C0_ADC))^2 +
      1.0 * (log(nca$AUC / fda$AUC21d_ADC))^2 +
      2.0 * (log(nca$t12 / fda$t_half_d))^2   # T½ x2 : signature clé TMDD
  }
  total
}

# ── Grille élargie : Vmax jusqu'à 2 mg/L/h, Km jusqu'à 300 µg/mL ─────────────
# La physique TMDD implique Km dans la gamme des conc. terminales (~10-100 µg/mL)
# et Vmax ~0.1-1.0 mg/L/h pour expliquer la variation CL × 1.5 entre doses
CL_grid <- exp(seq(log(1e-4),  log(2e-3),  length.out = 10))
VM_grid <- exp(seq(log(0.02),  log(2.0),   length.out = 10))
Km_grid <- exp(seq(log(1.0),   log(300.0), length.out = 10))

best_val <- 1e8; best_par <- c(1.115e-3, 0.060, 4.0)
for (cl in CL_grid) for (vm in VM_grid) for (km in Km_grid) {
  v <- wrss_tmdd(cl, vm, km)
  if (v < best_val) { best_val <- v; best_par <- c(cl, vm, km) }
}
cat(sprintf("Coarse: CL=%.3e Vmax=%.3f Km=%.1f  RMSE=%.1f%%\n",
            best_par[1], best_par[2], best_par[3], 100*sqrt(best_val/12)))

# Grille fine : ×0.4–×2.5 autour du meilleur
CL_f <- exp(seq(log(best_par[1]*0.4), log(best_par[1]*2.5), length.out = 8))
VM_f <- exp(seq(log(best_par[2]*0.3), log(best_par[2]*3.0), length.out = 8))
Km_f <- exp(seq(log(best_par[3]*0.3), log(best_par[3]*3.0), length.out = 8))

best_val2 <- best_val; best_par2 <- best_par
for (cl in CL_f) for (vm in VM_f) for (km in Km_f) {
  v <- wrss_tmdd(cl, vm, km)
  if (v < best_val2) { best_val2 <- v; best_par2 <- c(cl, vm, km) }
}
cat(sprintf("Fine:   CL=%.4e Vmax=%.4f Km=%.2f  RMSE=%.1f%%\n",
            best_par2[1], best_par2[2], best_par2[3], 100*sqrt(best_val2/12)))

# ── Raffinement Nelder-Mead 4D : CL_lin, Vmax_MM, Km_MM, V2_ADC ─────────────
# V2 libre : V2 fixé sur T½ linéaire ≠ T½ TMDD optimal
obj_tmdd <- function(theta) {
  CL_lin  <- exp(theta[1])
  Vmax_MM <- exp(theta[2])
  Km_MM   <- exp(theta[3])
  V2      <- exp(theta[4])
  if (any(c(CL_lin, Vmax_MM, Km_MM, V2) <= 0)) return(1e8)
  if (V2 < 0.01 || V2 > 0.60) return(1e8)   # garde-fou physique
  wrss_tmdd(CL_lin, Vmax_MM, Km_MM, V2)
}

theta0_4d <- c(log(best_par2), log(tdxd_nhp$V2_ADC))

opt <- optim(theta0_4d, obj_tmdd, method = "Nelder-Mead",
             control = list(maxit = 8000, reltol = 1e-10))

CL_lin_cal  <- exp(opt$par[1])
Vmax_MM_cal <- exp(opt$par[2])
Km_MM_cal   <- exp(opt$par[3])
V2_cal      <- exp(opt$par[4])
rmse_cal    <- 100 * sqrt(opt$value / 12)   # 12 = (C0+AUC+T½×2) × 3 doses

# Mise à jour des paramètres dans tdxd_nhp
tdxd_nhp$CL_lin  <- CL_lin_cal
tdxd_nhp$Vmax_MM <- Vmax_MM_cal
tdxd_nhp$Km_MM   <- Km_MM_cal
tdxd_nhp$V2_ADC  <- V2_cal

# ── Validation finale ─────────────────────────────────────────────────────────
cat("\n═══ VALIDATION FINALE — 2-cpt + TMDD ═══\n")
for (i in seq_along(doses_cal)) {
  fda <- fda_tk_nhp[[i]]
  sol <- sim_one_tmdd(doses_cal[i], CL_lin_cal, Vmax_MM_cal, Km_MM_cal, V2_cal)
  nca <- nca_tmdd(sol, fda)
  cat(sprintf("  %2d mg/kg: C0 %.1f/%.1f (%+.1f%%)  AUC %d/%d (%+.1f%%)  T½ %.2f/%.2f (%+.1f%%)\n",
              doses_cal[i],
              nca$C0,  fda$C0_ADC,      100*(nca$C0  / fda$C0_ADC      - 1),
              round(nca$AUC), round(fda$AUC21d_ADC), 100*(nca$AUC / fda$AUC21d_ADC - 1),
              nca$t12, fda$t_half_d,    100*(nca$t12 / fda$t_half_d    - 1)))
}
cat(sprintf("\n  RMSE(C0,AUC,T½) = %.1f%%\n", rmse_cal))
cat("\n  → Paramètres TMDD :\n")
cat(sprintf("    CL_lin  <- %.5e  # L/h\n",   CL_lin_cal))
cat(sprintf("    Vmax_MM <- %.5f  # mg/L/h\n", Vmax_MM_cal))
cat(sprintf("    Km_MM   <- %.2f      # mg/L (µg/mL)\n",   Km_MM_cal))
cat(sprintf("    V2_ADC  <- %.5f  # L  (calibré TMDD)\n\n", V2_cal))

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

# ── Graphiques individuels par dose ──────────────────────
for (i in seq_along(doses)) {
  s   <- sims[[i]]
  col <- dose_cols[i]
  d   <- doses[i]

  # ADC sérum — log scale
  plot(s$time_d, s$C_ADC1,
       type = "l", lwd = 2.5, col = col, log = "y",
       xlab = "Temps (jours)", ylab = "ADC [µg/mL]",
       main = sprintf("ADC sérum — %d mg/kg Q3W", d),
       ylim = c(1, max(s$C_ADC1) * 2))
  abline(v = dose_days, lty = 2, col = "grey60", lwd = 0.8)

  # Point C0 FDA (Day 1)
  fda  <- fda_tk_nhp[[i]]
  points(0.02, fda$C0_ADC, pch = 19, cex = 1.5, col = "black")
  legend("topright", c("Simulation", "C0 Table 7"),
         col = c(col, "black"), lwd = c(2.5, NA), pch = c(NA, 19),
         bty = "n", cex = 0.85)

  # DXd plasmatique
  plot(s$time_d, s$C_DXd_ngmL,
       type = "l", lwd = 2.5, col = col,
       xlab = "Temps (jours)", ylab = "DXd [ng/mL]",
       main = sprintf("DXd plasma — %d mg/kg Q3W", d))
  abline(v = dose_days, lty = 2, col = "grey60", lwd = 0.8)
  # Point C0_DXd FDA
  points(0.02, fda$C0_DXd_ng, pch = 17, cex = 1.5, col = "black")
  legend("topright", c("Simulation", "C0 Table 7"),
         col = c(col, "black"), lwd = c(2.5, NA), pch = c(NA, 17),
         bty = "n", cex = 0.85)
}

dev.off()
cat("  -> results_TDXD/NHP_PK_validation_Table7.pdf\n")

# ── Graphique superposé 3 doses ──────────────────────────
pdf("results_TDXD/NHP_PK_3doses_comparison.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1.5))

# ADC — 3 doses superposées (log)
y_min <- min(sapply(sims, function(s) min(s$C_ADC1[s$C_ADC1 > 0])))
y_max <- max(sapply(sims, function(s) max(s$C_ADC1)))
plot(sims[[1]]$time_d, sims[[1]]$C_ADC1,
     type = "n", log = "y",
     ylim = c(y_min, y_max * 3),
     xlab = "Temps (jours)", ylab = "ADC [µg/mL]",
     main = "ADC sérum — NHP Q3W × 3")
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$C_ADC1, col = dose_cols[i], lwd = 2)
abline(v = dose_days, lty = 2, col = "grey70")
# Points C0 Table 7
for (i in seq_along(doses))
  points(0.02, fda_tk_nhp[[i]]$C0_ADC, pch = 19, col = dose_cols[i], cex = 1.5)
legend("topright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.9)
text(0.5, y_max * 2.5, "● = C0 Table 7", cex = 0.8, col = "grey40")

# DXd — 3 doses superposées
y_max_dxd <- max(sapply(sims, function(s) max(s$C_DXd_ngmL)))
plot(sims[[1]]$time_d, sims[[1]]$C_DXd_ngmL,
     type = "n",
     ylim = c(0, y_max_dxd * 1.15),
     xlab = "Temps (jours)", ylab = "DXd [ng/mL]",
     main = "DXd plasma — NHP Q3W × 3")
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$C_DXd_ngmL, col = dose_cols[i], lwd = 2)
abline(v = dose_days, lty = 2, col = "grey70")
for (i in seq_along(doses))
  points(0.02, fda_tk_nhp[[i]]$C0_DXd_ng, pch = 17, col = dose_cols[i], cex = 1.5)
legend("topright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.9)

dev.off()
cat("  -> results_TDXD/NHP_PK_3doses_comparison.pdf\n")

# ════════════════════════════════════════════════════════
# Validation NCA : C0, AUC0-21d, T½ simulés vs Table 7
# ════════════════════════════════════════════════════════
cat("\n═══════════════════════════════════════════════════════════\n")
cat("  VALIDATION NCA — NHP Q3W × 3 (Day 1 simulé vs Table 7)\n")
cat("───────────────────────────────────────────────────────────\n")
cat(sprintf("  %-8s │ %-7s %-7s %-6s │ %-9s %-9s %-6s │ %-8s %-8s\n",
            "Dose", "C0_obs", "C0_sim", "ratio",
            "AUC_obs", "AUC_sim", "ratio",
            "T½_obs", "T½_sim"))
cat(sprintf("  %-8s │ %-7s %-7s %-6s │ %-9s %-9s %-6s │ %-8s %-8s\n",
            "mg/kg", "µg/mL","µg/mL","",
            "µg.d/mL","µg.d/mL","",
            "jours","jours"))
cat(sprintf("  %s\n", paste(rep("─", 85), collapse="")))

for (i in seq_along(doses)) {
  s    <- sims[[i]]
  fda  <- fda_tk_nhp[[i]]

  # C0 simulé = Cmax dans la première heure
  C0_sim  <- max(s$C_ADC1[s$time <= 24])

  # AUC0-21d simulé par règle trapézoïdale
  idx     <- s$time <= 504
  AUC_sim <- sum(diff(s$time[idx]) *
                 (s$C_ADC1[idx][-sum(idx)] + s$C_ADC1[idx][-1]) / 2) / 24

  # T½ simulé depuis la pente terminale (t > 21j après dose 1)
  idx_t    <- s$time >= 100 & s$time <= 480
  if (sum(idx_t) > 5) {
    lm_fit  <- lm(log(C_ADC1) ~ time, data = s[idx_t & s$C_ADC1 > 0, ])
    beta_sim <- abs(coef(lm_fit)[2])
    t_half_sim <- log(2) / beta_sim / 24   # h → j
  } else {
    t_half_sim <- NA
  }

  cat(sprintf("  %-8s │ %-7.1f %-7.1f %-6.3f │ %-9.0f %-9.0f %-6.3f │ %-8.2f %-8.2f\n",
              paste0(doses[i], " mg/kg"),
              fda$C0_ADC, C0_sim, C0_sim / fda$C0_ADC,
              fda$AUC21d_ADC, AUC_sim, AUC_sim / fda$AUC21d_ADC,
              fda$t_half_d, t_half_sim))
}
cat("═══════════════════════════════════════════════════════════\n")
cat("  Ratio ~1.0 = bonne concordance avec FDA Table 7\n")
cat("  T½ dose-dépendant modélisé via TMDD Michaelis-Menten (Vmax/Km calibrés optim)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# ════════════════════════════════════════════════════════
# Tableau comparatif PDF — Simulé vs FDA BLA Table 7
# ════════════════════════════════════════════════════════
nca_rows <- vector("list", length(doses))
for (i in seq_along(doses)) {
  s   <- sims[[i]]
  fda <- fda_tk_nhp[[i]]

  C0_sim  <- max(s$C_ADC1[s$time <= 24])
  idx     <- s$time <= 504
  AUC_sim <- sum(diff(s$time[idx]) *
                 (s$C_ADC1[idx][-sum(idx)] + s$C_ADC1[idx][-1]) / 2) / 24
  idx_t   <- s$time >= 100 & s$time <= 480
  lm_fit  <- lm(log(C_ADC1) ~ time, data = s[idx_t & s$C_ADC1 > 0, ])
  t12_sim <- log(2) / abs(coef(lm_fit)[2]) / 24

  nca_rows[[i]] <- list(
    dose    = doses[i],
    C0_obs  = fda$C0_ADC,
    C0_sim  = C0_sim,
    C0_r    = C0_sim / fda$C0_ADC,
    AUC_obs = fda$AUC21d_ADC,
    AUC_sim = AUC_sim,
    AUC_r   = AUC_sim / fda$AUC21d_ADC,
    t12_obs = fda$t_half_d,
    t12_sim = t12_sim,
    t12_r   = t12_sim / fda$t_half_d
  )
}

# ── helper : couleur ratio ────────────────────────────
ratio_col <- function(r) {
  err <- abs(r - 1)
  if      (err <= 0.10) "#2ca02c"   # vert   ≤ 10 %
  else if (err <= 0.25) "#ff7f0e"   # orange 10–25 %
  else                  "#d62728"   # rouge  > 25 %
}

pdf("results_TDXD/NHP_NCA_Table_FDA_vs_Sim.pdf", width = 14, height = 7)

# En-tête
par(mar = c(1, 1, 3, 1), bg = "white")
plot.new()
mtext("T-DXd NHP — Simulation vs FDA BLA 761139 Table 7 (Day 1, M+F mean, Q3W)",
      side = 3, line = 1, cex = 1.3, font = 2)
mtext("Doses : 3 / 10 / 30 mg/kg   |   Modèle : ADC 2-cpt + TMDD Michaelis-Menten, k_int = 0",
      side = 3, line = -0.2, cex = 0.85, col = "grey30")

# ── coordonnées tableau ───────────────────────────────
# colonnes : dose | C0_obs C0_sim C0_r | AUC_obs AUC_sim AUC_r | T½_obs T½_sim T½_r
col_x  <- c(0.05, 0.14, 0.22, 0.30,  0.40, 0.50, 0.58,  0.68, 0.77, 0.86)
row_y  <- c(0.83, 0.65, 0.45, 0.25)   # en-têtes + 3 doses

hdrs1  <- c("Dose\n(mg/kg)",
            "C0 FDA\n(µg/mL)", "C0 sim\n(µg/mL)", "ratio",
            "AUC FDA\n(µg·d/mL)", "AUC sim\n(µg·d/mL)", "ratio",
            "T½ FDA\n(jours)", "T½ sim\n(jours)", "ratio")

# Bandes de fond
rect(0.0, 0.18, 1.0, 0.90, col = "#f7f7f7", border = NA)
for (r in c(0.58, 0.38)) rect(0.0, r - 0.10, 1.0, r + 0.10, col = "white", border = NA)

# Séparateurs de groupes
segments(c(0.325, 0.615), 0.18, c(0.325, 0.615), 0.90,
         col = "#aaaaaa", lwd = 1.2, lty = 2)

# En-têtes de groupes
text(0.215, 0.93, "ADC sérum (DS-8201a)", font = 2, cex = 0.95)
text(0.490, 0.93, "AUC0-21j (DS-8201a)", font = 2, cex = 0.95)
text(0.775, 0.93, "Demi-vie terminale", font = 2, cex = 0.95)

# En-têtes de colonnes
for (j in seq_along(col_x))
  text(col_x[j], row_y[1], hdrs1[j], cex = 0.82, font = 2,
       col = ifelse(j == 1, "black", "grey20"))

# Données
for (i in seq_along(doses)) {
  r  <- nca_rows[[i]]
  yy <- row_y[i + 1]
  bg <- if (i %% 2 == 0) "#eaf3fb" else "white"
  rect(0.0, yy - 0.09, 1.0, yy + 0.09, col = bg, border = NA)

  # valeurs FDA
  text(col_x[1],  yy, sprintf("%d mg/kg", r$dose), cex = 0.9, font = 2)
  text(col_x[2],  yy, sprintf("%.1f",  r$C0_obs),  cex = 0.88)
  text(col_x[3],  yy, sprintf("%.1f",  r$C0_sim),  cex = 0.88)
  text(col_x[4],  yy, sprintf("%.3f",  r$C0_r),    cex = 0.88,
       col = ratio_col(r$C0_r), font = 2)

  text(col_x[5],  yy, sprintf("%.0f",  r$AUC_obs), cex = 0.88)
  text(col_x[6],  yy, sprintf("%.0f",  r$AUC_sim), cex = 0.88)
  text(col_x[7],  yy, sprintf("%.3f",  r$AUC_r),   cex = 0.88,
       col = ratio_col(r$AUC_r), font = 2)

  text(col_x[8],  yy, sprintf("%.2f",  r$t12_obs), cex = 0.88)
  text(col_x[9],  yy, sprintf("%.2f",  r$t12_sim), cex = 0.88)
  text(col_x[10], yy, sprintf("%.3f",  r$t12_r),   cex = 0.88,
       col = ratio_col(r$t12_r), font = 2)
}

# Bordure
rect(0.0, 0.18, 1.0, 0.90, col = NA, border = "#555555", lwd = 1.5)

# Légende couleurs
legend_x <- 0.05; legend_y <- 0.10
text(legend_x, legend_y + 0.04, "Légende ratio sim/FDA :", cex = 0.78, adj = 0)
rect(legend_x,       legend_y - 0.01, legend_x + 0.018, legend_y + 0.03,
     col = "#2ca02c", border = NA)
text(legend_x + 0.024, legend_y + 0.01, "≤ 10 %", cex = 0.75, adj = 0)
rect(legend_x + 0.09,  legend_y - 0.01, legend_x + 0.108, legend_y + 0.03,
     col = "#ff7f0e", border = NA)
text(legend_x + 0.114, legend_y + 0.01, "10–25 %", cex = 0.75, adj = 0)
rect(legend_x + 0.19,  legend_y - 0.01, legend_x + 0.208, legend_y + 0.03,
     col = "#d62728", border = NA)
text(legend_x + 0.214, legend_y + 0.01, "> 25 %", cex = 0.75, adj = 0)

text(0.62, legend_y + 0.01,
     "TMDD Michaelis-Menten : T½ dose-dépendant modélisé (Vmax/Km calibrés par optim)",
     cex = 0.75, col = "grey40", adj = 0)

dev.off()
cat("  -> results_TDXD/NHP_NCA_Table_FDA_vs_Sim.pdf\n")

# ════════════════════════════════════════════════════════
# Tableau CTCAE v5.0 — Grades anémie + neutropénie
# ════════════════════════════════════════════════════════
pdf("results_TDXD/NHP_CTCAE_Grades_Anemie_Neut.pdf", width = 13, height = 8)

# ── couleurs grades ───────────────────────────────────
gc <- c(G1 = "#ffffb2", G2 = "#fecc5c", G3 = "#fd8d3c", G4 = "#e31a1c")
gc_txt <- c(G1 = "black",  G2 = "black",   G3 = "white",   G4 = "white")

par(mar = c(1, 1, 3.5, 1), bg = "white")
plot.new()
mtext("CTCAE v5.0 — Grades hématologiques : Anémie & Neutropénie",
      side = 3, line = 1.8, cex = 1.35, font = 2)
mtext("Référence humain — applicable à l'interprétation des effets T-DXd (DS-8201a)",
      side = 3, line = 0.4, cex = 0.9, col = "grey30")

# ──────────────────────────────────────────────────────
# BLOC ANÉMIE  (hémoglobine)
# ──────────────────────────────────────────────────────
# En-tête bloc
rect(0.02, 0.71, 0.49, 0.90, col = "#2166ac", border = NA)
text(0.255, 0.805, "ANÉMIE  (Hémoglobine / Hgb)", col = "white", cex = 1.1, font = 2)

# Sous-en-têtes
sub_x <- c(0.055, 0.175, 0.355, 0.455)
rect(0.02, 0.63, 0.49, 0.71, col = "#d0e4f5", border = NA)
text(sub_x[1], 0.67, "Grade", font = 2, cex = 0.9)
text(sub_x[2], 0.67, "Hgb (g/dL)", font = 2, cex = 0.9)
text(sub_x[3], 0.67, "Définition clinique", font = 2, cex = 0.9)
text(sub_x[4], 0.67, "NHP ref.", font = 2, cex = 0.9)

anemie <- list(
  list(g="G1", hgb="LLN – 10.0",   def="Asymptomatique",                 nhp="LLN singe ≈ 12 g/dL"),
  list(g="G2", hgb="8.0 – <10.0",  def="Symptômes modérés / fatigue",    nhp=""),
  list(g="G3", hgb="<8.0",         def="Transfusion indiquée",           nhp=""),
  list(g="G4", hgb="<6.5 (typ.)",  def="Pronostic vital — urgence",      nhp="")
)
row_ys_a <- c(0.56, 0.47, 0.38, 0.29)

for (k in seq_along(anemie)) {
  a  <- anemie[[k]]; yy <- row_ys_a[k]
  bg <- if (k %% 2 == 0) "#f0f0f0" else "white"
  rect(0.02, yy - 0.045, 0.49, yy + 0.045, col = bg, border = NA)
  # badge grade
  rect(0.02, yy - 0.038, 0.09, yy + 0.038, col = gc[a$g], border = NA)
  text(0.055,   yy, a$g, col = gc_txt[a$g], font = 2, cex = 0.95)
  text(sub_x[2], yy, a$hgb, cex = 0.88)
  text(sub_x[3], yy, a$def, cex = 0.85)
  text(sub_x[4], yy, a$nhp, cex = 0.72, col = "grey40")
}
rect(0.02, 0.245, 0.49, 0.90, col = NA, border = "#2166ac", lwd = 1.5)

# Note LLN anémie
text(0.02, 0.21, "LLN humain : Hgb ~12.0 g/dL (femme) / 13.5 g/dL (homme)  |  NHP cynomolgus : 12–16 g/dL",
     cex = 0.72, col = "grey30", adj = 0)

# ──────────────────────────────────────────────────────
# BLOC NEUTROPÉNIE  (ANC)
# ──────────────────────────────────────────────────────
rect(0.51, 0.71, 0.98, 0.90, col = "#b2182b", border = NA)
text(0.745, 0.805,
     "NEUTROPÉNIE  (ANC — Absolute Neutrophil Count)",
     col = "white", cex = 1.05, font = 2)

sub_x2 <- c(0.545, 0.660, 0.840, 0.945)
rect(0.51, 0.63, 0.98, 0.71, col = "#fde0d9", border = NA)
text(sub_x2[1], 0.67, "Grade", font = 2, cex = 0.9)
text(sub_x2[2], 0.67, "ANC (×10⁹/L)", font = 2, cex = 0.9)
text(sub_x2[3], 0.67, "Définition clinique", font = 2, cex = 0.9)
text(sub_x2[4], 0.67, "NHP ref.", font = 2, cex = 0.9)

neutro <- list(
  list(g="G1", anc="LLN – 1.5",  def="Asymptomatique",                 nhp="LLN singe ≈ 1.0"),
  list(g="G2", anc="1.0 – <1.5", def="Risque infectieux modéré",        nhp=""),
  list(g="G3", anc="0.5 – <1.0", def="Prophylaxie G-CSF recommandée",  nhp=""),
  list(g="G4", anc="<0.5",       def="Neutropénie fébrile — urgence",  nhp="")
)
row_ys_n <- c(0.56, 0.47, 0.38, 0.29)

for (k in seq_along(neutro)) {
  n  <- neutro[[k]]; yy <- row_ys_n[k]
  bg <- if (k %% 2 == 0) "#f0f0f0" else "white"
  rect(0.51, yy - 0.045, 0.98, yy + 0.045, col = bg, border = NA)
  rect(0.51, yy - 0.038, 0.60, yy + 0.038, col = gc[n$g], border = NA)
  text(sub_x2[1], yy, n$g, col = gc_txt[n$g], font = 2, cex = 0.95)
  text(sub_x2[2], yy, n$anc, cex = 0.88)
  text(sub_x2[3], yy, n$def, cex = 0.85)
  text(sub_x2[4], yy, n$nhp, cex = 0.72, col = "grey40")
}
rect(0.51, 0.245, 0.98, 0.90, col = NA, border = "#b2182b", lwd = 1.5)

text(0.51, 0.21,
     "LLN humain : ANC ~1.8 ×10⁹/L  |  NHP cynomolgus : ANC ~1.0–8.0 ×10⁹/L (moy. ~3.0)",
     cex = 0.72, col = "grey30", adj = 0)

# ──────────────────────────────────────────────────────
# Pied de page
# ──────────────────────────────────────────────────────
# Bande légende grades
lx <- seq(0.02, 0.30, length.out = 4)
for (k in 1:4) {
  rect(lx[k], 0.04, lx[k] + 0.06, 0.10, col = gc[k], border = "grey50")
  text(lx[k] + 0.03, 0.07, paste0("Grade ", k),
       col = gc_txt[k], font = 2, cex = 0.78)
}
text(0.38, 0.07, "Source : NCI CTCAE v5.0 (2017)", cex = 0.75, col = "grey40", adj = 0)
text(0.38, 0.04, "T-DXd (DS-8201a) — effets hématologiques attendus : anémie & neutropénie (DXd topoisomérase I)",
     cex = 0.72, col = "grey30", adj = 0)

dev.off()
cat("  -> results_TDXD/NHP_CTCAE_Grades_Anemie_Neut.pdf\n")

cat("  Fichiers générés dans results_TDXD/ :\n")
cat("    -> NHP_PK_validation_Table7.pdf\n")
cat("    -> NHP_PK_3doses_comparison.pdf\n")
cat("    -> NHP_NCA_Table_FDA_vs_Sim.pdf\n")
cat("    -> NHP_CTCAE_Grades_Anemie_Neut.pdf\n")
