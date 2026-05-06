############################################################
# run_fgfr2_nhp.R
# Simulation PK/PD — FGFR2 — SINGE CYNOMOLGUS (NHP)
#
# Modèle (25 états) :
#   PK  : Inhibiteur FGFR2 — 2-cpt + TMDD | Damage ADN
#   PD  : Modèle Fornari 2019
#         MPP → CMP → Neut / Mono
#         MPP → MEP → Ret / RBC / Plt
#
# !! Paramètres PK à calibrer sur données FGFR2 NHP !!
# (structure issue du modèle T-DXd — PK et slopes adaptés FGFR2)
############################################################
library(deSolve)

source("parameters_fgfr2_nhp.R")
source("parameters_nhp.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("pkpd_fgfr2_nhp.R")

if (!dir.exists("results")) dir.create("results")

# ════════════════════════════════════════════════════════
# Paramètres consolidés FGFR2 + Fornari NHP
# ════════════════════════════════════════════════════════
build_pars <- function() {
  p <- init_pars
  for (nm in names(fgfr2_nhp)) p[[nm]] <- fgfr2_nhp[[nm]]
  p$k_dam_DXd  <- 0.017
  p$rate_fun   <- NULL
  p
}

# ── État initial (23 états) ──────────────────────────────
pd_state <- init_state[setdiff(names(init_state), c("C1", "C2", "Damage"))]
nhp_state0 <- c(
  C_ADC1 = 0, C_ADC2 = 0, Damage = 0,
  pd_state
)

# ════════════════════════════════════════════════════════
# Calibration TMDD — grille + Nelder-Mead sur FDA Table 7
# ════════════════════════════════════════════════════════
cat("\nCalibration TMDD (grille + Nelder-Mead) ...\n")
doses_cal <- c(3, 10, 30)

sim_one_tmdd <- function(dose_mgkg, CL_lin, Vmax_MM, Km_MM) {
  p <- build_pars()
  p$CL_lin  <- CL_lin
  p$Vmax_MM <- Vmax_MM
  p$Km_MM   <- Km_MM
  p$rate_fun <- make_nhp_infusion(dose_mgkg = dose_mgkg, BW_kg = 4.0,
                                   Tinfu_h = 0.5, n_cycles = 1)
  times <- c(seq(0, 2, by = 0.1), seq(3, 504, by = 1))
  sol <- tryCatch(
    suppressWarnings(as.data.frame(ode(
      y = nhp_state0, times = times, func = pkpd_nhp_ode,
      parms = p, method = "lsoda"))),
    error = function(e) NULL)
  if (is.null(sol) || any(is.nan(sol$C_ADC1)) || min(sol$C_ADC1) < -1e-6)
    return(NULL)
  sol$C_ADC1 <- pmax(sol$C_ADC1, 1e-12)
  sol
}

nca_tmdd <- function(sol) {
  if (is.null(sol)) return(list(C0=NA, AUC=NA, t12=NA, ok=FALSE))
  C0  <- max(sol$C_ADC1[sol$time <= 1])
  idx <- sol$time <= 504
  AUC <- sum(diff(sol$time[idx]) *
             (sol$C_ADC1[idx][-sum(idx)] + sol$C_ADC1[idx][-1]) / 2) / 24
  idt <- sol$time >= 100 & sol$time <= 480 & sol$C_ADC1 > 0
  if (sum(idt) < 5) return(list(C0=C0, AUC=AUC, t12=NA, ok=FALSE))
  lm_f <- tryCatch(lm(log(C_ADC1) ~ time, data = sol[idt, ]), error=function(e) NULL)
  if (is.null(lm_f) || coef(lm_f)[2] >= 0) return(list(C0=C0, AUC=AUC, t12=NA, ok=FALSE))
  list(C0=C0, AUC=AUC, t12=log(2)/(-coef(lm_f)[2])/24, ok=TRUE)
}

wrss_tmdd <- function(CL_lin, Vmax_MM, Km_MM) {
  total <- 0
  w_AUC <- c(2, 2, 1); w_t12 <- c(1, 1, 2)
  for (i in seq_along(doses_cal)) {
    fda <- fda_tk_nhp[[i]]
    nca <- nca_tmdd(sim_one_tmdd(doses_cal[i], CL_lin, Vmax_MM, Km_MM))
    if (!nca$ok || nca$C0<=0 || nca$AUC<=0 || nca$t12<=0) return(1e8)
    total <- total +
      (log(nca$C0/fda$C0_ADC))^2 +
      w_AUC[i]*(log(nca$AUC/fda$AUC21d_ADC))^2 +
      w_t12[i]*(log(nca$t12/fda$t_half_d))^2
  }
  total
}

CL_g <- exp(seq(log(5e-5), log(2e-3), length.out=7))
VM_g <- exp(seq(log(0.05),  log(6.0),  length.out=7))
Km_g <- exp(seq(log(50),    log(2000), length.out=9))
best_val <- 1e8; best_par <- c(1e-3, 0.2, 400)
for (cl in CL_g) for (vm in VM_g) for (km in Km_g) {
  v <- wrss_tmdd(cl, vm, km)
  if (v < best_val) { best_val <- v; best_par <- c(cl, vm, km) }
}
opt <- optim(log(best_par),
             function(th) wrss_tmdd(exp(th[1]), exp(th[2]), exp(th[3])),
             method = "Nelder-Mead",
             control = list(maxit = 5000, reltol = 1e-10))
CL_lin_cal  <- exp(opt$par[1])
Vmax_MM_cal <- exp(opt$par[2])
Km_MM_cal   <- exp(opt$par[3])
fgfr2_nhp$CL_lin  <- CL_lin_cal
fgfr2_nhp$Vmax_MM <- Vmax_MM_cal
fgfr2_nhp$Km_MM   <- Km_MM_cal

cat(sprintf("TMDD calibré : CL=%.3e  Vmax=%.4f  Km=%.1f  RMSE=%.1f%%\n",
            CL_lin_cal, Vmax_MM_cal, Km_MM_cal, 100*sqrt(opt$value/9)))

# ════════════════════════════════════════════════════════
# Simulation helper
# ════════════════════════════════════════════════════════
simulate_nhp <- function(dose_tdxd_mgkg, n_cycles = 3,
                          Tinfu_h = 0.5, BW_kg = 4.0) {
  p <- build_pars()
  p$CL_lin  <- CL_lin_cal
  p$Vmax_MM <- Vmax_MM_cal
  p$Km_MM   <- Km_MM_cal
  p$rate_fun <- make_nhp_infusion(
    dose_mgkg = dose_tdxd_mgkg, BW_kg = BW_kg,
    Tinfu_h = Tinfu_h, interval_h = fgfr2_nhp$interval_h,
    n_cycles = n_cycles)
  times <- seq(0, max(n_cycles * 21 * 24, 120 * 24), by = 1)
  sol <- as.data.frame(ode(
    y = nhp_state0, times = times,
    func = pkpd_nhp_ode, parms = p, method = "lsoda",
    hmax = Tinfu_h / 2))
  sol$time_d <- sol$time / 24
  sol
}

# ════════════════════════════════════════════════════════
# Simulations FGFR2 inhibiteur Q3W × 3 cycles
# ════════════════════════════════════════════════════════
cat("\nSimulations FGFR2 inhibiteur Q3W × 3 ...\n")
doses_fgfr2 <- c(3, 10, 30)
sims_fgfr2  <- lapply(doses_fgfr2, function(d) simulate_nhp(d, n_cycles = 3))

# ════════════════════════════════════════════════════════
# DONNÉES PK OBSERVÉES — à remplir manuellement
# time_h = temps en heures depuis la dose
# conc   = concentration FGFR2 mesurée (µg/mL = mg/L)
# ════════════════════════════════════════════════════════
pk_obs <- list(
  "3"  = data.frame(
    time_h = c(),   # ex: c(0.5, 1, 2, 4, 8, 24, 48, 168)
    conc   = c()    # ex: c(98, 95, 90, ...)
  ),
  "10" = data.frame(
    time_h = c(),
    conc   = c()
  ),
  "30" = data.frame(
    time_h = c(),
    conc   = c()
  )
)

# ════════════════════════════════════════════════════════
# Graphiques PK — validation FDA Table 7 + données observées
# ════════════════════════════════════════════════════════
dose_cols <- c("#2166ac", "#4dac26", "#d6604d")
dose_days <- c(0, 21, 42)

pdf("results/NHP_PK_validation_Table7.pdf", width = 10, height = 5)
par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1.5))
for (i in seq_along(doses_fgfr2)) {
  s <- sims_fgfr2[[i]]; col <- dose_cols[i]; d <- doses_fgfr2[i]
  fda <- fda_tk_nhp[[i]]

  plot(s$time_d, s$C_ADC1, type="l", lwd=2.5, col=col, log="y",
       xlab="Temps (j)", ylab="FGFR2 inhibiteur [µg/mL]",
       main=sprintf("FGFR2 — %d mg/kg Q3W", d),
       ylim=c(1, max(s$C_ADC1)*2))
  abline(v=dose_days, lty=2, col="grey60")
  points(0.02, fda$C0_ADC, pch=19, cex=1.5)
  obs <- pk_obs[[as.character(d)]]
  if (nrow(obs) > 0)
    points(obs$time_h / 24, obs$conc, pch=17, cex=1.4, col=col)
  has_obs <- nrow(obs) > 0
  legend("topright",
         c("Sim", "C0 FDA", if (has_obs) "Observé"),
         col   = c(col, "black", if (has_obs) col),
         lwd   = c(2.5, NA,     if (has_obs) NA),
         pch   = c(NA,  19,     if (has_obs) 17),
         bty   = "n", cex = 0.85)
}
dev.off()
cat("  -> results/NHP_PK_validation_Table7.pdf\n")

# ════════════════════════════════════════════════════════
# Graphiques PD — Fornari
# ════════════════════════════════════════════════════════
pdf("results/NHP_PD_Fornari_FGFR2.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4.5, 3, 1.5))

cell_info <- list(
  list(var="Neut",   label="Neutrophiles (10⁹/L)", base="Neut0"),
  list(var="Plt",    label="Plaquettes (10⁹/L)",   base="Plt0"),
  list(var="Ret",    label="Réticulocytes (10⁹/L)", base="Ret0"),
  list(var="RBC",    label="GR (10⁹/L)",            base="RBC0"),
  list(var="Mono",   label="Monocytes (10⁹/L)",     base="Mono0"),
  list(var="Damage", label="Damage (u.a.)",              base=NULL)
)

for (ci in cell_info) {
  y_vals <- lapply(sims_fgfr2, function(s) s[[ci$var]])
  y_max  <- max(unlist(y_vals), na.rm=TRUE) * 1.1
  y_min  <- min(unlist(y_vals), na.rm=TRUE) * 0.9
  plot(sims_fgfr2[[1]]$time_d, y_vals[[1]], type="n",
       ylim=c(y_min, y_max), xlab="Temps (j)", ylab=ci$label,
       main=ci$label)
  for (i in seq_along(doses_fgfr2))
    lines(sims_fgfr2[[i]]$time_d, y_vals[[i]], col=dose_cols[i], lwd=2)
  if (!is.null(ci$base))
    abline(h=init_pars[[ci$base]], lty=2, col="grey40")
  abline(v=dose_days, lty=3, col="grey70")
  legend("bottomright", paste0(doses_fgfr2, " mg/kg"),
         col=dose_cols, lwd=2, bty="n", cex=0.85)
}
dev.off()
cat("  -> results/NHP_PD_Fornari_FGFR2.pdf\n")

# ════════════════════════════════════════════════════════
# Validation NCA — table console
# ════════════════════════════════════════════════════════
cat("\n═══ VALIDATION NCA FDA Table 7 ═══\n")
cat(sprintf("  %-10s │ C0_obs  C0_sim  ratio │ AUC_obs AUC_sim ratio │ T½_obs T½_sim\n", "Dose"))
cat(sprintf("  %s\n", paste(rep("─",75), collapse="")))
for (i in seq_along(doses_fgfr2)) {
  s <- sims_fgfr2[[i]]; fda <- fda_tk_nhp[[i]]
  C0  <- max(s$C_ADC1[s$time <= 24])
  idx <- s$time <= 504
  AUC <- sum(diff(s$time[idx])*(s$C_ADC1[idx][-sum(idx)]+s$C_ADC1[idx][-1])/2)/24
  idt <- s$time>=100 & s$time<=480 & s$C_ADC1>0
  t12 <- if(sum(idt)>5) log(2)/abs(coef(lm(log(C_ADC1)~time,data=s[idt,]))[2])/24 else NA
  cat(sprintf("  %-10s │ %6.1f  %6.1f  %5.3f │ %7.0f %7.0f %5.3f │ %6.2f %6.2f\n",
              paste0(doses_fgfr2[i]," mg/kg"),
              fda$C0_ADC, C0, C0/fda$C0_ADC,
              fda$AUC21d_ADC, AUC, AUC/fda$AUC21d_ADC,
              fda$t_half_d, t12))
}
cat(paste(rep("═",57), collapse=""), "\n")

# ════════════════════════════════════════════════════════
# COMPTES CELLULAIRES — Jours 2, 8, 15, 22
# ════════════════════════════════════════════════════════
target_days <- c(2, 8, 15, 22)
cell_vars   <- c("Neut", "Mono", "Ret", "RBC", "Plt", "MPP", "CMP", "MEP")

closest_row <- function(sol, day) sol[which.min(abs(sol$time_d - day)), ]

records <- list()
for (i in seq_along(doses_fgfr2)) {
  s <- sims_fgfr2[[i]]
  for (d in target_days) {
    row <- closest_row(s, d)
    rec <- data.frame(Dose_mgkg = doses_fgfr2[i], Jour = d)
    for (v in cell_vars) rec[[v]] <- round(row[[v]], 3)
    records <- c(records, list(rec))
  }
}
cell_table <- do.call(rbind, records)

cat("\n════════════════════════════════════════════════════════════\n")
cat("  COMPTES CELLULAIRES AUX JOURS 2 / 8 / 15 / 22\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("  %-8s │ %-4s │ %8s %8s %8s %8s %8s │ %7s %7s %7s\n",
            "Dose","Jour","Neut","Mono","Ret","RBC","Plt","MPP","CMP","MEP"))
cat(sprintf("  %s\n", paste(rep("─", 88), collapse="")))

for (i in seq_len(nrow(cell_table))) {
  r <- cell_table[i, ]
  cat(sprintf("  %-8s │ J%-3d │ %8.3f %8.3f %8.1f %8.0f %8.1f │ %7.3f %7.3f %7.3f\n",
              paste0(r$Dose_mgkg," mg/kg"), r$Jour,
              r$Neut, r$Mono, r$Ret, r$RBC, r$Plt,
              r$MPP,  r$CMP,  r$MEP))
  if (i %% length(target_days) == 0)
    cat(sprintf("  %s\n", paste(rep("─", 88), collapse="")))
}
cat(sprintf("  %-8s │      │ %8.3f %8.3f %8.1f %8.0f %8.1f │ %7.3f %7.3f %7.3f\n",
            "Baseline", init_pars$Neut0, init_pars$Mono0, init_pars$Ret0,
            init_pars$RBC0, init_pars$Plt0,
            init_pars$MPP0, init_pars$CMP0, init_pars$MEP0))

# ════════════════════════════════════════════════════════
# NADIRS — valeur minimale + jour sur toute la simulation
# ════════════════════════════════════════════════════════
cat("\n════════════════════════════════════════════════════════════\n")
cat("  NADIRS (valeur min + jour) sur 3 cycles Q3W\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("  %-8s │ %14s %14s %14s %14s %14s\n",
            "Dose", "Neut", "Mono", "Ret", "RBC", "Plt"))
cat(sprintf("  %-8s │ %14s %14s %14s %14s %14s\n",
            "", "val  (jour)", "val  (jour)", "val  (jour)", "val  (jour)", "val  (jour)"))
cat(sprintf("  %s\n", paste(rep("─", 85), collapse="")))

nadir_records <- list()
for (i in seq_along(doses_fgfr2)) {
  s   <- sims_fgfr2[[i]]
  rec <- data.frame(Dose_mgkg = doses_fgfr2[i])
  row_parts <- c()
  for (v in c("Neut", "Mono", "Ret", "RBC", "Plt")) {
    idx      <- which.min(s[[v]])
    val      <- s[[v]][idx]
    day      <- s$time_d[idx]
    pct      <- round((val / init_pars[[paste0(v,"0")]] - 1) * 100, 1)
    rec[[paste0(v,"_nadir")]]    <- round(val, 3)
    rec[[paste0(v,"_nadir_day")]] <- round(day, 1)
    rec[[paste0(v,"_pct")]]      <- pct
    row_parts <- c(row_parts, sprintf("%6.2f (J%4.1f)", val, day))
  }
  nadir_records <- c(nadir_records, list(rec))
  cat(sprintf("  %-8s │ %14s %14s %14s %14s %14s\n",
              paste0(doses_fgfr2[i]," mg/kg"),
              row_parts[1], row_parts[2], row_parts[3],
              row_parts[4], row_parts[5]))
}
cat(sprintf("  %s\n", paste(rep("─", 85), collapse="")))
cat(sprintf("  %-8s │ %14s %14s %14s %14s %14s\n",
            "Baseline",
            sprintf("%6.1f  (J0)  ", init_pars$Neut0),
            sprintf("%6.2f  (J0)  ", init_pars$Mono0),
            sprintf("%6.1f  (J0)  ", init_pars$Ret0),
            sprintf("%6.0f  (J0)  ", init_pars$RBC0),
            sprintf("%6.1f  (J0)  ", init_pars$Plt0)))
cat("════════════════════════════════════════════════════════════\n")

nadir_table <- do.call(rbind, nadir_records)
write.csv(cell_table,   "results/cell_counts_J2_J8_J15_J22.csv", row.names = FALSE)
write.csv(nadir_table,  "results/nadirs.csv",                    row.names = FALSE)
cat("\n  -> results/cell_counts_J2_J8_J15_J22.csv\n")
cat("  -> results/nadirs.csv\n")
