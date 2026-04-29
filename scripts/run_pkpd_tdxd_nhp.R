############################################################
# run_pkpd_tdxd_nhp.R
# Modèle PKPD combiné NHP — T-DXd (DS-8201a) + Carboplatin
#
# ODE : 27 variables d'état
#   PK  : ADC 2-cpt + TMDD | DXd 1-cpt | Carboplatin 2-cpt
#   PD  : Modèle Fornari 2019 (moelle osseuse)
#         MPP → CMP → Neut/Mono
#         MPP → MEP → Ret/RBC/Plt
#
# Damage combiné :
#   carboplatin  → k_dam × C_libre_µM       (linéaire, Fornari)
#   DXd          → k_dam_DXd × Emax(DXd_ic) (saturé)
#
# Sources :
#   pkpd_model_FORNARI.R  — feedbacks et équations PD
#   parameters_tdxd_nhp.R — PK T-DXd + PK carboplatin NHP
#   parameters_nhp.R      — baselines hématologiques NHP
#   parameters_FORNARI_CORRECT.R — constantes dérivées (Eq. S4)
############################################################
library(deSolve)

source("parameters_tdxd_nhp.R")         # tdxd_nhp, carbo_nhp, fda_tk_nhp,
                                         # make_nhp_infusion, tdxd_nhp_state0
source("parameters_nhp.R")              # init_pars (baselines NHP)
source("parameters_FORNARI_CORRECT.R")  # constantes cinétiques dérivées

if (!dir.exists("results_TDXD")) dir.create("results_TDXD")

# ════════════════════════════════════════════════════════
# ODE COMBINÉE — 27 variables d'état
# ════════════════════════════════════════════════════════
# État : C_ADC1, C_ADC2, C_DXd, C_DXd_ic,
#        C1(carbo), C2(carbo), Damage,
#        MPP, CMP, MEP,
#        T1_Neut..Neut, T1_Mono..Mono,
#        T1_Ret..Ret, RBC,
#        T1_Plt..Plt

pkpd_nhp_combined <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    rate_in       <- if (!is.null(pars$rate_fun))       pars$rate_fun(time)       else 0
    rate_in_carbo <- if (!is.null(pars$rate_fun_carbo)) pars$rate_fun_carbo(time) else 0

    # ── T-DXd PK — ADC 2-cpt + TMDD Michaelis-Menten ────
    cycle_num <- max(1L, floor(time / interval_h) + 1L)
    krel_t    <- k_rel_c1 * cycle_num^krel_power *
                 ifelse(cycle_num > 1L, krel_factor, 1.0)

    dC_ADC1 <- rate_in / V1_ADC +
               (Q_ADC / V2_ADC) * C_ADC2 -
               (CL_lin / V1_ADC + Q_ADC / V1_ADC + k_int) * C_ADC1 -
               Vmax_MM * C_ADC1 / (Km_MM + C_ADC1)
    dC_ADC2 <- (Q_ADC / V1_ADC) * C_ADC1 - (Q_ADC / V2_ADC) * C_ADC2

    # ── DXd plasmatique et intracellulaire ────────────────
    release   <- mass_frac_DXd * krel_t * C_ADC1 * V1_ADC / V_DXd
    dC_DXd    <- release - (CL_DXd / V_DXd) * C_DXd -
                 k_inD * C_DXd + k_effD * C_DXd_ic * (V_ic / V_DXd)
    dC_DXd_ic <- k_inD * C_DXd * (V_DXd / V_ic) - k_effD * C_DXd_ic

    # ── Carboplatin — 2 compartiments (fit rxode2 NHP) ───
    # C1, C2 : noms Fornari conservés pour compatibilité
    Q_c <- if (is.na(Q_carbo) || is.null(Q_carbo)) 0 else Q_carbo
    dC1 <- rate_in_carbo / V1_carbo -
           (CL_carbo / V1_carbo + Q_c / V1_carbo) * C1 +
           (Q_c / V2_carbo) * C2
    dC2 <- (Q_c / V1_carbo) * C1 - (Q_c / V2_carbo) * C2

    # ── Damage combiné (Eq. 3 Fornari + DXd Emax) ────────
    fu_t         <- fu_inf + (fu0 - fu_inf) * exp(-k_bind * time)
    C_free_uM    <- fu_t * C1 * mgL_to_uM          # carboplatin libre (µM)
    E_DXd        <- (C_DXd_ic * mgL_to_uM_DXd) /
                    (IC50_DXd_uM + C_DXd_ic * mgL_to_uM_DXd)
    dDamage      <- k_dam * C_free_uM +
                    k_dam_DXd * E_DXd -
                    k_rep * Damage

    # ── Feedbacks Fornari (Eq. 11, 12, 13) ───────────────
    eps <- 1e-9
    r_stem   <- 0.5*(CMP0/max(CMP,eps)) + 0.5*(MEP0/max(MEP,eps))
    f_stem   <- pmin(pmax(r_stem,   0.2), 50)^gamma_stem

    r_matCMP  <- 0.5*(Neut0/max(Neut,eps)) + 0.5*(Mono0/max(Mono,eps))
    f_mat_CMP <- pmin(pmax(r_matCMP, 0.2), 50)^gamma_mat_CMP

    r_matMEP  <- RBC0 / max(RBC, eps)
    f_mat_MEP <- pmin(pmax(r_matMEP, 0.2), 50)^gamma_mat_MEP

    r_pRet     <- 0.5*(Ret0/max(Ret,eps)) + 0.5*(RBC0/max(RBC,eps))
    f_prol_Ret <- pmin(pmax(r_pRet, 0.2), 10)^gamma_prolTrans

    r_pPlt     <- Plt0 / max(Plt, eps)
    f_prol_Plt <- pmin(pmax(r_pPlt, 0.2), 10)^gamma_prolTrans

    # ── Progéniteurs (Eq. 1, 4) ──────────────────────────
    dMPP <- k_stem * f_stem +
            k_prol_MPP * (1 - Slope_MPP * Damage) * MPP -
            k_tr_CMP * f_mat_CMP * MPP -
            k_tr_MEP * f_mat_MEP * MPP

    dCMP <- k_prol_CMP * (1 - Slope_CMP * Damage) * CMP +
            k_tr_CMP * f_mat_CMP * MPP -
            (k_tr_Neut + k_tr_Mono) * CMP

    dMEP <- k_prol_MEP * (1 - Slope_MEP * Damage) * MEP +
            k_tr_MEP * f_mat_MEP * MPP -
            (k_tr_Ret + k_tr_Plt) * MEP

    # ── Neutrophiles ──────────────────────────────────────
    a_Neut   <- 3 / MTT_Neut
    dT1_Neut <- k_tr_Neut*CMP  - a_Neut*T1_Neut
    dT2_Neut <- a_Neut*T1_Neut - a_Neut*T2_Neut
    dT3_Neut <- a_Neut*T2_Neut - a_Neut*T3_Neut
    dNeut    <- a_Neut*T3_Neut - k_circ_Neut*Neut

    # ── Monocytes ─────────────────────────────────────────
    a_Mono   <- 3 / MTT_Mono
    dT1_Mono <- k_tr_Mono*CMP  - a_Mono*T1_Mono
    dT2_Mono <- a_Mono*T1_Mono - a_Mono*T2_Mono
    dT3_Mono <- a_Mono*T2_Mono - a_Mono*T3_Mono
    dMono    <- a_Mono*T3_Mono - k_circ_Mono*Mono

    # ── Réticulocytes / RBC ───────────────────────────────
    drug_ret <- delta_Ret * Slope_MEP * k_prol_Ret * Damage
    dT1_Ret  <- k_prol_Ret*f_prol_Ret*T1_Ret - drug_ret*T1_Ret +
                k_tr_Ret*MEP - a_Ret*T1_Ret
    dT2_Ret  <- k_prol_Ret*f_prol_Ret*T2_Ret - drug_ret*T2_Ret +
                a_Ret*T1_Ret - a_Ret*T2_Ret
    dT3_Ret  <- a_Ret*T2_Ret  - a_Ret*T3_Ret
    dRet     <- a_Ret*T3_Ret  - k_circ_Ret*Ret
    dRBC     <- k_circ_Ret*Ret - k_circ_RBC*RBC

    # ── Plaquettes ────────────────────────────────────────
    drug_plt <- delta_Plt * Slope_MEP * k_prol_Plt * Damage
    dT1_Plt  <- k_prol_Plt*f_prol_Plt*T1_Plt - drug_plt*T1_Plt +
                k_tr_Plt*MEP - a_Plt*T1_Plt
    dT2_Plt  <- k_prol_Plt*f_prol_Plt*T2_Plt - drug_plt*T2_Plt +
                a_Plt*T1_Plt - a_Plt*T2_Plt
    dT3_Plt  <- a_Plt*T2_Plt - a_Plt*T3_Plt
    dPlt     <- a_Plt*T3_Plt - k_circ_Plt*Plt

    list(c(dC_ADC1, dC_ADC2, dC_DXd, dC_DXd_ic,
           dC1, dC2, dDamage,
           dMPP, dCMP, dMEP,
           dT1_Neut, dT2_Neut, dT3_Neut, dNeut,
           dT1_Mono, dT2_Mono, dT3_Mono, dMono,
           dT1_Ret,  dT2_Ret,  dT3_Ret,  dRet, dRBC,
           dT1_Plt,  dT2_Plt,  dT3_Plt,  dPlt))
  })
}

# ── Paramètres consolidés (T-DXd + carboplatin + Fornari NHP) ─
build_pars <- function(base_pars = init_pars) {
  p <- base_pars
  # T-DXd
  for (nm in names(tdxd_nhp)) p[[nm]] <- tdxd_nhp[[nm]]
  # Carboplatin NHP (noms Fornari : CL, V1, Q, V2, fu0, fu_inf, k_bind)
  p$CL_carbo   <- carbo_nhp$CL
  p$V1_carbo   <- carbo_nhp$V1   # stocké séparément pour l'ODE
  p$V2_carbo   <- carbo_nhp$V2
  p$Q_carbo    <- carbo_nhp$Q
  p$fu0        <- carbo_nhp$fu0
  p$fu_inf     <- carbo_nhp$fu_inf
  p$k_bind     <- carbo_nhp$k_bind
  p$k_dam_DXd  <- 0.017          # effet DXd sur l'ADN (identique k_dam carbo)
  p$rate_fun       <- NULL
  p$rate_fun_carbo <- NULL
  p
}
base_pars <- build_pars()

# ── État initial complet (PK + Fornari PD) ────────────────────
nhp_state0 <- c(
  C_ADC1 = 0, C_ADC2 = 0, C_DXd = 0, C_DXd_ic = 0,
  C1 = 0, C2 = 0,
  Damage = 0,
  init_state[setdiff(names(init_state), c("C1", "C2", "Damage"))]
)

# ════════════════════════════════════════════════════════
# Calibration TMDD (PK seule — regarde uniquement C_ADC1)
# ════════════════════════════════════════════════════════
cat("\nCalibration TMDD (grille + Nelder-Mead) ...\n")
doses_cal <- c(3, 10, 30)

sim_one_tmdd <- function(dose_mgkg, CL_lin, Vmax_MM, Km_MM) {
  p <- build_pars()
  p$CL_lin  <- CL_lin
  p$Vmax_MM <- Vmax_MM
  p$Km_MM   <- Km_MM
  p$rate_fun <- make_nhp_infusion(dose_mgkg = dose_mgkg, BW_kg = 4.0,
                                   Tinfu_h = 0.5, interval_h = NULL, n_cycles = 1)
  times <- c(seq(0, 2, by = 0.1), seq(3, 504, by = 1))
  sol <- tryCatch(
    suppressWarnings(as.data.frame(ode(
      y = nhp_state0, times = times, func = pkpd_nhp_combined,
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
tdxd_nhp$CL_lin  <- CL_lin_cal
tdxd_nhp$Vmax_MM <- Vmax_MM_cal
tdxd_nhp$Km_MM   <- Km_MM_cal

cat(sprintf("TMDD calibré : CL=%.3e  Vmax=%.4f  Km=%.1f  RMSE=%.1f%%\n",
            CL_lin_cal, Vmax_MM_cal, Km_MM_cal,
            100*sqrt(opt$value/9)))

# ════════════════════════════════════════════════════════
# Simulation helper PKPD
# ════════════════════════════════════════════════════════
simulate_nhp <- function(dose_tdxd_mgkg, dose_carbo_mgkg = 0,
                          n_cycles = 3, Tinfu_h_tdxd = 0.5,
                          Tinfu_h_carbo = 1, BW_kg = 4.0) {
  p <- build_pars()
  p$CL_lin  <- CL_lin_cal
  p$Vmax_MM <- Vmax_MM_cal
  p$Km_MM   <- Km_MM_cal

  if (dose_tdxd_mgkg > 0)
    p$rate_fun <- make_nhp_infusion(
      dose_mgkg = dose_tdxd_mgkg, BW_kg = BW_kg,
      Tinfu_h = Tinfu_h_tdxd, interval_h = tdxd_nhp$interval_h,
      n_cycles = n_cycles)

  if (dose_carbo_mgkg > 0)
    p$rate_fun_carbo <- make_nhp_carbo_infusion(
      dose_mgkg = dose_carbo_mgkg, BW_kg = BW_kg,
      Tinfu_h = Tinfu_h_carbo, interval_h = tdxd_nhp$interval_h,
      n_cycles = n_cycles)

  times <- seq(0, n_cycles * 21 * 24, by = 1)

  sol <- as.data.frame(ode(
    y = nhp_state0, times = times,
    func = pkpd_nhp_combined, parms = p, method = "lsoda"))
  sol$time_d      <- sol$time / 24
  sol$C_DXd_ngmL  <- sol$C_DXd * 1e3
  sol
}

# ════════════════════════════════════════════════════════
# Simulations : T-DXd seul (validation FDA Table 7)
# ════════════════════════════════════════════════════════
cat("\nSimulations T-DXd seul Q3W × 3 ...\n")
doses_tdxd <- c(3, 10, 30)
sims_tdxd  <- lapply(doses_tdxd, function(d)
  simulate_nhp(dose_tdxd_mgkg = d, n_cycles = 3))

# ════════════════════════════════════════════════════════
# Graphiques PK — validation FDA Table 7
# ════════════════════════════════════════════════════════
dose_cols <- c("#2166ac", "#4dac26", "#d6604d")
dose_days <- c(0, 21, 42)

pdf("results_TDXD/NHP_PK_validation_Table7.pdf", width = 14, height = 9)
par(mfrow = c(2, 3), mar = c(4, 4.5, 3, 1.5))
for (i in seq_along(doses_tdxd)) {
  s <- sims_tdxd[[i]]; col <- dose_cols[i]; d <- doses_tdxd[i]
  fda <- fda_tk_nhp[[i]]

  plot(s$time_d, s$C_ADC1, type="l", lwd=2.5, col=col, log="y",
       xlab="Temps (j)", ylab="ADC [µg/mL]",
       main=sprintf("ADC — %d mg/kg Q3W", d),
       ylim=c(1, max(s$C_ADC1)*2))
  abline(v=dose_days, lty=2, col="grey60")
  points(0.02, fda$C0_ADC, pch=19, cex=1.5)
  legend("topright", c("Sim","C0 FDA"), col=c(col,"black"),
         lwd=c(2.5,NA), pch=c(NA,19), bty="n", cex=0.85)

  plot(s$time_d, s$C_DXd_ngmL, type="l", lwd=2.5, col=col,
       xlab="Temps (j)", ylab="DXd [ng/mL]",
       main=sprintf("DXd — %d mg/kg Q3W", d))
  abline(v=dose_days, lty=2, col="grey60")
  points(0.02, fda$C0_DXd_ng, pch=17, cex=1.5)
}
dev.off()
cat("  -> results_TDXD/NHP_PK_validation_Table7.pdf\n")

# ════════════════════════════════════════════════════════
# Graphiques PD — Fornari (Neut, Plt, Ret) T-DXd seul
# ════════════════════════════════════════════════════════
pdf("results_TDXD/NHP_PD_Fornari_TDXd.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4.5, 3, 1.5))

cell_info <- list(
  list(var="Neut", label="Neutrophiles (10⁹/L)", base="Neut0"),
  list(var="Plt",  label="Plaquettes (10⁹/L)",   base="Plt0"),
  list(var="Ret",  label="Réticulocytes (10⁹/L)", base="Ret0"),
  list(var="RBC",  label="GR (10⁹/L)",            base="RBC0"),
  list(var="Mono", label="Monocytes (10⁹/L)",     base="Mono0"),
  list(var="Damage", label="Damage (u.a.)",        base=NULL)
)

for (ci in cell_info) {
  y_vals <- lapply(sims_tdxd, function(s) s[[ci$var]])
  y_max  <- max(unlist(y_vals), na.rm=TRUE) * 1.1
  y_min  <- min(unlist(y_vals), na.rm=TRUE) * 0.9
  plot(sims_tdxd[[1]]$time_d, y_vals[[1]], type="n",
       ylim=c(y_min, y_max), xlab="Temps (j)", ylab=ci$label,
       main=ci$label)
  for (i in seq_along(doses_tdxd))
    lines(sims_tdxd[[i]]$time_d, y_vals[[i]], col=dose_cols[i], lwd=2)
  if (!is.null(ci$base))
    abline(h=init_pars[[ci$base]], lty=2, col="grey40")
  abline(v=dose_days, lty=3, col="grey70")
  legend("bottomright", paste0(doses_tdxd, " mg/kg"),
         col=dose_cols, lwd=2, bty="n", cex=0.85)
}
dev.off()
cat("  -> results_TDXD/NHP_PD_Fornari_TDXd.pdf\n")

# ════════════════════════════════════════════════════════
# Validation NCA — table console
# ════════════════════════════════════════════════════════
cat("\n═══ VALIDATION NCA FDA Table 7 ═══\n")
cat(sprintf("  %-10s │ C0_obs  C0_sim  ratio │ AUC_obs AUC_sim ratio │ T½_obs T½_sim\n",
            "Dose"))
cat(sprintf("  %s\n", paste(rep("─",75),collapse="")))
for (i in seq_along(doses_tdxd)) {
  s <- sims_tdxd[[i]]; fda <- fda_tk_nhp[[i]]
  C0  <- max(s$C_ADC1[s$time <= 24])
  idx <- s$time <= 504
  AUC <- sum(diff(s$time[idx])*(s$C_ADC1[idx][-sum(idx)]+s$C_ADC1[idx][-1])/2)/24
  idt <- s$time>=100 & s$time<=480 & s$C_ADC1>0
  t12 <- if(sum(idt)>5) log(2)/abs(coef(lm(log(C_ADC1)~time,data=s[idt,]))[2])/24 else NA
  cat(sprintf("  %-10s │ %6.1f  %6.1f  %5.3f │ %7.0f %7.0f %5.3f │ %6.2f %6.2f\n",
              paste0(doses_tdxd[i]," mg/kg"),
              fda$C0_ADC,C0,C0/fda$C0_ADC,
              fda$AUC21d_ADC,AUC,AUC/fda$AUC21d_ADC,
              fda$t_half_d,t12))
}
cat("═══════════════════════════════════════════════════════\n")
cat("\nFichiers générés dans results_TDXD/\n")
