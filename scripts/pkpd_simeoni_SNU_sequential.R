# =============================================================================
# Modele PK/PD TGI — Simeoni (2004) — Fc-silent FGFR2-huBPA-LP1
# Cohorte : SNU (xenogreffe souris)
#
# PK  : 2 compartiments IV bolus lineaire (parametres fixes)
#   dA1/dt = -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
#   dA2/dt =  (Q/V1)*A1 - (Q/V2)*A2
#   C1 = A1/V1
#
# PD  : Modele de Simeoni (compartiments de transit)
#   dx1/dt = [g(w)/w - k2*C1] * x1
#   dx2/dt =  k2*C1*x1 - k1*x2
#   dx3/dt =  k1*x2 - k1*x3
#   dx4/dt =  k1*x3 - k1*x4
#   g(w) = L0*w / (1 + (L0*w/L1)^psi)^(1/psi)   psi = 20
#
#   TSC (ug/L) : concentration seuil de cytotoxicite — estimee directement
#   k2 = L0/TSC (derive)
#   MTT = 4/k1  (jours)
#
# Ajustement SEQUENTIEL :
#   Etape A — Croissance  : L0 et L1 sur le groupe vehicule uniquement
#   Etape B — Efficacite  : k1 et TSC sur les 3 groupes traites (L0/L1 fixes)
#
#   Bornes TSC : [Cmax_1p2/10, Cmax_1p2]
#     Justification biologique :
#       1.2 mg/kg → effet observe (~350 vs ~1500 mm3 a j49) → TSC < Cmax_1p2
#       Borne basse = Cmax_1p2/10 pour autoriser solutions a TSC tres bas
#
# Schema posologique : dose unique IV bolus au jour 0
# Groupes : Vehicule | 1.2 mg/kg | 5.6 mg/kg | 11.8 mg/kg
# =============================================================================

library(deSolve)
library(DEoptim)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMETRES PK FIXES
# =============================================================================
# NOTE TMDD : Km et Vmax fixes (parametres PK, non identifiables depuis donnees PD).
#   Km = 1 ug/L  : C1 >> Km => regime zeroth-order (decline lineaire de C1).
#   Vmax = VMAX_FIXED ug/kg/j : fixe pour garantir que le drug reste actif a
#     5.6 mg/kg pendant toute la fenetre d'observation (50j).
#     Condition : (Dose_5p6 - TSC*V1) / Vmax > 50j
#     => Vmax < (5600 - TSC*V1) / 50 ~ 100 ug/kg/j pour TSC~3000.
#   Seuls k1 et TSC sont estimes en Etape B.

pk_data <- read.csv("pk_resultats_20260720_150418.csv", stringsAsFactors = FALSE)

# x1000 : mL/kg -> L/kg    x24 : /h -> /j
pk_fixed <- c(
  CL = (pk_data$CL[1] * 1000) * 24,
  V1 =  pk_data$V1[1] * 1000,
  V2 =  pk_data$V2[1] * 1000,
  Q  = (pk_data$Q[1]  * 1000) * 24
)

cat("=== Parametres PK fixes (souris, unites journalieres) ===\n")
cat(sprintf("CL = %.5e L/j/kg\n", pk_fixed["CL"]))
cat(sprintf("V1 = %.5e L/kg\n",   pk_fixed["V1"]))
cat(sprintf("V2 = %.5e L/kg\n",   pk_fixed["V2"]))
cat(sprintf("Q  = %.5e L/j/kg\n", pk_fixed["Q"]))

# Cmax (ng/mL = ug/L) issues du CSV
Cmax_1p2  <- pk_data$cmax[pk_data$Animal == "groupe 4"]
Cmax_5p6  <- pk_data$cmax[pk_data$Animal == "grp 5"]
Cmax_11p8 <- pk_data$cmax[pk_data$Animal == "grp 7"]

KM_FIXED   <- 1    # ug/L — fixe (C1 >> Km => regime zeroth-order)
VMAX_FIXED <- 50   # ug/kg/j — fixe : drug actif a 5.6mg/kg pendant >100j

cat(sprintf("\nCmax (ug/L) : 1.2 mg/kg = %.0f | 5.6 mg/kg = %.0f | 11.8 mg/kg = %.0f\n",
            Cmax_1p2, Cmax_5p6, Cmax_11p8))
cat(sprintf("Km   fixe = %.0f ug/L      (C1/Km @ 1.2mg/kg = %.0fx => regime zeroth-order)\n",
            KM_FIXED, Cmax_1p2 / KM_FIXED))
cat(sprintf("Vmax fixe = %.0f ug/kg/j  (drug actif a 5.6mg/kg pendant ~%.0fj)\n",
            VMAX_FIXED, (5600 - 3000 * pk_data$V1[1] * 1000) / VMAX_FIXED))

# =============================================================================
# 2. DONNEES TUMORALES — SNU
#
# Structure Excel — 2 blocs separes par une ligne d'en-tete "Time" :
#   Bloc 1 : moyennes  col 1=Temps, 2=Vehicule, 3=1.2mg/kg, 4=5.6mg/kg, 5=11.8mg/kg
#   Bloc 2 : SEM       (memes colonnes)
# =============================================================================

raw_snu <- suppressMessages(
  read_xlsx("Tumor_volum_SNU.xlsx", col_names = FALSE, .name_repair = "minimal")
)

header_rows <- which(grepl("Time", trimws(as.character(raw_snu[[1]])), ignore.case = TRUE))
stopifnot("2 blocs attendus (moyennes + SEM)" = length(header_rows) == 2)

extract_block <- function(raw, row_start, row_end) {
  blk <- raw[seq(row_start, row_end), ]
  blk[!is.na(suppressWarnings(as.numeric(as.character(blk[[1]])))), ]
}

b1 <- extract_block(raw_snu, header_rows[1] + 1, header_rows[2] - 1)
b2 <- extract_block(raw_snu, header_rows[2] + 1, nrow(raw_snu))

num_col <- function(df, j) suppressWarnings(as.numeric(as.character(df[[j]])))
fix_sem <- function(x) ifelse(is.na(x) | x <= 0, 1, x)

times_d   <- num_col(b1, 1)
tv_ctrl   <- num_col(b1, 2)
tv_d1p2   <- num_col(b1, 3)
tv_d5p6   <- num_col(b1, 4)
tv_d11p8  <- num_col(b1, 5)

sem_ctrl  <- fix_sem(num_col(b2, 2))
sem_d1p2  <- fix_sem(num_col(b2, 3))
sem_d5p6  <- fix_sem(num_col(b2, 4))
sem_d11p8 <- fix_sem(num_col(b2, 5))

get_tv0 <- function(tv) {
  v <- tv[times_d == 0]
  if (length(v) == 0 || is.na(v[1])) tv[!is.na(tv)][1] else v[1]
}

tv0_ctrl  <- get_tv0(tv_ctrl)
tv0_d1p2  <- get_tv0(tv_d1p2)
tv0_d5p6  <- get_tv0(tv_d5p6)
tv0_d11p8 <- get_tv0(tv_d11p8)

ok_ctrl  <- !is.na(tv_ctrl)
ok_d1p2  <- !is.na(tv_d1p2)
ok_d5p6  <- !is.na(tv_d5p6)
ok_d11p8 <- !is.na(tv_d11p8)

cat(sprintf("\nTV0 (j0) : Ctrl=%.0f | 1.2mg=%.0f | 5.6mg=%.0f | 11.8mg=%.0f mm3\n",
            tv0_ctrl, tv0_d1p2, tv0_d5p6, tv0_d11p8))

# =============================================================================
# 3. EQUATIONS DU MODELE DE SIMEONI  <- NE PAS MODIFIER
# =============================================================================

PSI <- 20

simeoni_rhs <- function(t, state, parms) {
  A1 <- max(state["A1"], 0); A2 <- max(state["A2"], 0)
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)

  CL <- parms["CL"]; V1 <- parms["V1"]
  V2 <- parms["V2"]; Q  <- parms["Q"]
  L0 <- parms["L0"]; L1 <- parms["L1"]
  k1 <- parms["k1"]; k2 <- parms["k2"]

  C1  <- A1 / V1
  dA1 <- -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
  dA2 <-  (Q/V1)*A1 - (Q/V2)*A2

  w  <- x1 + x2 + x3 + x4
  gw <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  dx1 <- (growth_rate - k2 * C1) * x1
  dx2 <- k2 * C1 * x1 - k1 * x2
  dx3 <- k1 * x2 - k1 * x3
  dx4 <- k1 * x3 - k1 * x4

  list(c(A1 = dA1, A2 = dA2, x1 = dx1, x2 = dx2, x3 = dx3, x4 = dx4))
}

simeoni_ctrl_rhs <- function(t, state, parms) {
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)
  L0 <- parms["L0"]; L1 <- parms["L1"]

  w  <- x1 + x2 + x3 + x4
  gw <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  list(c(x1 = growth_rate * x1, x2 = 0, x3 = 0, x4 = 0))
}

# =============================================================================
# 3B. ODE TMDD — PK Michaelis-Menten + PD Simeoni
#     Parametres supplementaires : Vmax [ug/kg/j], Km [ug/L] (fixe)
#     dA1/dt = -(CL/V1+Q/V1)*A1 + (Q/V2)*A2 - Vmax*C1/(Km+C1)
#     PD identique a simeoni_rhs  <- NE PAS MODIFIER simeoni_rhs/ctrl_rhs
# =============================================================================

simeoni_rhs_tmdd <- function(t, state, parms) {
  A1 <- max(state["A1"], 0); A2 <- max(state["A2"], 0)
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)

  CL   <- parms["CL"]; V1 <- parms["V1"]
  V2   <- parms["V2"]; Q  <- parms["Q"]
  L0   <- parms["L0"]; L1 <- parms["L1"]
  k1   <- parms["k1"]; k2 <- parms["k2"]
  Vmax <- parms["Vmax"]; Km <- parms["Km"]

  C1  <- A1 / V1
  dA1 <- -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2 - Vmax * C1 / (Km + C1)
  dA2 <-  (Q/V1)*A1 - (Q/V2)*A2

  w  <- x1 + x2 + x3 + x4
  gw <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  dx1 <- (growth_rate - k2 * C1) * x1
  dx2 <- k2 * C1 * x1 - k1 * x2
  dx3 <- k1 * x2 - k1 * x3
  dx4 <- k1 * x3 - k1 * x4

  list(c(A1 = dA1, A2 = dA2, x1 = dx1, x2 = dx2, x3 = dx3, x4 = dx4))
}

# =============================================================================
# 4. FONCTIONS DE SIMULATION  <- NE PAS MODIFIER
# =============================================================================

sim_ctrl_fn <- function(L0, L1, tv0, times_out) {
  t_all <- sort(unique(c(0, times_out)))
  out <- tryCatch(
    as.data.frame(lsoda(
      y     = c(x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
      times = t_all,
      func  = simeoni_ctrl_rhs,
      parms = c(L0 = L0, L1 = L1)
    )),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(times_out)))
  w_tot <- out$x1 + out$x2 + out$x3 + out$x4
  approx(out$time, w_tot, xout = times_out, rule = 2)$y
}

sim_treated <- function(dose_ugkg, tv0, params, times_out) {
  t_all  <- sort(unique(c(0, times_out)))
  params <- setNames(as.numeric(params),
                     c("CL", "V1", "V2", "Q", "L0", "L1", "k1", "k2"))
  out <- tryCatch(
    as.data.frame(lsoda(
      y        = c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
      times    = t_all,
      func     = simeoni_rhs,
      parms    = params,
      atol     = 1e-6, rtol = 1e-6,
      maxsteps = 50000
    )),
    error   = function(e) { message("lsoda error: ", conditionMessage(e)); NULL },
    warning = function(w) {
      suppressWarnings(
        as.data.frame(lsoda(
          y        = c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
          times    = t_all,
          func     = simeoni_rhs,
          parms    = params,
          atol     = 1e-4, rtol = 1e-4,
          maxsteps = 100000
        ))
      )
    }
  )
  if (is.null(out) || any(is.na(out$x1))) return(rep(NA_real_, length(times_out)))
  w_tot <- pmax(out$x1, 0) + pmax(out$x2, 0) + pmax(out$x3, 0) + pmax(out$x4, 0)
  approx(out$time, w_tot, xout = times_out, rule = 2)$y
}

# 4B. Simulation avec TMDD (Km fixe)
sim_treated_tmdd <- function(dose_ugkg, tv0, params, times_out) {
  t_all  <- sort(unique(c(0, times_out)))
  params <- setNames(as.numeric(params),
                     c("CL", "V1", "V2", "Q", "L0", "L1", "k1", "k2", "Vmax", "Km"))
  out <- tryCatch(
    as.data.frame(lsoda(
      y        = c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
      times    = t_all,
      func     = simeoni_rhs_tmdd,
      parms    = params,
      atol     = 1e-6, rtol = 1e-6,
      maxsteps = 50000
    )),
    error   = function(e) { message("lsoda error: ", conditionMessage(e)); NULL },
    warning = function(w) {
      suppressWarnings(
        as.data.frame(lsoda(
          y        = c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
          times    = t_all,
          func     = simeoni_rhs_tmdd,
          parms    = params,
          atol     = 1e-4, rtol = 1e-4,
          maxsteps = 100000
        ))
      )
    }
  )
  if (is.null(out) || any(is.na(out$x1))) return(rep(NA_real_, length(times_out)))
  w_tot <- pmax(out$x1, 0) + pmax(out$x2, 0) + pmax(out$x3, 0) + pmax(out$x4, 0)
  approx(out$time, w_tot, xout = times_out, rule = 2)$y
}

# =============================================================================
# 5A. ETAPE A — INTOUCHABLE
#     Croissance : L0 et L1 sur le groupe vehicule uniquement.
#     Bornes L0 : doublement dans [3, 14 j]
#     Bornes L1 : [TV_max observe, 200 x TV_max]
# =============================================================================

cat("\n===========================================================\n")
cat("ETAPE A — Ajustement vehicule : L0 et L1\n")
cat("===========================================================\n")

lm_ctrl <- lm(log(tv_ctrl[ok_ctrl]) ~ times_d[ok_ctrl])
L0_init <- max(coef(lm_ctrl)[2], 0.005)

lower_A <- c(log(0.05),
             log(max(tv_ctrl[ok_ctrl], na.rm = TRUE)))
upper_A <- c(log(0.25),
             log(max(tv_ctrl[ok_ctrl], na.rm = TRUE) * 200))

objective_vehicule <- function(logpar) {
  L0 <- exp(logpar[1])
  L1 <- exp(logpar[2])
  if (L0 <= 0 || L1 <= 0) return(1e12)

  pred <- sim_ctrl_fn(L0, L1, tv0_ctrl, times_d[ok_ctrl])
  if (any(is.na(pred))) return(1e12)
  pred <- pmax(pred, 0.1)

  sum(1 / sem_ctrl[ok_ctrl]^2 * (log(tv_ctrl[ok_ctrl]) - log(pred))^2,
      na.rm = TRUE)
}

set.seed(42)
fit_de_A <- DEoptim(
  fn      = objective_vehicule,
  lower   = lower_A,
  upper   = upper_A,
  control = DEoptim.control(
    NP      = 80,
    itermax = 600,
    F       = 0.8,
    CR      = 0.9,
    trace   = 100,
    reltol  = 1e-8,
    steptol = 150
  )
)

fit_local_A <- nlminb(
  start     = fit_de_A$optim$bestmem,
  objective = objective_vehicule,
  lower     = lower_A,
  upper     = upper_A,
  control   = list(eval.max = 2000, iter.max = 800,
                   rel.tol  = 1e-12, x.tol = 1e-12)
)

if (fit_local_A$objective < fit_de_A$optim$bestval) {
  best_par_A <- fit_local_A$par
  best_val_A <- fit_local_A$objective
  src_A      <- "nlminb"
} else {
  best_par_A <- as.numeric(fit_de_A$optim$bestmem)
  best_val_A <- fit_de_A$optim$bestval
  src_A      <- "DEoptim"
}

L0_est <- exp(best_par_A[1])
L1_est <- exp(best_par_A[2])

cat(sprintf("\n-> L0 = %.6f /j   [doublement = %.1f j]  (%s)\n",
            L0_est, log(2) / L0_est, src_A))
cat(sprintf("-> L1 = %.2f mm3/j  (%s)\n", L1_est, src_A))
cat(sprintf("   Objectif DEoptim : %.8f\n", fit_de_A$optim$bestval))
cat(sprintf("   Objectif retenu  : %.8f\n", best_val_A))

# =============================================================================
# 5B. ETAPE B — Efficacite : k1 et TSC (Km et Vmax fixes)
#
#   L0_est et L1_est fixes depuis l'Etape A.
#   Km   fixe = KM_FIXED   ug/L  : regime zeroth-order
#   Vmax fixe = VMAX_FIXED ug/kg/j : non identifiable depuis donnees PD,
#     fixe pour garantir suppression soutenue a 5.6 mg/kg sur 50j.
#   k2 derive : k2 = L0_est / TSC
#
#   Bornes TSC : [Cmax_1p2/10, Cmax_5p6]
#   Bornes k1  : MTT dans [3, 25 j]
# =============================================================================

cat("\n===========================================================\n")
cat("ETAPE B — Ajustement groupes traites : k1 et TSC (Km, Vmax fixes)\n")
cat("===========================================================\n")
cat(sprintf("L0 fixe = %.6f /j  |  L1 fixe = %.2f mm3/j\n", L0_est, L1_est))
cat(sprintf("Km   fixe = %.0f ug/L\n",    KM_FIXED))
cat(sprintf("Vmax fixe = %.0f ug/kg/j\n", VMAX_FIXED))

TSC_lower <- Cmax_1p2 / 10
TSC_upper <- Cmax_5p6
k1_lower  <- 4 / 25
k1_upper  <- 4 / 3

cat(sprintf("Bornes TSC : [%.0f, %.0f] ug/L\n", TSC_lower, TSC_upper))
cat(sprintf("Bornes k1  : [%.4f, %.4f] /j    (MTT dans [3, 25 j])\n", k1_lower, k1_upper))

lower_B <- c(log(k1_lower), log(TSC_lower))
upper_B <- c(log(k1_upper), log(TSC_upper))

objective_traites <- function(logpar) {
  k1  <- exp(logpar[1])
  TSC <- exp(logpar[2])
  if (k1 <= 0 || TSC <= 0) return(1e12)

  k2 <- L0_est / TSC

  params_all <- c(
    CL = unname(pk_fixed["CL"]), V1 = unname(pk_fixed["V1"]),
    V2 = unname(pk_fixed["V2"]), Q  = unname(pk_fixed["Q"]),
    L0 = L0_est, L1 = L1_est, k1 = k1, k2 = k2,
    Vmax = VMAX_FIXED, Km = KM_FIXED
  )

  pred_d1p2 <- sim_treated_tmdd(1200,  tv0_d1p2,  params_all, times_d[ok_d1p2])
  if (any(is.na(pred_d1p2))) return(1e12)
  pred_d1p2 <- pmax(pred_d1p2, 0.1)

  pred_d5p6 <- sim_treated_tmdd(5600,  tv0_d5p6,  params_all, times_d[ok_d5p6])
  if (any(is.na(pred_d5p6))) return(1e12)
  pred_d5p6 <- pmax(pred_d5p6, 0.1)

  pred_d11p8 <- sim_treated_tmdd(11800, tv0_d11p8, params_all, times_d[ok_d11p8])
  if (any(is.na(pred_d11p8))) return(1e12)
  pred_d11p8 <- pmax(pred_d11p8, 0.1)

  sum(1 / sem_d1p2[ok_d1p2]^2   * (log(tv_d1p2[ok_d1p2])   - log(pred_d1p2))^2,  na.rm = TRUE) +
  sum(1 / sem_d5p6[ok_d5p6]^2   * (log(tv_d5p6[ok_d5p6])   - log(pred_d5p6))^2,  na.rm = TRUE) +
  sum(1 / sem_d11p8[ok_d11p8]^2 * (log(tv_d11p8[ok_d11p8]) - log(pred_d11p8))^2, na.rm = TRUE)
}

# Test de sanite avant DEoptim
params_sanity <- c(
  CL = unname(pk_fixed["CL"]), V1 = unname(pk_fixed["V1"]),
  V2 = unname(pk_fixed["V2"]), Q  = unname(pk_fixed["Q"]),
  L0 = L0_est, L1 = L1_est,
  k1 = 4 / 15,
  k2 = L0_est / 3000,
  Vmax = VMAX_FIXED, Km = KM_FIXED
)
test_ok <- sim_treated_tmdd(1200, tv0_d1p2, params_sanity, times_d[ok_d1p2])
if (all(is.na(test_ok))) {
  stop("ERREUR : sim_treated_tmdd retourne NA avec des parametres raisonnables.")
}
cat(sprintf("  [Sanity check] sim_treated_tmdd 1.2 mg/kg : OK  (pred[1] = %.1f mm3)\n",
            test_ok[1]))

set.seed(42)
fit_de_B <- DEoptim(
  fn      = objective_traites,
  lower   = lower_B,
  upper   = upper_B,
  control = DEoptim.control(
    NP      = 80,
    itermax = 800,
    F       = 0.8,
    CR      = 0.9,
    trace   = 100,
    reltol  = 1e-8,
    steptol = 200
  )
)

fit_local_B <- nlminb(
  start     = fit_de_B$optim$bestmem,
  objective = objective_traites,
  lower     = lower_B,
  upper     = upper_B,
  control   = list(eval.max = 3000, iter.max = 1000,
                   rel.tol  = 1e-12, x.tol = 1e-12)
)

if (fit_local_B$objective < fit_de_B$optim$bestval) {
  best_par_B <- fit_local_B$par
  best_val_B <- fit_local_B$objective
  src_B      <- "nlminb"
} else {
  best_par_B <- as.numeric(fit_de_B$optim$bestmem)
  best_val_B <- fit_de_B$optim$bestval
  src_B      <- "DEoptim"
}

k1_est  <- exp(best_par_B[1])
TSC_est <- exp(best_par_B[2])
k2_est  <- L0_est / TSC_est

cat(sprintf("\n-> k1  = %.5f /j   [MTT = %.1f j]  (%s)\n", k1_est, 4 / k1_est, src_B))
cat(sprintf("-> TSC = %.0f ug/L  (%s)\n", TSC_est, src_B))
cat(sprintf("   k2 derive = %.3e L/ug/j\n", k2_est))
cat(sprintf("   Objectif DEoptim : %.8f\n", fit_de_B$optim$bestval))
cat(sprintf("   Objectif retenu  : %.8f\n", best_val_B))

# Verification position des parametres sur leurs bornes
tol <- 0.02
if (abs(log(TSC_est) - log(TSC_lower)) < tol * (log(TSC_upper) - log(TSC_lower)))
  cat("  [AVERT] TSC sur borne inferieure\n")
if (abs(log(TSC_est) - log(TSC_upper)) < tol * (log(TSC_upper) - log(TSC_lower)))
  cat("  [AVERT] TSC sur borne superieure\n")
if (abs(log(k1_est)  - log(k1_upper))  < tol * (log(k1_upper)  - log(k1_lower)))
  cat(sprintf("  [AVERT] k1 sur borne superieure (MTT=%.1fj)\n", 4/k1_est))
if (abs(log(k1_est)  - log(k1_lower))  < tol * (log(k1_upper)  - log(k1_lower)))
  cat(sprintf("  [AVERT] k1 sur borne inferieure (MTT=%.1fj)\n", 4/k1_est))

# =============================================================================
# 6. RESUME FINAL
# =============================================================================

mtt <- 4 / k1_est

cat("\n===========================================================\n")
cat("PARAMETRES FINAUX — SNU  (ajustement sequentiel, TMDD Km fixe)\n")
cat("===========================================================\n")
cat(sprintf("L0   = %.6f /j         [Etape A — %s]\n", L0_est,    src_A))
cat(sprintf("L1   = %.2f mm3/j      [Etape A — %s]\n", L1_est,    src_A))
cat(sprintf("k1   = %.5f /j         [Etape B — %s]\n", k1_est,    src_B))
cat(sprintf("TSC  = %.0f ug/L       [Etape B — %s]\n", TSC_est,   src_B))
cat(sprintf("Vmax = %.0f ug/kg/j    [fixe]\n",          VMAX_FIXED))
cat(sprintf("Km   = %.0f ug/L       [fixe]\n",          KM_FIXED))
cat("-----------------------------------------------------------\n")
cat(sprintf("MTT = %.1f j\n", mtt))
cat(sprintf("k2  = %.3e L/ug/j  (= L0/TSC)\n", k2_est))
cat(sprintf("Plausibilite : Cmax(1.2mg)=%.0f | TSC=%.0f | Cmax(5.6mg)=%.0f ug/L\n",
            Cmax_1p2, TSC_est, Cmax_5p6))
cat(sprintf("  Cmax/TSC @ 1.2 mg/kg : %.1f\n", Cmax_1p2  / TSC_est))
cat(sprintf("  Cmax/TSC @ 5.6 mg/kg : %.1f\n", Cmax_5p6  / TSC_est))
cat(sprintf("  Cmax/TSC @ 11.8 mg/kg: %.1f\n", Cmax_11p8 / TSC_est))

# =============================================================================
# 7. SIMULATION FINALE
# =============================================================================

params_best <- c(
  CL = unname(pk_fixed["CL"]), V1 = unname(pk_fixed["V1"]),
  V2 = unname(pk_fixed["V2"]), Q  = unname(pk_fixed["Q"]),
  L0 = L0_est, L1 = L1_est, k1 = k1_est, k2 = k2_est,
  Vmax = Vmax_est, Km = KM_FIXED
)

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim  <- sim_ctrl_fn(L0_est, L1_est, tv0_ctrl, times_sim)
pred_d1p2_sim  <- sim_treated_tmdd(1200,  tv0_d1p2,  params_best, times_sim)
pred_d5p6_sim  <- sim_treated_tmdd(5600,  tv0_d5p6,  params_best, times_sim)
pred_d11p8_sim <- sim_treated_tmdd(11800, tv0_d11p8, params_best, times_sim)

niv <- c("Vehicule", "1.2 mg/kg", "5.6 mg/kg", "11.8 mg/kg")

df_sim <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_sim,  Groupe = "Vehicule"),
  data.frame(jour = times_sim, TV = pred_d1p2_sim,  Groupe = "1.2 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d5p6_sim,  Groupe = "5.6 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d11p8_sim, Groupe = "11.8 mg/kg")
)
df_sim$Groupe <- factor(df_sim$Groupe, levels = niv)

df_obs <- rbind(
  data.frame(jour = times_d[ok_ctrl],  TV = tv_ctrl[ok_ctrl],    sem = sem_ctrl[ok_ctrl],   Groupe = "Vehicule"),
  data.frame(jour = times_d[ok_d1p2],  TV = tv_d1p2[ok_d1p2],   sem = sem_d1p2[ok_d1p2],  Groupe = "1.2 mg/kg"),
  data.frame(jour = times_d[ok_d5p6],  TV = tv_d5p6[ok_d5p6],   sem = sem_d5p6[ok_d5p6],  Groupe = "5.6 mg/kg"),
  data.frame(jour = times_d[ok_d11p8], TV = tv_d11p8[ok_d11p8], sem = sem_d11p8[ok_d11p8], Groupe = "11.8 mg/kg")
)
df_obs$Groupe <- factor(df_obs$Groupe, levels = niv)

# =============================================================================
# 8. GRAPHIQUE
# =============================================================================

cols <- c("Vehicule"    = "#888888",
          "1.2 mg/kg"  = "#74ADD1",
          "5.6 mg/kg"  = "#4393C3",
          "11.8 mg/kg" = "#2166AC")

ymax <- max(df_obs$TV + df_obs$sem, na.rm = TRUE)

subtitle_txt <- paste0(
  "A: L0=", round(L0_est, 4), "/j | L1=", format(round(L1_est), big.mark = " "), "mm3/j",
  "  ||  ",
  "B: k1=", round(k1_est, 3), "/j | TSC=", round(TSC_est, 0), "ug/L",
  " | Vmax=", VMAX_FIXED, "(fixe) | MTT=", round(mtt, 1), "j"
)

p_simeoni <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(data    = df_sim[df_sim$Groupe != "Vehicule", ],
            mapping = aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1) +
  geom_line(data    = df_sim[df_sim$Groupe == "Vehicule", ],
            mapping = aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1, linetype = "dashed") +
  geom_errorbar(data    = df_obs,
                mapping = aes(x = jour, ymin = TV - sem, ymax = TV + sem, color = Groupe),
                width = 0.8, linewidth = 0.5) +
  geom_point(data    = df_obs,
             mapping = aes(x = jour, y = TV, color = Groupe),
             size = 2.5) +
  annotate("point", x = 0,   y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text",  x = 1.5, y = -ymax * 0.06,
           label = "= Dose unique (j0)", hjust = 0, size = 3.2, color = "#CC0000") +
  scale_x_continuous(breaks = seq(0, max(times_d, na.rm = TRUE), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "PK/PD TGI - Simeoni (2004) - Fc-silent FGFR2-huBPA-LP1 - SNU",
    subtitle = subtitle_txt,
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm3)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 9, color = "grey50"),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

print(p_simeoni)
ggsave("scripts/plot_PKPD_simeoni_SNU.png", p_simeoni,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique sauvegarde : scripts/plot_PKPD_simeoni_SNU.png\n")

# =============================================================================
# 9. SAUVEGARDE
# =============================================================================

simeoni_results_SNU <- list(
  pk_fixed      = pk_fixed,
  L0            = L0_est,
  L1            = L1_est,
  k1            = k1_est,
  TSC_ugL       = TSC_est,
  Vmax_ugkgj    = VMAX_FIXED,
  Km_ugL        = KM_FIXED,
  k2            = k2_est,
  MTT_days      = mtt,
  obj_A_DEoptim = fit_de_A$optim$bestval,
  obj_A_final   = best_val_A,
  obj_B_DEoptim = fit_de_B$optim$bestval,
  obj_B_final   = best_val_B
)

save(simeoni_results_SNU, file = "scripts/resultats_PKPD_simeoni_SNU.RData")
cat("Resultats sauvegardes : scripts/resultats_PKPD_simeoni_SNU.RData\n")
