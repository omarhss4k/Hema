# =============================================================================
# PKPD — Fc-silent FGFR2-huBPA-LP1 — Simeoni 2004 — deSolve
# -----------------------------------------------------------------------------
# Groupes : Group 01 = Véhicule (contrôle)
#           Group 03 = 10 mg/kg
#           Group 04 =  3 mg/kg
# Observation : 49 jours
# Deux scénarios : (A) dose unique  (B) doses répétées toutes les 14j
# PK fixé depuis estimation 2-comp (méthode des résidus)
# PD estimé : k1, k2
# =============================================================================

library(deSolve)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS
# =============================================================================

load("scripts/resultats_PK2comp_FGFR2.RData")   # → pk2comp

CL <- unname(pk2comp$CL)   # L/h/kg
V1 <- unname(pk2comp$V1)   # L/kg
V2 <- unname(pk2comp$V2)   # L/kg
Q  <- unname(pk2comp$Q)    # L/h/kg

k10 <- CL / V1
k12 <- Q  / V1
k21 <- Q  / V2

cat("=== Paramètres PK fixés ===\n")
cat("CL =", round(CL, 6), "L/h/kg  =", round(CL*24, 4), "L/j/kg\n")
cat("V1 =", round(V1, 5), "L/kg\n")
cat("V2 =", round(V2, 5), "L/kg\n")
cat("Q  =", round(Q,  6), "L/h/kg\n")
cat("k10 =", round(k10, 6), "/h\n")
cat("k12 =", round(k12, 6), "/h\n")
cat("k21 =", round(k21, 6), "/h\n")

# =============================================================================
# 2. PARAMÈTRES PD DE CROISSANCE FIXÉS (convertis en heures)
# =============================================================================

l0 <- 0.146 / 24   # /h
l1 <- 0.334 / 24   # g/h
p  <- 20

# =============================================================================
# 3. DONNÉES TUMORALES
# =============================================================================

raw_pd  <- read_xlsx("data/TumorVolume_FGFR2.xlsx",
                     n_max     = 3,
                     col_types = c("text", rep("numeric", 15)))

time_d <- suppressWarnings(as.numeric(colnames(raw_pd)[-1]))
time_h <- time_d * 24

extract_group <- function(raw, pattern) {
  idx <- grep(pattern, raw[[1]], ignore.case = TRUE)
  if (length(idx) == 0) stop(paste("Groupe introuvable :", pattern))
  as.numeric(unlist(raw[idx[1], -1]))
}

w_ctrl <- extract_group(raw_pd, "Group 01")
w_d10  <- extract_group(raw_pd, "Group 03")
w_d3   <- extract_group(raw_pd, "Group 04")

mk <- function(t, w) {
  ok <- !is.na(w) & w > 0
  data.frame(t = t[ok], w = w[ok] / 1000)   # mm³ → g
}
dat_ctrl <- mk(time_h, w_ctrl)
dat_d10  <- mk(time_h, w_d10)
dat_d3   <- mk(time_h, w_d3)

cat("\nPoints par groupe (données converties mm³→g) :\n")
cat("  Contrôle  :", nrow(dat_ctrl),
    "  w_j0=", round(dat_ctrl$w[1],4), "g  w_fin=", round(tail(dat_ctrl$w,1),4), "g\n")
cat("  10 mg/kg  :", nrow(dat_d10),
    "  w_j0=", round(dat_d10$w[1],4),  "g  w_fin=", round(tail(dat_d10$w,1),4),  "g\n")
cat("   3 mg/kg  :", nrow(dat_d3),
    "  w_j0=", round(dat_d3$w[1],4),   "g  w_fin=", round(tail(dat_d3$w,1),4),   "g\n")

dose10 <- 10 * 1000   # µg/kg
dose3  <-  3 * 1000
dose0  <-  0

# =============================================================================
# 4. MASSE TUMORALE INITIALE w0
# =============================================================================

w0 <- mean(c(
  dat_ctrl$w[which.min(dat_ctrl$t)],
  dat_d10$w[which.min(dat_d10$t)],
  dat_d3$w[which.min(dat_d3$t)]
), na.rm = TRUE)

cat("\nw0 (moyenne premiers points) :", round(w0, 4), "g\n")

# =============================================================================
# 4b. CALIBRATION l0 / l1 SUR LE GROUPE CONTRÔLE
# -----------------------------------------------------------------------------
# Pour dose=0 : C1=0 → x2=x3=x4=0 toujours → w = x1
# Modèle réduit 1D : dw/dt = l0*w / (1+(l0/l1*w)^p)^(1/p)
# =============================================================================

growth_1d <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    list(c(l0 * w / (1 + (l0/l1 * w)^p)^(1/p)))
  })
}

obj_ctrl_growth <- function(logpar) {
  l0_f <- unname(exp(logpar[1]))
  l1_f <- unname(exp(logpar[2]))
  pars  <- c(l0=l0_f, l1=l1_f, p=p)
  times <- sort(unique(c(0, dat_ctrl$t)))

  out <- tryCatch(
    as.data.frame(lsoda(c(w=w0), times, growth_1d, pars,
                        rtol=1e-6, atol=1e-8)),
    error = function(e) NULL
  )
  if (is.null(out)) return(1e10)
  w_pred <- pmax(approx(out$time, out$w, xout=dat_ctrl$t, rule=2)$y, 1e-9)
  val <- mean((log(dat_ctrl$w) - log(w_pred))^2)
  if (!is.finite(val)) 1e10 else val
}

fit_ctrl_growth <- nlminb(
  start     = log(c(l0=l0, l1=l1)),
  objective = obj_ctrl_growth,
  control   = list(eval.max=2000, iter.max=1000, rel.tol=1e-12, x.tol=1e-12)
)
l0 <- unname(exp(fit_ctrl_growth$par[1]))
l1 <- unname(exp(fit_ctrl_growth$par[2]))

cat("\n=== l0/l1 calibrés sur le groupe contrôle ===\n")
cat("l0 =", round(l0*24, 5), "/j  (", round(l0, 7), "/h)  t½ expo =",
    round(log(2)/(l0*24), 1), "j\n")
cat("l1 =", round(l1*24, 5), "g/j  (", round(l1, 7), "g/h)\n")
cat("Objectif contrôle :", round(fit_ctrl_growth$objective, 6), "\n")

# =============================================================================
# 5. MODÈLE ODE — Simeoni 2004
# =============================================================================

pkpd_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    dA1 <- -(k10 + k12) * A1 + k21 * A2
    dA2 <-   k12 * A1   - k21 * A2
    C1  <- A1 / V1

    w      <- x1 + x2 + x3 + x4
    growth <- l0 * x1 / (1 + (l0/l1 * w)^p)^(1/p)

    dx1 <- growth - k2 * C1 * x1
    dx2 <- k2 * C1 * x1 - k1 * x2
    dx3 <- k1 * (x2 - x3)
    dx4 <- k1 * (x3 - x4)

    list(c(dA1, dA2, dx1, dx2, dx3, dx4))
  })
}

pars_fixed <- c(k10=k10, k12=k12, k21=k21, V1=V1,
                l0=l0, l1=l1, p=p)

# =============================================================================
# 6. FONCTIONS DE SIMULATION
# =============================================================================

# --- (A) Dose unique à t=0 ---
simulate_single <- function(dose, k1, k2, obs_times) {
  pars   <- c(pars_fixed, k1=k1, k2=k2)
  state0 <- c(A1=dose, A2=0, x1=w0, x2=0, x3=0, x4=0)
  times  <- sort(unique(c(0, obs_times)))

  out <- tryCatch(
    as.data.frame(lsoda(state0, times, pkpd_ode, pars,
                        rtol=1e-6, atol=1e-8)),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(obs_times)))
  w_total <- pmax(out$x1 + out$x2 + out$x3 + out$x4, 1e-9)
  approx(out$time, w_total, xout=obs_times, rule=2)$y
}

# --- (B) Doses répétées toutes les 14 jours ---
simulate_multi <- function(dose, k1, k2, obs_times,
                           n_doses=4, interval_h=14*24) {
  pars    <- c(pars_fixed, k1=k1, k2=k2)
  state0  <- c(A1=dose, A2=0, x1=w0, x2=0, x3=0, x4=0)
  t_doses <- seq(0, (n_doses-1)*interval_h, by=interval_h)
  t_end   <- max(obs_times)
  times   <- sort(unique(c(0, obs_times, t_doses)))

  # événements de dosage (impulsion sur A1 à chaque administration)
  events <- data.frame(
    var   = "A1",
    time  = t_doses[-1],    # première dose déjà dans state0
    value = dose,
    method = "add"
  )

  out <- tryCatch(
    as.data.frame(lsoda(state0, times, pkpd_ode, pars,
                        events = if (nrow(events) > 0) list(data=events) else NULL,
                        rtol=1e-6, atol=1e-8)),
    error = function(e) { message("lsoda error: ", e$message); NULL }
  )
  if (is.null(out)) return(rep(NA_real_, length(obs_times)))
  w_total <- pmax(out$x1 + out$x2 + out$x3 + out$x4, 1e-9)
  approx(out$time, w_total, xout=obs_times, rule=2)$y
}

# =============================================================================
# 7. FONCTION OBJECTIVE
# =============================================================================

make_objective <- function(sim_fn) {
  function(logpar) {
    k1 <- unname(exp(logpar[1]))
    k2 <- unname(exp(logpar[2]))

    pc  <- tryCatch(sim_fn(dose0,  k1, k2, dat_ctrl$t), error=function(e) NULL)
    p3  <- tryCatch(sim_fn(dose3,  k1, k2, dat_d3$t),   error=function(e) NULL)
    p10 <- tryCatch(sim_fn(dose10, k1, k2, dat_d10$t),  error=function(e) NULL)

    bad <- function(p) is.null(p) || any(is.na(p) | p <= 0)
    if (bad(pc) || bad(p3) || bad(p10)) return(1e10)

    val <- mean((log(dat_ctrl$w) - log(pc))^2)  +
           mean((log(dat_d3$w)   - log(p3))^2)  +
           mean((log(dat_d10$w)  - log(p10))^2)

    if (!is.finite(val)) 1e10 else val
  }
}

# =============================================================================
# 8. OPTIMISATION — multi-start nlminb
# =============================================================================

start_grid <- list(
  c(k1=0.5,  k2=1e-4),
  c(k1=0.1,  k2=1e-5),
  c(k1=1.0,  k2=1e-4),
  c(k1=0.5,  k2=1e-3),
  c(k1=0.2,  k2=5e-5)
)

run_optim <- function(obj_fn, label) {
  cat("\nOptimisation —", label, "...\n")
  best_obj <- Inf
  best_fit <- NULL
  for (s in start_grid) {
    fit_try <- tryCatch(
      nlminb(log(s), obj_fn,
             control=list(eval.max=3000, iter.max=1500,
                          rel.tol=1e-12, x.tol=1e-12)),
      error=function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(fit_try$objective) &&
        fit_try$objective < best_obj) {
      best_obj <- fit_try$objective
      best_fit <- fit_try
    }
  }
  if (is.null(best_fit)) stop(paste("Aucune convergence pour", label))
  cat("  k1 =", round(unname(exp(best_fit$par[1])), 6), "/h\n")
  cat("  k2 =", round(unname(exp(best_fit$par[2])), 8), "\n")
  cat("  Objectif :", round(best_fit$objective, 5), "\n")
  best_fit
}

fit_single <- run_optim(make_objective(simulate_single), "Dose unique")
fit_multi  <- run_optim(make_objective(simulate_multi),  "Doses répétées (14j)")

k1_single <- unname(exp(fit_single$par[1]));  k2_single <- unname(exp(fit_single$par[2]))
k1_multi  <- unname(exp(fit_multi$par[1]));   k2_multi  <- unname(exp(fit_multi$par[2]))

cat("\n=== Résultats PKPD ===\n")
cat("Dose unique      : k1 =", round(k1_single, 6), "/h  k2 =", round(k2_single, 8), "\n")
cat("Doses répétées   : k1 =", round(k1_multi,  6), "/h  k2 =", round(k2_multi,  8), "\n")

# =============================================================================
# 9. GRAPHIQUES
# =============================================================================

times_full <- seq(0, 49 * 24, by = 1)
lev        <- c("Contrôle", "3 mg/kg", "10 mg/kg")

make_df_sim <- function(sim_fn, k1, k2) {
  do.call(rbind, list(
    data.frame(t=times_full/24, w=sim_fn(dose0,  k1, k2, times_full), Groupe="Contrôle"),
    data.frame(t=times_full/24, w=sim_fn(dose3,  k1, k2, times_full), Groupe="3 mg/kg"),
    data.frame(t=times_full/24, w=sim_fn(dose10, k1, k2, times_full), Groupe="10 mg/kg")
  ))
}

df_obs <- rbind(
  data.frame(t=dat_ctrl$t/24, w=dat_ctrl$w, Groupe="Contrôle"),
  data.frame(t=dat_d3$t/24,   w=dat_d3$w,   Groupe="3 mg/kg"),
  data.frame(t=dat_d10$t/24,  w=dat_d10$w,  Groupe="10 mg/kg")
)
df_obs$Groupe <- factor(df_obs$Groupe, levels=lev)

plot_pkpd <- function(df_sim, title_suffix, dose_days) {
  df_sim$Groupe <- factor(df_sim$Groupe, levels=lev)
  ggplot() +
    geom_line(data=df_sim,  aes(x=t, y=w, color=Groupe), linewidth=1) +
    geom_point(data=df_obs, aes(x=t, y=w, color=Groupe), size=2.5) +
    geom_vline(xintercept=dose_days, linetype="dashed",
               color="grey60", linewidth=0.5) +
    annotate("text", x=dose_days, y=max(df_obs$w, na.rm=TRUE)*1.05,
             label=paste0("j", dose_days), size=2.8, color="grey40", hjust=0.5) +
    labs(
      title    = paste("PKPD Simeoni 2004 (deSolve) —", title_suffix),
      subtitle = paste0("Fc-silent FGFR2-huBPA-LP1 | ",
                        "PK fixé : CL=", round(CL*24, 4), " L/j/kg"),
      x = "Temps (jours)", y = "Masse tumorale (g)"
    ) +
    theme_bw(base_size=13)
}

df_sim_single <- make_df_sim(simulate_single, k1_single, k2_single)
plot_pkpd(df_sim_single, "Dose unique", dose_days=0)
ggsave("scripts/plot_PKPD_single_FGFR2.png", width=9, height=5, dpi=150)
cat("\nGraphique dose unique → scripts/plot_PKPD_single_FGFR2.png\n")

df_sim_multi <- make_df_sim(simulate_multi, k1_multi, k2_multi)
plot_pkpd(df_sim_multi, "Doses répétées (toutes les 14j)", dose_days=c(0,14,28,42))
ggsave("scripts/plot_PKPD_multi_FGFR2.png", width=9, height=5, dpi=150)
cat("Graphique doses répétées → scripts/plot_PKPD_multi_FGFR2.png\n")

# =============================================================================
# 10. TGI — Tumor Growth Inhibition à j49
# =============================================================================

tgi <- function(w_ctrl_end, w_treat_end, w0_val)
  round((1 - (w_treat_end - w0_val) / (w_ctrl_end - w0_val)) * 100, 1)

calc_tgi <- function(sim_fn, k1, k2, label) {
  t_seq <- seq(0, 49*24, by=24)
  wc  <- tail(sim_fn(dose0,  k1, k2, t_seq), 1)
  w3  <- tail(sim_fn(dose3,  k1, k2, t_seq), 1)
  w10 <- tail(sim_fn(dose10, k1, k2, t_seq), 1)
  cat("\nTGI à j49 —", label, ":\n")
  cat("   3 mg/kg  :", tgi(wc, w3,  w0), "%\n")
  cat("  10 mg/kg  :", tgi(wc, w10, w0), "%\n")
}

calc_tgi(simulate_single, k1_single, k2_single, "Dose unique")
calc_tgi(simulate_multi,  k1_multi,  k2_multi,  "Doses répétées")

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

save(k1_single, k2_single, k1_multi, k2_multi,
     pars_fixed, w0, dat_ctrl, dat_d3, dat_d10,
     file="scripts/resultats_PKPD_desolve_FGFR2.RData")
cat("\nRésultats → scripts/resultats_PKPD_desolve_FGFR2.RData\n")
