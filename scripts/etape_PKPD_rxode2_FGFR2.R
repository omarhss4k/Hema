# =============================================================================
# PKPD — Fc-silent FGFR2-huBPA-LP1 — Simeoni 2004 — rxode2
# Version corrigée :
#   - start_grid k2 recalibré pour doses en µg/kg
#   - approx rule=1 + vérification intégration complète
#   - TGI : gestion régression tumorale et dénominateur nul
#   - Stratégie : fit global (k1, k2) puis k2 libre par scénario (k1 fixé)
# =============================================================================

library(rxode2)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS
# =============================================================================

load("resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2

PK <- c(
  CL = unname(pk2comp_rxode2$CL),
  V1 = unname(pk2comp_rxode2$V1),
  V2 = unname(pk2comp_rxode2$V2),
  Q  = unname(pk2comp_rxode2$Q)
)

cat("=== Paramètres PK fixés ===\n")
cat("CL =", round(PK["CL"], 6), "L/h/kg  =", round(PK["CL"]*24, 4), "L/j/kg\n")
cat("V1 =", round(PK["V1"], 5), "L/kg\n")
cat("V2 =", round(PK["V2"], 5), "L/kg\n")
cat("Q  =", round(PK["Q"],  6), "L/h/kg\n")

# =============================================================================
# 2. PARAMÈTRES PD DE CROISSANCE FIXÉS
# =============================================================================

l0 <- 0.146 / 24   # /h
l1 <- 0.334 / 24   # g/h
p  <- 20

# C1_max pour 10 mg/kg → ordre de grandeur k2 attendu
C1_max_10 <- 10000 / PK["V1"]
k2_ref    <- l0 / C1_max_10   # effet ~l0 au pic : k2 * C1_max ≈ l0
cat("\nRepère k2 : C1_max(10 mg/kg) =", round(C1_max_10, 0),
    "µg/L  →  k2_ref ≈", formatC(k2_ref, format="e", digits=2), "\n")

# =============================================================================
# 3. DONNÉES TUMORALES (mm³ → g)
# =============================================================================

raw_pd  <- read_xlsx("TumorVolume_FGFR2.xlsx",
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

cat("\nPoints par groupe :\n")
cat("  Contrôle  :", nrow(dat_ctrl), "\n")
cat("  10 mg/kg  :", nrow(dat_d10),  "\n")
cat("   3 mg/kg  :", nrow(dat_d3),   "\n")

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
cat("\nw0 :", round(w0, 4), "g\n")

# =============================================================================
# 5. MODÈLE PKPD rxode2 — Simeoni 2004
# =============================================================================

pkpd_model <- suppressMessages(rxode2({
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
  d/dt(A2) <-  (Q/V1) * A1 - (Q/V2) * A2

  w      <- x1 + x2 + x3 + x4
  growth <- l0 * x1 / (1 + (l0/l1 * w)^p)^(1/p)

  d/dt(x1) <- growth - k2 * C1 * x1
  d/dt(x2) <- k2 * C1 * x1 - k1 * x2
  d/dt(x3) <- k1 * (x2 - x3)
  d/dt(x4) <- k1 * (x3 - x4)
}))

# =============================================================================
# 6. FONCTIONS DE SIMULATION — rule=1 + vérification intégration complète
# =============================================================================

inits_base <- c(A1 = 0, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)

.sim <- function(ev, params, obs_times) {
  out <- tryCatch(
    rxSolve(pkpd_model, as.list(params), ev, inits = inits_base),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(obs_times)))

  # Intégration incomplète → NA plutôt qu'extrapolation silencieuse
  if (max(out$time) < max(obs_times) - 1)
    return(rep(NA_real_, length(obs_times)))

  w_tot  <- pmax(out$x1 + out$x2 + out$x3 + out$x4, 1e-9)
  result <- approx(out$time, w_tot, xout = obs_times, rule = 1)$y
  if (any(is.na(result))) return(rep(NA_real_, length(obs_times)))
  result
}

simulate_single <- function(dose, params, obs_times, t_admin = 0) {
  ev <- eventTable()
  if (dose > 0)
    ev$add.dosing(dose = dose, nbr.doses = 1,
                  dosing.to = 1, start.time = t_admin)
  ev$add.sampling(sort(unique(c(0, obs_times))))
  .sim(ev, params, obs_times)
}

simulate_multi <- function(dose, params, obs_times,
                           t_admin = 0, n_doses = 4, interval_h = 14*24) {
  ev <- eventTable()
  if (dose > 0)
    ev$add.dosing(dose = dose, nbr.doses = n_doses,
                  dosing.interval = interval_h,
                  dosing.to = 1, start.time = t_admin)
  ev$add.sampling(sort(unique(c(0, obs_times))))
  .sim(ev, params, obs_times)
}

# =============================================================================
# 7. FONCTION OBJECTIVE GÉNÉRIQUE
# =============================================================================

params_fixed <- c(PK, l0 = l0, l1 = l1, p = p)

.residuals <- function(sim_fn, par) {
  pc  <- tryCatch(sim_fn(dose0,  par, dat_ctrl$t), error = function(e) NULL)
  p3  <- tryCatch(sim_fn(dose3,  par, dat_d3$t),   error = function(e) NULL)
  p10 <- tryCatch(sim_fn(dose10, par, dat_d10$t),  error = function(e) NULL)
  bad <- function(p) is.null(p) || any(is.na(p) | p <= 0)
  if (bad(pc) || bad(p3) || bad(p10)) return(1e10)
  val <- mean((log(dat_ctrl$w) - log(pc))^2)  +
         mean((log(dat_d3$w)   - log(p3))^2)  +
         mean((log(dat_d10$w)  - log(p10))^2)
  if (!is.finite(val)) 1e10 else val
}

# Objectif pour un scénario unique (k1 + k2 libres)
make_objective <- function(sim_fn) {
  function(logpar) {
    par <- c(params_fixed,
             k1 = unname(exp(logpar[1])),
             k2 = unname(exp(logpar[2])))
    .residuals(sim_fn, par)
  }
}

# Objectif k2 seul (k1 fixé)
make_objective_k2 <- function(sim_fn, k1_fixed) {
  function(logpar_k2) {
    par <- c(params_fixed,
             k1 = k1_fixed,
             k2 = unname(exp(logpar_k2[1])))
    .residuals(sim_fn, par)
  }
}

# Objectif global : mêmes k1/k2 doivent expliquer les deux scénarios
make_objective_global <- function() {
  function(logpar) {
    par <- c(params_fixed,
             k1 = unname(exp(logpar[1])),
             k2 = unname(exp(logpar[2])))
    rs  <- .residuals(simulate_single, par)
    rm  <- .residuals(simulate_multi,  par)
    if (rs >= 1e10 || rm >= 1e10) return(1e10)
    (rs + rm) / 2   # moyenne pour équipondérer les deux scénarios
  }
}

# =============================================================================
# 8. GRILLES DE DÉPART — k2 recalibré pour doses en µg/kg
#    k2_ref ≈ l0 / C1_max(10 mg/kg)  →  ordre 1e-8 à 1e-6
# =============================================================================

# C1_max ≈ 143 000 µg/L pour V1 ≈ 0.07 L/kg
# k2 * C1_max doit être ~ l0 (0.006/h) : k2 ≈ 4e-8
start_grid <- list(
  c(k1 = 0.5,  k2 = 1e-7),
  c(k1 = 0.1,  k2 = 1e-8),
  c(k1 = 1.0,  k2 = 5e-8),
  c(k1 = 0.2,  k2 = 1e-6),
  c(k1 = 0.5,  k2 = 1e-9)
)

run_optim <- function(obj_fn, starts, label) {
  cat("\nOptimisation —", label, "...\n")
  best_obj <- Inf;  best_fit <- NULL
  for (s in starts) {
    fit_try <- tryCatch(
      nlminb(log(s), obj_fn,
             control = list(eval.max = 3000, iter.max = 1500,
                            rel.tol = 1e-12, x.tol = 1e-12)),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(fit_try$objective) &&
        fit_try$objective < best_obj) {
      best_obj <- fit_try$objective;  best_fit <- fit_try
    }
  }
  if (is.null(best_fit)) stop(paste("Aucune convergence pour", label))
  best_fit
}

# =============================================================================
# 9. PHASE 1 — FIT GLOBAL (k1, k2 communs aux deux scénarios)
# =============================================================================

fit_global <- run_optim(make_objective_global(), start_grid, "Fit global")

k1_global <- unname(exp(fit_global$par[1]))
k2_global <- unname(exp(fit_global$par[2]))

cat("\n=== Fit global (k1, k2 communs) ===\n")
cat("  k1 =", round(k1_global, 6), "/h\n")
cat("  k2 =", formatC(k2_global, format="e", digits=3), "\n")
cat("  Objectif :", round(fit_global$objective, 5), "\n")

# =============================================================================
# 10. PHASE 2 — k2 LIBRE PAR SCÉNARIO (k1 = k1_global fixé)
#     Permet de détecter résistance ou accumulation
# =============================================================================

starts_k2 <- lapply(start_grid, function(s) s["k2"])

fit_k2_single <- run_optim(make_objective_k2(simulate_single, k1_global),
                           starts_k2, "k2 dose unique (k1 fixé)")
fit_k2_multi  <- run_optim(make_objective_k2(simulate_multi,  k1_global),
                           starts_k2, "k2 doses répétées (k1 fixé)")

k2_single <- unname(exp(fit_k2_single$par[1]))
k2_multi  <- unname(exp(fit_k2_multi$par[1]))

cat("\n=== Fits séparés k2 (k1 =", round(k1_global, 5), "/h fixé) ===\n")
cat("  k2 dose unique    :", formatC(k2_single, format="e", digits=3), "\n")
cat("  k2 doses répétées :", formatC(k2_multi,  format="e", digits=3), "\n")
ratio <- k2_multi / k2_single
cat("  Ratio k2_multi/k2_single :", round(ratio, 3),
    if (ratio > 1.5) "→ accumulation possible" else
    if (ratio < 0.67) "→ résistance possible"  else
    "→ pas de dérive notable", "\n")

# =============================================================================
# 11. GRAPHIQUES
# =============================================================================

times_full <- seq(0, 49 * 24, by = 1)
lev        <- c("Contrôle", "3 mg/kg", "10 mg/kg")
doses      <- list(dose0, dose3, dose10)
grps       <- c("Contrôle", "3 mg/kg", "10 mg/kg")

make_df_sim <- function(sim_fn, k1, k2) {
  par <- c(params_fixed, k1 = k1, k2 = k2)
  do.call(rbind, mapply(function(d, g)
    data.frame(t = times_full/24,
               w = sim_fn(d, par, times_full),
               Groupe = g),
    doses, grps, SIMPLIFY = FALSE))
}

df_obs <- rbind(
  data.frame(t = dat_ctrl$t/24, w = dat_ctrl$w, Groupe = "Contrôle"),
  data.frame(t = dat_d3$t/24,   w = dat_d3$w,   Groupe = "3 mg/kg"),
  data.frame(t = dat_d10$t/24,  w = dat_d10$w,  Groupe = "10 mg/kg")
)
df_obs$Groupe <- factor(df_obs$Groupe, levels = lev)

plot_pkpd <- function(df_sim, title_suffix, dose_days, k1, k2) {
  df_sim$Groupe <- factor(df_sim$Groupe, levels = lev)
  ggplot() +
    geom_line(data  = df_sim, aes(x=t, y=w, color=Groupe), linewidth=1) +
    geom_point(data = df_obs, aes(x=t, y=w, color=Groupe), size=2.5) +
    geom_vline(xintercept = dose_days, linetype="dashed",
               color="grey60", linewidth=0.5) +
    annotate("text", x = dose_days, y = max(df_obs$w, na.rm=TRUE)*1.05,
             label = paste0("j", dose_days), size=2.8,
             color="grey40", hjust=0.5) +
    labs(
      title    = paste("PKPD Simeoni 2004 (rxode2) —", title_suffix),
      subtitle = paste0("Fc-silent FGFR2-huBPA-LP1 | ",
                        "k1=", round(k1, 4), " /h  ",
                        "k2=", formatC(k2, format="e", digits=2)),
      x = "Temps (jours)", y = "Volume tumoral (g)"
    ) +
    theme_bw(base_size = 13)
}

# Graphique fit global
df_sim_global_s <- make_df_sim(simulate_single, k1_global, k2_global)
plot_pkpd(df_sim_global_s, "Dose unique — fit global",
          dose_days = 0, k1_global, k2_global)
ggsave("scripts/plot_PKPD_global_single_FGFR2.png", width=9, height=5, dpi=150)

df_sim_global_m <- make_df_sim(simulate_multi, k1_global, k2_global)
plot_pkpd(df_sim_global_m, "Doses répétées — fit global",
          dose_days = c(0,14,28,42), k1_global, k2_global)
ggsave("scripts/plot_PKPD_global_multi_FGFR2.png", width=9, height=5, dpi=150)

# Graphique fits séparés k2
df_sim_k2s <- make_df_sim(simulate_single, k1_global, k2_single)
plot_pkpd(df_sim_k2s, "Dose unique — k2 libre",
          dose_days = 0, k1_global, k2_single)
ggsave("scripts/plot_PKPD_k2sep_single_FGFR2.png", width=9, height=5, dpi=150)

df_sim_k2m <- make_df_sim(simulate_multi, k1_global, k2_multi)
plot_pkpd(df_sim_k2m, "Doses répétées — k2 libre",
          dose_days = c(0,14,28,42), k1_global, k2_multi)
ggsave("scripts/plot_PKPD_k2sep_multi_FGFR2.png", width=9, height=5, dpi=150)

cat("\nGraphiques → scripts/plot_PKPD_*_FGFR2.png\n")

# =============================================================================
# 12. TGI — gestion régression tumorale et dénominateur nul
# =============================================================================

tgi <- function(w_ctrl_end, w_treat_end, w0_val) {
  delta_ctrl <- w_ctrl_end - w0_val
  if (is.na(delta_ctrl) || delta_ctrl <= 0) return(NA_real_)  # contrôle n'a pas crû
  round((1 - (w_treat_end - w0_val) / delta_ctrl) * 100, 1)
  # Note : valeur > 100% = régression tumorale (W_treat < W0)
  #        valeur < 0%   = croissance malgré traitement
}

calc_tgi <- function(sim_fn, k1, k2, label) {
  par   <- c(params_fixed, k1 = k1, k2 = k2)
  t_seq <- seq(0, 49*24, by = 24)
  wc  <- tail(sim_fn(dose0,  par, t_seq), 1)
  w3  <- tail(sim_fn(dose3,  par, t_seq), 1)
  w10 <- tail(sim_fn(dose10, par, t_seq), 1)

  tgi3  <- tgi(wc, w3,  w0)
  tgi10 <- tgi(wc, w10, w0)

  cat("\nTGI à j49 —", label, ":\n")
  cat("   3 mg/kg  :", if (is.na(tgi3))  "NA" else
      if (tgi3 > 100) paste0(tgi3, "% (régression)") else paste0(tgi3, "%"), "\n")
  cat("  10 mg/kg  :", if (is.na(tgi10)) "NA" else
      if (tgi10 > 100) paste0(tgi10, "% (régression)") else paste0(tgi10, "%"), "\n")
}

cat("\n--- TGI fit global ---\n")
calc_tgi(simulate_single, k1_global, k2_global, "Dose unique")
calc_tgi(simulate_multi,  k1_global, k2_global, "Doses répétées")

cat("\n--- TGI k2 libre par scénario ---\n")
calc_tgi(simulate_single, k1_global, k2_single, "Dose unique")
calc_tgi(simulate_multi,  k1_global, k2_multi,  "Doses répétées")

# =============================================================================
# 13. SAUVEGARDE
# =============================================================================

save(k1_global, k2_global, k2_single, k2_multi,
     params_fixed, w0, dat_ctrl, dat_d3, dat_d10,
     file = "scripts/resultats_PKPD_rxode2_FGFR2.RData")
cat("\nRésultats → scripts/resultats_PKPD_rxode2_FGFR2.RData\n")
