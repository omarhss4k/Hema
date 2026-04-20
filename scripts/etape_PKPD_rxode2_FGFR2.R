# =============================================================================
# PKPD — Fc-silent FGFR2-huBPA-LP1 — Simeoni 2004 — rxode2
# -----------------------------------------------------------------------------
# Groupes : Group 01 = Véhicule (contrôle)
#           Group 03 = 10 mg/kg
#           Group 04 =  3 mg/kg
# Observation : 49 jours
# Deux scénarios : (A) dose unique  (B) doses répétées toutes les 14j
# PK fixé depuis estimation 2-comp rxode2
# PD estimé : k1, k2
# =============================================================================

library(rxode2)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2

PK <- c(
  CL = unname(pk2comp_rxode2$CL),   # L/h/kg
  V1 = unname(pk2comp_rxode2$V1),   # L/kg
  V2 = unname(pk2comp_rxode2$V2),   # L/kg
  Q  = unname(pk2comp_rxode2$Q)     # L/h/kg
)

cat("=== Paramètres PK fixés ===\n")
cat("CL =", round(PK["CL"], 6), "L/h/kg  =", round(PK["CL"]*24, 4), "L/j/kg\n")
cat("V1 =", round(PK["V1"], 5), "L/kg\n")
cat("V2 =", round(PK["V2"], 5), "L/kg\n")
cat("Q  =", round(PK["Q"],  6), "L/h/kg\n")

# =============================================================================
# 2. PARAMÈTRES PD DE CROISSANCE FIXÉS (stagiaire — convertis en heures)
# =============================================================================

l0 <- 0.146 / 24   # /h
l1 <- 0.334 / 24   # g/h
p  <- 20

# =============================================================================
# 3. DONNÉES TUMORALES
# -----------------------------------------------------------------------------
# Format Excel (wide) :
#   Ligne 1  : "Group" | 0 | 4 | 7 | 11 | 14 | 18 | 21 | 25 | 28 | ...  (jours)
#   Group 01 : Véhicule (ADC), 0 mg/kg
#   Group 03 : Fc-silent FGFR2-huBPA-LP1, 10 mg/kg
#   Group 04 : Fc-silent FGFR2-huBPA-LP1,  3 mg/kg
# =============================================================================

raw_pd  <- read_xlsx("data/TumorVolume_FGFR2.xlsx",
                     n_max     = 3,
                     col_types = c("text", rep("numeric", 15)))

# Noms de colonnes → temps en jours
time_d <- suppressWarnings(as.numeric(colnames(raw_pd)[-1]))
time_h <- time_d * 24

# Extraction robuste : correspondance partielle insensible à la casse
extract_group <- function(raw, pattern) {
  idx <- grep(pattern, raw[[1]], ignore.case = TRUE)
  if (length(idx) == 0) stop(paste("Groupe introuvable :", pattern))
  row <- raw[idx[1], -1]
  as.numeric(unlist(row))
}

w_ctrl <- extract_group(raw_pd, "Group 01")
w_d10  <- extract_group(raw_pd, "Group 03")
w_d3   <- extract_group(raw_pd, "Group 04")

mk <- function(t, w) {
  ok <- !is.na(w) & w > 0
  data.frame(t = t[ok], w = w[ok])
}
dat_ctrl <- mk(time_h, w_ctrl)
dat_d10  <- mk(time_h, w_d10)
dat_d3   <- mk(time_h, w_d3)

cat("\nPoints par groupe :\n")
cat("  Contrôle  :", nrow(dat_ctrl), "\n")
cat("  10 mg/kg  :", nrow(dat_d10),  "\n")
cat("   3 mg/kg  :", nrow(dat_d3),   "\n")

# Doses (µg/kg)
dose10 <- 10 * 1000
dose3  <-  3 * 1000
dose0  <-  0

# =============================================================================
# 4. MASSE TUMORALE INITIALE w0
#    Premier point d'observation (t = 0)
# =============================================================================

w0 <- mean(c(
  dat_ctrl$w[which.min(dat_ctrl$t)],
  dat_d10$w[which.min(dat_d10$t)],
  dat_d3$w[which.min(dat_d3$t)]
), na.rm = TRUE)

cat("\nw0 (moyenne premiers points) :", round(w0, 4), "g\n")

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
# 6. FONCTIONS DE SIMULATION
# =============================================================================

inits_base <- c(A1 = 0, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)

# --- (A) Dose unique ---
simulate_single <- function(dose, params, obs_times, t_admin = 0) {
  ev <- eventTable()
  if (dose > 0)
    ev$add.dosing(dose = dose, nbr.doses = 1,
                  dosing.to = 1, start.time = t_admin)
  ev$add.sampling(sort(unique(c(0, obs_times))))

  out <- tryCatch(
    rxSolve(pkpd_model, as.list(params), ev, inits = inits_base),
    error = function(e) { message("rxSolve error (single): ", e$message); NULL }
  )
  if (is.null(out)) return(rep(NA_real_, length(obs_times)))
  approx(out$time, out$w, xout = obs_times, rule = 2)$y
}

# --- (B) Doses répétées toutes les 14 jours ---
simulate_multi <- function(dose, params, obs_times,
                           t_admin = 0, n_doses = 4, interval_h = 14*24) {
  ev <- eventTable()
  if (dose > 0)
    ev$add.dosing(dose = dose, nbr.doses = n_doses,
                  dosing.interval = interval_h,
                  dosing.to = 1, start.time = t_admin)
  ev$add.sampling(sort(unique(c(0, obs_times))))

  out <- tryCatch(
    rxSolve(pkpd_model, as.list(params), ev, inits = inits_base),
    error = function(e) { message("rxSolve error (multi): ", e$message); NULL }
  )
  if (is.null(out)) return(rep(NA_real_, length(obs_times)))
  approx(out$time, out$w, xout = obs_times, rule = 2)$y
}

# =============================================================================
# 7. FONCTION OBJECTIVE GÉNÉRIQUE
# =============================================================================

params_fixed <- c(PK, l0 = l0, l1 = l1, p = p)

make_objective <- function(sim_fn) {
  function(logpar) {
    par <- c(params_fixed, k1 = exp(logpar[1]), k2 = exp(logpar[2]))

    pc  <- tryCatch(sim_fn(dose0,  par, dat_ctrl$t), error = function(e) NULL)
    p3  <- tryCatch(sim_fn(dose3,  par, dat_d3$t),   error = function(e) NULL)
    p10 <- tryCatch(sim_fn(dose10, par, dat_d10$t),  error = function(e) NULL)

    bad <- function(p) is.null(p) || any(is.na(p) | p <= 0)
    if (bad(pc) || bad(p3) || bad(p10)) return(1e10)

    val <- mean((log(dat_ctrl$w) - log(pc))^2)  +
           mean((log(dat_d3$w)   - log(p3))^2)   +
           mean((log(dat_d10$w)  - log(p10))^2)

    if (!is.finite(val)) return(1e10)
    val
  }
}

# =============================================================================
# 8. OPTIMISATION — multi-start nlminb pour les deux scénarios
# =============================================================================

start_grid <- list(
  c(k1 = 0.5,  k2 = 1e-4),
  c(k1 = 0.1,  k2 = 1e-5),
  c(k1 = 1.0,  k2 = 1e-4),
  c(k1 = 0.5,  k2 = 1e-3),
  c(k1 = 0.2,  k2 = 5e-5)
)

run_optim <- function(obj_fn, label) {
  cat("\nOptimisation —", label, "...\n")
  best_obj <- Inf
  best_fit <- NULL
  for (i in seq_along(start_grid)) {
    s <- start_grid[[i]]
    fit_try <- tryCatch(
      nlminb(log(s), obj_fn,
             control = list(eval.max = 3000, iter.max = 1500,
                            rel.tol = 1e-12, x.tol = 1e-12)),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(fit_try$objective) &&
        fit_try$objective < best_obj) {
      best_obj <- fit_try$objective
      best_fit <- fit_try
    }
  }
  if (is.null(best_fit)) stop(paste("Aucune convergence pour", label))
  cat("  k1 =", round(exp(best_fit$par[1]), 6), "/h\n")
  cat("  k2 =", round(exp(best_fit$par[2]), 8), "\n")
  cat("  Objectif :", round(best_fit$objective, 5), "\n")
  best_fit
}

fit_single <- run_optim(make_objective(simulate_single), "Dose unique")
fit_multi  <- run_optim(make_objective(simulate_multi),  "Doses répétées (14j)")

k1_single <- exp(fit_single$par[1]);  k2_single <- exp(fit_single$par[2])
k1_multi  <- exp(fit_multi$par[1]);   k2_multi  <- exp(fit_multi$par[2])

cat("\n=== Résultats PKPD ===\n")
cat("Dose unique      : k1 =", round(k1_single, 6), "/h  k2 =", round(k2_single, 8), "\n")
cat("Doses répétées   : k1 =", round(k1_multi,  6), "/h  k2 =", round(k2_multi,  8), "\n")

# =============================================================================
# 9. GRAPHIQUES
# =============================================================================

times_full  <- seq(0, 49 * 24, by = 1)
groups      <- c("Contrôle", "3 mg/kg", "10 mg/kg")
lev         <- c("Contrôle", "3 mg/kg", "10 mg/kg")

make_df_sim <- function(sim_fn, params_k1k2) {
  par <- c(params_fixed, k1 = params_k1k2[1], k2 = params_k1k2[2])
  do.call(rbind, mapply(function(dose, grp)
    data.frame(t = times_full / 24,
               w = sim_fn(dose, par, times_full),
               Groupe = grp),
    list(dose0, dose3, dose10), groups, SIMPLIFY = FALSE))
}

df_obs <- rbind(
  data.frame(t = dat_ctrl$t/24, w = dat_ctrl$w, Groupe = "Contrôle"),
  data.frame(t = dat_d3$t/24,   w = dat_d3$w,   Groupe = "3 mg/kg"),
  data.frame(t = dat_d10$t/24,  w = dat_d10$w,  Groupe = "10 mg/kg")
)
df_obs$Groupe <- factor(df_obs$Groupe, levels = lev)

plot_pkpd <- function(df_sim, df_obs, title_suffix, dose_days) {
  df_sim$Groupe <- factor(df_sim$Groupe, levels = lev)
  ggplot() +
    geom_line(data  = df_sim, aes(x=t, y=w, color=Groupe), linewidth=1) +
    geom_point(data = df_obs, aes(x=t, y=w, color=Groupe), size=2.5) +
    geom_vline(xintercept = dose_days, linetype="dashed",
               color="grey60", linewidth=0.5) +
    annotate("text", x = dose_days, y = max(df_obs$w, na.rm=TRUE)*1.05,
             label = paste0("j", dose_days), size=2.8, color="grey40", hjust=0.5) +
    labs(
      title    = paste("PKPD Simeoni 2004 (rxode2) —", title_suffix),
      subtitle = paste0("Fc-silent FGFR2-huBPA-LP1 | ",
                        "PK fixé : CL=", round(PK["CL"]*24,4), " L/j/kg"),
      x = "Temps (jours)", y = "Masse tumorale (g)"
    ) +
    theme_bw(base_size = 13)
}

# Graphique A — dose unique
df_sim_single <- make_df_sim(simulate_single, c(k1_single, k2_single))
plot_pkpd(df_sim_single, df_obs, "Dose unique", dose_days = 0)
ggsave("scripts/plot_PKPD_single_FGFR2.png", width=9, height=5, dpi=150)
cat("\nGraphique dose unique → scripts/plot_PKPD_single_FGFR2.png\n")

# Graphique B — doses répétées
df_sim_multi <- make_df_sim(simulate_multi, c(k1_multi, k2_multi))
plot_pkpd(df_sim_multi, df_obs, "Doses répétées (toutes les 14j)",
          dose_days = c(0, 14, 28, 42))
ggsave("scripts/plot_PKPD_multi_FGFR2.png", width=9, height=5, dpi=150)
cat("Graphique doses répétées → scripts/plot_PKPD_multi_FGFR2.png\n")

# =============================================================================
# 10. TGI — Tumor Growth Inhibition (au jour 49)
# =============================================================================

tgi <- function(w_ctrl_end, w_treat_end, w0_val) {
  round((1 - (w_treat_end - w0_val) / (w_ctrl_end - w0_val)) * 100, 1)
}

calc_tgi <- function(sim_fn, k1, k2, label) {
  par    <- c(params_fixed, k1=k1, k2=k2)
  t_seq  <- seq(0, 49*24, by = 24)   # un point par jour
  wc  <- tail(sim_fn(dose0,  par, t_seq), 1)
  w3  <- tail(sim_fn(dose3,  par, t_seq), 1)
  w10 <- tail(sim_fn(dose10, par, t_seq), 1)
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
     params_fixed, w0, dat_ctrl, dat_d3, dat_d10,
     file = "scripts/resultats_PKPD_rxode2_FGFR2.RData")
cat("\nRésultats → scripts/resultats_PKPD_rxode2_FGFR2.RData\n")
