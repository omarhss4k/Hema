# =============================================================================
# PKPD — Fc-silent FGFR2-huBPA-LP1 — Simeoni 2004 — Protocole GA6822
# Q2W × 3 doses (j0, j14, j28) | Simulation à 56 jours
# Identique à etape_PKPD_rxode2_FGFR2.R pour le fitting (données j0–j49) ;
# seules les projections graphiques et TGI sont à 56j / 3 doses.
# =============================================================================

library(deSolve)
library(DEoptim)
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
# 2. PARAMÈTRE PD FIXÉ — p (Hill) ; l0 et l1 estimés (section 5.5)
# =============================================================================

p <- 20   # coefficient de Hill — fixé (Simeoni 2004)

C1_max_10 <- 10000 / PK["V1"]
cat("\nRepère PK : C1_max(10 mg/kg) =", round(C1_max_10, 0), "µg/L\n")

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
dat_ctrl_all <- mk(time_h, w_ctrl)          # tous points contrôle (graphiques)
dat_ctrl     <- dat_ctrl_all[dat_ctrl_all$t <= 35 * 24, ]  # j≤35 seulement (censure j>35)
dat_d10      <- mk(time_h, w_d10)
dat_d3       <- mk(time_h, w_d3)

cat("\nPoints par groupe :\n")
cat("  Contrôle  :", nrow(dat_ctrl), "pts optimisation |", nrow(dat_ctrl_all), "pts total\n")
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
# 5. MODÈLE PKPD — Simeoni 2004 — ODE deSolve
# =============================================================================

pkpd_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    C1     <- A1 / V1
    w      <- x1 + x2 + x3 + x4
    growth <- l0 * x1 / (1 + (l0/l1 * w)^p)^(1/p)

    dA1 <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
    dA2 <-  (Q/V1) * A1 - (Q/V2) * A2
    dx1 <- growth - k2 * C1 * x1
    dx2 <- k2 * C1 * x1 - k1 * x2
    dx3 <- k1 * (x2 - x3)
    dx4 <- k1 * (x3 - x4)

    list(c(dA1, dA2, dx1, dx2, dx3, dx4))
  })
}

# =============================================================================
# 5.5 PRÉ-CALIBRATION l0, l1 DEPUIS LE GROUPE CONTRÔLE
#     dose = 0  →  C1 = 0  →  croissance pure  →  l0, l1 identifiables seuls
# =============================================================================

obj_growth <- function(logpar) {
  l0_t <- unname(exp(logpar[1]))
  l1_t <- unname(exp(logpar[2]))
  pars   <- c(as.list(PK), l0 = l0_t, l1 = l1_t, p = p, k1 = 1, k2 = 0)
  state0 <- c(A1 = 0, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)
  out <- tryCatch(
    as.data.frame(lsoda(state0, sort(unique(c(0, dat_ctrl$t))), pkpd_ode, pars,
                        rtol = 1e-6, atol = 1e-8)),
    error = function(e) NULL
  )
  if (is.null(out)) return(1e10)
  w_vec <- with(out, x1 + x2 + x3 + x4)
  if (any(!is.finite(w_vec))) return(1e10)
  pc  <- approx(out$time, pmax(w_vec, 1e-9), xout = dat_ctrl$t, rule = 2)$y
  val <- mean((log(dat_ctrl$w) - log(pc))^2)
  if (!is.finite(val)) 1e10 else val
}

cat("\nPré-calibration l0, l1 depuis le groupe contrôle (DEoptim)...\n")
fit_growth <- DEoptim(
  obj_growth,
  lower   = c(log(0.001/24), log(0.005/24)),
  upper   = c(log(0.5/24),   log(5/24)),
  control = DEoptim.control(NP = 20L, itermax = 500L,
                            F = 0.8, CR = 0.9, strategy = 2L,
                            trace = FALSE)
)

l0 <- unname(exp(fit_growth$optim$bestmem[1]))
l1 <- unname(exp(fit_growth$optim$bestmem[2]))

cat("  l0 =", round(l0 * 24, 4), "/j  (", round(l0, 7), "/h)\n")
cat("  l1 =", round(l1 * 24, 4), "g/j (", round(l1, 6), "g/h)\n")
cat("  Objectif contrôle :", round(fit_growth$optim$bestval, 5), "\n")

k2_ref <- l0 / C1_max_10
cat("  Repère k2 : k2_ref ≈", formatC(k2_ref, format = "e", digits = 2), "\n")

# =============================================================================
# 6. FONCTIONS DE SIMULATION — deSolve/lsoda + pénalité si échec
# =============================================================================

.PENALTY <- 1e6

.sim <- function(state0, pars_list, obs_times, events_df = NULL) {
  times <- sort(unique(c(0, obs_times)))
  out <- tryCatch(
    as.data.frame(lsoda(state0, times, pkpd_ode, pars_list,
                        events = if (!is.null(events_df)) list(data = events_df) else NULL,
                        rtol = 1e-6, atol = 1e-8)),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(.PENALTY, length(obs_times)))
  w_vec <- with(out, x1 + x2 + x3 + x4)
  if (any(!is.finite(w_vec)) || max(out$time) < max(obs_times) - 1)
    return(rep(.PENALTY, length(obs_times)))
  approx(out$time, pmax(w_vec, 1e-9), xout = obs_times, rule = 2)$y
}

simulate_single <- function(dose, params, obs_times, t_admin = 0) {
  state0    <- c(A1 = dose, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)
  pars_list <- as.list(params)
  .sim(state0, pars_list, obs_times, events_df = NULL)
}

# GA6822 : Q2W × 3 doses (j0, j14, j28) — n_doses = 3 par défaut
simulate_multi <- function(dose, params, obs_times,
                           t_admin = 0, n_doses = 3, interval_h = 14*24) {
  state0    <- c(A1 = dose, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)
  pars_list <- as.list(params)
  events_df <- if (n_doses > 1)
    data.frame(var    = "A1",
               time   = seq(t_admin + interval_h,
                            t_admin + (n_doses - 1) * interval_h,
                            by = interval_h),
               value  = dose,
               method = "add")
  else NULL
  .sim(state0, pars_list, obs_times, events_df)
}

# =============================================================================
# 7. FONCTION OBJECTIVE GÉNÉRIQUE
# =============================================================================

params_fixed <- c(PK, l0 = l0, l1 = l1, p = p)

STAB_FACTOR <- 50
.stable_k2 <- function(k2) k2 * C1_max_10 <= STAB_FACTOR * l0

.residuals <- function(sim_fn, par) {
  pc  <- tryCatch(sim_fn(dose0,  par, dat_ctrl$t), error = function(e) NULL)
  p3  <- tryCatch(sim_fn(dose3,  par, dat_d3$t),   error = function(e) NULL)
  p10 <- tryCatch(sim_fn(dose10, par, dat_d10$t),  error = function(e) NULL)

  bad <- function(p) is.null(p) || length(p) == 0 || any(!is.finite(p) | p <= 0)
  if (bad(pc) || bad(p3) || bad(p10)) return(1e10)

  val <- mean((log(dat_ctrl$w) - log(pc))^2)  +
         mean((log(dat_d3$w)   - log(p3))^2)  +
         mean((log(dat_d10$w)  - log(p10))^2)
  if (!is.finite(val)) 1e10 else val
}

make_objective_k2 <- function(sim_fn, k1_val) {
  function(logpar_k2) {
    k2 <- unname(exp(logpar_k2[1]))
    if (!.stable_k2(k2)) return(1e10)
    par <- c(params_fixed, k1 = k1_val, k2 = k2)
    .residuals(sim_fn, par)
  }
}

make_objective_k2_global <- function(k1_val) {
  function(logpar_k2) {
    k2 <- unname(exp(logpar_k2[1]))
    if (!.stable_k2(k2)) return(1e10)
    par <- c(params_fixed, k1 = k1_val, k2 = k2)
    rs  <- .residuals(simulate_single, par)
    rm  <- .residuals(simulate_multi,  par)
    if (rs >= 1e10 || rm >= 1e10) return(1e10)
    (rs + rm) / 2
  }
}

# =============================================================================
# 8. OPTIMISATION — DEoptim
# =============================================================================

LOG_LOWER_1 <- c(log(1e-9))
LOG_UPPER_1 <- c(log(0.1))

run_optim <- function(obj_fn, lower, upper, label, NP = NULL) {
  np  <- length(lower)
  NP  <- if (is.null(NP)) 10L * np else as.integer(NP)
  cat("\nDEoptim —", label, " (NP =", NP, ", itermax = 1000)...\n")
  out <- tryCatch(
    DEoptim(obj_fn, lower = lower, upper = upper,
            control = DEoptim.control(NP       = NP,
                                      itermax  = 1000L,
                                      F        = 0.8,
                                      CR       = 0.9,
                                      strategy = 2L,
                                      trace    = FALSE)),
    error = function(e) { message("DEoptim erreur : ", e$message); NULL }
  )
  if (is.null(out)) stop(paste("DEoptim échoué pour", label))
  out
}

# =============================================================================
# 9. PHASE 1 — k2 GLOBAL (k1 fixé)
# =============================================================================

k1_fixed <- 0.5 / 24   # /h (= 0.5 /j, transit moyen = 8 j)

fit_global <- run_optim(make_objective_k2_global(k1_fixed),
                        LOG_LOWER_1, LOG_UPPER_1, "Fit global k2 (k1 fixé)", NP = 10L)

k2_global <- unname(exp(fit_global$optim$bestmem[1]))

cat("\n=== Fit global k2 (k1 =", round(k1_fixed * 24, 4), "/j fixé) ===\n")
cat("  k2 =", formatC(k2_global, format="e", digits=3), "\n")
cat("  Objectif :", round(fit_global$optim$bestval, 5), "\n")

# =============================================================================
# 10. PHASE 2 — k2 LIBRE PAR SCÉNARIO
# =============================================================================

fit_k2_single <- run_optim(make_objective_k2(simulate_single, k1_fixed),
                           LOG_LOWER_1, LOG_UPPER_1,
                           "k2 dose unique (k1 fixé)", NP = 10L)
fit_k2_multi  <- run_optim(make_objective_k2(simulate_multi,  k1_fixed),
                           LOG_LOWER_1, LOG_UPPER_1,
                           "k2 Q2W×3 (k1 fixé)", NP = 10L)

k2_single <- unname(exp(fit_k2_single$optim$bestmem[1]))
k2_multi  <- unname(exp(fit_k2_multi$optim$bestmem[1]))

cat("\n=== Fits séparés k2 (k1 =", round(k1_fixed * 24, 4), "/j fixé) ===\n")
cat("  k2 dose unique :", formatC(k2_single, format="e", digits=3), "\n")
cat("  k2 Q2W×3       :", formatC(k2_multi,  format="e", digits=3), "\n")
ratio <- k2_multi / k2_single
cat("  Ratio k2_multi/k2_single :", round(ratio, 3),
    if (ratio > 1.5) "→ accumulation possible" else
    if (ratio < 0.67) "→ résistance possible"  else
    "→ pas de dérive notable", "\n")

# =============================================================================
# 11. GRAPHIQUES — simulation à 56 jours, doses Q2W (j0, j14, j28)
# =============================================================================

T_END     <- 56 * 24              # horizon de simulation (h)
DOSE_DAYS <- c(0, 14, 28)         # protocole GA6822 Q2W × 3

times_full <- seq(0, T_END, by = 1)
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

# Observations (tous points pour les graphiques, dont contrôle j>35)
df_obs <- rbind(
  data.frame(t = dat_ctrl_all$t/24, w = dat_ctrl_all$w, Groupe = "Contrôle"),
  data.frame(t = dat_d3$t/24,       w = dat_d3$w,       Groupe = "3 mg/kg"),
  data.frame(t = dat_d10$t/24,      w = dat_d10$w,      Groupe = "10 mg/kg")
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
      title    = paste("PKPD Simeoni 2004 — GA6822 Q2W×3 —", title_suffix),
      subtitle = paste0("Fc-silent FGFR2-huBPA-LP1 | ",
                        "k1=", round(k1 * 24, 3), " /j  ",
                        "k2=", formatC(k2, format="e", digits=2)),
      x = "Temps (jours)", y = "Volume tumoral (g)"
    ) +
    theme_bw(base_size = 13)
}

# Fit global (k2 commun)
df_sim_global_m <- make_df_sim(simulate_multi, k1_fixed, k2_global)
plot_pkpd(df_sim_global_m, "Doses répétées — fit global",
          dose_days = DOSE_DAYS, k1_fixed, k2_global)
ggsave("scripts/plot_GA6822_global_multi_FGFR2.png", width=9, height=5, dpi=150)

df_sim_global_s <- make_df_sim(simulate_single, k1_fixed, k2_global)
plot_pkpd(df_sim_global_s, "Dose unique — fit global",
          dose_days = 0, k1_fixed, k2_global)
ggsave("scripts/plot_GA6822_global_single_FGFR2.png", width=9, height=5, dpi=150)

# k2 libre par scénario
df_sim_k2m <- make_df_sim(simulate_multi, k1_fixed, k2_multi)
plot_pkpd(df_sim_k2m, "Q2W×3 — k2 libre",
          dose_days = DOSE_DAYS, k1_fixed, k2_multi)
ggsave("scripts/plot_GA6822_k2sep_multi_FGFR2.png", width=9, height=5, dpi=150)

df_sim_k2s <- make_df_sim(simulate_single, k1_fixed, k2_single)
plot_pkpd(df_sim_k2s, "Dose unique — k2 libre",
          dose_days = 0, k1_fixed, k2_single)
ggsave("scripts/plot_GA6822_k2sep_single_FGFR2.png", width=9, height=5, dpi=150)

cat("\nGraphiques → scripts/plot_GA6822_*_FGFR2.png\n")

# =============================================================================
# 12. TGI — calculé à j56
# =============================================================================

tgi <- function(w_ctrl_end, w_treat_end, w0_val) {
  delta_ctrl <- w_ctrl_end - w0_val
  if (is.na(delta_ctrl) || delta_ctrl <= 0) return(NA_real_)
  round((1 - (w_treat_end - w0_val) / delta_ctrl) * 100, 1)
}

calc_tgi <- function(sim_fn, k1, k2, label) {
  par   <- c(params_fixed, k1 = k1, k2 = k2)
  t_seq <- seq(0, T_END, by = 24)
  wc  <- tail(sim_fn(dose0,  par, t_seq), 1)
  w3  <- tail(sim_fn(dose3,  par, t_seq), 1)
  w10 <- tail(sim_fn(dose10, par, t_seq), 1)

  tgi3  <- tgi(wc, w3,  w0)
  tgi10 <- tgi(wc, w10, w0)

  cat("\nTGI à j56 —", label, ":\n")
  cat("   3 mg/kg  :", if (is.na(tgi3))  "NA" else
      if (tgi3 > 100) paste0(tgi3, "% (régression)") else paste0(tgi3, "%"), "\n")
  cat("  10 mg/kg  :", if (is.na(tgi10)) "NA" else
      if (tgi10 > 100) paste0(tgi10, "% (régression)") else paste0(tgi10, "%"), "\n")
}

cat("\n--- TGI j56 — fit global ---\n")
calc_tgi(simulate_single, k1_fixed, k2_global, "Dose unique")
calc_tgi(simulate_multi,  k1_fixed, k2_global, "Q2W×3")

cat("\n--- TGI j56 — k2 libre ---\n")
calc_tgi(simulate_single, k1_fixed, k2_single, "Dose unique")
calc_tgi(simulate_multi,  k1_fixed, k2_multi,  "Q2W×3")

# =============================================================================
# 13. SAUVEGARDE
# =============================================================================

save(k1_fixed, k2_global, k2_single, k2_multi,
     params_fixed, w0, dat_ctrl, dat_d3, dat_d10,
     file = "scripts/resultats_GA6822_FGFR2.RData")
cat("\nRésultats → scripts/resultats_GA6822_FGFR2.RData\n")
