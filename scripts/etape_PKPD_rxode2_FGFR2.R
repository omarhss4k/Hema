# =============================================================================
# PKPD — Fc-silent B/C huBPA-LP1 (FGFR2) — Simeoni 2004 — rxode2
# -----------------------------------------------------------------------------
# PK   : 2-compartiments, paramètres FIXÉS depuis l'estimation PK
# PD   : modèle Simeoni 2004 (l0, l1, p fixés ; k1, k2 à estimer)
# Unités : heures (tout le script)
# =============================================================================

library(rxode2)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS (issus de l'estimation 2-comp rxode2)
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2

PK <- c(
  CL = unname(pk2comp_rxode2$CL),   # L/h/kg
  V1 = unname(pk2comp_rxode2$V1),   # L/kg
  V2 = unname(pk2comp_rxode2$V2),   # L/kg
  Q  = unname(pk2comp_rxode2$Q)     # L/h/kg
)

cat("=== Paramètres PK fixés ===\n")
cat("CL =", round(PK["CL"], 6), "L/h/kg  →", round(PK["CL"]*24, 4), "L/j/kg\n")
cat("V1 =", round(PK["V1"], 5), "L/kg\n")
cat("V2 =", round(PK["V2"], 5), "L/kg\n")
cat("Q  =", round(PK["Q"],  6), "L/h/kg\n")

# =============================================================================
# 2. PARAMÈTRES PD DE CROISSANCE FIXÉS (stagiaire)
#    Convertis en heures (données originales en jours)
# =============================================================================

l0 <- 0.146 / 24   # /h  (taux croissance exponentielle)
l1 <- 0.334 / 24   # g/h (taux croissance linéaire plateau)
p  <- 20            # exposant de transition (sans unité)

cat("\n=== Paramètres PD de croissance (fixés) ===\n")
cat("l0 =", round(l0, 6), "/h  (=", 0.146, "/j)\n")
cat("l1 =", round(l1, 6), "g/h (=", 0.334, "g/j)\n")
cat("p  =", p, "\n")

# =============================================================================
# 3. DONNÉES TUMORALES
# -----------------------------------------------------------------------------
# FORMAT ATTENDU pour le fichier Excel :
#   Colonne 1 : temps (jours)
#   Colonne 2 : volume tumoral groupe contrôle (g ou mm³)
#   Colonne 3 : volume tumoral dose 1 mg/kg
#   Colonne 4 : volume tumoral dose 5 mg/kg
#   Colonne 5 : volume tumoral dose 10 mg/kg
#
# Remplacer le chemin ci-dessous par le vrai fichier :
# =============================================================================

raw_pd <- read_xlsx("data/TumorVolume_FGFR2.xlsx")   # ← À RENSEIGNER

time_d  <- as.numeric(raw_pd[[1]])   # jours
w_ctrl  <- as.numeric(raw_pd[[2]])   # contrôle (véhicule)
w_d1    <- as.numeric(raw_pd[[3]])   # 1  mg/kg
w_d5    <- as.numeric(raw_pd[[4]])   # 5  mg/kg
w_d10   <- as.numeric(raw_pd[[5]])   # 10 mg/kg

# Conversion jours → heures
time_h <- time_d * 24

# Doses (µg/kg)
dose10 <- 10 * 1000
dose5  <-  5 * 1000
dose1  <-  1 * 1000
dose0  <-  0            # contrôle

# Jour d'administration → heures
t_admin_h <- 13 * 24   # ← Ajuster si nécessaire

# Nettoyage par groupe
mk <- function(t, w) {
  ok <- !is.na(w) & w > 0
  data.frame(t = t[ok], w = w[ok])
}
dat_ctrl <- mk(time_h, w_ctrl)
dat_d1   <- mk(time_h, w_d1)
dat_d5   <- mk(time_h, w_d5)
dat_d10  <- mk(time_h, w_d10)

cat("\nPoints par groupe :\n")
cat("  Contrôle :", nrow(dat_ctrl), "\n")
cat("  1 mg/kg  :", nrow(dat_d1),  "\n")
cat("  5 mg/kg  :", nrow(dat_d5),  "\n")
cat("  10 mg/kg :", nrow(dat_d10), "\n")

# =============================================================================
# 4. ESTIMATION DE w0 — masse tumorale initiale
#    Régression log-linéaire sur les points pré-traitement
# =============================================================================

est_w0 <- function(df, t_admin) {
  pre <- df[df$t < t_admin, ]
  if (nrow(pre) < 2) {
    cat("  Moins de 2 points pré-traitement — w0 = premier point observé\n")
    return(df$w[which.min(df$t)])
  }
  fit <- lm(log(w) ~ t, data = pre)
  exp(predict(fit, newdata = data.frame(t = 0)))
}

w0_ctrl <- est_w0(dat_ctrl, t_admin_h)
w0_d1   <- est_w0(dat_d1,   t_admin_h)
w0_d5   <- est_w0(dat_d5,   t_admin_h)
w0_d10  <- est_w0(dat_d10,  t_admin_h)
w0      <- mean(c(w0_ctrl, w0_d1, w0_d5, w0_d10), na.rm = TRUE)

cat("\nw0 estimés :", round(c(w0_ctrl, w0_d1, w0_d5, w0_d10), 4), "\n")
cat("w0 commun  :", round(w0, 4), "g\n")

# =============================================================================
# 5. MODÈLE PKPD rxode2 — Simeoni 2004
# -----------------------------------------------------------------------------
# États :
#   A1, A2 : quantités PK (µg/kg)
#   x1     : cellules proliférantes (g)
#   x2-x4  : compartiments de transit (cellules endommagées)
#   C1 = A1/V1 : concentration (observée)
#   w  = x1+x2+x3+x4 : masse tumorale totale (g)
# =============================================================================

pkpd_model <- rxode2({
  # --- Pharmacocinétique ---
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
  d/dt(A2) <-  (Q/V1) * A1 - (Q/V2) * A2

  # --- Pharmacodynamique (Simeoni 2004) ---
  w      <- x1 + x2 + x3 + x4
  growth <- l0 * x1 / (1 + (l0/l1 * w)^p)^(1/p)

  d/dt(x1) <- growth - k2 * C1 * x1
  d/dt(x2) <- k2 * C1 * x1 - k1 * x2
  d/dt(x3) <- k1 * (x2 - x3)
  d/dt(x4) <- k1 * (x3 - x4)
})

# =============================================================================
# 6. SIMULATION D'UN GROUPE
#    dose = 0 pour le contrôle (pas d'événement de dosage)
# =============================================================================

simulate_group <- function(dose, params, obs_times, w0_val, t_admin) {
  inits <- c(A1 = 0, A2 = 0,
             x1 = w0_val, x2 = 0, x3 = 0, x4 = 0)

  ev <- eventTable()
  if (dose > 0)
    ev$add.dosing(dose = dose, nbr.doses = 1,
                  dosing.to = 1, start.time = t_admin)
  ev$add.sampling(sort(unique(c(0, obs_times))))

  out <- tryCatch(
    rxSolve(pkpd_model, params, ev, inits = inits),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(obs_times)))

  w_tot <- out$x1 + out$x2 + out$x3 + out$x4
  approx(out$time, w_tot, xout = obs_times, rule = 2)$y
}

# =============================================================================
# 7. FONCTION OBJECTIVE — résidus log, poids égaux par groupe
#    Paramètres libres : k1, k2  (log-transformés)
#    Paramètres fixés  : CL, V1, V2, Q, l0, l1, p, w0
# =============================================================================

params_fixed <- c(PK, l0 = l0, l1 = l1, p = p)

objective <- function(logpar) {
  k1 <- exp(logpar[1])
  k2 <- exp(logpar[2])
  par <- c(params_fixed, k1 = k1, k2 = k2)

  pc  <- tryCatch(simulate_group(dose0,  par, dat_ctrl$t, w0, t_admin_h),
                  error = function(e) NULL)
  p1  <- tryCatch(simulate_group(dose1,  par, dat_d1$t,   w0, t_admin_h),
                  error = function(e) NULL)
  p5  <- tryCatch(simulate_group(dose5,  par, dat_d5$t,   w0, t_admin_h),
                  error = function(e) NULL)
  p10 <- tryCatch(simulate_group(dose10, par, dat_d10$t,  w0, t_admin_h),
                  error = function(e) NULL)

  bad <- function(p, d) is.null(p) || any(is.na(p) | p <= 0)
  if (bad(pc,  dat_ctrl)) return(1e10)
  if (bad(p1,  dat_d1))   return(1e10)
  if (bad(p5,  dat_d5))   return(1e10)
  if (bad(p10, dat_d10))  return(1e10)

  mean((log(dat_ctrl$w) - log(pc))^2)  +
  mean((log(dat_d1$w)   - log(p1))^2)  +
  mean((log(dat_d5$w)   - log(p5))^2)  +
  mean((log(dat_d10$w)  - log(p10))^2)
}

# =============================================================================
# 8. OPTIMISATION — multi-start nlminb (k1, k2)
#    Plages biologiquement raisonnables :
#    k1 : 0.01–5 /h  (transit ≈ taux de mort cellulaire)
#    k2 : 1e-7–0.01  (sensibilité au médicament)
# =============================================================================

start_grid <- list(
  c(k1 = 0.5,  k2 = 1e-4),
  c(k1 = 0.1,  k2 = 1e-5),
  c(k1 = 1.0,  k2 = 1e-4),
  c(k1 = 0.5,  k2 = 1e-3),
  c(k1 = 0.2,  k2 = 5e-5)
)

cat("\nOptimisation nlminb (multi-start) ...\n")
best_obj <- Inf
best_fit <- NULL

for (i in seq_along(start_grid)) {
  s <- start_grid[[i]]
  cat("  Départ", i, ": k1 =", s["k1"], " k2 =", s["k2"], "...\n")

  fit_try <- tryCatch(
    nlminb(log(s), objective,
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

if (is.null(best_fit)) stop("Aucun point de départ n'a convergé.")

k1_est <- exp(best_fit$par[1])
k2_est <- exp(best_fit$par[2])

cat("\n=== Paramètres PKPD estimés ===\n")
cat("k1 =", round(k1_est, 6), "/h\n")
cat("k2 =", round(k2_est, 8), "\n")
cat("Objectif :", round(best_fit$objective, 5), "\n")

# =============================================================================
# 9. SIMULATION FINALE ET GRAPHIQUE
# =============================================================================

best_params <- c(params_fixed, k1 = k1_est, k2 = k2_est)
times_full  <- seq(0, max(time_h) * 1.05, by = 1)

sim <- list(
  ctrl = simulate_group(dose0,  best_params, times_full, w0, t_admin_h),
  d1   = simulate_group(dose1,  best_params, times_full, w0, t_admin_h),
  d5   = simulate_group(dose5,  best_params, times_full, w0, t_admin_h),
  d10  = simulate_group(dose10, best_params, times_full, w0, t_admin_h)
)

groups <- c("Contrôle", "1 mg/kg", "5 mg/kg", "10 mg/kg")

df_sim <- do.call(rbind, mapply(function(s, g)
  data.frame(t = times_full / 24, w = s, Groupe = g),
  sim, groups, SIMPLIFY = FALSE))

df_obs <- rbind(
  data.frame(t = dat_ctrl$t / 24, w = dat_ctrl$w, Groupe = "Contrôle"),
  data.frame(t = dat_d1$t   / 24, w = dat_d1$w,   Groupe = "1 mg/kg"),
  data.frame(t = dat_d5$t   / 24, w = dat_d5$w,   Groupe = "5 mg/kg"),
  data.frame(t = dat_d10$t  / 24, w = dat_d10$w,  Groupe = "10 mg/kg")
)

lev <- c("Contrôle", "1 mg/kg", "5 mg/kg", "10 mg/kg")
df_sim$Groupe <- factor(df_sim$Groupe, levels = lev)
df_obs$Groupe <- factor(df_obs$Groupe, levels = lev)

ggplot() +
  geom_line(data  = df_sim, aes(x = t, y = w, color = Groupe), linewidth = 1) +
  geom_point(data = df_obs, aes(x = t, y = w, color = Groupe), size = 2.5) +
  geom_vline(xintercept = t_admin_h / 24, linetype = "dashed", color = "grey50") +
  annotate("text", x = t_admin_h/24 + 0.3,
           y = max(df_obs$w, na.rm = TRUE) * 0.95,
           label = "Administration", hjust = 0, size = 3.5, color = "grey40") +
  labs(
    title    = "PKPD Simeoni 2004 (rxode2) — Fc-silent B/C huBPA-LP1 (FGFR2)",
    subtitle = paste0(
      "k1 = ", round(k1_est, 5), " /h",
      "   k2 = ", round(k2_est, 7),
      "   (PK fixé : CL=", round(PK["CL"]*24, 4),
      " L/j/kg, V1=", round(PK["V1"], 4), " L/kg)"
    ),
    x = "Temps (jours)",
    y = "Masse tumorale (g)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PKPD_rxode2_FGFR2.png", width = 8, height = 5, dpi = 150)
cat("\nGraphique → scripts/plot_PKPD_rxode2_FGFR2.png\n")

# =============================================================================
# 10. SAUVEGARDE
# =============================================================================

save(k1_est, k2_est, best_params, w0,
     times_full, sim, dat_ctrl, dat_d1, dat_d5, dat_d10,
     file = "scripts/resultats_PKPD_rxode2_FGFR2.RData")
cat("Résultats → scripts/resultats_PKPD_rxode2_FGFR2.RData\n")
