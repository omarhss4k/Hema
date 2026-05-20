# =============================================================================
# Modèle PK/PD TGI — Simeoni (2004) — Fc-silent FGFR2-huBPA-LP1
#
# PK  : 2 compartiments IV bolus (paramètres fixés depuis fit précédent)
# PD  : Modèle de Simeoni — compartiments de transit
#
#   dx1/dt = [g(w)/w - k2·C] · x1           ← cellules proliférantes
#   dx2/dt =  k2·C · x1 - k1·x2             ← cellules endommagées (transit 1)
#   dx3/dt =  k1·x2 - k1·x3                 ← transit 2
#   dx4/dt =  k1·x3 - k1·x4                 ← transit 3 → mort
#   w      =  x1 + x2 + x3 + x4             ← volume tumoral total
#
#   g(w)   = λ0·w / (1 + (λ0·w / λ1)^ψ)^(1/ψ)   ψ = 20 (fixé)
#              ↳ exponentiel si w << λ1/λ0
#              ↳ linéaire (≈ λ1) si w >> λ1/λ0
#
#   TSC    ≈ λ0 / k2   (Tumor Static Concentration — phase exponentielle)
#   MTT    = 4 / k1    (Mean Transit Time, 3 compartiments de transit + x1)
#
# Paramètres PD estimés : λ0, λ1, k1, k2  (log-espace)
# Unités : temps en jours | dose en µg/kg | concentration en µg/L
#
# Schéma : Q2W × 4 (j0, j14, j28, j42)
# Groupes : Contrôle (Group 01) | 3 mg/kg (Group 04) | 10 mg/kg (Group 03)
# =============================================================================

library(deSolve)
library(DEoptim)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS (conversion h → j)
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2

pk_fixed <- c(
  CL = unname(pk2comp_rxode2$CL) * 24,
  V1 = unname(pk2comp_rxode2$V1),
  V2 = unname(pk2comp_rxode2$V2),
  Q  = unname(pk2comp_rxode2$Q)  * 24
)

cat("=== Paramètres PK fixés (jours) ===\n")
cat(sprintf("CL = %.5f L/j/kg\n", pk_fixed["CL"]))
cat(sprintf("V1 = %.5f L/kg\n",   pk_fixed["V1"]))
cat(sprintf("V2 = %.5f L/kg\n",   pk_fixed["V2"]))
cat(sprintf("Q  = %.5f L/j/kg\n", pk_fixed["Q"]))

# =============================================================================
# 2. DONNÉES TUMORALES
# =============================================================================

raw_all <- suppressMessages(
  read_xlsx("TumorVolume_FGFR2.xlsx",
            col_names = FALSE, .name_repair = "minimal")
)
raw_all    <- raw_all[, colSums(!is.na(raw_all)) > 0]
header_idx <- which(raw_all[[1]] == "Group")

times_all <- suppressWarnings(
  as.numeric(as.character(unlist(raw_all[header_idx[1], -1])))
)
valid   <- !is.na(times_all)
times_d <- times_all[valid]

extract_block <- function(start, end) {
  blk <- raw_all[start:end, ]
  blk <- blk[!is.na(blk[[1]]) & blk[[1]] != "", ]
  dat <- blk[, c(TRUE, valid)]
  for (j in 2:ncol(dat)) storage.mode(dat[[j]]) <- "numeric"
  dat
}

blk_mean <- extract_block(header_idx[1] + 1, header_idx[2] - 1)
blk_sem  <- extract_block(header_idx[2] + 1, nrow(raw_all))

get_row <- function(blk, pat)
  as.numeric(unlist(blk[grep(pat, blk[[1]], ignore.case = TRUE)[1], -1]))

tv_ctrl  <- get_row(blk_mean, "Group 01")
tv_iso   <- get_row(blk_mean, "Group 02")
tv_d3    <- get_row(blk_mean, "Group 04")
tv_d10   <- get_row(blk_mean, "Group 03")

sem_ctrl <- pmax(get_row(blk_sem, "Group 01"), 1)
sem_iso  <- pmax(get_row(blk_sem, "Group 02"), 1)
sem_d3   <- pmax(get_row(blk_sem, "Group 04"), 1)
sem_d10  <- pmax(get_row(blk_sem, "Group 03"), 1)

get_tv0 <- function(tv) {
  v <- tv[times_d == 0]
  if (length(v) == 0 || is.na(v[1])) tv[!is.na(tv)][1] else v[1]
}
tv0_ctrl <- get_tv0(tv_ctrl)
tv0_iso  <- get_tv0(tv_iso)
tv0_d3   <- get_tv0(tv_d3)
tv0_d10  <- get_tv0(tv_d10)

ok_ctrl <- !is.na(tv_ctrl)
ok_iso  <- !is.na(tv_iso)
ok_d3   <- !is.na(tv_d3)
ok_d10  <- !is.na(tv_d10)

# Troncature du contrôle à j≤39 : après j39 les animaux atteignent l'endpoint
# humanitaire et sont sacrifiés → les TV post-j39 (≈142 mm³) sont des artefacts
# de censure qui biaiseraient λ0 vers le bas.
# ok_ctrl_fit : masque restreint pour l'estimation
# ok_ctrl     : masque complet conservé pour l'affichage
CTRL_CENSOR_DAY <- 39
ok_ctrl_fit <- ok_ctrl & times_d <= CTRL_CENSOR_DAY

cat(sprintf("\nTV initiale (j0) : Ctrl=%.0f | Isotype=%.0f | 3mg=%.0f | 10mg=%.0f mm³\n",
            tv0_ctrl, tv0_iso, tv0_d3, tv0_d10))
cat(sprintf("Contrôle : %d points valides, %d retenus pour le fit (≤ j%d)\n",
            sum(ok_ctrl), sum(ok_ctrl_fit), CTRL_CENSOR_DAY))

# =============================================================================
# 3. ÉQUATIONS DU MODÈLE (deSolve)
# =============================================================================

PSI <- 20   # exposant de la fonction de croissance (fixé, Simeoni 2004)

simeoni_rhs <- function(t, state, parms) {
  # Clamp à 0 : empêche x1<0 de s'auto-amplifier (instabilité numérique)
  A1 <- max(state["A1"], 0); A2 <- max(state["A2"], 0)
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)

  CL <- parms["CL"]; V1 <- parms["V1"]
  V2 <- parms["V2"]; Q  <- parms["Q"]
  L0 <- parms["L0"]; L1 <- parms["L1"]
  k1 <- parms["k1"]; k2 <- parms["k2"]

  # PK
  C1  <- A1 / V1
  dA1 <- -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
  dA2 <-  (Q/V1)*A1 - (Q/V2)*A2

  # PD — fonction de croissance de Simeoni
  w   <- x1 + x2 + x3 + x4
  gw  <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)

  growth_rate <- if (w > 1e-12) gw / w else L0

  dx1 <- (growth_rate - k2 * C1) * x1
  dx2 <- k2 * C1 * x1 - k1 * x2
  dx3 <- k1 * x2 - k1 * x3
  dx4 <- k1 * x3 - k1 * x4

  list(c(A1 = dA1, A2 = dA2,
         x1 = dx1, x2 = dx2, x3 = dx3, x4 = dx4))
}

# Contrôle : même modèle sans drogue (A1=A2=0 fixes, k2 ignoré)
simeoni_ctrl_rhs <- function(t, state, parms) {
  x1 <- max(state["x1"], 0)
  x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0)
  x4 <- max(state["x4"], 0)
  L0 <- parms["L0"]; L1 <- parms["L1"]

  w  <- x1 + x2 + x3 + x4
  gw <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  dx1 <- growth_rate * x1
  dx2 <- 0; dx3 <- 0; dx4 <- 0

  list(c(x1 = dx1, x2 = dx2, x3 = dx3, x4 = dx4))
}

# =============================================================================
# 4. FONCTIONS DE SIMULATION
# =============================================================================

dose_days <- c(0, 14, 28, 42)

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

# Redémarre l'intégrateur à chaque dose : plus robuste que lsoda events
sim_treated <- function(dose_ugkg, tv0, params, times_out) {
  t_end  <- max(times_out)
  breaks <- c(dose_days, t_end + 1)   # limites de chaque intervalle dose

  state <- c(A1 = dose_ugkg, A2 = 0,
             x1 = tv0,       x2 = 0, x3 = 0, x4 = 0)

  t_all <- numeric(0)
  w_all <- numeric(0)

  for (i in seq_along(dose_days)) {
    t_start <- dose_days[i]
    t_stop  <- min(breaks[i + 1], t_end)
    if (t_start >= t_end) break

    # Bolus additif à chaque dose (sauf la 1re déjà dans les CI)
    if (i > 1) state["A1"] <- state["A1"] + dose_ugkg

    t_seg <- sort(unique(c(t_start,
                           times_out[times_out > t_start & times_out <= t_stop],
                           t_stop)))
    if (length(t_seg) < 2) t_seg <- c(t_start, t_stop)

    seg <- tryCatch(
      as.data.frame(lsoda(
        y     = state,
        times = t_seg,
        func  = simeoni_rhs,
        parms = params,
        atol  = 1e-6, rtol = 1e-6
      )),
      error = function(e) NULL
    )
    if (is.null(seg) || any(is.na(seg$x1))) return(rep(NA_real_, length(times_out)))

    w_seg <- pmax(seg$x1, 0) + pmax(seg$x2, 0) +
             pmax(seg$x3, 0) + pmax(seg$x4, 0)

    # Évite les doublons au raccord entre segments
    keep   <- if (length(t_all) > 0) seg$time > tail(t_all, 1) else rep(TRUE, nrow(seg))
    t_all  <- c(t_all, seg$time[keep])
    w_all  <- c(w_all, w_seg[keep])

    # État final → conditions initiales du segment suivant
    last  <- seg[nrow(seg), ]
    state <- c(A1 = last$A1, A2 = last$A2,
               x1 = last$x1, x2 = last$x2, x3 = last$x3, x4 = last$x4)
  }

  if (length(t_all) == 0) return(rep(NA_real_, length(times_out)))
  approx(t_all, w_all, xout = times_out, rule = 2)$y
}

# =============================================================================
# 5. FONCTION OBJECTIVE
#    Paramètres libres : L0, L1, k1, k2 (log-espace)
#    Pondération par 1/SEM² sur log(TV)
# =============================================================================

objective_simeoni <- function(logpar) {
  par <- exp(logpar)
  # L0 et L1 fixés hors optimisation — seuls k1 et k2 sont libres
  L0 <- L0_fixed; L1 <- L1_fixed
  k1 <- unname(par[1]); k2 <- unname(par[2])
  if (any(par <= 0)) return(1e12)

  params_all <- c(pk_fixed, L0 = L0, L1 = L1, k1 = k1, k2 = k2)

  # — Contrôle (tronqué à j≤CTRL_CENSOR_DAY) —
  t_c       <- times_d[ok_ctrl_fit]
  pred_ctrl <- sim_ctrl_fn(L0, L1, tv0_ctrl, t_c)
  obs_ctrl  <- tv_ctrl[ok_ctrl_fit]
  w_ctrl    <- 1 / sem_ctrl[ok_ctrl_fit]^2
  if (any(is.na(pred_ctrl))) return(1e12)
  pred_ctrl <- pmax(pred_ctrl, 0.1)

  # — 3 mg/kg —
  t_3     <- times_d[ok_d3]
  pred_d3 <- sim_treated(3000, tv0_d3, params_all, t_3)
  obs_d3  <- tv_d3[ok_d3]
  w_d3    <- 1 / sem_d3[ok_d3]^2
  if (any(is.na(pred_d3))) return(1e12)
  pred_d3 <- pmax(pred_d3, 0.1)

  # — 10 mg/kg —
  t_10     <- times_d[ok_d10]
  pred_d10 <- sim_treated(10000, tv0_d10, params_all, t_10)
  obs_d10  <- tv_d10[ok_d10]
  w_d10    <- 1 / sem_d10[ok_d10]^2
  if (any(is.na(pred_d10))) return(1e12)
  pred_d10 <- pmax(pred_d10, 0.1)

  sum(w_ctrl * (log(obs_ctrl) - log(pred_ctrl))^2, na.rm = TRUE) +
  sum(w_d3   * (log(obs_d3)   - log(pred_d3))^2,   na.rm = TRUE) +
  sum(w_d10  * (log(obs_d10)  - log(pred_d10))^2,  na.rm = TRUE)
}

# =============================================================================
# 6. VALEURS INITIALES
#
# Stratégie : λ0 ET L1 fixés depuis le contrôle. Seuls k1 et k2 sont libres.
#
# λ0 : régression log-linéaire sur la phase exponentielle (j0–j21)
# L1 : estimé depuis le plateau du contrôle (j21–j39)
#      Dans Simeoni, la transition expo→linéaire se produit à w ≈ L1/λ0
#      → L1 = λ0 × TV_plateau
# =============================================================================

# λ0 FIXÉ : phase exponentielle j0–j21
L0_EXPO_DAY  <- 21
ok_ctrl_expo <- ok_ctrl & times_d <= L0_EXPO_DAY
lm_ctrl      <- lm(log(tv_ctrl[ok_ctrl_expo]) ~ times_d[ok_ctrl_expo])
L0_fixed     <- max(coef(lm_ctrl)[2], 0.005)

# L1 FIXÉ : médiane du contrôle sur la phase de plateau (j21–j39)
ok_ctrl_plateau <- ok_ctrl & times_d >= L0_EXPO_DAY & times_d <= CTRL_CENSOR_DAY
TV_plateau      <- median(tv_ctrl[ok_ctrl_plateau], na.rm = TRUE)
L1_fixed        <- L0_fixed * TV_plateau

# k1 : valeur initiale MTT ~ 7 jours
k1_init <- 4 / 7

# k2 : TSC ≈ Cmax(10 mg/kg)
Cmax_10mgkg <- 10000 / as.numeric(pk_fixed["V1"])
k2_init     <- L0_fixed / Cmax_10mgkg

init_pd <- c(k1 = k1_init, k2 = k2_init)

cat("\n=== Paramètres fixés depuis le contrôle ===\n")
cat(sprintf("λ0 (L0_fixed) = %.5f /j     (t½ = %.1f j — régression j0–j%d)\n",
            L0_fixed, log(2)/L0_fixed, L0_EXPO_DAY))
cat(sprintf("λ1 (L1_fixed) = %.2f mm³/j  (= λ0 × TV_plateau, TV_plateau = %.0f mm³)\n",
            L1_fixed, TV_plateau))
cat("\n=== Valeurs initiales (paramètres libres : k1, k2) ===\n")
cat(sprintf("k1 = %.4f /j   (transit rate, MTT = %.1f j)\n", k1_init, 4/k1_init))
cat(sprintf("k2 = %.2e L/µg/j\n", k2_init))
cat(sprintf("TSC initiale ≈ %.1f µg/L  (= L0/k2)\n", L0_fixed / k2_init))

# =============================================================================
# 6b. TEST SIMULATION AUX VALEURS INITIALES (diagnostic)
# =============================================================================

params_test <- c(pk_fixed,
                 L0 = unname(L0_fixed), L1 = unname(L1_fixed),
                 k1 = unname(k1_init),  k2 = unname(k2_init))

cat("\n=== TEST simulation aux valeurs initiales ===\n")
tc_test <- sim_ctrl_fn(L0_fixed, L1_fixed, tv0_ctrl, times_d[ok_ctrl_fit])
cat(sprintf("Contrôle : %s\n",
    if (any(is.na(tc_test))) "ECHEC (NA)" else
    paste(round(head(tc_test, 4)), collapse = " | ")))

t3_test <- sim_treated(3000, tv0_d3, params_test, times_d[ok_d3])
cat(sprintf("3 mg/kg  : %s\n",
    if (any(is.na(t3_test))) "ECHEC (NA)" else
    paste(round(head(t3_test, 4)), collapse = " | ")))

t10_test <- sim_treated(10000, tv0_d10, params_test, times_d[ok_d10])
cat(sprintf("10 mg/kg : %s\n",
    if (any(is.na(t10_test))) "ECHEC (NA)" else
    paste(round(head(t10_test, 4)), collapse = " | ")))

obj0 <- objective_simeoni(log(init_pd))
cat(sprintf("Objectif initial : %.4f %s\n", obj0,
    if (obj0 >= 1e11) "— ECHEC (1e12 = toutes simulations NA)" else "— OK"))

# =============================================================================
# 7. OPTIMISATION — DEoptim (évolution différentielle, optimiseur global)
#
#    2 paramètres libres : k1, k2  (λ0 et L1 fixés depuis le contrôle)
#    Bornes physiologiques en log-espace :
#      k1  ∈ [0.05, 4.0]   /j        (MTT = 4/k1 ∈ [1, 80] j)
#      k2  ∈ [1e-8, 1e-3]  L/µg/j
# =============================================================================

lower_log <- c(log(0.05), log(1e-8))
upper_log <- c(log(4.0),  log(1e-3))

cat("\nOptimisation DEoptim en cours (λ0 et L1 fixés, 2 paramètres libres)...\n")

set.seed(42)
fit_de <- DEoptim(
  fn      = objective_simeoni,
  lower   = lower_log,
  upper   = upper_log,
  control = DEoptim.control(
    NP      = 60,      # population : 30 × nb_params (2)
    itermax = 600,
    F       = 0.8,
    CR      = 0.9,
    trace   = 100,
    reltol  = 1e-8,
    steptol = 150
  )
)

# Affinage local
fit_local <- nlminb(
  start     = fit_de$optim$bestmem,
  objective = objective_simeoni,
  lower     = lower_log,
  upper     = upper_log,
  control   = list(eval.max = 3000, iter.max = 1000,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

best_log       <- fit_local$par
best_pd        <- exp(best_log)
names(best_pd) <- c("k1", "k2")

# =============================================================================
# 8. PARAMÈTRES DÉRIVÉS
# =============================================================================

# TSC (phase exponentielle) : L0_fixed / k2
tsc <- L0_fixed / unname(best_pd["k2"])

# MTT : 4 / k1
mtt <- 4 / unname(best_pd["k1"])

cat("\n=== Paramètres PD estimés (Simeoni + DEoptim) ===\n")
cat(sprintf("λ0 (L0) = %.6f /j         (FIXÉ — régression contrôle j0–j%d)\n", L0_fixed, L0_EXPO_DAY))
cat(sprintf("λ1 (L1) = %.2f mm³/j      (FIXÉ — plateau contrôle)\n",          L1_fixed))
cat(sprintf("k1      = %.5f /j         (transit rate)\n",                best_pd["k1"]))
cat(sprintf("k2      = %.2e L/µg/j    (potence drogue)\n",               best_pd["k2"]))
cat("─────────────────────────────────────────────────────────\n")
cat(sprintf("TSC   = %.1f µg/L   (= λ0/k2, concentration statique tumorale)\n", tsc))
cat(sprintf("MTT   = %.1f j      (= 4/k1, délai de réponse)\n", mtt))
cat(sprintf("Objectif DEoptim = %.6f\n", fit_de$optim$bestval))
cat(sprintf("Objectif final   = %.6f\n", fit_local$objective))

# =============================================================================
# 9. SIMULATION FINALE
# =============================================================================

params_best <- c(pk_fixed,
                 L0 = L0_fixed,
                 L1 = L1_fixed,
                 k1 = unname(best_pd["k1"]),
                 k2 = unname(best_pd["k2"]))

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim <- sim_ctrl_fn(L0_fixed, L1_fixed, tv0_ctrl, times_sim)
pred_iso_sim  <- sim_ctrl_fn(L0_fixed, L1_fixed, tv0_iso,  times_sim)
pred_d3_sim   <- sim_treated(3000,  tv0_d3,  params_best, times_sim)
pred_d10_sim  <- sim_treated(10000, tv0_d10, params_best, times_sim)

df_sim <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_sim, Groupe = "Contrôle"),
  data.frame(jour = times_sim, TV = pred_iso_sim,  Groupe = "Isotype 10mg/kg"),
  data.frame(jour = times_sim, TV = pred_d3_sim,   Groupe = "3 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d10_sim,  Groupe = "10 mg/kg")
)

# Sépare les points contrôle valides (fit) des points censurés (post-j39)
ok_ctrl_cens <- ok_ctrl & times_d > CTRL_CENSOR_DAY

df_obs <- rbind(
  data.frame(jour = times_d[ok_ctrl_fit],  TV = tv_ctrl[ok_ctrl_fit],
             sem  = sem_ctrl[ok_ctrl_fit],   Groupe = "Contrôle",        censored = FALSE),
  data.frame(jour = times_d[ok_ctrl_cens], TV = tv_ctrl[ok_ctrl_cens],
             sem  = sem_ctrl[ok_ctrl_cens],  Groupe = "Contrôle",        censored = TRUE),
  data.frame(jour = times_d[ok_iso],       TV = tv_iso[ok_iso],
             sem  = sem_iso[ok_iso],         Groupe = "Isotype 10mg/kg", censored = FALSE),
  data.frame(jour = times_d[ok_d3],        TV = tv_d3[ok_d3],
             sem  = sem_d3[ok_d3],           Groupe = "3 mg/kg",         censored = FALSE),
  data.frame(jour = times_d[ok_d10],       TV = tv_d10[ok_d10],
             sem  = sem_d10[ok_d10],         Groupe = "10 mg/kg",        censored = FALSE)
)

niv <- c("Contrôle", "Isotype 10mg/kg", "3 mg/kg", "10 mg/kg")
df_sim$Groupe <- factor(df_sim$Groupe, levels = niv)
df_obs$Groupe <- factor(df_obs$Groupe, levels = niv)

# =============================================================================
# 10. GRAPHIQUE
# =============================================================================

cols <- c("Contrôle"       = "#888888",
          "Isotype 10mg/kg"= "#CC4444",
          "3 mg/kg"        = "#4393C3",
          "10 mg/kg"       = "#2166AC")
ymax <- max(df_obs$TV + df_obs$sem, na.rm = TRUE)

subtitle_txt <- paste0(
  "λ0 = ", round(L0_fixed, 4), " /j (fixé)  |  λ1 = ", round(L1_fixed, 0), " mm³/j (fixé)  |  ",
  "k2 = ", formatC(best_pd["k2"], digits = 3, format = "e"), " L/µg/j  |  ",
  "k1 = ", round(best_pd["k1"], 3), " /j  |  ",
  "TSC = ", round(tsc, 0), " µg/L  |  ",
  "MTT = ", round(mtt, 1), " j"
)

df_obs_fit  <- df_obs[!df_obs$censored, ]
df_obs_cens <- df_obs[df_obs$censored,  ]

p_simeoni <- ggplot() +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = CTRL_CENSOR_DAY, linetype = "dotted",
             color = "#888888", linewidth = 0.6) +
  annotate("text", x = CTRL_CENSOR_DAY + 0.5, y = ymax * 0.55,
           label = sprintf("Endpoint\ncontrôle\n(j%d)", CTRL_CENSOR_DAY),
           hjust = 0, size = 2.8, color = "#888888") +
  geom_line(data = df_sim[df_sim$Groupe != "Contrôle", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1) +
  geom_line(data = df_sim[df_sim$Groupe == "Contrôle", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1, linetype = "dashed") +
  geom_errorbar(data = df_obs_fit,
                aes(x = jour, ymin = TV - sem, ymax = TV + sem, color = Groupe),
                width = 0.8, linewidth = 0.5) +
  geom_point(data = df_obs_fit,
             aes(x = jour, y = TV, color = Groupe),
             size = 2.5) +
  # Points censurés contrôle : symbole creux pour signaler l'exclusion du fit
  geom_point(data = df_obs_cens,
             aes(x = jour, y = TV, color = Groupe),
             size = 2.5, shape = 1, stroke = 0.8) +
  annotate("point", x = dose_days, y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text",  x = max(dose_days) + 1.5, y = -ymax * 0.06,
           label = "= Traitement", hjust = 0, size = 3.2, color = "#CC0000") +
  scale_x_continuous(breaks = seq(0, max(times_d, na.rm = TRUE), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "Modèle PK/PD TGI — Simeoni (2004) — Fc-silent FGFR2-huBPA-LP1",
    subtitle = subtitle_txt,
    caption  = sprintf("Contrôle tronqué à j≤%d (endpoint humanitaire) — points creux = exclus du fit",
                       CTRL_CENSOR_DAY),
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm³)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 9, color = "grey50"),
    plot.caption    = element_text(size = 8,  color = "grey50", hjust = 0),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

print(p_simeoni)
ggsave("scripts/plot_PKPD_simeoni_FGFR2.png", p_simeoni,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_PKPD_simeoni_FGFR2.png\n")

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

simeoni_results <- list(
  pk_fixed          = pk_fixed,
  L0                = L0_fixed,
  L0_fixed          = TRUE,
  L1                = L1_fixed,
  L1_fixed          = TRUE,
  k1                = unname(best_pd["k1"]),
  k2                = unname(best_pd["k2"]),
  TSC_ugL           = tsc,
  MTT_days          = mtt,
  ctrl_censor_day   = CTRL_CENSOR_DAY,
  objective_deoptim = fit_de$optim$bestval,
  objective_final   = fit_local$objective
)

save(simeoni_results, file = "scripts/resultats_PKPD_simeoni_FGFR2.RData")
cat("Résultats → scripts/resultats_PKPD_simeoni_FGFR2.RData\n")
