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
# Paramètres PD estimés (libres) : λ1, k1, k2  | λ0 fixé (log-espace)
# Unités : temps en jours | dose en µg/kg | concentration en µg/L
#
# Schéma : Q2W × 3 (j0, j14, j28) — VARIANTE 3 doses
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

cat(sprintf("\nTV initiale (j0) : Ctrl=%.0f | Isotype=%.0f | 3mg=%.0f | 10mg=%.0f mm³\n",
            tv0_ctrl, tv0_iso, tv0_d3, tv0_d10))
cat(sprintf("Contrôle : %d points valides\n", sum(ok_ctrl)))

# =============================================================================
# 3. ÉQUATIONS DU MODÈLE (deSolve)
# =============================================================================

PSI <- 20   # exposant de la fonction de croissance (fixé, Simeoni 2004)

simeoni_rhs <- function(t, state, parms) {
  # Indexation positionnelle + na.rm : robuste aux NaN/Inf explorés par lsoda
  A1 <- max(state[1], 0, na.rm = TRUE); A2 <- max(state[2], 0, na.rm = TRUE)
  x1 <- max(state[3], 0, na.rm = TRUE); x2 <- max(state[4], 0, na.rm = TRUE)
  x3 <- max(state[5], 0, na.rm = TRUE); x4 <- max(state[6], 0, na.rm = TRUE)

  CL <- parms["CL"]; V1 <- parms["V1"]
  V2 <- parms["V2"]; Q  <- parms["Q"]
  L0 <- parms["L0"]; L1 <- parms["L1"]
  k1 <- parms["k1"]; k2 <- parms["k2"]

  # PK
  C1  <- A1 / V1
  dA1 <- -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
  dA2 <-  (Q/V1)*A1 - (Q/V2)*A2

  # PD — fonction de croissance de Simeoni
  w <- x1 + x2 + x3 + x4
  if (is.finite(w) && w > 1e-12) {
    ratio <- (L0 * w / L1)^PSI
    growth_rate <- if (is.finite(ratio)) L0 / (1 + ratio)^(1/PSI)
                   else as.numeric(L1) / w
  } else {
    growth_rate <- as.numeric(L0)
  }

  dx1 <- (growth_rate - k2 * C1) * x1
  dx2 <- k2 * C1 * x1 - k1 * x2
  dx3 <- k1 * x2 - k1 * x3
  dx4 <- k1 * x3 - k1 * x4

  list(c(dA1, dA2, dx1, dx2, dx3, dx4))
}

# Contrôle : même modèle sans drogue (A1=A2=0, k2 ignoré)
simeoni_ctrl_rhs <- function(t, state, parms) {
  x1 <- max(state[1], 0, na.rm = TRUE)
  x2 <- max(state[2], 0, na.rm = TRUE)
  x3 <- max(state[3], 0, na.rm = TRUE)
  x4 <- max(state[4], 0, na.rm = TRUE)
  L0 <- parms["L0"]; L1 <- parms["L1"]

  w <- x1 + x2 + x3 + x4
  if (is.finite(w) && w > 1e-12) {
    ratio <- (L0 * w / L1)^PSI
    growth_rate <- if (is.finite(ratio)) L0 / (1 + ratio)^(1/PSI)
                   else as.numeric(L1) / w
  } else {
    growth_rate <- as.numeric(L0)
  }

  dx1 <- growth_rate * x1
  dx2 <- 0; dx3 <- 0; dx4 <- 0

  list(c(dx1, dx2, dx3, dx4))
}

# =============================================================================
# 4. FONCTIONS DE SIMULATION
# =============================================================================

dose_days <- c(0, 14, 28)   # Q2W × 3 doses

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
        atol  = 1e-2, rtol = 1e-2
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
# 5. FONCTIONS OBJECTIVES — FIT EN 2 ÉTAPES
# =============================================================================

L1_FIXED <- 1e6   # mm³/j — non identifiable en phase exponentielle

obj_ctrl <- function(log_L0) {
  L0   <- exp(log_L0)
  pred <- sim_ctrl_fn(L0, L1_FIXED, tv0_ctrl, times_d[ok_ctrl])
  if (any(is.na(pred))) return(1e12)
  obs  <- tv_ctrl[ok_ctrl]
  wt   <- 1 / sem_ctrl[ok_ctrl]^2
  pred <- pmax(pred, 0.1)
  sum(wt * (log(obs) - log(pred))^2, na.rm = TRUE)
}

obj_treated <- function(log_k2) {
  k2 <- exp(log_k2)
  params <- c(pk_fixed, L0 = L0_ctrl, L1 = L1_FIXED, k1 = K1_FIXED, k2 = k2)

  pred_d3 <- sim_treated(3000,  tv0_d3,  params, times_d[ok_d3])
  if (any(is.na(pred_d3))) return(1e12)
  pred_d3 <- pmax(pred_d3, 0.1)
  obs_d3  <- tv_d3[ok_d3];  wt_d3  <- 1 / sem_d3[ok_d3]^2

  pred_d10 <- sim_treated(10000, tv0_d10, params, times_d[ok_d10])
  if (any(is.na(pred_d10))) return(1e12)
  pred_d10 <- pmax(pred_d10, 0.1)
  obs_d10  <- tv_d10[ok_d10]; wt_d10 <- 1 / sem_d10[ok_d10]^2

  sum(wt_d3  * (log(obs_d3)  - log(pred_d3))^2,  na.rm = TRUE) +
  sum(wt_d10 * (log(obs_d10) - log(pred_d10))^2, na.rm = TRUE)
}

# =============================================================================
# 6. ÉTAPE 1 — FIT CONTRÔLE : estimation de λ0
# =============================================================================

cat("\n══════════════════════════════════════════════════════════\n")
cat("ÉTAPE 1 — Fit contrôle seul (λ0 libre, λ1 fixé à 1e6)\n")
cat("══════════════════════════════════════════════════════════\n")

lm_ctrl  <- lm(log(tv_ctrl[ok_ctrl]) ~ times_d[ok_ctrl])
L0_init  <- max(coef(lm_ctrl)[2], 0.005)
cat(sprintf("λ0 initial (régression log-linéaire) : %.5f /j\n", L0_init))

fit_ctrl <- nlminb(
  start     = log(L0_init),
  objective = obj_ctrl,
  lower     = log(0.001),
  upper     = log(1.0),
  control   = list(eval.max = 2000, iter.max = 1000, rel.tol = 1e-12)
)

L0_ctrl <- exp(fit_ctrl$par)
cat(sprintf("λ0 (L0) = %.5f /j   (t½ = %.1f j)\n", L0_ctrl, log(2)/L0_ctrl))
cat(sprintf("λ1 (L1) = %.0f mm³/j  (FIXÉ)\n", L1_FIXED))
cat(sprintf("Objectif contrôle = %.8f\n", fit_ctrl$objective))

# =============================================================================
# 7. ÉTAPE 2 — FIT GROUPES TRAITÉS : estimation de k1 et k2
# =============================================================================

cat("\n══════════════════════════════════════════════════════════\n")
cat("ÉTAPE 2 — Fit traités seuls (k2 libre ; λ0, λ1, k1 gelés)\n")
cat("══════════════════════════════════════════════════════════\n")

K1_FIXED <- 4 / 14
cat(sprintf("k1 (K1_FIXED) = %.4f /j  (MTT = %.0f j — FIXÉ, non identifiable)\n",
            K1_FIXED, 4/K1_FIXED))

Cmax_10 <- 10000 / as.numeric(pk_fixed["V1"])
k2_init <- L0_ctrl / Cmax_10
cat(sprintf("k2 init = %.2e L/µg/j (TSC ≈ %.0f µg/L)\n", k2_init, L0_ctrl/k2_init))

lower_k2 <- log(1e-8)
upper_k2 <- log(1e-3)

set.seed(42)
de_treated <- DEoptim(
  fn      = obj_treated,
  lower   = lower_k2,
  upper   = upper_k2,
  control = DEoptim.control(
    NP      = 40,
    itermax = 400,
    F       = 0.8, CR = 0.9,
    trace   = 100,
    reltol  = 1e-8, steptol = 150
  )
)

loc_treated <- nlminb(
  start     = de_treated$optim$bestmem,
  objective = obj_treated,
  lower     = lower_k2,
  upper     = upper_k2,
  control   = list(eval.max = 2000, iter.max = 1000, rel.tol = 1e-12)
)

if (!is.finite(loc_treated$objective) ||
    loc_treated$objective > de_treated$optim$bestval * 10) {
  best_log_k2 <- de_treated$optim$bestmem
  best_obj_k  <- de_treated$optim$bestval
} else {
  best_log_k2 <- loc_treated$par
  best_obj_k  <- loc_treated$objective
}

k1_best <- K1_FIXED
k2_best <- exp(unname(best_log_k2))

# =============================================================================
# 8. PARAMÈTRES DÉRIVÉS
# =============================================================================

tsc <- L0_ctrl / k2_best
mtt <- 4 / k1_best

cat("\n=== Paramètres PD finaux (Simeoni — fit 2 étapes) ===\n")
cat(sprintf("λ0 (L0) = %.6f /j         (étape 1 — fit contrôle)\n", L0_ctrl))
cat(sprintf("λ1 (L1) = %.0f mm³/j      (FIXÉ)\n",                    L1_FIXED))
cat(sprintf("k1      = %.5f /j         (FIXÉ — MTT = %.0f j)\n",      k1_best, 4/k1_best))
cat(sprintf("k2      = %.2e L/µg/j    (étape 2 — fit traités)\n",    k2_best))
cat("─────────────────────────────────────────────────────────\n")
cat(sprintf("TSC = %.1f µg/L  (= λ0/k2)\n", tsc))
cat(sprintf("MTT = %.1f j     (= 4/k1)\n",  mtt))
cat(sprintf("Objectif contrôle = %.8f\n", fit_ctrl$objective))
cat(sprintf("Objectif traités  = %.8f\n", best_obj_k))

# =============================================================================
# 9. SIMULATION FINALE
# =============================================================================

params_best <- c(pk_fixed,
                 L0 = L0_ctrl,
                 L1 = L1_FIXED,
                 k1 = k1_best,
                 k2 = k2_best)

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim <- sim_ctrl_fn(L0_ctrl, L1_FIXED, tv0_ctrl, times_sim)
pred_iso_sim  <- sim_ctrl_fn(L0_ctrl, L1_FIXED, tv0_iso,  times_sim)
pred_d3_sim   <- sim_treated(3000,  tv0_d3,  params_best, times_sim)
pred_d10_sim  <- sim_treated(10000, tv0_d10, params_best, times_sim)

df_sim <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_sim, Groupe = "Contrôle"),
  data.frame(jour = times_sim, TV = pred_iso_sim,  Groupe = "Isotype 10mg/kg"),
  data.frame(jour = times_sim, TV = pred_d3_sim,   Groupe = "3 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d10_sim,  Groupe = "10 mg/kg")
)

df_obs <- rbind(
  data.frame(jour = times_d[ok_ctrl], TV = tv_ctrl[ok_ctrl],
             sem  = sem_ctrl[ok_ctrl], Groupe = "Contrôle"),
  data.frame(jour = times_d[ok_iso],  TV = tv_iso[ok_iso],
             sem  = sem_iso[ok_iso],   Groupe = "Isotype 10mg/kg"),
  data.frame(jour = times_d[ok_d3],   TV = tv_d3[ok_d3],
             sem  = sem_d3[ok_d3],     Groupe = "3 mg/kg"),
  data.frame(jour = times_d[ok_d10],  TV = tv_d10[ok_d10],
             sem  = sem_d10[ok_d10],   Groupe = "10 mg/kg")
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
  "λ0 = ", round(L0_ctrl, 4), " /j (fit ctrl)  |  ",
  "k2 = ", formatC(k2_best, digits = 3, format = "e"), " L/µg/j  |  ",
  "k1 = ", round(k1_best, 3), " /j  |  ",
  "TSC = ", round(tsc, 0), " µg/L  |  ",
  "MTT = ", round(mtt, 1), " j"
)

p_simeoni <- ggplot() +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(data = df_sim[df_sim$Groupe != "Contrôle", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1) +
  geom_line(data = df_sim[df_sim$Groupe == "Contrôle", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1, linetype = "dashed") +
  geom_errorbar(data = df_obs,
                aes(x = jour, ymin = TV - sem, ymax = TV + sem, color = Groupe),
                width = 0.8, linewidth = 0.5) +
  geom_point(data = df_obs,
             aes(x = jour, y = TV, color = Groupe),
             size = 2.5) +
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
ggsave("scripts/plot_PKPD_simeoni_FGFR2_3doses.png", p_simeoni,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_PKPD_simeoni_FGFR2_3doses.png\n")

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

simeoni_results <- list(
  pk_fixed           = pk_fixed,
  L0                 = L0_ctrl,
  L0_source          = "fit contrôle étape 1",
  L1                 = L1_FIXED,
  L1_fixed           = TRUE,
  k1                 = k1_best,
  k1_fixed           = TRUE,
  k2                 = k2_best,
  TSC_ugL            = tsc,
  MTT_days           = mtt,
  objective_ctrl     = fit_ctrl$objective,
  objective_traites  = best_obj_k
)

save(simeoni_results, file = "scripts/resultats_PKPD_simeoni_FGFR2_3doses.RData")
cat("Résultats → scripts/resultats_PKPD_simeoni_FGFR2_3doses.RData\n")
