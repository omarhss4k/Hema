# =============================================================================
# Modèle PK/PD TGI — Simeoni (2004) — Fc-silent FGFR2-huBPA-LP1
# Cohorte : SNU (xénogreffe souris)
#
# PK  : 2 compartiments IV bolus (paramètres fixés — données souris)
# PD  : Modèle de Simeoni — compartiments de transit
#
#   dx1/dt = [g(w)/w - k2·C1] · x1
#   dx2/dt =  k2·C1·x1 - k1·x2
#   dx3/dt =  k1·x2 - k1·x3
#   dx4/dt =  k1·x3 - k1·x4
#   w = x1+x2+x3+x4
#   g(w) = λ0·w / (1 + (λ0·w/λ1)^ψ)^(1/ψ)   ψ=20
#
#   TSC = λ0/k2   (µg/L — doit être entre Cmax(1.2mg/kg) et Cmax(5.6mg/kg))
#   MTT = 4/k1    (jours)
#
# Ajustement SÉQUENTIEL :
#   Étape A : L0 et L1 sur le groupe véhicule uniquement
#   Étape B : k1 et k2 (TSC) sur les 3 groupes traités, L0/L1 fixés
#
# Paramètres PD estimés (4) : λ0, λ1, k1, k2
#
# Schéma : dose unique (j0)
# Groupes : Vehicule | 1.2 mg/kg (Group 4) | 5.6 mg/kg (Group 5) | 11.8 mg/kg (Group 7)
# =============================================================================

library(deSolve)
library(DEoptim)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS (Importation depuis le CSV)
# =============================================================================

pk_data <- read.csv("pk_resultats_20260611_110140.csv", stringsAsFactors = FALSE)

# *1000 : conversion mL/kg → L/kg (ou équivalent unité)
# *24   : conversion par heure → par jour
pk_fixed <- c(
  CL = (pk_data$CL[1] * 1000) * 24,
  V1 = pk_data$V1[1] * 1000,
  V2 = pk_data$V2[1] * 1000,
  Q  = (pk_data$Q[1]  * 1000) * 24
)

cat("=== Paramètres PK fixés (souris, jours) ===\n")
cat(sprintf("CL = %.5e L/j/kg\n", pk_fixed["CL"]))
cat(sprintf("V1 = %.5e L/kg\n",   pk_fixed["V1"]))
cat(sprintf("V2 = %.5e L/kg\n",   pk_fixed["V2"]))
cat(sprintf("Q  = %.5e L/j/kg\n", pk_fixed["Q"]))

# Extraction des Cmax depuis le CSV (en ng/mL = µg/L)
Cmax_1p2  <- pk_data$cmax[pk_data$Animal == "A"]
Cmax_5p6  <- pk_data$cmax[pk_data$Animal == "C"]
Cmax_11p8 <- pk_data$cmax[pk_data$Animal == "D"]

cat(sprintf("Cmax lues (µg/L) : 1.2mg=%.0f | 5.9mg=%.0f | 11.8mg=%.0f\n",
            Cmax_1p2, Cmax_5p6, Cmax_11p8))

# =============================================================================
# 2. DONNÉES TUMORALES — SNU
# =============================================================================

raw_snu <- suppressMessages(
  read_xlsx("Tumor_volum_SNU.xlsx", col_names = FALSE, .name_repair = "minimal")
)

header_rows <- which(grepl("Time", trimws(as.character(raw_snu[[1]])), ignore.case = TRUE))
stopifnot("2 blocs attendus (moyennes + SEM)" = length(header_rows) == 2)

b1 <- raw_snu[seq(header_rows[1] + 1, header_rows[2] - 1), ]
b1 <- b1[!is.na(suppressWarnings(as.numeric(as.character(b1[[1]])))), ]

b2 <- raw_snu[seq(header_rows[2] + 1, nrow(raw_snu)), ]
b2 <- b2[!is.na(suppressWarnings(as.numeric(as.character(b2[[1]])))), ]

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

cat(sprintf("\nTV initiale (j0) : Ctrl=%.0f | 1.2mg=%.0f | 5.6mg=%.0f | 11.8mg=%.0f mm³\n",
            tv0_ctrl, tv0_d1p2, tv0_d5p6, tv0_d11p8))

# =============================================================================
# 3. ÉQUATIONS DU MODÈLE DE SIMEONI  (inchangées)
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
# 4. FONCTIONS DE SIMULATION  (inchangées)
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
  t_all <- sort(unique(c(0, times_out)))
  out <- tryCatch(
    as.data.frame(lsoda(
      y     = c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
      times = t_all,
      func  = simeoni_rhs,
      parms = params,
      atol  = 1e-6, rtol = 1e-6
    )),
    error = function(e) NULL
  )
  if (is.null(out) || any(is.na(out$x1))) return(rep(NA_real_, length(times_out)))
  w_tot <- pmax(out$x1, 0) + pmax(out$x2, 0) + pmax(out$x3, 0) + pmax(out$x4, 0)
  approx(out$time, w_tot, xout = times_out, rule = 2)$y
}

# =============================================================================
# 5A. ÉTAPE A — Ajustement croissance : L0 et L1 sur le groupe véhicule
#
#   Bornes L0 : temps de doublement ∈ [3, 14 j]  →  L0 ∈ [ln2/14, ln2/3]
#   Bornes L1 : volume maximum plausible ∈ [TV_max, 50×TV_max]
# =============================================================================

cat("\n=== ÉTAPE A : Ajustement groupe Véhicule (L0, L1) ===\n")

# Valeurs initiales
lm_ctrl <- lm(log(tv_ctrl[ok_ctrl]) ~ times_d[ok_ctrl])
L0_init <- max(coef(lm_ctrl)[2], 0.005)
L1_init <- max(tv_ctrl[ok_ctrl], na.rm = TRUE) * L0_init * 50

# Bornes biologiquement contraintes
L0_lower <- log(0.05)   # doublement ≈ 14 j
L0_upper <- log(0.25)   # doublement ≈  3 j
L1_lower <- log(max(tv_ctrl[ok_ctrl], na.rm = TRUE))
L1_upper <- log(max(tv_ctrl[ok_ctrl], na.rm = TRUE) * 200)

objective_vehicule <- function(logpar) {
  L0 <- exp(logpar[1])
  L1 <- exp(logpar[2])
  if (L0 <= 0 || L1 <= 0) return(1e12)

  pred <- sim_ctrl_fn(L0, L1, tv0_ctrl, times_d[ok_ctrl])
  if (any(is.na(pred))) return(1e12)
  pred <- pmax(pred, 0.1)

  sum(1/sem_ctrl[ok_ctrl]^2 * (log(tv_ctrl[ok_ctrl]) - log(pred))^2,
      na.rm = TRUE)
}

set.seed(42)
fit_de_A <- DEoptim(
  fn      = objective_vehicule,
  lower   = c(L0_lower, L1_lower),
  upper   = c(L0_upper, L1_upper),
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
  lower     = c(L0_lower, L1_lower),
  upper     = c(L0_upper, L1_upper),
  control   = list(eval.max = 2000, iter.max = 800,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

L0_est <- exp(fit_local_A$par[1])
L1_est <- exp(fit_local_A$par[2])

cat(sprintf("→ L0 (λ0)  = %.6f /j   (doublement = %.1f j)\n",
            L0_est, log(2)/L0_est))
cat(sprintf("→ L1 (λ1)  = %.2f mm³/j\n", L1_est))
cat(sprintf("   Objectif DEoptim = %.6f\n", fit_de_A$optim$bestval))
cat(sprintf("   Objectif final   = %.6f\n", fit_local_A$objective))

# =============================================================================
# 5B. ÉTAPE B — Ajustement efficacité : k1 et k2 sur groupes traités
#
#   L0 et L1 fixés à partir de l'Étape A
#   Bornes TSC : [Cmax_1p2/3, Cmax_5p6]  →  borne inf élargie pour faible dose
#   Bornes k1  : MTT ∈ [2, 21 j]  →  k1 ∈ [4/21, 4/2]
# =============================================================================

cat("\n=== ÉTAPE B : Ajustement groupes traités (k1, k2) ===\n")
cat(sprintf("L0 fixé = %.6f /j | L1 fixé = %.2f mm³/j\n", L0_est, L1_est))

# Bornes k2 dérivées des bornes TSC
tsc_lower <- Cmax_1p2 / 3     # élargi pour capter légère régression à faible dose
tsc_upper <- Cmax_5p6
k2_lower  <- L0_est / tsc_upper
k2_upper  <- L0_est / tsc_lower

cat(sprintf("Bornes TSC : [%.0f, %.0f] µg/L\n", tsc_lower, tsc_upper))
cat(sprintf("→ Bornes k2 : [%.2e, %.2e] L/µg/j\n", k2_lower, k2_upper))

objective_traites <- function(logpar) {
  k1 <- exp(logpar[1])
  k2 <- exp(logpar[2])
  if (k1 <= 0 || k2 <= 0) return(1e12)

  params_all <- c(pk_fixed, L0 = L0_est, L1 = L1_est, k1 = k1, k2 = k2)

  pred_d1p2 <- sim_treated(1200,  tv0_d1p2,  params_all, times_d[ok_d1p2])
  if (any(is.na(pred_d1p2))) return(1e12)
  pred_d1p2 <- pmax(pred_d1p2, 0.1)

  pred_d5p6 <- sim_treated(5600,  tv0_d5p6,  params_all, times_d[ok_d5p6])
  if (any(is.na(pred_d5p6))) return(1e12)
  pred_d5p6 <- pmax(pred_d5p6, 0.1)

  pred_d11p8 <- sim_treated(11800, tv0_d11p8, params_all, times_d[ok_d11p8])
  if (any(is.na(pred_d11p8))) return(1e12)
  pred_d11p8 <- pmax(pred_d11p8, 0.1)

  sum(1/sem_d1p2[ok_d1p2]^2   * (log(tv_d1p2[ok_d1p2])   - log(pred_d1p2))^2,  na.rm = TRUE) +
  sum(1/sem_d5p6[ok_d5p6]^2   * (log(tv_d5p6[ok_d5p6])   - log(pred_d5p6))^2,  na.rm = TRUE) +
  sum(1/sem_d11p8[ok_d11p8]^2 * (log(tv_d11p8[ok_d11p8]) - log(pred_d11p8))^2, na.rm = TRUE)
}

k1_init <- 4 / 7
lower_B <- c(log(4/21),    log(k2_lower))
upper_B <- c(log(4/2),     log(k2_upper))

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
                   rel.tol = 1e-12, x.tol = 1e-12)
)

k1_est <- exp(fit_local_B$par[1])
k2_est <- exp(fit_local_B$par[2])

cat(sprintf("→ k1 = %.5f /j   (MTT = %.1f j)\n", k1_est, 4/k1_est))
cat(sprintf("→ k2 = %.3e L/µg/j\n", k2_est))
cat(sprintf("   Objectif DEoptim = %.6f\n", fit_de_B$optim$bestval))
cat(sprintf("   Objectif final   = %.6f\n", fit_local_B$objective))

# =============================================================================
# 6. PARAMÈTRES DÉRIVÉS FINAUX
# =============================================================================

tsc <- L0_est / k2_est
mtt <- 4 / k1_est

cat("\n=== Paramètres PD estimés — SNU (ajustement séquentiel) ===\n")
cat(sprintf("λ0 (L0) = %.6f /j         [Étape A]\n", L0_est))
cat(sprintf("λ1 (L1) = %.2f mm³/j      [Étape A]\n", L1_est))
cat(sprintf("k1      = %.5f /j         [Étape B]\n", k1_est))
cat(sprintf("k2      = %.3e L/µg/j    [Étape B]\n",  k2_est))
cat("──────────────────────────────────────────────────────────\n")
cat(sprintf("TSC = %.0f µg/L  (Cmax_1.2=%.0f | Cmax_5.6=%.0f)\n",
            tsc, Cmax_1p2, Cmax_5p6))
cat(sprintf("MTT = %.1f j\n", mtt))

# =============================================================================
# 7. SIMULATION FINALE
# =============================================================================

params_best <- c(pk_fixed,
                 L0 = L0_est, L1 = L1_est,
                 k1 = k1_est, k2 = k2_est)

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim  <- sim_ctrl_fn(L0_est, L1_est, tv0_ctrl,  times_sim)
pred_d1p2_sim  <- sim_treated(1200,  tv0_d1p2,  params_best, times_sim)
pred_d5p6_sim  <- sim_treated(5600,  tv0_d5p6,  params_best, times_sim)
pred_d11p8_sim <- sim_treated(11800, tv0_d11p8, params_best, times_sim)

df_sim <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_sim,  Groupe = "Vehicule"),
  data.frame(jour = times_sim, TV = pred_d1p2_sim,  Groupe = "1.2 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d5p6_sim,  Groupe = "5.6 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d11p8_sim, Groupe = "11.8 mg/kg")
)

df_obs <- rbind(
  data.frame(jour = times_d[ok_ctrl],  TV = tv_ctrl[ok_ctrl],    sem = sem_ctrl[ok_ctrl],   Groupe = "Vehicule"),
  data.frame(jour = times_d[ok_d1p2],  TV = tv_d1p2[ok_d1p2],   sem = sem_d1p2[ok_d1p2],  Groupe = "1.2 mg/kg"),
  data.frame(jour = times_d[ok_d5p6],  TV = tv_d5p6[ok_d5p6],   sem = sem_d5p6[ok_d5p6],  Groupe = "5.6 mg/kg"),
  data.frame(jour = times_d[ok_d11p8], TV = tv_d11p8[ok_d11p8], sem = sem_d11p8[ok_d11p8], Groupe = "11.8 mg/kg")
)

niv <- c("Vehicule", "1.2 mg/kg", "5.6 mg/kg", "11.8 mg/kg")
df_sim$Groupe <- factor(df_sim$Groupe, levels = niv)
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
  "Étape A — λ0 = ", round(L0_est, 4), " /j  |  ",
  "λ1 = ", round(L1_est, 0), " mm³/j  |  ",
  "Étape B — k2 = ", formatC(k2_est, digits = 3, format = "e"), " L/µg/j  |  ",
  "k1 = ", round(k1_est, 3), " /j  |  ",
  "TSC = ", round(tsc, 0), " µg/L  |  MTT = ", round(mtt, 1), " j"
)

p_simeoni <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(data = df_sim[df_sim$Groupe != "Vehicule", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1) +
  geom_line(data = df_sim[df_sim$Groupe == "Vehicule", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1, linetype = "dashed") +
  geom_errorbar(data = df_obs,
                aes(x = jour, ymin = TV - sem, ymax = TV + sem, color = Groupe),
                width = 0.8, linewidth = 0.5) +
  geom_point(data = df_obs,
             aes(x = jour, y = TV, color = Groupe),
             size = 2.5) +
  annotate("point", x = 0, y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text",  x = 1.5, y = -ymax * 0.06,
           label = "= Dose unique (j0)", hjust = 0, size = 3.2, color = "#CC0000") +
  scale_x_continuous(breaks = seq(0, max(times_d, na.rm = TRUE), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "Modèle PK/PD TGI — Simeoni (2004) — Fc-silent FGFR2-huBPA-LP1 — SNU",
    subtitle = subtitle_txt,
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm³)",
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
cat("\nGraphique → scripts/plot_PKPD_simeoni_SNU.png\n")

# =============================================================================
# 9. SAUVEGARDE
# =============================================================================

simeoni_results_SNU <- list(
  pk_fixed              = pk_fixed,
  L0                    = L0_est,
  L1                    = L1_est,
  k1                    = k1_est,
  k2                    = k2_est,
  TSC_ugL               = tsc,
  MTT_days              = mtt,
  objective_A_deoptim   = fit_de_A$optim$bestval,
  objective_A_final     = fit_local_A$objective,
  objective_B_deoptim   = fit_de_B$optim$bestval,
  objective_B_final     = fit_local_B$objective
)

save(simeoni_results_SNU, file = "scripts/resultats_PKPD_simeoni_SNU.RData")
cat("Résultats → scripts/resultats_PKPD_simeoni_SNU.RData\n")
