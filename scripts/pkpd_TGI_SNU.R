# =============================================================================
# Modèle PK/PD TGI — Fc-silent FGFR2-huBPA-LP1
# Cohorte : SNU (xénogreffe souris)
#
# PK  : 2 compartiments IV bolus (paramètres fixés — données souris)
#         dA1/dt = -(CL/V1 + Q/V1)·A1 + (Q/V2)·A2
#         dA2/dt =  (Q/V1)·A1 - (Q/V2)·A2
#         C1 = A1/V1
#
# PD  : Modèle Emax simple
#         dTV/dt = kg · TV − ki · [C1/(EC50 + C1)] · TV
#
#   TSC = kg · EC50 / (ki − kg)   [µg/L, si ki > kg]
#   TSC doit être entre Cmax(1.2 mg/kg) et Cmax(5.6 mg/kg)
#
# Paramètres PD estimés (3) : kg, ki, EC50
#
# Schéma : dose unique (j0)
# Groupes : Vehicule | 1.2 mg/kg (Group 4) | 5.6 mg/kg (Group 5) | 11.8 mg/kg (Group 7)
#
# Structure Excel (Tumor_volum_SNU.xlsx) — 2 blocs séparés par une ligne vide :
#   Bloc 1 : moyennes  — colonnes A=Temps, B=Vehicule, C=Group4, D=Group5, F=Group7
#   Bloc 2 : SEM      — mêmes colonnes
# =============================================================================

library(deSolve)
library(DEoptim)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS (données souris)
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2

pk_fixed <- c(
  CL = unname(pk2comp_rxode2$CL) * 24,
  V1 = unname(pk2comp_rxode2$V1),
  V2 = unname(pk2comp_rxode2$V2),
  Q  = unname(pk2comp_rxode2$Q)  * 24
)

cat("=== Paramètres PK fixés (souris, jours) ===\n")
cat(sprintf("CL = %.5f L/j/kg\n", pk_fixed["CL"]))
cat(sprintf("V1 = %.5f L/kg\n",   pk_fixed["V1"]))
cat(sprintf("V2 = %.5f L/kg\n",   pk_fixed["V2"]))
cat(sprintf("Q  = %.5f L/j/kg\n", pk_fixed["Q"]))

V1 <- as.numeric(pk_fixed["V1"])

# Cmax théorique IV bolus (µg/L) = dose(µg/kg) / V1(L/kg)
Cmax_1p2  <- 1200  / V1
Cmax_5p6  <- 5600  / V1
Cmax_11p8 <- 11800 / V1

cat(sprintf("\nCmax(1.2  mg/kg) = %.0f µg/L\n", Cmax_1p2))
cat(sprintf("Cmax(5.6  mg/kg) = %.0f µg/L\n", Cmax_5p6))
cat(sprintf("Cmax(11.8 mg/kg) = %.0f µg/L\n", Cmax_11p8))

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
tv_d1p2   <- num_col(b1, 3)   # 1.2 mg/kg
tv_d5p6   <- num_col(b1, 4)   # 5.6 mg/kg
tv_d11p8  <- num_col(b1, 6)   # 11.8 mg/kg  (col 5 vide)

sem_ctrl  <- fix_sem(num_col(b2, 2))
sem_d1p2  <- fix_sem(num_col(b2, 3))
sem_d5p6  <- fix_sem(num_col(b2, 4))
sem_d11p8 <- fix_sem(num_col(b2, 6))

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
# 3. ÉQUATIONS DU MODÈLE PK/PD (2-cpt PK + Emax PD)
# =============================================================================

CL_pk <- as.numeric(pk_fixed["CL"])
V2_pk <- as.numeric(pk_fixed["V2"])
Q_pk  <- as.numeric(pk_fixed["Q"])

pkpd_rhs <- function(t, state, parms) {
  A1 <- max(state["A1"], 0)
  A2 <- max(state["A2"], 0)
  TV <- max(state["TV"], 0)

  kg   <- parms["kg"]
  ki   <- parms["ki"]
  EC50 <- parms["EC50"]

  C1  <- A1 / V1
  Inh <- C1 / (EC50 + C1)

  dA1 <- -(CL_pk/V1 + Q_pk/V1) * A1 + (Q_pk/V2_pk) * A2
  dA2 <-  (Q_pk/V1) * A1 - (Q_pk/V2_pk) * A2
  dTV <-  kg * TV - ki * Inh * TV

  list(c(A1 = dA1, A2 = dA2, TV = dTV))
}

# =============================================================================
# 4. FONCTIONS DE SIMULATION
# =============================================================================

sim_ctrl_fn <- function(kg, tv0, times_out) tv0 * exp(kg * times_out)

sim_treated <- function(dose_ugkg, tv0, params, times_out) {
  t_all <- sort(unique(c(0, times_out)))
  out <- tryCatch(
    as.data.frame(lsoda(
      y     = c(A1 = dose_ugkg, A2 = 0, TV = tv0),
      times = t_all,
      func  = pkpd_rhs,
      parms = params
    )),
    error = function(e) NULL
  )
  if (is.null(out) || any(is.na(out$TV)) || any(out$TV < 0))
    return(rep(NA_real_, length(times_out)))
  approx(out$time, out$TV, xout = times_out, rule = 2)$y
}

# =============================================================================
# 5. FONCTION OBJECTIVE
#    Paramètres : kg, ki, EC50  (log-espace)
# =============================================================================

objective_pkpd <- function(logpar) {
  par  <- exp(logpar)
  kg   <- par[1]; ki <- par[2]; EC50 <- par[3]
  if (any(par <= 0)) return(1e12)

  params_all <- c(kg = kg, ki = ki, EC50 = EC50)

  pc <- sim_ctrl_fn(kg, tv0_ctrl, times_d[ok_ctrl])
  if (any(pc <= 0 | is.na(pc))) return(1e12)

  p1 <- sim_treated(1200,  tv0_d1p2,  params_all, times_d[ok_d1p2])
  if (any(is.na(p1) | p1 <= 0)) return(1e12)

  p5 <- sim_treated(5600,  tv0_d5p6,  params_all, times_d[ok_d5p6])
  if (any(is.na(p5) | p5 <= 0)) return(1e12)

  p11 <- sim_treated(11800, tv0_d11p8, params_all, times_d[ok_d11p8])
  if (any(is.na(p11) | p11 <= 0)) return(1e12)

  sum(1/sem_ctrl[ok_ctrl]^2    * (log(tv_ctrl[ok_ctrl])    - log(pc))^2,  na.rm = TRUE) +
  sum(1/sem_d1p2[ok_d1p2]^2   * (log(tv_d1p2[ok_d1p2])   - log(p1))^2,  na.rm = TRUE) +
  sum(1/sem_d5p6[ok_d5p6]^2   * (log(tv_d5p6[ok_d5p6])   - log(p5))^2,  na.rm = TRUE) +
  sum(1/sem_d11p8[ok_d11p8]^2 * (log(tv_d11p8[ok_d11p8]) - log(p11))^2, na.rm = TRUE)
}

# =============================================================================
# 6. VALEURS INITIALES ET BORNES
#
#   EC50 ciblé entre Cmax(1.2) et Cmax(5.6) → initialisation à la moyenne géométrique
#   Bornes dynamiques : [Cmax(1.2)/10 , Cmax(5.6)*10]
# =============================================================================

lm_ctrl  <- lm(log(tv_ctrl[ok_ctrl]) ~ times_d[ok_ctrl])
kg_init  <- max(coef(lm_ctrl)[2], 0.005)
ki_init  <- kg_init * 3

# EC50 (µg/L) — moyenne géométrique entre Cmax(1.2) et Cmax(5.6)
EC50_init    <- sqrt(Cmax_1p2 * Cmax_5p6)
EC50_lower   <- Cmax_1p2  / 10
EC50_upper   <- Cmax_5p6  * 10

cat("\n=== Valeurs initiales ===\n")
cat(sprintf("kg       = %.5f /j  (t½ tumeur = %.1f j)\n", kg_init, log(2)/kg_init))
cat(sprintf("ki       = %.5f /j\n", ki_init))
cat(sprintf("EC50     = %.0f µg/L  (entre Cmax(1.2)=%.0f et Cmax(5.6)=%.0f)\n",
            EC50_init, Cmax_1p2, Cmax_5p6))

lower_log <- c(log(0.005), log(0.01), log(EC50_lower))
upper_log <- c(log(1.0),   log(10.0), log(EC50_upper))

# =============================================================================
# 7. OPTIMISATION — DEoptim + affinage nlminb
# =============================================================================

cat("\nOptimisation DEoptim en cours...\n")
set.seed(42)
fit_de <- DEoptim(
  fn      = objective_pkpd,
  lower   = lower_log,
  upper   = upper_log,
  control = DEoptim.control(
    NP      = 120,
    itermax = 600,
    F       = 0.8,
    CR      = 0.9,
    trace   = 100,
    reltol  = 1e-8,
    steptol = 150
  )
)

fit_local <- nlminb(
  start     = fit_de$optim$bestmem,
  objective = objective_pkpd,
  lower     = lower_log,
  upper     = upper_log,
  control   = list(eval.max = 3000, iter.max = 1000,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

best_pd        <- exp(fit_local$par)
names(best_pd) <- c("kg", "ki", "EC50")

# =============================================================================
# 8. PARAMÈTRES DÉRIVÉS
# =============================================================================

ratio_ki_kg <- unname(best_pd["ki"] / best_pd["kg"])

tsc_conc <- if (best_pd["ki"] > best_pd["kg"])
  unname(best_pd["kg"] * best_pd["EC50"] / (best_pd["ki"] - best_pd["kg"]))
else NA_real_

tsc_mgkg <- if (!is.na(tsc_conc)) (tsc_conc * V1) / 1000 else NA_real_

cat("\n=== Paramètres estimés — SNU ===\n")
cat(sprintf("kg   = %.6f /j   (t½ tumeur = %.1f j)\n", best_pd["kg"], log(2)/best_pd["kg"]))
cat(sprintf("ki   = %.6f /j\n", best_pd["ki"]))
cat(sprintf("EC50 = %.0f µg/L = %.2f mg/kg (via V1)\n",
            best_pd["EC50"], best_pd["EC50"] * V1 / 1000))
cat("──────────────────────────────────────────────────────────\n")
cat(sprintf("ki/kg = %.2f  →  %s\n",
            ratio_ki_kg,
            ifelse(ratio_ki_kg > 1,
                   "régression tumorale atteignable (ki > kg)",
                   "ralentissement sans régression (ki ≤ kg)")))
if (!is.na(tsc_mgkg))
  cat(sprintf("TSC   = %.0f µg/L = %.2f mg/kg\n", tsc_conc, tsc_mgkg))
cat(sprintf("Objectif DEoptim = %.6f\n", fit_de$optim$bestval))
cat(sprintf("Objectif final   = %.6f\n", fit_local$objective))

# =============================================================================
# 9. SIMULATION FINALE
# =============================================================================

params_best <- c(kg = best_pd["kg"], ki = best_pd["ki"], EC50 = best_pd["EC50"])

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim  <- sim_ctrl_fn(best_pd["kg"], tv0_ctrl,  times_sim)
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
# 10. GRAPHIQUE
# =============================================================================

cols <- c("Vehicule"    = "#888888",
          "1.2 mg/kg"  = "#74ADD1",
          "5.6 mg/kg"  = "#4393C3",
          "11.8 mg/kg" = "#2166AC")

ymax <- max(df_obs$TV + df_obs$sem, na.rm = TRUE)

subtitle_txt <- paste0(
  "kg = ", round(best_pd["kg"], 4), " /j  |  ",
  "ki = ", round(best_pd["ki"], 4), " /j  |  ",
  "EC50 = ", round(best_pd["EC50"], 0), " µg/L",
  if (!is.na(tsc_mgkg)) paste0("  |  TSC = ", round(tsc_mgkg, 2), " mg/kg") else ""
)

p_pkpd <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(data  = df_sim,
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1) +
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
    title    = "Modèle PK/PD TGI — Fc-silent FGFR2-huBPA-LP1 — SNU",
    subtitle = subtitle_txt,
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm³)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 8, color = "grey50"),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

print(p_pkpd)
ggsave("scripts/plot_PKPD_TGI_SNU.png", p_pkpd,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_PKPD_TGI_SNU.png\n")

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

pkpd_results_SNU <- list(
  kg           = unname(best_pd["kg"]),
  ki           = unname(best_pd["ki"]),
  EC50_ugL     = unname(best_pd["EC50"]),
  EC50_mgkg    = unname(best_pd["EC50"]) * V1 / 1000,
  ratio_ki_kg  = ratio_ki_kg,
  TSC_ugL      = tsc_conc,
  TSC_mgkg     = tsc_mgkg,
  t_half_tumor = log(2) / unname(best_pd["kg"]),
  pk_fixed     = pk_fixed,
  objective    = fit_local$objective,
  convergence  = fit_local$convergence
)

save(pkpd_results_SNU, file = "scripts/resultats_PKPD_TGI_SNU.RData")
cat("Résultats → scripts/resultats_PKPD_TGI_SNU.RData\n")
