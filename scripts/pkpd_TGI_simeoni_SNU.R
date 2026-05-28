# =============================================================================
# Modèle PK/PD TGI — Simeoni (2004) — Fc-silent FGFR2-huBPA-LP1
# Cohorte : SNU (xénogreffe)
#
# PK  : 1 compartiment IV bolus — paramètres estimés depuis les données souris
#         dA1/dt = -ke · A1      A1 en µg/kg
#
# PD  : Modèle de Simeoni — compartiments de transit
#         dx1/dt = [g(w)/w - k2e·A1] · x1      k2e = k2/V1_souris
#         dx2/dt =  k2e·A1·x1 - k1·x2
#         dx3/dt =  k1·x2 - k1·x3
#         dx4/dt =  k1·x3 - k1·x4
#         w = x1+x2+x3+x4
#         g(w) = λ0·w / (1 + (λ0·w/λ1)^ψ)^(1/ψ)   ψ=20
#
# Paramètres estimés (5) : λ0, λ1, k1, k2e, ke
#
#   TSC_dose = λ0/k2e  (µg/kg — dose statique tumorale)
#   MTT      = 4/k1    (jours)
#   t½_drug  = ln2/ke  (jours)
#
# POURQUOI 1-cpt estimé plutôt que PK singe fixée ?
#   Avec la PK singe (V1 ~ 0.033 L/kg), Cmax(5.6 mg/kg) >> TSC, mais
#   Cmax(1.2 mg/kg) >> TSC aussi → le modèle ne pouvait pas expliquer
#   pourquoi 1.2 mg/kg ne fait pas régresser alors que 5.6 mg/kg oui.
#   En estimant ke et k2e depuis les données souris, la TSC_dose peut
#   se placer entre 1.2 et 5.6 mg/kg, ce qui est cohérent avec les données.
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
# 1. DONNÉES TUMORALES — SNU
#
# Colonnes attendues dans chaque bloc :
#   1 (A) = Temps (jours)
#   2 (B) = Vehicule
#   3 (C) = Group 4 — 1.2 mg/kg
#   4 (D) = Group 5 — 5.6 mg/kg
#   5 (E) = colonne vide (espacement)
#   6 (F) = Group 7 — 11.8 mg/kg
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

cat(sprintf("TV initiale (j0) : Ctrl=%.0f | 1.2mg=%.0f | 5.6mg=%.0f | 11.8mg=%.0f mm³\n",
            tv0_ctrl, tv0_d1p2, tv0_d5p6, tv0_d11p8))

# =============================================================================
# 2. ÉQUATIONS DU MODÈLE
# =============================================================================

PSI <- 20

simeoni_rhs <- function(t, state, parms) {
  A1 <- max(state["A1"], 0)
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)

  ke  <- parms["ke"];  k2e <- parms["k2e"]
  L0  <- parms["L0"];  L1  <- parms["L1"];  k1 <- parms["k1"]

  dA1 <- -ke * A1

  w  <- x1 + x2 + x3 + x4
  gw <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  dx1 <- (growth_rate - k2e * A1) * x1
  dx2 <- k2e * A1 * x1 - k1 * x2
  dx3 <- k1 * x2 - k1 * x3
  dx4 <- k1 * x3 - k1 * x4

  list(c(A1 = dA1, x1 = dx1, x2 = dx2, x3 = dx3, x4 = dx4))
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
# 3. FONCTIONS DE SIMULATION
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
      y     = c(A1 = dose_ugkg, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
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
# 4. FONCTION OBJECTIVE
#    Paramètres : L0, L1, k1, k2e, ke  (log-espace)
# =============================================================================

objective_simeoni <- function(logpar) {
  par <- exp(logpar)
  L0  <- unname(par[1]); L1  <- unname(par[2])
  k1  <- unname(par[3]); k2e <- unname(par[4]); ke <- unname(par[5])
  if (any(par <= 0)) return(1e12)

  params_all <- c(ke = ke, k2e = k2e, L0 = L0, L1 = L1, k1 = k1)

  pred_ctrl <- sim_ctrl_fn(L0, L1, tv0_ctrl, times_d[ok_ctrl])
  if (any(is.na(pred_ctrl))) return(1e12)
  pred_ctrl <- pmax(pred_ctrl, 0.1)

  pred_d1p2 <- sim_treated(1200,  tv0_d1p2,  params_all, times_d[ok_d1p2])
  if (any(is.na(pred_d1p2))) return(1e12)
  pred_d1p2 <- pmax(pred_d1p2, 0.1)

  pred_d5p6 <- sim_treated(5600,  tv0_d5p6,  params_all, times_d[ok_d5p6])
  if (any(is.na(pred_d5p6))) return(1e12)
  pred_d5p6 <- pmax(pred_d5p6, 0.1)

  pred_d11p8 <- sim_treated(11800, tv0_d11p8, params_all, times_d[ok_d11p8])
  if (any(is.na(pred_d11p8))) return(1e12)
  pred_d11p8 <- pmax(pred_d11p8, 0.1)

  sum(1/sem_ctrl[ok_ctrl]^2    * (log(tv_ctrl[ok_ctrl])    - log(pred_ctrl))^2,   na.rm = TRUE) +
  sum(1/sem_d1p2[ok_d1p2]^2   * (log(tv_d1p2[ok_d1p2])   - log(pred_d1p2))^2,  na.rm = TRUE) +
  sum(1/sem_d5p6[ok_d5p6]^2   * (log(tv_d5p6[ok_d5p6])   - log(pred_d5p6))^2,  na.rm = TRUE) +
  sum(1/sem_d11p8[ok_d11p8]^2 * (log(tv_d11p8[ok_d11p8]) - log(pred_d11p8))^2, na.rm = TRUE)
}

# =============================================================================
# 5. VALEURS INITIALES
# =============================================================================

lm_ctrl <- lm(log(tv_ctrl[ok_ctrl]) ~ times_d[ok_ctrl])
L0_init  <- max(coef(lm_ctrl)[2], 0.005)
L1_init  <- max(tv_ctrl[ok_ctrl], na.rm = TRUE) * L0_init * 50
k1_init  <- 4 / 7
ke_init  <- log(2) / 7          # t½ ~ 7 j pour un mAb chez la souris
k2e_init <- L0_init / 2500      # TSC_dose cible ~ 2500 µg/kg (entre 1.2 et 5.6 mg/kg)

init_pd <- c(L0 = L0_init, L1 = L1_init, k1 = k1_init, k2e = k2e_init, ke = ke_init)

cat("\n=== Valeurs initiales ===\n")
cat(sprintf("L0   = %.5f /j   (λ0)\n",                          init_pd["L0"]))
cat(sprintf("L1   = %.2f mm³/j\n",                               init_pd["L1"]))
cat(sprintf("k1   = %.4f /j   (MTT = %.1f j)\n",                 init_pd["k1"], 4/init_pd["k1"]))
cat(sprintf("k2e  = %.2e /j/(µg/kg)\n",                          init_pd["k2e"]))
cat(sprintf("ke   = %.4f /j   (t½ = %.1f j)\n",                  init_pd["ke"],  log(2)/init_pd["ke"]))
cat(sprintf("TSC_dose initiale ≈ %.0f µg/kg = %.1f mg/kg\n",
            init_pd["L0"]/init_pd["k2e"], init_pd["L0"]/init_pd["k2e"]/1000))

# =============================================================================
# 6. TEST DIAGNOSTIC
# =============================================================================

cat("\n=== TEST aux valeurs initiales ===\n")
params_test <- c(ke = ke_init, k2e = k2e_init, L0 = L0_init, L1 = L1_init, k1 = k1_init)

tc <- sim_ctrl_fn(L0_init, L1_init, tv0_ctrl, times_d[ok_ctrl])
cat(sprintf("Contrôle  : %s\n", if (any(is.na(tc))) "ECHEC" else paste(round(head(tc,4)), collapse=" | ")))

t1 <- sim_treated(1200,  tv0_d1p2,  params_test, times_d[ok_d1p2])
cat(sprintf("1.2 mg/kg : %s\n", if (any(is.na(t1))) "ECHEC" else paste(round(head(t1,4)), collapse=" | ")))

t5 <- sim_treated(5600,  tv0_d5p6,  params_test, times_d[ok_d5p6])
cat(sprintf("5.6 mg/kg : %s\n", if (any(is.na(t5))) "ECHEC" else paste(round(head(t5,4)), collapse=" | ")))

t11 <- sim_treated(11800, tv0_d11p8, params_test, times_d[ok_d11p8])
cat(sprintf("11.8mg/kg : %s\n", if (any(is.na(t11))) "ECHEC" else paste(round(head(t11,4)), collapse=" | ")))

obj0 <- objective_simeoni(log(unname(init_pd)))
cat(sprintf("Objectif initial : %.4f %s\n", obj0,
            if (obj0 >= 1e11) "— ECHEC" else "— OK"))

# =============================================================================
# 7. OPTIMISATION — DEoptim + affinage nlminb
#
#   Bornes (log-espace) :
#     L0   ∈ [0.005, 1.0]   /j
#     L1   ∈ [10, 1e6]      mm³/j
#     k1   ∈ [0.2, 4.0]     /j       (MTT ∈ [1, 20] j)
#     k2e  ∈ [1e-9, 1e-3]   /j/(µg/kg)
#     ke   ∈ [0.03, 2.0]    /j       (t½ ∈ 8h – 23j)
# =============================================================================

lower_log <- c(log(0.005), log(10),   log(0.2), log(1e-9), log(0.03))
upper_log <- c(log(1.0),   log(1e6),  log(4.0), log(1e-3), log(2.0))

cat("\nOptimisation DEoptim en cours...\n")
set.seed(42)
fit_de <- DEoptim(
  fn      = objective_simeoni,
  lower   = lower_log,
  upper   = upper_log,
  control = DEoptim.control(
    NP      = 150,    # 30 × nb_params (5)
    itermax = 800,
    F       = 0.8,
    CR      = 0.9,
    trace   = 100,
    reltol  = 1e-8,
    steptol = 200
  )
)

fit_local <- nlminb(
  start     = fit_de$optim$bestmem,
  objective = objective_simeoni,
  lower     = lower_log,
  upper     = upper_log,
  control   = list(eval.max = 3000, iter.max = 1000,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

best_pd        <- exp(fit_local$par)
names(best_pd) <- c("L0", "L1", "k1", "k2e", "ke")

# =============================================================================
# 8. PARAMÈTRES DÉRIVÉS
# =============================================================================

tsc_dose <- unname(best_pd["L0"] / best_pd["k2e"])   # µg/kg
tsc_mgkg <- tsc_dose / 1000                           # mg/kg
mtt      <- 4 / unname(best_pd["k1"])
t_half   <- log(2) / unname(best_pd["ke"])

cat("\n=== Paramètres estimés — SNU (PK 1-cpt souris) ===\n")
cat(sprintf("λ0 (L0) = %.6f /j\n",       best_pd["L0"]))
cat(sprintf("λ1 (L1) = %.2f mm³/j\n",    best_pd["L1"]))
cat(sprintf("k1      = %.5f /j\n",        best_pd["k1"]))
cat(sprintf("k2e     = %.3e /j/(µg/kg)\n", best_pd["k2e"]))
cat(sprintf("ke      = %.5f /j\n",        best_pd["ke"]))
cat("──────────────────────────────────────────────────\n")
cat(sprintf("TSC_dose = %.0f µg/kg = %.2f mg/kg\n", tsc_dose, tsc_mgkg))
cat(sprintf("MTT      = %.1f j\n", mtt))
cat(sprintf("t½_drug  = %.1f j\n", t_half))
cat(sprintf("Objectif DEoptim = %.6f\n", fit_de$optim$bestval))
cat(sprintf("Objectif final   = %.6f\n", fit_local$objective))

# =============================================================================
# 9. SIMULATION FINALE
# =============================================================================

params_best <- c(ke  = unname(best_pd["ke"]),
                 k2e = unname(best_pd["k2e"]),
                 L0  = unname(best_pd["L0"]),
                 L1  = unname(best_pd["L1"]),
                 k1  = unname(best_pd["k1"]))

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim  <- sim_ctrl_fn(unname(best_pd["L0"]), unname(best_pd["L1"]), tv0_ctrl,  times_sim)
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
  "λ0 = ", round(best_pd["L0"], 4), " /j  |  ",
  "k2e = ", formatC(best_pd["k2e"], digits = 3, format = "e"), " /j/(µg/kg)  |  ",
  "k1 = ", round(best_pd["k1"], 3), " /j  |  ",
  "ke = ", round(best_pd["ke"], 3), " /j  |  ",
  "TSC = ", round(tsc_mgkg, 2), " mg/kg  |  ",
  "MTT = ", round(mtt, 1), " j  |  ",
  "t½ = ", round(t_half, 1), " j"
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
    plot.subtitle   = element_text(size = 8, color = "grey50"),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

print(p_simeoni)
ggsave("scripts/plot_PKPD_simeoni_SNU.png", p_simeoni,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_PKPD_simeoni_SNU.png\n")

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

simeoni_results_SNU <- list(
  L0                = unname(best_pd["L0"]),
  L1                = unname(best_pd["L1"]),
  k1                = unname(best_pd["k1"]),
  k2e               = unname(best_pd["k2e"]),
  ke                = unname(best_pd["ke"]),
  TSC_dose_ugkg     = tsc_dose,
  TSC_mgkg          = tsc_mgkg,
  MTT_days          = mtt,
  t_half_drug_days  = t_half,
  objective_deoptim = fit_de$optim$bestval,
  objective_final   = fit_local$objective
)

save(simeoni_results_SNU, file = "scripts/resultats_PKPD_simeoni_SNU.RData")
cat("Résultats → scripts/resultats_PKPD_simeoni_SNU.RData\n")
