# =============================================================================
# Modèle PK/PD TGI — Fc-silent FGFR2-huBPA-LP1
#
# PK  : 2 compartiments IV bolus (paramètres fixés depuis fit précédent)
# PD  : dTV/dt = kg·TV − ke·[C1/(EC50+C1)]·TV
#
# Unités : temps en jours | dose en µg/kg | concentration en µg/L
#
# Schéma : Q2W × 4 (j0, j14, j28, j42)
# Groupes : Contrôle (Group 01) | 3 mg/kg (Group 04) | 10 mg/kg (Group 03)
# =============================================================================

library(rxode2)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. PARAMÈTRES PK FIXÉS (conversion h → j)
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2

pk_fixed <- c(
  CL = unname(pk2comp_rxode2$CL) * 24,   # L/h/kg → L/j/kg
  V1 = unname(pk2comp_rxode2$V1),         # L/kg
  V2 = unname(pk2comp_rxode2$V2),         # L/kg
  Q  = unname(pk2comp_rxode2$Q)  * 24    # L/h/kg → L/j/kg
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
tv_d3    <- get_row(blk_mean, "Group 04")
tv_d10   <- get_row(blk_mean, "Group 03")

# SEM plancher à 1 pour éviter poids infinis
sem_ctrl <- pmax(get_row(blk_sem, "Group 01"), 1)
sem_d3   <- pmax(get_row(blk_sem, "Group 04"), 1)
sem_d10  <- pmax(get_row(blk_sem, "Group 03"), 1)

# Volumes initiaux (j0)
get_tv0 <- function(tv) {
  v <- tv[times_d == 0]
  if (length(v) == 0 || is.na(v[1])) tv[!is.na(tv)][1] else v[1]
}
tv0_ctrl <- get_tv0(tv_ctrl)
tv0_d3   <- get_tv0(tv_d3)
tv0_d10  <- get_tv0(tv_d10)

# Masques : points non-NA
ok_ctrl <- !is.na(tv_ctrl)
ok_d3   <- !is.na(tv_d3)
ok_d10  <- !is.na(tv_d10)

cat(sprintf("\nTV initiale (j0) : Ctrl=%.0f | 3mg=%.0f | 10mg=%.0f mm³\n",
            tv0_ctrl, tv0_d3, tv0_d10))

# =============================================================================
# 3. MODÈLE rxode2 — PK/PD couplé
#
#   A1, A2 : quantités PK (µg/kg)
#   C1     : concentration centrale = A1/V1 (µg/L)
#   Inh    : fraction d'inhibition via modèle Emax (Hill n=1)
#   TV     : volume tumoral (mm³)
# =============================================================================

mod_pkpd <- rxode2({
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
  d/dt(A2) <-  (Q/V1) * A1 - (Q/V2) * A2
  Inh      <-  C1 / (EC50 + C1)
  d/dt(TV) <-  kg * TV - ke * Inh * TV
})

# =============================================================================
# 4. FONCTIONS DE SIMULATION
# =============================================================================

dose_days <- c(0, 14, 28, 42)   # Q2W × 4

# Contrôle : solution analytique (pas de drug)
sim_ctrl_fn <- function(kg, tv0, times_out) {
  tv0 * exp(kg * times_out)
}

# Groupe traité : premier bolus dans les conditions initiales, doses 2-4 via eventTable
# (évite l'ambiguïté rxode2 sur l'ordre dose/sampling à t=0)
sim_treated <- function(dose_ugkg, tv0, params, times_out) {
  ev <- eventTable()
  ev$add.dosing(dose = dose_ugkg, dosing.to = 1,
                nbr.doses = length(dose_days) - 1L,
                dosing.interval = 14,
                start.time = dose_days[2])          # j14, j28, j42
  ev$add.sampling(sort(unique(c(0, times_out))))

  out <- tryCatch(
    rxSolve(mod_pkpd, params, ev,
            inits = c(A1 = dose_ugkg, A2 = 0, TV = tv0)),   # bolus j0 en CI
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(times_out)))
  approx(out$time, out$TV, xout = times_out, rule = 2)$y
}

# =============================================================================
# 5. FONCTION OBJECTIVE
#    Paramètres libres : kg, ke, EC50 (optimisés en log-espace)
#    PK fixé
#    Pondération par 1/SEM² sur log(TV)
# =============================================================================

objective_pkpd <- function(logpar) {
  par        <- exp(logpar)
  kg  <- par[1]; ke <- par[2]; EC50 <- par[3]
  if (any(par <= 0)) return(1e12)

  params_all <- c(pk_fixed, kg = kg, ke = ke, EC50 = EC50)

  # — Contrôle —
  t_c       <- times_d[ok_ctrl]
  pred_ctrl <- sim_ctrl_fn(kg, tv0_ctrl, t_c)
  obs_ctrl  <- tv_ctrl[ok_ctrl]
  w_ctrl    <- 1 / sem_ctrl[ok_ctrl]^2
  if (any(pred_ctrl <= 0 | is.na(pred_ctrl))) return(1e12)

  # — 3 mg/kg —
  t_3     <- times_d[ok_d3]
  pred_d3 <- sim_treated(3000, tv0_d3, params_all, t_3)
  obs_d3  <- tv_d3[ok_d3]
  w_d3    <- 1 / sem_d3[ok_d3]^2
  if (any(is.na(pred_d3) | pred_d3 <= 0)) return(1e12)

  # — 10 mg/kg —
  t_10     <- times_d[ok_d10]
  pred_d10 <- sim_treated(10000, tv0_d10, params_all, t_10)
  obs_d10  <- tv_d10[ok_d10]
  w_d10    <- 1 / sem_d10[ok_d10]^2
  if (any(is.na(pred_d10) | pred_d10 <= 0)) return(1e12)

  sum(w_ctrl * (log(obs_ctrl) - log(pred_ctrl))^2, na.rm = TRUE) +
  sum(w_d3   * (log(obs_d3)   - log(pred_d3))^2,   na.rm = TRUE) +
  sum(w_d10  * (log(obs_d10)  - log(pred_d10))^2,  na.rm = TRUE)
}

# =============================================================================
# 6. VALEURS INITIALES PD
# =============================================================================

# kg : régression log-linéaire sur le contrôle
lm_ctrl <- lm(log(tv_ctrl[ok_ctrl]) ~ times_d[ok_ctrl])
kg_init <- max(coef(lm_ctrl)[2], 0.005)

# EC50 : 10 % du Cmax estimé à 10 mg/kg (= dose/V1)
# as.numeric() retire le nom hérité de pk_fixed["V1"] pour éviter
# que c(..., EC50 = ...) reçoive un élément nommé "V1" au lieu de "EC50"
EC50_init <- 10000 / as.numeric(pk_fixed["V1"]) * 0.10

# ke : doit être > kg pour qu'une régression soit possible
ke_init <- kg_init * 3

init_pd <- c(kg = kg_init, ke = ke_init, EC50 = EC50_init)

cat("\n=== Valeurs initiales PD ===\n")
cat(sprintf("kg   = %.5f /j  (t½ = %.1f j)\n", init_pd["kg"], log(2)/init_pd["kg"]))
cat(sprintf("ke   = %.5f /j\n",                  init_pd["ke"]))
cat(sprintf("EC50 = %.1f µg/L\n",                init_pd["EC50"]))

# =============================================================================
# 7. OPTIMISATION — nlminb (log-espace)
# =============================================================================

cat("\nOptimisation nlminb en cours...\n")

fit_pd <- nlminb(
  start     = log(init_pd),
  objective = objective_pkpd,
  control   = list(eval.max = 5000, iter.max = 2000,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

best_pd        <- exp(fit_pd$par)
names(best_pd) <- c("kg", "ke", "EC50")

# =============================================================================
# 8. PARAMÈTRES DÉRIVÉS
# =============================================================================

# ke/kg > 1 → inhibition maximale > croissance → régression possible
ratio_ke_kg <- unname(best_pd["ke"] / best_pd["kg"])

# Tumor Static Concentration (TSC) : C tel que kg = ke·C/(EC50+C)
# → TSC = kg·EC50 / (ke − kg)  [valide seulement si ke > kg]
tsc <- if (best_pd["ke"] > best_pd["kg"])
  unname(best_pd["kg"] * best_pd["EC50"] / (best_pd["ke"] - best_pd["kg"]))
else NA_real_

cat("\n=== Paramètres PD estimés ===\n")
cat(sprintf("kg   = %.6f /j   (t½ tumeur sans traitement = %.1f j)\n",
            best_pd["kg"], log(2) / best_pd["kg"]))
cat(sprintf("ke   = %.6f /j\n", best_pd["ke"]))
cat(sprintf("EC50 = %.2f µg/L\n", best_pd["EC50"]))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("ke/kg = %.2f  →  %s\n",
            ratio_ke_kg,
            ifelse(ratio_ke_kg > 1,
                   "régression tumorale atteignable (ke > kg)",
                   "ralentissement sans régression (ke ≤ kg)")))
if (!is.na(tsc))
  cat(sprintf("TSC   = %.1f µg/L  (concentration statique tumorale)\n", tsc))
cat(sprintf("Objectif final = %.6f\n", fit_pd$objective))
cat(sprintf("Convergence    : %s (code %d)\n",
            fit_pd$message, fit_pd$convergence))

# =============================================================================
# 9. SIMULATION FINALE
# =============================================================================

params_best <- c(pk_fixed,
                 kg   = best_pd["kg"],
                 ke   = best_pd["ke"],
                 EC50 = best_pd["EC50"])

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim <- sim_ctrl_fn(best_pd["kg"], tv0_ctrl, times_sim)
pred_d3_sim   <- sim_treated(3000,  tv0_d3,  params_best, times_sim)
pred_d10_sim  <- sim_treated(10000, tv0_d10, params_best, times_sim)

df_sim <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_sim, Groupe = "Contrôle"),
  data.frame(jour = times_sim, TV = pred_d3_sim,   Groupe = "3 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d10_sim,  Groupe = "10 mg/kg")
)

df_obs <- rbind(
  data.frame(jour = times_d[ok_ctrl], TV = tv_ctrl[ok_ctrl],
             sem  = sem_ctrl[ok_ctrl],  Groupe = "Contrôle"),
  data.frame(jour = times_d[ok_d3],   TV = tv_d3[ok_d3],
             sem  = sem_d3[ok_d3],      Groupe = "3 mg/kg"),
  data.frame(jour = times_d[ok_d10],  TV = tv_d10[ok_d10],
             sem  = sem_d10[ok_d10],    Groupe = "10 mg/kg")
)

niv <- c("Contrôle", "3 mg/kg", "10 mg/kg")
df_sim$Groupe <- factor(df_sim$Groupe, levels = niv)
df_obs$Groupe <- factor(df_obs$Groupe, levels = niv)

# =============================================================================
# 10. GRAPHIQUE
# =============================================================================

cols <- c("Contrôle" = "#888888", "3 mg/kg" = "#4393C3", "10 mg/kg" = "#2166AC")
ymax <- max(df_obs$TV + df_obs$sem, na.rm = TRUE)

subtitle_txt <- paste0(
  "kg = ", round(best_pd["kg"], 4),  " /j  |  ",
  "ke = ", round(best_pd["ke"], 4),  " /j  |  ",
  "EC50 = ", round(best_pd["EC50"], 0), " µg/L  |  ",
  "ke/kg = ", round(ratio_ke_kg, 2),
  if (!is.na(tsc)) paste0("  |  TSC = ", round(tsc, 0), " µg/L") else ""
)

p_pkpd <- ggplot() +
  geom_vline(xintercept = dose_days, linetype = "dashed",
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
  annotate("point", x = dose_days, y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text",  x = max(dose_days) + 1.5, y = -ymax * 0.06,
           label = "= Traitement", hjust = 0, size = 3.2, color = "#CC0000") +
  scale_x_continuous(breaks = seq(0, max(times_d, na.rm = TRUE), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "Modèle PK/PD TGI — Fc-silent FGFR2-huBPA-LP1",
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

print(p_pkpd)
ggsave("scripts/plot_PKPD_TGI_FGFR2.png", p_pkpd,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_PKPD_TGI_FGFR2.png\n")

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

pkpd_results <- list(
  pk_fixed    = pk_fixed,
  kg          = unname(best_pd["kg"]),
  ke          = unname(best_pd["ke"]),
  EC50        = unname(best_pd["EC50"]),
  ratio_ke_kg = ratio_ke_kg,
  TSC_ugL     = tsc,
  t_half_tumor = log(2) / unname(best_pd["kg"]),
  objective   = fit_pd$objective,
  convergence = fit_pd$convergence
)

save(pkpd_results, file = "scripts/resultats_PKPD_TGI_FGFR2.RData")
cat("Résultats → scripts/resultats_PKPD_TGI_FGFR2.RData\n")
