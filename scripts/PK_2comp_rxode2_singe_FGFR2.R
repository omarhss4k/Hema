# =============================================================================
# PK 2-compartiments — Fc-silent B/C huBPA-LP1 (FGFR2) — Singe
# rxode2 + nlminb  |  IV bolus unique  |  4 groupes : 3, 10, 20, 30 mg/kg
#
# Structure Excel (PK_singe_FGFR2.xlsx) :
#   4 blocs séparés par des lignes d'en-tête "Time (h)"
#   Bloc 1 : 3 mg/kg  (animaux 1001, 1002)
#   Bloc 2 : 10 mg/kg (animaux 2001, 2002)
#   Bloc 3 : 30 mg/kg (animaux 3002, 3101)
#   Bloc 4 : 20 mg/kg (animaux 4001, 4002)
#
# Unités : temps en heures | dose en µg/kg | concentration en ng/mL (= µg/L)
# =============================================================================

library(rxode2)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. LECTURE BRUTE DU FICHIER MULTI-BLOCS
# =============================================================================

raw <- suppressMessages(
  read_xlsx("PK_singe_FGFR2.xlsx",
            sheet      = "Sheet1",
            col_names  = FALSE,
            .name_repair = "minimal")
)

# Repérer toutes les lignes d'en-tête ("Time (h)")
header_rows <- which(raw[[1]] == "Time (h)")
cat("Lignes 'Time (h)' trouvées :", header_rows, "\n")
stopifnot("4 blocs attendus" = length(header_rows) == 4)

# Bornes de chaque bloc (données = lignes entre header+1 et prochain header-1)
block_ends <- c(header_rows[-1] - 1, nrow(raw))
blocks_range <- mapply(
  function(h, e) c(h + 1, e),
  header_rows, block_ends,
  SIMPLIFY = FALSE
)

# =============================================================================
# 2. EXTRACTION D'UN BLOC → data.frame (temps + animaux)
# =============================================================================

parse_block <- function(header_row, data_rows) {
  # Lire l'en-tête pour récupérer les IDs animaux
  hdr <- as.character(unlist(raw[header_row, ]))
  # Colonnes avec un ID animal (non "Time (h)", non NA)
  animal_cols <- which(!is.na(hdr) & hdr != "Time (h)")
  animal_ids  <- hdr[animal_cols]

  # Extraire les données de ce bloc
  blk <- raw[data_rows[1]:data_rows[2], ]

  # Colonne temps
  t_raw <- suppressWarnings(as.numeric(as.character(blk[[1]])))

  # Concentrations par animal (BLQ → NA)
  conc_list <- lapply(animal_cols, function(j) {
    x <- as.character(blk[[j]])
    x[toupper(trimws(x)) == "BLQ"] <- NA
    suppressWarnings(as.numeric(x))
  })

  df <- data.frame(time = t_raw)
  for (k in seq_along(animal_ids)) {
    df[[animal_ids[k]]] <- conc_list[[k]]
  }

  # Supprimer les lignes sans temps valide
  df[!is.na(df$time), ]
}

# Lire les 4 blocs
blk1 <- parse_block(header_rows[1], blocks_range[[1]])  # 3  mg/kg
blk2 <- parse_block(header_rows[2], blocks_range[[2]])  # 10 mg/kg
blk3 <- parse_block(header_rows[3], blocks_range[[3]])  # 30 mg/kg
blk4 <- parse_block(header_rows[4], blocks_range[[4]])  # 20 mg/kg

cat("\nBloc 1 (3  mg/kg) :", nrow(blk1), "lignes, colonnes:", names(blk1), "\n")
cat("Bloc 2 (10 mg/kg) :", nrow(blk2), "lignes, colonnes:", names(blk2), "\n")
cat("Bloc 3 (30 mg/kg) :", nrow(blk3), "lignes, colonnes:", names(blk3), "\n")
cat("Bloc 4 (20 mg/kg) :", nrow(blk4), "lignes, colonnes:", names(blk4), "\n")

# =============================================================================
# 3. MOYENNE GÉOMÉTRIQUE PAR GROUPE DE DOSE + NETTOYAGE
# =============================================================================

# Moyenne géométrique (NA ignorés) → robuste pour les PK log-normales
geom_mean <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  exp(mean(log(x)))
}

agg_group <- function(blk, dose_mgkg) {
  # Colonnes de concentration = tout sauf "time"
  conc_cols <- setdiff(names(blk), "time")
  conc_mat  <- as.matrix(blk[, conc_cols])
  avg       <- apply(conc_mat, 1, geom_mean)
  ok        <- !is.na(avg) & avg > 0 & blk$time > 0   # exclure t=0 et BLQ
  data.frame(
    t    = blk$time[ok],
    C    = avg[ok],
    dose = dose_mgkg * 1000,   # mg/kg → µg/kg
    Dose = paste0(dose_mgkg, " mg/kg")
  )
}

d3  <- agg_group(blk1,  3)
d10 <- agg_group(blk2, 10)
d30 <- agg_group(blk3, 30)
d20 <- agg_group(blk4, 20)

cat("\nPoints retenus par groupe :\n")
cat("  3  mg/kg :", nrow(d3),  "points\n")
cat("  10 mg/kg :", nrow(d10), "points\n")
cat("  20 mg/kg :", nrow(d20), "points\n")
cat("  30 mg/kg :", nrow(d30), "points\n")

# Jeu complet pour le graphique (brutes individuelles + BLQ)
make_obs_df <- function(blk, dose_mgkg) {
  conc_cols <- setdiff(names(blk), "time")
  conc_mat  <- as.matrix(blk[, conc_cols])
  avg       <- apply(conc_mat, 1, geom_mean)
  blq       <- apply(is.na(as.matrix(blk[, conc_cols])) |
                       apply(as.matrix(blk[, conc_cols]), 2,
                             function(x) toupper(trimws(as.character(x))) == "BLQ"),
                     1, all)
  data.frame(
    t    = blk$time,
    C    = avg,
    BLQ  = blq,
    Dose = paste0(dose_mgkg, " mg/kg")
  )
}

df_obs_all <- do.call(rbind, list(
  make_obs_df(blk1,  3),
  make_obs_df(blk2, 10),
  make_obs_df(blk4, 20),
  make_obs_df(blk3, 30)
))

# =============================================================================
# 4. MODÈLE rxode2 — 2 compartiments IV bolus
#
#   dA1/dt = -(CL/V1 + Q/V1)·A1 + (Q/V2)·A2
#   dA2/dt =  (Q/V1)·A1 - (Q/V2)·A2
#   C1 = A1/V1  (µg/kg / L/kg = µg/L = ng/mL)
# =============================================================================

mod2comp <- rxode2({
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
  d/dt(A2) <-  (Q/V1) * A1 - (Q/V2) * A2
})

sim_dose <- function(dose_ugkg, params, times) {
  ev <- eventTable()
  ev$add.dosing(dose = dose_ugkg, nbr.doses = 1, dosing.to = 1)
  ev$add.sampling(sort(unique(c(0, times))))

  out <- tryCatch(rxSolve(mod2comp, params, ev), error = function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, length(times)))
  approx(out$time, out$C1, xout = times, rule = 2)$y
}

# =============================================================================
# 5. VALEURS INITIALES — mAb IV chez le singe (typiques)
# =============================================================================

init_params <- c(
  CL = 0.005,   # L/h/kg  (~0.12 L/j/kg, typique IgG primate)
  V1 = 0.05,    # L/kg    (compartiment central ≈ volume plasmatique)
  V2 = 0.04,    # L/kg    (compartiment périphérique)
  Q  = 0.003    # L/h/kg  (clairance intercompartimentale)
)

cat("\n=== Valeurs initiales ===\n")
cat(sprintf("CL = %.4f L/h/kg\n", init_params["CL"]))
cat(sprintf("V1 = %.4f L/kg\n",   init_params["V1"]))
cat(sprintf("V2 = %.4f L/kg\n",   init_params["V2"]))
cat(sprintf("Q  = %.4f L/h/kg\n", init_params["Q"]))

# =============================================================================
# 6. FONCTION OBJECTIVE — résidus log², poids égaux, 4 groupes
# =============================================================================

objective <- function(logpar) {
  par        <- exp(logpar)
  names(par) <- c("CL", "V1", "V2", "Q")
  if (any(par <= 0)) return(1e10)

  p3  <- tryCatch(sim_dose(d3$dose[1],  par, d3$t),  error = function(e) NULL)
  p10 <- tryCatch(sim_dose(d10$dose[1], par, d10$t), error = function(e) NULL)
  p20 <- tryCatch(sim_dose(d20$dose[1], par, d20$t), error = function(e) NULL)
  p30 <- tryCatch(sim_dose(d30$dose[1], par, d30$t), error = function(e) NULL)

  if (is.null(p3)  || any(is.na(p3)  | p3  <= 0)) return(1e10)
  if (is.null(p10) || any(is.na(p10) | p10 <= 0)) return(1e10)
  if (is.null(p20) || any(is.na(p20) | p20 <= 0)) return(1e10)
  if (is.null(p30) || any(is.na(p30) | p30 <= 0)) return(1e10)

  mean((log(d3$C)  - log(p3))^2)  +
  mean((log(d10$C) - log(p10))^2) +
  mean((log(d20$C) - log(p20))^2) +
  mean((log(d30$C) - log(p30))^2)
}

# =============================================================================
# 7. OPTIMISATION — nlminb en log-espace
# =============================================================================

cat("\nOptimisation nlminb en cours...\n")

fit <- nlminb(
  start     = log(init_params),
  objective = objective,
  control   = list(eval.max = 3000, iter.max = 1500,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

best_par        <- exp(fit$par)
names(best_par) <- c("CL", "V1", "V2", "Q")

# =============================================================================
# 8. PARAMÈTRES DÉRIVÉS
# =============================================================================

k10 <- best_par["CL"] / best_par["V1"]
k12 <- best_par["Q"]  / best_par["V1"]
k21 <- best_par["Q"]  / best_par["V2"]
Vss <- best_par["V1"] + best_par["V2"]

sum_k <- k10 + k12 + k21
disc  <- sqrt((k10 + k12 - k21)^2 + 4 * k12 * k21)
alpha <- (sum_k + disc) / 2
beta  <- (sum_k - disc) / 2

t_half_alpha <- log(2) / alpha
t_half_beta  <- log(2) / beta

cat("\n=== Paramètres PK 2-compartiments — Singe (rxode2 + nlminb) ===\n")
cat("--- Paramètres macro ---\n")
cat(sprintf("CL  = %.6f L/h/kg   = %.4f L/j/kg\n",
            best_par["CL"], best_par["CL"] * 24))
cat(sprintf("V1  = %.5f L/kg    (compartiment central)\n",  best_par["V1"]))
cat(sprintf("V2  = %.5f L/kg    (compartiment périphérique)\n", best_par["V2"]))
cat(sprintf("Vss = %.5f L/kg    (volume de distribution à l'état stationnaire)\n", Vss))
cat(sprintf("Q   = %.6f L/h/kg  = %.4f L/j/kg\n",
            best_par["Q"], best_par["Q"] * 24))
cat("--- Constantes de vitesse ---\n")
cat(sprintf("k10 = %.6f /h  (élimination)\n",  k10))
cat(sprintf("k12 = %.6f /h  (C1 → C2)\n",      k12))
cat(sprintf("k21 = %.6f /h  (C2 → C1)\n",      k21))
cat("--- Demi-vies ---\n")
cat(sprintf("α   = %.6f /h\n",  alpha))
cat(sprintf("β   = %.7f /h\n",  beta))
cat(sprintf("t½α = %.2f h  (phase de distribution)\n", t_half_alpha))
cat(sprintf("t½β = %.1f h  = %.2f jours  (phase d'élimination)\n",
            t_half_beta, t_half_beta / 24))
cat("--- Ajustement ---\n")
cat(sprintf("Objectif final : %.6f  (somme résidus² log)\n", fit$objective))
cat(sprintf("Convergence    : %s (code %d)\n", fit$message, fit$convergence))

# =============================================================================
# 9. SIMULATION FINALE
# =============================================================================

t_max  <- max(d3$t, d10$t, d20$t, d30$t) * 1.05
t_sim  <- seq(0, t_max, by = 1)

sim_all <- do.call(rbind, lapply(
  list(c(d3$dose[1],   "3 mg/kg"),
       c(d10$dose[1],  "10 mg/kg"),
       c(d20$dose[1],  "20 mg/kg"),
       c(d30$dose[1],  "30 mg/kg")),
  function(x) data.frame(
    t    = t_sim,
    C    = sim_dose(as.numeric(x[1]), best_par, t_sim),
    Dose = x[2]
  )
))

niv <- c("3 mg/kg", "10 mg/kg", "20 mg/kg", "30 mg/kg")
sim_all$Dose      <- factor(sim_all$Dose,      levels = niv)
df_obs_all$Dose   <- factor(df_obs_all$Dose,   levels = niv)

# Données pour le graphique (uniquement points valides, t > 0)
df_pts <- df_obs_all[df_obs_all$t > 0 & !df_obs_all$BLQ & !is.na(df_obs_all$C), ]

# =============================================================================
# 10. GRAPHIQUE
# =============================================================================

cols <- c("3 mg/kg"  = "#74ADD1",
          "10 mg/kg" = "#4393C3",
          "20 mg/kg" = "#2166AC",
          "30 mg/kg" = "#053061")

p <- ggplot() +
  geom_line(data  = sim_all,
            aes(x = t, y = C, color = Dose),
            linewidth = 1) +
  geom_point(data = df_pts,
             aes(x = t, y = C, color = Dose),
             size = 2.5) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_color_manual(values = cols) +
  labs(
    title    = "PK 2-compartiments (rxode2) — Fc-silent B/C huBPA-LP1 (FGFR2) — Singe",
    subtitle = paste0(
      "V1 = ", round(best_par["V1"], 4), " L/kg",
      "   V2 = ", round(best_par["V2"], 4), " L/kg",
      "   CL = ", round(best_par["CL"] * 24, 4), " L/j/kg",
      "   t½α = ", round(t_half_alpha, 1), " h",
      "   t½β = ", round(t_half_beta / 24, 1), " j"
    ),
    x     = "Temps (heures)",
    y     = "Concentration (ng/mL, échelle log)",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 10, color = "grey50")
  )

print(p)
ggsave("scripts/plot_PK2comp_rxode2_singe_FGFR2.png",
       p, width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_PK2comp_rxode2_singe_FGFR2.png\n")

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

pk2comp_singe <- list(
  espece       = "singe",
  CL           = unname(best_par["CL"]),
  V1           = unname(best_par["V1"]),
  V2           = unname(best_par["V2"]),
  Q            = unname(best_par["Q"]),
  k10          = unname(k10),
  k12          = unname(k12),
  k21          = unname(k21),
  Vss          = unname(Vss),
  alpha        = unname(alpha),
  beta         = unname(beta),
  t_half_alpha = unname(t_half_alpha),
  t_half_beta  = unname(t_half_beta),
  objective    = fit$objective,
  convergence  = fit$convergence
)

save(pk2comp_singe,
     file = "scripts/resultats_PK2comp_rxode2_singe_FGFR2.RData")
cat("Résultats → scripts/resultats_PK2comp_rxode2_singe_FGFR2.RData\n")
