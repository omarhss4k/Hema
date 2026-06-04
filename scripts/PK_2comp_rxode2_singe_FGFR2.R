# =============================================================================
# PK 2-compartiments — hBPA-LP1 (FGFR2) — Singe
# rxode2 + nlminb  |  IV bolus unique
#
# Structure Excel (PK_singe_FGFR2.xlsx) :
#   Multi-blocs verticaux séparés par des lignes "Time (h)"
#   Col 1 = temps (h), col 2 = concentration (ng/mL) — 1 colonne par groupe
#   Temps : 0.5, 4, 24, 48, 168, 336, 504, 648 h
#   Groupes : 1, 3, 3.5, 5, 10, 11.8 mg/kg (ordre libre dans le fichier)
#
# Unités : temps en heures | dose en µg/kg | concentration en ng/mL (= µg/L)
# =============================================================================

library(rxode2)
library(ggplot2)
library(readxl)

options(encoding = "UTF-8")
if (.Platform$OS.type == "unix") Sys.setlocale("LC_ALL", "C.UTF-8")

# =============================================================================
# 1. LECTURE BRUTE — multi-blocs verticaux
# =============================================================================

raw <- suppressMessages(
  read_xlsx("PK_singe_FGFR2.xlsx",
            sheet        = "Sheet1",
            col_names    = FALSE,
            .name_repair = "minimal")
)

# Repérer les lignes d'en-tête "Time (h)" dans la colonne 1
header_rows <- grep("Time.*\\(h\\)", trimws(as.character(raw[[1]])),
                    ignore.case = TRUE, perl = TRUE)
cat("Blocs 'Time (h)' trouvés :", length(header_rows), "→", header_rows, "\n")
stopifnot("Au moins 1 bloc attendu" = length(header_rows) >= 1)

# Bornes de chaque bloc
block_ends   <- c(header_rows[-1] - 1, nrow(raw))
blocks_range <- mapply(function(h, e) c(h + 1, e),
                       header_rows, block_ends, SIMPLIFY = FALSE)

# =============================================================================
# 2. PARSING — 1 colonne temps + 1 colonne concentration par bloc
# =============================================================================

extract_dose_mgkg <- function(s) {
  m <- regmatches(s, regexpr("[0-9]+(?:\\.[0-9]+)?(?=\\s*mg/kg)", s, perl = TRUE))
  if (length(m) == 0) return(NA_real_)
  as.numeric(m)
}

to_num <- function(x) {
  x <- trimws(as.character(x))
  x[toupper(x) == "BLQ" | x == ""] <- NA
  suppressWarnings(as.numeric(x))
}

parse_bloc <- function(header_row, data_rows) {
  hdr       <- trimws(as.character(unlist(raw[header_row, ])))
  dose_mgkg <- extract_dose_mgkg(hdr[2])
  if (is.na(dose_mgkg)) return(NULL)

  blk  <- raw[data_rows[1]:data_rows[2], ]
  t    <- suppressWarnings(as.numeric(as.character(blk[[1]])))
  conc <- to_num(blk[[2]])

  ok_fit <- !is.na(t) & t > 0 & !is.na(conc) & conc > 0
  ok_obs <- !is.na(t)

  list(
    fit = data.frame(t = t[ok_fit], C = conc[ok_fit],
                     dose_ugkg = dose_mgkg * 1000,
                     Dose = paste0(dose_mgkg, " mg/kg")),
    obs = data.frame(t = t[ok_obs], C = conc[ok_obs],
                     BLQ = is.na(conc[ok_obs]),
                     Dose = paste0(dose_mgkg, " mg/kg"))
  )
}

blocs <- lapply(seq_along(header_rows), function(i)
  parse_bloc(header_rows[i], blocks_range[[i]]))
blocs <- Filter(Negate(is.null), blocs)

# =============================================================================
# 3. NETTOYAGE — tri par dose, vérification
# =============================================================================

dose_vals  <- sapply(blocs, function(b) b$fit$dose_ugkg[1])
blocs      <- blocs[order(dose_vals)]

fit_groups <- lapply(blocs, function(b) b$fit)
df_obs_all <- do.call(rbind, lapply(blocs, function(b) b$obs))

cat("\nPoints retenus par groupe :\n")
for (g in fit_groups)
  cat(sprintf("  %-12s : %d points | dose = %g µg/kg\n",
              g$Dose[1], nrow(g), g$dose_ugkg[1]))

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
  CL = 0.001,   # L/h/kg
  V1 = 0.033,   # L/kg
  V2 = 0.025,   # L/kg
  Q  = 0.002    # L/h/kg
)

cat("\n=== Valeurs initiales ===\n")
cat(sprintf("CL = %.4f L/h/kg\n", init_params["CL"]))
cat(sprintf("V1 = %.4f L/kg\n",   init_params["V1"]))
cat(sprintf("V2 = %.4f L/kg\n",   init_params["V2"]))
cat(sprintf("Q  = %.4f L/h/kg\n", init_params["Q"]))

# =============================================================================
# 6. FONCTION OBJECTIVE — résidus log², poids égaux, tous groupes
# =============================================================================

objective <- function(logpar) {
  par        <- exp(logpar)
  names(par) <- c("CL", "V1", "V2", "Q")
  if (any(par <= 0)) return(1e10)

  resid <- sapply(fit_groups, function(g) {
    pred <- tryCatch(sim_dose(g$dose_ugkg[1], par, g$t), error = function(e) NULL)
    if (is.null(pred) || any(is.na(pred) | pred <= 0)) return(NA_real_)
    mean((log(g$C) - log(pred))^2)
  })

  if (any(is.na(resid))) return(1e10)
  sum(resid)
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

t_max <- max(sapply(fit_groups, function(g) max(g$t))) * 1.05
t_sim <- seq(0, t_max, by = 1)
niv   <- sapply(blocs, function(b) b$fit$Dose[1])   # trié par dose

sim_all <- do.call(rbind, lapply(blocs, function(b) {
  data.frame(
    t    = t_sim,
    C    = sim_dose(b$fit$dose_ugkg[1], best_par, t_sim),
    Dose = b$fit$Dose[1]
  )
}))

sim_all$Dose    <- factor(sim_all$Dose,    levels = niv)
df_obs_all$Dose <- factor(df_obs_all$Dose, levels = niv)

df_pts <- df_obs_all[df_obs_all$t > 0 & !df_obs_all$BLQ & !is.na(df_obs_all$C), ]

# =============================================================================
# 10. GRAPHIQUE
# =============================================================================

n_grp <- length(niv)
blues <- colorRampPalette(c("#C6DBEF", "#08306B"))(n_grp)
cols  <- setNames(blues, niv)

p <- ggplot() +
  geom_line(data  = sim_all,
            aes(x = t, y = C, color = Dose),
            linewidth = 1) +
  geom_point(data = df_pts,
             aes(x = t, y = C, color = Dose),
             size = 2.5) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
  scale_color_manual(values = cols) +
  labs(
    title    = "PK 2-compartiments (rxode2) — hBPA-LP1 (FGFR2) — Singe",
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
