############################################################
# calibrate_neut.R — Calibration 1D : Slope_CMP
# Stratégie :
#   Ajuster Slope_CMP pour que le nadir de Neut (cycle 1 et 2)
#   corresponde aux données digitalisées Fornari 2019 Figure 4.
#
# Cibles (données Figure 4 digitalisées) :
#   nadir1_Neut ≈ 2.52 × 10⁹/L  (j12.4, cycle 1)
#   nadir2_Neut ≈ 2.32 × 10⁹/L  (j34.4, cycle 2)
#
# Paramètre : Slope_CMP contrôle la cytotoxicité sur les CMP
#   → affecte directement les neutrophiles (et monocytes secondairement)
#   → n'affecte pas les plaquettes (MEP) ni les érythrocytes
############################################################

if (interactive()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("pkpd_model_FORNARI.R")
source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")

# ── Dose Calvert (AUC=5, GFR=125 mL/min — patient ref. Fornari S11) ──
AUC_target   <- 5
GFR_mLmin    <- 125
dose_calvert <- AUC_target * (GFR_mLmin + 25)
times_hu     <- seq(0, 63 * 24, by = 1)

# ── Fonction de simulation silencieuse ──────────────────────────────
run_sim_neut <- function(slope_cmp) {
  p <- init_pars
  p$Slope_CMP <- slope_cmp
  p$rate_fun  <- make_repeated_infusion(
    dose_mg    = dose_calvert,
    Tinfu_h    = 1,
    interval_h = 21 * 24,
    n_cycles   = 2
  )
  invisible(capture.output(
    sim <- simulate_all(times_hu, p, init_state)
  ))
  # Cycle 1 : j7–j21  (168 h – 504 h)
  c1 <- sim[sim$time >= 168  & sim$time <= 504,  ]
  # Cycle 2 : j28–j42 (672 h – 1008 h)
  c2 <- sim[sim$time >= 672  & sim$time <= 1008, ]
  n1 <- min(c1$Neut, na.rm = TRUE)
  n2 <- min(c2$Neut, na.rm = TRUE)
  d1 <- c1$days[which.min(c1$Neut)]
  d2 <- c2$days[which.min(c2$Neut)]
  list(nadir1 = n1, nadir2 = n2, day1 = d1, day2 = d2)
}

# ── Bornes Slope_CMP ──────────────────────────────────────────────
sc_min <- 0.30   # en dessous : nadir trop élevé
sc_max <- 3.00   # au dessus  : kill dominant (Slope_CMP × Damage > 1)

# ── Tolérances ±15 % autour des valeurs digitalisées ──────────────
t1_lo <- 2.52 * 0.85 ; t1_hi <- 2.52 * 1.15   # [2.14, 2.90]
t2_lo <- 2.32 * 0.85 ; t2_hi <- 2.32 * 1.15   # [1.97, 2.67]

# ── Point de départ : valeur courante ─────────────────────────────
sc <- max(sc_min, min(sc_max, init_pars$Slope_CMP))

cat(sprintf("\n╔══════════════════════════════════════════════════╗\n"))
cat(sprintf("║  CALIBRATION NEUTROPHILES — Slope_CMP            ║\n"))
cat(sprintf("╚══════════════════════════════════════════════════╝\n"))
cat(sprintf("  Objectif nadir1 = 2.52  plage [2.14 – 2.90]  (donnée j12.4)\n"))
cat(sprintf("  Objectif nadir2 = 2.32  plage [1.97 – 2.67]  (donnée j34.4)\n\n"))
cat(sprintf("  Départ : Slope_CMP = %.3f\n\n", sc))

MAX_ITER <- 80
ok <- FALSE
n1 <- NA ; n2 <- NA

for (iter in 1:MAX_ITER) {

  res <- run_sim_neut(sc)
  n1  <- res$nadir1 ; n2 <- res$nadir2

  ok1 <- n1 >= t1_lo & n1 <= t1_hi
  ok2 <- n2 >= t2_lo & n2 <= t2_hi

  cat(sprintf("Iter %3d | sc=%.3f | n1=%5.3f [%s] n2=%5.3f [%s]",
              iter, sc,
              n1, ifelse(ok1, "OK", "!!"),
              n2, ifelse(ok2, "OK", "!!")))

  if (ok1 && ok2) {
    cat("  *** CONVERGÉ ***\n")
    ok <- TRUE
    break
  }

  new_sc <- sc

  # ── Priorité 1 : nadir1 hors plage → ajuster Slope_CMP ──────────
  # nadir trop haut → augmenter Slope_CMP (plus de kill)
  # nadir trop bas  → diminuer Slope_CMP
  if (!ok1) {
    if (n1 > t1_hi) new_sc <- sc + 0.05
    else            new_sc <- sc - 0.05
  }

  # ── Priorité 2 : nadir1 OK, nadir2 hors plage ────────────────────
  # Le nadir du cycle 2 dépend aussi de Slope_CMP (en l'absence
  # d'un paramètre de récupération distinct). On ajuste finement.
  if (ok1 && !ok2) {
    if (n2 > t2_hi) new_sc <- sc + 0.02   # ajustement plus fin
    else            new_sc <- sc - 0.02
  }

  # ── Bornes ────────────────────────────────────────────────────────
  new_sc <- max(sc_min, min(sc_max, new_sc))

  if (abs(new_sc - sc) < 1e-9) {
    cat("  [paramètre bloqué aux bornes, arrêt]\n")
    break
  }

  cat(sprintf("  → sc %.3f→%.3f\n", sc, new_sc))
  sc <- new_sc
}

# ── Résumé final ──────────────────────────────────────────────────
cat(sprintf("\n══════════════════════════════════════════════════\n"))
cat(sprintf("  RÉSULTATS FINAUX\n"))
cat(sprintf("══════════════════════════════════════════════════\n"))
cat(sprintf("  Slope_CMP   = %.3f\n", sc))
cat(sprintf("  Neut_nadir1 = %.3f  (cible 2.52, plage 2.14–2.90) %s\n",
            n1, ifelse(n1 >= t1_lo & n1 <= t1_hi, "✓", "✗")))
cat(sprintf("  Neut_nadir2 = %.3f  (cible 2.32, plage 1.97–2.67) %s\n",
            n2, ifelse(n2 >= t2_lo & n2 <= t2_hi, "✓", "✗")))
cat(sprintf("  Convergence : %s\n", ifelse(ok, "OUI ✓", "NON ✗")))
cat(sprintf("══════════════════════════════════════════════════\n\n"))

# ── Mise à jour de parameters_human.R si convergé ──────────────────
if (ok) {
  lines <- readLines("parameters_human.R")
  lines <- gsub(
    "^(init_pars\\$Slope_CMP\\s*<-).*",
    sprintf("\\1 %.3f  # calibration auto : nadir_Neut1~%.2f nadir_Neut2~%.2f", sc, n1, n2),
    lines, perl = TRUE
  )
  writeLines(lines, "parameters_human.R")
  cat(sprintf("  ✓ parameters_human.R mis à jour : Slope_CMP = %.3f\n\n", sc))
} else {
  cat("  ⚠ Pas de mise à jour (convergence non atteinte)\n\n")
  cat(sprintf("  Valeur courante : Slope_CMP = %.3f\n", sc))
  cat(sprintf("  → Derniers nadirs : n1=%.3f  n2=%.3f\n\n", n1, n2))
}
