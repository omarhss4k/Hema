############################################################
# calibrate_plt.R — Calibration 2D : delta_Plt + gamma_prolTrans
# Stratégie :
#   Phase 1 : delta_Plt pour mettre nadir1 dans [144.5, 195.5]
#   Phase 2 : gamma_prolTrans pour mettre nadir2 dans [126.9, 171.5]
#   (réajustement delta_Plt si nadir1 sort de sa plage)
#
# Cibles tirées des données digitalisées Fornari 2019 Figure 4 :
#   nadir1 ≈ 170 × 10⁹/L (Plt, j14.5, cycle 1)
#   nadir2 ≈ 149 × 10⁹/L (Plt, j35.7, cycle 2)
# Paramètres CL cohérents avec AUC=5 (GFR=125, CL=9.0 L/h)
############################################################

if (interactive()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("pkpd_model_FORNARI.R")
source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")

# ── Dose Calvert (AUC=5, GFR=125 mL/min — normal renal function, Fornari 2019 S11) ──
AUC_target   <- 5
GFR_mLmin    <- 125
dose_calvert <- AUC_target * (GFR_mLmin + 25)
times_hu     <- seq(0, 63 * 24, by = 1)

# ── Fonction de simulation silencieuse (supprime les cat() de simulate_all) ──
run_sim_quiet <- function(delta_plt, gamma_pt) {
  p <- init_pars
  p$delta_Plt       <- delta_plt
  p$gamma_prolTrans <- gamma_pt
  p$rate_fun <- make_repeated_infusion(
    dose_mg    = dose_calvert,
    Tinfu_h    = 1,
    interval_h = 21 * 24,
    n_cycles   = 2
  )
  invisible(capture.output(
    sim <- simulate_all(times_hu, p, init_state)
  ))
  c1 <- sim[sim$time >= 168  & sim$time <= 504,  ]
  c2 <- sim[sim$time >= 672  & sim$time <= 1008, ]
  n1 <- min(c1$Plt, na.rm=TRUE)
  n2 <- min(c2$Plt, na.rm=TRUE)
  d1 <- c1$days[which.min(c1$Plt)]
  d2 <- c2$days[which.min(c2$Plt)]
  list(nadir1=n1, nadir2=n2, day1=d1, day2=d2)
}

# ── Bornes ──
dp_min  <- 0.40 ; dp_max  <- 2.00   # élargi pour couvrir l'espace après fix CL
gpt_min <- 0.30 ; gpt_max <- 0.80

# ── Tolérances (cibles issues des données Figure 4 digitalisées) ──
t1_lo <- 144.5 ; t1_hi <- 195.5   # 170 ± 15%  (nadir cycle 1 observé : 169.7)
t2_lo <- 126.9 ; t2_hi <- 171.5   # 149 ± 15%  (nadir cycle 2 observé : 149.4)

# ── Paramètres initiaux ──
dp  <- max(dp_min, min(dp_max, init_pars$delta_Plt))
gpt <- max(gpt_min, min(gpt_max, init_pars$gamma_prolTrans))

cat(sprintf("\n╔══════════════════════════════════════════════════╗\n"))
cat(sprintf("║  CALIBRATION PLAQUETTES — STRATÉGIE 2D          ║\n"))
cat(sprintf("╚══════════════════════════════════════════════════╝\n"))
cat(sprintf("  Objectif nadir1 = 170  plage [144.5 – 195.5]  (donnée j14.5)\n"))
cat(sprintf("  Objectif nadir2 = 149  plage [126.9 – 171.5]  (donnée j35.7)\n\n"))
cat(sprintf("  Départ : delta_Plt=%.2f  gamma_prolTrans=%.2f\n\n", dp, gpt))

MAX_ITER <- 100
ok <- FALSE
n1 <- NA ; n2 <- NA

for (iter in 1:MAX_ITER) {

  res <- run_sim_quiet(dp, gpt)
  n1 <- res$nadir1 ; n2 <- res$nadir2

  ok1 <- n1 >= t1_lo & n1 <= t1_hi
  ok2 <- n2 >= t2_lo & n2 <= t2_hi

  cat(sprintf("Iter %3d | dp=%.2f gpt=%.2f | n1=%6.1f [%s] n2=%6.1f [%s]",
              iter, dp, gpt,
              n1, ifelse(ok1, "OK", "!!"),
              n2, ifelse(ok2, "OK", "!!")))

  if (ok1 && ok2) {
    cat("  *** CONVERGÉ ***\n")
    ok <- TRUE
    break
  }

  new_dp  <- dp
  new_gpt <- gpt

  # ── Priorité 1 : nadir1 hors plage → ajuster delta_Plt ──
  if (!ok1) {
    if (n1 > t1_hi) new_dp <- dp + 0.05   # trop haut → plus de kill
    else            new_dp <- dp - 0.05   # trop bas  → moins de kill
  }

  # ── Priorité 2 : nadir2 hors plage avec nadir1 OK ──
  if (ok1 && !ok2) {
    if (n2 > t2_hi) {
      if (gpt > gpt_min) {
        # diminuer gamma_prolTrans : récupération plus lente → cycle 2 plus profond
        new_gpt <- gpt - 0.05
      } else {
        # gamma_prolTrans à son minimum → augmenter delta_Plt pour creuser les 2 nadirs
        # (l'asymétrie déjà créée par gpt=min fait que nadir2 descend plus vite que nadir1)
        new_dp <- dp + 0.05
      }
    } else {
      # nadir2 trop bas → augmenter gamma_prolTrans
      if (gpt < gpt_max) new_gpt <- gpt + 0.05
      else               new_dp  <- dp - 0.05
    }
  }

  # ── Cas mixte : les deux hors plage dans la même direction ──
  if (!ok1 && !ok2) {
    # On ajuste d'abord delta_Plt (déjà fait ci-dessus)
    # Si on est à la borne de delta_Plt, on essaie gamma_prolTrans
    if (new_dp == dp) {  # delta_Plt bloqué par ses bornes
      if (n2 > t2_hi) new_gpt <- gpt - 0.05
      else            new_gpt <- gpt + 0.05
    }
  }

  # ── Appliquer les bornes ──
  new_dp  <- max(dp_min,  min(dp_max,  new_dp))
  new_gpt <- max(gpt_min, min(gpt_max, new_gpt))

  if (abs(new_dp - dp) < 1e-9 && abs(new_gpt - gpt) < 1e-9) {
    cat("  [paramètres bloqués aux bornes, arrêt]\n")
    break
  }

  cat(sprintf("  → dp %.2f→%.2f  gpt %.2f→%.2f\n", dp, new_dp, gpt, new_gpt))

  dp  <- new_dp
  gpt <- new_gpt
}

# ── Résumé final ──
cat(sprintf("\n══════════════════════════════════════════════════\n"))
cat(sprintf("  RÉSULTATS FINAUX\n"))
cat(sprintf("══════════════════════════════════════════════════\n"))
cat(sprintf("  delta_Plt       = %.2f\n", dp))
cat(sprintf("  gamma_prolTrans = %.2f\n", gpt))
cat(sprintf("  Plt_nadir1      = %.1f  (cible 170, plage 144.5–195.5) %s\n",
            n1, ifelse(n1>=t1_lo & n1<=t1_hi, "✓", "✗")))
cat(sprintf("  Plt_nadir2      = %.1f  (cible 149, plage 126.9–171.5) %s\n",
            n2, ifelse(n2>=t2_lo & n2<=t2_hi, "✓", "✗")))
cat(sprintf("  Asymétrie       = %.1f  (nadir1 - nadir2)\n", n1-n2))
cat(sprintf("  Convergence     : %s\n", ifelse(ok, "OUI ✓", "NON ✗")))
cat(sprintf("══════════════════════════════════════════════════\n\n"))

# ── Mise à jour de parameters_human.R si convergé ──
if (ok) {
  lines <- readLines("parameters_human.R")
  lines <- gsub(
    "^(init_pars\\$delta_Plt\\s*<-).*",
    sprintf("\\1 %.2f  # calibration auto : nadir1~%.0f nadir2~%.0f", dp, n1, n2),
    lines, perl=TRUE
  )
  lines <- gsub(
    "^(init_pars\\$gamma_prolTrans\\s*<-).*",
    sprintf("\\1 %.2f", gpt),
    lines, perl=TRUE
  )
  writeLines(lines, "parameters_human.R")
  cat(sprintf("  ✓ parameters_human.R mis à jour :\n"))
  cat(sprintf("      delta_Plt       = %.2f\n", dp))
  cat(sprintf("      gamma_prolTrans = %.2f\n\n", gpt))
} else {
  cat("  ⚠ Pas de mise à jour (convergence non atteinte)\n\n")
}
