############################################################
# calibrate_slope_mep.R
# Scan Slope_MEP pour T-DXd rat
#
# Objectif :
#   MEP@20 mg/kg < 10% (sous seuil histopathologique n=4)
#   MEP@60 mg/kg > 20% (au-dessus du seuil de détection)
#   FDA BLA 761139 : érythroblastes↓ seulement à ≥60 mg/kg
############################################################
library(deSolve)

source("../scripts/parameters_rat.R")
source("../scripts/parameters_FORNARI_CORRECT.R")
source("parameters_tdxd_rat.R")
source("pkpd_tdxd_rat.R")

# ── Paramètres de base ────────────────────────────────────
pars_fornari_pd <- init_pars
pars_fornari_pd[["k_dam"]] <- NULL
pars_fornari_pd[["k_rep"]] <- NULL

pars_base <- c(pars_fornari_pd, tdxd_pars)

state_pd    <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_full <- c(tdxd_state0, state_pd)

MEP0_ref <- pars_base$MEP0
Ret0_ref <- pars_base$Ret0

# ── Fonction : simuler et extraire nadirs MEP et Ret ─────
sim_nadirs <- function(dose_mgkg, Slope_MEP_val) {
  pars_i <- pars_base
  pars_i$Slope_MEP <- Slope_MEP_val
  pars_i$rate_fun  <- make_tdxd_infusion(
    dose_mgkg  = dose_mgkg, BW_kg = 0.25,
    Tinfu_h    = 0.5,
    interval_h = 21 * 24,
    n_cycles   = 3
  )
  times <- seq(0, 63 * 24, by = 2)

  # Silencieux
  out <- suppressMessages(suppressWarnings(
    as.data.frame(lsoda(
      y     = state0_full,
      times = times,
      func  = pkpd_tdxd_fornari,
      parms = pars_i,
      rtol  = 1e-6, atol = 1e-8,
      maxsteps = 300000
    ))
  ))
  list(
    MEP_nadir_pct = 100 * (min(out$MEP, na.rm=TRUE) - MEP0_ref) / MEP0_ref,
    Ret_nadir_pct = 100 * (min(out$Ret, na.rm=TRUE) - Ret0_ref) / Ret0_ref
  )
}

# ── Scan de Slope_MEP ────────────────────────────────────
cat("═══════════════════════════════════════════════════════════\n")
cat("  SCAN Slope_MEP — Calibration T-DXd Rat\n")
cat("  Cible : MEP@20 < -10%  |  MEP@60 > -20%\n")
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("  %-10s  %-14s  %-14s  %-14s  %-14s  %s\n",
            "Slope_MEP", "MEP@20", "MEP@60", "Ret@20", "Ret@60", "Status"))
cat(sprintf("  %-10s  %-14s  %-14s  %-14s  %-14s  %s\n",
            "─────────", "──────────────","──────────────",
            "──────────────","──────────────","───────"))

slope_vals <- c(2.19, 1.8, 1.56, 1.4, 1.2, 1.0, 0.85)

results <- data.frame(
  Slope_MEP    = slope_vals,
  MEP20_pct    = NA_real_,
  MEP60_pct    = NA_real_,
  Ret20_pct    = NA_real_,
  Ret60_pct    = NA_real_
)

for (i in seq_along(slope_vals)) {
  s <- slope_vals[i]
  r20 <- sim_nadirs(20,  s)
  r60 <- sim_nadirs(60,  s)

  results$MEP20_pct[i] <- r20$MEP_nadir_pct
  results$MEP60_pct[i] <- r60$MEP_nadir_pct
  results$Ret20_pct[i] <- r20$Ret_nadir_pct
  results$Ret60_pct[i] <- r60$Ret_nadir_pct

  ok20 <- r20$MEP_nadir_pct > -10
  ok60 <- r60$MEP_nadir_pct < -20
  status <- if (ok20 && ok60) "OK  ✓" else if (ok20) "MEP60 trop faible" else if (ok60) "MEP20 trop fort" else "les deux hors cible"

  cat(sprintf("  %-10.2f  %-14s  %-14s  %-14s  %-14s  %s\n",
              s,
              sprintf("%+.1f%%", r20$MEP_nadir_pct),
              sprintf("%+.1f%%", r60$MEP_nadir_pct),
              sprintf("%+.1f%%", r20$Ret_nadir_pct),
              sprintf("%+.1f%%", r60$Ret_nadir_pct),
              status))
}

cat("═══════════════════════════════════════════════════════════\n")

# ── Sélection de la valeur calibrée ──────────────────────
valid <- results[results$MEP20_pct > -10 & results$MEP60_pct < -20, ]

if (nrow(valid) > 0) {
  # Prendre la valeur la plus élevée (la plus conservatrice)
  best <- valid[which.max(valid$Slope_MEP), ]
  cat(sprintf("\n  Valeur retenue : Slope_MEP_tdxd_rat = %.2f\n", best$Slope_MEP))
  cat(sprintf("    MEP@20 = %+.1f%%  (cible < -10%%, FDA : pas d'effet)\n", best$MEP20_pct))
  cat(sprintf("    MEP@60 = %+.1f%%  (cible < -20%%, FDA : érythroblastes↓)\n", best$MEP60_pct))
} else {
  cat("\n  ATTENTION : aucune valeur ne satisfait les deux critères\n")
  cat("  → interpolation linéaire nécessaire\n")
  # Interpolation
  df_ok20 <- results[results$MEP20_pct > -10, ]
  if (nrow(df_ok20) > 0) {
    cat(sprintf("  → Valeurs satisfaisant MEP@20 < -10%% : Slope_MEP ≤ %.2f\n",
                max(df_ok20$Slope_MEP)))
  }
}
