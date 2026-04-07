############################################################
# run_pkpd_tdxd_rat.R
# Simulation complète PK/PD T-DXd — RAT
#
# Sources :
#   ../scripts/parameters_rat.R             (baselines PD, rat)
#   ../scripts/parameters_FORNARI_CORRECT.R (constantes PD dérivées)
#   parameters_tdxd_rat.R                   (PK T-DXd)
#   pkpd_tdxd_rat.R                         (ODE fusionné)
#   ../scripts/plots.R                      (graphiques Fornari)
#
# Scénarios :
#   1. Dose unique 5 mg/kg IV Q3W × 3 cycles
#   2. Dose unique 10 mg/kg IV Q3W × 3 cycles
############################################################
library(deSolve)

# ── Chargement des paramètres ────────────────────────────
source("../scripts/parameters_rat.R")             # init_pars, init_state
source("../scripts/parameters_FORNARI_CORRECT.R") # init_pars enrichi
source("parameters_tdxd_rat.R")                   # tdxd_pars, tdxd_state0
source("pkpd_tdxd_rat.R")                         # pkpd_tdxd_fornari()
source("../scripts/plots.R")                      # save_all_cells()

if (!dir.exists("results_PKPD")) dir.create("results_PKPD")

# ── Fusion des paramètres ────────────────────────────────
# Base : paramètres PD Fornari (init_pars, après FORNARI_CORRECT)
# Ajout : paramètres PK T-DXd (tdxd_pars)
# Règle : tdxd_pars a la priorité sur les conflits (k_dam, k_rep)
pars_fornari_pd <- init_pars
pars_fornari_pd[["k_dam"]] <- NULL
pars_fornari_pd[["k_rep"]] <- NULL

pars_full <- c(pars_fornari_pd, tdxd_pars)

# ── Fusion des états initiaux ────────────────────────────
# PK T-DXd : C_ADC1, C_ADC2, C_DXd, C_DXd_ic, Damage (= 0 à t=0)
# PD Fornari : tous sauf C1, C2, Damage (carboplatin → retirés)
state_pd <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_full <- c(tdxd_state0, state_pd)

cat("═══════════════════════════════════════════════════════════\n")
cat("  MODÈLE FUSIONNÉ : T-DXd PK + Fornari PD\n")
cat(sprintf("  %d états  |  %d paramètres\n",
            length(state0_full), length(pars_full)))
cat(sprintf("  Krel : %.4f (C1) → %.4f (C2) → %.4f h⁻¹ (C3)  [Yin 2020]\n",
            pars_full$k_rel_c1,
            pars_full$k_rel_c1 * 2^pars_full$krel_power * pars_full$krel_factor,
            pars_full$k_rel_c1 * 3^pars_full$krel_power * pars_full$krel_factor))
cat("═══════════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════
# Scénario 1 — 5 mg/kg Q3W × 3 cycles
# ══════════════════════════════════════════════════════════
cat("=== Scénario 1 : T-DXd 5 mg/kg Q3W × 3 cycles ===\n")

pars_s1          <- pars_full
pars_s1$rate_fun <- make_tdxd_infusion(
  dose_mgkg  = 5, BW_kg = 0.25,
  Tinfu_h    = 0.5,
  interval_h = 21 * 24,
  n_cycles   = 3
)

times_s1     <- seq(0, 63 * 24, by = 1)   # 63 jours, pas = 1h
dose_days_s1 <- c(0, 21, 42)

sim_s1 <- simulate_pkpd_tdxd(times_s1, pars_s1, state0_full)

save_all_cells(
  sim        = sim_s1,
  pars       = pars_s1,
  file       = "results_PKPD/PKPD_5mgkg_Q3Wx3_cells.pdf",
  titre      = "T-DXd 5 mg/kg Q3W × 3 — PD hématologique (Rat)",
  dose_days  = dose_days_s1
)

# Graphique PK associé
pdf("results_PKPD/PKPD_5mgkg_Q3Wx3_pk.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4.2, 3, 1))

# ADC
plot(sim_s1$time_d, sim_s1$C_ADC1,
     type = "l", lwd = 2, col = "#2166ac",
     xlab = "Temps (jours)", ylab = "ADC [mg/L]",
     main = "ADC sérum")
abline(v = dose_days_s1, lty = 2, col = "grey60")

# DXd intracell.
plot(sim_s1$time_d, sim_s1$C_DXd_ic_uM,
     type = "l", lwd = 2, col = "#d6604d",
     xlab = "Temps (jours)", ylab = "DXd intracell. [µM]",
     main = "DXd intracellulaire")
abline(h = pars_full$IC50_DXd_uM, lty = 3, col = "#d6604d")
abline(v = dose_days_s1, lty = 2, col = "grey60")

# Damage
plot(sim_s1$time_d, sim_s1$Damage,
     type = "l", lwd = 2, col = "#1a9641",
     xlab = "Temps (jours)", ylab = "Damage (γH2AX norm.)",
     main = "Dommages ADN → driver PD")
abline(v = dose_days_s1, lty = 2, col = "grey60")

dev.off()
cat("  -> results_PKPD/PKPD_5mgkg_Q3Wx3_pk.pdf\n\n")

# ══════════════════════════════════════════════════════════
# Scénario 2 — 10 mg/kg Q3W × 3 cycles
# ══════════════════════════════════════════════════════════
cat("=== Scénario 2 : T-DXd 10 mg/kg Q3W × 3 cycles ===\n")

pars_s2          <- pars_full
pars_s2$rate_fun <- make_tdxd_infusion(
  dose_mgkg  = 10, BW_kg = 0.25,
  Tinfu_h    = 0.5,
  interval_h = 21 * 24,
  n_cycles   = 3
)

sim_s2 <- simulate_pkpd_tdxd(times_s1, pars_s2, state0_full)

save_all_cells(
  sim        = sim_s2,
  pars       = pars_s2,
  file       = "results_PKPD/PKPD_10mgkg_Q3Wx3_cells.pdf",
  titre      = "T-DXd 10 mg/kg Q3W × 3 — PD hématologique (Rat)",
  dose_days  = dose_days_s1
)

# ══════════════════════════════════════════════════════════
# Résumé
# ══════════════════════════════════════════════════════════
cat("\n═══════════════════════════════════════════════════════════\n")
cat("  COMPARAISON NADIRS\n")
cat("─────────────────────────────────────────────────────────\n")
cat(sprintf("  %-10s  %-12s  %-12s  %-12s  %-12s\n",
            "Dose", "Neut nadir", "Plt nadir", "Ret nadir", "Damage max"))
cat(sprintf("  %-10s  %-12s  %-12s  %-12s  %-12s\n",
            "baseline",
            sprintf("%.2f",  pars_full$Neut0),
            sprintf("%.0f",  pars_full$Plt0),
            sprintf("%.0f",  pars_full$Ret0), "—"))
cat(sprintf("  %-10s  %-12s  %-12s  %-12s  %-12s\n",
            "5 mg/kg",
            sprintf("%.4f", min(sim_s1$Neut)),
            sprintf("%.1f",  min(sim_s1$Plt)),
            sprintf("%.1f",  min(sim_s1$Ret)),
            sprintf("%.4f",  max(sim_s1$Damage))))
cat(sprintf("  %-10s  %-12s  %-12s  %-12s  %-12s\n",
            "10 mg/kg",
            sprintf("%.4f", min(sim_s2$Neut)),
            sprintf("%.1f",  min(sim_s2$Plt)),
            sprintf("%.1f",  min(sim_s2$Ret)),
            sprintf("%.4f",  max(sim_s2$Damage))))
cat("═══════════════════════════════════════════════════════════\n")

# ── Repères FDA BLA 761139 — toxicologie rat Q3W × 3 ─────
# DS-8201a ne lie pas HER2 rat → effets purement DXd-dépendants
cat("\n  Repères toxicologie rat (FDA BLA 761139, Q3W × 3) :\n")
cat("  ├ NOAEL hémato  : ~10 mg/kg     (pas d'effet à 10 mg/kg)\n")
cat("  ├ Réticulocytes ↓ à ≥ 20 mg/kg\n")
cat("  ├ Érythroblastes ↓ à ≥ 60 mg/kg\n")
cat("  └ Myélocytes    ↓ à ≥ 197 mg/kg\n")
cat("  -> Ret nadir 5mg/kg =", sprintf("%.1f", min(sim_s1$Ret)),
    "/ baseline", sprintf("%.0f", pars_full$Ret0), "\n")
cat("  -> Ret nadir 10mg/kg =", sprintf("%.1f", min(sim_s2$Ret)),
    "/ baseline", sprintf("%.0f", pars_full$Ret0), "\n")
cat("  (Attendre ↓Ret à 20 mg/kg → simuler ce scénario si calibration requise)\n")

cat("\n═══════════════════════════════════════════════════════════\n")
cat("  Fichiers générés dans results_PKPD/ :\n")
cat("    -> PKPD_5mgkg_Q3Wx3_cells.pdf\n")
cat("    -> PKPD_5mgkg_Q3Wx3_pk.pdf\n")
cat("    -> PKPD_10mgkg_Q3Wx3_cells.pdf\n")
cat("═══════════════════════════════════════════════════════════\n")
