############################################################
# run_pkpd_tdxd_rat.R
# Simulation complète PK/PD T-DXd -- RAT
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

# -- Chargement des paramètres ----------------------------
source("../etape1_fornari_carboplatin_rat/parameters_rat.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("parameters_tdxd_rat.R")
source("pkpd_tdxd_rat.R")
source("../shared/plots.R")

if (!dir.exists("results")) dir.create("results")

# -- Fusion des paramètres --------------------------------
# Base : paramètres PD Fornari (init_pars, après FORNARI_CORRECT)
# Ajout : paramètres PK T-DXd (tdxd_pars)
# Règle : tdxd_pars a la priorité sur les conflits (k_dam, k_rep)
pars_fornari_pd <- init_pars
pars_fornari_pd[["k_dam"]] <- NULL
pars_fornari_pd[["k_rep"]] <- NULL

pars_full <- c(pars_fornari_pd, tdxd_pars)

# -- Override Slope_MEP calibré T-DXd --------------------
# Slope_MEP carboplatin (2.19) produit MEP@20=-14% → au-dessus du seuil
# histopathologique n=4. Calibré à 1.00 : MEP@20=-9.8%, MEP@60=-23.7%
pars_full$Slope_MEP <- Slope_MEP_tdxd_rat   # 1.00, depuis parameters_tdxd_rat.R

# -- Fusion des états initiaux ----------------------------
# PK T-DXd : C_ADC1, C_ADC2, C_DXd, C_DXd_ic, Damage (= 0 à t=0)
# PD Fornari : tous sauf C1, C2, Damage (carboplatin → retirés)
state_pd <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_full <- c(tdxd_state0, state_pd)

cat("═══════════════════════════════════════════════════════════\n")
cat("  MODÈLE FUSIONNÉ : T-DXd PK + Fornari PD\n")
cat(sprintf("  %d états  |  %d paramètres\n",
            length(state0_full), length(pars_full)))
cat(sprintf("  Krel : %.4f (C1) → %.4f (C2) → %.4f h- (C3)  [Yin 2020]\n",
            pars_full$k_rel_c1,
            pars_full$k_rel_c1 * 2^pars_full$krel_power * pars_full$krel_factor,
            pars_full$k_rel_c1 * 3^pars_full$krel_power * pars_full$krel_factor))
cat("═══════════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════
# Scénario 1 -- 5 mg/kg Q3W × 3 cycles
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
  file       = "results/PKPD_5mgkg_Q3Wx3_cells.pdf",
  titre      = "T-DXd 5 mg/kg Q3W × 3 -- PD hématologique (Rat)",
  dose_days  = dose_days_s1
)

# Graphique PK associé
pdf("results/PKPD_5mgkg_Q3Wx3_pk.pdf", width = 12, height = 4)
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
cat("  -> results/PKPD_5mgkg_Q3Wx3_pk.pdf\n\n")

# ══════════════════════════════════════════════════════════
# Scénario 2 -- 10 mg/kg Q3W × 3 cycles
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
  file       = "results/PKPD_10mgkg_Q3Wx3_cells.pdf",
  titre      = "T-DXd 10 mg/kg Q3W × 3 -- PD hématologique (Rat)",
  dose_days  = dose_days_s1
)

# ══════════════════════════════════════════════════════════
# Scénario 3 -- 20 mg/kg Q3W × 3 cycles (FDA seuil Ret)
# ══════════════════════════════════════════════════════════
cat("=== Scénario 3 : T-DXd 20 mg/kg Q3W × 3 cycles (FDA seuil Ret) ===\n")

pars_s3          <- pars_full
pars_s3$rate_fun <- make_tdxd_infusion(
  dose_mgkg  = 20, BW_kg = 0.25,
  Tinfu_h    = 0.5,
  interval_h = 21 * 24,
  n_cycles   = 3
)

sim_s3 <- simulate_pkpd_tdxd(times_s1, pars_s3, state0_full)

# ══════════════════════════════════════════════════════════
# Scénario 4 -- 60 mg/kg Q3W × 3 cycles (FDA seuil Ery/Plt)
# ══════════════════════════════════════════════════════════
cat("=== Scénario 4 : T-DXd 60 mg/kg Q3W × 3 cycles (FDA seuil Ery) ===\n")

pars_s4          <- pars_full
pars_s4$rate_fun <- make_tdxd_infusion(
  dose_mgkg  = 60, BW_kg = 0.25,
  Tinfu_h    = 0.5,
  interval_h = 21 * 24,
  n_cycles   = 3
)

sim_s4 <- simulate_pkpd_tdxd(times_s1, pars_s4, state0_full)

# ══════════════════════════════════════════════════════════
# Scénario 5 -- 197 mg/kg Q3W × 3 cycles (FDA seuil Neut/Mye)
# ══════════════════════════════════════════════════════════
cat("=== Scénario 5 : T-DXd 197 mg/kg Q3W × 3 cycles (FDA seuil Mye) ===\n")

pars_s5          <- pars_full
pars_s5$rate_fun <- make_tdxd_infusion(
  dose_mgkg  = 197, BW_kg = 0.25,
  Tinfu_h    = 0.5,
  interval_h = 21 * 24,
  n_cycles   = 3
)

sim_s5 <- simulate_pkpd_tdxd(times_s1, pars_s5, state0_full)

# ══════════════════════════════════════════════════════════
# Graphique dose-réponse -- validation FDA BLA 761139
# ══════════════════════════════════════════════════════════
pdf("results/PKPD_FDA_validation_dose_response.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1))

doses_sim  <- c(5, 10, 20, 60, 197)
sims_list  <- list(sim_s1, sim_s2, sim_s3, sim_s4, sim_s5)

# -- Couleurs par dose --
dose_cols  <- c("#2166ac", "#4dac26", "#d6604d", "#f4a582", "#b2182b")

# Fonctions utilitaires locales
pct_nadir <- function(sim, col, ref) 100 * (min(sim[[col]], na.rm=TRUE) - ref) / ref

Ret0   <- pars_full$Ret0
MEP0   <- pars_full$MEP0
Neut0  <- pars_full$Neut0
Plt0   <- pars_full$Plt0

ret_pcts  <- sapply(sims_list,  function(s) pct_nadir(s, "Ret",  Ret0))
mep_pcts  <- sapply(sims_list,  function(s) pct_nadir(s, "MEP",  MEP0))
neut_pcts <- sapply(sims_list,  function(s) pct_nadir(s, "Neut", Neut0))
plt_pcts  <- sapply(sims_list,  function(s) pct_nadir(s, "Plt",  Plt0))

# -- 1. Dose-réponse Réticulocytes --
barplot(ret_pcts, names.arg = paste0(doses_sim, "\nmg/kg"),
        col = dose_cols, ylim = c(min(ret_pcts) * 1.2, 5),
        ylab = "Nadir Ret (%Δ baseline)", main = "Réticulocytes (sang)",
        border = NA)
abline(h = -20, lty = 2, col = "red", lwd = 1.5)
text(0.5, -21, "Seuil détection -20%", col = "red", adj = 0, cex = 0.75)
abline(h = 0, col = "grey40")
# Annotation FDA
text(0.7, -15, "FDA: ↓ à ≥20 mg/kg", col = "#d6604d", cex = 0.8, font = 2)

# -- 2. Dose-réponse Érythroblastes (MEP) --
barplot(mep_pcts, names.arg = paste0(doses_sim, "\nmg/kg"),
        col = dose_cols, ylim = c(min(mep_pcts) * 1.2, 5),
        ylab = "Nadir MEP (%Δ baseline)", main = "Érythroblastes / MEP (moelle)",
        border = NA)
abline(h = -20, lty = 2, col = "red", lwd = 1.5)
abline(h = -10, lty = 3, col = "orange", lwd = 1.5)
text(0.5, -21, "Seuil détection -20%", col = "red",    adj = 0, cex = 0.75)
text(0.5, -11, "Seuil detect. -10%",   col = "orange", adj = 0, cex = 0.75)
abline(h = 0, col = "grey40")
text(0.7, -5, "FDA: ↓ à ≥60 mg/kg", col = "#d6604d", cex = 0.8, font = 2)

# -- 3. Dose-réponse Neutrophiles (Myélocytes proxy) --
barplot(neut_pcts, names.arg = paste0(doses_sim, "\nmg/kg"),
        col = dose_cols, ylim = c(min(neut_pcts) * 1.2, 5),
        ylab = "Nadir Neut (%Δ baseline)", main = "Neutrophiles / Myélocytes (moelle)",
        border = NA)
abline(h = -20, lty = 2, col = "red", lwd = 1.5)
text(0.5, -21, "Seuil détection -20%", col = "red", adj = 0, cex = 0.75)
abline(h = 0, col = "grey40")
text(0.7, -5, "FDA: ↓ à ≥197 mg/kg", col = "#d6604d", cex = 0.8, font = 2)

# -- 4. Réticulocytes -- profil temporel multi-doses --
all_times_d <- sim_s1$time_d
y_min <- min(sapply(sims_list, function(s) min(s$Ret, na.rm=TRUE)))
y_max <- Ret0 * 1.15

plot(all_times_d, sim_s1$Ret, type = "n",
     xlim = range(all_times_d), ylim = c(y_min, y_max),
     xlab = "Temps (jours)", ylab = "Ret [×10⁶/kg]",
     main = "Réticulocytes -- profil temporel")
for (i in seq_along(sims_list)) {
  lines(sims_list[[i]]$time_d, sims_list[[i]]$Ret,
        col = dose_cols[i], lwd = 1.5)
}
abline(v = c(0, 21, 42), lty = 2, col = "grey70")
abline(h = Ret0, lty = 3, col = "grey40")
legend("topright", paste0(doses_sim, " mg/kg"), col = dose_cols,
       lwd = 1.5, bty = "n", cex = 0.8)

# -- 5. MEP -- profil temporel multi-doses --
y_min_mep <- min(sapply(sims_list, function(s) min(s$MEP, na.rm=TRUE)))
plot(all_times_d, sim_s1$MEP, type = "n",
     xlim = range(all_times_d), ylim = c(y_min_mep, MEP0 * 1.15),
     xlab = "Temps (jours)", ylab = "MEP [×10⁶/kg]",
     main = "Érythroblastes / MEP -- profil temporel")
for (i in seq_along(sims_list)) {
  lines(sims_list[[i]]$time_d, sims_list[[i]]$MEP,
        col = dose_cols[i], lwd = 1.5)
}
abline(v = c(0, 21, 42), lty = 2, col = "grey70")
abline(h = MEP0, lty = 3, col = "grey40")
legend("topright", paste0(doses_sim, " mg/kg"), col = dose_cols,
       lwd = 1.5, bty = "n", cex = 0.8)

# -- 6. Neutrophiles -- profil temporel multi-doses --
y_min_neut <- min(sapply(sims_list, function(s) min(s$Neut, na.rm=TRUE)))
plot(all_times_d, sim_s1$Neut, type = "n",
     xlim = range(all_times_d), ylim = c(y_min_neut, Neut0 * 1.15),
     xlab = "Temps (jours)", ylab = "Neut [×10⁶/kg]",
     main = "Neutrophiles -- profil temporel")
for (i in seq_along(sims_list)) {
  lines(sims_list[[i]]$time_d, sims_list[[i]]$Neut,
        col = dose_cols[i], lwd = 1.5)
}
abline(v = c(0, 21, 42), lty = 2, col = "grey70")
abline(h = Neut0, lty = 3, col = "grey40")
legend("topright", paste0(doses_sim, " mg/kg"), col = dose_cols,
       lwd = 1.5, bty = "n", cex = 0.8)

dev.off()
cat("  -> results_PKPD/PKPD_FDA_validation_dose_response.pdf\n")

# ══════════════════════════════════════════════════════════
# Résumé
# ══════════════════════════════════════════════════════════
cat("\n═══════════════════════════════════════════════════════════\n")
cat("  VALIDATION FDA BLA 761139 -- DOSE-RÉPONSE RAT Q3W×3\n")
cat("  (Slope_MEP calibré = 1.00 | Slope_MEP carbo = 2.19)\n")
cat("---------------------------------------------------------\n")
cat(sprintf("  %-8s  %-9s  %-9s  %-9s  %-9s  %-9s  %s\n",
            "Dose", "Ret%", "MEP%", "Neut%", "Plt%", "Dmg max", "FDA observé"))
cat(sprintf("  %-8s  %-9s  %-9s  %-9s  %-9s  %-9s  %s\n",
            "--------", "---------","---------","---------","---------","---------","------------------"))

fda_obs <- c(
  "NOAEL hémato",
  "NOAEL hémato",
  "Ret↓",
  "Ret↓ + Ery↓ + Plt↑",
  "Ret↓ + Ery↓ + Mye↓"
)

for (i in seq_along(doses_sim)) {
  sim_i <- sims_list[[i]]
  d     <- doses_sim[i]

  ret_p  <- 100 * (min(sim_i$Ret,  na.rm=TRUE) - Ret0)  / Ret0
  mep_p  <- 100 * (min(sim_i$MEP,  na.rm=TRUE) - MEP0)  / MEP0
  neut_p <- 100 * (min(sim_i$Neut, na.rm=TRUE) - Neut0) / Neut0
  plt_p  <- 100 * (min(sim_i$Plt,  na.rm=TRUE) - Plt0)  / Plt0
  dmg_mx <- max(sim_i$Damage, na.rm=TRUE)

  cat(sprintf("  %-8s  %-9s  %-9s  %-9s  %-9s  %-9s  %s\n",
              paste0(d, " mg/kg"),
              sprintf("%+.1f%%", ret_p),
              sprintf("%+.1f%%", mep_p),
              sprintf("%+.1f%%", neut_p),
              sprintf("%+.1f%%", plt_p),
              sprintf("%.4f",   dmg_mx),
              fda_obs[i]))
}

cat("---------------------------------------------------------\n")
cat("  Seuils FDA (détection histopathologique, n=4) :\n")
cat("  Ret  : seuil ≥20 mg/kg → modèle prédit ↓ à 20 mg/kg ✓\n")
cat("  MEP  : seuil ≥60 mg/kg → MEP@20<10%, MEP@60>20%    ✓\n")
cat("  Neut : seuil ≥197 mg/kg → Neut@60<20%, Neut@197>20%?\n")
cat(sprintf("           Neut@60=%.1f%%  Neut@197=%.1f%%\n",
            100*(min(sim_s4$Neut,na.rm=TRUE)-Neut0)/Neut0,
            100*(min(sim_s5$Neut,na.rm=TRUE)-Neut0)/Neut0))
cat("═══════════════════════════════════════════════════════════\n")

cat("\n═══════════════════════════════════════════════════════════\n")
cat("  Fichiers générés dans results_PKPD/ :\n")
cat("    -> PKPD_5mgkg_Q3Wx3_cells.pdf\n")
cat("    -> PKPD_5mgkg_Q3Wx3_pk.pdf\n")
cat("    -> PKPD_10mgkg_Q3Wx3_cells.pdf\n")
cat("    -> PKPD_FDA_validation_dose_response.pdf  (NEW)\n")
cat("═══════════════════════════════════════════════════════════\n")
