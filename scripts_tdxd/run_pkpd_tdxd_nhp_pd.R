############################################################
# run_pkpd_tdxd_nhp_pd.R
# Simulation PK/PD — T-DXd — SINGE CYNOMOLGUS (NHP)
#
# Protocole FDA BLA 761139 Table 7 :
#   "3-Month Intermittent IV Dose Toxicity Study"
#   Doses Q3W × 5 cycles : 3, 10, 30 mg/kg
#
# Résultats générés dans results_TDXD/ :
#   NHP_PKPD_PK_profiles.pdf   — ADC, DXd, Damage par dose
#   NHP_PKPD_cells.pdf         — profils hémato 8 lignées
#   NHP_PKPD_dose_response.pdf — nadir relatif vs dose
############################################################
library(deSolve)

# ── Chargement paramètres ─────────────────────────────
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("parameters_tdxd_nhp.R")    # PK NHP (tdxd_nhp, make_nhp_infusion)
source("parameters_tdxd_nhp_pd.R") # PD NHP (nhp_pd, nhp_pkpd_state0)
source("pkpd_tdxd_nhp.R")          # ODE pkpd_nhp_ode

if (!dir.exists("results_TDXD")) dir.create("results_TDXD")

# ════════════════════════════════════════════════════════
# Fusion paramètres PK + PD dans un seul vecteur
# ════════════════════════════════════════════════════════
make_nhp_pkpd_pars <- function(dose_mgkg, BW_kg = 4.0,
                                Tinfu_h = 0.5, n_cycles = 5) {
  pars <- c(as.list(tdxd_nhp), as.list(nhp_pd))

  # k_int = 0 (TMDD via MM, pas d'internalisation linéaire)
  pars$k_int <- 0

  # Fonction de perfusion
  pars$rate_fun <- make_nhp_infusion(
    dose_mgkg  = dose_mgkg,
    BW_kg      = BW_kg,
    Tinfu_h    = Tinfu_h,
    interval_h = tdxd_nhp$interval_h,
    n_cycles   = n_cycles
  )
  pars
}

# ════════════════════════════════════════════════════════
# Simulation helper
# ════════════════════════════════════════════════════════
simulate_nhp_pkpd <- function(dose_mgkg, n_cycles = 5,
                               BW_kg = 4.0, Tinfu_h = 0.5,
                               dt_h = 2) {
  pars  <- make_nhp_pkpd_pars(dose_mgkg, BW_kg, Tinfu_h, n_cycles)

  # Inclure les discontinuités de perfusion pour éviter que lsoda
  # ne saute par-dessus la fenêtre Tinfu_h = 0.5h avec dt = 2h
  dose_starts <- seq(0, by = tdxd_nhp$interval_h, length.out = n_cycles)
  disc_times  <- c(dose_starts, dose_starts + Tinfu_h + 1e-6)
  times <- sort(unique(c(seq(0, n_cycles * 21 * 24, by = dt_h), disc_times)))

  sol <- tryCatch(
    suppressWarnings(as.data.frame(ode(
      y      = nhp_pkpd_state0,
      times  = times,
      func   = pkpd_nhp_ode,
      parms  = pars,
      method = "lsoda",
      hmax   = Tinfu_h / 2   # pas interne max < Tinfu_h : le solveur ne peut pas
                              # sauter par-dessus la fenêtre de perfusion de 0.5 h
    ))),
    error = function(e) { cat("ODE error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(sol)) return(NULL)

  sol$time_d        <- sol$time / 24
  sol$C_DXd_ngmL    <- sol$C_DXd    * 1e3
  sol$C_DXd_ic_uM   <- sol$C_DXd_ic * (1000 / tdxd_nhp$MW_DXd) / nhp_pd$IC50_DXd_uM
  sol$Neut_pct      <- sol$Neut / nhp_pd$Neut0 * 100
  sol$Ret_pct       <- sol$Ret  / nhp_pd$Ret0  * 100
  sol$RBC_pct       <- sol$RBC  / nhp_pd$RBC0  * 100
  sol$Plt_pct       <- sol$Plt  / nhp_pd$Plt0  * 100
  sol$MEP_pct       <- sol$MEP  / nhp_pd$MEP0  * 100
  sol$CMP_pct       <- sol$CMP  / nhp_pd$CMP0  * 100
  sol
}

# ════════════════════════════════════════════════════════
# Simulation 3 doses
# ════════════════════════════════════════════════════════
doses     <- c(3, 10, 30)
dose_cols <- c("#2166ac", "#4dac26", "#d6604d")
dose_days <- (0:4) * 21   # 5 cycles Q3W

cat("\n── Simulation Q3W × 5 cycles ─────────────────────────────\n")
sims <- lapply(doses, function(d) {
  cat(sprintf("  %d mg/kg ... ", d))
  s <- simulate_nhp_pkpd(d, n_cycles = 5)
  if (!is.null(s)) cat("OK\n") else cat("ERREUR\n")
  s
})
names(sims) <- paste0(doses, "_mgkg")

# ════════════════════════════════════════════════════════
# FIGURE 1 — Profils PK : ADC, DXd plasma, Damage
# ════════════════════════════════════════════════════════
pdf("results_TDXD/NHP_PKPD_PK_profiles.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4.5, 3.5, 1.5))

# ── ADC sérum (log) ──
ylim_adc <- range(sapply(sims, function(s) range(s$C_ADC1[s$C_ADC1 > 0])),
                  na.rm = TRUE)
plot(NA, xlim = c(0, 105), ylim = ylim_adc, log = "y",
     xlab = "Temps (jours)", ylab = "ADC sérum [µg/mL]",
     main = "ADC DS-8201a — NHP Q3W × 5")
abline(v = dose_days, lty = 2, col = "grey70", lwd = 0.8)
for (i in seq_along(doses)) {
  s <- sims[[i]]
  lines(s$time_d, s$C_ADC1, col = dose_cols[i], lwd = 2)
}
legend("topright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.85)

# ── DXd plasma ──
ylim_dxd <- range(sapply(sims, function(s) range(s$C_DXd_ngmL)),
                  na.rm = TRUE)
plot(NA, xlim = c(0, 105), ylim = c(0, ylim_dxd[2] * 1.1),
     xlab = "Temps (jours)", ylab = "DXd plasma [ng/mL]",
     main = "DXd plasma — NHP Q3W × 5")
abline(v = dose_days, lty = 2, col = "grey70", lwd = 0.8)
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$C_DXd_ngmL, col = dose_cols[i], lwd = 2)
legend("topright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.85)

# ── DXd intracellulaire (ratio IC50) ──
plot(NA, xlim = c(0, 105),
     ylim = c(0, max(sapply(sims, function(s) max(s$C_DXd_ic_uM))) * 1.1),
     xlab = "Temps (jours)", ylab = "DXd_ic / IC50 [%]",
     main = "DXd intracellulaire (ratio IC50)")
abline(v = dose_days, lty = 2, col = "grey70", lwd = 0.8)
abline(h = 100, lty = 3, col = "red", lwd = 1)
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$C_DXd_ic_uM * 100,
        col = dose_cols[i], lwd = 2)
legend("topright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.85)
text(80, 100, "IC50", col = "red", cex = 0.8, adj = c(0, -0.3))

# ── Damage ADN ──
plot(NA, xlim = c(0, 105),
     ylim = c(0, max(sapply(sims, function(s) max(s$Damage))) * 1.1),
     xlab = "Temps (jours)", ylab = "Damage ADN (normalisé)",
     main = "Dommage ADN γH2AX — NHP")
abline(v = dose_days, lty = 2, col = "grey70", lwd = 0.8)
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$Damage, col = dose_cols[i], lwd = 2)
legend("topright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.85)

# ── MEP (progéniteurs érythro-megakaryocytaires) ──
plot(NA, xlim = c(0, 105), ylim = c(50, 110),
     xlab = "Temps (jours)", ylab = "MEP (% baseline)",
     main = "Progéniteurs MEP — NHP")
abline(v = dose_days, lty = 2, col = "grey70", lwd = 0.8)
abline(h = 100, lty = 1, col = "grey50", lwd = 0.8)
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$MEP_pct, col = dose_cols[i], lwd = 2)
legend("bottomright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.85)

# ── CMP (progéniteurs myéloïdes) ──
plot(NA, xlim = c(0, 105), ylim = c(50, 110),
     xlab = "Temps (jours)", ylab = "CMP (% baseline)",
     main = "Progéniteurs CMP — NHP")
abline(v = dose_days, lty = 2, col = "grey70", lwd = 0.8)
abline(h = 100, lty = 1, col = "grey50", lwd = 0.8)
for (i in seq_along(doses))
  lines(sims[[i]]$time_d, sims[[i]]$CMP_pct, col = dose_cols[i], lwd = 2)
legend("bottomright", paste0(doses, " mg/kg"), col = dose_cols,
       lwd = 2, bty = "n", cex = 0.85)

dev.off()
cat("  -> results_TDXD/NHP_PKPD_PK_profiles.pdf\n")

# ════════════════════════════════════════════════════════
# FIGURE 2 — Profils hématologiques (8 cellules)
# ════════════════════════════════════════════════════════
cell_cfg <- list(
  list(var="Neut_pct", lab="Neutrophiles (% baseline)", base=nhp_pd$Neut0,
       unit="×10⁹/L", seuil=100/nhp_pd$Neut0*1.0),   # LLN NHP ~1.0
  list(var="Ret_pct",  lab="Réticulocytes (% baseline)", base=nhp_pd$Ret0,  unit="×10⁹/L", seuil=NA),
  list(var="RBC_pct",  lab="GR / RBC (% baseline)",      base=nhp_pd$RBC0,  unit="×10⁹/L", seuil=NA),
  list(var="Plt_pct",  lab="Plaquettes (% baseline)",    base=nhp_pd$Plt0,  unit="×10⁹/L", seuil=NA),
  list(var="MEP_pct",  lab="MEP (% baseline)",           base=nhp_pd$MEP0,  unit="×10⁶/kg", seuil=NA),
  list(var="CMP_pct",  lab="CMP (% baseline)",           base=nhp_pd$CMP0,  unit="×10⁶/kg", seuil=NA),
  list(var="Mono",     lab="Monocytes [×10⁹/L]",        base=nhp_pd$Mono0, unit="×10⁹/L", seuil=NA),
  list(var="Neut",     lab="Neutrophiles [×10⁹/L]",     base=nhp_pd$Neut0, unit="×10⁹/L", seuil=1.0)
)

pdf("results_TDXD/NHP_PKPD_cells.pdf", width = 16, height = 12)
par(mfrow = c(2, 4), mar = c(4, 4.5, 3.5, 1.5))

for (cfg in cell_cfg) {
  is_pct <- grepl("%", cfg$lab)
  ylim_vals <- if (is_pct) {
    all_y <- unlist(lapply(sims, function(s) s[[cfg$var]]))
    c(min(60, min(all_y, na.rm = TRUE) * 0.95),
      max(110, max(all_y, na.rm = TRUE) * 1.05))
  } else {
    all_y <- unlist(lapply(sims, function(s) s[[cfg$var]]))
    c(min(all_y, na.rm = TRUE) * 0.90,
      max(all_y, na.rm = TRUE) * 1.10)
  }

  plot(NA, xlim = c(0, 105), ylim = ylim_vals,
       xlab = "Temps (jours)", ylab = cfg$lab,
       main = gsub(" \\(.*", "", cfg$lab))
  abline(v = dose_days, lty = 2, col = "grey70", lwd = 0.8)

  if (is_pct) abline(h = 100, lty = 1, col = "grey50", lwd = 0.8)

  if (!is.na(cfg$seuil)) {
    lln_pct <- if (is_pct) cfg$seuil else cfg$seuil
    abline(h = lln_pct, lty = 3, col = "#d62728", lwd = 1.2)
    text(100, lln_pct, "LLN NHP", col = "#d62728", cex = 0.7, adj = c(1, -0.3))
  }

  for (i in seq_along(doses))
    lines(sims[[i]]$time_d, sims[[i]][[cfg$var]], col = dose_cols[i], lwd = 2)

  legend("bottomright", paste0(doses, " mg/kg"), col = dose_cols,
         lwd = 2, bty = "n", cex = 0.75)
}

dev.off()
cat("  -> results_TDXD/NHP_PKPD_cells.pdf\n")

# ════════════════════════════════════════════════════════
# FIGURE 3 — Dose-réponse : nadir relatif par cellule
# ════════════════════════════════════════════════════════
get_nadir <- function(sol, var, baseline) {
  vals <- sol[[var]]
  min(vals, na.rm = TRUE) / baseline * 100
}

nadir_tbl <- data.frame(
  dose = doses,
  Neut = sapply(sims, get_nadir, "Neut", nhp_pd$Neut0),
  Ret  = sapply(sims, get_nadir, "Ret",  nhp_pd$Ret0),
  RBC  = sapply(sims, get_nadir, "RBC",  nhp_pd$RBC0),
  Plt  = sapply(sims, get_nadir, "Plt",  nhp_pd$Plt0),
  MEP  = sapply(sims, get_nadir, "MEP",  nhp_pd$MEP0),
  CMP  = sapply(sims, get_nadir, "CMP",  nhp_pd$CMP0)
)

cat("\n═══ NADIR RELATIF (% baseline) ═══\n")
print(round(nadir_tbl, 1))

pdf("results_TDXD/NHP_PKPD_dose_response.pdf", width = 10, height = 6)
par(mar = c(5, 5, 3.5, 8), xpd = FALSE)

cell_vars  <- c("MEP", "CMP", "Ret", "Plt", "Neut", "RBC")
cell_cols  <- c("#e41a1c", "#ff7f00", "#4daf4a",
                "#984ea3", "#377eb8", "#a65628")
cell_names <- c("MEP", "CMP", "Réticulocytes", "Plaquettes",
                "Neutrophiles", "GR")

plot(NA, xlim = c(2, 35), ylim = c(50, 105),
     xlab = "Dose (mg/kg)", ylab = "Nadir (% baseline)",
     main = "Dose-réponse hématologique — NHP Q3W × 5 cycles",
     log = "x", xaxt = "n")
axis(1, at = doses, labels = paste0(doses, "\nmg/kg"))
abline(h = 100, lty = 1, col = "grey70")
abline(h = c(80, 90), lty = 3, col = "grey80")

for (j in seq_along(cell_vars)) {
  y <- nadir_tbl[[cell_vars[j]]]
  lines(doses, y, col = cell_cols[j], lwd = 2, type = "b", pch = 16, cex = 1.2)
}

par(xpd = TRUE)
legend(x = 38, y = 105,
       legend = cell_names, col = cell_cols,
       lwd = 2, pch = 16, bty = "n", cex = 0.85,
       title = "Lignée")
par(xpd = FALSE)

dev.off()
cat("  -> results_TDXD/NHP_PKPD_dose_response.pdf\n")

# ════════════════════════════════════════════════════════
# RÉSUMÉ CONSOLE
# ════════════════════════════════════════════════════════
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  BILAN PK/PD NHP — T-DXd Q3W × 5 cycles\n")
cat("───────────────────────────────────────────────────────────────\n")
cat(sprintf("  %-8s │ %-8s %-8s │ %-6s %-6s %-6s %-6s\n",
            "Dose", "Cmax ADC", "Cmax DXd", "MEP", "Ret", "Neut", "Plt"))
cat(sprintf("  %-8s │ %-8s %-8s │ %-6s %-6s %-6s %-6s\n",
            "mg/kg", "µg/mL", "ng/mL", "%nadr", "%nadr", "%nadr", "%nadr"))
cat(sprintf("  %s\n", strrep("─", 65)))
for (i in seq_along(doses)) {
  s    <- sims[[i]]
  cmax_adc <- max(s$C_ADC1[s$time <= 24])
  cmax_dxd <- max(s$C_DXd_ngmL[s$time <= 24])
  cat(sprintf("  %-8s │ %-8.1f %-8.3f │ %-6.1f %-6.1f %-6.1f %-6.1f\n",
              paste0(doses[i], " mg/kg"),
              cmax_adc, cmax_dxd,
              nadir_tbl$MEP[i], nadir_tbl$Ret[i],
              nadir_tbl$Neut[i], nadir_tbl$Plt[i]))
}
cat("═══════════════════════════════════════════════════════════════\n")
cat("  NOTE : Effets PD faibles dus aux faibles concentrations\n")
cat(sprintf("  DXd_ic NHP (Krel=%.2e h⁻¹) vs rat (Krel=1.19e-3 h⁻¹)\n",
            tdxd_nhp$k_rel_c1))
cat("  Cohérent avec utilisation du singe comme espèce de sécurité.\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("  Fichiers générés dans results_TDXD/ :\n")
cat("    -> NHP_PKPD_PK_profiles.pdf\n")
cat("    -> NHP_PKPD_cells.pdf\n")
cat("    -> NHP_PKPD_dose_response.pdf\n")
