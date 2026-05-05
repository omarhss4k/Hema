############################################################
# plot_rat_validation.R
# Figure de présentation — Validation rat vs FDA BLA 761139
# 3 doses clés : 20, 60, 197 mg/kg Q3W × 3 cycles
############################################################
library(deSolve)

source("../etape1_fornari_carboplatin_rat/parameters_rat.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("parameters_tdxd_rat.R")
source("pkpd_tdxd_rat.R")

# ── Paramètres ────────────────────────────────────────────
pars_fornari_pd <- init_pars
pars_fornari_pd[["k_dam"]] <- NULL
pars_fornari_pd[["k_rep"]] <- NULL
pars_full <- c(pars_fornari_pd, tdxd_pars)
pars_full$Slope_MEP <- Slope_MEP_tdxd_rat   # 1.00

state_pd    <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_full <- c(tdxd_state0, state_pd)
times       <- seq(0, 63 * 24, by = 2)   # 63 jours = 3 cycles Q3W

# ── Simulations 3 doses FDA ───────────────────────────────
doses <- c(20, 60, 197)
cols  <- c("#4dac26", "#f4a582", "#b2182b")   # vert → orange → rouge

cat("Simulation des 3 doses FDA (20 / 60 / 197 mg/kg)...\n")
sims <- lapply(doses, function(d) {
  cat(sprintf("  %d mg/kg...\n", d))
  p <- pars_full
  p$rate_fun <- make_tdxd_infusion(d, 0.25, 0.5, 21*24, 3)
  suppressMessages(suppressWarnings(simulate_pkpd_tdxd(times, p, state0_full)))
})

Ret0  <- pars_full$Ret0;  MEP0 <- pars_full$MEP0;  Neut0 <- pars_full$Neut0
pct <- function(s, v, r) round(100*(min(s[[v]],na.rm=T)-r)/r, 1)

ret_pct  <- sapply(sims, pct, "Ret",  Ret0)
mep_pct  <- sapply(sims, pct, "MEP",  MEP0)
neut_pct <- sapply(sims, pct, "Neut", Neut0)

cat("\nRésultats :\n")
for (i in 1:3) cat(sprintf("  %3d mg/kg — Ret=%+.1f%%  MEP=%+.1f%%  Neut=%+.1f%%\n",
                            doses[i], ret_pct[i], mep_pct[i], neut_pct[i]))

# ══════════════════════════════════════════════════════════
# FIGURE 1 — Résumé validation : nadirs vs seuils FDA
# ══════════════════════════════════════════════════════════
pdf("results_PKPD/rat_validation_FDA_summary.pdf", width = 10, height = 6)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 1.5), bg = "white", oma = c(0, 0, 3, 0))

plot_bar <- function(pcts, title, fda_seuil_txt, fda_dose_txt,
                     seuil1 = -20, seuil2 = NULL) {
  y_min <- min(pcts) * 1.35
  bp <- barplot(pcts,
                names.arg = paste0(doses, "\nmg/kg"),
                col    = cols,
                ylim   = c(y_min, 8),
                ylab   = "%Δ nadir vs baseline",
                main   = title,
                border = "white",
                cex.names = 1.05, cex.axis = 0.95, cex.lab = 1.0,
                cex.main  = 1.1)
  abline(h = 0, col = "grey50", lwd = 1)

  # Seuil de détection histologique
  abline(h = seuil1, lty = 2, col = "#b2182b", lwd = 2)
  text(0.2, seuil1 - 1.5,
       sprintf("Seuil détection\n(%+d%%)", seuil1),
       col = "#b2182b", cex = 0.78, adj = 0)

  if (!is.null(seuil2)) {
    abline(h = seuil2, lty = 3, col = "#f4a582", lwd = 1.5)
    text(0.2, seuil2 + 1, sprintf("%+d%%", seuil2),
         col = "#f4a582", cex = 0.75, adj = 0)
  }

  # Valeurs au-dessus/dessous des barres
  for (i in seq_along(pcts)) {
    v <- pcts[i]
    vjust <- if (v < 0) v - 2.5 else v + 1.5
    text(bp[i], vjust, sprintf("%.1f%%", v),
         cex = 0.88, col = cols[i], font = 2)
  }

  # Annotation FDA
  mtext(fda_seuil_txt, side = 3, line = 0.1, cex = 0.78, col = "#555555")
  mtext(fda_dose_txt,  side = 1, line = 3.8, cex = 0.82, col = "#b2182b", font = 2)
}

plot_bar(ret_pct,
         "Réticulocytes (sang)",
         "FDA : chute détectable à ≥ 20 mg/kg",
         "→ modèle reproduit le seuil ✓")

plot_bar(mep_pct,
         "Érythroblastes / MEP (moelle)",
         "FDA : chute détectable à ≥ 60 mg/kg",
         "→ modèle reproduit le seuil ✓",
         seuil2 = -10)

plot_bar(neut_pct,
         "Neutrophiles (sang)",
         "FDA : chute détectable à ≥ 197 mg/kg",
         "→ modèle reproduit le seuil ✓")

mtext("Validation rat — T-DXd Q3W × 3 cycles  |  FDA BLA 761139",
      outer = TRUE, cex = 1.15, font = 2, line = 1.2)

dev.off()
cat("  -> results_PKPD/rat_validation_FDA_summary.pdf\n")

# ══════════════════════════════════════════════════════════
# FIGURE 2 — Profils temporels (Ret + MEP + Neut)
# ══════════════════════════════════════════════════════════
pdf("results_PKPD/rat_profils_temporels.pdf", width = 12, height = 5)
par(mfrow = c(1, 3), mar = c(4, 4.5, 3.5, 1), bg = "white",
    oma = c(0, 0, 2.5, 0))

dose_days <- c(0, 21, 42)

plot_profile <- function(var, ref, ylabel, title) {
  y_min <- min(sapply(sims, function(s) min(s[[var]], na.rm=TRUE))) * 0.88
  times_d <- sims[[1]]$time_d
  plot(times_d, sims[[1]][[var]], type = "n",
       xlim = c(0, 63), ylim = c(y_min, ref * 1.12),
       xlab = "Temps (jours)", ylab = ylabel, main = title,
       cex.lab = 1.0, cex.main = 1.05)
  abline(v = dose_days, lty = 3, col = "grey75", lwd = 1.2)
  abline(h = ref, lty = 3, col = "grey40", lwd = 1.2)
  for (i in 1:3) lines(sims[[i]]$time_d, sims[[i]][[var]],
                       col = cols[i], lwd = 2.2)
  legend("bottomright", paste0(doses, " mg/kg"),
         col = cols, lwd = 2.2, bty = "n", cex = 0.88)
  text(64, ref * 1.04, "Baseline", col = "grey40", cex = 0.78, adj = 1)
  for (d in dose_days) text(d + 1, y_min * 1.03,
                             paste0("C", which(dose_days == d) + 0),
                             col = "grey60", cex = 0.75)
}

plot_profile("Ret",  Ret0,  "Ret [×10⁶/kg]",  "Réticulocytes")
plot_profile("MEP",  MEP0,  "MEP [×10⁶/kg]",  "Érythroblastes (MEP)")
plot_profile("Neut", Neut0, "Neut [×10⁶/kg]", "Neutrophiles")

mtext("Profils temporels — T-DXd rat (20 / 60 / 197 mg/kg Q3W × 3)",
      outer = TRUE, cex = 1.1, font = 2, line = 0.8)

dev.off()
cat("  -> results_PKPD/rat_profils_temporels.pdf\n")
