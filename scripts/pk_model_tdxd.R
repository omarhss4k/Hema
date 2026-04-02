############################################################
# pk_model_tdxd.R
# Modèle PK T-DXd — Chaîne : ADC sérum → Libération DXd
#
# États :
#   C_ADC1  : ADC compartiment central       [mg/L]
#   C_ADC2  : ADC compartiment périphérique  [mg/L]
#   C_DXd   : DXd plasma (payload libre)     [mg/L]
#
# Équations :
#   dC_ADC1 = rate_in/V1 - (CL_ADC/V1 + k_int)*C_ADC1
#             - Q/V1*C_ADC1 + Q/V2*C_ADC2
#   dC_ADC2 = Q/V1*C_ADC1 - Q/V2*C_ADC2
#   dC_DXd  = mass_frac_DXd * k_int * C_ADC1 * V1/V_DXd
#             - CL_DXd/V_DXd * C_DXd
#
# Note sur les unités :
#   mass_frac_DXd × k_int × C_ADC1 × V1  →  mg DXd/h
#   divisé par V_DXd                      →  mg/L/h = (mg/L)/h
#   C_DXd_uM = C_DXd [mg/L] × 1000 / MW_DXd [g/mol]
############################################################
library(deSolve)

# ── ODE ──────────────────────────────────────────────────
pk_tdxd_ode <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    # Taux de perfusion externe
    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0

    # ADC — 2 compartiments avec fuite par internalisation
    dC_ADC1 <- rate_in / V1_ADC +
               (Q_ADC / V2_ADC) * C_ADC2 -
               (CL_ADC / V1_ADC + Q_ADC / V1_ADC + k_int) * C_ADC1

    dC_ADC2 <- (Q_ADC / V1_ADC) * C_ADC1 -
               (Q_ADC / V2_ADC) * C_ADC2

    # DXd — libéré par l'internalisation de l'ADC
    # Source : mass_frac_DXd × k_int × (C_ADC1 × V1_ADC) / V_DXd
    release_DXd <- mass_frac_DXd * k_int * C_ADC1 * V1_ADC / V_DXd

    dC_DXd <- release_DXd - (CL_DXd / V_DXd) * C_DXd

    list(c(dC_ADC1, dC_ADC2, dC_DXd))
  })
}

# ── Wrapper simulation ───────────────────────────────────
simulate_pk_tdxd <- function(times, pars, state0,
                             rtol = 1e-8, atol = 1e-10) {
  out <- as.data.frame(lsoda(
    y        = state0,
    times    = sort(unique(times)),
    func     = pk_tdxd_ode,
    parms    = pars,
    rtol     = rtol,
    atol     = atol,
    maxsteps = 500000
  ))

  # Colonnes dérivées
  out$time_h  <- out$time
  out$time_d  <- out$time / 24

  # Conversion DXd → µM
  out$C_DXd_uM <- out$C_DXd * pars$mgL_to_uM_DXd

  # Rapport DXd_uM / IC50 (indice d'exposition pharmacologique)
  out$DXd_over_IC50 <- out$C_DXd_uM / pars$IC50_DXd_uM

  # Diagnostics
  Cmax_ADC  <- max(out$C_ADC1,     na.rm = TRUE)
  Cmax_DXd  <- max(out$C_DXd_uM,  na.rm = TRUE)
  Tmax_DXd  <- out$time_h[which.max(out$C_DXd_uM)]
  AUC_ADC   <- tryCatch(
    sum(diff(out$time_h) * (head(out$C_ADC1, -1) + tail(out$C_ADC1, -1)) / 2),
    error = function(e) NA_real_
  )
  AUC_DXd   <- tryCatch(
    sum(diff(out$time_h) * (head(out$C_DXd_uM, -1) + tail(out$C_DXd_uM, -1)) / 2),
    error = function(e) NA_real_
  )

  cat(sprintf("  ✓ ADC  : Cmax=%.2f mg/L  AUC=%.1f mg·h/L\n",
              Cmax_ADC, AUC_ADC))
  cat(sprintf("  ✓ DXd  : Cmax=%.4f µM  Tmax=%.1f h  AUC=%.3f µM·h\n",
              Cmax_DXd, Tmax_DXd, AUC_DXd))
  cat(sprintf("  ✓ DXd/IC50 max = %.3f  (%s)\n",
              Cmax_DXd / pars$IC50_DXd_uM,
              ifelse(Cmax_DXd > pars$IC50_DXd_uM, "exposition > IC50", "exposition < IC50")))

  out
}

# ── Graphiques PK ────────────────────────────────────────
plot_pk_tdxd <- function(sim, pars, titre = "T-DXd PK — Rat",
                         dose_times_h = 0, file = NULL) {
  if (!is.null(file)) pdf(file, width = 10, height = 8)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # ADC central (linéaire)
  plot(sim$time_h, sim$C_ADC1, type = "l", lwd = 2, col = "#2166ac",
       xlab = "Temps (h)", ylab = "ADC [mg/L]",
       main = "ADC — Compartiment central")
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # ADC central (log)
  idx <- sim$C_ADC1 > 0
  plot(sim$time_h[idx], log10(sim$C_ADC1[idx]), type = "l", lwd = 2, col = "#2166ac",
       xlab = "Temps (h)", ylab = "log10 ADC [mg/L]",
       main = "ADC — Échelle log")
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # DXd (µM)
  plot(sim$time_h, sim$C_DXd_uM, type = "l", lwd = 2, col = "#d6604d",
       xlab = "Temps (h)", ylab = "DXd [µM]",
       main = "DXd payload — Plasma")
  abline(h = pars$IC50_DXd_uM, lty = 3, col = "#d6604d", lwd = 1.5)
  text(max(sim$time_h) * 0.6, pars$IC50_DXd_uM * 1.15,
       sprintf("IC50 = %.2f µM", pars$IC50_DXd_uM), col = "#d6604d", cex = 0.85)
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # DXd / IC50
  plot(sim$time_h, sim$DXd_over_IC50, type = "l", lwd = 2, col = "#4dac26",
       xlab = "Temps (h)", ylab = "DXd / IC50",
       main = "Indice d'exposition DXd/IC50")
  abline(h = 1, lty = 3, col = "grey40")
  abline(v = dose_times_h, lty = 2, col = "grey60")

  mtext(titre, outer = TRUE, line = -1.5, font = 2, cex = 1.1)

  if (!is.null(file)) {
    dev.off()
    cat(sprintf("  -> Graphiques sauvegardés : %s\n", file))
  }
}
