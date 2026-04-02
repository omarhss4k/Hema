############################################################
# pk_model_tdxd.R
# Modèle PK/Damage T-DXd — Chaîne complète :
#   ADC sérum → Libération DXd → DXd intracell. → Dommage ADN
#
# États (5) :
#   C_ADC1   : ADC compartiment central          [mg/L]
#   C_ADC2   : ADC compartiment périphérique     [mg/L]
#   C_DXd    : DXd plasma                        [mg/L]
#   C_DXd_ic : DXd intracellulaire (moelle)      [mg/L]
#   Damage   : Dommages ADN normalisés (γH2AX)   [sans unité]
#
# Équations :
#   dC_ADC1   = rate_in/V1 - (CL_ADC/V1 + Q/V1 + k_int)*C_ADC1
#               + Q/V2*C_ADC2
#   dC_ADC2   = Q/V1*C_ADC1 - Q/V2*C_ADC2
#   dC_DXd    = mass_frac_DXd*k_int*C_ADC1*V1/V_DXd
#               - CL_DXd/V_DXd*C_DXd
#               + k_effD*C_DXd_ic*(V_ic/V_DXd)
#               - k_inD*C_DXd
#   dC_DXd_ic = k_inD*C_DXd*(V_DXd/V_ic) - k_effD*C_DXd_ic
#   dDamage   = k_dam * E_drug - k_rep * Damage
#               où E_drug = Emax (C_DXd_ic_uM / (IC50 + C_DXd_ic_uM))
#
# Connexion avec Fornari :
#   Ce "Damage" remplace directement le "Damage" carboplatin dans
#   pkpd_model_FORNARI.R (termes Slope_MPP/CMP/MEP × Damage).
############################################################
library(deSolve)

# ── ODE ──────────────────────────────────────────────────
pk_tdxd_ode <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    # Taux de perfusion externe [mg/h]
    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0

    # ── ADC — 2 compartiments avec fuite par internalisation ──
    dC_ADC1 <- rate_in / V1_ADC +
               (Q_ADC / V2_ADC) * C_ADC2 -
               (CL_ADC / V1_ADC + Q_ADC / V1_ADC + k_int) * C_ADC1

    dC_ADC2 <- (Q_ADC / V1_ADC) * C_ADC1 -
               (Q_ADC / V2_ADC) * C_ADC2

    # ── DXd plasma — alimenté par internalisation + efflux cellulaire ──
    # Source ADC : mass_frac × k_int × (C_ADC1 × V1) / V_DXd
    release_ADC <- mass_frac_DXd * k_int * C_ADC1 * V1_ADC / V_DXd

    # Échanges membranaires (Vasalou 2024)
    flux_in_cell  <- k_inD  * C_DXd                    # plasma → cellule
    flux_out_cell <- k_effD * C_DXd_ic * (V_ic / V_DXd) # cellule → plasma

    dC_DXd <- release_ADC - (CL_DXd / V_DXd) * C_DXd -
              flux_in_cell + flux_out_cell

    # ── DXd intracellulaire (moelle osseuse) ──
    # Conservation de masse : flux entrant rapporté au volume V_ic
    dC_DXd_ic <- k_inD * C_DXd * (V_DXd / V_ic) - k_effD * C_DXd_ic

    # ── Dommages ADN (γH2AX, modèle Emax) ──
    C_DXd_ic_uM <- C_DXd_ic * mgL_to_uM_DXd             # [µM]
    E_drug       <- C_DXd_ic_uM / (IC50_DXd_uM + C_DXd_ic_uM)  # Emax [0,1]

    dDamage <- k_dam * E_drug - k_rep * Damage

    list(c(dC_ADC1, dC_ADC2, dC_DXd, dC_DXd_ic, dDamage))
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

  out$time_h       <- out$time
  out$time_d       <- out$time / 24
  out$C_DXd_uM     <- out$C_DXd    * pars$mgL_to_uM_DXd
  out$C_DXd_ic_uM  <- out$C_DXd_ic * pars$mgL_to_uM_DXd
  out$accum_ratio  <- ifelse(out$C_DXd > 0,
                             out$C_DXd_ic / out$C_DXd, 0)

  # Diagnostics
  Cmax_ADC     <- max(out$C_ADC1,       na.rm = TRUE)
  Cmax_DXd_pl  <- max(out$C_DXd_uM,    na.rm = TRUE)
  Cmax_DXd_ic  <- max(out$C_DXd_ic_uM, na.rm = TRUE)
  Tmax_DXd_ic  <- out$time_h[which.max(out$C_DXd_ic_uM)]
  Damage_max   <- max(out$Damage,       na.rm = TRUE)

  cat(sprintf("  ✓ ADC        : Cmax = %.2f mg/L\n", Cmax_ADC))
  cat(sprintf("  ✓ DXd plasma : Cmax = %.5f µM  (%.2f%% IC50)\n",
              Cmax_DXd_pl, 100 * Cmax_DXd_pl / pars$IC50_DXd_uM))
  cat(sprintf("  ✓ DXd intra  : Cmax = %.4f µM  Tmax = %.1f h  (%.1f%% IC50)\n",
              Cmax_DXd_ic, Tmax_DXd_ic,
              100 * Cmax_DXd_ic / pars$IC50_DXd_uM))
  cat(sprintf("  ✓ Ratio ic/plasma (éq) = %.1fx\n",
              pars$V_DXd / pars$V_ic))
  cat(sprintf("  ✓ Damage max = %.4f  %s\n",
              Damage_max,
              ifelse(Damage_max * 2.05 > 1,
                     "<<< Slope_MPP × Damage > 1 : nadir potentiellement profond",
                     "OK (< 1/Slope)")))
  out
}

# ── Graphiques PK/Damage ─────────────────────────────────
plot_pk_tdxd <- function(sim, pars, titre = "T-DXd PK — Rat",
                         dose_times_h = 0, file = NULL) {
  if (!is.null(file)) pdf(file, width = 12, height = 9)

  par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1))

  # 1 — ADC central (log)
  idx <- sim$C_ADC1 > 1e-8
  plot(sim$time_h[idx], log10(sim$C_ADC1[idx]),
       type = "l", lwd = 2, col = "#2166ac",
       xlab = "Temps (h)", ylab = "log10 ADC [mg/L]",
       main = "ADC sérum — Échelle log")
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # 2 — DXd plasma (µM)
  plot(sim$time_h, sim$C_DXd_uM,
       type = "l", lwd = 2, col = "#4393c3",
       xlab = "Temps (h)", ylab = "DXd plasma [µM]",
       main = "DXd plasma")
  abline(h = pars$IC50_DXd_uM, lty = 3, col = "grey40")
  text(max(sim$time_h) * 0.55, pars$IC50_DXd_uM * 1.2,
       sprintf("IC50 = %.2f µM", pars$IC50_DXd_uM), col = "grey40", cex = 0.8)
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # 3 — DXd intracellulaire (µM)
  plot(sim$time_h, sim$C_DXd_ic_uM,
       type = "l", lwd = 2, col = "#d6604d",
       xlab = "Temps (h)", ylab = "DXd intracell. [µM]",
       main = "DXd intracellulaire (moelle)")
  abline(h = pars$IC50_DXd_uM, lty = 3, col = "#d6604d", lwd = 1.5)
  text(max(sim$time_h) * 0.55, pars$IC50_DXd_uM * 1.15,
       sprintf("IC50 = %.2f µM", pars$IC50_DXd_uM), col = "#d6604d", cex = 0.8)
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # 4 — Ratio ic/plasma
  plot(sim$time_h, sim$accum_ratio,
       type = "l", lwd = 2, col = "#762a83",
       xlab = "Temps (h)", ylab = "C_ic / C_plasma",
       main = "Ratio accumulation intra/plasma")
  abline(h = pars$V_DXd / pars$V_ic, lty = 3, col = "grey50")
  text(max(sim$time_h) * 0.55, pars$V_DXd / pars$V_ic * 1.05,
       sprintf("Éq. = %.0fx", pars$V_DXd / pars$V_ic), col = "grey50", cex = 0.8)
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # 5 — Damage
  plot(sim$time_h, sim$Damage,
       type = "l", lwd = 2, col = "#1a9641",
       xlab = "Temps (h)", ylab = "Damage (γH2AX norm.)",
       main = "Dommages ADN")
  abline(v = dose_times_h, lty = 2, col = "grey60")

  # 6 — E_drug (Emax)
  E_drug <- sim$C_DXd_ic_uM / (pars$IC50_DXd_uM + sim$C_DXd_ic_uM)
  plot(sim$time_h, E_drug,
       type = "l", lwd = 2, col = "#f4a582",
       xlab = "Temps (h)", ylab = "E_drug = Cic / (IC50 + Cic)",
       main = "Effet Emax (driver Damage)", ylim = c(0, max(E_drug) * 1.2))
  abline(h = 0.5, lty = 3, col = "grey50")
  text(max(sim$time_h) * 0.55, 0.52, "E = 0.5", col = "grey50", cex = 0.8)
  abline(v = dose_times_h, lty = 2, col = "grey60")

  mtext(titre, outer = TRUE, line = -1.2, font = 2, cex = 1.1)

  if (!is.null(file)) {
    dev.off()
    cat(sprintf("  -> Graphiques sauvegardés : %s\n", file))
  }
}
