############################################################
# pkpd_model_FORNARI.R
# ODE système complet — Équations 1–9 (Fornari 2019)
# Feedbacks Eq. 11, 12, 13
############################################################
library(deSolve)

# ── Wrapper simulation ──────────────────────────────────
simulate_all <- function(times, pars, state0,
                         rtol = 1e-7, atol = 1e-9) {
  out <- as.data.frame(lsoda(
    y        = state0,
    times    = sort(unique(times)),
    func     = pkpd_fornari,
    parms    = pars,
    rtol     = rtol,
    atol     = atol,
    maxsteps = 200000
  ))
  out$days <- out$time / 24
  cat(sprintf("  ✓ MPP∈[%.2f,%.0f]  Neut_nadir=%.4f  Plt_nadir=%.1f  Ret_nadir=%.3f\n",
              min(out$MPP, na.rm=TRUE), max(out$MPP, na.rm=TRUE),
              min(out$Neut, na.rm=TRUE),
              min(out$Plt,  na.rm=TRUE),
              min(out$Ret,  na.rm=TRUE)))

  # ── Diagnostic PK/Damage ─────────────────────────────────
  # Si Slope × Damage_max > 1, le kill domine la proliferation
  # → nadir trop profond (signe que les params PK generent trop d'exposition)
  D_max  <- max(out$Damage, na.rm = TRUE)
  C1_max <- max(out$C1,     na.rm = TRUE)
  cat(sprintf("  ✓ PK/PD diag : C1_max=%.2f mg/L (%.1f uM)  Damage_max=%.4f\n",
              C1_max, C1_max * 1000/371.25, D_max))
  for (sl_name in c("Slope_MPP", "Slope_CMP", "Slope_MEP")) {
    sl_val <- pars[[sl_name]]
    if (!is.null(sl_val)) {
      prod <- sl_val * D_max
      cat(sprintf("    %s=%.3f  x  Damage_max=%.4f  =  %.3f  %s\n",
                  sl_name, sl_val, D_max, prod,
                  ifelse(prod > 1, "<<< KILL DOMINANT : nadir trop profond", "OK")))
    }
  }
  out
}

# ── ODE ──────────────────────────────────────────────────
pkpd_fornari <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    # ── PK (2 compartiments) ──
    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0
    dC1 <- rate_in/V1 - (CL/V1)*C1 - (Q/V1)*C1 + (Q/V2)*C2
    dC2 <- (Q/V1)*C1 - (Q/V2)*C2

    # ── Platine libre → Damage (Eq. 3) ──
    fu_t      <- fu_inf + (fu0 - fu_inf) * exp(-k_bind * time)
    C_free_uM <- fu_t * C1 * mgL_to_uM
    dDamage   <- k_dam * C_free_uM - k_rep * Damage

    # ── Feedbacks (Eq. 11, 12, 13) ──
    eps <- 1e-9

    # fdbk_stem : descendance MPP (Eq. 11)
    r_stem   <- 0.5*(CMP0/max(CMP, eps)) + 0.5*(MEP0/max(MEP, eps))
    f_stem   <- pmin(pmax(r_stem,   0.2), 50)^gamma_stem

    # fdbk_mat_CMP : cellules circulantes blanches (Eq. 12)
    r_matCMP  <- 0.5*(Neut0/max(Neut, eps)) + 0.5*(Mono0/max(Mono, eps))
    f_mat_CMP <- pmin(pmax(r_matCMP, 0.2), 50)^gamma_mat_CMP

    # fdbk_mat_MEP : RBC (Eq. 12)
    r_matMEP  <- RBC0 / max(RBC, eps)
    f_mat_MEP <- pmin(pmax(r_matMEP, 0.2), 50)^gamma_mat_MEP

    # fdbk_prol_Ret : Ret + RBC (Eq. 13)
    r_pRet     <- 0.5*(Ret0/max(Ret, eps)) + 0.5*(RBC0/max(RBC, eps))
    f_prol_Ret <- pmin(pmax(r_pRet, 0.2), 10)^gamma_prolTrans

    # fdbk_prol_Plt : Plt (Eq. 13)
    r_pPlt     <- Plt0 / max(Plt, eps)
    f_prol_Plt <- pmin(pmax(r_pPlt, 0.2), 10)^gamma_prolTrans

    # ── Progeniteurs (Eq. 1, 4) ──
    # Modèle linéaire Fornari : (1 - Slope × Damage)
    # k_dam réduit pour que Slope × Damage_max < 1 au niveau typique
    dMPP <- k_stem * f_stem +
            k_prol_MPP * (1 - Slope_MPP * Damage) * MPP -
            k_tr_CMP * f_mat_CMP * MPP -
            k_tr_MEP * f_mat_MEP * MPP

    dCMP <- k_prol_CMP * (1 - Slope_CMP * Damage) * CMP +
            k_tr_CMP * f_mat_CMP * MPP -
            (k_tr_Neut + k_tr_Mono) * CMP

    dMEP <- k_prol_MEP * (1 - Slope_MEP * Damage) * MEP +
            k_tr_MEP * f_mat_MEP * MPP -
            (k_tr_Ret + k_tr_Plt) * MEP

    # ── Neutrophiles — transit non prolifératif (Eq. 5, Fornari 2019) ──
    # Pas d'effet drogue direct sur les compartiments transit (contrairement à Ret/Plt)
    a_Neut   <- 3 / MTT_Neut
    dT1_Neut <- k_tr_Neut*CMP   - a_Neut*T1_Neut
    dT2_Neut <- a_Neut*T1_Neut  - a_Neut*T2_Neut
    dT3_Neut <- a_Neut*T2_Neut  - a_Neut*T3_Neut
    dNeut    <- a_Neut*T3_Neut  - k_circ_Neut*Neut

    # ── Monocytes — transit non prolifératif (Eq. 5) ──
    a_Mono   <- 3 / MTT_Mono
    dT1_Mono <- k_tr_Mono*CMP     - a_Mono*T1_Mono
    dT2_Mono <- a_Mono*T1_Mono   - a_Mono*T2_Mono
    dT3_Mono <- a_Mono*T2_Mono   - a_Mono*T3_Mono
    dMono    <- a_Mono*T3_Mono   - k_circ_Mono*Mono

    # ── Réticulocytes — T1/T2 prolifératifs (Eq. 7) ──
    # Cap Emax : drug_ret plafonné à k_prol_Ret × f_prol_Ret
    # → taux de prolifération net ≥ 0 (arrest complet au max, pas de mort active)
    drug_ret <- min(delta_Ret * Slope_MEP * k_prol_Ret * Damage,
                    k_prol_Ret * f_prol_Ret)

    dT1_Ret <- k_prol_Ret * f_prol_Ret * T1_Ret -
               drug_ret * T1_Ret +
               k_tr_Ret * MEP - a_Ret * T1_Ret

    dT2_Ret <- k_prol_Ret * f_prol_Ret * T2_Ret -
               drug_ret * T2_Ret +
               a_Ret * T1_Ret - a_Ret * T2_Ret

    dT3_Ret <- a_Ret * T2_Ret - a_Ret * T3_Ret
    dRet    <- a_Ret * T3_Ret - k_circ_Ret * Ret
    dRBC    <- k_circ_Ret * Ret - k_circ_RBC * RBC

    # ── Plaquettes — T1/T2 prolifératifs (Eq. 7) ──
    # Formule correcte (papier Fornari Eq. 7) : drug_plt inclut k_prol_Plt
    # → drug_plt en /h cohérent avec k_prol_Plt (correction du bug e2cc5ec)
    # → drug_ret inclut déjà k_prol_Ret → symétrie maintenant respectée
    drug_plt <- delta_Plt * Slope_MEP * k_prol_Plt * Damage

    dT1_Plt <- k_prol_Plt * f_prol_Plt * T1_Plt -
               drug_plt * T1_Plt +
               k_tr_Plt * MEP - a_Plt * T1_Plt

    dT2_Plt <- k_prol_Plt * f_prol_Plt * T2_Plt -
               drug_plt * T2_Plt +
               a_Plt * T1_Plt - a_Plt * T2_Plt

    dT3_Plt <- a_Plt * T2_Plt - a_Plt * T3_Plt
    dPlt    <- a_Plt * T3_Plt - k_circ_Plt * Plt

    list(c(dC1, dC2, dDamage,
           dMPP, dCMP, dMEP,
           dT1_Neut, dT2_Neut, dT3_Neut, dNeut,
           dT1_Mono, dT2_Mono, dT3_Mono, dMono,
           dT1_Ret,  dT2_Ret,  dT3_Ret,  dRet, dRBC,
           dT1_Plt,  dT2_Plt,  dT3_Plt,  dPlt))
  })
}
