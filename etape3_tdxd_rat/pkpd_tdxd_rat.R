############################################################
# pkpd_tdxd_rat.R
# ODE fusionné : T-DXd PK + Fornari PD
#
# STATES (25) :
#   -- PK T-DXd (5) --------------------------------------
#   C_ADC1    ADC compartiment central          [mg/L]
#   C_ADC2    ADC compartiment périphérique     [mg/L]
#   C_DXd     DXd plasma                        [mg/L]
#   C_DXd_ic  DXd intracellulaire (moelle)      [mg/L]
#   Damage    Dommages ADN (γH2AX normalisés)   [-]
#   -- PD Fornari (20) -----------------------------------
#   MPP       Progéniteurs multipotents         [×10⁶/kg]
#   CMP       Progéniteurs myéloïdes            [×10⁶/kg]
#   MEP       Progéniteurs érythro/mégacary.    [×10⁶/kg]
#   T1-T3_Neut, Neut    Transit + Neutrophiles
#   T1-T3_Mono, Mono    Transit + Monocytes
#   T1-T3_Ret,  Ret, RBC Transit + Rétic. + GR
#   T1-T3_Plt,  Plt     Transit + Plaquettes
#
# INTERFACE : seul "Damage" connecte PK et PD
#   PK → Damage via Emax(C_DXd_ic, IC50)
#   Damage → PD via Slope_MPP/CMP/MEP (Fornari Eq. 1,4)
############################################################
library(deSolve)

# -- ODE fusionné ------------------------------------------
pkpd_tdxd_fornari <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    # ════════════════════════════════════════════════════
    # BLOC 1 -- PK T-DXd
    # ════════════════════════════════════════════════════

    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0

    # -- Krel temps-dépendant (Yin 2020) ------------------
    # Krel(cycle) = k_rel_c1 × cycle^krel_power × (krel_factor si cycle > 1)
    cycle_num <- max(1L, floor(time / interval_h) + 1L)
    krel_t    <- k_rel_c1 * cycle_num^krel_power *
                 ifelse(cycle_num > 1L, krel_factor, 1.0)

    # ADC sérum -- 2 compartiments + fuite par internalisation (k_int)
    dC_ADC1 <- rate_in / V1_ADC +
               (Q_ADC / V2_ADC) * C_ADC2 -
               (CL_ADC / V1_ADC + Q_ADC / V1_ADC + k_int) * C_ADC1

    dC_ADC2 <- (Q_ADC / V1_ADC) * C_ADC1 -
               (Q_ADC / V2_ADC) * C_ADC2

    # DXd plasma -- libéré via Krel(t) (empirique, Yin 2020) + efflux cellulaire
    release_ADC   <- mass_frac_DXd * krel_t * C_ADC1 * V1_ADC / V_DXd
    flux_in_cell  <- k_inD  * C_DXd
    flux_out_cell <- k_effD * C_DXd_ic * (V_ic / V_DXd)

    dC_DXd <- release_ADC - (CL_DXd / V_DXd) * C_DXd -
              flux_in_cell + flux_out_cell

    # DXd intracellulaire (moelle osseuse)
    dC_DXd_ic <- k_inD * C_DXd * (V_DXd / V_ic) - k_effD * C_DXd_ic

    # Dommages ADN -- Emax
    # Mode 1 (défaut) : driver = C_DXd_ic [µM]  (Yin 2020, DXd intracell.)
    # Mode 2 (use_ADC_driver=TRUE) : driver = C_ADC1 [µg/mL]  (PFB-10 HPC assay)
    use_adc <- if (!is.null(pars$use_ADC_driver)) pars$use_ADC_driver else FALSE
    if (use_adc) {
      E_drug <- C_ADC1 / (IC50_ADC_ugmL + C_ADC1)
    } else {
      C_DXd_ic_uM <- C_DXd_ic * mgL_to_uM_DXd
      E_drug       <- C_DXd_ic_uM / (IC50_DXd_uM + C_DXd_ic_uM)
    }
    dDamage <- k_dam * E_drug - k_rep * Damage

    # Damage effectif pour le kill (seuil D0 = réparation cellulaire)
    # D_kill = max(0, Damage - D0) :
    #   - si Damage < D0 (trough) → Kill = 0 → CMP récupère entre cycles
    #   - si Damage > D0 (pic)    → Kill proportionnel à l'excès
    # D0 = 0 (défaut, rat DXd-driver) ou D0 = Damage_threshold (humain ADC-driver)
    D0     <- if (!is.null(pars$Damage_threshold)) pars$Damage_threshold else 0
    D_kill <- pmax(0, Damage - D0)

    # -- Kill term : Emax saturant (si ED50_kill défini) ou linéaire (Fornari) --
    # Emax : kill_X = Emax_X × D_kill / (ED50_kill + D_kill)  ∈ [0, Emax_X] ≤ 1
    # → kill ne peut jamais dépasser 1 → pas de prolifération négative
    # → IIV sur Emax_CMP crée distribution de sensibilité → grade 0 possible
    # Linéaire (défaut rat) : kill_X = Slope_X × D_kill  (backward compat.)
    if (!is.null(pars$ED50_kill)) {
      ED50k    <- pars$ED50_kill
      Emx_MPP  <- if (!is.null(pars$Emax_MPP_kill)) pars$Emax_MPP_kill else min(1, Slope_MPP * ED50k)
      Emx_CMP  <- if (!is.null(pars$Emax_CMP_kill)) pars$Emax_CMP_kill else min(1, Slope_CMP * ED50k)
      Emx_MEP  <- if (!is.null(pars$Emax_MEP_kill)) pars$Emax_MEP_kill else min(1, Slope_MEP * ED50k)
      kill_MPP <- pmin(1, Emx_MPP * D_kill / (ED50k + D_kill))
      kill_CMP <- pmin(1, Emx_CMP * D_kill / (ED50k + D_kill))
      kill_MEP <- pmin(1, Emx_MEP * D_kill / (ED50k + D_kill))
    } else {
      # Linéaire Fornari (original) : kill = Slope × D_kill
      # NOTE : Slope_CMP_tdxd_human >> 1 → kill > 1 empiriquement nécessaire
      # pour reproduire G3-4 ~20% FDA (CMP_ss = MPP_input/(Slope×k_prol+k_out))
      kill_MPP <- Slope_MPP * D_kill
      kill_CMP <- Slope_CMP * D_kill
      kill_MEP <- Slope_MEP * D_kill
      kill_MEP <- Slope_MEP * D_kill
    }

    # ════════════════════════════════════════════════════
    # BLOC 2 -- PD Fornari (Equations 1-9 + feedbacks 11-13)
    # Identique à pkpd_model_FORNARI.R -- seul "Damage" change
    # ════════════════════════════════════════════════════
    eps <- 1e-9

    # Feedbacks (Eq. 11, 12, 13)
    r_stem    <- 0.5*(CMP0/max(CMP, eps)) + 0.5*(MEP0/max(MEP, eps))
    f_stem    <- pmin(pmax(r_stem,   0.2), 50)^gamma_stem

    r_matCMP  <- 0.5*(Neut0/max(Neut, eps)) + 0.5*(Mono0/max(Mono, eps))
    f_mat_CMP <- pmin(pmax(r_matCMP, 0.2), 50)^gamma_mat_CMP

    r_matMEP  <- RBC0 / max(RBC, eps)
    f_mat_MEP <- pmin(pmax(r_matMEP, 0.2), 50)^gamma_mat_MEP

    r_pRet    <- 0.5*(Ret0/max(Ret, eps)) + 0.5*(RBC0/max(RBC, eps))
    f_prol_Ret <- pmin(pmax(r_pRet, 0.2), 10)^gamma_prolTrans

    r_pPlt     <- Plt0 / max(Plt, eps)
    f_prol_Plt <- pmin(pmax(r_pPlt, 0.2), 10)^gamma_prolTrans

    # Progéniteurs (Eq. 1, 4) -- kill via kill_MPP/CMP/MEP (Emax ou linéaire)
    dMPP0_kdd <- if (!is.null(pars$k_depl_direct)) pars$k_depl_direct else 0
    dMPP <- k_stem * f_stem +
            k_prol_MPP * (1 - kill_MPP) * MPP -
            k_tr_CMP   * f_mat_CMP * MPP -
            k_tr_MEP   * f_mat_MEP * MPP -
            (if (!is.null(pars$k_depl_direct_MPP)) pars$k_depl_direct_MPP else dMPP0_kdd) * pmax(0, Damage - (if (!is.null(pars$Damage_threshold)) pars$Damage_threshold else 0)) * MPP

    # Deplation DIRECTE optionnelle des progeniteurs (guardee, OFF par defaut) :
    # permet d'atteindre le G4 myeloide observe quand le blocage de prolif seul
    # plafonne a G3. Actif seulement si pars$k_depl_direct est defini.
    kdd      <- if (!is.null(pars$k_depl_direct)) pars$k_depl_direct else 0
    depl_CMP <- kdd * D_kill
    depl_MPP <- (if (!is.null(pars$k_depl_direct_MPP)) pars$k_depl_direct_MPP else kdd) * D_kill

    dCMP <- k_prol_CMP * (1 - kill_CMP) * CMP +
            k_tr_CMP   * f_mat_CMP * MPP -
            (k_tr_Neut + k_tr_Mono) * CMP -
            depl_CMP * CMP

    dMEP <- k_prol_MEP * (1 - kill_MEP) * MEP +
            k_tr_MEP   * f_mat_MEP * MPP -
            (k_tr_Ret + k_tr_Plt) * MEP

    # Neutrophiles -- transit non prolifératif (Eq. 5)
    a_Neut   <- 3 / MTT_Neut
    dT1_Neut <- k_tr_Neut * CMP  - a_Neut * T1_Neut
    dT2_Neut <- a_Neut * T1_Neut - a_Neut * T2_Neut
    dT3_Neut <- a_Neut * T2_Neut - a_Neut * T3_Neut
    dNeut    <- a_Neut * T3_Neut - k_circ_Neut * Neut

    # Monocytes -- transit non prolifératif (Eq. 5)
    a_Mono   <- 3 / MTT_Mono
    dT1_Mono <- k_tr_Mono * CMP  - a_Mono * T1_Mono
    dT2_Mono <- a_Mono * T1_Mono - a_Mono * T2_Mono
    dT3_Mono <- a_Mono * T2_Mono - a_Mono * T3_Mono
    dMono    <- a_Mono * T3_Mono - k_circ_Mono * Mono

    # Réticulocytes -- T1/T2 prolifératifs (Eq. 7)
    drug_ret <- delta_Ret * Slope_MEP * k_prol_Ret * Damage

    dT1_Ret <- k_prol_Ret * f_prol_Ret * T1_Ret - drug_ret * T1_Ret +
               k_tr_Ret  * MEP - a_Ret * T1_Ret
    dT2_Ret <- k_prol_Ret * f_prol_Ret * T2_Ret - drug_ret * T2_Ret +
               a_Ret * T1_Ret - a_Ret * T2_Ret
    dT3_Ret <- a_Ret * T2_Ret - a_Ret * T3_Ret
    dRet    <- a_Ret * T3_Ret - k_circ_Ret * Ret
    dRBC    <- k_circ_Ret * Ret - k_circ_RBC * RBC

    # Plaquettes -- T1/T2 prolifératifs (Eq. 7)
    drug_plt <- delta_Plt * Slope_MEP * k_prol_Plt * Damage

    dT1_Plt <- k_prol_Plt * f_prol_Plt * T1_Plt - drug_plt * T1_Plt +
               k_tr_Plt  * MEP - a_Plt * T1_Plt
    dT2_Plt <- k_prol_Plt * f_prol_Plt * T2_Plt - drug_plt * T2_Plt +
               a_Plt * T1_Plt - a_Plt * T2_Plt
    dT3_Plt <- a_Plt * T2_Plt - a_Plt * T3_Plt
    # Kill direct sur plaquettes circulantes (bypass feedback transit)
    # Slope_Plt_direct > 0 pour sous-groupe susceptible (3.4% FDA pooled G3-4)
    Slope_Plt_d <- if (!is.null(pars$Slope_Plt_direct)) pars$Slope_Plt_direct else 0
    dPlt    <- a_Plt * T3_Plt - k_circ_Plt * Plt - Slope_Plt_d * D_kill * Plt

    list(c(
      # PK
      dC_ADC1, dC_ADC2, dC_DXd, dC_DXd_ic, dDamage,
      # PD
      dMPP, dCMP, dMEP,
      dT1_Neut, dT2_Neut, dT3_Neut, dNeut,
      dT1_Mono, dT2_Mono, dT3_Mono, dMono,
      dT1_Ret,  dT2_Ret,  dT3_Ret,  dRet, dRBC,
      dT1_Plt,  dT2_Plt,  dT3_Plt,  dPlt
    ))
  })
}

# -- Wrapper simulation ------------------------------------
simulate_pkpd_tdxd <- function(times, pars, state0,
                               rtol = 1e-7, atol = 1e-9) {
  out <- as.data.frame(lsoda(
    y        = state0,
    times    = sort(unique(times)),
    func     = pkpd_tdxd_fornari,
    parms    = pars,
    rtol     = rtol,
    atol     = atol,
    maxsteps = 500000
  ))

  out$time_h       <- out$time
  out$time_d       <- out$time / 24
  out$days         <- out$time / 24   # requis par plots.R (save_all_cells)
  out$C_DXd_uM     <- out$C_DXd    * pars$mgL_to_uM_DXd
  out$C_DXd_ic_uM  <- out$C_DXd_ic * pars$mgL_to_uM_DXd

  # Diagnostics PK
  cat(sprintf("  ✓ ADC Cmax       = %.2f mg/L\n",
              max(out$C_ADC1,      na.rm = TRUE)))
  cat(sprintf("  ✓ DXd_ic Cmax    = %.4f µM  (%.1f%% IC50)\n",
              max(out$C_DXd_ic_uM, na.rm = TRUE),
              100 * max(out$C_DXd_ic_uM, na.rm = TRUE) / pars$IC50_DXd_uM))
  cat(sprintf("  ✓ Damage max     = %.4f\n",
              max(out$Damage,      na.rm = TRUE)))

  # AUCtau cycle 1 (trapèze)
  # NOTE : si k_int>0 (humain), CL_eff = CL_ADC + k_int×V1 >> CL_ADC_FDA
  #        → AUCtau modèle < FDA (735 µg·d/mL SS) car k_int non découplé de CL_ADC
  #        Pour rat (k_int=0) : AUCtau modèle ≈ Dose/CL_ADC (valide vs FDA TK)
  ivl <- if (!is.null(pars$interval_h)) pars$interval_h else 504
  c1  <- out[out$time_h <= ivl + 0.5, ]    # légère tolérance numérique
  if (nrow(c1) > 1) {
    dt       <- diff(c1$time_h)
    auc_adc1 <- sum(dt * (head(c1$C_ADC1, -1) + tail(c1$C_ADC1, -1)) / 2) / 24
    auc_dxd1 <- sum(dt * (head(c1$C_DXd,  -1) + tail(c1$C_DXd,  -1)) / 2) * 1000 / 24
    k_int_val <- if (!is.null(pars$k_int)) pars$k_int else 0
    if (k_int_val > 0) {
      cat(sprintf("  ✓ AUCtau ADC C1  = %.0f µg·d/mL  [k_int actif → CL_eff≈%.3fL/h > CL_FDA]\n",
                  auc_adc1, pars$CL_ADC + k_int_val * pars$V1_ADC))
      cat(sprintf("  ✓ AUCtau DXd C1  = %.1f ng·d/mL  (FDA SS=28 non comparable : k_int inclus dans CL_FDA)\n",
                  auc_dxd1))
    } else {
      racc <- 1.3
      cat(sprintf("  ✓ AUCtau ADC C1  = %.0f µg·d/mL  (extrapolé SS×%.1f = %.0f µg·d/mL)\n",
                  auc_adc1, racc, auc_adc1 * racc))
      cat(sprintf("  ✓ AUCtau DXd C1  = %.1f ng·d/mL  (extrapolé SS×%.1f = %.1f ng·d/mL)\n",
                  auc_dxd1, racc, auc_dxd1 * racc))
    }
  }

  # Diagnostics PD -- nadirs
  cat(sprintf("  ✓ Neut nadir     = %.4f  (baseline = %.2f)\n",
              min(out$Neut, na.rm = TRUE), pars$Neut0))
  cat(sprintf("  ✓ Plt  nadir     = %.1f   (baseline = %.0f)\n",
              min(out$Plt,  na.rm = TRUE), pars$Plt0))
  cat(sprintf("  ✓ Ret  nadir     = %.1f   (baseline = %.0f)\n",
              min(out$Ret,  na.rm = TRUE), pars$Ret0))

  # Alerte kill dominant
  D_max <- max(out$Damage, na.rm = TRUE)
  for (sl in c("Slope_MPP", "Slope_CMP", "Slope_MEP")) {
    val <- pars[[sl]]
    if (!is.null(val) && val * D_max > 1)
      cat(sprintf("  !! %s × Damage_max = %.2f > 1 : nadir très profond\n",
                  sl, val * D_max))
  }

  out
}
