############################################################
# pkpd_fgfr2_nhp.R
# ODE PK/PD — T-DXd (DS-8201a) — SINGE CYNOMOLGUS (NHP)
#
# Modèle (23 états) :
#   PK (3) : ADC 2-cpt + TMDD Michaelis-Menten + Damage
#   PD (20): Modèle Fornari 2019 adapté NHP
#            MPP, CMP, MEP (progéniteurs)
#            Transit + circulants : Neut, Mono, Ret, RBC, Plt
#
# Driver PD : concentration ADC plasmatique (C_ADC1) → Damage
#   dDamage/dt = k_dam × C_ADC1 − k_rep × Damage
############################################################

pkpd_nhp_ode <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    # ── Perfusion IV ────────────────────────────────────
    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0

    # ════════════════════════════════════════════════════
    # PK — ADC 2-cpt + TMDD Michaelis-Menten
    # ════════════════════════════════════════════════════
    dC_ADC1 <- rate_in / V1_ADC +
               (Q_ADC / V2_ADC) * C_ADC2 -
               (CL_lin / V1_ADC + Q_ADC / V1_ADC + k_int) * C_ADC1 -
               Vmax_MM * C_ADC1 / (Km_MM + C_ADC1)

    dC_ADC2 <- (Q_ADC / V1_ADC) * C_ADC1 -
               (Q_ADC / V2_ADC) * C_ADC2

    # ════════════════════════════════════════════════════
    # DOMMAGE ADN — driver : C_ADC1 en µM (Fornari-like)
    # dDamage/dt = k_dam × C_ADC1(µM) − k_rep × Damage
    # ════════════════════════════════════════════════════
    C_ADC1_uM <- C_ADC1 * mgL_to_uM_ADC   # mg/L → µM  (÷ MW_ADC × 1000)
    dDamage   <- k_dam * C_ADC1_uM - k_rep * Damage

    # ════════════════════════════════════════════════════
    # PD — Modèle Fornari (kill linéaire sur Damage)
    # ════════════════════════════════════════════════════
    D_kill <- max(0, Damage)

    kill_MPP <- Slope_MPP * D_kill
    kill_CMP <- Slope_CMP * D_kill
    kill_MEP <- Slope_MEP * D_kill

    # ── Feedbacks (Eq. 11-13 Fornari) ───────────────────
    r_stem    <- 0.5 * (CMP0 / max(CMP, 1e-6)) + 0.5 * (MEP0 / max(MEP, 1e-6))
    f_stem    <- (max(min(r_stem,    50), 0.2))^gamma_stem

    r_matCMP  <- 0.5 * (Neut0 / max(Neut, 1e-6)) + 0.5 * (Mono0 / max(Mono, 1e-6))
    f_mat_CMP <- (max(min(r_matCMP,  50), 0.2))^gamma_mat_CMP

    r_matMEP  <- RBC0 / max(RBC, 1e-6)
    f_mat_MEP <- (max(min(r_matMEP,  50), 0.2))^gamma_mat_MEP

    # ── Progéniteurs ────────────────────────────────────
    dMPP <- k_stem * f_stem +
            k_prol_MPP * (1 - kill_MPP) * MPP -
            (k_tr_CMP * f_mat_CMP + k_tr_MEP * f_mat_MEP) * MPP

    dCMP <- k_prol_CMP * (1 - kill_CMP) * CMP +
            k_tr_CMP * f_mat_CMP * MPP -
            (k_tr_Neut + k_tr_Mono) * CMP

    dMEP <- k_prol_MEP * (1 - kill_MEP) * MEP +
            k_tr_MEP * f_mat_MEP * MPP -
            (k_tr_Ret + k_tr_Plt) * MEP

    # ── Transit Neutrophiles ─────────────────────────────
    a_Neut   <- 3 / MTT_Neut
    dT1_Neut <- k_tr_Neut * CMP  - a_Neut * T1_Neut
    dT2_Neut <- a_Neut * T1_Neut - a_Neut * T2_Neut
    dT3_Neut <- a_Neut * T2_Neut - a_Neut * T3_Neut
    dNeut    <- a_Neut * T3_Neut - k_circ_Neut * Neut

    # ── Transit Monocytes ────────────────────────────────
    a_Mono   <- 3 / MTT_Mono
    dT1_Mono <- k_tr_Mono * CMP  - a_Mono * T1_Mono
    dT2_Mono <- a_Mono * T1_Mono - a_Mono * T2_Mono
    dT3_Mono <- a_Mono * T2_Mono - a_Mono * T3_Mono
    dMono    <- a_Mono * T3_Mono - k_circ_Mono * Mono

    # ── Transit Réticulocytes (prolifératifs) ────────────
    f_prol_Ret <- (max(min(RBC0 / max(RBC, 1e-6), 50), 0.2))^gamma_prolTrans
    drug_ret   <- delta_Ret * Slope_MEP * k_prol_Ret * D_kill

    dT1_Ret <- k_prol_Ret * f_prol_Ret * T1_Ret - drug_ret * T1_Ret +
               k_tr_Ret * MEP - a_Ret * T1_Ret
    dT2_Ret <- k_prol_Ret * f_prol_Ret * T2_Ret - drug_ret * T2_Ret +
               a_Ret * T1_Ret - a_Ret * T2_Ret
    dT3_Ret <- a_Ret * T2_Ret - a_Ret * T3_Ret
    dRet    <- a_Ret * T3_Ret - k_circ_Ret * Ret
    dRBC    <- k_circ_Ret * Ret - k_circ_RBC * RBC

    # ── Transit Plaquettes (prolifératifs) ───────────────
    f_prol_Plt <- (max(min(Plt0 / max(Plt, 1e-6), 50), 0.2))^gamma_prolTrans
    drug_plt   <- delta_Plt * Slope_MEP * k_prol_Plt * D_kill

    dT1_Plt <- k_prol_Plt * f_prol_Plt * T1_Plt - drug_plt * T1_Plt +
               k_tr_Plt * MEP - a_Plt * T1_Plt
    dT2_Plt <- k_prol_Plt * f_prol_Plt * T2_Plt - drug_plt * T2_Plt +
               a_Plt * T1_Plt - a_Plt * T2_Plt
    dT3_Plt <- a_Plt * T2_Plt - a_Plt * T3_Plt
    dPlt    <- a_Plt * T3_Plt - k_circ_Plt * Plt

    list(c(
      dC_ADC1, dC_ADC2, dDamage,
      dMPP, dCMP, dMEP,
      dT1_Neut, dT2_Neut, dT3_Neut, dNeut,
      dT1_Mono, dT2_Mono, dT3_Mono, dMono,
      dT1_Ret,  dT2_Ret,  dT3_Ret,  dRet, dRBC,
      dT1_Plt,  dT2_Plt,  dT3_Plt,  dPlt
    ))
  })
}
