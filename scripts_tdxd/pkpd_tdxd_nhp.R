############################################################
# pkpd_tdxd_nhp.R
# ODE PK/PD — T-DXd (DS-8201a) — SINGE CYNOMOLGUS (NHP)
#
# Modèle (25 états) :
#   PK (5) : ADC 2-cpt + TMDD Michaelis-Menten + DXd 1-cpt
#            + DXd intracellulaire + Dommage ADN
#   PD (20): Modèle Fornari 2019 adapté NHP
#            MPP, CMP, MEP (progéniteurs)
#            Transit + circulants : Neut, Mono, Ret, RBC, Plt
#
# Différences vs pkpd_tdxd_rat.R :
#   → Élimination ADC : CL_lin + TMDD Michaelis-Menten (Vmax/Km)
#     au lieu de CL_lin seul (rat exprime peu HER2)
#   → Paramètres PD : baselines et MTT NHP-spécifiques
#
# Prérequis :
#   source("parameters_tdxd_nhp.R")     # PK NHP
#   source("parameters_tdxd_nhp_pd.R")  # PD NHP
############################################################

pkpd_nhp_ode <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    # ── Perfusion IV ────────────────────────────────────
    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0

    # ── Krel cycle-dépendant (Yin 2020) ─────────────────
    cycle_num <- max(1L, floor(time / interval_h) + 1L)
    krel_t    <- k_rel_c1 * cycle_num^krel_power *
                 ifelse(cycle_num > 1L, krel_factor, 1.0)

    # ════════════════════════════════════════════════════
    # PK — ADC 2-cpt + TMDD Michaelis-Menten
    # ════════════════════════════════════════════════════
    dC_ADC1 <- rate_in / V1_ADC +
               (Q_ADC / V2_ADC) * C_ADC2 -
               (CL_lin / V1_ADC + Q_ADC / V1_ADC + k_int) * C_ADC1 -
               Vmax_MM * C_ADC1 / (Km_MM + C_ADC1)

    dC_ADC2 <- (Q_ADC / V1_ADC) * C_ADC1 -
               (Q_ADC / V2_ADC) * C_ADC2

    # ── PK DXd plasma ───────────────────────────────────
    release  <- mass_frac_DXd * krel_t * C_ADC1 * V1_ADC / V_DXd
    flux_in  <- k_inD  * C_DXd
    flux_out <- k_effD * C_DXd_ic * (V_ic / V_DXd)

    dC_DXd    <- release - (CL_DXd / V_DXd) * C_DXd - flux_in + flux_out

    # ── PK DXd intracellulaire (moelle osseuse) ─────────
    dC_DXd_ic <- k_inD * C_DXd * (V_DXd / V_ic) - k_effD * C_DXd_ic

    # ════════════════════════════════════════════════════
    # DOMMAGE ADN (γH2AX normalisé)
    # ════════════════════════════════════════════════════
    C_DXd_ic_uM <- C_DXd_ic * mgL_to_uM_DXd
    E_drug      <- C_DXd_ic_uM / (IC50_DXd_uM + C_DXd_ic_uM)

    dDamage <- k_dam * E_drug - k_rep * Damage

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
    a_Neut  <- 3 / MTT_Neut
    dT1_Neut <- k_tr_Neut * CMP        - a_Neut * T1_Neut
    dT2_Neut <- a_Neut * T1_Neut       - a_Neut * T2_Neut
    dT3_Neut <- a_Neut * T2_Neut       - a_Neut * T3_Neut
    dNeut    <- a_Neut * T3_Neut       - k_circ_Neut * Neut

    # ── Transit Monocytes ────────────────────────────────
    a_Mono  <- 3 / MTT_Mono
    dT1_Mono <- k_tr_Mono * CMP        - a_Mono * T1_Mono
    dT2_Mono <- a_Mono * T1_Mono       - a_Mono * T2_Mono
    dT3_Mono <- a_Mono * T2_Mono       - a_Mono * T3_Mono
    dMono    <- a_Mono * T3_Mono       - k_circ_Mono * Mono

    # ── Transit Réticulocytes (prolifératifs) ────────────
    f_prol_Ret <- (max(min(RBC0 / max(RBC, 1e-6), 50), 0.2))^gamma_prolTrans
    drug_ret   <- delta_Ret * Slope_MEP * k_prol_Ret * D_kill

    dT1_Ret <- k_prol_Ret * f_prol_Ret * T1_Ret - drug_ret * T1_Ret +
               k_tr_Ret * MEP - a_Ret * T1_Ret
    dT2_Ret <- k_prol_Ret * f_prol_Ret * T2_Ret - drug_ret * T2_Ret +
               a_Ret * T1_Ret - a_Ret * T2_Ret
    dT3_Ret <- a_Ret * T2_Ret                   - a_Ret * T3_Ret
    dRet    <- a_Ret * T3_Ret                   - k_circ_Ret * Ret
    dRBC    <- k_circ_Ret * Ret                 - k_circ_RBC * RBC

    # ── Transit Plaquettes (prolifératifs) ───────────────
    f_prol_Plt <- (max(min(Plt0 / max(Plt, 1e-6), 50), 0.2))^gamma_prolTrans
    drug_plt   <- delta_Plt * Slope_MEP * k_prol_Plt * D_kill

    dT1_Plt <- k_prol_Plt * f_prol_Plt * T1_Plt - drug_plt * T1_Plt +
               k_tr_Plt * MEP - a_Plt * T1_Plt
    dT2_Plt <- k_prol_Plt * f_prol_Plt * T2_Plt - drug_plt * T2_Plt +
               a_Plt * T1_Plt - a_Plt * T2_Plt
    dT3_Plt <- a_Plt * T2_Plt                   - a_Plt * T3_Plt
    dPlt    <- a_Plt * T3_Plt                   - k_circ_Plt * Plt

    list(c(
      dC_ADC1, dC_ADC2, dC_DXd, dC_DXd_ic, dDamage,
      dMPP, dCMP, dMEP,
      dT1_Neut, dT2_Neut, dT3_Neut, dNeut,
      dT1_Mono, dT2_Mono, dT3_Mono, dMono,
      dT1_Ret,  dT2_Ret,  dT3_Ret,  dRet, dRBC,
      dT1_Plt,  dT2_Plt,  dT3_Plt,  dPlt
    ))
  })
}
