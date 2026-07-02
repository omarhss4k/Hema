############################################################
# pkpd_emax_ic50_test.R  --  VARIANTE EXPERIMENTALE (test de faisabilite)
#
# But : tester si l'IC50 CFU utilisee DIRECTEMENT comme EC50 d'un effet
#       Emax par compartiment reproduit les grades FDA de l'etape 5.
#
# Mecanisme (vs Slope*Damage lineaire de l'ODE validee) :
#   Cnm    = C_ADC1[ug/mL] * 1e6 / MW_ADC          (conversion -> nM)
#   E_X    = Cnm / (IC50_X_nM + Cnm)               (Emax=1, EC50 = IC50 CFU)
#   deplet = kmax_X * E_X * X   (DEPLETION DIRECTE, pas inhibition de prolif)
#
# -> l'IC50 mesuree devient l'EC50 de l'effet sur ce compartiment.
# -> kmax_X (taux de mort max, 1/h) porte l'AMPLITUDE (non fournie par l'IC50).
#
# Base : pkpd_tdxd_fornari (2-cmt humain, 25 etats), structure PD INCHANGEE.
############################################################

pkpd_emax_ic50 <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0
    cycle_num <- max(1L, floor(time / interval_h) + 1L)
    krel_t    <- k_rel_c1 * cycle_num^krel_power * ifelse(cycle_num > 1L, krel_factor, 1.0)

    # PK ADC 2-cmt
    dC_ADC1 <- rate_in / V1_ADC + (Q_ADC/V2_ADC)*C_ADC2 -
               (CL_ADC/V1_ADC + Q_ADC/V1_ADC + k_int) * C_ADC1
    dC_ADC2 <- (Q_ADC/V1_ADC)*C_ADC1 - (Q_ADC/V2_ADC)*C_ADC2
    release_ADC   <- mass_frac_DXd * krel_t * C_ADC1 * V1_ADC / V_DXd
    dC_DXd    <- release_ADC - (CL_DXd/V_DXd)*C_DXd - k_inD*C_DXd + k_effD*C_DXd_ic*(V_ic/V_DXd)
    dC_DXd_ic <- k_inD*C_DXd*(V_DXd/V_ic) - k_effD*C_DXd_ic
    dDamage   <- 0   # inutilise dans cette variante (garde l'etat pour compat.)

    # ---- EFFET Emax-IC50 par compartiment (deplation directe) --------------
    Cnm   <- C_ADC1 * 1e6 / MW_ADC              # ug/mL -> nM
    E_CMP <- Cnm / (IC50_CMP_nM + Cnm)
    E_MEP <- Cnm / (IC50_MEP_nM + Cnm)
    E_MPP <- Cnm / (IC50_MPP_nM + Cnm)
    depl_MPP <- kmax_MPP * E_MPP
    depl_CMP <- kmax_CMP * E_CMP
    depl_MEP <- kmax_MEP * E_MEP

    # ---- PD Fornari (structure INCHANGEE) ----------------------------------
    eps <- 1e-9
    r_stem    <- 0.5*(CMP0/max(CMP,eps)) + 0.5*(MEP0/max(MEP,eps))
    f_stem    <- pmin(pmax(r_stem,0.2),50)^gamma_stem
    r_matCMP  <- 0.5*(Neut0/max(Neut,eps)) + 0.5*(Mono0/max(Mono,eps))
    f_mat_CMP <- pmin(pmax(r_matCMP,0.2),50)^gamma_mat_CMP
    r_matMEP  <- RBC0/max(RBC,eps)
    f_mat_MEP <- pmin(pmax(r_matMEP,0.2),50)^gamma_mat_MEP
    r_pRet    <- 0.5*(Ret0/max(Ret,eps)) + 0.5*(RBC0/max(RBC,eps))
    f_prol_Ret <- pmin(pmax(r_pRet,0.2),10)^gamma_prolTrans
    r_pPlt     <- Plt0/max(Plt,eps)
    f_prol_Plt <- pmin(pmax(r_pPlt,0.2),10)^gamma_prolTrans

    # Progeniteurs : prolif NORMALE + deplation directe (-depl_X * X)
    dMPP <- k_stem*f_stem + k_prol_MPP*MPP - k_tr_CMP*f_mat_CMP*MPP -
            k_tr_MEP*f_mat_MEP*MPP - depl_MPP*MPP
    dCMP <- k_prol_CMP*CMP + k_tr_CMP*f_mat_CMP*MPP -
            (k_tr_Neut+k_tr_Mono)*CMP - depl_CMP*CMP
    dMEP <- k_prol_MEP*MEP + k_tr_MEP*f_mat_MEP*MPP -
            (k_tr_Ret+k_tr_Plt)*MEP - depl_MEP*MEP

    a_Neut <- 3/MTT_Neut
    dT1_Neut <- k_tr_Neut*CMP - a_Neut*T1_Neut
    dT2_Neut <- a_Neut*T1_Neut - a_Neut*T2_Neut
    dT3_Neut <- a_Neut*T2_Neut - a_Neut*T3_Neut
    dNeut    <- a_Neut*T3_Neut - k_circ_Neut*Neut
    a_Mono <- 3/MTT_Mono
    dT1_Mono <- k_tr_Mono*CMP - a_Mono*T1_Mono
    dT2_Mono <- a_Mono*T1_Mono - a_Mono*T2_Mono
    dT3_Mono <- a_Mono*T2_Mono - a_Mono*T3_Mono
    dMono    <- a_Mono*T3_Mono - k_circ_Mono*Mono
    # Deplation des precurseurs erythro/mega PROLIFERANTS (T1/T2), meme EC50 (E_MEP)
    #   -> reproduit l'atteinte des progeniteurs en division (cf drug_ret etape5),
    #      necessaire car la deplation du MEP seul est tamponnee (RBC vie longue).
    depl_retp <- kmax_ret_prol * E_MEP
    depl_pltp <- kmax_plt_prol * E_MEP
    dT1_Ret <- k_prol_Ret*f_prol_Ret*T1_Ret - depl_retp*T1_Ret + k_tr_Ret*MEP - a_Ret*T1_Ret
    dT2_Ret <- k_prol_Ret*f_prol_Ret*T2_Ret - depl_retp*T2_Ret + a_Ret*T1_Ret - a_Ret*T2_Ret
    dT3_Ret <- a_Ret*T2_Ret - a_Ret*T3_Ret
    dRet    <- a_Ret*T3_Ret - k_circ_Ret*Ret
    dRBC    <- k_circ_Ret*Ret - k_circ_RBC*RBC
    dT1_Plt <- k_prol_Plt*f_prol_Plt*T1_Plt - depl_pltp*T1_Plt + k_tr_Plt*MEP - a_Plt*T1_Plt
    dT2_Plt <- k_prol_Plt*f_prol_Plt*T2_Plt - depl_pltp*T2_Plt + a_Plt*T1_Plt - a_Plt*T2_Plt
    dT3_Plt <- a_Plt*T2_Plt - a_Plt*T3_Plt
    dPlt    <- a_Plt*T3_Plt - k_circ_Plt*Plt

    list(c(dC_ADC1, dC_ADC2, dC_DXd, dC_DXd_ic, dDamage,
           dMPP, dCMP, dMEP,
           dT1_Neut,dT2_Neut,dT3_Neut,dNeut,
           dT1_Mono,dT2_Mono,dT3_Mono,dMono,
           dT1_Ret,dT2_Ret,dT3_Ret,dRet,dRBC,
           dT1_Plt,dT2_Plt,dT3_Plt,dPlt))
  })
}
