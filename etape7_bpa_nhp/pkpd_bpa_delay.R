############################################################
# pkpd_bpa_delay.R
# ODE Fornari PD + PK ADC 2-cmt + SIGNAL TOXIQUE RETARDE
#
# Motivation : la toxicite du BPA (ADC-STING) est DECOUPLEE de la PK
#   (chute a j18-27 pour une dose unique a j0). Le modele Fornari
#   standard (kill ~ Damage instantane) ne peut pas la reproduire.
#
# Solution : le Damage (rapide, suit la PK) alimente une chaine de
#   transit lente Sig1 -> Sig2. Le kill depend de Sig2 (retarde de
#   ~2/k_delay). k_delay controle le delai (a estimer).
#
# ETATS (27) : PK(5) + Damage + Sig1 + Sig2 ... non :
#   C_ADC1,C_ADC2,C_DXd,C_DXd_ic,Damage (5) + Sig1,Sig2 (2) + PD(20) = 27
############################################################
library(deSolve)

pkpd_bpa_delay <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0
    cycle_num <- max(1L, floor(time / interval_h) + 1L)
    krel_t    <- k_rel_c1 * cycle_num^krel_power * ifelse(cycle_num > 1L, krel_factor, 1.0)

    # -- PK ADC 2-cmt --
    dC_ADC1 <- rate_in / V1_ADC + (Q_ADC/V2_ADC)*C_ADC2 -
               (CL_ADC/V1_ADC + Q_ADC/V1_ADC + k_int) * C_ADC1
    dC_ADC2 <- (Q_ADC/V1_ADC)*C_ADC1 - (Q_ADC/V2_ADC)*C_ADC2
    release_ADC <- mass_frac_DXd * krel_t * C_ADC1 * V1_ADC / V_DXd
    dC_DXd    <- release_ADC - (CL_DXd/V_DXd)*C_DXd - k_inD*C_DXd + k_effD*C_DXd_ic*(V_ic/V_DXd)
    dC_DXd_ic <- k_inD*C_DXd*(V_DXd/V_ic) - k_effD*C_DXd_ic

    # -- Effet instantane du drug (suit la PK) --
    use_adc <- if (!is.null(pars$use_ADC_driver)) pars$use_ADC_driver else FALSE
    if (use_adc) E_drug <- C_ADC1/(IC50_ADC_ugmL + C_ADC1)
    else { C_DXd_ic_uM <- C_DXd_ic*mgL_to_uM_DXd; E_drug <- C_DXd_ic_uM/(IC50_DXd_uM + C_DXd_ic_uM) }
    dDamage <- k_dam * E_drug - k_rep * Damage

    # -- SIGNAL TOXIQUE RETARDE : chaine de transit lente --
    kd <- if (!is.null(pars$k_delay)) pars$k_delay else 0.01   # /h
    dSig1 <- kd * (Damage - Sig1)
    dSig2 <- kd * (Sig1   - Sig2)

    # -- Kill pilote par Sig2 (retarde), lineaire --
    D0     <- if (!is.null(pars$Damage_threshold)) pars$Damage_threshold else 0
    D_kill <- pmax(0, Sig2 - D0)
    kill_MPP <- Slope_MPP * D_kill
    kill_CMP <- Slope_CMP * D_kill
    kill_MEP <- Slope_MEP * D_kill

    # -- PD Fornari (structure inchangee) --
    eps <- 1e-9
    r_stem    <- 0.5*(CMP0/max(CMP,eps)) + 0.5*(MEP0/max(MEP,eps)); f_stem<-pmin(pmax(r_stem,0.2),50)^gamma_stem
    r_matCMP  <- 0.5*(Neut0/max(Neut,eps)) + 0.5*(Mono0/max(Mono,eps)); f_mat_CMP<-pmin(pmax(r_matCMP,0.2),50)^gamma_mat_CMP
    r_matMEP  <- RBC0/max(RBC,eps); f_mat_MEP<-pmin(pmax(r_matMEP,0.2),50)^gamma_mat_MEP
    r_pRet    <- 0.5*(Ret0/max(Ret,eps))+0.5*(RBC0/max(RBC,eps)); f_prol_Ret<-pmin(pmax(r_pRet,0.2),10)^gamma_prolTrans
    r_pPlt    <- Plt0/max(Plt,eps); f_prol_Plt<-pmin(pmax(r_pPlt,0.2),10)^gamma_prolTrans

    dMPP <- k_stem*f_stem + k_prol_MPP*(1-kill_MPP)*MPP - k_tr_CMP*f_mat_CMP*MPP - k_tr_MEP*f_mat_MEP*MPP
    dCMP <- k_prol_CMP*(1-kill_CMP)*CMP + k_tr_CMP*f_mat_CMP*MPP - (k_tr_Neut+k_tr_Mono)*CMP
    dMEP <- k_prol_MEP*(1-kill_MEP)*MEP + k_tr_MEP*f_mat_MEP*MPP - (k_tr_Ret+k_tr_Plt)*MEP
    a_Neut<-3/MTT_Neut
    dT1_Neut<-k_tr_Neut*CMP-a_Neut*T1_Neut; dT2_Neut<-a_Neut*T1_Neut-a_Neut*T2_Neut
    dT3_Neut<-a_Neut*T2_Neut-a_Neut*T3_Neut; dNeut<-a_Neut*T3_Neut-k_circ_Neut*Neut
    a_Mono<-3/MTT_Mono
    dT1_Mono<-k_tr_Mono*CMP-a_Mono*T1_Mono; dT2_Mono<-a_Mono*T1_Mono-a_Mono*T2_Mono
    dT3_Mono<-a_Mono*T2_Mono-a_Mono*T3_Mono; dMono<-a_Mono*T3_Mono-k_circ_Mono*Mono
    drug_ret <- delta_Ret*Slope_MEP*k_prol_Ret*D_kill
    dT1_Ret<-k_prol_Ret*f_prol_Ret*T1_Ret-drug_ret*T1_Ret+k_tr_Ret*MEP-a_Ret*T1_Ret
    dT2_Ret<-k_prol_Ret*f_prol_Ret*T2_Ret-drug_ret*T2_Ret+a_Ret*T1_Ret-a_Ret*T2_Ret
    dT3_Ret<-a_Ret*T2_Ret-a_Ret*T3_Ret; dRet<-a_Ret*T3_Ret-k_circ_Ret*Ret; dRBC<-k_circ_Ret*Ret-k_circ_RBC*RBC
    drug_plt <- delta_Plt*Slope_MEP*k_prol_Plt*D_kill
    dT1_Plt<-k_prol_Plt*f_prol_Plt*T1_Plt-drug_plt*T1_Plt+k_tr_Plt*MEP-a_Plt*T1_Plt
    dT2_Plt<-k_prol_Plt*f_prol_Plt*T2_Plt-drug_plt*T2_Plt+a_Plt*T1_Plt-a_Plt*T2_Plt
    dT3_Plt<-a_Plt*T2_Plt-a_Plt*T3_Plt; dPlt<-a_Plt*T3_Plt-k_circ_Plt*Plt

    list(c(dC_ADC1,dC_ADC2,dC_DXd,dC_DXd_ic,dDamage, dSig1,dSig2,
           dMPP,dCMP,dMEP,
           dT1_Neut,dT2_Neut,dT3_Neut,dNeut, dT1_Mono,dT2_Mono,dT3_Mono,dMono,
           dT1_Ret,dT2_Ret,dT3_Ret,dRet,dRBC, dT1_Plt,dT2_Plt,dT3_Plt,dPlt))
  })
}
