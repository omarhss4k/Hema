############################################################
# pkpd_bpa_immuno.R
# Modele IMMUNO-PD pour ADC-STING (BPA)
# Inspire de l'architecture "mediateur soluble" de Chen & D'Argenio
# (AAPS J 2020), en substituant au feedback G-CSF (homeostatique,
# negatif) un SIGNAL STING auto-entretenu (positif, retarde).
#
# Signal S :
#   dS/dt = k_trig * E_adc            (declenchement bref par l'ADC)
#         + k_amp  * S * (1 - S/Smax) (AUTO-AMPLIFICATION, persiste sans ADC)
#         - k_off  * S                (resolution lente)
#   E_adc = C_ADC1 / (EC50_ADC + C_ADC1)
#
# Deux effets de S (les deux directions) :
#   Stimulation (rebond)  : relargage matures x (1 + Emax_stim*S/(EC50s+S))
#     -> EC50s petit -> sature tot -> domine quand S est modere (precoce)
#   Suppression (kill)    : kill_X = Slope_X * S   (progeniteurs)
#     -> lineaire -> domine quand S est grand (tardif)
#
# -> rebond precoce PUIS effondrement retarde emergent de la dynamique de S.
#
# ETATS (24) : C_ADC1,C_ADC2, S, MPP,CMP,MEP, T1-3_Neut,Neut,
#   T1-3_Mono,Mono, T1-3_Ret,Ret,RBC, T1-3_Plt,Plt
############################################################
library(deSolve)

pkpd_bpa_immuno <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0

    # -- PK ADC 2-cmt --
    dC_ADC1 <- rate_in/V1_ADC + (Q_ADC/V2_ADC)*C_ADC2 -
               (CL_ADC/V1_ADC + Q_ADC/V1_ADC)*C_ADC1
    dC_ADC2 <- (Q_ADC/V1_ADC)*C_ADC1 - (Q_ADC/V2_ADC)*C_ADC2

    # -- DEUX MEDIATEURS immunitaires (aigu + chronique) --
    E_adc <- C_ADC1 / (EC50_ADC + C_ADC1)
    # Sa : reponse AIGUE -- montee/descente rapides (cytokines inflammatoires)
    dSa <- k_a*E_adc - koff_a*Sa
    # Sc : reponse CHRONIQUE -- declenchee par Sa, auto-amplifiee, tardive
    dSc <- k_c*Sa*(1 - Sc/Smax) - koff_c*Sc

    stim <- Emax_stim * Sa/(EC50_stim + Sa)   # rebond -> pilote par l'AIGU (Sa)
    depl_CMP <- k_depl_CMP * Sc                # effondrement -> pilote par le CHRONIQUE (Sc)
    depl_MEP <- k_depl_MEP * Sc
    depl_MPP <- k_depl_MPP * Sc

    # -- PD (structure Fornari, feedbacks conserves) --
    eps <- 1e-9
    r_stem   <- 0.5*(CMP0/max(CMP,eps))+0.5*(MEP0/max(MEP,eps)); f_stem<-pmin(pmax(r_stem,0.2),50)^gamma_stem
    r_matCMP <- 0.5*(Neut0/max(Neut,eps))+0.5*(Mono0/max(Mono,eps)); f_mat_CMP<-pmin(pmax(r_matCMP,0.2),50)^gamma_mat_CMP
    r_matMEP <- RBC0/max(RBC,eps); f_mat_MEP<-pmin(pmax(r_matMEP,0.2),50)^gamma_mat_MEP
    r_pRet   <- 0.5*(Ret0/max(Ret,eps))+0.5*(RBC0/max(RBC,eps)); f_prol_Ret<-pmin(pmax(r_pRet,0.2),10)^gamma_prolTrans
    r_pPlt   <- Plt0/max(Plt,eps); f_prol_Plt<-pmin(pmax(r_pPlt,0.2),10)^gamma_prolTrans

    dMPP <- k_stem*f_stem + k_prol_MPP*MPP - k_tr_CMP*f_mat_CMP*MPP - k_tr_MEP*f_mat_MEP*MPP - depl_MPP*MPP
    dCMP <- k_prol_CMP*CMP + k_tr_CMP*f_mat_CMP*MPP - (k_tr_Neut+k_tr_Mono)*CMP - depl_CMP*CMP
    dMEP <- k_prol_MEP*MEP + k_tr_MEP*f_mat_MEP*MPP - (k_tr_Ret+k_tr_Plt)*MEP - depl_MEP*MEP

    # Neutrophiles : relargage terminal stimule par S (rebond = demargination)
    a_Neut<-3/MTT_Neut
    dT1_Neut<-k_tr_Neut*CMP-a_Neut*T1_Neut; dT2_Neut<-a_Neut*T1_Neut-a_Neut*T2_Neut
    dT3_Neut<-a_Neut*T2_Neut-a_Neut*(1+stim)*T3_Neut
    dNeut   <-a_Neut*(1+stim)*T3_Neut - k_circ_Neut*Neut
    a_Mono<-3/MTT_Mono
    dT1_Mono<-k_tr_Mono*CMP-a_Mono*T1_Mono; dT2_Mono<-a_Mono*T1_Mono-a_Mono*T2_Mono
    dT3_Mono<-a_Mono*T2_Mono-a_Mono*T3_Mono; dMono<-a_Mono*T3_Mono-k_circ_Mono*Mono

    # Reticulocytes : relargage stimule par S (reticulocytose de stress)
    dT1_Ret<-k_prol_Ret*f_prol_Ret*T1_Ret+k_tr_Ret*MEP-a_Ret*T1_Ret
    dT2_Ret<-k_prol_Ret*f_prol_Ret*T2_Ret+a_Ret*T1_Ret-a_Ret*T2_Ret
    dT3_Ret<-a_Ret*T2_Ret-a_Ret*(1+stim)*T3_Ret
    dRet   <-a_Ret*(1+stim)*T3_Ret - k_circ_Ret*Ret
    dRBC   <-k_circ_Ret*Ret - k_circ_RBC*RBC

    dT1_Plt<-k_prol_Plt*f_prol_Plt*T1_Plt+k_tr_Plt*MEP-a_Plt*T1_Plt
    dT2_Plt<-k_prol_Plt*f_prol_Plt*T2_Plt+a_Plt*T1_Plt-a_Plt*T2_Plt
    dT3_Plt<-a_Plt*T2_Plt-a_Plt*T3_Plt; dPlt<-a_Plt*T3_Plt-k_circ_Plt*Plt

    list(c(dC_ADC1,dC_ADC2,dSa,dSc, dMPP,dCMP,dMEP,
           dT1_Neut,dT2_Neut,dT3_Neut,dNeut, dT1_Mono,dT2_Mono,dT3_Mono,dMono,
           dT1_Ret,dT2_Ret,dT3_Ret,dRet,dRBC, dT1_Plt,dT2_Plt,dT3_Plt,dPlt),
         Sa=Sa, Sc=Sc, E_adc=E_adc)
  })
}
