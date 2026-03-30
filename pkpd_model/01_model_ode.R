# =============================================================================
# 01_model_ode.R
# Fornari 2019 - Modèle PKPD Carboplatine / Toxicité Hématologique
#
# Système d'équations différentielles complet :
#   - PK 2 compartiments (platine libre)
#   - Dommage ADN (kdam / krep)
#   - Progéniteurs : MPP, CMP, MEP
#   - Compartiments de transit (non-prolifératifs) : Neut, Mono (3 x chacun)
#   - Compartiments de transit (prolifératifs)     : Ret, Plt  (2 prolif + 1 maturation)
#   - Cellules circulantes : Neut, Mono, Ret, Plt, RBC
#
# Référence : Fornari et al. (2019) CPT:PSP
# Unités    : temps en heures, concentrations en µM, cellules en 10^9/L
# =============================================================================

#' Système ODE du modèle PKPD Fornari
#'
#' @param t  Temps courant (h)
#' @param y  Vecteur d'état (25 composantes, voir ci-dessous)
#' @param p  Liste de paramètres (issue de build_parameters())
#' @param dose_fun Fonction dose(t) → quantité de platine injectée (µmol) dans Cen
#'
#' @return Liste avec dydt (vecteur 25)
#'
#' Index du vecteur d'état :
#'  PK         : 1=Cen,   2=Per
#'  Dommage    : 3=Damage
#'  Progén.    : 4=MPP,   5=CMP,   6=MEP
#'  Transit Neut : 7=T1N, 8=T2N,  9=T3N
#'  Transit Mono : 10=T1Mo,11=T2Mo,12=T3Mo
#'  Transit Ret  : 13=T1R, 14=T2R, 15=T3R
#'  Transit Plt  : 16=T1P, 17=T2P, 18=T3P
#'  Circulants   : 19=Neut,20=Mono,21=Ret,22=Plt
#'  RBC          : 23=RBC
#'  (réservé IIV): 24-25 non utilisés ici
pkpd_ode <- function(t, y, p, dose_fun = NULL) {

  # --- Déballage du vecteur d'état ---
  Cen    <- max(y[1],  0)
  Per    <- max(y[2],  0)
  Damage <- max(y[3],  0)

  MPP    <- max(y[4],  1e-6)
  CMP    <- max(y[5],  1e-6)
  MEP    <- max(y[6],  1e-6)

  T1N    <- max(y[7],  0)
  T2N    <- max(y[8],  0)
  T3N    <- max(y[9],  0)

  T1Mo   <- max(y[10], 0)
  T2Mo   <- max(y[11], 0)
  T3Mo   <- max(y[12], 0)

  T1R    <- max(y[13], 0)
  T2R    <- max(y[14], 0)
  T3R    <- max(y[15], 0)

  T1P    <- max(y[16], 0)
  T2P    <- max(y[17], 0)
  T3P    <- max(y[18], 0)

  Neut   <- max(y[19], 1e-6)
  Mono   <- max(y[20], 1e-6)
  Ret    <- max(y[21], 1e-6)
  Plt    <- max(y[22], 1e-6)
  RBC    <- max(y[23], 1e-6)

  # --- Concentration plasmatique de platine libre (µM) ---
  Cp <- Cen / p$VCen   # µmol / L = µM

  # --- Terme de bolus IV (dose en µmol/h injectée dans Cen) ---
  bolus_rate <- if (!is.null(dose_fun)) dose_fun(t) else 0

  # =========================================================================
  # 1. PHARMACOCINÉTIQUE  (modèle 2 compartiments, platine libre)
  #    Équation S1 — Fornari 2019
  # =========================================================================
  dCen <- bolus_rate - (p$CL / p$VCen) * Cen - (p$QC / p$VCen) * Cen + (p$QC / p$VPer) * Per
  dPer <- (p$QC / p$VCen) * Cen - (p$QC / p$VPer) * Per

  # =========================================================================
  # 2. DOMMAGE ADN
  #    dDamage/dt = kdam * Cp - krep * Damage   (Équation 3)
  # =========================================================================
  dDamage <- p$kdam * Cp - p$krep * Damage

  # =========================================================================
  # 3. FONCTIONS DE FEEDBACK  (Équations 11-13)
  # =========================================================================

  # a) Feedback sur le flux entrant dans MPP (depuis les HSC)
  #    fdbkstem = (0.5 * CMP0/CMP + 0.5 * MEP0/MEP)^γstem
  fdbk_stem <- (0.5 * p$CMP0 / CMP + 0.5 * p$MEP0 / MEP)^p$gamma_stem

  # b) Feedback sur la maturation MPP → CMP (signaux myeloïdes)
  #    fdbkmatCMP = (0.5 * Neut0/Neut + 0.5 * Mono0/Mono)^γmatCMP
  fdbk_mat_CMP <- (0.5 * p$Neut0 / Neut + 0.5 * p$Mono0 / Mono)^p$gamma_mat_CMP

  # c) Feedback sur la maturation MPP → MEP (signaux érythroïdes)
  #    fdbkmatMEP = (RBC0/RBC)^γmatMEP
  fdbk_mat_MEP <- (p$RBC0 / RBC)^p$gamma_mat_MEP

  # d) Feedback prolifération Ret (érythropoïétine)
  #    fdbkprolRet = (0.5 * Ret0/Ret + 0.5 * RBC0/RBC)^γprolTrans
  fdbk_prol_Ret <- (0.5 * p$Ret0 / Ret + 0.5 * p$RBC0 / RBC)^p$gamma_prol_trans

  # e) Feedback prolifération Plt (thrombopoïétine)
  #    fdbkprolPlt = (Plt0/Plt)^γprolTrans
  fdbk_prol_Plt <- (p$Plt0 / Plt)^p$gamma_prol_trans

  # =========================================================================
  # 4. TAUX DE MATURATION  (Équations 6, 8)
  #    a = 3/MTT  (3 compartiments de transit)
  # =========================================================================
  aNeut <- 3 / p$MTTNeut
  aMono <- 3 / p$MTTMono
  aRet  <- 3 / p$MTTRet
  aPlt  <- 3 / p$MTTPlt

  # =========================================================================
  # 5. EFFETS DU MÉDICAMENT SUR LA PROLIFÉRATION
  #    Inhibition linéaire : (1 - Slope * Damage)
  #    Clampé à 0 pour éviter des valeurs négatives
  # =========================================================================
  drug_eff_MPP <- max(0, 1 - p$SlopeMPP * Damage)
  drug_eff_CMP <- max(0, 1 - p$SlopeCMP * Damage)
  drug_eff_MEP <- max(0, 1 - p$SlopeMEP * Damage)

  # =========================================================================
  # 6. PROGÉNITEURS  (Équations 1, 4)
  # =========================================================================

  # MPP (Multipotent Progenitor)
  # dMPP/dt = kstem*fdbkstem + kprolMPP*(1 - SlopeMPP*Damage)*MPP
  #           - ktrCMP*fdbkmatCMP*MPP - ktrMEP*fdbkmatMEP*MPP
  dMPP <- p$kstem * fdbk_stem +
          p$kprolMPP * drug_eff_MPP * MPP -
          p$ktrCMP * fdbk_mat_CMP * MPP -
          p$ktrMEP * fdbk_mat_MEP * MPP

  # CMP (Common Myeloid Progenitor)
  # dCMP/dt = kprolCMP*(1-SlopeCMP*Damage)*CMP + ktrCMP*fdbkmatCMP*MPP
  #           - ktrNeut*CMP - ktrMono*CMP
  dCMP <- p$kprolCMP * drug_eff_CMP * CMP +
          p$ktrCMP * fdbk_mat_CMP * MPP -
          p$ktrNeut * CMP -
          p$ktrMono * CMP

  # MEP (Megakaryocyte-Erythrocyte Progenitor)
  # dMEP/dt = kprolMEP*(1-SlopeMEP*Damage)*MEP + ktrMEP*fdbkmatMEP*MPP
  #           - ktrRet*MEP - ktrPlt*MEP
  dMEP <- p$kprolMEP * drug_eff_MEP * MEP +
          p$ktrMEP * fdbk_mat_MEP * MPP -
          p$ktrRet * MEP -
          p$ktrPlt * MEP

  # =========================================================================
  # 7. TRANSIT NON-PROLIFÉRATIF : Neutrophiles  (Équation 5)
  # =========================================================================
  dT1N <- p$ktrNeut * CMP - aNeut * T1N
  dT2N <- aNeut * T1N     - aNeut * T2N
  dT3N <- aNeut * T2N     - aNeut * T3N

  # =========================================================================
  # 8. TRANSIT NON-PROLIFÉRATIF : Monocytes  (Équation 5)
  # =========================================================================
  dT1Mo <- p$ktrMono * CMP - aMono * T1Mo
  dT2Mo <- aMono * T1Mo    - aMono * T2Mo
  dT3Mo <- aMono * T2Mo    - aMono * T3Mo

  # =========================================================================
  # 9. TRANSIT PROLIFÉRATIF : Réticulocytes  (Équation 7)
  #    Les compartiments T1R et T2R prolifèrent et sont sensibles au médicament
  #    via le terme δRet * SlopeMEP * Damage
  # =========================================================================
  drug_kill_Ret <- p$delta_Ret * p$SlopeMEP * Damage

  dT1R <- p$kprolRet * fdbk_prol_Ret * T1R -
          drug_kill_Ret * p$kprolRet * T1R +
          p$ktrRet * MEP -
          aRet * T1R

  dT2R <- p$kprolRet * fdbk_prol_Ret * T2R -
          drug_kill_Ret * p$kprolRet * T2R +
          aRet * T1R -
          aRet * T2R

  dT3R <- aRet * T2R - aRet * T3R

  # =========================================================================
  # 10. TRANSIT PROLIFÉRATIF : Plaquettes  (Équation 7)
  # =========================================================================
  drug_kill_Plt <- p$delta_Plt * p$SlopeMEP * Damage

  dT1P <- p$kprolPlt * fdbk_prol_Plt * T1P -
          drug_kill_Plt * p$kprolPlt * T1P +
          p$ktrPlt * MEP -
          aPlt * T1P

  dT2P <- p$kprolPlt * fdbk_prol_Plt * T2P -
          drug_kill_Plt * p$kprolPlt * T2P +
          aPlt * T1P -
          aPlt * T2P

  dT3P <- aPlt * T2P - aPlt * T3P

  # =========================================================================
  # 11. CELLULES CIRCULANTES  (Équation 9)
  # =========================================================================
  dNeut <- aNeut * T3N  - p$kcircNeut * Neut
  dMono <- aMono * T3Mo - p$kcircMono * Mono
  dRet  <- aRet  * T3R  - p$kcircRet  * Ret
  dPlt  <- aPlt  * T3P  - p$kcircPlt  * Plt

  # RBC : production depuis les réticulocytes
  dRBC  <- p$kcircRet * Ret - p$kcircRBC * RBC

  # =========================================================================
  # Assemblage du vecteur de dérivées
  # =========================================================================
  dydt <- c(
    dCen,                         # 1
    dPer,                         # 2
    dDamage,                      # 3
    dMPP, dCMP, dMEP,             # 4-6
    dT1N, dT2N, dT3N,             # 7-9
    dT1Mo, dT2Mo, dT3Mo,          # 10-12
    dT1R, dT2R, dT3R,             # 13-15
    dT1P, dT2P, dT3P,             # 16-18
    dNeut, dMono, dRet, dPlt,     # 19-22
    dRBC                          # 23
  )

  list(dydt)
}


#' Noms des variables d'état (pour débogage / plots)
state_names <- c(
  "Cen", "Per", "Damage",
  "MPP", "CMP", "MEP",
  "T1Neut", "T2Neut", "T3Neut",
  "T1Mono", "T2Mono", "T3Mono",
  "T1Ret",  "T2Ret",  "T3Ret",
  "T1Plt",  "T2Plt",  "T3Plt",
  "Neut", "Mono", "Ret", "Plt", "RBC"
)
