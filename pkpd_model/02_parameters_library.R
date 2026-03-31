# =============================================================================
# 02_parameters_library.R
# Dictionnaire de paramètres — Modèle PKPD Carboplatine (Fornari 2019)
#
# Contenu :
#   1. Paramètres fixes humains (physiologie + drug)
#   2. Valeurs basales des cellules (état initial)
#   3. Paramètres dérivés à l'état stationnaire (calculés analytiquement)
#   4. Structure IIV (Sigma Table S4 — log-additive residual variability)
#   5. Fonction principale : build_parameters(species, overrides)
#
# Références :
#   - Fornari et al. (2019) CPT:PSP — Table 1, Table 2, Table S2, Table S3, Table S4
#   - De Carlo et al. (2025) Br J Clin Pharmacol — IIV Friberg (reference secondaire)
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 1 — PARAMÈTRES FIXES RAT  (Table 1 & 2, Fornari 2019)
# -----------------------------------------------------------------------------
.rat_physiology <- list(

  # Taux de clairance circulatoire (1/h)
  kcircNeut = 0.17,
  kcircMono = 0.06,
  kcircPlt  = 0.01,
  kcircRBC  = 0.0007,

  # Valeurs basales circulantes (10^9 cells/L)
  Neut0 = 1.32,
  Mono0 = 0.16,
  Plt0  = 967,
  Ret0  = 235,
  RBC0  = 8161,

  # Valeurs basales progéniteurs (10^9 cells/L)
  MPP0  = 69,
  CMP0  = 214,
  MEP0  = 156,

  # Temps de transit moyen (h)
  MTTNeut = 61.5,
  MTTMono = 62.2,
  MTTPlt  = 74,
  MTTRet  = 75,

  # Paramètres de reparamétrisation (sans unité)
  # lambda1 = T2Ret0/T1Ret0,  lambda2 = T2Plt0/T1Plt0
  # lambda3 = (ktrNeut+ktrMono)/kprolCMP, etc.
  lambda1 = 2.0,
  lambda2 = 2.0,
  lambda3 = 1.8,
  lambda4 = 1.8,
  lambda5 = 1.8,

  # Exposants de feedback (sans unité)
  gamma_stem      = 0.07,
  gamma_mat_CMP   = 0.60,
  gamma_mat_MEP   = 0.30,
  gamma_prol_trans = 0.70,

  # Dommage ADN (1/h)
  kdam = 0.017,
  krep = 0.017,

  # Paramètres drogue (rat) — Table 2
  SlopeMPP  = 2.05,   # (1/µM)
  SlopeCMP  = 1.47,   # (1/µM)
  SlopeMEP  = 2.19,   # (1/µM)
  delta_Ret = 2.80,   # sans unité
  delta_Plt = 0.54,   # sans unité

  # PK Carboplatine rat — Table S3
  VCen = 0.26,   # L
  VPer = 0.30,   # L
  CL   = 0.42,   # L/h
  QC   = 0.06    # L/h
)

# -----------------------------------------------------------------------------
# SECTION 2 — PARAMÈTRES HUMAINS (Table 1 — Fornari 2019 + allométrie)
# -----------------------------------------------------------------------------
.human_physiology <- list(

  # Taux de clairance circulatoire (1/h)
  kcircNeut = 0.10,
  kcircMono = 0.04,
  kcircPlt  = 0.0052,
  kcircRBC  = 0.00037,

  # Valeurs basales circulantes (10^9 cells/L) — Tableau 1 humain
  Neut0 = 4.50,    # médiane population adulte saine
  Mono0 = 0.50,
  Plt0  = 345.0,
  Ret0  = 77.5,    # médiane (40-115)
  RBC0  = 5000.0,  # médiane (4100-5900)

  # Valeurs basales progéniteurs humains (allométrie + proportions Table 1)
  # MPP = α * M^(3/4),  log10(α)=7.46, M=70 kg → MPP ≈ 1.3 × 10^9/L
  MPP0  = 1.3,
  CMP0  = 20.9,    # CMP0/MPP0 = 16
  MEP0  = 14.3,    # MEP0/MPP0 = 11

  # Temps de transit moyen humains (h) — Table 1
  MTTNeut = 210.0,
  MTTMono = 121.5,
  MTTPlt  = 168.0,
  MTTRet  = 66.0,

  # Reparamétrisation (identiques rat — fixes)
  lambda1 = 2.0,
  lambda2 = 2.0,
  lambda3 = 1.8,
  lambda4 = 1.8,
  lambda5 = 1.8,

  # Exposants de feedback (identiques rat — fixes)
  gamma_stem       = 0.07,
  gamma_mat_CMP    = 0.60,
  gamma_mat_MEP    = 0.30,
  gamma_prol_trans = 0.70,

  # Dommage ADN (1/h) — identiques rat
  kdam = 0.017,
  krep = 0.017,

  # Paramètres drogue humains — Table S2 (scaling IC50 rat→human)
  # SlopeH = SlopeR * (IC50_R / IC50_H)
  # IC50_CMP: rat=2.83, human=4.04 → facteur = 2.83/4.04 = 0.701
  # IC50_MEP: rat=0.86, human=1.56 → facteur = 0.86/1.56 = 0.551
  SlopeMPP  = 0.79,   # (1/µM)
  SlopeCMP  = 0.57,   # (1/µM) — 1.47 * 0.701 ≈ 1.03 → valeur Fornari Table S2
  SlopeMEP  = 0.66,   # (1/µM) — 2.19 * 0.551 ≈ 1.21 → valeur Fornari Table S2
  delta_Ret = 2.80,   # inchangé
  delta_Plt = 0.54,   # inchangé

  # PK Carboplatine humain — modèle Zandvliet 2008 / De Carlo 2025
  # Paramètres population (fixés; variabilité gérée dans simulation_core)
  VCen = 15.5,    # L  (Zandvliet 2008)
  VPer = 9.86,    # L
  CL   = NA,      # L/h — calculé dynamiquement via Calvert/Cockcroft-Gault
  QC   = 3.46     # L/h
)

# -----------------------------------------------------------------------------
# SECTION 3 — CALCUL DES PARAMÈTRES À L'ÉTAT STATIONNAIRE
#             (Équations S3, S4, S5-S7 — Fornari 2019)
#
#  Logique : à l'état stationnaire (dX/dt = 0, Damage = 0),
#  on résout algébriquement tous les taux de prolifération et de transition.
# -----------------------------------------------------------------------------

#' Calcule les paramètres dynamiques à partir des valeurs basales
#'
#' @param p  Liste partielle de paramètres (physio + drug)
#' @return   Liste complète avec kprol*, ktr*, kstem ajoutés
.compute_steady_state <- function(p) {

  # Taux de maturation
  aNeut <- 3 / p$MTTNeut
  aMono <- 3 / p$MTTMono
  aRet  <- 3 / p$MTTRet
  aPlt  <- 3 / p$MTTPlt

  # kcircRet : conservation de masse à l'état stationnaire
  # kcircRet * Ret0 = kcircRBC * RBC0
  p$kcircRet <- p$kcircRBC * p$RBC0 / p$Ret0

  # --- Transit Neut & Mono (non-prolifératif, 3 compartiments) ---
  # T1=T2=T3= (kcircNeut/aNeut)*Neut0  (Équation S3)
  T1N0 <- (p$kcircNeut / aNeut) * p$Neut0
  T1Mo0 <- (p$kcircMono / aMono) * p$Mono0

  # Flux entrants depuis CMP
  # ktrNeut * CMP0 = aNeut * T1N0  →  ktrNeut = aNeut * T1N0 / CMP0
  p$ktrNeut <- aNeut * T1N0 / p$CMP0
  p$ktrMono <- aMono * T1Mo0 / p$CMP0

  # Taux de prolifération CMP (Équation S4)
  # A l'état stationnaire : kprolCMP * CMP0 = ktrNeut*CMP0 + ktrMono*CMP0
  # Avec lambda3 : ktrNeut+ktrMono = lambda3 * kprolCMP
  p$kprolCMP <- (p$ktrNeut + p$ktrMono) / p$lambda3

  # --- Transit Ret (prolifératif : T1R, T2R ; non-prolifératif : T3R) ---
  # Équation S3 : T3R0 = T2R0 = (kcircRet/aRet)*Ret0
  T3R0 <- (p$kcircRet / aRet) * p$Ret0
  T2R0 <- T3R0
  T1R0 <- T2R0 / p$lambda1     # lambda1 = T2R0/T1R0

  # kprolRet — bilan SS de T2Ret (Équation S4, Fornari) :
  #   kprolRet*T2R0 + aRet*T1R0 = aRet*T2R0
  #   → kprolRet = aRet*(1 - 1/lambda1) = aRet*(lambda1-1)/lambda1
  p$kprolRet <- aRet * (p$lambda1 - 1) / p$lambda1

  # ktrRet — bilan SS de T1Ret :
  #   kprolRet*T1R0 + ktrRet*MEP0 = aRet*T1R0
  #   → ktrRet = (aRet - kprolRet)*T1R0 / MEP0
  #            = aRet * T2R0 / (lambda1^2 * MEP0)
  p$ktrRet <- (aRet - p$kprolRet) * T1R0 / p$MEP0

  # --- Transit Plt (prolifératif : T1P, T2P ; non-prolifératif : T3P) ---
  T3P0 <- (p$kcircPlt / aPlt) * p$Plt0
  T2P0 <- T3P0
  T1P0 <- T2P0 / p$lambda2

  # kprolPlt — bilan SS de T2Plt :
  #   kprolPlt = aPlt*(lambda2-1)/lambda2
  p$kprolPlt <- aPlt * (p$lambda2 - 1) / p$lambda2

  # ktrPlt — bilan SS de T1Plt :
  #   ktrPlt = (aPlt - kprolPlt)*T1P0 / MEP0
  p$ktrPlt <- (aPlt - p$kprolPlt) * T1P0 / p$MEP0

  # Taux de prolifération MEP (Équation S4)
  # lambda4 = (ktrRet + ktrPlt) / kprolMEP
  p$kprolMEP <- (p$ktrRet + p$ktrPlt) / p$lambda4

  # --- Transitions MPP → CMP et MPP → MEP ---
  # lambda5 = (ktrCMP + ktrMEP) / kprolMPP
  # Conservation à l'état stationnaire :
  # ktrCMP * MPP0 = kprolCMP * CMP0  (flux entrant = flux de prolifération CMP)
  p$ktrCMP <- p$kprolCMP * p$CMP0 / p$MPP0
  p$ktrMEP <- p$kprolMEP * p$MEP0 / p$MPP0

  # Taux de prolifération MPP
  p$kprolMPP <- (p$ktrCMP + p$ktrMEP) / p$lambda5

  # Flux source depuis HSC (kstem)
  # À l'état stationnaire : kstem + kprolMPP*MPP0 = ktrCMP*MPP0 + ktrMEP*MPP0
  # → kstem = (ktrCMP + ktrMEP - kprolMPP) * MPP0
  p$kstem <- (p$ktrCMP + p$ktrMEP - p$kprolMPP) * p$MPP0
  p$kstem <- max(p$kstem, 1e-8)

  # --- Stockage des états initiaux de transit (nécessaires pour CI) ---
  p$T1N0  <- T1N0
  p$T1Mo0 <- T1Mo0
  p$T1R0  <- T1R0
  p$T2R0  <- T2R0
  p$T3R0  <- T3R0
  p$T1P0  <- T1P0
  p$T2P0  <- T2P0
  p$T3P0  <- T3P0

  p
}

# -----------------------------------------------------------------------------
# SECTION 4 — VARIABILITÉ INTER-INDIVIDUELLE
#             Table S4 — Fornari 2019 (sigma log-additive)
#             + IIV De Carlo 2025 (omega log-normal pour PK)
# -----------------------------------------------------------------------------

#' Sigmas de l'erreur résiduelle log-additive (Table S4 — Fornari 2019)
#' Utilisés pour générer la variabilité entre patients dans la population
IIV_sigma <- list(
  sigma_MPP  = 0.33,
  sigma_CMP  = 0.19,
  sigma_MEP  = 0.33,
  sigma_Neut = 0.17,
  sigma_Mono = 0.17,
  sigma_Plt  = 0.17,
  sigma_Ret  = 0.36,
  sigma_RBC  = 0.06
)

#' CV% des paramètres physiologiques (Table 1 — Fornari 2019, rat)
#' Utilisés pour la variabilité des paramètres individuels (population)
IIV_params_CV <- list(
  cv_kcircNeut = 0.45,
  cv_kcircMono = 0.45,
  cv_kcircPlt  = 0.46,
  cv_kcircRBC  = 0.45,
  cv_Neut0     = 0.026,
  cv_Mono0     = 0.026,
  cv_Plt0      = 0.025,
  cv_Ret0      = 0.027,
  cv_RBC0      = 0.022,
  cv_MPP0      = 0.046,
  cv_CMP0      = 0.041,
  cv_MEP0      = 0.052,
  cv_MTTNeut   = 0.56,
  cv_MTTMono   = 0.55,
  cv_MTTPlt    = 0.056,
  cv_MTTRet    = 0.048,
  cv_SlopeMPP  = 0.118,
  cv_SlopeCMP  = 0.086,
  cv_SlopeMEP  = 0.103,
  cv_delta_Ret = 0.21,
  cv_delta_Plt = 0.304
)

#' IIV PK humain (omega log-normal — De Carlo 2025 / Zandvliet 2008)
IIV_pk_omega <- list(
  omega_CL  = 0.13,   # Zandvliet 2008
  omega_VCen = 0.54,
  omega_VPer = 0.31,
  omega_QC   = 0.46
)

# -----------------------------------------------------------------------------
# SECTION 5 — CALCUL DE LA CLAIRANCE HUMAINE (Calvert + Cockcroft-Gault)
# -----------------------------------------------------------------------------

#' Calcule la clairance du carboplatine (L/h) selon Cockcroft-Gault + Zandvliet
#'
#' @param age     Âge (ans)
#' @param weight  Poids (kg)
#' @param creat   Créatinine sérique (µmol/L)
#' @param sex     "M" ou "F"
#' @return CL (L/h)
calc_carboplatin_CL <- function(age = 60, weight = 70, creat = 80, sex = "M",
                                GFR = NULL) {
  # Si GFR fourni directement, on l'utilise (plus simple et cohérent avec Calvert)
  if (!is.null(GFR)) {
    # Calvert : CL_carbo (mL/min) = GFR + 25
    CL <- (GFR + 25) * 60 / 1000   # L/h
    return(CL)
  }

  # Cockcroft-Gault : SCr en µmol/L → CLcr en mL/min
  # CLcr (mL/min) = (140-age) × weight × sex_factor / (0.815 × SCr_µmol/L)
  # (0.815 = 72 / 88.4, conversion depuis la formule originale en mg/dL)
  sex_factor <- if (sex == "F") 0.85 else 1.0
  CLcg_mLmin <- (140 - age) * weight * sex_factor / (0.815 * creat)

  # Conversion mL/min → L/h, puis modèle Zandvliet 2008 :
  # CL_carbo (L/h) = CLcg (L/h) × theta1 + theta2
  # theta1 = 0.76, theta2 = 1.5 L/h (clairance non-rénale)
  CLcg_Lh <- CLcg_mLmin * 60 / 1000
  CL <- CLcg_Lh * 0.76 + 1.5
  CL
}

#' Calcule la dose de carboplatine via la formule de Calvert
#'
#' @param AUC_target AUC cible (mg.min/mL, typiquement 5)
#' @param GFR        Débit de filtration glomérulaire (mL/min, typiquement 125)
#' @return Dose (mg)
calc_calvert_dose <- function(AUC_target = 5, GFR = 125) {
  Dose_mg <- AUC_target * (GFR + 25)
  Dose_mg
}

#' Convertit une dose en mg vers µmol (PM carboplatine = 371.25 g/mol)
#'
#' @param dose_mg  Dose en mg
#' @return Dose en µmol
mg_to_umol <- function(dose_mg) {
  dose_mg / 371.25 * 1000
}

# -----------------------------------------------------------------------------
# SECTION 6 — FONCTION PRINCIPALE
# -----------------------------------------------------------------------------

#' Construit la liste complète de paramètres pour un individu
#'
#' @param species   "human" (défaut) ou "rat"
#' @param overrides Liste nommée de paramètres à remplacer (optionnel)
#' @param patient   Liste de covariables patient (age, weight, creat, sex)
#'                  utilisée pour calculer CL humain si non fourni
#'
#' @return Liste de paramètres complète (physiologie + SS dérivés + PK)
build_parameters <- function(species = "human",
                             overrides = list(),
                             patient = list(age = 60, weight = 70,
                                            creat = 80, sex = "M")) {

  # 1. Sélection de la base
  if (species == "rat") {
    p <- .rat_physiology
  } else {
    p <- .human_physiology
    # Calcul de CL humain via Cockcroft-Gault si non spécifié
    if (is.na(p$CL)) {
      p$CL <- calc_carboplatin_CL(
        age    = patient$age,
        weight = patient$weight,
        creat  = patient$creat,
        sex    = patient$sex,
        GFR    = patient$GFR   # NULL si non fourni → utilise CG
      )
    }
  }

  # 2. Application des overrides (IIV ou tests de sensibilité)
  for (nm in names(overrides)) {
    p[[nm]] <- overrides[[nm]]
  }

  # 3. Calcul des paramètres à l'état stationnaire
  p <- .compute_steady_state(p)

  # 4. Métadonnées
  p$species <- species

  p
}
