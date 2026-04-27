############################################################
# parameters_tdxd_nhp_pd.R
# Paramètres PD — T-DXd — SINGE CYNOMOLGUS (NHP)
#
# Méthode :
#   1. Baselines hématologiques : valeurs de référence cynomolgus
#      (Chamanza et al. 2010, SNBL reference ranges)
#   2. MTT (Mean Transit Times) : allométrie depuis rat BW^0.25
#      facteur = (BW_nhp/BW_rat)^0.25 = (4/0.25)^0.25 = 2.0
#   3. k_circ : spécifiques NHP (durée de vie cellules circulantes)
#   4. Paramètres dérivés : Équation S4 exacte (Fornari 2019)
#   5. Slopes T-DXd : point de départ = rat calibré FDA
#      (driver DXd intracellulaire, nécessite calibration si données PD NHP)
#
# Sources :
#   FDA BLA 761139 Multi-Discipline Review (2019), Table 7
#   Fornari et al. 2019, CPT:PSP (modèle hématopoïèse)
#   Chamanza et al. 2010 (valeurs référence cynomolgus)
#   SNBL cynomolgus monkey historical hematology reference ranges
#
# Unités : temps [h] | concentrations circulantes [×10⁹/L]
#          progéniteurs [×10⁶/kg] | taux [h⁻¹]
############################################################

nhp_pd <- list()

# ════════════════════════════════════════════════════════
# 1. PARAMÈTRES ALLOMÉTRIQUES
# ════════════════════════════════════════════════════════
BW_rat_pd <- 0.25   # kg
BW_nhp_pd <- 4.0    # kg

# Facteur MTT : temps biologiques ∝ BW^0.25
allo_MTT <- (BW_nhp_pd / BW_rat_pd)^0.25   # = 16^0.25 = 2.0

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PD T-DXd — SINGE CYNOMOLGUS (Fornari adapté)              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("  Allométrie MTT : (%.1f/%.2f)^0.25 = %.2f\n\n",
            BW_nhp_pd, BW_rat_pd, allo_MTT))

# ════════════════════════════════════════════════════════
# 2. PARAMÈTRES λ (identiques rat — Fornari Table 1)
# ════════════════════════════════════════════════════════
nhp_pd$lambda1 <- 2.0   # amplification transit Ret
nhp_pd$lambda2 <- 2.0   # amplification transit Plt
nhp_pd$lambda3 <- 1.8   # prolifération CMP
nhp_pd$lambda4 <- 1.8   # prolifération MEP
nhp_pd$lambda5 <- 1.8   # prolifération MPP

# ════════════════════════════════════════════════════════
# 3. BASELINES HÉMATOLOGIQUES NHP CYNOMOLGUS
#    Référence : Chamanza 2010, SNBL ranges
#    Unités circulants : ×10⁹/L | progéniteurs : ×10⁶/kg
# ════════════════════════════════════════════════════════

# Progéniteurs moelle osseuse — densité/kg similaire au rat
# (cellularité moelle comparable entre espèces)
nhp_pd$MPP0 <- 69.0    # ×10⁶/kg  (rat Fornari Table 1)
nhp_pd$CMP0 <- 214.0   # ×10⁶/kg
nhp_pd$MEP0 <- 156.0   # ×10⁶/kg

# Cellules circulantes — valeurs référence cynomolgus
nhp_pd$Neut0 <- 2.50    # ×10⁹/L  (plage 1.0–8.0, moy ~3.5)
nhp_pd$Mono0 <- 0.35    # ×10⁹/L  (plage 0.1–0.7)
nhp_pd$Ret0  <- 30.0    # ×10⁹/L  (plage 15–50)
nhp_pd$RBC0  <- 5200.0  # ×10⁹/L  (= 5.2 ×10¹²/L, plage 4.5–6.5)
nhp_pd$Plt0  <- 400.0   # ×10⁹/L  (plage 200–600)

cat("── Baselines NHP (réf. cynomolgus) ──\n")
cat(sprintf("  MPP=%.0f  CMP=%.0f  MEP=%.0f  [×10⁶/kg]\n",
            nhp_pd$MPP0, nhp_pd$CMP0, nhp_pd$MEP0))
cat(sprintf("  Neut=%.2f  Mono=%.2f  Ret=%.0f  RBC=%.0f  Plt=%.0f  [×10⁹/L]\n\n",
            nhp_pd$Neut0, nhp_pd$Mono0, nhp_pd$Ret0, nhp_pd$RBC0, nhp_pd$Plt0))

# ════════════════════════════════════════════════════════
# 4. MTT ALLOMÉTRIQUES (rat × 2.0)
#    Valeurs rat (Fornari Table 1) : Neut=61.5h, Mono=62.2h,
#    Ret=75h, Plt=74h
# ════════════════════════════════════════════════════════
nhp_pd$MTT_Neut <- 61.5 * allo_MTT    # = 123.0 h = 5.1 j
nhp_pd$MTT_Mono <- 62.2 * allo_MTT    # = 124.4 h = 5.2 j
nhp_pd$MTT_Ret  <- 75.0 * allo_MTT    # = 150.0 h = 6.25 j
nhp_pd$MTT_Plt  <- 74.0 * allo_MTT    # = 148.0 h = 6.17 j

cat("── MTT allométriques NHP (rat × 2.0) ──\n")
cat(sprintf("  MTT_Neut=%.1fh (%.1fj)  MTT_Mono=%.1fh (%.1fj)\n",
            nhp_pd$MTT_Neut, nhp_pd$MTT_Neut/24,
            nhp_pd$MTT_Mono, nhp_pd$MTT_Mono/24))
cat(sprintf("  MTT_Ret =%.1fh (%.1fj)  MTT_Plt =%.1fh (%.1fj)\n\n",
            nhp_pd$MTT_Ret, nhp_pd$MTT_Ret/24,
            nhp_pd$MTT_Plt, nhp_pd$MTT_Plt/24))

# ════════════════════════════════════════════════════════
# 5. k_CIRC — durées de vie NHP spécifiques
#    Neutrophiles  : ~6h sang → k = 0.17 /h  (identique rat)
#    Monocytes     : ~17h sang → k = 0.06 /h (identique rat)
#    GR (RBC)      : ~100 j cynomolgus vs 60j rat
#                    k_circ_RBC = 1/(100×24) = 4.17e-4 /h
#    Plaquettes    : ~8.5 j cynomolgus vs 4.2j rat
#                    k_circ_Plt = 1/(8.5×24) = 4.90e-3 /h
# ════════════════════════════════════════════════════════
nhp_pd$k_circ_Neut <- 0.17               # /h  (même rat)
nhp_pd$k_circ_Mono <- 0.06               # /h  (même rat)
nhp_pd$k_circ_RBC  <- 1 / (100 * 24)    # /h  ≈ 4.17e-4 (T½ 100j)
nhp_pd$k_circ_Plt  <- 1 / (8.5 * 24)    # /h  ≈ 4.90e-3 (T½ 8.5j)

cat("── k_circ NHP (durées de vie spécifiques) ──\n")
cat(sprintf("  k_circ_Neut=%.4f /h  (T½=%.1fh=%.1fj)\n",
            nhp_pd$k_circ_Neut, 1/nhp_pd$k_circ_Neut, 1/nhp_pd$k_circ_Neut/24))
cat(sprintf("  k_circ_Mono=%.4f /h  (T½=%.1fh=%.1fj)\n",
            nhp_pd$k_circ_Mono, 1/nhp_pd$k_circ_Mono, 1/nhp_pd$k_circ_Mono/24))
cat(sprintf("  k_circ_RBC =%.2e /h  (T½=%.0fj)\n",
            nhp_pd$k_circ_RBC, 1/nhp_pd$k_circ_RBC/24))
cat(sprintf("  k_circ_Plt =%.2e /h  (T½=%.1fj)\n\n",
            nhp_pd$k_circ_Plt, 1/nhp_pd$k_circ_Plt/24))

# ════════════════════════════════════════════════════════
# 6. PARAMÈTRES DÉRIVÉS — Équation S4 exacte (Fornari 2019)
# ════════════════════════════════════════════════════════

# k_circ_Ret dérivé (Eq S4)
nhp_pd$k_circ_Ret <- nhp_pd$k_circ_RBC * nhp_pd$RBC0 / nhp_pd$Ret0

# Constantes a (maturation transit)
nhp_pd$a_Ret <- 3 / nhp_pd$MTT_Ret
nhp_pd$a_Plt <- 3 / nhp_pd$MTT_Plt
a_Neut_nhp   <- 3 / nhp_pd$MTT_Neut
a_Mono_nhp   <- 3 / nhp_pd$MTT_Mono

# k_tr_Neut, k_tr_Mono (Eq S4)
nhp_pd$k_tr_Neut <- nhp_pd$k_circ_Neut * nhp_pd$Neut0 / nhp_pd$CMP0
nhp_pd$k_tr_Mono <- nhp_pd$k_circ_Mono * nhp_pd$Mono0 / nhp_pd$CMP0

# k_tr_Ret, k_tr_Plt (Eq S4, ÷λ²)
nhp_pd$k_tr_Ret <- nhp_pd$k_circ_Ret * nhp_pd$Ret0 /
                   (nhp_pd$lambda1^2 * nhp_pd$MEP0)
nhp_pd$k_tr_Plt <- nhp_pd$k_circ_Plt * nhp_pd$Plt0 /
                   (nhp_pd$lambda2^2 * nhp_pd$MEP0)

# k_prol_Ret, k_prol_Plt (Eq S4)
nhp_pd$k_prol_Ret <- nhp_pd$a_Ret * (1 - 1 / nhp_pd$lambda1)
nhp_pd$k_prol_Plt <- nhp_pd$a_Plt * (1 - 1 / nhp_pd$lambda2)

# k_prol_CMP, k_prol_MEP (Eq S4)
nhp_pd$k_prol_CMP <- (nhp_pd$k_tr_Neut + nhp_pd$k_tr_Mono) / nhp_pd$lambda3
nhp_pd$k_prol_MEP <- (nhp_pd$k_tr_Ret  + nhp_pd$k_tr_Plt)  / nhp_pd$lambda4

# k_tr_CMP, k_tr_MEP (Eq S4)
nhp_pd$k_tr_CMP <- (nhp_pd$k_tr_Neut + nhp_pd$k_tr_Mono - nhp_pd$k_prol_CMP) *
                   nhp_pd$CMP0 / nhp_pd$MPP0
nhp_pd$k_tr_MEP <- (nhp_pd$k_tr_Ret  + nhp_pd$k_tr_Plt  - nhp_pd$k_prol_MEP) *
                   nhp_pd$MEP0 / nhp_pd$MPP0

# k_prol_MPP (Eq S4)
nhp_pd$k_prol_MPP <- (nhp_pd$k_tr_CMP + nhp_pd$k_tr_MEP) / nhp_pd$lambda5

# k_stem (Eq S4)
nhp_pd$k_stem <- (nhp_pd$k_tr_CMP + nhp_pd$k_tr_MEP - nhp_pd$k_prol_MPP) *
                 nhp_pd$MPP0

cat("── Paramètres dérivés Eq S4 ──\n")
cat(sprintf("  k_circ_Ret = %.5f /h  (Eq S4 : k_RBC×RBC0/Ret0)\n", nhp_pd$k_circ_Ret))
cat(sprintf("  k_tr_Neut  = %.6f /h\n", nhp_pd$k_tr_Neut))
cat(sprintf("  k_tr_Mono  = %.6f /h\n", nhp_pd$k_tr_Mono))
cat(sprintf("  k_tr_Ret   = %.6f /h\n", nhp_pd$k_tr_Ret))
cat(sprintf("  k_tr_Plt   = %.6f /h\n", nhp_pd$k_tr_Plt))
cat(sprintf("  k_prol_CMP = %.6f /h\n", nhp_pd$k_prol_CMP))
cat(sprintf("  k_prol_MEP = %.6f /h\n", nhp_pd$k_prol_MEP))
cat(sprintf("  k_tr_CMP   = %.6f /h\n", nhp_pd$k_tr_CMP))
cat(sprintf("  k_tr_MEP   = %.6f /h\n", nhp_pd$k_tr_MEP))
cat(sprintf("  k_prol_MPP = %.6f /h\n", nhp_pd$k_prol_MPP))
cat(sprintf("  k_stem     = %.6f\n\n",  nhp_pd$k_stem))

# ════════════════════════════════════════════════════════
# 7. PARAMÈTRES FEEDBACK (Fornari Table 1 — identiques rat)
# ════════════════════════════════════════════════════════
nhp_pd$gamma_stem      <- 0.07
nhp_pd$gamma_mat_CMP   <- 0.60
nhp_pd$gamma_mat_MEP   <- 0.30
nhp_pd$gamma_prolTrans <- 0.70

# ════════════════════════════════════════════════════════
# 8. PARAMÈTRES EFFET MÉDICAMENT — T-DXd NHP
#    Driver : DXd intracellulaire (identique rat, Eq damage)
#    Slope_MEP = 1.00 (calibré FDA rat, point de départ NHP)
#    NOTE : Les concentrations DXd_ic NHP sont ~10× inférieures
#    au rat (Krel_nhp/Krel_rat ≈ 0.09). Les effets PD attendus
#    sont minimes, cohérent avec utilisation du singe comme espèce
#    de sécurité pour T-DXd.
# ════════════════════════════════════════════════════════
nhp_pd$k_dam       <- 0.017   # /h (identique rat/humain)
nhp_pd$k_rep       <- 0.017   # /h
nhp_pd$IC50_DXd_uM <- 0.31    # µM (topo I, Yin 2020)

# Slopes (point de départ = rat calibré FDA)
nhp_pd$Slope_MPP <- 0.15
nhp_pd$Slope_CMP <- 1.47
nhp_pd$Slope_MEP <- 1.00   # calibré FDA rat Table 6

# Multiplicateurs effet sur lignées prolifératives transit
nhp_pd$delta_Ret <- 2.8
nhp_pd$delta_Plt <- 0.54

cat("── Effet médicament T-DXd NHP ──\n")
cat(sprintf("  Driver    : DXd intracellulaire (C_DXd_ic [µM])\n"))
cat(sprintf("  IC50_DXd  : %.2f µM  |  k_dam/k_rep : %.3f /h\n",
            nhp_pd$IC50_DXd_uM, nhp_pd$k_dam))
cat(sprintf("  Slope_MPP : %.2f  |  Slope_CMP : %.2f  |  Slope_MEP : %.2f\n",
            nhp_pd$Slope_MPP, nhp_pd$Slope_CMP, nhp_pd$Slope_MEP))
cat(sprintf("  delta_Ret : %.2f  |  delta_Plt : %.2f\n\n",
            nhp_pd$delta_Ret, nhp_pd$delta_Plt))
cat("  NOTE: DXd_ic NHP << IC50 (Krel_NHP ≈ 0.09×Krel_rat)\n")
cat("  → Effets PD faibles attendus = cohérent avec la biologie\n\n")

# ════════════════════════════════════════════════════════
# 9. ÉTATS INITIAUX TRANSIT (Équation S3)
# ════════════════════════════════════════════════════════

# Transit Neut et Mono (compartiments non-prolifératifs)
T_Neut0 <- nhp_pd$k_circ_Neut * nhp_pd$Neut0 / a_Neut_nhp
T_Mono0 <- nhp_pd$k_circ_Mono * nhp_pd$Mono0 / a_Mono_nhp

# Transit Ret et Plt (compartiments prolifératifs, Eq S3)
T2_Ret0 <- nhp_pd$k_circ_Ret * nhp_pd$Ret0 / nhp_pd$a_Ret
T1_Ret0 <- T2_Ret0 / nhp_pd$lambda1
T2_Plt0 <- nhp_pd$k_circ_Plt * nhp_pd$Plt0 / nhp_pd$a_Plt
T1_Plt0 <- T2_Plt0 / nhp_pd$lambda2

cat("── États initiaux transit (Eq S3) ──\n")
cat(sprintf("  T_Neut = %.4f  |  T_Mono = %.4f\n",   T_Neut0, T_Mono0))
cat(sprintf("  T1_Ret = %.4f  |  T2=T3_Ret = %.4f\n", T1_Ret0, T2_Ret0))
cat(sprintf("  T1_Plt = %.4f  |  T2=T3_Plt = %.4f\n\n", T1_Plt0, T2_Plt0))

# ════════════════════════════════════════════════════════
# 10. ÉTAT INITIAL COMPLET (PK + PD)
# ════════════════════════════════════════════════════════
nhp_pkpd_state0 <- c(
  # PK ADC (2-cpt) + DXd + Damage
  C_ADC1   = 0,
  C_ADC2   = 0,
  C_DXd    = 0,
  C_DXd_ic = 0,
  Damage   = 0,
  # PD — Progéniteurs
  MPP = nhp_pd$MPP0,
  CMP = nhp_pd$CMP0,
  MEP = nhp_pd$MEP0,
  # PD — Neutrophiles (transit + circulants)
  T1_Neut = T_Neut0, T2_Neut = T_Neut0, T3_Neut = T_Neut0,
  Neut    = nhp_pd$Neut0,
  # PD — Monocytes
  T1_Mono = T_Mono0, T2_Mono = T_Mono0, T3_Mono = T_Mono0,
  Mono    = nhp_pd$Mono0,
  # PD — Réticulocytes + GR
  T1_Ret = T1_Ret0, T2_Ret = T2_Ret0, T3_Ret = T2_Ret0,
  Ret    = nhp_pd$Ret0,
  RBC    = nhp_pd$RBC0,
  # PD — Plaquettes
  T1_Plt = T1_Plt0, T2_Plt = T2_Plt0, T3_Plt = T2_Plt0,
  Plt    = nhp_pd$Plt0
)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PARAMÈTRES PD NHP INITIALISÉS ✓  (25 états PK+PD)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
