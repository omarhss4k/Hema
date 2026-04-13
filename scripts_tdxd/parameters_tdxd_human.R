############################################################
# parameters_tdxd_human.R
# Paramètres PK/PD — T-DXd (Trastuzumab deruxtecan) — HUMAIN
#
# PK  : Yin et al. 2020 (PopPK 2-compartiments)
#       k_int : Vasalou et al. 2024 (internalisation HER2)
# PD  : Fornari 2019 (baselines humains, Table 1)
#       Slopes carboplatin humain (Table 2) — point de départ
#
# Dose approuvée : 5.4 mg/kg Q3W IV (90 min)
# Validation     : FDA BLA 761139 Clinical Pharmacology Review
############################################################

# ── Baselines PD humains (Fornari 2019 Table 1) ──────────
source("../scripts/parameters_human.R")
source("../scripts/parameters_FORNARI_CORRECT.R")

# ── Paramètres PK T-DXd humain ───────────────────────────
# Source : Yin et al. 2020, PopPK estimés sur patients HER2+
# k_int  : Vasalou 2024 (t½ internalisation = 46h)
# k_rel  : Yin 2020, cycle-dépendant

tdxd_pars_hu <- list(

  # ── PK ADC — 2 compartiments ───────────────────────────
  CL_ADC      = 0.421  / 24,    # 0.01754 L/h  (Yin 2020)
  V1_ADC      = 2.77,            # L            (Yin 2020)
  V2_ADC      = 5.16,            # L            (Yin 2020)
  Q_ADC       = 0.174  / 24,    # 0.00725 L/h  (Yin 2020)

  # ── PK DXd payload — 1 compartiment ────────────────────
  V_DXd       = 17 * 1.73,      # 29.41 L      (Yin 2020)
  CL_DXd      = 19.2,           # L/h          (Yin 2020)

  # ── Mécanistiques ───────────────────────────────────────
  k_int       = log(2) / 46,    # 0.01507 h⁻¹  (Vasalou 2024, t½=46h HER2)

  # k_rel_c1 calibré depuis FDA BLA 761139 (Clinical Pharm Review) :
  #   Cmax_DXd = 4.4 ng/mL = 0.0044 mg/L à C_ADC1_Cmax = 122 mg/L
  #   Pseudo-SS : C_DXd = mass_frac × k_rel × C_ADC1 × V1 / CL_DXd
  #   → k_rel = 0.0044 × 19.2 / (0.03885 × 122 × 2.77) = 0.00644 h⁻¹
  #   (Yin 2020 = 0.0159 h⁻¹ : formulation ODE différente → non comparable directement)
  k_rel_c1    = 0.00644,        # h⁻¹  calibré FDA BLA clinical PK

  krel_power  = -0.137,         # exposant cycle (Yin 2020)
  krel_factor = 0.830,          # réduction cycle > 1
  interval_h  = 21 * 24,        # 504 h = Q3W

  # ── Échanges DXd membranaires ───────────────────────────
  k_inD       = 0.7,            # h⁻¹  (entrée intracell.)
  k_effD      = 0.7,            # h⁻¹  (efflux)

  # ── Compartiment intracellulaire moelle ─────────────────
  # Moelle ~ 1.5% BW × 70kg × 70% intracell. ≈ 0.74 L
  V_ic        = 0.003 * (70 / 0.25),   # 0.84 L  (rat × 280)

  # ── Propriétés ADC ──────────────────────────────────────
  DAR           = 8,
  MW_ADC        = 148000,
  MW_DXd        = 718.8,
  mass_frac_DXd = 8 * 718.8 / 148000,  # 0.03885

  # ── Dommages ADN ────────────────────────────────────────
  k_dam         = 0.017,        # h⁻¹  (Fornari)
  k_rep         = 0.017,        # h⁻¹  (Fornari)

  # ── IC50 ADC (T-DXd entier) — assay PFB-10 sur HPC ──────
  # Source : PFB-10 Colony Forming Unit assay, progéniteurs humains
  # Driver : C_ADC1 [µg/mL = mg/L]  (pas le DXd libre)
  # C_Avg clinique = 700 µg/mL·jour / 21 jours = 33.3 µg/mL
  IC50_ADC_ery_ugmL  = 27.3,   # µg/mL — Total Erythroid (MEP/Ret)
  IC50_ADC_neut_ugmL = 28.1,   # µg/mL — CFU-GM (CMP/Neut)
  # Moyenne géométrique pour Damage partagé :
  IC50_ADC_ugmL      = sqrt(27.3 * 28.1),   # 27.70 µg/mL
  # Flag : utiliser C_ADC1 comme driver (pas C_DXd_ic)
  use_ADC_driver     = TRUE,
  IC50_DXd_uM        = 0.31    # conservé mais inactif si use_ADC_driver=TRUE
)

tdxd_pars_hu$mgL_to_uM_DXd <- 1000 / tdxd_pars_hu$MW_DXd

# ── Cibles de validation clinique FDA BLA 761139 ─────────
# 5.4 mg/kg Q3W, géométrique moyen cycle 1
tdxd_hu_targets <- list(
  dose_mgkg       = 5.4,
  BW_kg           = 70,
  Tinfu_h         = 1.5,
  Cmax_ADC_mgL    = 122,        # µg/mL  (FDA Clinical Pharm Review)
  Cmax_DXd_ngmL   = 4.4,        # ng/mL  (DXd libre plasmatique)
  Cmax_DXd_mgL    = 4.4e-3,
  tol_Cmax        = 0.30        # ±30%
)

# ── Effets cliniques attendus (FDA BLA 761139, Safety Review) ─
# 5.4 mg/kg Q3W (DESTINY-Breast01, n=184)
# Grade 3-4 neutropenia : ~20% patients → Neut < 1.0 × 10⁹/L
# Grade 1-2 neutropenia : ~35% patients → 1.0 < Neut < 2.0
# Grade 3-4 anemia      : ~10%
# → Le modèle doit prédire une chute modérée de Neut (~20-40%)
#   cohérente avec Grade 1-2 en médiane

# ── État initial ─────────────────────────────────────────
tdxd_hu_state0 <- c(
  C_ADC1    = 0,
  C_ADC2    = 0,
  C_DXd     = 0,
  C_DXd_ic  = 0,
  Damage    = 0
)

# ── Résumé ───────────────────────────────────────────────
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  PK/PD T-DXd — HUMAIN (Yin 2020 + Vasalou 2024)       ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat("── ADC (2-compartiments) ──\n")
cat(sprintf("  CL_ADC = %.5f L/h  V1 = %.2f L  V2 = %.2f L\n",
            tdxd_pars_hu$CL_ADC, tdxd_pars_hu$V1_ADC, tdxd_pars_hu$V2_ADC))
cat("── DXd payload ──\n")
cat(sprintf("  CL_DXd = %.1f L/h  V_DXd = %.2f L\n",
            tdxd_pars_hu$CL_DXd, tdxd_pars_hu$V_DXd))
cat("── Mécanistiques ──\n")
cat(sprintf("  k_int = %.5f h⁻¹  Krel_C1 = %.4f h⁻¹\n",
            tdxd_pars_hu$k_int, tdxd_pars_hu$k_rel_c1))
cat(sprintf("  V_ic  = %.3f L  (ratio V_DXd/V_ic = %.0f)\n",
            tdxd_pars_hu$V_ic,
            tdxd_pars_hu$V_DXd / tdxd_pars_hu$V_ic))
cat("── Baselines PD (Fornari 2019, humain) ──\n")
cat(sprintf("  Neut0=%.1f  Ret0=%.0f  MEP0=%.1f  Plt0=%.0f  [×10⁹/L]\n",
            init_pars$Neut0, init_pars$Ret0,
            init_pars$MEP0,  init_pars$Plt0))
cat(sprintf("  Slope_MEP=%.2f  Slope_CMP=%.2f  Slope_MPP=%.2f\n",
            init_pars$Slope_MEP, init_pars$Slope_CMP, init_pars$Slope_MPP))
cat(sprintf("  Driver toxicité : C_ADC1 [µg/mL]  (assay PFB-10)\n"))
cat(sprintf("  IC50_ADC : ery=%.1f µg/mL  neut=%.1f µg/mL  shared=%.2f µg/mL\n",
            tdxd_pars_hu$IC50_ADC_ery_ugmL,
            tdxd_pars_hu$IC50_ADC_neut_ugmL,
            tdxd_pars_hu$IC50_ADC_ugmL))
cat(sprintf("  C_Avg_clinique = 33.3 µg/mL → E_drug_avg = %.3f\n\n",
            33.3 / (tdxd_pars_hu$IC50_ADC_ugmL + 33.3)))
