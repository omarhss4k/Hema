############################################################
# parameters_tdxd_human.R
# Paramètres PK/PD — T-DXd (Trastuzumab deruxtecan) — HUMAIN
#
# PK ADC : FDA BLA 761139 p.89/93 (= Yin et al. 2020 popPK 2-cpt)
#          Vc=2.77L, CL=0.42 L/j, T½_ADC=5.7j (NCA)
# PK DXd : FDA BLA 761139 p.93 — T½_DXd=5.8j (NCA apparent)
#          CL_DXd dérivé de T½ FDA (remplace Yin 2020 CL=19.2 L/h)
# k_int  : Vasalou et al. 2024 (internalisation HER2, t½=46h)
# PD     : Fornari 2019 (baselines humains, Table 1)
#
# Dose approuvée : 5.4 mg/kg Q3W IV (90 min)
# Driver E-R     : DXd Cavg (FDA BLA 761139 p.95)
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
  # CL_DXd recalculé depuis FDA BLA 761139 (p.93) NCA T½_DXd = 5.8j (apparent)
  # T½_app = ln(2) × V_DXd / CL_DXd  →  CL_DXd = ln(2) × 29.41 / (5.8×24) = 0.1464 L/h
  # (Yin 2020 CL_DXd = 19.2 L/h donne T½ = 1.06h — incompatible avec T½ FDA apparent)
  CL_DXd      = log(2) * (17 * 1.73) / (5.8 * 24), # 0.1464 L/h  (FDA BLA 761139 T½=5.8j)

  # ── Mécanistiques ───────────────────────────────────────
  k_int       = log(2) / 46,    # 0.01507 h⁻¹  (Vasalou 2024, t½=46h HER2)

  # k_rel_c1 recalibré avec nouveau CL_DXd (FDA T½=5.8j) :
  #   Pseudo-SS : C_DXd = mass_frac × k_rel × C_ADC1 × V1 / CL_DXd
  #   Cible : Cmax_DXd = 4.4 ng/mL = 0.0044 mg/L (FDA BLA p.89)
  #   → k_rel = 0.0044 × CL_DXd / (0.03885 × 122 × 2.77)
  #   → k_rel = 0.0044 × 0.1464 / (0.03885 × 122 × 2.77) = 4.91e-05 h⁻¹
  #   Ratio k_rel/CL_DXd conservé → Cmax_DXd inchangé (= 4.4 ng/mL)
  k_rel_c1    = 0.0044 * (log(2) * (17*1.73) / (5.8*24)) /
                (8 * 718.8/148000 * 122 * 2.77),  # ≈ 4.91e-05 h⁻¹  (FDA BLA 761139)

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
  # Conservé pour référence mais inactif (use_ADC_driver = FALSE)
  IC50_ADC_ery_ugmL  = 27.3,
  IC50_ADC_neut_ugmL = 28.1,
  IC50_ADC_ugmL      = sqrt(27.3 * 28.1),   # 27.70 µg/mL

  # ── Driver toxicité ─────────────────────────────────────
  # FDA BLA 761139 (p.95) : DXd Cavg est le prédicteur E-R (statistique)
  # Avec CL_DXd corrigé (T½=5.8j FDA) : C_DXd_ic Cavg ~ 0.058 µM
  #   → E_drug ~ 16%  (vs ~3% avec T½=1.06h)
  # use_ADC_driver=FALSE désormais cohérent avec FDA (DXd driver actif)
  use_ADC_driver     = FALSE,
  IC50_DXd_uM        = 0.31    # µM — DXd libre (Yin 2020 / FDA BLA 761139)
)

tdxd_pars_hu$mgL_to_uM_DXd <- 1000 / tdxd_pars_hu$MW_DXd

# ── Slope_CMP calibré pour T-DXd humain ──────────────────
# Calibration depuis DESTINY-Breast01 (FDA BLA 761139, n=184) :
#   Driver : C_DXd_ic [µM]  (FDA p.95 — DXd Cavg prédicteur E-R)
#   CL_DXd corrigé (T½=5.8j) → C_DXd_ic Cavg ~ 0.058 µM → E_drug ~ 16%
#   Cible : G3-4 neutropenie ~16-20%
#   Slope_CMP = 12.0 → RECALIBRATION NECESSAIRE après changement de driver
#   (ancien driver C_ADC1 avec E_drug_avg ~ 54% → Slope × 0.54 = 6.5 equiv.)
Slope_CMP_tdxd_human <- 12.0   # A recalibrer — run calibrate_slope_cmp_human.R

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
cat("── DXd payload (FDA BLA 761139 T½=5.8j) ──\n")
cat(sprintf("  CL_DXd = %.4f L/h  V_DXd = %.2f L  T½ = %.1f h (%.1f j)\n",
            tdxd_pars_hu$CL_DXd, tdxd_pars_hu$V_DXd,
            log(2) * tdxd_pars_hu$V_DXd / tdxd_pars_hu$CL_DXd,
            log(2) * tdxd_pars_hu$V_DXd / tdxd_pars_hu$CL_DXd / 24))
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
cat(sprintf("  Driver toxicité : C_DXd_ic [µM]  (FDA BLA 761139 p.95)\n"))
cat(sprintf("  IC50_DXd = %.2f µM\n", tdxd_pars_hu$IC50_DXd_uM))
Cavg_DXd_ic_approx <- 35 * 1.19e-3 * 1000 / tdxd_pars_hu$MW_DXd
cat(sprintf("  C_DXd_ic Cavg approx = %.4f µM → E_drug_avg ~ %.1f%%\n",
            Cavg_DXd_ic_approx,
            100 * Cavg_DXd_ic_approx / (tdxd_pars_hu$IC50_DXd_uM + Cavg_DXd_ic_approx)))
