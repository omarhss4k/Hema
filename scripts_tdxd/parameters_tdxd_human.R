############################################################
# parameters_tdxd_human.R
# Paramètres PK/PD — T-DXd (Trastuzumab deruxtecan) — HUMAIN
#
# PK ADC : FDA BLA 761139 p.89/93 (= Yin et al. 2020 popPK 2-cpt)
#          Vc=2.77L, CL=0.42 L/j, T½_ADC=5.7j (NCA)
# PK DXd : Yin et al. 2020 (CL_DXd=19.2 L/h, V=29.4L → T½_int=1.06h)
#          NOTE: T½_DXd=5.8j FDA (p.93) est APPARENT (flip-flop NCA)
#          → taux limité par libération ADC, pas par CL intrinsèque de DXd
#          → CL_DXd=19.2 L/h est COMPATIBLE avec les données FDA
# k_int  : Vasalou et al. 2024 (internalisation HER2, t½=46h)
# PD     : Fornari 2019 (baselines humains, Table 1)
#
# Dose approuvée : 5.4 mg/kg Q3W IV (90 min)
# Driver E-R     : C_ADC1 surrogate (FDA p.95 : DXd Cavg, compromis)
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
  # NOTE: T½_int_DXd = ln(2)×29.41/19.2 = 1.06h (intrinsèque)
  # T½_DXd = 5.8j (FDA NCA p.93) est apparent = flip-flop (limité par libération ADC)

  # ── Mécanistiques ───────────────────────────────────────
  k_int       = log(2) / 46,    # 0.01507 h⁻¹  (Vasalou 2024, t½=46h HER2)

  # k_rel_c1 calibré depuis FDA BLA 761139 (Clinical Pharm Review) :
  #   Pseudo-SS (T½_DXd=1.06h → équilibre rapide) : C_DXd ≈ mass_frac×k_rel×C_ADC1×V1/CL_DXd
  #   Cmax_DXd = 4.4 ng/mL à C_ADC1_Cmax = 122 mg/L
  #   → k_rel = 0.0044 × 19.2 / (0.03885 × 122 × 2.77) = 0.00644 h⁻¹
  k_rel_c1    = 0.00644,        # h⁻¹  calibré FDA BLA clinical PK (Yin 2020)

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
  # Mais T½_DXd modèle = 1.06h → Cavg_DXd_ic ~ 0.015 µM << IC50=0.31 µM → E_drug~5%
  # → toxicité nulle avec DXd driver pur
  # Compromis : use_ADC_driver=TRUE (C_ADC1 surrogate, T½=23j > Q3W → E_drug persistant)
  use_ADC_driver     = TRUE,
  IC50_DXd_uM        = 0.31,   # µM — conservé (inactif si use_ADC_driver=TRUE)

  # ── Seuil de réparation ADN (D0) ────────────────────────
  # D_kill = max(0, Damage - D0) : Kill=0 entre cycles (trough ≈ 0 << D0)
  # → CMP récupère entre cycles → grade 0 possible
  Damage_threshold   = 0.05    # D0 seuil (rat: 0 par défaut dans ODE)
)

tdxd_pars_hu$mgL_to_uM_DXd <- 1000 / tdxd_pars_hu$MW_DXd

# ── Modèle de mélange TRIMODAL (DESTINY-Breast01) ───────────────────────────
# FDA BLA 761139 pooled (N=234) : G0=71% G1=7% G2=7% G3=13% G4=3%
# → 3 sous-groupes indépendants :
#   66% résistants  → G0          (Slope très faible)
#   14% modérés     → G1-G2       (Slope intermédiaire, à calibrer)
#   20% sensibles   → G3-4 ~80%   (Slope élevé)

p_sensitive_tdxd     <- 0.20   # G3-4 = 20% × 80% ≈ 16% (pooled FDA 16.2%)
p_moderate_tdxd      <- 0.14   # G1+G2 = 7%+7% = 14% (FDA pooled)
# p_resistant = 1 - p_sensitive - p_moderate = 0.66

Slope_resist_tdxd    <- 0.10   # kill_max ≈ 7%  → G0 garanti
Slope_moderate_tdxd  <- 12.0   # kill intermédiaire → G1-G2
Slope_sensitive_tdxd <- 59.1   # kill fort → G3-4 (nécessaire pour séparation binaire)

# Compatibilité backward
Slope_CMP_tdxd_human <- Slope_sensitive_tdxd

# ── Mixture TRIMODALE MEP — Anémie (DESTINY-Breast01) ───────────────────────
# FDA BLA 761139 (n=184) : G0=30% G1=37% G2=24% G3=8% G4=1%
# → 3 sous-groupes indépendants (tirages indép. de la mixture CMP) :
#   30% résistants MEP → G0          (Slope=0)
#   37% légers MEP     → G1          (Slope_light, à calibrer)
#   33% sensibles MEP  → G2-G3       (Slope_sensitive, inchangé)

p_MEP_light_tdxd         <- 0.33   # G1 cible (iter.8)
p_sensitive_mep_tdxd     <- 0.33   # G2-G3 cible (inchangé)
# p_MEP_resistant = 1 - 0.33 - 0.33 = 0.34 → G0 cible (iter.8)

Slope_MEP_resist_tdxd    <- 0.01   # résistants MEP : quasi-nul → G0 (<0.5% chute RBC)
Slope_MEP_light_tdxd     <- 0.18   # légers → centre G1 (ω élargi, Slope légèrement baissé)
Slope_MEP_sensitive_tdxd <- 1.25   # sensibles → G2-G3 (inchangé)
omega_Slope_MEP_sensitive <- 0.50  # IIV spread G2/G3
omega_Slope_MEP_light     <- 0.20  # IIV élargi (biologique), vs 0.10 artificiellement bas

# ── Cibles de validation clinique FDA BLA 761139 ─────────
# 5.4 mg/kg Q3W, popPK géométrique moyen steady-state
# Source : FDA BLA 761139 Multi-Discipline Review, p.89 (Table 11)
tdxd_hu_targets <- list(
  dose_mgkg          = 5.4,
  BW_kg              = 70,
  Tinfu_h            = 1.5,
  # PK ADC (fam-trastuzumab deruxtecan)
  Cmax_ADC_mgL       = 122,       # µg/mL  popPK SS geometric mean (CV=20%)
  AUCtau_ADC_ugdmL   = 735,       # µg·d/mL popPK SS geometric mean (CV=31%)
  Racc_ADC           = 1.3,       # accumulation ratio SS/C1 (Q3W)
  # PK DXd (MAAA-1181a, free plasma)
  Cmax_DXd_ngmL      = 4.4,       # ng/mL  popPK SS geometric mean (CV=40%)
  Cmax_DXd_mgL       = 4.4e-3,
  AUCtau_DXd_ngdmL   = 28,        # ng·d/mL popPK SS geometric mean (CV=37%)
  # Tolerance
  tol_Cmax           = 0.30       # ±30%
)

# ── Effets cliniques attendus (FDA BLA 761139, Safety Review) ─
# 5.4 mg/kg Q3W
# Pooled 5.4 mg/kg (N=234, Table 44 FDA p.189-190) :
#   Neutropénie G3-4 : 16.2%  (lab any grade : 61.5%)
#   Anémie      G3-4 :  7.3%  (lab any grade : 69.7%)
#   Thrombocy.  G3-4 :  3.4%  (lab any grade : 36.8%)
# DESTINY-Breast01 (U201, N=184) — cibles primaires modèle :
#   Neutropénie G3-4 : ~20%   Tout grade : ~29%
#   Anémie      G3-4 :  ~9%   Tout grade : ~70%
#   Thrombocytopénie G3-4 : ~3-4%
# → E-R : DXd Cavg (cycle-averaged) = prédicteur hématologie (FDA p.94)

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
cat(sprintf("  CL_DXd = %.1f L/h  V_DXd = %.2f L  T½_int = %.2f h\n",
            tdxd_pars_hu$CL_DXd, tdxd_pars_hu$V_DXd,
            log(2) * tdxd_pars_hu$V_DXd / tdxd_pars_hu$CL_DXd))
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
cat(sprintf("  Driver toxicité : C_ADC1 [µg/mL]  (assay PFB-10, surrogate)\n"))
cat(sprintf("  IC50_ADC : ery=%.1f µg/mL  neut=%.1f µg/mL  shared=%.2f µg/mL\n",
            tdxd_pars_hu$IC50_ADC_ery_ugmL,
            tdxd_pars_hu$IC50_ADC_neut_ugmL,
            tdxd_pars_hu$IC50_ADC_ugmL))
cat(sprintf("  C_Avg_clinique = 33.3 µg/mL → E_drug_avg = %.3f\n\n",
            33.3 / (tdxd_pars_hu$IC50_ADC_ugmL + 33.3)))
