############################################################
# parameters_tdxd_rat.R
# Paramètres PK — T-DXd (Trastuzumab deruxtecan) — RAT
#
# Sources :
#   ADC serum : Yin et al. (2020), 2-compartiments humain
#               → transposition allométrique rat (250g)
#   DXd payload : Yin et al. (2020), scalé allométriquement
#   Mécanistiques : Vasalou et al. (2024)
#
# Allométrie : CL ∝ BW^0.75 | V ∝ BW^1.0
# Temps en HEURES | Concentrations en mg/L (= µg/mL)
############################################################

tdxd_pars <- list()

# ── Poids corporels ──────────────────────────────────────
tdxd_pars$BW_human <- 70     # kg
tdxd_pars$BW_rat   <- 0.25   # kg

# Facteurs allométriques (rat / humain)
scale_CL <- (tdxd_pars$BW_rat / tdxd_pars$BW_human)^0.75  # ~0.01460
scale_V  <- (tdxd_pars$BW_rat / tdxd_pars$BW_human)^1.0   # ~0.003571

# ── Poids moléculaires ───────────────────────────────────
tdxd_pars$MW_ADC  <- 148000  # g/mol  (T-DXd ~148 kDa)
tdxd_pars$MW_DXd  <- 718.8   # g/mol  (exatecan derivative)

# ── Propriétés ADC ───────────────────────────────────────
tdxd_pars$DAR <- 8   # Drug-to-Antibody Ratio (molaire)

# Fraction massique payload par molécule ADC
# = DAR × (MW_DXd / MW_ADC)
tdxd_pars$mass_frac_DXd <- tdxd_pars$DAR * tdxd_pars$MW_DXd / tdxd_pars$MW_ADC
# ≈ 8 × 718.8/148000 ≈ 0.03885

# ── PK ADC — 2 compartiments (Yin 2020 → rat) ───────────
# Valeurs humaines
CL_ADC_human_Lday <- 0.421   # L/jour
V1_ADC_human_L    <- 2.77    # L
V2_ADC_human_L    <- 5.16    # L
Q_ADC_human_Lday  <- 0.174   # L/jour

# Transposition en heures puis rat
tdxd_pars$CL_ADC <- (CL_ADC_human_Lday / 24) * scale_CL  # L/h, rat
tdxd_pars$V1_ADC <- V1_ADC_human_L            * scale_V   # L,   rat
tdxd_pars$V2_ADC <- V2_ADC_human_L            * scale_V   # L,   rat
tdxd_pars$Q_ADC  <- (Q_ADC_human_Lday  / 24) * scale_CL  # L/h, rat

# ── PK DXd (payload libre) — 1 compartiment (Yin 2020 → rat) ──
# BSA humaine standard ≈ 1.73 m² → V_DXd_human ≈ 17 × 1.73 = 29.4 L
V_DXd_human_L    <- 17 * 1.73   # L
CL_DXd_human_Lh  <- 19.2        # L/h

tdxd_pars$V_DXd  <- V_DXd_human_L   * scale_V   # L,   rat
tdxd_pars$CL_DXd <- CL_DXd_human_Lh * scale_CL  # L/h, rat

# ── Constantes mécanistiques ─────────────────────────────
# Internalisation de l'ADC dans la cellule (t½ = 46h, Vasalou 2024)
tdxd_pars$k_int  <- log(2) / 46     # h⁻¹ ≈ 0.01507

# Libération du payload — Krel TEMPS-DÉPENDANT (Yin 2020, Eq. finale)
# Krel(cycle) = k_rel_c1 × cycle^(-0.137) × (0.830 si cycle > 1)
# → −25% au cycle 2, −29% au cycle 3, −39% au cycle 10
tdxd_pars$k_rel_c1    <- 0.0159     # h⁻¹  (valeur cycle 1, Yin 2020)
tdxd_pars$krel_power  <- -0.137     # exposant puissance par cycle
tdxd_pars$krel_factor <- 0.830      # réduction additionnelle cycles > 1
tdxd_pars$interval_h  <- 21 * 24   # h    (Q3W = 504h — à ajuster si autre schéma)

# Échanges membranaires DXd (Vasalou 2024)
tdxd_pars$k_inD  <- 0.7             # h⁻¹  (entrée intracellulaire)
tdxd_pars$k_effD <- 0.7             # h⁻¹  (efflux)

# ── Paramètre PD (pour connexion future) ─────────────────
tdxd_pars$IC50_DXd_uM <- 0.31       # µM (topoisomérase I)

# ── Compartiment intracellulaire DXd ─────────────────────
# V_ic : volume intracellulaire des cellules de moelle osseuse (rat)
#   Moelle osseuse ~ 1.5% du BW → 0.25 × 0.015 = 3.75 mL
#   Fraction intracellulaire ~ 70% → ~2.6 mL ≈ 0.0026 L
#   → ratio V_DXd/V_ic ≈ 35 → accumulation intracell ~35×
tdxd_pars$V_ic <- 0.003      # L  (volume intracell. moelle, rat)

# Ratio d'accumulation à l'équilibre :
# C_DXd_ic / C_DXd = (k_inD/k_effD) × (V_DXd/V_ic)
# = 1 × (V_DXd/V_ic) car k_inD = k_effD

# ── Dommage ADN (γH2AX, modèle Fornari adapté DXd) ───────
# k_dam : taux de formation des dommages (proportionnel à E_drug)
# k_rep : taux de réparation ADN (t½ réparation ≈ 41h)
tdxd_pars$k_dam <- 0.017     # h⁻¹
tdxd_pars$k_rep <- 0.017     # h⁻¹

# ── Calibration directe FDA BLA 761139 (rat) ─────────────
# DS-8201a ne lie PAS HER2 rat → pas de TMDD
#   → CL_ADC réelle < prédiction allométrique (pas de clairance cible-médiée)
#   → Krel_rat << Krel_humain (libération DXd = clivage passif linker seulement,
#      sans endocytose récepteur-médiée)
#
# Source : Table 6, 6-Week Intermittent IV Toxicity Study of DS-8201a in Rats
#   20 mg/kg Q3W×3 : C0_ADC = 439 µg/mL,  AUC0-21d = 1776 µg·d/mL,  DXd_C0 = 0.819 ng/mL
#   60 mg/kg Q3W×3 : C0_ADC = 1400 µg/mL, AUC0-21d = 4903 µg·d/mL,  DXd_C0 = 2.49 ng/mL

# Sauvegarde des valeurs Yin 2020 / Vasalou (humain) avant remplacement
k_rel_c1_yin2020 <- tdxd_pars$k_rel_c1   # 0.0159 h⁻¹ — conservé pour tdxd_pars_human
k_int_vasalou    <- tdxd_pars$k_int       # 0.01507 h⁻¹ — conservé pour tdxd_pars_human

# ── Données FDA Table 6 (males, Day 1) ───────────────────
# C0 [µg/mL = mg/L], AUC0-21d [µg.d/mL = mg.d/L], T½ [jours]
fda_tk_rat <- list(
  list(dose_mgkg = 20, C0 = 439,  AUC21d = 1776, t_half_d = 8.07),
  list(dose_mgkg = 60, C0 = 1400, AUC21d = 4903, t_half_d = 8.52)
)

V1_v <- CL_v <- V2_v <- c()

for (d in fda_tk_rat) {
  dose_mg <- d$dose_mgkg * tdxd_pars$BW_rat

  # V1 [L] = dose [mg] / C0 [mg/L]
  V1 <- dose_mg / d$C0

  # Correction AUC0-21d → AUC0-inf : fraction éliminée = 1 - exp(-β × 504h)
  beta    <- log(2) / (d$t_half_d * 24)    # h⁻¹
  frac    <- 1 - exp(-beta * 504)
  AUC_inf <- d$AUC21d / frac               # mg.d/L

  # CL [L/h] = dose [mg] / (AUC [mg.d/L] × 24 h/d)
  CL  <- dose_mg / (AUC_inf * 24)

  # V_ss [L] = CL / β ; V2 = V_ss - V1
  Vss <- CL / beta
  V2  <- Vss - V1

  V1_v <- c(V1_v, V1); CL_v <- c(CL_v, CL); V2_v <- c(V2_v, V2)
}

V1_ADC_fda_rat  <- mean(V1_v)   # 0.01105 L  (vs allom. 0.00989)
CL_ADC_fda_rat  <- mean(CL_v)   # 1.01e-4 L/h (vs allom. 2.56e-4)
V2_ADC_fda_rat  <- mean(V2_v)   # 0.01803 L  (vs allom. 0.01843 — quasi-identique)

# Application V1, V2, CL avant le calcul de Krel (Krel dépend de V1)
tdxd_pars$V1_ADC <- V1_ADC_fda_rat
tdxd_pars$V2_ADC <- V2_ADC_fda_rat
tdxd_pars$CL_ADC <- CL_ADC_fda_rat
tdxd_pars$k_int  <- 0   # pas d'internalisation récepteur-médiée sans HER2 rat

# Krel rat calibrée [h⁻¹] — pseudo-équilibre DXd :
#   C_DXd_ss = Krel × C_ADC × mass_frac × V1_ADC / CL_DXd
#   → Krel = C_DXd [mg/L] × CL_DXd / (mass_frac × C_ADC [mg/L] × V1_ADC [L])
# NOTE : utilise tdxd_pars$V1_ADC déjà mis à jour ci-dessus
krel_fda_20  <- (0.819e-3) * tdxd_pars$CL_DXd /
                (tdxd_pars$mass_frac_DXd * 439  * tdxd_pars$V1_ADC)
krel_fda_60  <- (2.49e-3)  * tdxd_pars$CL_DXd /
                (tdxd_pars$mass_frac_DXd * 1400 * tdxd_pars$V1_ADC)
krel_fda_rat <- mean(c(krel_fda_20, krel_fda_60))
# → ~0.00119 h⁻¹  (vs Yin2020 0.01590 : facteur ~0.075)

tdxd_pars$k_rel_c1 <- krel_fda_rat   # h⁻¹ — calibré FDA rat

# ── Conversion de concentration ──────────────────────────
# C_DXd [mg/L] → C_DXd [µM] : × 1000 / MW_DXd
tdxd_pars$mgL_to_uM_DXd <- 1000 / tdxd_pars$MW_DXd  # µM per mg/L

# ── Paramètres PK humains (non scalés) ──────────────────
# Utilisés pour validation vs données FDA BLA 761139
# Source : Yin et al. 2020, PopPK DS-8201a ; FDA BLA 761139
tdxd_pars_human <- list(
  CL_ADC       = CL_ADC_human_Lday / 24,  # 0.01754 L/h
  V1_ADC       = V1_ADC_human_L,           # 2.77 L
  V2_ADC       = V2_ADC_human_L,           # 5.16 L
  Q_ADC        = Q_ADC_human_Lday  / 24,  # 0.00725 L/h
  V_DXd        = V_DXd_human_L,            # 29.41 L
  CL_DXd       = CL_DXd_human_Lh,         # 19.2 L/h
  k_int        = k_int_vasalou,         # Vasalou 2024, humain (0.01507 h⁻¹)
  k_rel_c1     = k_rel_c1_yin2020,     # Yin 2020, humain (0.0159 h⁻¹, non modifié)
  krel_power   = tdxd_pars$krel_power,
  krel_factor  = tdxd_pars$krel_factor,
  interval_h   = tdxd_pars$interval_h,
  k_inD        = tdxd_pars$k_inD,
  k_effD       = tdxd_pars$k_effD,
  V_ic         = tdxd_pars$V_ic * (70 / 0.25),  # V_ic scalé humain : 0.003 × 280 = 0.84 L
                                                 # moelle ~ 1.5% × 70 kg × 70% intracell.
  mass_frac_DXd = tdxd_pars$mass_frac_DXd,
  DAR          = tdxd_pars$DAR,
  MW_ADC       = tdxd_pars$MW_ADC,
  MW_DXd       = tdxd_pars$MW_DXd,
  mgL_to_uM_DXd = tdxd_pars$mgL_to_uM_DXd,
  IC50_DXd_uM  = tdxd_pars$IC50_DXd_uM,
  k_dam        = tdxd_pars$k_dam,
  k_rep        = tdxd_pars$k_rep
)

# ── Cibles de validation — FDA BLA 761139 (ENHERTU, 2019) ─
# Expositions géométriques moyennes à 5.4 mg/kg Q3W, régime approuvé
# Source : Clinical Pharmacology Review, BLA 761139 ; Yin et al. 2020
# ADC  : µg/mL = mg/L
# DXd  : ng/mL (plasma, DXd libre)
tdxd_fda_targets <- list(
  # ── Dose ────────────────────────────────────────────────
  dose_mgkg    = 5.4,
  BW_kg        = 70,
  Tinfu_h      = 1.5,     # perfusion 90 min (protocole clinique)
  # ── ADC sérum ───────────────────────────────────────────
  Cmax_ADC_mgL = 122,     # µg/mL ≈ mg/L  (géométrique, cycle 1, SS~cycle 3)
  AUCtau_ADC_mgLh = 20000, # µg·h/mL approx.  (Dose/CL = 378/0.01754 ≈ 21550 h·mg/L)
  # ── DXd plasma ──────────────────────────────────────────
  Cmax_DXd_ngmL = 4.4,   # ng/mL  (DXd libre plasmatique)
  Cmax_DXd_mgL  = 4.4e-3, # ng/mL → mg/L
  Cmax_DXd_uM   = (4.4e-3) * (1000 / 718.8),  # ≈ 0.00612 µM
  # ── Toxicologie rat — données hématologiques (FDA BLA Table 6) ──────────
  # DS-8201a ne lie pas HER2 rat → effets DXd-dépendants (clivage passif)
  # Étude 6 semaines, Q3W × 3 doses
  # ≥ 20 mg/kg : ↓ réticulocytes (sang)
  # ≥ 60 mg/kg : ↓ érythroblastes (BM), ↓ leucocytes/lymphocytes/neutrophiles,
  #              ↑ plaquettes (thrombocytose réactionnelle)
  # 197 mg/kg  : ↓ myélocytes (BM), ↓ monocytes, + effets non-hémato
  rat_NOAEL_mgkg         = 10,   # mg/kg (doses < 20 : pas d'effet hémato)
  rat_ret_threshold_mgkg = 20,   # ≥ 20 mg/kg → ↓ Ret
  rat_ery_threshold_mgkg = 60,   # ≥ 60 mg/kg → ↓ érythroblastes BM
  rat_plt_reactif_mgkg   = 60,   # ≥ 60 mg/kg → ↑ Plt (thrombocytose réactionnelle)
  rat_mye_threshold_mgkg = 197   # 197 mg/kg → ↓ myélocytes BM
)

# Tolérance de validation (±30% pour Cmax, ±40% pour AUC)
tdxd_fda_targets$tol_Cmax <- 0.30
tdxd_fda_targets$tol_AUC  <- 0.40

# ── Calibration PD — Slope_MEP T-DXd rat ────────────────
# Slope_MEP carboplatin (parameters_rat.R) = 2.19
# Pour T-DXd : Damage_max plus élevé (accumulation intracell DXd)
# → Slope_MEP = 1.00 calibré pour satisfaire les seuils FDA BLA :
#   MEP@20 mg/kg = -9.8%  < 10% (sous seuil histopathologique, n=4)
#   MEP@60 mg/kg = -23.7% > 20% (détectable, FDA : érythroblastes↓ à ≥60 mg/kg)
#   (scan calibrate_slope_mep.R : 2.19→1.80→1.40→1.00)
Slope_MEP_tdxd_rat <- 1.00

# ── Fonction d'administration IV (perfusion courte) ──────
make_tdxd_infusion <- function(dose_mgkg, BW_kg = 0.25,
                               Tinfu_h = 0.5, interval_h = NULL,
                               n_cycles = 1) {
  dose_mg  <- dose_mgkg * BW_kg
  rate_mgh <- dose_mg / Tinfu_h                     # mg/h

  if (is.null(interval_h) || n_cycles == 1) {
    # Dose unique
    function(t) if (t >= 0 & t < Tinfu_h) rate_mgh else 0
  } else {
    # Doses répétées
    t_starts <- seq(0, by = interval_h, length.out = n_cycles)
    function(t) {
      if (any(t >= t_starts & t < (t_starts + Tinfu_h))) rate_mgh else 0
    }
  }
}

# ── État initial ─────────────────────────────────────────
tdxd_state0 <- c(
  C_ADC1    = 0,   # ADC compartiment central        [mg/L]
  C_ADC2    = 0,   # ADC compartiment périphérique   [mg/L]
  C_DXd     = 0,   # DXd plasma                      [mg/L]
  C_DXd_ic  = 0,   # DXd intracellulaire (moelle)    [mg/L]
  Damage    = 0    # Dommages ADN (γH2AX normalisés) [sans unité]
)

# ── Résumé des paramètres ────────────────────────────────
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  PK T-DXd — RAT (calibré FDA BLA 761139 + Yin 2020)    ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("Allométrie : scale_CL=%.5f  scale_V=%.6f\n", scale_CL, scale_V))
cat(sprintf("Calibration FDA rat (Table 6) : V1×%.2f  CL×%.3f  Krel×%.3f  k_int=0\n\n",
            V1_ADC_fda_rat / (V1_ADC_human_L * scale_V),
            CL_ADC_fda_rat / (CL_ADC_human_Lday/24 * scale_CL),
            krel_fda_rat / k_rel_c1_yin2020))
cat("── ADC (2-compartiments) ──\n")
cat(sprintf("  CL_ADC = %.4e L/h  [FDA]  (allom.=%.4e)\n",
            tdxd_pars$CL_ADC, CL_ADC_human_Lday/24 * scale_CL))
cat(sprintf("  V1_ADC = %.5f L    [FDA]  (allom.=%.5f)\n",
            tdxd_pars$V1_ADC, V1_ADC_human_L * scale_V))
cat(sprintf("  V2_ADC = %.5f L    [FDA]  (allom.=%.5f)\n",
            tdxd_pars$V2_ADC, V2_ADC_human_L * scale_V))
cat(sprintf("  Q_ADC  = %.4e L/h  (allom., human: %.5f L/h)\n\n",
            tdxd_pars$Q_ADC, Q_ADC_human_Lday/24))
cat("── DXd payload (1-compartiment) ──\n")
cat(sprintf("  CL_DXd = %.5f L/h  (human: %.1f L/h)\n",
            tdxd_pars$CL_DXd, CL_DXd_human_Lh))
cat(sprintf("  V_DXd  = %.5f L    (human: %.1f L)\n\n",
            tdxd_pars$V_DXd, V_DXd_human_L))
cat("── Constantes mécanistiques ──\n")
cat(sprintf("  k_int  = %.5f h⁻¹  [RAT=0, pas de liaison HER2]  (Vasalou humain = %.5f)\n",
            tdxd_pars$k_int, k_int_vasalou))
cat(sprintf("  Krel   = %.6f × Cycle^(%.3f) × (%.3f si Cycle>1)  [FDA rat calibré]\n",
            tdxd_pars$k_rel_c1, tdxd_pars$krel_power, tdxd_pars$krel_factor))
cat(sprintf("         (Yin2020 humain = %.4f h⁻¹ → facteur rat/humain = %.3f)\n",
            k_rel_c1_yin2020, tdxd_pars$k_rel_c1 / k_rel_c1_yin2020))
cat(sprintf("         Cycle1=%.6f  Cycle2=%.6f  Cycle3=%.6f h⁻¹\n",
            tdxd_pars$k_rel_c1,
            tdxd_pars$k_rel_c1 * 2^tdxd_pars$krel_power * tdxd_pars$krel_factor,
            tdxd_pars$k_rel_c1 * 3^tdxd_pars$krel_power * tdxd_pars$krel_factor))
cat(sprintf("  k_inD  = k_effD = %.1f h⁻¹\n", tdxd_pars$k_inD))
cat(sprintf("  DAR    = %d  |  mass_frac_DXd = %.5f\n",
            tdxd_pars$DAR, tdxd_pars$mass_frac_DXd))
cat(sprintf("  IC50_DXd = %.2f µM\n",    tdxd_pars$IC50_DXd_uM))
cat(sprintf("  V_ic     = %.4f L  (ratio V_DXd/V_ic = %.0f → accum. ~%.0fx)\n",
            tdxd_pars$V_ic,
            tdxd_pars$V_DXd / tdxd_pars$V_ic,
            tdxd_pars$V_DXd / tdxd_pars$V_ic))
cat(sprintf("  k_dam    = %.4f h⁻¹  k_rep = %.4f h⁻¹\n\n",
            tdxd_pars$k_dam, tdxd_pars$k_rep))
