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

# ── Conversion de concentration ──────────────────────────
# C_DXd [mg/L] → C_DXd [µM] : × 1000 / MW_DXd
tdxd_pars$mgL_to_uM_DXd <- 1000 / tdxd_pars$MW_DXd  # µM per mg/L

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
cat("║  PK T-DXd — RAT (allométrie depuis Yin 2020)            ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("Allométrie : scale_CL=%.5f  scale_V=%.6f\n\n", scale_CL, scale_V))
cat("── ADC (2-compartiments) ──\n")
cat(sprintf("  CL_ADC = %.4e L/h  (human: %.5f L/h)\n",
            tdxd_pars$CL_ADC, CL_ADC_human_Lday/24))
cat(sprintf("  V1_ADC = %.5f L    (human: %.2f L)\n",
            tdxd_pars$V1_ADC, V1_ADC_human_L))
cat(sprintf("  V2_ADC = %.5f L    (human: %.2f L)\n",
            tdxd_pars$V2_ADC, V2_ADC_human_L))
cat(sprintf("  Q_ADC  = %.4e L/h  (human: %.5f L/h)\n\n",
            tdxd_pars$Q_ADC, Q_ADC_human_Lday/24))
cat("── DXd payload (1-compartiment) ──\n")
cat(sprintf("  CL_DXd = %.5f L/h  (human: %.1f L/h)\n",
            tdxd_pars$CL_DXd, CL_DXd_human_Lh))
cat(sprintf("  V_DXd  = %.5f L    (human: %.1f L)\n\n",
            tdxd_pars$V_DXd, V_DXd_human_L))
cat("── Constantes mécanistiques ──\n")
cat(sprintf("  k_int  = %.5f h⁻¹  (t½ internalisation = %.1f h)\n",
            tdxd_pars$k_int, log(2)/tdxd_pars$k_int))
cat(sprintf("  Krel   = %.4f × Cycle^(%.3f) × (%.3f si Cycle>1)  [Yin 2020, temps-dep.]\n",
            tdxd_pars$k_rel_c1, tdxd_pars$krel_power, tdxd_pars$krel_factor))
cat(sprintf("         Cycle1=%.4f  Cycle2=%.4f  Cycle3=%.4f h⁻¹\n",
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
