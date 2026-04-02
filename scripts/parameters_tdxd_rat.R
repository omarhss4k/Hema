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

# Libération du payload intracellulaire (Vasalou 2024)
tdxd_pars$k_rel  <- 0.0159          # h⁻¹

# Échanges membranaires DXd (Vasalou 2024)
tdxd_pars$k_inD  <- 0.7             # h⁻¹  (entrée intracellulaire)
tdxd_pars$k_effD <- 0.7             # h⁻¹  (efflux)

# ── Paramètre PD (pour connexion future) ─────────────────
tdxd_pars$IC50_DXd_uM <- 0.31       # µM (topoisomérase I)

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

# ── État initial (PK seul) ───────────────────────────────
tdxd_state0 <- c(
  C_ADC1 = 0,   # ADC compartiment central    [mg/L]
  C_ADC2 = 0,   # ADC compartiment périphérique [mg/L]
  C_DXd  = 0    # DXd plasma                  [mg/L]
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
cat(sprintf("  k_rel  = %.4f h⁻¹\n", tdxd_pars$k_rel))
cat(sprintf("  k_inD  = k_effD = %.1f h⁻¹\n", tdxd_pars$k_inD))
cat(sprintf("  DAR    = %d  |  mass_frac_DXd = %.5f\n",
            tdxd_pars$DAR, tdxd_pars$mass_frac_DXd))
cat(sprintf("  IC50_DXd = %.2f µM\n\n", tdxd_pars$IC50_DXd_uM))
