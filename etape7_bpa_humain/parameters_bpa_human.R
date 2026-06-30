############################################################
# parameters_bpa_human.R
# Parametres PK/PD -- anti-FLT3_chBPA-STINGa20 ADC -- HUMAIN
#
# PK ADC : placeholder T-DXd (Yin 2020) -- a remplacer par PK BPA
#           lorsque les donnees seront disponibles
# Slopes  : calibres depuis IC50 CFU humaines via ratio T-DXd
#           Slope_BPA = Slope_tdxd x (IC50_tdxd_nM / IC50_BPA_nM)
# IC50 reference T-DXd : 27.3 ug/mL (erythroid), 28.1 ug/mL (myeloid)
#                     -> 184 nM (erythroid), 190 nM (myeloid)  [MW=148000]
#
# IC50 BPA (CFU humaines) :
#   Simu 1 (sous-traitant) : IC50_MEP=155 nM, IC50_CMP=7.8 nM
#   Simu 2 (interne)       : IC50_MEP=155 nM, IC50_CMP=0.015 nM
############################################################

source("../etape2_carboplatin_humain/parameters_human.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("../etape5_tdxd_humain/parameters_tdxd_human.R")

# -- IC50 reference T-DXd converties en nM (MW_ADC = 148000 g/mol) -----------
IC50_MEP_tdxd_nM <- 27.3e3 / 148000 * 1e6   # 184.5 nM
IC50_CMP_tdxd_nM <- 28.1e3 / 148000 * 1e6   # 189.9 nM

# -- IC50 BPA (CFU humaines, en nM) -------------------------------------------
IC50_MEP_BPA_nM  <- 155     # erythroid (sous-traitant -- interne = ND)
IC50_CMP_BPA_s1  <- 7.8     # myeloid Simu 1 (sous-traitant)
IC50_CMP_BPA_s2  <- 0.015   # myeloid Simu 2 (interne)

# -- Slopes BPA calcules par ratio IC50 ---------------------------------------
# Slope_BPA = Slope_tdxd_human * (IC50_tdxd_nM / IC50_BPA_nM)
# Rappel : Slope_sensitive_tdxd (humain) ~ 59 000 (calibre FDA BLA)
#          Slope_MEP_sensitive_tdxd      ~ 1 000
#
# Simu 1 (sous-traitant) :
#   Slope_CMP = 59000 * (190/7.8)   ~ 1 440 000  (24x plus potent que T-DXd)
#   Slope_MEP = 1000  * (184/155)   ~     1 190
Slope_MEP_BPA_s1 <- Slope_MEP_sensitive_tdxd * (IC50_MEP_tdxd_nM / IC50_MEP_BPA_nM)
Slope_CMP_BPA_s1 <- Slope_sensitive_tdxd     * (IC50_CMP_tdxd_nM / IC50_CMP_BPA_s1)

# Simu 2 (interne) :
#   Slope_CMP = 59000 * (190/0.015) ~ 747 000 000 (EXTREME : 12 667x plus potent)
#   -> kill_CMP = 1 garanti des la premiere dose pour 100% des patients
Slope_MEP_BPA_s2 <- Slope_MEP_sensitive_tdxd * (IC50_MEP_tdxd_nM / IC50_MEP_BPA_nM)
Slope_CMP_BPA_s2 <- Slope_sensitive_tdxd     * (IC50_CMP_tdxd_nM / IC50_CMP_BPA_s2)

# -- PK BPA : parametres souris issus du Shiny ---------------------------------
#
#   1. Shiny PK -> charger donnees BPA souris -> Modelisation 2-cmt
#   2. Lire CL, V1, V2, Q dans le tableau Parametres
#   3. Les unites Shiny sont L/kg/h (CL, Q) et L/kg (V1, V2) -- souris
#   4. Renseigner ci-dessous, l'allometrie est appliquee automatiquement

BW_MOUSE_KG <- 0.025   # 25 g (poids moyen souris)

# ---- Valeurs souris issues du Shiny [L/kg/h et L/kg] ----
BPA_CL_mouse_per_kg  <- tdxd_pars_hu$CL_ADC / BW_MOUSE_KG   # REMPLACER [L/kg/h]
BPA_V1_mouse_per_kg  <- tdxd_pars_hu$V1_ADC / BW_MOUSE_KG   # REMPLACER [L/kg]
BPA_V2_mouse_per_kg  <- tdxd_pars_hu$V2_ADC / BW_MOUSE_KG   # REMPLACER [L/kg]
BPA_Q_mouse_per_kg   <- tdxd_pars_hu$Q_ADC  / BW_MOUSE_KG   # REMPLACER [L/kg/h]

# Exemple :
# BPA_CL_mouse_per_kg <- 0.85   # L/kg/h
# BPA_V1_mouse_per_kg <- 0.077  # L/kg
# BPA_V2_mouse_per_kg <- 0.15   # L/kg
# BPA_Q_mouse_per_kg  <- 0.12   # L/kg/h

# ---- Allometrie souris -> humain ------------------------------------------
# Regles standard (Boxenbaum 1982) :
#   CL, Q  : echelle BW^0.75  (metabolisme energetique)
#   V1, V2 : echelle BW^1.0   (volumes de distribution proportionnels au poids)
#
# Formule :
#   CL_human [L/h] = CL_mouse [L/kg/h] x BW_mouse x (BW_human/BW_mouse)^0.75
#   V_human  [L]   = V_mouse  [L/kg]   x BW_human
#
# Interpretation :
#   k10 = CL/V -> k10_human = k10_mouse x (BW_mouse/BW_human)^0.25
#   Avec BW_mouse=0.025 kg, BW_human=70 kg :
#   k10_human = k10_mouse x (0.025/70)^0.25 = k10_mouse x 0.133
#   -> T½ humain ~ 7.5x plus long que souris a meme dose mg/kg
#   -> AUC humain ~ 7.5x plus elevee -> toxicite potentiellement plus marquee

allo_CL <- BW_MOUSE_KG * (BW_KG / BW_MOUSE_KG)^0.75   # facteur allometrique CL/Q
allo_V  <- BW_KG                                         # facteur allometrique V

BPA_CL_hu <- BPA_CL_mouse_per_kg * allo_CL   # [L/h]
BPA_V1_hu <- BPA_V1_mouse_per_kg * allo_V    # [L]
BPA_V2_hu <- BPA_V2_mouse_per_kg * allo_V    # [L]
BPA_Q_hu  <- BPA_Q_mouse_per_kg  * allo_CL   # [L/h]

cat(sprintf("Allometrie souris->humain (BW_mouse=%.3f kg -> BW_human=%.0f kg) :\n",
            BW_MOUSE_KG, BW_KG))
cat(sprintf("  facteur CL/Q = x%.1f  |  facteur V = x%.0f\n", allo_CL, allo_V))
cat(sprintf("  CL  : %.4f L/kg/h (souris) -> %.4f L/h (humain)\n",
            BPA_CL_mouse_per_kg, BPA_CL_hu))
cat(sprintf("  V1  : %.4f L/kg   (souris) -> %.3f L   (humain)\n",
            BPA_V1_mouse_per_kg, BPA_V1_hu))
cat(sprintf("  V2  : %.4f L/kg   (souris) -> %.3f L   (humain)\n",
            BPA_V2_mouse_per_kg, BPA_V2_hu))
cat(sprintf("  Q   : %.4f L/kg/h (souris) -> %.4f L/h (humain)\n",
            BPA_Q_mouse_per_kg, BPA_Q_hu))
cat(sprintf("  T½ alpha humain estimee : %.1f h\n",
            log(2) / ((BPA_CL_hu/BPA_V1_hu + BPA_Q_hu/BPA_V1_hu + BPA_Q_hu/BPA_V2_hu +
            sqrt((BPA_CL_hu/BPA_V1_hu + BPA_Q_hu/BPA_V1_hu + BPA_Q_hu/BPA_V2_hu)^2 -
                 4*BPA_CL_hu/BPA_V1_hu * BPA_Q_hu/BPA_V2_hu)) / 2)))

bpa_pars        <- tdxd_pars_hu    # herite payload, k_int, k_rel, etc.
bpa_pars$CL_ADC <- BPA_CL_hu
bpa_pars$V1_ADC <- BPA_V1_hu
bpa_pars$V2_ADC <- BPA_V2_hu
bpa_pars$Q_ADC  <- BPA_Q_hu

# -- IIV (log-normal) ---------------------------------------------------------
omega_CL        <- 0.35   # Yin 2020
omega_V1        <- 0.20   # Yin 2020
omega_Slope_CMP <- 0.33   # Fornari Table S4
omega_Slope_MEP <- 0.33   # Fornari Table S4

# -- Protocole ----------------------------------------------------------------
# A definir selon le composé BPA -- placeholder Q3W identique T-DXd
DOSE_MGKG  <- 5.4
TINFU_H    <- 1.5
INTERVAL_H <- 21 * 24
N_CYCLES   <- 6
BW_KG      <- 70

# -- Resume ------------------------------------------------------------------
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║  PK/PD BPA ADC -- HUMAIN (prospectif)               ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")
cat(sprintf("IC50 MEP BPA     = %.1f nM\n", IC50_MEP_BPA_nM))
cat(sprintf("IC50 CMP BPA s1  = %.2f nM  (sous-traitant)\n", IC50_CMP_BPA_s1))
cat(sprintf("IC50 CMP BPA s2  = %.4f nM  (interne)\n", IC50_CMP_BPA_s2))
cat("\n-- Slopes calcules --\n")
cat(sprintf("Slope_MEP s1 = s2 = %.3f\n", Slope_MEP_BPA_s1))
cat(sprintf("Slope_CMP s1      = %.1f\n", Slope_CMP_BPA_s1))
cat(sprintf("Slope_CMP s2      = %.0f  (toxicite myeloide extreme attendue)\n",
            Slope_CMP_BPA_s2))
