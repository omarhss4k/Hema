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

# -- PK BPA : coller ici les parametres exportes depuis le Shiny ---------------
#
#   1. Ouvrir le Shiny PK (shiny_pk/app.R)
#   2. Charger les donnees PK BPA (concentration vs temps)
#   3. Onglet "Modelisation" -> ajuster le modele 2 compartiments
#   4. Lire les valeurs CL, V1, V2, Q dans le tableau "Parametres"
#   5. Remplacer les valeurs ci-dessous :
#
#   UNITES ATTENDUES : CL en L/h, V en L  (BW=70 kg)
#   Si le Shiny sort en L/kg/h et L/kg, multiplier par BW_KG=70

BPA_CL  <- tdxd_pars_hu$CL_ADC   # REMPLACER par valeur Shiny  [L/h]
BPA_V1  <- tdxd_pars_hu$V1_ADC   # REMPLACER par valeur Shiny  [L]
BPA_V2  <- tdxd_pars_hu$V2_ADC   # REMPLACER par valeur Shiny  [L]
BPA_Q   <- tdxd_pars_hu$Q_ADC    # REMPLACER par valeur Shiny  [L/h]

# Exemple une fois les donnees disponibles :
# BPA_CL <- 0.012    # L/h
# BPA_V1 <- 3.5      # L
# BPA_V2 <- 6.2      # L
# BPA_Q  <- 0.008    # L/h

bpa_pars        <- tdxd_pars_hu   # herite tous les autres params (payload, k_int, etc.)
bpa_pars$CL_ADC <- BPA_CL         # <- seuls ces 4 parametres changent
bpa_pars$V1_ADC <- BPA_V1
bpa_pars$V2_ADC <- BPA_V2
bpa_pars$Q_ADC  <- BPA_Q

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
