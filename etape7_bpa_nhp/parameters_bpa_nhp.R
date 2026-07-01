############################################################
# parameters_bpa_nhp.R
# Parametres PK/PD -- anti-FLT3_chBPA-STINGa20 ADC -- SINGE (cynomolgus)
#
# PK ADC : allometrie souris (Shiny, 1-cmt) -> cynomolgus (4 kg)
#          CL, Q  : BW^0.75  (Boxenbaum 1982)
#          V1, V2 : BW^1.0
#
# PD     : parametres Fornari humain (NHP ~ humain physiologiquement)
# Slopes : identiques humain (IC50 mesurees sur CFU humaines)
############################################################

source("../etape2_carboplatin_humain/parameters_human.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("../etape5_tdxd_humain/parameters_tdxd_human.R")

# -- IC50 reference T-DXd [nM] ------------------------------------------------
# C[nM] = C[ug/mL] * 1e6 / MW[g/mol]  (ex: 27.3 ug/mL / 148000 g/mol * 1e6 = 184.5 nM)
IC50_MEP_tdxd_nM <- 27.3 / 148000 * 1e6   # 184.5 nM
IC50_CMP_tdxd_nM <- 28.1 / 148000 * 1e6   # 189.9 nM

# -- IC50 BPA (CFU humaines, en nM) -------------------------------------------
IC50_MEP_BPA_nM  <- 155
IC50_CMP_BPA_s1  <- 7.8
IC50_CMP_BPA_s2  <- 0.015

# -- Slopes BPA ---------------------------------------------------------------
Slope_MEP_BPA_s1 <- Slope_MEP_sensitive_tdxd * (IC50_MEP_tdxd_nM / IC50_MEP_BPA_nM)
Slope_CMP_BPA_s1 <- Slope_sensitive_tdxd     * (IC50_CMP_tdxd_nM / IC50_CMP_BPA_s1)
Slope_MEP_BPA_s2 <- Slope_MEP_sensitive_tdxd * (IC50_MEP_tdxd_nM / IC50_MEP_BPA_nM)
Slope_CMP_BPA_s2 <- Slope_sensitive_tdxd     * (IC50_CMP_tdxd_nM / IC50_CMP_BPA_s2)

# -- PK BPA : valeurs souris Shiny (1-cmt) ------------------------------------
BW_MOUSE_KG <- 0.025
BW_NHP_KG   <- 4.0     # cynomolgus moyen

BPA_CL_mouse_per_kg <- 0.000422   # L/kg/h  (Shiny, %RSE=22.7%)
BPA_V1_mouse_per_kg <- 0.0787     # L/kg    (Shiny, %RSE=34.7%)

# -- Allometrie souris -> singe -----------------------------------------------
allo_CL_nhp <- BW_MOUSE_KG * (BW_NHP_KG / BW_MOUSE_KG)^0.75
allo_V_nhp  <- BW_NHP_KG

BPA_CL_nhp <- BPA_CL_mouse_per_kg * allo_CL_nhp   # [L/h]
BPA_V1_nhp <- BPA_V1_mouse_per_kg * allo_V_nhp    # [L]

cat(sprintf("Allometrie souris->singe (BW_mouse=%.3f kg -> BW_NHP=%.0f kg) :\n",
            BW_MOUSE_KG, BW_NHP_KG))
cat(sprintf("  facteur CL = x%.2f  |  facteur V = x%.0f\n", allo_CL_nhp, allo_V_nhp))
cat(sprintf("  CL  : %.6f L/kg/h (souris) -> %.6f L/h (singe)\n",
            BPA_CL_mouse_per_kg, BPA_CL_nhp))
cat(sprintf("  V1  : %.4f L/kg   (souris) -> %.4f L   (singe)\n",
            BPA_V1_mouse_per_kg, BPA_V1_nhp))
cat(sprintf("  T½ singe : %.1f h  (%.1f jours)\n",
            log(2) / (BPA_CL_nhp / BPA_V1_nhp),
            log(2) / (BPA_CL_nhp / BPA_V1_nhp) / 24))

# -- Parametres PK NHP --------------------------------------------------------
bpa_pars_nhp        <- tdxd_pars_hu
bpa_pars_nhp$CL_ADC <- BPA_CL_nhp
bpa_pars_nhp$V1_ADC <- BPA_V1_nhp
bpa_pars_nhp$V2_ADC <- BPA_V1_nhp   # 1-cmt effectif
bpa_pars_nhp$Q_ADC  <- 0

# -- IIV ----------------------------------------------------------------------
omega_CL        <- 0.35
omega_V1        <- 0.20
omega_Slope_CMP <- 0.33
omega_Slope_MEP <- 0.33

# -- Protocole NHP ------------------------------------------------------------
# Etude de toxicite preclinique typique : Q3W x 4 cycles
BW_KG      <- BW_NHP_KG
DOSE_MGKG  <- 5.4
TINFU_H    <- 1.5
INTERVAL_H <- 21 * 24
N_CYCLES   <- 4

# -- Resume -------------------------------------------------------------------
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║  PK/PD BPA ADC -- SINGE CYNOMOLGUS (prospectif)     ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")
cat(sprintf("IC50 MEP BPA     = %.1f nM\n", IC50_MEP_BPA_nM))
cat(sprintf("IC50 CMP BPA s1  = %.2f nM  (sous-traitant)\n", IC50_CMP_BPA_s1))
cat(sprintf("IC50 CMP BPA s2  = %.4f nM  (interne)\n", IC50_CMP_BPA_s2))
cat(sprintf("Slope_MEP s1=s2  = %.3f\n", Slope_MEP_BPA_s1))
cat(sprintf("Slope_CMP s1     = %.1f\n", Slope_CMP_BPA_s1))
cat(sprintf("Slope_CMP s2     = %.0f\n", Slope_CMP_BPA_s2))
cat(sprintf("Protocole        : %.1f mg/kg Q3W x%d (NHP, BW=%.0f kg)\n",
            DOSE_MGKG, N_CYCLES, BW_KG))
