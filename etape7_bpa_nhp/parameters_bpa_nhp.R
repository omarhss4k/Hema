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

# -- Neutraliser la decroissance de relargage cycle-dependante T-DXd ----------
# krel_power=-0.137 et krel_factor=0.830 (Yin 2020) sont calibres sur la
# cinetique de clivage du linker T-DXd (GGFG) -- non pertinents pour BPA.
# Sans cette neutralisation, Krel(cycle) decroit -> chaque cycle libere
# moins de payload -> toxicite artificiellement attenuee des le cycle 2.
bpa_pars_nhp$krel_power  <- 0
bpa_pars_nhp$krel_factor <- 1.0

# -- Neutraliser k_int (internalisation mediee par la cible T-DXd = HER2) -----
# k_int=0.01507 h-1 (Vasalou 2024) represente la degradation de l'ADC via
# liaison/internalisation HER2 -- sans rapport avec la cible BPA (FLT3),
# pour laquelle aucune donnee d'internalisation n'est disponible.
# BUG MAJEUR : ce terme s'ajoute directement dans l'ODE
#   dC_ADC1 = rate_in/V1 - (CL_ADC/V1 + k_int) * C1
# et domine totalement la clearance allometrique du BPA (T½ reelle chute
# de ~19j a ~2j chez le singe), ecrasant Damage a quasi-zero des le cycle 2
# et masquant toute toxicite cumulative entre cycles Q3W.
bpa_pars_nhp$k_int <- 0

# -- Plafonner le kill via modele Emax (au lieu du lineaire Slope*D_kill) -----
# Le modele lineaire par defaut (pkpd_tdxd_rat.R, branche sans ED50_kill) n'a
# AUCUN plafond : kill_CMP = Slope_CMP * D_kill peut depasser 1 des que Damage
# franchit le seuil D0, ce qui rend (1-kill_CMP) negatif -> ablation totale
# instantanee. Avec Slope_CMP~1439, la moindre dose franchissant D0 produit
# une "falaise" tout-ou-rien (ex. observe : 0.1 mg/kg = G0 partout, 0.25 mg/kg
# = G4 chez 77% des animaux), non realiste biologiquement.
#
# En fixant ED50_kill, l'ODE bascule sur la branche Emax deja codee :
#   Emax_X   = min(1, Slope_X * ED50_kill)   (calcule automatiquement)
#   kill_X   = Emax_X * D_kill / (ED50_kill + D_kill)   (borne a Emax_X <= 1)
# Slope_CMP (1439) sature Emax_CMP=1 quelle que soit la dose testee (ablation
# totale possible a dose suffisante), tandis que Slope_MEP (1.19) plus faible
# donne Emax_MEP=min(1, 1.19*ED50_kill) < 1 -- preserve le ratio de puissance
# CMP/MEP calibre via IC50 tout en bornant chaque compartiment a un effet
# maximal plausible.
# ED50_kill=0.10 choisi ~ ordre de grandeur du D_kill observe autour de
# 0.4-0.5 mg/kg (milieu de la gamme de dose testee) -- HYPOTHESE a
# recalibrer des que des donnees dose-toxicite BPA reelles seront dispo.
bpa_pars_nhp$ED50_kill <- 0.10

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
