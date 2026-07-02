############################################################
# bpa_build_params.R
# BUILDER UNIQUE et TRACÉ des paramètres PK/PD du BPA ADC
# (anti-FLT3_chBPA-STINGa20), pour SINGE et HUMAIN.
#
# OBJECTIF : supprimer l'heritage silencieux "bpa_pars <- tdxd_pars_hu".
# Ici, CHAQUE parametre issu du T-DXd recoit un STATUT explicite :
#
#   [REMPLACE]   valeur BPA propre (PK Shiny + allometrie, Slopes ratio IC50)
#   [CORRIGE]    parametre T-DXd non pertinent pour BPA -> neutralise (bug corrige)
#   [INACTIF]    calcule par l'ODE mais SANS effet sur le PD tant que
#                use_ADC_driver=TRUE (sous-systeme DXd deconnecte du Damage)
#   [HYPOTHESE]  valeur T-DXd conservee FAUTE DE DONNEE BPA -> a recalibrer
#
# Usage :
#   source("../shared/bpa_build_params.R")
#   p <- build_bpa_pars(BW_target_kg = 4,  IC50_CMP_nM = 7.8)    # singe s1
#   p <- build_bpa_pars(BW_target_kg = 70, IC50_CMP_nM = 0.015)  # humain s2
#
# Prerequis charges AVANT : parameters_tdxd_human.R (fournit tdxd_pars_hu),
#   parameters_FORNARI_CORRECT.R (Slope_*_sensitive_tdxd, init_pars/state).
############################################################

# -- Donnees PK souris (Shiny, modele 1-cmt retenu) ---------------------------
BPA_CL_mouse_per_kg <- 0.000422   # L/kg/h  (%RSE=22.7%)
BPA_V1_mouse_per_kg <- 0.0787     # L/kg    (%RSE=34.7%)
BW_MOUSE_KG         <- 0.025

# -- IC50 reference T-DXd (nM) : conversion ug/mL CORRECTE (pas de x1000) ------
#    C[nM] = C[ug/mL] * 1e6 / MW[g/mol]
IC50_MEP_tdxd_nM_REF <- 27.3 / 148000 * 1e6   # 184.5 nM
IC50_CMP_tdxd_nM_REF <- 28.1 / 148000 * 1e6   # 189.9 nM
IC50_MEP_BPA_nM_REF  <- 155                    # erythroide (commun s1/s2)

build_bpa_pars <- function(BW_target_kg,
                           IC50_CMP_nM,
                           ED50_kill   = 0.10,
                           verbose     = TRUE) {

  # ══════════════════════════════════════════════════════════
  # [REMPLACE] PK ADC : allometrie souris -> cible (Boxenbaum 1982)
  #   CL, Q  : BW^0.75   |   V1, V2 : BW^1.0
  # ══════════════════════════════════════════════════════════
  allo_CL <- BW_MOUSE_KG * (BW_target_kg / BW_MOUSE_KG)^0.75
  allo_V  <- BW_target_kg

  CL_ADC <- BPA_CL_mouse_per_kg * allo_CL   # [L/h]
  V1_ADC <- BPA_V1_mouse_per_kg * allo_V    # [L]
  # modele 1-cmt effectif : V2=V1, Q~0 (2-cmt non identifiable, n=2 souris)
  V2_ADC <- V1_ADC
  Q_ADC  <- 0

  # ══════════════════════════════════════════════════════════
  # [REMPLACE] Slopes BPA : ratio de puissance IC50 vs T-DXd
  #   Slope_BPA = Slope_tdxd * (IC50_tdxd / IC50_BPA)
  # ══════════════════════════════════════════════════════════
  Slope_MEP <- Slope_MEP_sensitive_tdxd * (IC50_MEP_tdxd_nM_REF / IC50_MEP_BPA_nM_REF)
  Slope_CMP <- Slope_sensitive_tdxd     * (IC50_CMP_tdxd_nM_REF / IC50_CMP_nM)
  Slope_MPP <- init_pars$Slope_MPP      # herite Fornari (progeniteur amont)

  # ══════════════════════════════════════════════════════════
  # Construction de la liste, en repartant de tdxd_pars_hu MAIS
  # en ecrasant explicitement chaque champ selon son statut.
  # ══════════════════════════════════════════════════════════
  p <- tdxd_pars_hu

  # [REMPLACE] -- PK propre BPA
  p$CL_ADC <- CL_ADC;  p$V1_ADC <- V1_ADC
  p$V2_ADC <- V2_ADC;  p$Q_ADC  <- Q_ADC

  # [REMPLACE] -- Slopes propres BPA
  p$Slope_MEP <- Slope_MEP
  p$Slope_CMP <- Slope_CMP
  p$Slope_MPP <- Slope_MPP

  # [CORRIGE] -- k_int : internalisation mediee HER2 (Vasalou) -> nulle pour FLT3
  p$k_int <- 0

  # [CORRIGE] -- decroissance de relargage cycle-dependante (linker GGFG T-DXd)
  p$krel_power  <- 0
  p$krel_factor <- 1.0

  # [AJOUTE/HYPOTHESE] -- kill Emax borne (evite l'ablation "falaise" du modele
  #   lineaire non plafonne). ED50_kill NON calibre sur donnees BPA -> hypothese.
  p$ED50_kill <- ED50_kill

  # ---- Parametres HERITES conserves (documentation du statut) ----
  #
  # [HYPOTHESE ACTIVE] IC50_ADC_ugmL = 27.70 (T-DXd, assay PFB-10)
  #   -> pilote E_drug = C_ADC1/(IC50_ADC_ugmL + C_ADC1) car use_ADC_driver=TRUE.
  #   C'est l'echelle concentration->Damage. Valeur T-DXd faute d'IC50 ADC-entier
  #   BPA en ug/mL ; la potence differentielle BPA est portee par les Slopes.
  #   >>> Principale hypothese structurelle a challenger. <<<
  #   (laisse tel quel : p$IC50_ADC_ugmL inchange)
  #
  # [HYPOTHESE] k_dam = k_rep = 0.017 h-1 (cinetique dommage/reparation ADN,
  #   Fornari generique -- pas de donnee BPA). Damage_threshold = 0.05 (seuil D0
  #   calibre T-DXd). (laisses tels quels)
  #
  # [HERITE structurel] use_ADC_driver = TRUE (choix T-DXd : driver = C_ADC1).
  #
  # [INACTIF tant que use_ADC_driver=TRUE] tout le sous-systeme DXd est calcule
  #   par l'ODE mais N'AFFECTE PAS le Damage :
  #     k_rel_c1, V_DXd, CL_DXd, k_inD, k_effD, V_ic, IC50_DXd_uM,
  #     DAR, MW_ADC, MW_DXd, mass_frac_DXd
  #   -> valeurs T-DXd sans impact sur la prediction ; a re-specifier seulement
  #      si on bascule sur le driver DXd intracellulaire (use_ADC_driver=FALSE).

  # -- Tracabilite -----------------------------------------------------------
  if (verbose) {
    t_half <- log(2) / (CL_ADC / V1_ADC)
    cat(sprintf("\n[build_bpa_pars] BW_cible=%.0f kg | IC50_CMP=%.4g nM | ED50_kill=%.3g\n",
                BW_target_kg, IC50_CMP_nM, ED50_kill))
    cat(sprintf("  [REMPLACE] CL=%.5f L/h  V1=%.4f L  (T1/2=%.1f h = %.1f j)\n",
                CL_ADC, V1_ADC, t_half, t_half/24))
    cat(sprintf("  [REMPLACE] Slope_CMP=%.1f  Slope_MEP=%.2f\n", Slope_CMP, Slope_MEP))
    cat(sprintf("  [CORRIGE ] k_int=0  krel_power=0  krel_factor=1.0\n"))
    cat(sprintf("  [HYPOTHESE] ED50_kill=%.3g  IC50_ADC_ugmL=%.2f (T-DXd)  D0=%.3g\n",
                ED50_kill, p$IC50_ADC_ugmL, p$Damage_threshold))
    cat(sprintf("  [INACTIF ] sous-systeme DXd (k_rel,V_DXd,CL_DXd,...) : use_ADC_driver=%s\n",
                p$use_ADC_driver))
  }

  p
}
