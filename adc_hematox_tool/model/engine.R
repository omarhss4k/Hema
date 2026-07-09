############################################################
# engine.R -- moteur de l'outil ADC-HematoTox
# Charge le modele et fournit les fonctions de haut niveau :
#   load_model(), build_pars(), simulate(), ctcae_grade()
# Ne pas modifier pour un usage courant (voir config/compounds.R).
############################################################

# -- charge tout le modele (ODE + parametres biologiques) --
# Utilise source(..., chdir=TRUE) : le repertoire courant est change
# TEMPORAIREMENT le temps de chaque source (pour resoudre les sous-source
# en chemins locaux) puis restaure automatiquement -> pas de setwd() manuel,
# pas de risque de laisser la session RStudio dans un mauvais repertoire.
load_model <- function(model_dir = "model") {
  suppressMessages(suppressWarnings({
    library(deSolve)
    source(file.path(model_dir, "parameters_baseline.R"),  chdir = TRUE)
    source(file.path(model_dir, "parameters_fornari.R"),   chdir = TRUE)
    source(file.path(model_dir, "parameters_adc.R"),       chdir = TRUE)
    source(file.path(model_dir, "parameters_infusion.R"),  chdir = TRUE)
    source(file.path(model_dir, "model_ode.R"),            chdir = TRUE)
  }))
  invisible(TRUE)
}

# -- construit le jeu de parametres pour un composE + une IC50 + la calibration --
# cp    : liste compound (PK, IC50, BW...) depuis config/compounds.R
# ic50  : IC50 myeloide du composE (nM)
# calib : liste CALIB (PD calibre + ancres de potency)
# ic50_ref : IC50 myeloide de la REFERENCE (ancre de translation)
build_pars <- function(cp, ic50, calib, ic50_ref) {
  p <- c(init_pars, tdxd_pars_hu)
  p$CL_ADC <- cp$CL; p$V1_ADC <- cp$V1; p$Q_ADC <- cp$Q; p$V2_ADC <- cp$V2
  p$k_int <- 0; p$krel_power <- 0; p$krel_factor <- 1
  p$use_ADC_driver <- TRUE; p$Damage_threshold <- 0
  # potency translatee : IC50_ADC = echelle_ref x (IC50_composE / IC50_ref)
  p$IC50_ADC_ugmL <- calib$IC50ADC_scale * (ic50 / ic50_ref)
  p$k_rep <- calib$k_rep; p$ED50_kill <- calib$ED50_kill
  p$Emax_CMP_kill <- calib$Emax_CMP; p$Emax_MPP_kill <- 1; p$Emax_MEP_kill <- calib$Emax_MEP
  p$Slope_CMP <- 1; p$Slope_MPP <- 1; p$Slope_MEP <- 1e-6
  p$k_depl_direct <- calib$k_depl_direct; p$k_depl_direct_MPP <- calib$k_depl_direct * 0.3
  p
}

# -- simule une dose (mg/kg), rend la trajectoire (jours) --
simulate <- function(p, dose, BW_kg, Tinfu_h = 1.5, tmax_day = 40) {
  st <- c(tdxd_hu_state0, init_state[!names(init_state) %in% c("C1", "C2", "Damage")])
  p$rate_fun <- make_tdxd_infusion(dose_mgkg = dose, BW_kg = BW_kg,
                  Tinfu_h = Tinfu_h, interval_h = 21 * 24, n_cycles = 1)
  o <- as.data.frame(lsoda(y = st, times = seq(0, tmax_day * 24, by = 6),
         func = pkpd_tdxd_fornari, parms = p,
         rtol = 1e-4, atol = 1e-6, maxsteps = 5e5, hmax = 1))
  o$td <- o$time / 24
  o
}

# -- grade CTCAE neutropenie a partir du fold nadir et de la baseline absolue --
ctcae_grade <- function(fold, baseline_abs) {
  x <- fold * baseline_abs
  if (x < 0.5) "G4" else if (x < 1) "G3" else if (x < 1.5) "G2" else if (x < 2) "G1" else "G0"
}
