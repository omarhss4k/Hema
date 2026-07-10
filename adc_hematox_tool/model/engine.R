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

# -- allometrie (Boxenbaum) : scale une PK d'une espece source vers cible --
#   CL, Q  ∝ BW^0.75      V1, V2 ∝ BW^1.0
# cp doit contenir pk_bw (BW ou la PK a ete mesuree, kg) ; BW_kg = cible.
# Alerte si le V scalE est physiquement implausible pour la cible.
allometric_scale <- function(cp) {
  if (is.null(cp$pk_bw) || isFALSE(cp$allometry)) return(cp)
  r <- cp$BW_kg / cp$pk_bw
  fCL <- r^0.75; fV <- r^1.0
  cat(sprintf("[allometrie] %s : BW %.3g -> %.3g kg (x%.0f) | CLx%.2f  Vx%.0f\n",
              cp$name, cp$pk_bw, cp$BW_kg, r, fCL, fV))
  cp$CL <- cp$CL * fCL; cp$Q <- cp$Q * fCL
  cp$V1 <- cp$V1 * fV; cp$V2 <- cp$V2 * fV
  vplasma <- 0.06 * cp$BW_kg   # ~60 mL/kg de volume sanguin
  if (cp$V1 > 3 * vplasma)
    cat(sprintf("[!] V1 scalE = %.3f L >> volume plasmatique attendu (~%.3f L) : verifier BW source et unites de concentration.\n",
                cp$V1, vplasma))
  cp
}

# -- construit le jeu de parametres pour un composE + une IC50 + la calibration --
# cp    : liste compound (PK, IC50, BW...) depuis config/compounds.R
# ic50  : IC50 myeloide du composE (nM)
# calib : liste CALIB (PD calibre + ancres de potency)
# ic50_ref : IC50 myeloide de la REFERENCE (ancre de translation)
build_pars <- function(cp, ic50, calib, ic50_ref) {
  cp <- allometric_scale(cp)          # scale PK si allometrie demandee
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

# -- lignes de seuils CTCAE (neutropenie) sur un graphe en FOLD --
# baseline_abs : baseline absolue (x10^9/L) ; ymax/xmax : cadre du plot.
# Trace les seuils G1..G4 convertis en fold, seulement ceux qui rentrent.
draw_grade_lines <- function(baseline_abs, ymax, xmax) {
  seuils <- c(G1 = 2.0, G2 = 1.5, G3 = 1.0, G4 = 0.5)   # x10^9/L
  couleurs <- c(G1 = "#f6c000", G2 = "#f08000", G3 = "#d84040", G4 = "#a01010")
  for (g in names(seuils)) {
    yf <- seuils[[g]] / baseline_abs
    if (yf <= ymax) {
      abline(h = yf, lty = 3, col = couleurs[[g]], lwd = 1)
      text(xmax, yf, g, col = couleurs[[g]], cex = 0.7, adj = c(1, -0.3))
    }
  }
}
