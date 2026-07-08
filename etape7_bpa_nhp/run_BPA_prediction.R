############################################################
# run_BPA_prediction.R
# PREDICTION de la myelotoxicite du BPA (ADC) a partir de sa
# PK + IC50, SANS donnees in vivo BPA.
#
# Principe (translation) :
#  - le PD est CALIBRE sur l'ADC cytotoxique de reference
#    (Fornari IC50 + deplation directe, cf. run_ADCcyto_G4.R) ;
#  - on branche la PK propre du BPA ;
#  - la POTENCY est translatee via le rapport des IC50 in vitro :
#       IC50_ADC_ugmL(BPA) = IC50_ADC_ugmL(ref) x (IC50_BPA / IC50_ref)
#
# Sort neutrophiles + monocytes (lignees myeloides) et grades CTCAE.
#
# ⚠️ CAVEAT : si le BPA agit AUSSI par voie immunitaire (STING),
#    cette prediction cytotoxique est un PLANCHER (la vraie
#    toxicite pourrait etre superieure).
#
# Lancer depuis etape7_bpa_nhp/ :  Rscript run_BPA_prediction.R
#
# CONFIDENTIEL -- usage interne.
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("../etape5_tdxd_humain/parameters_tdxd_human.R")
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")   # make_tdxd_infusion
  source("../etape3_tdxd_rat/pkpd_tdxd_rat.R")          # pkpd_tdxd_fornari
}))

## ═══════ [1] PK du BPA (2-cmt, depuis ton Shiny) ═══════
BW_KG   <- 2.5
CL_BPA  <- 0.000467305265769336   # L/h
V1_BPA  <- 0.0586062942310418     # L
V2_BPA  <- 0.0424186702071202     # L
Q_BPA   <- 0.00102316504406689    # L/h
TINFU_H <- 1.5
DOSES   <- c(0.3, 1, 3)           # mg/kg

## ═══════ [2] POTENCY : reference + IC50 du BPA ═══════
# Reference = ADC cytotoxique calibre
IC50_REF      <- 0.008615   # nM  (IC50 myelo ADC de reference)
IC50ADC_REF   <- 27.70      # ug/mL (echelle de potency du modele calibre)
# IC50 myeloide du BPA -- A TRANCHER (2 valeurs candidates) :
IC50_BPA_LIST <- c(0.015, 7.8)   # nM

## ═══════ [3] PD CALIBRE (fige, = run_ADCcyto_G4.R) ═══════
K_REP <- 0.002; ED50_KILL <- 0.1; EMAX_CMP <- 0.982
K_DEPL_DIR <- 0.25; EMAX_MEP <- 1e-4

mkp <- function(ic50_bpa) {
  p <- c(init_pars, tdxd_pars_hu)
  p$CL_ADC <- CL_BPA; p$V1_ADC <- V1_BPA; p$Q_ADC <- Q_BPA; p$V2_ADC <- V2_BPA
  p$k_int <- 0; p$krel_power <- 0; p$krel_factor <- 1
  p$use_ADC_driver <- TRUE; p$Damage_threshold <- 0
  # translation de la potency par le rapport d'IC50 in vitro :
  p$IC50_ADC_ugmL <- IC50ADC_REF * (ic50_bpa / IC50_REF)
  p$k_rep <- K_REP; p$ED50_kill <- ED50_KILL
  p$Emax_CMP_kill <- EMAX_CMP; p$Emax_MPP_kill <- 1; p$Emax_MEP_kill <- EMAX_MEP
  p$Slope_CMP <- 1; p$Slope_MPP <- 1; p$Slope_MEP <- 1e-6
  p$k_depl_direct <- K_DEPL_DIR; p$k_depl_direct_MPP <- K_DEPL_DIR * 0.3
  p
}

## ═══════ [4] SIMULATION ═══════
sim <- function(ic50_bpa, dose) {
  p  <- mkp(ic50_bpa)
  st <- c(tdxd_hu_state0, init_state[!names(init_state) %in% c("C1", "C2", "Damage")])
  p$rate_fun <- make_tdxd_infusion(dose_mgkg = dose, BW_kg = BW_KG,
                  Tinfu_h = TINFU_H, interval_h = 21 * 24, n_cycles = 1)
  o <- as.data.frame(lsoda(y = st, times = seq(0, 40 * 24, by = 6),
         func = pkpd_tdxd_fornari, parms = p,
         rtol = 1e-4, atol = 1e-6, maxsteps = 5e5, hmax = 1))
  o$td <- o$time / 24
  o
}

## ═══════ [5] FIGURE (fourchette IC50) ═══════
cols <- c("0.3" = "#2166ac", "1" = "#1a9641", "3" = "#b2182b")
dir.create("results", showWarnings = FALSE)
pdf("results/BPA_prediction_bracket.pdf", width = 13, height = 5.6)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3.5, 1))
for (cc in list(c("Neut", init_pars$Neut0, "Neutrophiles"),
                c("Mono", init_pars$Mono0, "Monocytes"))) {
  L <- cc[1]; b <- as.numeric(cc[2])
  plot(NA, xlim = c(0, 40), ylim = c(0, 1.15), xlab = "Jour",
       ylab = "fold vs baseline", main = paste("BPA predit —", cc[3]))
  abline(h = 1, lty = 2, col = "grey60")
  abline(h = 0.21, lty = 3, col = "red"); text(39, 0.27, "G4", col = "red", cex = .7, adj = 1)
  for (d in DOSES) {
    lines(sim(IC50_BPA_LIST[1], d)$td, sim(IC50_BPA_LIST[1], d)[[L]] / b, col = cols[as.character(d)], lwd = 2.6, lty = 1)
    lines(sim(IC50_BPA_LIST[2], d)$td, sim(IC50_BPA_LIST[2], d)[[L]] / b, col = cols[as.character(d)], lwd = 2, lty = 3)
  }
  legend("bottomleft",
    c(paste0(names(cols), " mg/kg"), "",
      paste0("IC50=", IC50_BPA_LIST[1], " (puissant)"),
      paste0("IC50=", IC50_BPA_LIST[2], " (faible)")),
    col = c(cols, NA, "grey30", "grey30"),
    lwd = c(2.6, 2.6, 2.6, NA, 2.6, 2), lty = c(1, 1, 1, NA, 1, 3), bty = "n", cex = .78)
}
dev.off()

## ═══════ [6] RESUME CONSOLE (grades CTCAE) ═══════
cn <- function(x) if (x < 0.5) "G4" else if (x < 1) "G3" else if (x < 1.5) "G2" else if (x < 2) "G1" else "G0"
BASE_NEUT <- 2.4; BASE_MONO <- 0.8
for (ic in IC50_BPA_LIST) {
  cat(sprintf("\n### BPA IC50_myelo = %s nM  (IC50_ADC = %.0f ug/mL) ###\n",
      ic, IC50ADC_REF * (ic / IC50_REF)))
  for (d in DOSES) {
    o <- sim(ic, d)
    nN <- min(o$Neut) / init_pars$Neut0; nM <- min(o$Mono) / init_pars$Mono0
    cat(sprintf(" dose %-4s mg/kg : Neut %.2f@j%.0f (%s) | Mono %.2f@j%.0f (%s)\n",
      d, nN, o$td[which.min(o$Neut)], cn(nN * BASE_NEUT),
      nM, o$td[which.min(o$Mono)], cn(nM * BASE_MONO)))
  }
}
cat("\n-> results/BPA_prediction_bracket.pdf\n")
cat("⚠️ Prediction CYTOTOXIQUE = plancher si le BPA agit aussi par voie STING.\n")
