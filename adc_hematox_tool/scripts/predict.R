############################################################
# scripts/predict.R  --  PREDIT la myelotoxicite d'un nouveau composE
# a partir de sa PK + IC50, sans donnees in vivo.
#
# Lancer depuis la racine de l'outil :
#     Rscript scripts/predict.R
#
# Sorties : results/prediction_<composE>.pdf + grades CTCAE console.
############################################################
source("config/compounds.R")
source("model/engine.R")
load_model("model")

## -- calibration : lue depuis results/ si presente, sinon config --
CAL_FILE <- "results/calibration.rds"
if (file.exists(CAL_FILE)) {
  calib <- readRDS(CAL_FILE)
  cat(sprintf("[calibration lue depuis %s]\n", CAL_FILE))
} else {
  calib <- CALIB
  cat("[calibration par defaut depuis config/compounds.R -- lance scripts/calibrate.R pour re-caler]\n")
}
IC50_REF <- REFERENCE$IC50_myelo

cp <- NEW_COMPOUND
dir.create("results", showWarnings = FALSE)

## -- FIGURE (une courbe par dose x IC50) --
cols <- c("#2166ac", "#1a9641", "#b2182b", "#762a83", "#e08214")
ltys <- 1:4   # une texture par valeur d'IC50
pdf(sprintf("results/prediction_%s.pdf", cp$name), width = 13, height = 5.6)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3.5, 1))
for (cc in list(c("Neut", init_pars$Neut0, "Neutrophiles"),
                c("Mono", init_pars$Mono0, "Monocytes"))) {
  L <- cc[1]; b <- as.numeric(cc[2])
  plot(NA, xlim = c(0, 40), ylim = c(0, 1.15), xlab = "Jour",
       ylab = "fold vs baseline", main = paste(cp$name, "predit —", cc[3]))
  abline(h = 1, lty = 2, col = "grey60")
  abline(h = 0.21, lty = 3, col = "red"); text(39, 0.27, "G4", col = "red", cex = .7, adj = 1)
  for (di in seq_along(cp$doses)) {
    for (ii in seq_along(cp$IC50_myelo)) {
      p <- build_pars(cp, cp$IC50_myelo[ii], calib, IC50_REF)
      o <- simulate(p, cp$doses[di], cp$BW_kg, cp$Tinfu_h)
      lines(o$td, o[[L]] / b, col = cols[di], lwd = 2.4, lty = ltys[ii])
    }
  }
  legend("bottomleft",
    c(paste0(cp$doses, " mg/kg"), "", paste0("IC50=", cp$IC50_myelo, " nM")),
    col = c(cols[seq_along(cp$doses)], NA, rep("grey30", length(cp$IC50_myelo))),
    lwd = 2, lty = c(rep(1, length(cp$doses)), NA, ltys[seq_along(cp$IC50_myelo)]),
    bty = "n", cex = .75)
}
dev.off()

## -- RESUME CONSOLE (grades CTCAE) --
cat(sprintf("\n===== PREDICTION : %s =====\n", cp$name))
for (ic in cp$IC50_myelo) {
  cat(sprintf("\n### IC50_myelo = %s nM  (IC50_ADC = %.0f ug/mL) ###\n",
      ic, calib$IC50ADC_scale * (ic / IC50_REF)))
  for (d in cp$doses) {
    p <- build_pars(cp, ic, calib, IC50_REF)
    o <- simulate(p, d, cp$BW_kg, cp$Tinfu_h)
    nN <- min(o$Neut) / init_pars$Neut0; nM <- min(o$Mono) / init_pars$Mono0
    cat(sprintf(" dose %-4s mg/kg : Neut %.2f@j%.0f (%s) | Mono %.2f@j%.0f (%s)\n",
      d, nN, o$td[which.min(o$Neut)], ctcae_grade(nN, cp$baseline_neut),
      nM, o$td[which.min(o$Mono)], ctcae_grade(nM, cp$baseline_mono)))
  }
}
cat(sprintf("\n-> results/prediction_%s.pdf\n", cp$name))
if (!is.null(cp$note)) cat(sprintf("NOTE : %s\n", cp$note))
