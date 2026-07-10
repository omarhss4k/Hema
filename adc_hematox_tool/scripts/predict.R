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
  base_abs <- if (L == "Neut") cp$baseline_neut else cp$baseline_mono
  draw_grade_lines(base_abs, ymax = 1.15, xmax = 40)
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

## -- RESUME CONSOLE + EXPORT CSV (grades CTCAE) --
cat(sprintf("\n===== PREDICTION : %s =====\n", cp$name))
tab <- data.frame()   # accumulateur pour l'export CSV
for (ic in cp$IC50_myelo) {
  cat(sprintf("\n### IC50_myelo = %s nM  (IC50_ADC = %.0f ug/mL) ###\n",
      ic, calib$IC50ADC_scale * (ic / IC50_REF)))
  for (d in cp$doses) {
    p <- build_pars(cp, ic, calib, IC50_REF)
    o <- simulate(p, d, cp$BW_kg, cp$Tinfu_h)
    nN <- min(o$Neut) / init_pars$Neut0; jN <- o$td[which.min(o$Neut)]
    nM <- min(o$Mono) / init_pars$Mono0; jM <- o$td[which.min(o$Mono)]
    gN <- ctcae_grade(nN, cp$baseline_neut); gM <- ctcae_grade(nM, cp$baseline_mono)
    cat(sprintf(" dose %-4s mg/kg : Neut %.2f@j%.0f (%s) | Mono %.2f@j%.0f (%s)\n",
      d, nN, jN, gN, nM, jM, gM))
    tab <- rbind(tab, data.frame(
      compose = cp$name, IC50_nM = ic, dose_mgkg = d,
      neut_nadir_fold = round(nN, 3), neut_nadir_jour = round(jN),
      neut_absolu = round(nN * cp$baseline_neut, 3), neut_grade = gN,
      mono_nadir_fold = round(nM, 3), mono_nadir_jour = round(jM),
      mono_absolu = round(nM * cp$baseline_mono, 3), mono_grade = gM))
  }
}
csv_out <- sprintf("results/prediction_%s_grades.csv", cp$name)
write.csv(tab, csv_out, row.names = FALSE)
cat(sprintf("\n-> results/prediction_%s.pdf\n", cp$name))
cat(sprintf("-> %s\n", csv_out))
if (!is.null(cp$note)) cat(sprintf("NOTE : %s\n", cp$note))
