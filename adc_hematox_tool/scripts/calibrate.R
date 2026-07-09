############################################################
# scripts/calibrate.R  --  Cale/valide le modele sur l'ADC de REFERENCE
# (PK + IC50 + profils observes) et SAUVEGARDE la calibration, qui sera
# lue automatiquement par scripts/predict.R.
#
# Lancer depuis la racine de l'outil :
#     Rscript scripts/calibrate.R
#
# Necessite les profils de reference dans data/ (voir README, confidentiel).
# Sorties : results/calibration_fit.pdf + results/calibration.rds
############################################################
source("config/compounds.R")
source("model/engine.R")
load_model("model")

ref <- REFERENCE
dir.create("results", showWarnings = FALSE)

## -- sauvegarde de la calibration (alimente predict.R) --
saveRDS(CALIB, "results/calibration.rds")
cat("-> results/calibration.rds  (lue par scripts/predict.R)\n")

## -- validation graphique si les profils observes sont disponibles --
## (jamais de quit() : sourcE dans une session RStudio, quit() fermerait
##  la session et ferait "disparaitre" le projet.)
if (!file.exists(ref$profiles)) {
  cat(sprintf("[!] %s absent : calibration sauvegardee, mais pas de figure de validation.\n", ref$profiles))
  cat("    (Depose tes profils de reference dans data/ pour visualiser le calage.)\n")
} else {
  obs <- read.csv(ref$profiles, stringsAsFactors = FALSE); obs$jr <- round(obs$time_day)
  bi <- function(aid, L) { s <- obs[obs$Animal_Id == aid, ]; s <- s[order(s$time_day), ]
    v <- s[[L]][!is.na(s[[L]])]; if (!length(v)) NA else v[1] }

  cols <- c("0.3" = "#2166ac", "1" = "#1a9641", "3" = "#b2182b")
  pdf("results/calibration_fit.pdf", width = 13, height = 5.4)
  par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3, 1))
  for (cc in list(c("Neut", init_pars$Neut0, "Neutrophiles"),
                  c("Mono", init_pars$Mono0, "Monocytes"))) {
    L <- cc[1]; b <- as.numeric(cc[2])
    plot(NA, xlim = c(0, 32), ylim = c(0, 3.5), xlab = "Jour",
         ylab = "fold vs baseline", main = paste(cc[3], "—", ref$name, "(calage)"))
    abline(h = 1, lty = 2, col = "grey60"); abline(h = 0.21, lty = 3, col = "red")
    for (d in ref$doses) {
      p <- build_pars(ref, ref$IC50_myelo, CALIB, ref$IC50_myelo)  # ref sur elle-meme -> ratio 1
      o <- simulate(p, d, ref$BW_kg, ref$Tinfu_h, tmax_day = 35)
      lines(o$td, o[[L]] / b, col = cols[as.character(d)], lwd = 2.5)
      od <- obs[obs$Dose_mgkg == d, ]
      for (k in 1:nrow(od)) { bb <- bi(od$Animal_Id[k], L); fo <- od[[L]][k] / bb
        term <- (od$jr[k] > 22 & fo < 0.05)
        points(od$jr[k], fo, col = cols[as.character(d)], pch = ifelse(term, 1, 19), cex = 1.2, lwd = 1.4) }
    }
    legend("topright", paste0(names(cols), " mg/kg"), col = cols, lwd = 2.5, pch = 19, bty = "n", cex = .8)
  }
  dev.off()
  cat("-> results/calibration_fit.pdf\n")
}
