############################################################
# scripts/run_all.R  --  Enchaine le pipeline complet :
#   1) calibrate.R  (cale sur l'ADC de reference -> calibration.rds)
#   2) predict.R    (predit le nouveau composE -> PDF + CSV + grades)
#
# Lancer depuis la racine de l'outil :
#     Rscript scripts/run_all.R
#
# Chaque etape s'execute dans son propre processus R (isole et propre).
# Si les profils de reference sont absents, l'etape 1 sauvegarde quand
# meme la calibration par defaut et l'etape 2 continue.
############################################################
rscript <- file.path(R.home("bin"), "Rscript")

cat("========== [1/2] CALIBRATION ==========\n")
s1 <- system2(rscript, "scripts/calibrate.R")
if (s1 != 0) cat("[!] calibrate.R a renvoye un code non nul -- on continue avec la calibration disponible.\n")

cat("\n========== [2/2] PREDICTION ==========\n")
s2 <- system2(rscript, "scripts/predict.R")
if (s2 != 0) stop("predict.R a echoue.")

cat("\n========== TERMINE ==========\n")
cat("Sorties dans results/ : calibration_fit.pdf, prediction_*.pdf, prediction_*_grades.csv\n")
