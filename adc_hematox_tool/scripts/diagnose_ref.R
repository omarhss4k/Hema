############################################################
# scripts/diagnose_ref.R
# Diagnostic OBJECTIF des profils de reference : y a-t-il un vrai signal
# de depletion (nadir dose-dependant) ou pas ?
#
# Lit data/reference_profiles.csv (genere par import_fgfr2.R) et sort, par
# dose et par lignee, le nadir en fold vs baseline, le jour, le grade, et
# un verdict. Aucune donnee brute a transmettre : renvoie juste ce resume.
#
#     Rscript scripts/diagnose_ref.R
############################################################
f <- "data/reference_profiles.csv"
if (!file.exists(f)) stop("Lance d'abord scripts/import_fgfr2.R -> ", f)
d <- read.csv(f, stringsAsFactors = FALSE)

# baseline par animal = 1re valeur (jour le plus precoce, pre-dose inclus)
baseline <- function(aid, L) {
  s <- d[d$Animal_Id == aid, ]; s <- s[order(s$time_day), ]
  v <- s[[L]][!is.na(s[[L]])]; if (length(v)) v[1] else NA
}
d$foldN <- mapply(function(a, x) x / baseline(a, "Neut"), d$Animal_Id, d$Neut)
d$foldM <- mapply(function(a, x) x / baseline(a, "Mono"), d$Animal_Id, d$Mono)

cat("=================== DIAGNOSTIC REFERENCE ===================\n")
cat(sprintf("Fichier : %s | %d lignes | doses : %s mg/kg\n\n", f, nrow(d),
            paste(sort(unique(d$Dose_mgkg)), collapse = ", ")))

verdicts <- c()
for (L in c("N", "M")) {
  lab <- if (L == "N") "NEUTROPHILES" else "MONOCYTES"
  fcol <- if (L == "N") "foldN" else "foldM"
  cat(sprintf("---- %s (fold vs baseline) ----\n", lab))
  worst <- 1
  for (dose in sort(unique(d$Dose_mgkg))) {
    sub <- d[d$Dose_mgkg == dose & is.finite(d[[fcol]]) & d$time_day > 0, ]
    if (!nrow(sub)) next
    i <- which.min(sub[[fcol]]); nadir <- sub[[fcol]][i]; jour <- sub$time_day[i]
    worst <- min(worst, nadir)
    cat(sprintf("  %5s mg/kg : nadir %.2f x baseline  @ jour %s  (min sur %d pts)\n",
                dose, nadir, jour, nrow(sub)))
  }
  # verdict : un vrai signal myelotoxique descend nettement sous 1 (< ~0.5 = G3+)
  v <- if (worst >= 0.8) "PAS de signal (reste >= 0.8x : non-myelotoxique)"
       else if (worst >= 0.5) "signal FAIBLE (nadir 0.5-0.8x)"
       else "signal NET de depletion (nadir < 0.5x)"
  cat(sprintf("  => nadir le plus bas toutes doses = %.2f  -> %s\n\n", worst, v))
  verdicts <- c(verdicts, sprintf("%s: %s", lab, v))
}
cat("=================== VERDICT ===================\n")
for (v in verdicts) cat(" -", v, "\n")
cat("\nSi 'PAS de signal' : FGFR2 confirme non-myelotoxique -> mauvais ancrage\n",
    "de calibration de depletion (rien a caler). Voir recommandation.\n", sep = "")
