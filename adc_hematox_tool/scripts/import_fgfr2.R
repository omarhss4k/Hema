############################################################
# scripts/import_fgfr2.R
# Convertit le fichier hemato FGFR2 brut au format attendu par calibrate.R.
#
#   Source : data/nhp_hema_data.csv
#     colonnes : Animal_Id, Dose_mgkg, time_day,
#                Neut_1e3_uL, Mono_1e3_uL, Plt_1e3_uL, RBC_1e6_uL, Ret_1e9_L
#   Sortie : data/reference_profiles.csv
#     colonnes : Animal_Id, Dose_mgkg, time_day, Neut, Mono
#
# Unites : Neut_1e3_uL et Mono_1e3_uL sont en 10^3/uL = 10^9/L
#          -> identiques aux unites du modele, copie directe (aucune conversion).
#
# Lancer depuis la racine de l'outil :
#     Rscript scripts/import_fgfr2.R
############################################################

src <- "data/nhp_hema_data.csv"
if (!file.exists(src)) stop(sprintf("Fichier source absent : %s", src))

raw <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)

# -- selection tolerante des colonnes (accepte plusieurs orthographes) --
pick <- function(df, ...) {
  for (n in c(...)) if (n %in% names(df)) return(df[[n]])
  NULL
}

out <- data.frame(
  Animal_Id = pick(raw, "Animal_Id", "Animal_ID", "animal"),
  Dose_mgkg = pick(raw, "Dose_mgkg", "dose_mgkg", "Dose"),
  time_day  = pick(raw, "time_day", "jour", "day"),
  Neut      = pick(raw, "Neut_1e3_uL", "Neut"),
  Mono      = pick(raw, "Mono_1e3_uL", "Mono"),
  stringsAsFactors = FALSE
)

# garde les lignes ayant au moins une mesure
out <- out[!(is.na(out$Neut) & is.na(out$Mono)), ]

write.csv(out, "data/reference_profiles.csv", row.names = FALSE)
cat(sprintf("-> data/reference_profiles.csv  (%d lignes, doses : %s mg/kg)\n",
            nrow(out), paste(sort(unique(out$Dose_mgkg)), collapse = ", ")))

# baselines moyennes (1re valeur pre-dose de chaque animal) -> utiles pour la config
bl <- function(L) {
  v <- sapply(unique(out$Animal_Id), function(a) {
    s <- out[out$Animal_Id == a, ]; s <- s[order(s$time_day), ]
    w <- s[[L]][!is.na(s[[L]])]; if (length(w)) w[1] else NA
  })
  mean(v, na.rm = TRUE)
}
cat(sprintf("Baselines moyennes observees : Neut = %.2f | Mono = %.2f  (x10^9/L)\n",
            bl("Neut"), bl("Mono")))
cat("   -> reporte ces valeurs dans config/compounds.R (baseline_neut / baseline_mono).\n")
