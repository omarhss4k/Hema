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

# -- lecture robuste : detecte le separateur (Excel FR = ';' + decimales ',') --
read_flex <- function(path) {
  l1 <- readLines(path, n = 1L, warn = FALSE)
  n <- function(ch) lengths(regmatches(l1, gregexpr(ch, l1, fixed = TRUE)))
  sep <- if (n(";") > n(",")) ";" else if (n("\t") > n(",")) "\t" else ","
  dec <- if (sep == ";") "," else "."   # FR : ';' -> decimales ','
  df <- read.csv(path, sep = sep, dec = dec, stringsAsFactors = FALSE,
                 check.names = FALSE, fill = TRUE, strip.white = TRUE)
  df <- df[, !grepl("^\\s*$|^X$|^X\\.", names(df)), drop = FALSE]  # vire colonnes vides
  cat(sprintf("[lecture] separateur='%s' decimale='%s' | colonnes : %s\n",
              ifelse(sep=="\t","TAB",sep), dec, paste(names(df), collapse = ", ")))
  df
}
raw <- read_flex(src)

# -- selection tolerante des colonnes (accepte plusieurs orthographes) --
pick <- function(df, ...) {
  for (n in c(...)) if (n %in% names(df)) return(df[[n]])
  NULL
}

# conversion numerique tolerante (gere une eventuelle decimale ',')
num <- function(x) as.numeric(gsub(",", ".", gsub("\\s", "", as.character(x))))

out <- data.frame(
  Animal_Id = pick(raw, "Animal_Id", "Animal_ID", "animal"),
  Dose_mgkg = num(pick(raw, "Dose_mgkg", "dose_mgkg", "Dose")),
  time_day  = num(pick(raw, "time_day", "jour", "day")),
  Neut      = num(pick(raw, "Neut_1e3_uL", "Neut")),
  Mono      = num(pick(raw, "Mono_1e3_uL", "Mono")),
  stringsAsFactors = FALSE
)
if (is.null(out$Neut) && is.null(out$Mono))
  stop("Colonnes Neut/Mono introuvables. Colonnes vues : ",
       paste(names(raw), collapse = ", "))

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
