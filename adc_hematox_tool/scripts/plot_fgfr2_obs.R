############################################################
# scripts/plot_fgfr2_obs.R
# Graphes des profils hemato OBSERVES (FGFR2) -- toutes les lignees
# presentes dans le fichier (Neut, Mono, Plt, RBC, Ret).
#
# Source : data/nhp_hema_data.csv (Excel FR ';' + decimales ',' accepte)
# Sorties : results/fgfr2_obs_fold.pdf   (fold vs baseline)
#           results/fgfr2_obs_absolu.pdf (valeurs absolues)
#
#     Rscript scripts/plot_fgfr2_obs.R
############################################################
src <- "data/nhp_hema_data.csv"
if (!file.exists(src)) stop("Fichier absent : ", src, " (depose ton hemato FGFR2)")
dir.create("results", showWarnings = FALSE)

# -- lecture robuste (separateur/decimale auto) --
read_flex <- function(path) {
  l1 <- readLines(path, n = 1L, warn = FALSE)
  n <- function(ch) lengths(regmatches(l1, gregexpr(ch, l1, fixed = TRUE)))
  sep <- if (n(";") > n(",")) ";" else if (n("\t") > n(",")) "\t" else ","
  dec <- if (sep == ";") "," else "."
  df <- read.csv(path, sep = sep, dec = dec, stringsAsFactors = FALSE,
                 check.names = FALSE, fill = TRUE, strip.white = TRUE)
  df[, !grepl("^\\s*$|^X$|^X\\.", names(df)), drop = FALSE]
}
raw <- read_flex(src)
num <- function(x) as.numeric(gsub(",", ".", gsub("\\s", "", as.character(x))))
pick <- function(...) { for (n in c(...)) if (n %in% names(raw)) return(raw[[n]]); NULL }

d <- data.frame(
  Animal = pick("Animal_Id", "Animal_ID", "animal"),
  dose   = num(pick("Dose_mgkg", "dose_mgkg", "Dose")),
  jour   = num(pick("time_day", "jour", "day")),
  stringsAsFactors = FALSE
)
# lignees disponibles (nom affiche -> colonnes possibles)
LIN <- list(
  Neutrophiles = c("Neut_1e3_uL", "Neut"),
  Monocytes    = c("Mono_1e3_uL", "Mono"),
  Plaquettes   = c("Plt_1e3_uL", "Plt"),
  "Glob. rouges" = c("RBC_1e6_uL", "RBC"),
  Reticulocytes = c("Ret_1e9_L", "Ret"))
for (L in names(LIN)) { v <- pick(LIN[[L]][1], LIN[[L]][2]); if (!is.null(v)) d[[L]] <- num(v) }
present <- names(LIN)[names(LIN) %in% names(d)]
cat("Lignees tracees :", paste(present, collapse = ", "), "\n")
cat("Doses :", paste(sort(unique(d$dose)), collapse = ", "), "mg/kg\n")

doses <- sort(unique(d$dose))
pal <- grDevices::colorRampPalette(c("#2166ac", "#1a9641", "#f0a000", "#b2182b"))
cols <- setNames(pal(length(doses)), as.character(doses))

baseline <- function(aid, L) { s <- d[d$Animal == aid, ]; s <- s[order(s$jour), ]
  v <- s[[L]][!is.na(s[[L]])]; if (length(v)) v[1] else NA }

# grille de panneaux
np <- length(present); nc <- min(3, np); nr <- ceiling(np / nc)

draw <- function(mode) {         # mode = "fold" ou "abs"
  fn <- sprintf("results/fgfr2_obs_%s.pdf", ifelse(mode == "fold", "fold", "absolu"))
  pdf(fn, width = 4.6 * nc, height = 4.0 * nr)
  par(mfrow = c(nr, nc), mar = c(4, 4.3, 2.6, 1))
  for (L in present) {
    yv <- d[[L]]
    if (mode == "fold") yv <- mapply(function(a, x) x / baseline(a, L), d$Animal, d[[L]])
    ylab <- if (mode == "fold") "fold vs baseline" else "valeur (unites source)"
    ylim <- if (mode == "fold") c(0, max(2, max(yv, na.rm = TRUE) * 1.05)) else c(0, max(yv, na.rm = TRUE) * 1.05)
    plot(NA, xlim = range(d$jour, na.rm = TRUE), ylim = ylim,
         xlab = "Jour", ylab = ylab, main = L)
    if (mode == "fold") abline(h = 1, lty = 2, col = "grey60")
    for (a in unique(d$Animal)) {
      s <- d[d$Animal == a, ]; s <- s[order(s$jour), ]
      yy <- if (mode == "fold") s[[L]] / baseline(a, L) else s[[L]]
      cc <- cols[as.character(s$dose[1])]
      lines(s$jour, yy, col = cc, lwd = 1.6); points(s$jour, yy, col = cc, pch = 19, cex = 0.9)
    }
    if (L == present[1])
      legend("topright", paste0(doses, " mg/kg"), col = cols, lwd = 2, pch = 19, bty = "n", cex = .85)
  }
  dev.off(); cat("->", fn, "\n")
}
draw("fold"); draw("abs")
