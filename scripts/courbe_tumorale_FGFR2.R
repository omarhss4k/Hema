# =============================================================================
# Courbe de croissance tumorale — Fc-silent FGFR2-huBPA-LP1
# Lit la structure complète du fichier Excel :
#   - données individuelles  → mean ± SEM calculés
#   - données déjà moyennées → tracé sans barres d'erreur
# =============================================================================

library(readxl)
library(ggplot2)

# =============================================================================
# 1. LECTURE COMPLÈTE DU FICHIER
# =============================================================================

raw_full <- read_xlsx("TumorVolume_FGFR2.xlsx",
                      col_types  = "guess",
                      .name_repair = "minimal")

# Supprimer colonnes entièrement vides
raw_full <- raw_full[, colSums(!is.na(raw_full)) > 0]

cat("Structure du fichier :", nrow(raw_full), "lignes ×", ncol(raw_full), "colonnes\n")
cat("Première colonne :\n")
print(raw_full[[1]])

# Extraire les temps (en jours) depuis les noms de colonnes (col 2+)
time_d_all <- suppressWarnings(as.numeric(colnames(raw_full)[-1]))
valid_t    <- !is.na(time_d_all)
time_d     <- time_d_all[valid_t]
raw_data   <- raw_full[, c(TRUE, valid_t)]

cat("\nTemps (jours) :", time_d, "\n")

# =============================================================================
# 2. EXTRACTION PAR GROUPE
# =============================================================================

labels <- raw_data[[1]]

get_mat <- function(pattern) {
  idx <- grep(pattern, labels, ignore.case = TRUE)
  if (length(idx) == 0) stop(paste("Groupe introuvable :", pattern))
  m <- as.matrix(raw_data[idx, -1])
  storage.mode(m) <- "numeric"
  m
}

mat_ctrl <- get_mat("Group 01")
mat_d3   <- get_mat("Group 04")
mat_d10  <- get_mat("Group 03")

cat("\nLignes par groupe — Contrôle:", nrow(mat_ctrl),
    "/ 3 mg/kg:", nrow(mat_d3),
    "/ 10 mg/kg:", nrow(mat_d10), "\n")

has_individuals <- max(nrow(mat_ctrl), nrow(mat_d3), nrow(mat_d10)) > 1
if (has_individuals) cat("→ Données individuelles détectées : SEM calculée\n") else
  cat("→ Données moyennées uniquement : SEM non disponible\n")

# =============================================================================
# 3. CALCUL MEAN ± SEM
# =============================================================================

summarise_mat <- function(mat, groupe) {
  mn  <- colMeans(mat, na.rm = TRUE)
  sem <- if (nrow(mat) > 1)
    apply(mat, 2, function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
  else
    rep(NA_real_, ncol(mat))
  data.frame(jour = time_d, mean = mn, sem = sem, Groupe = groupe)
}

df_plot <- rbind(
  summarise_mat(mat_ctrl, "Contrôle"),
  summarise_mat(mat_d3,   "3 mg/kg"),
  summarise_mat(mat_d10,  "10 mg/kg")
)
df_plot <- df_plot[!is.na(df_plot$mean), ]
df_plot$Groupe <- factor(df_plot$Groupe,
                         levels = c("Contrôle", "3 mg/kg", "10 mg/kg"))

# =============================================================================
# 4. GRAPHIQUE
# =============================================================================

dose_days <- c(0, 14, 28, 42)   # Q2W × 4 — ajuster si protocole différent
ymax      <- max(df_plot$mean + ifelse(is.na(df_plot$sem), 0, df_plot$sem),
                 na.rm = TRUE)
cols <- c("Contrôle" = "#888888", "3 mg/kg" = "#4393C3", "10 mg/kg" = "#2166AC")

p <- ggplot(df_plot, aes(x = jour, y = mean, color = Groupe, group = Groupe)) +
  # Lignes de doses
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  # Courbes
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  # Barres d'erreur (si SEM disponible)
  { if (has_individuals)
      geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                    width = 0.8, linewidth = 0.5) } +
  # Triangles de traitement (sous l'axe)
  annotate("point", x = dose_days, y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text", x = max(dose_days) + 1, y = -ymax * 0.06,
           label = "= Traitement", hjust = 0, size = 3, color = "#CC0000") +
  # Axes
  scale_x_continuous(breaks = seq(0, max(df_plot$jour), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "Fc-silent FGFR2-huBPA-LP1 — Courbe de croissance tumorale",
    subtitle = if (has_individuals) "mean ± SEM" else "Données moyennées (SEM non disponible)",
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm³)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "right",
    plot.subtitle    = element_text(size = 10, color = "grey50"),
    plot.margin      = margin(t = 5, r = 10, b = 20, l = 5)
  )

print(p)
ggsave("scripts/plot_courbe_tumorale_FGFR2.png", p, width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_courbe_tumorale_FGFR2.png\n")
