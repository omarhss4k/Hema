# =============================================================================
# Courbe de croissance tumorale — Fc-silent FGFR2-huBPA-LP1
# Structure Excel :
#   Bloc 1 : moyennes (Group | j0 | j4 | ...)
#   <ligne vide>
#   Bloc 2 : SEM     (Group | j0 | j4 | ...)
# =============================================================================

library(readxl)
library(ggplot2)

# =============================================================================
# 1. LECTURE BRUTE (toutes lignes, sans en-tête automatique)
# =============================================================================

raw_all <- suppressMessages(
  read_xlsx("TumorVolume_FGFR2.xlsx",
            col_names    = FALSE,
            .name_repair = "minimal")
)

# Repérer les lignes "Group" (en-têtes de chaque bloc)
header_idx <- which(raw_all[[1]] == "Group")
cat("Lignes 'Group' trouvées :", header_idx, "\n")
stopifnot("Deux blocs attendus (means + SEM)" = length(header_idx) == 2)

# Extraire les temps (en jours) depuis le premier en-tête
times <- suppressWarnings(as.numeric(as.character(unlist(raw_all[header_idx[1], -1]))))
valid <- !is.na(times)
times <- times[valid]
cat("Temps (jours) :", times, "\n")

# =============================================================================
# 2. EXTRACTION DES DEUX BLOCS
# =============================================================================

extract_block <- function(start_row, end_row) {
  blk <- raw_all[start_row:end_row, ]
  blk <- blk[!is.na(blk[[1]]) & blk[[1]] != "", ]   # suppr. lignes vides
  # Garder seulement les colonnes avec des temps valides
  dat <- blk[, c(TRUE, valid)]
  storage.mode(dat[[1]]) <- "character"
  for (j in 2:ncol(dat)) storage.mode(dat[[j]]) <- "numeric"
  dat
}

# Bloc moyennes : lignes après header_idx[1] jusqu'avant header_idx[2]
blk_mean <- extract_block(header_idx[1] + 1, header_idx[2] - 1)
# Bloc SEM : lignes après header_idx[2] jusqu'à la fin
blk_sem  <- extract_block(header_idx[2] + 1, nrow(raw_all))

cat("\nGroupes (moyennes) :", blk_mean[[1]], "\n")
cat("Groupes (SEM)      :", blk_sem[[1]],  "\n")

# =============================================================================
# 3. ASSEMBLAGE PAR GROUPE
# =============================================================================

patterns <- c("Group 01", "Group 03", "Group 04")
labels   <- c("Contrôle", "10 mg/kg", "3 mg/kg")

get_row <- function(blk, pattern)
  as.numeric(unlist(blk[grep(pattern, blk[[1]], ignore.case = TRUE)[1], -1]))

df_list <- lapply(seq_along(patterns), function(i) {
  mn  <- get_row(blk_mean, patterns[i])
  sem <- get_row(blk_sem,  patterns[i])
  ok  <- !is.na(mn)
  data.frame(
    jour   = times[ok],
    mean   = mn[ok],
    sem    = sem[ok],
    Groupe = labels[i]
  )
})

df_plot <- do.call(rbind, df_list)
df_plot$Groupe <- factor(df_plot$Groupe,
                         levels = c("Contrôle", "3 mg/kg", "10 mg/kg"))

cat("\nAperçu des données :\n")
print(head(df_plot, 10))

# =============================================================================
# 4. GRAPHIQUE
# =============================================================================

dose_days <- c(0, 14, 28, 42)    # Q2W × 4 — adapter si protocole différent
ymax  <- max(df_plot$mean + df_plot$sem, na.rm = TRUE)
cols  <- c("Contrôle" = "#888888", "3 mg/kg" = "#4393C3", "10 mg/kg" = "#2166AC")

p <- ggplot(df_plot, aes(x = jour, y = mean, color = Groupe, group = Groupe)) +
  # Lignes verticales des doses
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  # Courbes + points + barres SEM
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.8, linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  # Triangles de traitement sous l'axe
  annotate("point", x = dose_days, y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text",  x = max(dose_days) + 1.5, y = -ymax * 0.06,
           label = "= Traitement", hjust = 0, size = 3.2, color = "#CC0000") +
  # Échelles
  scale_x_continuous(breaks = seq(0, max(df_plot$jour), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "Fc-silent FGFR2-huBPA-LP1 — Courbe de croissance tumorale",
    subtitle = "mean ± SEM",
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm³)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 10, color = "grey50"),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

print(p)
ggsave("scripts/plot_courbe_tumorale_FGFR2.png", p,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_courbe_tumorale_FGFR2.png\n")

# =============================================================================
# 5. SAUVEGARDE DES SEM (pour pondération PKPD)
# =============================================================================

sem_ctrl <- get_row(blk_sem, "Group 01")
sem_d10  <- get_row(blk_sem, "Group 03")
sem_d3   <- get_row(blk_sem, "Group 04")

dat_sem <- list(
  time   = times,
  ctrl   = sem_ctrl,
  d3     = sem_d3,
  d10    = sem_d10
)

save(dat_sem, file = "scripts/sem_tumorale_FGFR2.RData")
cat("SEM sauvegardées → scripts/sem_tumorale_FGFR2.RData\n")
