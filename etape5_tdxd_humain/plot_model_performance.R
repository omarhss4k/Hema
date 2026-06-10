############################################################
# plot_model_performance.R
# Évaluation de la performance prédictive du modèle PK/PD
# IC95% Wilson, RMSE/MAE, graphique de calibration
#
# Prérequis : run_pkpd_tdxd_human_population.R
#             (→ results/population_results.rds)
############################################################

library(ggplot2)
library(ggrepel)

# -- Chargement résultats modèle ---------------------------
rds_path <- "results/population_results.rds"

if (!file.exists(rds_path)) {
  stop(paste0(
    "Fichier manquant : ", rds_path, "\n",
    "Lancez d'abord run_pkpd_tdxd_human_population.R"
  ))
}

results <- readRDS(rds_path)
n       <- nrow(results)
cat(sprintf("\n═══ Performance prédictive -- N = %d patients ═══\n\n", n))

# -- Données FDA de référence (BLA 761139, DESTINY-Breast01, n=184) -
fda_neut   <- c(G0=71, G1=7,  G2=7,  G3=13, G4=3)
fda_anemia <- c(G0=30, G1=37, G2=24, G3=8,  G4=1)
fda_plt    <- c(G0=63, G1=30, G2=4,  G3=2,  G4=1)

fda_tout_grade <- c(
  Neut   = 100 - fda_neut["G0"],
  Anemia = 100 - fda_anemia["G0"],
  Plt    = 100 - fda_plt["G0"]
)
fda_g34 <- c(
  Neut   = fda_neut["G3"]   + fda_neut["G4"],
  Anemia = fda_anemia["G3"] + fda_anemia["G4"],
  Plt    = fda_plt["G3"]    + fda_plt["G4"]
)

# -- Calcul des proportions modèle ------------------------
grade_order <- c("G0", "G1", "G2", "G3", "G4")

tab_neut   <- table(factor(results$Grade_Neut,   levels = grade_order))
tab_anemia <- table(factor(results$Grade_Anemia, levels = grade_order))
tab_plt    <- table(factor(results$Grade_Plt,    levels = grade_order))

pct_n <- setNames(as.numeric(100 * tab_neut   / n), grade_order)
pct_a <- setNames(as.numeric(100 * tab_anemia / n), grade_order)
pct_p <- setNames(as.numeric(100 * tab_plt    / n), grade_order)

mod_tout_grade <- c(
  Neut   = sum(pct_n[c("G1","G2","G3","G4")]),
  Anemia = sum(pct_a[c("G1","G2","G3","G4")]),
  Plt    = sum(pct_p[c("G1","G2","G3","G4")])
)
mod_g34 <- c(
  Neut   = sum(pct_n[c("G3","G4")]),
  Anemia = sum(pct_a[c("G3","G4")]),
  Plt    = sum(pct_p[c("G3","G4")])
)

# -- IC95% Wilson -----------------------------------------
# Formule : p ± 1.96 * sqrt(p*(1-p)/n)  avec p en proportion [0,1]
wilson_ci <- function(pct_vec, n) {
  p    <- pct_vec / 100
  marg <- 1.96 * sqrt(p * (1 - p) / n) * 100
  data.frame(
    p   = pct_vec,
    lo  = pmax(0, pct_vec - marg),
    hi  = pmin(100, pct_vec + marg)
  )
}

ci_neut_tg  <- wilson_ci(mod_tout_grade["Neut"],   n)
ci_neut_g34 <- wilson_ci(mod_g34["Neut"],          n)
ci_an_tg    <- wilson_ci(mod_tout_grade["Anemia"], n)
ci_an_g34   <- wilson_ci(mod_g34["Anemia"],        n)
ci_plt_tg   <- wilson_ci(mod_tout_grade["Plt"],    n)
ci_plt_g34  <- wilson_ci(mod_g34["Plt"],           n)

# -- Tableau IC95% -----------------------------------------
ic95_df <- data.frame(
  Toxicite  = c("Neutropénie",      "Neutropénie",
                "Anémie",           "Anémie",
                "Thrombocytopénie", "Thrombocytopénie"),
  Metrique  = rep(c("Tout grade", "G3-4"), 3),
  Pred_pct  = c(mod_tout_grade["Neut"],   mod_g34["Neut"],
                mod_tout_grade["Anemia"], mod_g34["Anemia"],
                mod_tout_grade["Plt"],    mod_g34["Plt"]),
  IC95_lo   = c(ci_neut_tg$lo,  ci_neut_g34$lo,
                ci_an_tg$lo,    ci_an_g34$lo,
                ci_plt_tg$lo,   ci_plt_g34$lo),
  IC95_hi   = c(ci_neut_tg$hi,  ci_neut_g34$hi,
                ci_an_tg$hi,    ci_an_g34$hi,
                ci_plt_tg$hi,   ci_plt_g34$hi),
  FDA_pct   = c(fda_tout_grade["Neut"],   fda_g34["Neut"],
                fda_tout_grade["Anemia"], fda_g34["Anemia"],
                fda_tout_grade["Plt"],    fda_g34["Plt"]),
  stringsAsFactors = FALSE,
  row.names = NULL
)

cat("IC95% Wilson (méthode binomiale approchée, N=300) :\n")
cat("-------------------------------------------------------------\n")
fmt <- "  %-20s %-12s : %5.1f%%  [%5.1f%% - %5.1f%%]  (FDA: %5.1f%%)\n"
for (i in seq_len(nrow(ic95_df))) {
  cat(sprintf(fmt,
              ic95_df$Toxicite[i],
              ic95_df$Metrique[i],
              ic95_df$Pred_pct[i],
              ic95_df$IC95_lo[i],
              ic95_df$IC95_hi[i],
              ic95_df$FDA_pct[i]))
}

# Sauvegarde CSV
dir.create("results", showWarnings = FALSE)
write.csv(ic95_df, "results/grade_IC95.csv", row.names = FALSE)
cat("\n  -> results/grade_IC95.csv\n")

# -- Métriques de performance RMSE / MAE ------------------
pred_vec <- ic95_df$Pred_pct
obs_vec  <- ic95_df$FDA_pct

rmse_val <- sqrt(mean((pred_vec - obs_vec)^2))
mae_val  <- mean(abs(pred_vec - obs_vec))

cat(sprintf("\nMétriques de performance (6 comparaisons prédites vs FDA) :\n"))
cat(sprintf("  RMSE = %.2f points de pourcentage\n", rmse_val))
cat(sprintf("  MAE  = %.2f points de pourcentage\n", mae_val))

# -- Graphique de calibration ------------------------------
cal_df <- ic95_df
cal_df$label <- paste0(
  ifelse(cal_df$Toxicite == "Neutropénie",      "Neut",
  ifelse(cal_df$Toxicite == "Anémie",           "Anémie",
                                                 "Thrombo")),
  " ", cal_df$Metrique
)

tox_colors <- c(
  "Neutropénie"      = "#2c7bb6",
  "Anémie"           = "#d7191c",
  "Thrombocytopénie" = "#1a9641"
)
metric_shapes <- c("Tout grade" = 16, "G3-4" = 17)  # rond = tout grade, triangle = G3-4

# Annotation RMSE/MAE (coin supérieur gauche)
ann_text <- sprintf("RMSE = %.2f pp\nMAE  = %.2f pp", rmse_val, mae_val)

p_cal <- ggplot(cal_df, aes(x = FDA_pct, y = Pred_pct,
                             color = Toxicite, shape = Metrique)) +

  # Diagonale y = x
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey60", linewidth = 0.8) +

  # IC95% verticaux (sur axe Prédictions)
  geom_errorbar(aes(ymin = IC95_lo, ymax = IC95_hi),
                width = 0.8, linewidth = 0.7, alpha = 0.6) +

  # Points
  geom_point(size = 4, stroke = 0.8) +

  # Labels
  geom_text_repel(aes(label = label),
                  size = 3.2, fontface = "italic",
                  box.padding = 0.4, point.padding = 0.3,
                  max.overlaps = 20, show.legend = FALSE) +

  # Annotation RMSE/MAE
  annotate("text",
           x = -Inf, y = Inf,
           label = ann_text,
           hjust = -0.1, vjust = 1.3,
           size = 3.5, family = "mono",
           color = "grey30") +

  scale_color_manual(values = tox_colors, name = "Toxicité") +
  scale_shape_manual(values = metric_shapes, name = "Métrique") +

  scale_x_continuous(limits = c(0, 80),
                     labels = function(x) paste0(x, "%")) +
  scale_y_continuous(limits = c(0, 80),
                     labels = function(x) paste0(x, "%")) +

  labs(
    x        = "Observé FDA (DESTINY-Breast01, n=184)",
    y        = "Prédit modèle (N=300) ± IC95% Wilson",
    title    = "Calibration : Modèle PK/PD vs FDA",
    subtitle = sprintf(
      "6 comparaisons (3 toxicités × 2 métriques) -- RMSE=%.2f pp, MAE=%.2f pp",
      rmse_val, mae_val)
  ) +

  theme_classic(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    panel.grid.major = element_line(color = "grey92"),
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.title       = element_text(face = "bold")
  )

# -- Export -----------------------------------------------
dir.create("results_PKPD_human", showWarnings = FALSE)

out_pdf <- "results_PKPD_human/calibration_plot.pdf"
out_png <- "results_PKPD_human/calibration_plot.png"

ggsave(out_pdf, plot = p_cal, width = 7, height = 6.5, device = cairo_pdf)
ggsave(out_png, plot = p_cal, width = 7, height = 6.5, dpi = 300)

cat(sprintf("\n  -> %s\n", out_pdf))
cat(sprintf("  -> %s\n\n", out_png))
