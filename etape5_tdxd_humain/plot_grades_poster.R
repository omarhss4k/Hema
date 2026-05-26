############################################################
# plot_grades_poster.R
# Figure poster : grades CTCAE — T-DXd 5.4 mg/kg Q3W × 6
# Neutropénie | Anémie | Thrombocytopénie
# Modèle semi-mécaniste (N=300) vs FDA DESTINY-Breast01
#
# Ajoute lignes de référence FDA :
#   --- FDA tout grade (G≥1)
#   ··· FDA G3-4
#
# Prérequis : run_pkpd_tdxd_human_population.R (→ results/population_results.rds)
############################################################

library(ggplot2)

# ── Thème poster ─────────────────────────────────────────
theme_poster_grades <- theme_classic(base_size = 15) +
  theme(
    strip.background   = element_rect(fill = "#2c3e50", color = NA),
    strip.text         = element_text(color = "white", face = "bold", size = 15),
    axis.text          = element_text(size = 13),
    axis.title         = element_text(size = 14, face = "bold"),
    legend.title       = element_text(size = 13, face = "bold"),
    legend.text        = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey88"),
    plot.title         = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(size = 12, hjust = 0.5, color = "grey40"),
    legend.position    = "bottom"
  )

# ── Définitions grades ────────────────────────────────────
grade_order     <- c("G0", "G1", "G2", "G3", "G4")
grade_order_rev <- rev(grade_order)   # G4 au bas du stacked bar

grade_cols <- c(
  G0 = "#4575b4",   # bleu foncé  — pas de toxicité
  G1 = "#91bfdb",   # bleu clair
  G2 = "#fee090",   # jaune
  G3 = "#fc8d59",   # orange
  G4 = "#d73027"    # rouge
)

grade_labels <- c(
  G0 = "G0 (aucune)",
  G1 = "G1 (légère)",
  G2 = "G2 (modérée)",
  G3 = "G3 (sévère)",
  G4 = "G4 (critique)"
)

# ── Données FDA (BLA 761139) ──────────────────────────────
# DESTINY-Breast01, n=184 (U201)
fda_neut   <- c(G0=71, G1=7,  G2=7,  G3=13, G4=3)
fda_anemia <- c(G0=30, G1=37, G2=24, G3=8,  G4=1)
# FDA pooled 5.4 mg/kg (N=234, Table 44 BLA 761139)
fda_plt    <- c(G0=63, G1=30, G2=4,  G3=2,  G4=1)

tox_levels <- c("Neutropénie", "Anémie", "Thrombocytopénie")

# ── Chargement résultats modèle ───────────────────────────
rds_path <- "results/population_results.rds"

if (!file.exists(rds_path)) {
  stop(paste0(
    "Fichier manquant : ", rds_path, "\n",
    "Lancez d'abord run_pkpd_tdxd_human_population.R"
  ))
}

results <- readRDS(rds_path)
n_ok    <- nrow(results)
cat(sprintf("Résultats chargés : N = %d patients\n", n_ok))

tab_neut   <- table(factor(results$Grade_Neut,   levels = grade_order))
tab_anemia <- table(factor(results$Grade_Anemia, levels = grade_order))
tab_plt    <- table(factor(results$Grade_Plt,    levels = grade_order))

pct_n <- setNames(as.numeric(round(100 * tab_neut   / n_ok, 1)), grade_order)
pct_a <- setNames(as.numeric(round(100 * tab_anemia / n_ok, 1)), grade_order)
pct_p <- setNames(as.numeric(round(100 * tab_plt    / n_ok, 1)), grade_order)

# ── Assemblage long format ────────────────────────────────
make_long <- function(pct_mod, pct_fda, tox_name) {
  rbind(
    data.frame(source = "Modèle", grade = grade_order,
               pct = as.numeric(pct_mod),         tox = tox_name,
               stringsAsFactors = FALSE),
    data.frame(source = "FDA",    grade = grade_order,
               pct = as.numeric(pct_fda),         tox = tox_name,
               stringsAsFactors = FALSE)
  )
}

df <- rbind(
  make_long(pct_n, fda_neut,   "Neutropénie"),
  make_long(pct_a, fda_anemia, "Anémie"),
  make_long(pct_p, fda_plt,    "Thrombocytopénie")
)

df$grade  <- factor(df$grade,  levels = grade_order_rev)  # G4 bas → G0 haut
df$source <- factor(df$source, levels = c("Modèle", "FDA"))
df$tox    <- factor(df$tox,    levels = tox_levels)

# Calcul positions ymax/ymid pour labels (ordre G4 → G0 dans chaque groupe)
df <- df[order(df$tox, df$source, df$grade), ]
df$ymax <- ave(df$pct, paste(df$tox, df$source), FUN = cumsum)
df$ymid <- df$ymax - df$pct / 2

df_label <- df[df$pct >= 6, ]

# ── Références FDA par facette ────────────────────────────
# Dans le stacked bar (G4 bas, G0 haut) :
#   tout grade  = G1+G2+G3+G4 = 100 - G0  (hauteur depuis le bas jusqu'au bas de G0)
#   G3-4        = G4+G3                    (hauteur des 2 premières couches depuis le bas)
fda_refs <- data.frame(
  tox       = factor(tox_levels, levels = tox_levels),
  toutgrade = c(100 - fda_neut["G0"],
                100 - fda_anemia["G0"],
                100 - fda_plt["G0"]),
  g34       = c(fda_neut["G4"]   + fda_neut["G3"],
                fda_anemia["G4"] + fda_anemia["G3"],
                fda_plt["G4"]    + fda_plt["G3"]),
  row.names = NULL
)

# ── Annotations texte pour les lignes de référence ────────
# Positionnées à x = 2.5 (juste à droite de la barre FDA, x=2)
# avec coord_cartesian(xlim = c(0.4, 3.2)) pour l'espace
ann_tg <- data.frame(
  tox   = fda_refs$tox,
  y     = fda_refs$toutgrade,
  label = paste0("FDA tout grade: ", round(fda_refs$toutgrade), "%"),
  stringsAsFactors = FALSE
)
ann_g34 <- data.frame(
  tox   = fda_refs$tox,
  y     = fda_refs$g34,
  label = paste0("FDA G3-4: ", round(fda_refs$g34), "%"),
  stringsAsFactors = FALSE
)
ann_tg$tox  <- factor(ann_tg$tox,  levels = tox_levels)
ann_g34$tox <- factor(ann_g34$tox, levels = tox_levels)

# ── Figure ────────────────────────────────────────────────
p <- ggplot(df, aes(x = source, y = pct, fill = grade)) +

  # Barres empilées (Modèle | FDA)
  geom_col(width = 0.65, color = "white", linewidth = 0.4) +

  # Labels % sur les segments significatifs (≥6%)
  geom_text(data = df_label,
            aes(x = source, y = ymid,
                label = sprintf("%.0f%%", pct)),
            size = 3.8, fontface = "bold", color = "white",
            inherit.aes = FALSE) +

  # Ligne de référence FDA tout grade (tiretée, rouge)
  geom_hline(data = fda_refs,
             aes(yintercept = toutgrade),
             linetype = "dashed", color = "#e74c3c",
             linewidth = 1.0, inherit.aes = FALSE) +

  # Ligne de référence FDA G3-4 (pointillée, rouge foncé)
  geom_hline(data = fda_refs,
             aes(yintercept = g34),
             linetype = "dotted", color = "#9b2335",
             linewidth = 1.2, inherit.aes = FALSE) +

  # Label pour ligne tout grade
  geom_text(data = ann_tg,
            aes(y = y, label = label),
            x = 2.5, hjust = 0, vjust = -0.45,
            color = "#e74c3c", size = 3.5, fontface = "italic",
            inherit.aes = FALSE) +

  # Label pour ligne G3-4
  geom_text(data = ann_g34,
            aes(y = y, label = label),
            x = 2.5, hjust = 0, vjust = -0.45,
            color = "#9b2335", size = 3.5, fontface = "italic",
            inherit.aes = FALSE) +

  # Palette grades (G4 → G0, de bas en haut)
  scale_fill_manual(
    values = grade_cols,
    breaks = grade_order,           # G0→G4 dans la légende
    labels = grade_labels,
    name   = "Grade CTCAE v5",
    guide  = guide_legend(reverse = TRUE, nrow = 1,
                          override.aes = list(color = NA))
  ) +

  scale_y_continuous(
    limits = c(0, 106),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +

  # Espace à droite pour les annotations texte
  coord_cartesian(xlim = c(0.45, 3.2), clip = "off") +

  facet_wrap(~tox, nrow = 1) +

  labs(
    x        = NULL,
    y        = "Patients (%)",
    title    = "Hématotoxicité T-DXd 5.4 mg/kg Q3W × 6 cycles — Modèle vs FDA",
    subtitle = sprintf(
      "Modèle semi-mécaniste (N = %d) vs DESTINY-Breast01 (n = 184, BLA 761139)",
      n_ok)
  ) +

  theme_poster_grades

# ── Export ───────────────────────────────────────────────
dir.create("results_PKPD_human", showWarnings = FALSE)

out_pdf <- "results_PKPD_human/poster_grades_tdxd.pdf"
out_png <- "results_PKPD_human/poster_grades_tdxd.png"

ggsave(out_pdf, plot = p, width = 14, height = 6.5, device = cairo_pdf)
ggsave(out_png, plot = p, width = 14, height = 6.5, dpi = 300)

cat(sprintf("  -> %s\n", out_pdf))
cat(sprintf("  -> %s\n", out_png))

# ── Résumé comparatif console ─────────────────────────────
cat("\n═══ Synthèse Modèle vs FDA ═══\n")
fmt <- "  %-18s : tout grade %5.1f%%  |  G3-4 %5.1f%%\n"
cat(sprintf(fmt, "Neut  — Modèle",
            sum(pct_n[c("G1","G2","G3","G4")]), sum(pct_n[c("G3","G4")])))
cat(sprintf(fmt, "Neut  — FDA",
            100 - fda_neut["G0"], fda_neut["G3"] + fda_neut["G4"]))
cat(sprintf(fmt, "Anémie — Modèle",
            sum(pct_a[c("G1","G2","G3","G4")]), sum(pct_a[c("G3","G4")])))
cat(sprintf(fmt, "Anémie — FDA",
            100 - fda_anemia["G0"], fda_anemia["G3"] + fda_anemia["G4"]))
cat(sprintf(fmt, "Thrombo — Modèle",
            sum(pct_p[c("G1","G2","G3","G4")]), sum(pct_p[c("G3","G4")])))
cat(sprintf(fmt, "Thrombo — FDA",
            100 - fda_plt["G0"], fda_plt["G3"] + fda_plt["G4"]))
cat("\n")

# Lancer analyses de performance
source("plot_model_performance.R")
source("sensitivity_analysis.R")
