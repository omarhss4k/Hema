############################################################
# 67_plots_pd_poster_nhp.R
# Poster figure — PD predicted profiles ONLY (no observed data)
# FGFR2 inhibitor NHP  |  Q3W × 3 cycles + recovery
#
# 3 figures :
#   A. poster_PD_pred_4panels.pdf/png  — 4 cell types, facet × dose
#   B. poster_PD_pred_alldoses.pdf/png — 4 cell types, lines per dose
#   C. poster_PD_pred_nadir.pdf/png    — nadir barplot % vs baseline
#
# Prérequis : 68_run_fgfr2_nhp.R sourcé (objets sims_fgfr2, doses_fgfr2, init_pars)
############################################################

library(ggplot2)
library(dplyr)
library(tidyr)

if (!exists("sims_fgfr2") || !exists("doses_fgfr2")) {
  cat("Lancement de 68_run_fgfr2_nhp.R...\n")
  source("68_run_fgfr2_nhp.R")
}

if (!dir.exists("results")) dir.create("results")

# ── Palette & thème ──────────────────────────────────────
dose_cols <- c("4 mg/kg"  = "#2166ac",
               "13 mg/kg" = "#4dac26",
               "26 mg/kg" = "#f4a582",
               "39 mg/kg" = "#d6604d")

dose_days <- c(0, 21, 42)

theme_poster_pred <- theme_classic(base_size = 15) +
  theme(
    plot.title        = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle     = element_text(size = 12, hjust = 0.5, color = "grey40"),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12),
    legend.title      = element_text(size = 13, face = "bold"),
    legend.text       = element_text(size = 12),
    legend.position   = "bottom",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    strip.background  = element_rect(fill = "#2c3e50", color = NA),
    strip.text        = element_text(color = "white", face = "bold", size = 13),
    plot.margin       = margin(10, 15, 10, 10)
  )

# ── Long format simulé ───────────────────────────────────
sim_long <- bind_rows(lapply(seq_along(doses_fgfr2), function(i) {
  sims_fgfr2[[i]] %>%
    select(time_d, Neut, Plt, Ret, RBC) %>%
    mutate(Dose = factor(paste0(doses_fgfr2[i], " mg/kg"),
                         levels = paste0(doses_fgfr2, " mg/kg")))
})) %>%
  pivot_longer(cols = c(Neut, Plt, Ret, RBC),
               names_to  = "Cell",
               values_to = "Value") %>%
  mutate(Cell = factor(Cell,
                       levels = c("Neut", "Plt", "Ret", "RBC"),
                       labels = c("Neutrophils (10⁹/L)", "Platelets (10⁹/L)",
                                  "Reticulocytes (10⁹/L)", "RBC (10¹²/L)")))

# Baselines
base_df <- data.frame(
  Cell     = factor(c("Neutrophils (10⁹/L)", "Platelets (10⁹/L)",
                      "Reticulocytes (10⁹/L)", "RBC (10¹²/L)"),
                    levels = c("Neutrophils (10⁹/L)", "Platelets (10⁹/L)",
                               "Reticulocytes (10⁹/L)", "RBC (10¹²/L)")),
  baseline = c(init_pars$Neut0, init_pars$Plt0,
               init_pars$Ret0,  init_pars$RBC0)
)

# Nadirs
nadir_df <- sim_long %>%
  group_by(Dose, Cell) %>%
  slice_min(Value, n = 1, with_ties = FALSE) %>%
  ungroup()

xbreaks <- seq(0, 120, by = 21)
xlabels <- paste0("D", xbreaks)

# ════════════════════════════════════════════════════════
# Figure A — 4 panels × 4 doses (facet Cell, color Dose)
# ════════════════════════════════════════════════════════
pA <- ggplot(sim_long, aes(x = time_d, y = Value, color = Dose)) +

  # Dose days
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.45) +

  # Baseline par panel
  geom_hline(data = base_df, aes(yintercept = baseline),
             linetype = "dotted", color = "grey35",
             linewidth = 0.8, inherit.aes = FALSE) +

  # Courbes prédites
  geom_line(linewidth = 1.15) +

  # Marque nadir (triangle bas)
  geom_point(data = nadir_df,
             aes(x = time_d, y = Value, color = Dose),
             shape = 25, size = 3.5, fill = "white",
             stroke = 1.8, inherit.aes = FALSE) +

  scale_color_manual(values = dose_cols, name = "Dose") +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +

  facet_wrap(~ Cell, scales = "free_y", ncol = 2) +

  labs(
    title    = "Predicted Hematological Profiles — FGFR2 Inhibitor NHP",
    subtitle = "Q3W × 3 cycles + recovery  |  ▽ predicted nadir  |  ··· species baseline  |  --- dose day",
    x = "Time (days)", y = NULL
  ) +
  theme_poster_pred +
  guides(color = guide_legend(nrow = 1,
                              override.aes = list(linewidth = 2.5)))

ggsave("results/poster_PD_pred_4panels.pdf",
       pA, width = 14, height = 10, device = cairo_pdf)
ggsave("results/poster_PD_pred_4panels.png",
       pA, width = 14, height = 10, dpi = 300)
cat("  -> results/poster_PD_pred_4panels.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure B — toutes doses superposées, facet cell type
# Vue "overlay" : montre clairement la dose-réponse
# ════════════════════════════════════════════════════════

# Normalisation % baseline pour comparaison visuelle
sim_pct <- sim_long %>%
  left_join(base_df, by = "Cell") %>%
  mutate(Pct = Value / baseline * 100)

nadir_pct <- sim_pct %>%
  group_by(Dose, Cell) %>%
  slice_min(Pct, n = 1, with_ties = FALSE) %>%
  ungroup()

pB <- ggplot(sim_pct, aes(x = time_d, y = Pct, color = Dose)) +

  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey75", linewidth = 0.4) +

  geom_hline(yintercept = 100, linetype = "dotted",
             color = "grey35", linewidth = 0.8) +

  # Zone nadir grade 2 (−25% = 75%)
  geom_hline(yintercept = 75, linetype = "longdash",
             color = "#e74c3c", linewidth = 0.7, alpha = 0.7) +

  geom_line(linewidth = 1.2) +

  geom_point(data = nadir_pct,
             aes(x = time_d, y = Pct, color = Dose),
             shape = 25, size = 3.5, fill = "white",
             stroke = 1.8, inherit.aes = FALSE) +

  # Annotation grade 2
  annotate("text", x = 5, y = 73,
           label = "–25% (Grade 2 threshold)", color = "#e74c3c",
           size = 3.5, hjust = 0, fontface = "italic") +

  scale_color_manual(values = dose_cols, name = "Dose") +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.05, 0.05))) +

  facet_wrap(~ Cell, ncol = 2) +

  labs(
    title    = "Dose-Dependent Hematological Suppression — FGFR2 Inhibitor NHP",
    subtitle = "% of baseline  |  ▽ predicted nadir  |  ··· baseline (100%)  |  — Grade 2 threshold (−25%)",
    x = "Time (days)", y = "% of baseline"
  ) +
  theme_poster_pred +
  guides(color = guide_legend(nrow = 1,
                              override.aes = list(linewidth = 2.5)))

ggsave("results/poster_PD_pred_alldoses.pdf",
       pB, width = 14, height = 10, device = cairo_pdf)
ggsave("results/poster_PD_pred_alldoses.png",
       pB, width = 14, height = 10, dpi = 300)
cat("  -> results/poster_PD_pred_alldoses.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure C — Barplot nadir % vs baseline (4 doses × 4 lignées)
# ════════════════════════════════════════════════════════
nadir_bar <- nadir_pct %>%
  mutate(PctChange = Pct - 100)

# Seuils de grade CTCAE v5 (% vs baseline approximatif)
pC <- ggplot(nadir_bar,
             aes(x = Dose, y = PctChange, fill = Dose)) +

  geom_col(width = 0.65, color = "white", linewidth = 0.3) +

  geom_hline(yintercept = -25, linetype = "dashed",
             color = "#e67e22", linewidth = 0.9) +
  geom_hline(yintercept = -50, linetype = "longdash",
             color = "#e74c3c", linewidth = 0.9) +

  annotate("text", x = Inf, y = -23,
           label = "Grade 2 (−25%)", color = "#e67e22",
           size = 3.5, hjust = 1.05, fontface = "italic") +
  annotate("text", x = Inf, y = -48,
           label = "Grade 3 (−50%)", color = "#e74c3c",
           size = 3.5, hjust = 1.05, fontface = "italic") +

  # Labels valeurs
  geom_text(aes(label = sprintf("%.0f%%", PctChange),
                vjust = ifelse(PctChange < 0, 1.35, -0.4)),
            size = 3.8, fontface = "bold", color = "white") +

  scale_fill_manual(values = dose_cols, name = "Dose") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.12, 0.08))) +

  facet_wrap(~ Cell, ncol = 2) +

  labs(
    title    = "Predicted Hematological Nadir — % Change from Baseline",
    subtitle = "FGFR2 inhibitor NHP  |  Q3W × 3 cycles  |  CTCAE v5 severity thresholds",
    x = NULL, y = "% change vs baseline"
  ) +
  theme_poster_pred +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11)) +
  guides(fill = guide_legend(nrow = 1))

ggsave("results/poster_PD_pred_nadir.pdf",
       pC, width = 13, height = 10, device = cairo_pdf)
ggsave("results/poster_PD_pred_nadir.png",
       pC, width = 13, height = 10, dpi = 300)
cat("  -> results/poster_PD_pred_nadir.pdf / .png\n")

# ── Console récap ─────────────────────────────────────────
cat("\n═══ Nadir summary (predicted, % vs baseline) ═══\n")
nadir_bar %>%
  select(Dose, Cell, Nadir_pct = Pct, Change_pct = PctChange) %>%
  arrange(Cell, Dose) %>%
  { cat(format(., digits = 3)); cat("\n") }
