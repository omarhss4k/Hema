############################################################
# plots_tdxd_nhp.R
# Graphiques publication-ready — ggplot2 uniquement
# Prérequis : run_pkpd_tdxd_nhp.R déjà sourcé
############################################################
library(ggplot2)
library(dplyr)
library(tidyr)

if (!dir.exists("results_TDXD")) dir.create("results_TDXD")

# ── Palette ─────────────────────────────────────────────
dose_cols <- c("3 mg/kg"  = "#2166ac",
               "10 mg/kg" = "#4dac26",
               "30 mg/kg" = "#d6604d")
dose_days <- c(0, 21, 42)

theme_poster <- theme_bw(base_size = 15) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    axis.title       = element_text(size = 14),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 13),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#dce8f5"),
    strip.text       = element_text(face = "bold", size = 13)
  )

# ── Long format ──────────────────────────────────────────
sim_long <- bind_rows(lapply(seq_along(doses_tdxd), function(i) {
  sims_tdxd[[i]] %>%
    select(time_d, Neut, Plt, Ret, RBC, Damage) %>%
    mutate(Dose = paste0(doses_tdxd[i], " mg/kg"))
})) %>%
  pivot_longer(cols = c(Neut, Plt, Ret, RBC, Damage),
               names_to = "Cellule", values_to = "Valeur") %>%
  mutate(
    Dose    = factor(Dose, levels = paste0(doses_tdxd, " mg/kg")),
    Cellule = factor(Cellule, levels = c("Neut", "Plt", "Ret", "RBC", "Damage"),
                     labels = c("Neutrophiles (10⁹/L)", "Plaquettes (10⁹/L)",
                                "Réticulocytes (10⁹/L)", "GR (10⁹/L)",
                                "Damage (u.a.)"))
  )

# Baselines par panel
base_df <- data.frame(
  Cellule = c("Neutrophiles (10⁹/L)", "Plaquettes (10⁹/L)",
              "Réticulocytes (10⁹/L)", "GR (10⁹/L)"),
  baseline = c(init_pars$Neut0, init_pars$Plt0,
               init_pars$Ret0,  init_pars$RBC0)
)

# Nadirs par panel + dose
nadir_df <- sim_long %>%
  filter(Cellule != "Damage (u.a.)") %>%
  group_by(Dose, Cellule) %>%
  slice_min(Valeur, n = 1) %>%
  ungroup()

# ════════════════════════════════════════════════════════
# Figure 1 — 4 panels cellulaires (facet_wrap)
# ════════════════════════════════════════════════════════
df_cells <- sim_long %>% filter(Cellule != "Damage (u.a.)")

p_cells <- ggplot(df_cells, aes(x = time_d, y = Valeur, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_hline(data = base_df, aes(yintercept = baseline),
             linetype = "dotted", color = "grey40",
             linewidth = 0.8, inherit.aes = FALSE) +
  geom_line(linewidth = 1.1) +
  geom_point(data = nadir_df, shape = 25, size = 3,
             fill = "white", stroke = 1.5) +
  scale_color_manual(values = dose_cols) +
  scale_x_continuous(breaks = seq(0, 63, by = 21),
                     labels = paste0("J", seq(0, 63, by = 21))) +
  facet_wrap(~Cellule, scales = "free_y", ncol = 2) +
  labs(title    = "Predicted Hematological Profiles — T-DXd Q3W × 3 cycles",
       subtitle = "▽ nadir  |  ··· baseline  |  --- dose day",
       x = "Temps (jours)", y = NULL) +
  theme_poster

ggsave("results_TDXD/poster_PD_4panels.pdf",
       p_cells, width = 13, height = 10, dpi = 300)
ggsave("results_TDXD/poster_PD_4panels.png",
       p_cells, width = 13, height = 10, dpi = 300)
cat("  -> results_TDXD/poster_PD_4panels.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure 2 — Damage seul
# ════════════════════════════════════════════════════════
p_damage <- ggplot(sim_long %>% filter(Cellule == "Damage (u.a.)"),
                   aes(x = time_d, y = Valeur, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = dose_cols) +
  scale_x_continuous(breaks = seq(0, 63, by = 21),
                     labels = paste0("J", seq(0, 63, by = 21))) +
  labs(title = "DNA Damage — model driver",
       x = "Temps (jours)", y = "Damage (u.a.)") +
  theme_poster

ggsave("results_TDXD/poster_Damage.pdf", p_damage, width = 7, height = 5, dpi = 300)
ggsave("results_TDXD/poster_Damage.png", p_damage, width = 7, height = 5, dpi = 300)
cat("  -> results_TDXD/poster_Damage.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure 3 — Barplot nadir % changement
# ════════════════════════════════════════════════════════
base_vals <- c("Neutrophiles (10⁹/L)" = init_pars$Neut0,
               "Plaquettes (10⁹/L)"   = init_pars$Plt0,
               "Réticulocytes (10⁹/L)"= init_pars$Ret0,
               "GR (10⁹/L)"           = init_pars$RBC0)

nadir_pct <- nadir_df %>%
  mutate(Pct = (Valeur / base_vals[as.character(Cellule)] - 1) * 100)

p_bar <- ggplot(nadir_pct, aes(x = Cellule, y = Pct, fill = Dose)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept =   0, linewidth = 0.5) +
  geom_hline(yintercept = -25, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text", x = 4.4, y = -23,
           label = "Grade 2 (-25%)", color = "firebrick",
           size = 4, hjust = 1) +
  scale_fill_manual(values = dose_cols) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Nadir — % change from baseline",
       x = NULL, y = "% change vs baseline") +
  theme_poster

ggsave("results_TDXD/poster_nadir_barplot.pdf", p_bar, width = 8, height = 5, dpi = 300)
ggsave("results_TDXD/poster_nadir_barplot.png", p_bar, width = 8, height = 5, dpi = 300)
cat("  -> results_TDXD/poster_nadir_barplot.pdf / .png\n")

cat("\nTous les graphiques generes dans results_TDXD/\n")
