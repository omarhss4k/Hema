############################################################
# plots_tdxd_nhp.R
# Graphiques publication-ready pour le poster
# Prérequis : run_pkpd_tdxd_nhp.R déjà sourcé
#             (sims_tdxd, init_pars, doses_tdxd disponibles)
############################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

if (!dir.exists("results_TDXD")) dir.create("results_TDXD")

# ── Palette Pierre Fabre ────────────────────────────────
dose_cols  <- c("3 mg/kg"  = "#2166ac",
                "10 mg/kg" = "#4dac26",
                "30 mg/kg" = "#d6604d")

dose_days  <- c(0, 21, 42)   # jours de dose Q3W

# ── Thème poster ────────────────────────────────────────
theme_poster <- theme_bw(base_size = 16) +
  theme(
    plot.title       = element_text(face = "bold", size = 17),
    axis.title       = element_text(size = 15),
    axis.text        = element_text(size = 13),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 13),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#e8f0f7"),
    strip.text       = element_text(face = "bold", size = 14)
  )

# ════════════════════════════════════════════════════════
# Mise en forme des simulations → long format
# ════════════════════════════════════════════════════════
sim_long <- bind_rows(lapply(seq_along(doses_tdxd), function(i) {
  sims_tdxd[[i]] %>%
    select(time_d, Neut, Mono, Ret, RBC, Plt, Damage) %>%
    mutate(Dose = paste0(doses_tdxd[i], " mg/kg"))
})) %>%
  mutate(Dose = factor(Dose, levels = paste0(doses_tdxd, " mg/kg")))

# Baselines
baselines <- c(
  Neut   = init_pars$Neut0,
  Mono   = init_pars$Mono0,
  Ret    = init_pars$Ret0,
  RBC    = init_pars$RBC0,
  Plt    = init_pars$Plt0
)

# ════════════════════════════════════════════════════════
# Figure principale : 4 panels (Neut, Plt, Ret, RBC)
# ════════════════════════════════════════════════════════
make_panel <- function(var, ylabel, baseline_val, ymin_frac = 0.7) {

  df <- sim_long %>% select(time_d, Dose, value = all_of(var))

  nadir_df <- df %>%
    group_by(Dose) %>%
    slice_min(value, n = 1) %>%
    ungroup()

  ggplot(df, aes(x = time_d, y = value, color = Dose)) +
    # Lignes de dose
    geom_vline(xintercept = dose_days, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    # Baseline
    geom_hline(yintercept = baseline_val, linetype = "dotted",
               color = "grey40", linewidth = 0.8) +
    # Courbes simulées
    geom_line(linewidth = 1.2) +
    # Points nadir
    geom_point(data = nadir_df, aes(x = time_d, y = value),
               shape = 25, size = 3, fill = "white", stroke = 1.5) +
    scale_color_manual(values = dose_cols) +
    scale_x_continuous(breaks = seq(0, 63, by = 21),
                       labels = paste0("J", seq(0, 63, by = 21))) +
    labs(x = "Temps (jours)", y = ylabel) +
    coord_cartesian(ylim = c(baseline_val * ymin_frac,
                             baseline_val * 1.05)) +
    theme_poster
}

p_neut <- make_panel("Neut", "Neutrophiles (10⁹/L)", baselines["Neut"], 0.75)
p_plt  <- make_panel("Plt",  "Plaquettes (10⁹/L)",   baselines["Plt"],  0.80)
p_ret  <- make_panel("Ret",  "Réticulocytes (10⁹/L)", baselines["Ret"], 0.55)
p_rbc  <- make_panel("RBC",  "GR (10⁹/L)",            baselines["RBC"], 0.95)

# ── Assemblage gridExtra ─────────────────────────────────
# Légende commune extraite d'un seul panel
get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}
legend_common <- get_legend(p_neut)

panels_noleg <- lapply(list(p_neut, p_plt, p_ret, p_rbc),
                       function(p) p + theme(legend.position = "none"))

fig_pd <- gridExtra::arrangeGrob(
  grobs  = c(panels_noleg, list(legend_common)),
  layout_matrix = rbind(c(1,2), c(3,4), c(5,5)),
  heights = c(4, 4, 0.8),
  top    = "Predicted Hematological Profiles — T-DXd Q3W × 3 cycles"
)

ggsave("results_TDXD/poster_PD_4panels.pdf",
       fig_pd, width = 14, height = 11, dpi = 300)
ggsave("results_TDXD/poster_PD_4panels.png",
       fig_pd, width = 14, height = 11, dpi = 300)
cat("  -> results_TDXD/poster_PD_4panels.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure Damage — pour montrer le driver PD
# ════════════════════════════════════════════════════════
p_damage <- ggplot(sim_long, aes(x = time_d, y = Damage, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.5) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = dose_cols) +
  scale_x_continuous(breaks = seq(0, 63, by = 21),
                     labels = paste0("J", seq(0, 63, by = 21))) +
  labs(title = "DNA Damage (model driver)",
       x = "Temps (jours)", y = "Damage (u.a.)") +
  theme_poster

ggsave("results_TDXD/poster_Damage.pdf",
       p_damage, width = 7, height = 5, dpi = 300)
ggsave("results_TDXD/poster_Damage.png",
       p_damage, width = 7, height = 5, dpi = 300)
cat("  -> results_TDXD/poster_Damage.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure nadir — barplot % changement vs baseline
# ════════════════════════════════════════════════════════
nadir_pct <- bind_rows(lapply(seq_along(doses_tdxd), function(i) {
  s <- sims_tdxd[[i]]
  data.frame(
    Dose     = paste0(doses_tdxd[i], " mg/kg"),
    Cellule  = c("Neut", "Plt", "Ret", "RBC"),
    Pct_chg  = c(
      (min(s$Neut) / init_pars$Neut0 - 1) * 100,
      (min(s$Plt)  / init_pars$Plt0  - 1) * 100,
      (min(s$Ret)  / init_pars$Ret0  - 1) * 100,
      (min(s$RBC)  / init_pars$RBC0  - 1) * 100
    )
  )
})) %>%
  mutate(Dose    = factor(Dose, levels = paste0(doses_tdxd, " mg/kg")),
         Cellule = factor(Cellule, levels = c("Neut", "Plt", "Ret", "RBC")))

p_nadir_bar <- ggplot(nadir_pct, aes(x = Cellule, y = Pct_chg, fill = Dose)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_hline(yintercept = -25, linetype = "dashed",
             color = "firebrick", linewidth = 0.7) +
  annotate("text", x = 4.5, y = -23, label = "Grade 2 (−25%)",
           color = "firebrick", size = 4, hjust = 1) +
  scale_fill_manual(values = dose_cols) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Nadir — % change from baseline",
       x = NULL, y = "% change vs baseline") +
  theme_poster

ggsave("results_TDXD/poster_nadir_barplot.pdf",
       p_nadir_bar, width = 8, height = 5, dpi = 300)
ggsave("results_TDXD/poster_nadir_barplot.png",
       p_nadir_bar, width = 8, height = 5, dpi = 300)
cat("  -> results_TDXD/poster_nadir_barplot.pdf / .png\n")

cat("\nTous les graphiques générés dans results_TDXD/\n")
