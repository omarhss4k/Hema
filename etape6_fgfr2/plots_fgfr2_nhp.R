############################################################
# plots_fgfr2_nhp.R
# Graphiques publication-ready — ggplot2 uniquement
# Prérequis : run_fgfr2_nhp.R déjà sourcé
############################################################
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Auto-source si objets manquants ─────────────────────
if (!exists("sims_fgfr2") || !exists("doses_fgfr2")) {
  cat("Lancement de run_fgfr2_nhp.R...\n")
  source("run_fgfr2_nhp.R")
}

if (!dir.exists("results")) dir.create("results")

# ════════════════════════════════════════════════════════
# DONNÉES OBSERVÉES — à remplir manuellement
# Colonnes obligatoires : dose_mgkg, jour, Neut, Plt, Ret, RBC
# Ret = réticulocytes ABSOLUS (10⁹/L)  [= %RET/100 × RBC]
# Mettre NA si non mesuré
# ════════════════════════════════════════════════════════
obs_data <- data.frame(
  dose_mgkg = c(),   # ex: c(3, 3, 10, 10, 30, 30)
  jour      = c(),   # ex: c(2, 8,  2,  8,  2,  8)
  Neut      = c(),
  Plt       = c(),
  Ret       = c(),
  RBC       = c()
)

# ── Palette ─────────────────────────────────────────────
dose_cols <- c("4 mg/kg"  = "#2166ac",
               "13 mg/kg" = "#4dac26",
               "26 mg/kg" = "#f4a582",
               "39 mg/kg" = "#d6604d")
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
sim_long <- bind_rows(lapply(seq_along(doses_fgfr2), function(i) {
  sims_fgfr2[[i]] %>%
    select(time_d, Neut, Plt, Ret, RBC, Damage) %>%
    mutate(Dose = paste0(doses_fgfr2[i], " mg/kg"))
})) %>%
  pivot_longer(cols = c(Neut, Plt, Ret, RBC, Damage),
               names_to = "Cellule", values_to = "Valeur") %>%
  mutate(
    Dose    = factor(Dose, levels = paste0(doses_fgfr2, " mg/kg")),
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

# Observed data en long format (si rempli)
obs_long <- if (nrow(obs_data) > 0) {
  obs_data %>%
    mutate(Dose = paste0(dose_mgkg, " mg/kg")) %>%
    pivot_longer(cols = c(Neut, Plt, Ret, RBC),
                 names_to = "Cellule", values_to = "Valeur") %>%
    mutate(
      Dose    = factor(Dose, levels = paste0(doses_fgfr2, " mg/kg")),
      Cellule = factor(Cellule,
                       levels = c("Neut", "Plt", "Ret", "RBC"),
                       labels = c("Neutrophiles (10⁹/L)", "Plaquettes (10⁹/L)",
                                  "Réticulocytes (10⁹/L)", "GR (10⁹/L)"))
    ) %>%
    filter(!is.na(Valeur))
} else NULL

p_cells <- ggplot(df_cells, aes(x = time_d, y = Valeur, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_hline(data = base_df, aes(yintercept = baseline),
             linetype = "dotted", color = "grey40",
             linewidth = 0.8, inherit.aes = FALSE) +
  geom_line(linewidth = 1.1) +
  geom_point(data = nadir_df, shape = 25, size = 3,
             fill = "white", stroke = 1.5) +
  { if (!is.null(obs_long))
      geom_point(data = obs_long, aes(x = jour, y = Valeur, color = Dose),
                 shape = 16, size = 2.5, inherit.aes = FALSE)
  } +
  scale_color_manual(values = dose_cols) +
  scale_x_continuous(breaks = seq(0, 120, by = 21),
                     labels = paste0("J", seq(0, 120, by = 21))) +
  facet_wrap(~Cellule, scales = "free_y", ncol = 2) +
  labs(title    = "Predicted Hematological Profiles — FGFR2 inhibitor Q3W × 3 cycles (+ recovery)",
       subtitle = "▽ nadir  |  ··· baseline  |  --- dose day",
       x = "Temps (jours)", y = NULL) +
  theme_poster

ggsave("results/poster_PD_4panels.pdf",
       p_cells, width = 13, height = 10, dpi = 300)
ggsave("results/poster_PD_4panels.png",
       p_cells, width = 13, height = 10, dpi = 300)
cat("  -> results/poster_PD_4panels.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure 2 — Damage seul
# ════════════════════════════════════════════════════════
p_damage <- ggplot(sim_long %>% filter(Cellule == "Damage (u.a.)"),
                   aes(x = time_d, y = Valeur, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = dose_cols) +
  scale_x_continuous(breaks = seq(0, 120, by = 21),
                     labels = paste0("J", seq(0, 120, by = 21))) +
  labs(title = "DNA Damage — model driver",
       x = "Temps (jours)", y = "Damage (u.a.)") +
  theme_poster

ggsave("results/poster_Damage.pdf", p_damage, width = 7, height = 5, dpi = 300)
ggsave("results/poster_Damage.png", p_damage, width = 7, height = 5, dpi = 300)
cat("  -> results/poster_Damage.pdf / .png\n")

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

ggsave("results/poster_nadir_barplot.pdf", p_bar, width = 8, height = 5, dpi = 300)
ggsave("results/poster_nadir_barplot.png", p_bar, width = 8, height = 5, dpi = 300)
cat("  -> results/poster_nadir_barplot.pdf / .png\n")

cat("\nTous les graphiques generes dans results/\n")

# ════════════════════════════════════════════════════════
# TABLEAU POSTER — Nadirs & validation PK
# Export PNG prêt à intégrer dans le poster
# ════════════════════════════════════════════════════════
library(gridExtra)
library(grid)

# ── 1. Tableau des nadirs ────────────────────────────────
cells_info <- list(
  list(var = "Neut", label = "Neutrophiles\n(10⁹/L)", base = init_pars$Neut0),
  list(var = "Plt",  label = "Plaquettes\n(10⁹/L)",   base = init_pars$Plt0),
  list(var = "Ret",  label = "Réticulocytes\n(10⁹/L)", base = init_pars$Ret0),
  list(var = "RBC",  label = "GR\n(10⁹/L)",            base = init_pars$RBC0)
)

nadir_tbl <- do.call(rbind, lapply(seq_along(doses_fgfr2), function(i) {
  s <- sims_fgfr2[[i]]
  row <- data.frame(Dose = paste0(doses_fgfr2[i], " mg/kg"), stringsAsFactors = FALSE)
  for (ci in cells_info) {
    idx <- which.min(s[[ci$var]])
    val <- round(s[[ci$var]][idx], 1)
    day <- round(s$time_d[idx], 0)
    pct <- round((s[[ci$var]][idx] / ci$base - 1) * 100, 0)
    row[[ci$label]] <- sprintf("%.1f  (J%d  |  %+d%%)", val, day, pct)
  }
  row
}))

# Ligne baseline
base_row <- data.frame(Dose = "Baseline", stringsAsFactors = FALSE)
for (ci in cells_info) base_row[[ci$label]] <- round(ci$base, 1)
nadir_tbl <- rbind(nadir_tbl, base_row)

tt_nad <- ttheme_minimal(
  core = list(
    fg_params  = list(cex = 0.88, col = c(rep("black", 3), "#555555")),
    bg_params  = list(fill = c("#dce8f5", "#edf4ff", "#f5f9ff", "#f7f7f7"))
  ),
  colhead = list(
    fg_params = list(cex = 0.90, fontface = "bold", col = "white"),
    bg_params = list(fill = "#1a3a5c")
  ),
  rowhead = list(fg_params = list(cex = 0.88))
)

tbl_nad_grob <- tableGrob(nadir_tbl, rows = NULL, theme = tt_nad)

title_nad <- textGrob(
  "Nadirs hématologiques prédits — FGFR2 inhibiteur NHP  (Q3W × 3 cycles)",
  gp = gpar(fontsize = 11, fontface = "bold", col = "#1a3a5c")
)
note_nad <- textGrob(
  "Format : valeur nadir  (jour du nadir  |  % vs baseline)",
  gp = gpar(fontsize = 8.5, col = "grey45", fontface = "italic")
)

png("results/poster_table_nadirs.png",
    width = 1100, height = 300, res = 150)
grid.draw(arrangeGrob(title_nad, tbl_nad_grob, note_nad,
                      heights = c(0.13, 0.74, 0.08),
                      padding = unit(4, "mm")))
dev.off()
cat("  -> results/poster_table_nadirs.png\n")

# ── 2. Tableau validation PK vs FDA Table 7 ─────────────
pk_tbl <- do.call(rbind, lapply(seq_along(doses_fgfr2), function(i) {
  s   <- sims_fgfr2[[i]]
  fda <- if (i <= length(fda_tk_nhp)) fda_tk_nhp[[i]] else NULL
  C0_sim  <- round(max(s$C_ADC1[s$time <= 24]), 0)
  idx_auc <- s$time <= 504
  AUC_sim <- round(sum(diff(s$time[idx_auc]) *
                   (s$C_ADC1[idx_auc][-sum(idx_auc)] +
                    s$C_ADC1[idx_auc][-1]) / 2) / 24, 0)
  idt <- s$time >= 100 & s$time <= 480 & s$C_ADC1 > 0
  t12_sim <- if (sum(idt) > 5)
    round(log(2) / abs(coef(lm(log(C_ADC1) ~ time, data = s[idt,]))[2]) / 24, 1)
  else NA_real_

  data.frame(
    Dose         = paste0(doses_fgfr2[i], " mg/kg"),
    `C0 obs`     = if (!is.null(fda)) round(fda$C0_ADC, 0)    else NA_real_,
    `C0 sim`     = C0_sim,
    `AUC21 obs`  = if (!is.null(fda)) round(fda$AUC21d_ADC,0) else NA_real_,
    `AUC21 sim`  = AUC_sim,
    `t1/2 obs`   = if (!is.null(fda)) round(fda$t_half_d, 1)  else NA_real_,
    `t1/2 sim`   = t12_sim,
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
}))

tt_pk <- ttheme_minimal(
  core = list(
    fg_params = list(cex = 0.88),
    bg_params = list(fill = c("#dce8f5", "#edf4ff", "#f5f9ff"))
  ),
  colhead = list(
    fg_params = list(cex = 0.90, fontface = "bold", col = "white"),
    bg_params = list(fill = "#1a3a5c")
  )
)

tbl_pk_grob <- tableGrob(pk_tbl, rows = NULL, theme = tt_pk)

title_pk <- textGrob(
  "Validation PK FGFR2 inhibiteur NHP — vs données précliniques",
  gp = gpar(fontsize = 11, fontface = "bold", col = "#1a3a5c")
)
note_pk <- textGrob(
  "C₀ (µg/mL)  |  AUC₂₁ (µg·h/mL)  |  t½ (jours)   — obs = données précliniques",
  gp = gpar(fontsize = 8.5, col = "grey45", fontface = "italic")
)

png("results/poster_table_pk_validation.png",
    width = 900, height = 250, res = 150)
grid.draw(arrangeGrob(title_pk, tbl_pk_grob, note_pk,
                      heights = c(0.13, 0.74, 0.08),
                      padding = unit(4, "mm")))
dev.off()
cat("  -> results/poster_table_pk_validation.png\n")

# ════════════════════════════════════════════════════════
# Figure 5 — Décalage cinétique PK → hématotoxicité
# ────────────────────────────────────────────────────────
# 4 panneaux (un par dose) : PK normalisée (% Cmax) +
# Neut / Plt / Ret normalisés (% baseline) sur le même axe.
# L'œil voit directement le retard cinétique.
# ════════════════════════════════════════════════════════

norm_long <- bind_rows(lapply(seq_along(doses_fgfr2), function(i) {
  s    <- sims_fgfr2[[i]]
  dose <- paste0(doses_fgfr2[i], " mg/kg")
  Cmax <- max(s$C_ADC1)
  bind_rows(
    data.frame(time_d = s$time_d, Dose = dose,
               Variable = "FGFR2 inhib. (% Cmax)",
               Valeur   = s$C_ADC1 / Cmax * 100),
    data.frame(time_d = s$time_d, Dose = dose,
               Variable = "Neutrophiles (% baseline)",
               Valeur   = s$Neut / init_pars$Neut0 * 100),
    data.frame(time_d = s$time_d, Dose = dose,
               Variable = "Plaquettes (% baseline)",
               Valeur   = s$Plt  / init_pars$Plt0  * 100),
    data.frame(time_d = s$time_d, Dose = dose,
               Variable = "Réticulocytes (% baseline)",
               Valeur   = s$Ret  / init_pars$Ret0  * 100)
  )
})) %>%
  mutate(
    Dose = factor(Dose, levels = paste0(doses_fgfr2, " mg/kg")),
    Variable = factor(Variable, levels = c(
      "FGFR2 inhib. (% Cmax)",
      "Neutrophiles (% baseline)",
      "Plaquettes (% baseline)",
      "Réticulocytes (% baseline)"
    ))
  )

var_cols_kin <- c(
  "FGFR2 inhib. (% Cmax)"     = "#333333",
  "Neutrophiles (% baseline)"  = "#2166ac",
  "Plaquettes (% baseline)"    = "#4dac26",
  "Réticulocytes (% baseline)" = "#d6604d"
)
var_lty_kin <- c("FGFR2 inhib. (% Cmax)" = "solid",
                 "Neutrophiles (% baseline)"  = "solid",
                 "Plaquettes (% baseline)"    = "dashed",
                 "Réticulocytes (% baseline)" = "dotdash")

p_kinetics <- ggplot(norm_long,
                     aes(x = time_d, y = Valeur,
                         color = Variable, linetype = Variable)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey80", linewidth = 0.35) +
  geom_hline(yintercept = 100, linetype = "dotted",
             color = "grey50", linewidth = 0.55) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(values = var_cols_kin) +
  scale_linetype_manual(values = var_lty_kin) +
  scale_x_continuous(breaks = seq(0, 120, by = 21),
                     labels = paste0("J", seq(0, 120, by = 21))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, NA)) +
  facet_wrap(~ Dose, ncol = 2) +
  labs(
    title    = "Décalage cinétique PK → hématotoxicité — FGFR2 inhibiteur NHP",
    subtitle = "FGFR2 inhib. : % du Cmax  |  cellules : % de la valeur basale  |  --- jour de dose",
    x = "Temps (jours)", y = "% (normalisé)"
  ) +
  theme_poster +
  theme(legend.position = "bottom",
        legend.title    = element_blank())

ggsave("results/poster_PK_PD_kinetics.pdf",
       p_kinetics, width = 13, height = 10, dpi = 300)
ggsave("results/poster_PK_PD_kinetics.png",
       p_kinetics, width = 13, height = 10, dpi = 300)
cat("  -> results/poster_PK_PD_kinetics.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure 6 — 2 rangées synchronisées : PK (log) / Neut (lin)
# Vue toutes doses — montre clairement le retard de nadir
# ════════════════════════════════════════════════════════
pk_long_all <- bind_rows(lapply(seq_along(doses_fgfr2), function(i) {
  sims_fgfr2[[i]] %>%
    select(time_d, C_ADC1) %>%
    mutate(Dose = factor(paste0(doses_fgfr2[i], " mg/kg"),
                         levels = paste0(doses_fgfr2, " mg/kg")))
}))

neut_long_all <- bind_rows(lapply(seq_along(doses_fgfr2), function(i) {
  sims_fgfr2[[i]] %>%
    select(time_d, Neut, Plt, Ret) %>%
    mutate(Dose = factor(paste0(doses_fgfr2[i], " mg/kg"),
                         levels = paste0(doses_fgfr2, " mg/kg")))
})) %>%
  pivot_longer(cols = c(Neut, Plt, Ret),
               names_to = "Cellule", values_to = "Valeur") %>%
  mutate(Cellule = factor(Cellule,
                          levels = c("Neut", "Plt", "Ret"),
                          labels = c("Neutrophiles", "Plaquettes", "Réticulocytes")))

xbreaks <- seq(0, 120, by = 21)
xlabels <- paste0("J", xbreaks)

p_row1 <- ggplot(pk_long_all, aes(x = time_d, y = C_ADC1, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey75", linewidth = 0.35) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = dose_cols) +
  scale_y_log10() +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  labs(title = "PK — FGFR2 inhibiteur (échelle log)",
       x = NULL, y = "Concentration (µg/mL)") +
  theme_poster +
  theme(legend.position = "none",
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(5, 5, 0, 5))

p_row2 <- ggplot(neut_long_all,
                 aes(x = time_d, y = Valeur, color = Dose, linetype = Cellule)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey75", linewidth = 0.35) +
  geom_hline(data = data.frame(
    Cellule = factor(c("Neutrophiles","Plaquettes","Réticulocytes")),
    base    = c(init_pars$Neut0, init_pars$Plt0, init_pars$Ret0)),
    aes(yintercept = base), color = "grey40", linetype = "dotted",
    linewidth = 0.55, inherit.aes = FALSE) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = dose_cols) +
  scale_linetype_manual(values = c("Neutrophiles"  = "solid",
                                   "Plaquettes"    = "dashed",
                                   "Réticulocytes" = "dotdash")) +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  facet_wrap(~ Dose, ncol = 4) +
  labs(title    = "Réponse hématologique (nadir retardé vs PK)",
       subtitle = "··· baseline  |  --- jour de dose",
       x = "Temps (jours)", y = "Cellules (10⁹/L)",
       linetype = "Lignée") +
  theme_poster +
  theme(legend.position = "bottom",
        plot.margin     = margin(0, 5, 5, 5))

# Assembler avec gridExtra
g_combined2rows <- gridExtra::arrangeGrob(p_row1, p_row2,
                                          nrow = 2, heights = c(1, 2.2))
ggsave("results/poster_PK_PD_2rows.pdf",
       g_combined2rows, width = 14, height = 10, dpi = 300)
ggsave("results/poster_PK_PD_2rows.png",
       g_combined2rows, width = 14, height = 10, dpi = 300)
cat("  -> results/poster_PK_PD_2rows.pdf / .png\n")

# ════════════════════════════════════════════════════════
# CONCLUSIONS POSTER — texte révisé
# Copier-coller dans l'outil de mise en page du poster
# ════════════════════════════════════════════════════════
# • A semi-mechanistic PK/PD model was developed for the FGFR2 inhibitor
#   in NHP by coupling a 2-compartment TMDD PK model with the Fornari 2019
#   hematopoietic progenitor framework
#
# • PK parameters were calibrated and validated against preclinical NHP data
#   across 3 dose levels (3, 10, 30 mg/kg Q3W × 3 cycles)
#
# • The model predicts dose-dependent hematological nadirs:
#   reticulocytes at Day ~8, neutrophils and platelets at Day ~15 post-dose
#
# • At 30 mg/kg, deepest suppression is predicted for reticulocytes and
#   neutrophils — consistent with dose-limiting toxicity observed in NHP
#
# • This mechanistic framework provides a quantitative basis for
#   preclinical-to-clinical translation of FGFR2 inhibitor hematotoxicity
#   and supports rational dose optimization in oncology
