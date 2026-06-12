############################################################
# plots_fgfr2_nhp.R
# Graphiques publication-ready -- ggplot2 uniquement
# Prérequis : run_fgfr2_nhp.R déjà sourcé
############################################################
library(ggplot2)
library(dplyr)
library(tidyr)

# -- Auto-source si objets manquants ---------------------
if (!exists("sims_fgfr2") || !exists("doses_fgfr2")) {
  cat("Lancement de run_fgfr2_nhp.R...\n")
  source("run_fgfr2_nhp.R")
}

if (!dir.exists("results")) dir.create("results")

# ════════════════════════════════════════════════════════
# DONNÉES OBSERVÉES -- chargées depuis nhp_hema_data.csv
# Remplir le CSV avec les valeurs réelles (cf. data_nhp.R)
# ════════════════════════════════════════════════════════
source("data_nhp.R")

# -- Palette ---------------------------------------------
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

# -- Long format ------------------------------------------
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
                     labels = c("Neutrophils (10⁹/L)", "Platelets (10⁹/L)",
                                "Reticulocytes (10⁹/L)", "RBC (10⁹/L)",
                                "Damage (a.u.)"))
  )

# Baselines par panel
base_df <- data.frame(
  Cellule = c("Neutrophils (10⁹/L)", "Platelets (10⁹/L)",
              "Reticulocytes (10⁹/L)", "RBC (10⁹/L)"),
  baseline = c(init_pars$Neut0, init_pars$Plt0,
               init_pars$Ret0,  init_pars$RBC0)
)

# Nadirs par panel + dose
nadir_df <- sim_long %>%
  filter(Cellule != "Damage (a.u.)") %>%
  group_by(Dose, Cellule) %>%
  slice_min(Valeur, n = 1) %>%
  ungroup()

# ════════════════════════════════════════════════════════
# Figure 1 -- 4 panels cellulaires (facet_wrap)
# ════════════════════════════════════════════════════════
df_cells <- sim_long %>% filter(Cellule != "Damage (a.u.)")

# Observed data en long format (si rempli)
obs_has_data <- nrow(obs_data) > 0 &&
                any(!is.na(obs_data[, c("Neut","Plt","RBC","Ret")]))

obs_long <- if (obs_has_data) {
  obs_data %>%
    mutate(Dose = paste0(dose_mgkg, " mg/kg")) %>%
    pivot_longer(cols = c(Neut, Plt, Ret, RBC),
                 names_to = "Cellule", values_to = "Valeur") %>%
    mutate(
      Dose      = factor(Dose, levels = paste0(doses_fgfr2, " mg/kg")),
      Animal_Id = factor(Animal_Id),
      Cellule   = factor(Cellule,
                         levels = c("Neut", "Plt", "Ret", "RBC"),
                         labels = c("Neutrophils (10⁹/L)", "Platelets (10⁹/L)",
                                    "Reticulocytes (10⁹/L)", "RBC (10⁹/L)"))
    ) %>%
    filter(!is.na(Valeur))
} else NULL

# Données observées retirées pour confidentialité (Pierre Fabre)
# Le modèle reproduit bien les tendances — voir §3.3.3 pour discussion
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
  scale_x_continuous(breaks = seq(0, 120, by = 21),
                     labels = paste0("D", seq(0, 120, by = 21))) +
  facet_wrap(~Cellule, scales = "free_y", ncol = 2) +
  labs(title    = "Profils hématologiques prédits — composé en développement interne, Q3W × 3 cycles",
       subtitle = "▽ nadir  |  ··· baseline  |  --- jour de dose  |  données observées non reproduites (confidentielles)",
       x = "Temps (jours)", y = NULL) +
  theme_poster

ggsave("results/poster_PD_4panels.pdf",
       p_cells, width = 13, height = 10, dpi = 300)
ggsave("results/poster_PD_4panels.png",
       p_cells, width = 13, height = 10, dpi = 300)
cat("  -> results/poster_PD_4panels.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure 2 -- Damage seul
# ════════════════════════════════════════════════════════
p_damage <- ggplot(sim_long %>% filter(Cellule == "Damage (a.u.)"),
                   aes(x = time_d, y = Valeur, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = dose_cols) +
  scale_x_continuous(breaks = seq(0, 120, by = 21),
                     labels = paste0("D", seq(0, 120, by = 21))) +
  labs(title = "DNA Damage -- model driver",
       x = "Time (days)", y = "Damage (a.u.)") +
  theme_poster

ggsave("results/poster_Damage.pdf", p_damage, width = 7, height = 5, dpi = 300)
ggsave("results/poster_Damage.png", p_damage, width = 7, height = 5, dpi = 300)
cat("  -> results/poster_Damage.pdf / .png\n")

# ════════════════════════════════════════════════════════
# Figure 3 -- Barplot nadir % changement
# ════════════════════════════════════════════════════════
base_vals <- c("Neutrophils (10⁹/L)" = init_pars$Neut0,
               "Platelets (10⁹/L)"   = init_pars$Plt0,
               "Reticulocytes (10⁹/L)"= init_pars$Ret0,
               "RBC (10⁹/L)"           = init_pars$RBC0)

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
  labs(title = "Nadir -- % change from baseline",
       x = NULL, y = "% change vs baseline") +
  theme_poster

ggsave("results/poster_nadir_barplot.pdf", p_bar, width = 8, height = 5, dpi = 300)
ggsave("results/poster_nadir_barplot.png", p_bar, width = 8, height = 5, dpi = 300)
cat("  -> results/poster_nadir_barplot.pdf / .png\n")

cat("\nTous les graphiques generes dans results/\n")

# ════════════════════════════════════════════════════════
# TABLEAU POSTER -- Nadirs & validation PK
# Export PNG prêt à intégrer dans le poster
# ════════════════════════════════════════════════════════
library(gridExtra)
library(grid)

# -- 1. Tableau des nadirs --------------------------------
cells_info <- list(
  list(var = "Neut", label = "Neutrophils\n(10⁹/L)",   base = init_pars$Neut0),
  list(var = "Plt",  label = "Platelets\n(10⁹/L)",     base = init_pars$Plt0),
  list(var = "Ret",  label = "Reticulocytes\n(10⁹/L)", base = init_pars$Ret0),
  list(var = "RBC",  label = "RBC\n(10⁹/L)",           base = init_pars$RBC0)
)

nadir_tbl <- do.call(rbind, lapply(seq_along(doses_fgfr2), function(i) {
  s <- sims_fgfr2[[i]]
  row <- data.frame(Dose = paste0(doses_fgfr2[i], " mg/kg"), stringsAsFactors = FALSE)
  for (ci in cells_info) {
    idx <- which.min(s[[ci$var]])
    val <- round(s[[ci$var]][idx], 1)
    day <- round(s$time_d[idx], 0)
    pct <- round((s[[ci$var]][idx] / ci$base - 1) * 100, 0)
    row[[ci$label]] <- sprintf("%.1f  (D%d  |  %+d%%)", val, day, pct)
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
  "Predicted Hematological Nadirs -- FGFR2 inhibitor NHP  (Q3W × 3 cycles)",
  gp = gpar(fontsize = 11, fontface = "bold", col = "#1a3a5c")
)
note_nad <- textGrob(
  "Format: nadir value  (day of nadir  |  % vs baseline)",
  gp = gpar(fontsize = 8.5, col = "grey45", fontface = "italic")
)

png("results/poster_table_nadirs.png",
    width = 1100, height = 300, res = 150)
grid.draw(arrangeGrob(title_nad, tbl_nad_grob, note_nad,
                      heights = c(0.13, 0.74, 0.08),
                      padding = unit(4, "mm")))
dev.off()
cat("  -> results/poster_table_nadirs.png\n")

# -- 2. Tableau validation PK vs FDA Table 7 -------------
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
  "PK Validation -- FGFR2 inhibitor NHP vs preclinical data",
  gp = gpar(fontsize = 11, fontface = "bold", col = "#1a3a5c")
)
note_pk <- textGrob(
  "C₀ (µg/mL)  |  AUC₂₁ (µg·h/mL)  |  t½ (days)   -- obs = preclinical data",
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
# Figure 5 -- Décalage cinétique PK → hématotoxicité
# --------------------------------------------------------
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
               Variable = "Neutrophils (% baseline)",
               Valeur   = s$Neut / init_pars$Neut0 * 100),
    data.frame(time_d = s$time_d, Dose = dose,
               Variable = "Platelets (% baseline)",
               Valeur   = s$Plt  / init_pars$Plt0  * 100),
    data.frame(time_d = s$time_d, Dose = dose,
               Variable = "Reticulocytes (% baseline)",
               Valeur   = s$Ret  / init_pars$Ret0  * 100)
  )
})) %>%
  mutate(
    Dose = factor(Dose, levels = paste0(doses_fgfr2, " mg/kg")),
    Variable = factor(Variable, levels = c(
      "FGFR2 inhib. (% Cmax)",
      "Neutrophils (% baseline)",
      "Platelets (% baseline)",
      "Reticulocytes (% baseline)"
    ))
  )

var_cols_kin <- c(
  "FGFR2 inhib. (% Cmax)"     = "#333333",
  "Neutrophils (% baseline)"  = "#2166ac",
  "Platelets (% baseline)"    = "#4dac26",
  "Reticulocytes (% baseline)" = "#d6604d"
)
var_lty_kin <- c("FGFR2 inhib. (% Cmax)" = "solid",
                 "Neutrophils (% baseline)"  = "solid",
                 "Platelets (% baseline)"    = "dashed",
                 "Reticulocytes (% baseline)" = "dotdash")

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
                     labels = paste0("D", seq(0, 120, by = 21))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, NA)) +
  facet_wrap(~ Dose, ncol = 2) +
  labs(
    title    = "PK → Hematotoxicity Kinetic Shift -- FGFR2 inhibitor NHP",
    subtitle = "FGFR2 inhib.: % of Cmax  |  cells: % of baseline  |  --- dose day",
    x = "Time (days)", y = "% (normalized)"
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
# Figure 6 -- 2 rangées synchronisées : PK (log) / Neut (lin)
# Vue toutes doses -- montre clairement le retard de nadir
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
                          labels = c("Neutrophils", "Platelets", "Reticulocytes")))

xbreaks <- seq(0, 120, by = 21)
xlabels <- paste0("D", xbreaks)

p_row1 <- ggplot(pk_long_all, aes(x = time_d, y = C_ADC1, color = Dose)) +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey75", linewidth = 0.35) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = dose_cols) +
  scale_y_log10() +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  labs(title = "PK -- FGFR2 inhibitor (log scale)",
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
    Cellule = factor(c("Neutrophils","Platelets","Reticulocytes")),
    base    = c(init_pars$Neut0, init_pars$Plt0, init_pars$Ret0)),
    aes(yintercept = base), color = "grey40", linetype = "dotted",
    linewidth = 0.55, inherit.aes = FALSE) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = dose_cols) +
  scale_linetype_manual(values = c("Neutrophils"   = "solid",
                                   "Platelets"     = "dashed",
                                   "Reticulocytes" = "dotdash")) +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  facet_wrap(~ Dose, ncol = 4) +
  labs(title    = "Hematological Response (delayed nadir vs PK)",
       subtitle = "··· baseline  |  --- dose day",
       x = "Time (days)", y = "Cells (10⁹/L)",
       linetype = "Cell type") +
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
# Figure 7 -- Profils PREDITS individuels (% baseline propre)
# Une simulation par animal (baseline = valeur J-3 individuelle)
# 8 panneaux, 4 lignees, J-3 = 100%
# Points observes superposes
# ════════════════════════════════════════════════════════
if (obs_has_data) {

  # -- Helper : re-derive Eq. S4 pour une baseline individuelle --
  derive_pars_individual <- function(ip_base, neut0, plt0, rbc0, ret0) {
    ip <- ip_base
    ip$Neut0 <- neut0
    ip$Plt0  <- plt0
    ip$RBC0  <- rbc0
    ip$Ret0  <- ret0
    # Eq. S4
    ip$k_circ_Ret <- ip$k_circ_RBC * ip$RBC0  / ip$Ret0
    ip$a_Ret      <- 3 / ip$MTT_Ret
    ip$a_Plt      <- 3 / ip$MTT_Plt
    ip$k_tr_Neut  <- ip$k_circ_Neut * ip$Neut0 / ip$CMP0
    ip$k_tr_Mono  <- ip$k_circ_Mono * ip$Mono0 / ip$CMP0
    ip$k_tr_Ret   <- ip$k_circ_Ret  * ip$Ret0  / (ip$lambda1^2 * ip$MEP0)
    ip$k_tr_Plt   <- ip$k_circ_Plt  * ip$Plt0  / (ip$lambda2^2 * ip$MEP0)
    ip$k_prol_Ret <- ip$a_Ret * (1 - 1/ip$lambda1)
    ip$k_prol_Plt <- ip$a_Plt * (1 - 1/ip$lambda2)
    ip$k_prol_CMP <- (ip$k_tr_Neut + ip$k_tr_Mono) / ip$lambda3
    ip$k_prol_MEP <- (ip$k_tr_Ret  + ip$k_tr_Plt)  / ip$lambda4
    ip$k_tr_CMP   <- (ip$k_tr_Neut + ip$k_tr_Mono - ip$k_prol_CMP) * ip$CMP0 / ip$MPP0
    ip$k_tr_MEP   <- (ip$k_tr_Ret  + ip$k_tr_Plt  - ip$k_prol_MEP) * ip$MEP0 / ip$MPP0
    ip$k_prol_MPP <- (ip$k_tr_CMP  + ip$k_tr_MEP)  / ip$lambda5
    ip$k_stem     <- (ip$k_tr_CMP  + ip$k_tr_MEP - ip$k_prol_MPP) * ip$MPP0
    # Transit states initiaux (Eq. S3)
    a_Neut <- 3 / ip$MTT_Neut
    a_Mono <- 3 / ip$MTT_Mono
    T_Neut <- ip$k_circ_Neut * ip$Neut0 / a_Neut
    T_Mono <- ip$k_circ_Mono * ip$Mono0 / a_Mono
    T2_Ret <- ip$k_circ_Ret  * ip$Ret0  / ip$a_Ret
    T1_Ret <- T2_Ret / ip$lambda1
    T2_Plt <- ip$k_circ_Plt  * ip$Plt0  / ip$a_Plt
    T1_Plt <- T2_Plt / ip$lambda2
    attr(ip, "state0_pd") <- c(
      MPP = ip$MPP0, CMP = ip$CMP0, MEP = ip$MEP0,
      T1_Neut=T_Neut, T2_Neut=T_Neut, T3_Neut=T_Neut, Neut=ip$Neut0,
      T1_Mono=T_Mono, T2_Mono=T_Mono, T3_Mono=T_Mono, Mono=ip$Mono0,
      T1_Ret=T1_Ret,  T2_Ret=T2_Ret,  T3_Ret=T2_Ret,  Ret=ip$Ret0,
      RBC=ip$RBC0,
      T1_Plt=T1_Plt,  T2_Plt=T2_Plt,  T3_Plt=T2_Plt,  Plt=ip$Plt0
    )
    ip
  }

  # -- Simulation individuelle par animal ----------------
  animal_info <- list(
    list(id="1001", dose=4),  list(id="1002", dose=4),
    list(id="2001", dose=13), list(id="2002", dose=13),
    list(id="4001", dose=26), list(id="4002", dose=26),
    list(id="3101", dose=39), list(id="3002", dose=39)
  )

  animal_order <- c("1001 (4 mg/kg)",  "1002 (4 mg/kg)",
                    "2001 (13 mg/kg)", "2002 (13 mg/kg)",
                    "4001 (26 mg/kg)", "4002 (26 mg/kg)",
                    "3101 (39 mg/kg)", "3002 (39 mg/kg)")

  cat("\nSimulations individuelles (8 animaux)...\n")
  sims_ind <- lapply(animal_info, function(a) {
    base <- obs_data[obs_data$Animal_Id == a$id & obs_data$jour == -3, ]
    if (nrow(base) == 0 || any(is.na(c(base$Neut, base$Plt, base$RBC, base$Ret)))) {
      cat(sprintf("  %s : baseline J-3 manquante -- ignore\n", a$id))
      return(NULL)
    }
    ip <- derive_pars_individual(init_pars,
                                 neut0 = base$Neut, plt0 = base$Plt,
                                 rbc0  = base$RBC,  ret0 = base$Ret)
    # Ajouter params FGFR2
    for (nm in names(fgfr2_nhp)) ip[[nm]] <- fgfr2_nhp[[nm]]
    ip$k_dam_DXd <- 0.017
    ip$rate_fun  <- make_nhp_infusion(
      dose_mgkg  = a$dose, BW_kg = 4.0,
      Tinfu_h    = 0.5,
      interval_h = fgfr2_nhp$interval_h,
      n_cycles   = 3)
    state0_pd <- attr(ip, "state0_pd")
    state0 <- c(C_ADC1=0, C_ADC2=0, Damage=0, state0_pd)
    times  <- seq(0, 90*24, by = 1)
    sol <- tryCatch(
      as.data.frame(ode(y=state0, times=times,
                        func=pkpd_nhp_ode, parms=ip,
                        method="lsoda", hmax=0.25)),
      error = function(e) { cat(sprintf("  %s erreur: %s\n", a$id, e$message)); NULL }
    )
    if (is.null(sol)) return(NULL)
    sol$time_d    <- sol$time / 24
    sol$Animal_Id <- a$id
    sol$dose_mgkg <- a$dose
    sol$Neut_base <- base$Neut
    sol$Plt_base  <- base$Plt
    sol$RBC_base  <- base$RBC
    sol$Ret_base  <- base$Ret
    cat(sprintf("  %s (%d mg/kg) OK\n", a$id, a$dose))
    sol
  })
  names(sims_ind) <- sapply(animal_info, `[[`, "id")

  # -- Long format normalise ------------------------------
  ind_pred_long <- bind_rows(lapply(sims_ind, function(s) {
    if (is.null(s)) return(NULL)
    lbl <- paste0(s$Animal_Id[1], " (", s$dose_mgkg[1], " mg/kg)")
    data.frame(
      Animal_label = factor(lbl, levels = animal_order),
      time_d  = s$time_d,
      Neutrophils   = s$Neut / s$Neut_base[1] * 100,
      Platelets     = s$Plt  / s$Plt_base[1]  * 100,
      RBC           = s$RBC  / s$RBC_base[1]  * 100,
      Reticulocytes = s$Ret  / s$Ret_base[1]  * 100
    )
  })) %>%
    pivot_longer(cols = c(Neutrophils, Platelets, RBC, Reticulocytes),
                 names_to = "Cellule", values_to = "Pct") %>%
    mutate(Cellule = factor(Cellule,
                            levels = c("Neutrophils","Platelets",
                                       "RBC","Reticulocytes")))

  # -- Observations en % baseline individuelle ------------
  baseline_ind <- obs_data %>%
    filter(jour == -3) %>%
    select(Animal_Id, Neut_base=Neut, Plt_base=Plt,
           RBC_base=RBC, Ret_base=Ret)

  obs_ind_pct <- obs_data %>%
    left_join(baseline_ind, by="Animal_Id") %>%
    mutate(
      Neutrophils   = Neut / Neut_base * 100,
      Platelets     = Plt  / Plt_base  * 100,
      RBC           = RBC  / RBC_base  * 100,
      Reticulocytes = Ret  / Ret_base  * 100,
      Animal_label  = factor(
        paste0(Animal_Id, " (", dose_mgkg, " mg/kg)"),
        levels = animal_order)
    ) %>%
    select(Animal_label, jour,
           Neutrophils, Platelets, RBC, Reticulocytes) %>%
    pivot_longer(cols = c(Neutrophils, Platelets, RBC, Reticulocytes),
                 names_to = "Cellule", values_to = "Pct") %>%
    mutate(Cellule = factor(Cellule,
                            levels = c("Neutrophils","Platelets",
                                       "RBC","Reticulocytes"))) %>%
    filter(!is.na(Pct))

  # -- Palette -------------------------------------------
  cell_cols_ind <- c("Neutrophils"   = "#2166ac",
                     "Platelets"     = "#4dac26",
                     "RBC"           = "#d6604d",
                     "Reticulocytes" = "#984ea3")

  # -- Figure --------------------------------------------
  p_ind_pred <- ggplot(ind_pred_long,
                       aes(x = time_d, y = Pct,
                           color = Cellule, group = Cellule)) +
    geom_hline(yintercept = 100, linetype = "dashed",
               color = "grey45", linewidth = 0.6) +
    geom_vline(xintercept = dose_days, linetype = "dotted",
               color = "grey75", linewidth = 0.35) +
    geom_line(linewidth = 1.1) +
    geom_point(data = obs_ind_pct,
               aes(x = jour, y = Pct, color = Cellule),
               shape = 19, size = 2.4, inherit.aes = FALSE) +
    scale_color_manual(values = cell_cols_ind, name = "Cell type") +
    scale_x_continuous(breaks = c(-3, 0, 2, 8, 12, 15, 22, 25),
                       labels = c("D-3","D0","D2","D8","D12","D15","D22","D25")) +
    coord_cartesian(xlim = c(-3, 25)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    facet_wrap(~ Animal_label, ncol = 4) +
    labs(
      title    = "Individual Predicted Profiles -- FGFR2 inhibitor NHP",
      subtitle = "Line = model prediction  |  Points = observations  |  % of individual day-3 baseline  |  --- 100%  |  ··· dose day",
      x = "Time (days)", y = "% of individual baseline"
    ) +
    theme_poster +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  ggsave("results/poster_PD_predicted_individual.pdf",
         p_ind_pred, width = 15, height = 9, dpi = 300)
  ggsave("results/poster_PD_predicted_individual.png",
         p_ind_pred, width = 15, height = 9, dpi = 300)
  cat("  -> results/poster_PD_predicted_individual.pdf / .png\n")

  # ════════════════════════════════════════════════════════
  # Figure 8 -- Profils PREDITS individuels SANS observé
  # Identique à Fig 7 mais sans geom_point observations
  # ════════════════════════════════════════════════════════
  p_ind_pred_noobs <- ggplot(ind_pred_long,
                             aes(x = time_d, y = Pct,
                                 color = Cellule, group = Cellule)) +
    geom_hline(yintercept = 100, linetype = "dashed",
               color = "grey45", linewidth = 0.6) +
    geom_vline(xintercept = dose_days, linetype = "dotted",
               color = "grey75", linewidth = 0.35) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = cell_cols_ind, name = "Cell type") +
    scale_x_continuous(breaks = c(0, 21, 42, 56, 70, 84),
                       labels = c("D0","D21","D42","D56","D70","D84")) +
    coord_cartesian(xlim = c(0, 90)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    facet_wrap(~ Animal_label, ncol = 4) +
    labs(
      title    = "Individual Predicted Profiles -- FGFR2 inhibitor NHP",
      subtitle = "Model prediction  |  % of individual day-3 baseline  |  --- 100%  |  ··· dose day (D0/D21/D42)",
      x = "Time (days)", y = "% of individual baseline"
    ) +
    theme_poster +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  ggsave("results/poster_PD_predicted_individual_noobs.pdf",
         p_ind_pred_noobs, width = 15, height = 9, dpi = 300)
  ggsave("results/poster_PD_predicted_individual_noobs.png",
         p_ind_pred_noobs, width = 15, height = 9, dpi = 300)
  cat("  -> results/poster_PD_predicted_individual_noobs.pdf / .png\n")

} else {
  cat("  [Figure 7 ignoree : nhp_hema_data.csv vide]\n")
}

# ════════════════════════════════════════════════════════
# CONCLUSIONS POSTER -- texte révisé
# Copier-coller dans l'outil de mise en page du poster
# ════════════════════════════════════════════════════════
# * A semi-mechanistic PK/PD model was developed for the FGFR2 inhibitor
#   in NHP by coupling a 2-compartment TMDD PK model with the Fornari 2019
#   hematopoietic progenitor framework
#
# * PK parameters were calibrated and validated against preclinical NHP data
#   across 3 dose levels (3, 10, 30 mg/kg Q3W × 3 cycles)
#
# * The model predicts dose-dependent hematological nadirs:
#   reticulocytes at Day ~8, neutrophils and platelets at Day ~15 post-dose
#
# * At 30 mg/kg, deepest suppression is predicted for reticulocytes and
#   neutrophils -- consistent with dose-limiting toxicity observed in NHP
#
# * This mechanistic framework provides a quantitative basis for
#   preclinical-to-clinical translation of FGFR2 inhibitor hematotoxicity
#   and supports rational dose optimization in oncology
