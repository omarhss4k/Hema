############################################################
# run_sensitivity_baselines.R
# Analyse de sensibilité : impact des baselines sur les nadirs
# Varie Neut0, Plt0, Ret0, RBC0 (min / médiane / max) — une à la fois
############################################################
source("pkpd_model_FORNARI.R")
source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")
source("plots_grades.R")

library(ggplot2)
library(gridExtra)
library(grid)

if (!dir.exists("results_HUMAN")) dir.create("results_HUMAN")

# ── Dose Calvert (GFR=78 médian) ──────────────────────
AUC_target   <- 5
GFR_mLmin    <- 78
dose_calvert <- AUC_target * (GFR_mLmin + 25)
times_hu     <- seq(0, 63 * 24, by = 1)
dose_days_hu <- c(0, 21)

# ── Ranges des baselines (littérature, Table 1 Fornari) ──
baseline_ranges <- list(
  Neut0 = c(min = 2.0,    med = 4.5,    max = 7.0),
  Plt0  = c(min = 150.0,  med = 250.0,  max = 400.0),  # range normal adulte (×10⁹/L)
  Ret0  = c(min = 40.0,   med = 70.0,   max = 115.0),
  RBC0  = c(min = 4100.0, med = 5000.0, max = 5900.0)
)

# ── Fonction : simuler et extraire les nadirs ──────────
run_one <- function(pars_mod, dose) {
  pars_mod$rate_fun <- make_repeated_infusion(
    dose_mg = dose, Tinfu_h = 1,
    interval_h = 21 * 24, n_cycles = 2
  )
  sim <- simulate_all(times_hu, pars_mod, init_state)
  list(
    nadir_neut  = min(sim$Neut, na.rm = TRUE),
    grade_neut  = assign_grade_neut(min(sim$Neut, na.rm = TRUE)),
    nadir_plt   = min(sim$Plt,  na.rm = TRUE),
    grade_plt   = assign_grade_plt( min(sim$Plt,  na.rm = TRUE)),
    sim         = sim
  )
}

# ── Référence (valeurs nominales) ─────────────────────
cat("=== Référence (baselines nominales) ===\n")
ref <- run_one(init_pars, dose_calvert)
cat(sprintf("  Neut nadir = %.3f  → Grade %d\n", ref$nadir_neut, ref$grade_neut))
cat(sprintf("  Plt  nadir = %.1f  → Grade %d\n", ref$nadir_plt,  ref$grade_plt))

# ── Boucle sensibilité ────────────────────────────────
results <- list()

for (param in names(baseline_ranges)) {
  for (lvl in c("min", "med", "max")) {
    val     <- baseline_ranges[[param]][lvl]
    pars_i  <- init_pars
    pars_i[[param]] <- val

    # Recalibrer init_state pour ce paramètre
    state_i         <- init_state
    state_i[param]  <- val

    # Pour Neut0 / Plt0 : rééquilibrer les compartiments transit
    if (param == "Neut0") {
      a <- 3 / pars_i$MTT_Neut
      T <- pars_i$k_circ_Neut * val / a
      state_i["T1_Neut"] <- T
      state_i["T2_Neut"] <- T
      state_i["T3_Neut"] <- T
    }
    if (param == "Plt0") {
      a  <- 3 / pars_i$MTT_Plt
      T  <- pars_i$k_circ_Plt * val / a
      T1 <- T / pars_i$lambda2
      state_i["T1_Plt"] <- T1
      state_i["T2_Plt"] <- T
      state_i["T3_Plt"] <- T
    }
    if (param == "Ret0") {
      a <- 3 / pars_i$MTT_Ret
      T <- pars_i$k_circ_RBC * val / a
      state_i["T1_Ret"] <- T
      state_i["T2_Ret"] <- T
      state_i["T3_Ret"] <- T
    }

    pars_i$rate_fun <- make_repeated_infusion(
      dose_mg = dose_calvert, Tinfu_h = 1,
      interval_h = 21 * 24, n_cycles = 2
    )

    sim_i <- simulate_all(times_hu, pars_i, state_i)
    nadir_n <- min(sim_i$Neut, na.rm = TRUE)
    nadir_p <- min(sim_i$Plt,  na.rm = TRUE)
    grade_n <- assign_grade_neut(nadir_n)
    grade_p <- assign_grade_plt(nadir_p)

    results[[paste(param, lvl, sep = "_")]] <- list(
      param      = param,
      level      = lvl,
      value      = val,
      nadir_neut = nadir_n,
      grade_neut = grade_n,
      nadir_plt  = nadir_p,
      grade_plt  = grade_p,
      sim        = sim_i
    )
    cat(sprintf("  %s=%g : Neut nadir=%.3f (G%d) | Plt nadir=%.1f (G%d)\n",
                param, val, nadir_n, grade_n, nadir_p, grade_p))
  }
}

# ── Tableau résumé ────────────────────────────────────
df_res <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    param      = r$param,
    level      = r$level,
    value      = r$value,
    nadir_neut = round(r$nadir_neut, 3),
    grade_neut = r$grade_neut,
    nadir_plt  = round(r$nadir_plt,  1),
    grade_plt  = r$grade_plt,
    stringsAsFactors = FALSE
  )
}))

cat("\n=== Tableau de sensibilité ===\n")
print(df_res, row.names = FALSE)
write.csv(df_res, "results_HUMAN/sensitivity_baselines.csv", row.names = FALSE)
cat("\n✓ Tableau sauvegardé : results_HUMAN/sensitivity_baselines.csv\n")

# ── Figure : courbes Neut par niveau de Neut0 ─────────
make_sensitivity_plot <- function(cell, ylabel, thresholds, color_pal) {
  df_lines <- NULL
  for (param in names(baseline_ranges)) {
    for (lvl in c("min", "med", "max")) {
      key <- paste(param, lvl, sep = "_")
      r   <- results[[key]]
      df_lines <- rbind(df_lines, data.frame(
        days  = r$sim$days,
        value = r$sim[[cell]],
        param = param,
        level = lvl
      ))
    }
  }

  # Ajouter la référence
  df_lines <- rbind(df_lines, data.frame(
    days  = ref$sim$days,
    value = ref$sim[[cell]],
    param = "REF",
    level = "med"
  ))

  df_lines$label <- paste0(df_lines$param, " (", df_lines$level, ")")
  df_lines$label[df_lines$param == "REF"] <- "REF (nominal)"

  level_colors <- c(min = "#2166ac", med = "#4dac26", max = "#d73027")

  ymin_plot <- min(thresholds) * 0.1
  ymax_plot <- max(df_lines$value, na.rm = TRUE) * 3

  p <- ggplot(df_lines[df_lines$param != "REF", ],
              aes(x = days, y = value,
                  color = level, linetype = param, group = label)) +

    # Bandes de grade
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=ymin_plot, ymax=thresholds[4],
             fill=GRADE_COLORS["G4"], alpha=0.10) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=thresholds[4], ymax=thresholds[3],
             fill=GRADE_COLORS["G3"], alpha=0.09) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=thresholds[3], ymax=thresholds[2],
             fill=GRADE_COLORS["G2"], alpha=0.08) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=thresholds[2], ymax=thresholds[1],
             fill=GRADE_COLORS["G1"], alpha=0.07) +

    geom_line(linewidth = 0.7, alpha = 0.75) +

    # Courbe REF en noir épais
    geom_line(data = df_lines[df_lines$param == "REF", ],
              aes(x = days, y = value),
              color = "black", linewidth = 1.2,
              linetype = "solid", inherit.aes = FALSE) +

    # Lignes de seuil
    geom_hline(yintercept = thresholds,
               linetype = "dashed", linewidth = 0.3,
               color = c(GRADE_COLORS["G1"], GRADE_COLORS["G2"],
                         GRADE_COLORS["G3"], GRADE_COLORS["G4"])) +

    # Doses
    geom_vline(xintercept = dose_days_hu,
               linetype = "dotted", color = "grey55",
               linewidth = 0.35, alpha = 0.6) +

    scale_color_manual(values = level_colors, name = "Niveau") +
    scale_linetype_manual(
      values = c(Neut0 = "solid", Plt0 = "dashed",
                 Ret0  = "dotdash", RBC0 = "twodash"),
      name = "Baseline variée"
    ) +
    scale_y_log10(
      limits = c(ymin_plot, ymax_plot),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(title = ylabel,
         x = "Time (d)",
         y = expression(10^9~cells~L^{-1})) +
    theme_bw(base_size = 9.5) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major   = element_line(color = "grey92"),
          plot.title         = element_text(face = "bold", size = 9, hjust = 0.5),
          legend.position    = "right",
          legend.key.width   = unit(1.5, "lines"),
          axis.title         = element_text(size = 7.5),
          axis.text          = element_text(size = 7))
}

p_neut <- make_sensitivity_plot("Neut", "Neutrophils",
                                 as.numeric(NEUT_THRESHOLDS), NULL)
p_plt  <- make_sensitivity_plot("Plt",  "Platelets",
                                 as.numeric(PLT_THRESHOLDS),  NULL)

pdf("results_HUMAN/sensitivity_baselines.pdf", width = 14, height = 6)
grid.arrange(
  p_neut, p_plt, ncol = 2,
  top = textGrob(
    sprintf("Sensibilité aux baselines — AUC=5, GFR=78, Q21D×2\nNoir = référence nominale | Couleur = min/med/max | Style = paramètre varié"),
    gp = gpar(fontface = "bold", fontsize = 10)
  )
)
dev.off()

cat("✓ Figure sauvegardée : results_HUMAN/sensitivity_baselines.pdf\n")
