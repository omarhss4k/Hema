############################################################
# 56_analyse_sensibilite.R
# Analyse de sensibilité — paramètres PK/PD → nadirs / grades
#
# Utilise les données existantes : results/53_population_tdxd_humain.rds
# (pas de re-simulation)
#
# Prérequis : 53_run_population_tdxd_humain.R
#             (→ results/53_population_tdxd_humain.rds)
############################################################

library(ggplot2)

# ── Chargement résultats modèle ───────────────────────────
rds_path <- "results/53_population_tdxd_humain.rds"

if (!file.exists(rds_path)) {
  stop(paste0(
    "Fichier manquant : ", rds_path, "\n",
    "Lancez d'abord 53_run_population_tdxd_humain.R"
  ))
}

results <- readRDS(rds_path)
n       <- nrow(results)
cat(sprintf("\n═══ Analyse de sensibilité — N = %d patients ═══\n\n", n))

# ── Paramètres d'entrée et sorties ────────────────────────
params <- c("Slope_CMP", "Slope_MEP", "CL_ADC", "V1_ADC", "Slope_Plt_direct")
nadirs <- c("Neut_nadir", "RBC_nadir", "Plt_nadir")

# Vérifier que toutes les colonnes existent
all_cols <- c(params, nadirs, "Grade_Neut")
missing_cols <- setdiff(all_cols, colnames(results))
if (length(missing_cols) > 0) {
  warning(paste("Colonnes manquantes dans 53_population_tdxd_humain.rds :",
                paste(missing_cols, collapse = ", ")))
  # Filtrer pour ne garder que les colonnes disponibles
  params <- intersect(params, colnames(results))
  nadirs <- intersect(nadirs, colnames(results))
}

# ── a) Corrélations de Spearman ───────────────────────────
cat("a) Corrélations de Spearman (paramètre → nadir) :\n")
cat("─────────────────────────────────────────────────────────────\n")

spearman_list <- list()

for (param in params) {
  row_vals <- numeric(length(nadirs))
  names(row_vals) <- nadirs
  for (nadir in nadirs) {
    ct <- cor.test(results[[param]], results[[nadir]],
                   method = "spearman", exact = FALSE)
    row_vals[nadir] <- ct$estimate
  }
  spearman_list[[param]] <- row_vals
}

spearman_df <- do.call(rbind, spearman_list)
spearman_df <- as.data.frame(spearman_df)

# Affichage console
header_fmt <- "  %-22s  %12s  %12s  %12s\n"
cat(sprintf(header_fmt, "Paramètre", "Neut_nadir", "RBC_nadir", "Plt_nadir"))
cat(sprintf("  %s\n", paste(rep("-", 62), collapse = "")))

for (param in params) {
  row <- spearman_df[param, ]
  neut_val  <- if ("Neut_nadir"  %in% names(row)) sprintf("%+.3f", row[["Neut_nadir"]])  else "  N/A"
  rbc_val   <- if ("RBC_nadir"   %in% names(row)) sprintf("%+.3f", row[["RBC_nadir"]])   else "  N/A"
  plt_val   <- if ("Plt_nadir"   %in% names(row)) sprintf("%+.3f", row[["Plt_nadir"]])   else "  N/A"
  cat(sprintf("  %-22s  %12s  %12s  %12s\n", param, neut_val, rbc_val, plt_val))
}

# ── b) Régression logistique P(G3-4 neutropénie) ─────────
cat("\nb) Régression logistique P(Grade G3-4 Neutropénie) :\n")
cat("─────────────────────────────────────────────────────────────\n")

if ("Grade_Neut" %in% colnames(results)) {
  results$G34_Neut <- as.integer(results$Grade_Neut %in% c("G3", "G4"))

  log_params <- intersect(c("Slope_CMP", "Slope_MEP", "CL_ADC", "V1_ADC"),
                          colnames(results))

  if (length(log_params) >= 1) {
    formula_str <- paste("G34_Neut ~", paste(log_params, collapse = " + "))
    logit_formula <- as.formula(formula_str)

    # Standardisation des prédicteurs pour stabilité numérique
    results_std <- results
    for (p in log_params) {
      sd_p <- sd(results[[p]], na.rm = TRUE)
      if (sd_p > 0) results_std[[p]] <- scale(results[[p]])[,1]
    }

    tryCatch({
      fit <- glm(logit_formula, data = results_std, family = binomial)
      cat(sprintf("  Formule : %s\n\n", formula_str))
      cat("  Odds Ratios (IC95%) :\n")

      coefs    <- coef(fit)
      ci_logit <- suppressMessages(confint(fit))
      or_df    <- data.frame(
        Parametre = names(coefs),
        OR        = exp(coefs),
        CI_lo     = exp(ci_logit[, 1]),
        CI_hi     = exp(ci_logit[, 2]),
        p_value   = coef(summary(fit))[, 4],
        stringsAsFactors = FALSE
      )

      or_fmt <- "  %-22s  OR = %6.3f  [%6.3f – %6.3f]  p = %.4f\n"
      for (i in seq_len(nrow(or_df))) {
        cat(sprintf(or_fmt,
                    or_df$Parametre[i],
                    or_df$OR[i],
                    or_df$CI_lo[i],
                    or_df$CI_hi[i],
                    or_df$p_value[i]))
      }
    }, error = function(e) {
      cat(sprintf("  Erreur régression logistique : %s\n", e$message))
    })
  } else {
    cat("  Paramètres insuffisants pour la régression logistique.\n")
  }
} else {
  cat("  Colonne Grade_Neut absente — régression ignorée.\n")
}

# ── c) Tornado plot — corrélations Neut_nadir ────────────
cat("\nc) Tornado plot — corrélations de Spearman avec Neut_nadir :\n")

if ("Neut_nadir" %in% nadirs) {
  neut_cors <- sapply(params, function(p) {
    spearman_df[p, "Neut_nadir"]
  })

  tornado_df <- data.frame(
    Parametre   = params,
    Spearman_r  = as.numeric(neut_cors),
    stringsAsFactors = FALSE
  )
  # Tri par valeur absolue décroissante
  tornado_df <- tornado_df[order(abs(tornado_df$Spearman_r)), ]
  tornado_df$Parametre <- factor(tornado_df$Parametre,
                                  levels = tornado_df$Parametre)
  tornado_df$Direction <- ifelse(tornado_df$Spearman_r < 0, "Négatif", "Positif")

  p_tornado <- ggplot(tornado_df,
                      aes(x = Spearman_r, y = Parametre, fill = Direction)) +

    geom_col(width = 0.6, color = "white", linewidth = 0.4) +

    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.7) +

    geom_text(aes(label  = sprintf("r = %+.3f", Spearman_r),
                  hjust  = ifelse(Spearman_r >= 0, -0.15, 1.15)),
              size = 3.6, fontface = "bold", color = "grey20") +

    scale_fill_manual(
      values = c("Négatif" = "#d73027", "Positif" = "#4575b4"),
      name   = "Direction"
    ) +

    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      labels = function(x) sprintf("%.2f", x)
    ) +

    labs(
      x        = "Corrélation de Spearman (r)",
      y        = NULL,
      title    = "Analyse de sensibilité — Paramètres vs Nadir neutrophilique",
      subtitle = sprintf(
        "Corrélations de Spearman entre paramètres PK/PD et Neut_nadir (N=%d)",
        n)
    ) +

    theme_classic(base_size = 13) +
    theme(
      axis.text.y      = element_text(size = 12, face = "bold"),
      axis.text.x      = element_text(size = 11),
      axis.title.x     = element_text(size = 12, face = "bold"),
      legend.position  = "bottom",
      panel.grid.major.x = element_line(color = "grey90"),
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "grey40")
    )

  # ── d) Export ─────────────────────────────────────────
  dir.create("results_PKPD_human", showWarnings = FALSE)

  out_pdf <- "results_PKPD_human/sensitivity_tornado.pdf"
  out_png <- "results_PKPD_human/sensitivity_tornado.png"

  ggsave(out_pdf, plot = p_tornado, width = 8, height = 5, device = cairo_pdf)
  ggsave(out_png, plot = p_tornado, width = 8, height = 5, dpi = 300)

  cat(sprintf("  -> %s\n", out_pdf))
  cat(sprintf("  -> %s\n", out_png))
} else {
  cat("  Neut_nadir non disponible — tornado plot ignoré.\n")
}

# ── e) Résumé console des paramètres les plus influents ──
cat("\ne) Résumé — Paramètres les plus influents sur Neut_nadir :\n")
cat("─────────────────────────────────────────────────────────────\n")

if ("Neut_nadir" %in% nadirs) {
  sorted_cors <- sort(abs(neut_cors), decreasing = TRUE)
  for (i in seq_along(sorted_cors)) {
    param_name <- names(sorted_cors)[i]
    r_val      <- neut_cors[param_name]
    direction  <- ifelse(r_val < 0, "↓ nadir (plus toxique)", "↑ nadir (protecteur)")
    cat(sprintf("  %d. %-22s  r = %+.3f  →  %s\n",
                i, param_name, r_val, direction))
  }

  top_param <- names(sorted_cors)[1]
  cat(sprintf(
    "\n  Paramètre le plus déterminant pour la neutropénie : %s (|r| = %.3f)\n\n",
    top_param, sorted_cors[1]))
} else {
  cat("  Neut_nadir non disponible.\n")
}
