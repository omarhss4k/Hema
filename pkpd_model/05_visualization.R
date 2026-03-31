# =============================================================================
# 05_visualization.R
# Fonctions de visualisation ggplot2 — Modèle PKPD Carboplatine (Fornari 2019)
#
# Plots disponibles :
#   1. plot_population_dynamics()  — profils médiane + IC90 par lignée
#   2. plot_pk_profile()           — concentration plasmatique PK
#   3. plot_damage_profile()       — dommage ADN au cours du temps
#   4. plot_grade_distribution()   — diagramme en barres des grades NCI-CTCAE
#   5. plot_nadir_boxplot()        — boîtes à moustaches des nadirs
#   6. plot_individual_vs_pop()    — superposition individu + population
#   7. save_all_plots()            — sauvegarde groupée PNG/PDF
#
# Dépendances :
#   - ggplot2, dplyr, tidyr, scales, patchwork (CRAN)
#   - 03_simulation_core.R (pour summarize_population)
#   - 04_scoring_toxicity.R (pour grade_distribution, nadir_summary)
# =============================================================================

# --- Chargement des dépendances (auto-installation si absent) ---
.required_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork")
for (pkg in .required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installation de '%s'...", pkg))
    install.packages(pkg, repos = getOption("repos"),
                     lib = .libPaths()[1], quiet = TRUE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, lib.loc = .libPaths())
  )
}

# =============================================================================
# THÈME GLOBAL
# =============================================================================

.pkpd_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background  = element_rect(fill = "#2C3E50", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold", size = base_size - 1),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey90"),
      axis.title        = element_text(size = base_size),
      axis.text         = element_text(size = base_size - 1),
      legend.position   = "bottom",
      legend.key.size   = unit(0.5, "cm"),
      plot.title        = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle     = element_text(size = base_size, hjust = 0, colour = "grey40")
    )
}

# Labels lisibles pour les variables d'état
.var_labels <- c(
  Neut   = "Neutrophiles (10⁹/L)",
  Mono   = "Monocytes (10⁹/L)",
  Plt    = "Plaquettes (10⁹/L)",
  Ret    = "Réticulocytes (10⁹/L)",
  RBC    = "GR (10⁹/L)",
  Damage = "Dommage ADN (u.a.)",
  Cen    = "Platine libre central (µmol)",
  Per    = "Platine libre péri. (µmol)",
  MPP    = "MPP (10⁹/L)",
  CMP    = "CMP (10⁹/L)",
  MEP    = "MEP (10⁹/L)"
)

# Couleurs par lignée
.lineage_colors <- c(
  Neut   = "#2166AC",
  Mono   = "#4DAC26",
  Plt    = "#D01C8B",
  Ret    = "#F1A340",
  RBC    = "#B2182B",
  Damage = "#762A83",
  MPP    = "#35978F",
  CMP    = "#BF812D",
  MEP    = "#80CDC1"
)

# =============================================================================
# SECTION 1 — PROFILS POPULATIONNELS (médiane + IC 5-95%)
# =============================================================================

#' Trace les profils temporels d'une population (médiane ± IC90)
#'
#' @param pop_summary  data.frame issu de summarize_population()
#' @param variables    Sous-ensemble de variables à tracer
#'                     (NULL = toutes les variables présentes)
#' @param thresholds   Liste nommée de seuils horizontaux (ex: list(Neut=0.5))
#'                     Ces seuils sont tracés en rouge pointillé (grade 4)
#' @param ncol         Nombre de colonnes dans le facet_wrap
#' @param title        Titre du graphique
#' @param dose_times_d Vecteur de temps de dose (jours) pour lignes verticales
#'
#' @return Objet ggplot2
plot_population_dynamics <- function(pop_summary,
                                     variables    = NULL,
                                     thresholds   = list(
                                       Neut = 0.5,   # Grade 4
                                       Plt  = 25     # Grade 4
                                     ),
                                     ncol          = 2,
                                     title         = "Dynamique hématologique — Population (N=1000)",
                                     dose_times_d  = 0) {

  df <- pop_summary
  if (!is.null(variables)) df <- df[df$variable %in% variables, ]

  # Labels lisibles
  df$var_label <- .var_labels[df$variable]
  df$var_label[is.na(df$var_label)] <- df$variable[is.na(df$var_label)]
  df$var_label <- factor(df$var_label)

  # Couleurs
  df$color_key <- .lineage_colors[df$variable]
  df$color_key[is.na(df$color_key)] <- "steelblue"

  # Seuils de toxicité
  threshold_df <- data.frame(
    variable  = names(thresholds),
    threshold = unlist(thresholds),
    stringsAsFactors = FALSE
  )
  threshold_df$var_label <- .var_labels[threshold_df$variable]
  threshold_df$var_label[is.na(threshold_df$var_label)] <- threshold_df$variable[is.na(threshold_df$var_label)]

  p <- ggplot(df, aes(x = time_d)) +
    # Ruban IC 5-95%
    geom_ribbon(aes(ymin = p05, ymax = p95, fill = variable), alpha = 0.20) +
    # Ruban IQR 25-75%
    geom_ribbon(aes(ymin = p25, ymax = p75, fill = variable), alpha = 0.35) +
    # Médiane
    geom_line(aes(y = p50, colour = variable), linewidth = 0.9) +
    # Lignes de dose
    geom_vline(xintercept = dose_times_d, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    # Seuils de toxicité (grade 4)
    geom_hline(data = threshold_df,
               aes(yintercept = threshold),
               linetype = "dotted", colour = "#CC0000", linewidth = 0.7) +
    facet_wrap(~ var_label, scales = "free_y", ncol = ncol) +
    scale_colour_manual(values = .lineage_colors, guide = "none") +
    scale_fill_manual(values   = .lineage_colors, guide = "none") +
    labs(
      title    = title,
      subtitle = "Médiane (trait) | IQR 25-75% (ombre foncée) | IC 5-95% (ombre claire) | seuil Grade 4 (pointillé rouge)",
      x        = "Temps (jours)",
      y        = "Concentration cellulaire"
    ) +
    .pkpd_theme()

  p
}

# =============================================================================
# SECTION 2 — PROFIL PK INDIVIDUEL OU POPULATIONNEL
# =============================================================================

#' Trace le profil PK (concentration plasmatique de platine libre)
#'
#' @param sim_df    data.frame individuel (simulate_patient()) OU
#'                 pop_summary filtré sur variable == "Cen"
#' @param mode      "individual" ou "population"
#' @param title     Titre
#'
#' @return ggplot2
plot_pk_profile <- function(sim_df, mode = "individual",
                            title = "Profil PK — Platine libre plasmatique") {

  if (mode == "individual") {
    # Concentration plasmatique µM = Cen / VCen (VCen stocké dans p, non dispo ici)
    # On trace directement Cen (µmol dans compartiment central)
    p <- ggplot(sim_df, aes(x = time_d, y = Cen)) +
      geom_line(colour = "#2166AC", linewidth = 1) +
      labs(title = title, x = "Temps (jours)",
           y = "Carboplatine — compartiment central (µmol)") +
      .pkpd_theme()
  } else {
    p <- ggplot(sim_df, aes(x = time_d)) +
      geom_ribbon(aes(ymin = p05, ymax = p95), fill = "#2166AC", alpha = 0.2) +
      geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#2166AC", alpha = 0.35) +
      geom_line(aes(y = p50), colour = "#2166AC", linewidth = 1) +
      labs(title = title, x = "Temps (jours)",
           y = "Carboplatine — compartiment central (µmol)") +
      .pkpd_theme()
  }
  p
}

# =============================================================================
# SECTION 3 — PROFIL DE DOMMAGE ADN
# =============================================================================

#' Trace le profil de dommage ADN
#'
#' @param sim_df    data.frame individuel
#'
#' @return ggplot2
plot_damage_profile <- function(sim_df,
                                title = "Dommage ADN au cours du temps") {
  ggplot(sim_df, aes(x = time_d, y = Damage)) +
    geom_area(fill = "#762A83", alpha = 0.3) +
    geom_line(colour = "#762A83", linewidth = 1) +
    labs(title = title, x = "Temps (jours)",
         y = "Dommage ADN (u.a.)") +
    .pkpd_theme()
}

# =============================================================================
# SECTION 4 — DISTRIBUTION DES GRADES NCI-CTCAE
# =============================================================================

#' Diagramme en barres empilées de la distribution des grades
#'
#' @param grade_dist_df  data.frame issu de grade_distribution()
#' @param title          Titre
#'
#' @return ggplot2
plot_grade_distribution <- function(grade_dist_df,
                                    title = "Distribution des grades NCI-CTCAE au nadir") {

  # Palettes NCI-CTCAE standards
  grade_colors <- c(
    "0" = "#4CAF50",   # vert
    "1" = "#CDDC39",   # vert-jaune
    "2" = "#FFC107",   # orange
    "3" = "#FF5722",   # rouge-orange
    "4" = "#B71C1C"    # rouge foncé
  )

  df <- grade_dist_df
  df$grade_f <- factor(df$grade, levels = 0:4,
                       labels = paste0("Grade ", 0:4))

  # Labels de lignée
  df$lineage_label <- .var_labels[df$lineage]
  df$lineage_label[is.na(df$lineage_label)] <- df$lineage[is.na(df$lineage_label)]

  ggplot(df, aes(x = lineage_label, y = pct, fill = grade_f)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct >= 3, paste0(pct, "%"), "")),
              position = position_stack(vjust = 0.5),
              size = 3.2, colour = "white", fontface = "bold") +
    scale_fill_manual(values = grade_colors, name = "Grade NCI-CTCAE") +
    scale_y_continuous(labels = percent_format(scale = 1),
                       limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title    = title,
      subtitle = "Pourcentage de patients par grade au nadir (N=1000)",
      x        = NULL,
      y        = "Patients (%)"
    ) +
    .pkpd_theme() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# =============================================================================
# SECTION 5 — BOÎTES À MOUSTACHES DES NADIRS
# =============================================================================

#' Boîtes à moustaches des valeurs de nadir par lignée
#'
#' @param scores_df  data.frame issu de score_population()
#' @param title      Titre
#'
#' @return ggplot2
plot_nadir_boxplot <- function(scores_df,
                               title = "Distribution des nadirs par lignée") {

  df <- scores_df
  df$lineage_label <- .var_labels[df$lineage]
  df$lineage_label[is.na(df$lineage_label)] <- df$lineage[is.na(df$lineage_label)]

  ggplot(df, aes(x = lineage_label, y = nadir_value, fill = lineage)) +
    geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.8, linewidth = 0.5) +
    scale_fill_manual(values = .lineage_colors, guide = "none") +
    labs(
      title    = title,
      subtitle = "Boîtes : IQR ; moustaches : 1.5×IQR ; points : valeurs extrêmes",
      x        = NULL,
      y        = "Valeur au nadir (10⁹/L)"
    ) +
    facet_wrap(~ lineage_label, scales = "free_y", ncol = 2) +
    .pkpd_theme() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# =============================================================================
# SECTION 6 — SUPERPOSITION INDIVIDU + POPULATION
# =============================================================================

#' Superpose un profil individuel sur les bandes populationnelles
#'
#' @param pop_summary  data.frame issu de summarize_population()
#' @param ind_sim      data.frame individuel (simulate_patient())
#' @param variable     Variable à tracer (ex: "Plt")
#' @param title        Titre
#'
#' @return ggplot2
plot_individual_vs_pop <- function(pop_summary, ind_sim,
                                   variable = "Plt",
                                   title    = NULL) {

  pop <- pop_summary[pop_summary$variable == variable, ]
  ind <- ind_sim[, c("time_d", variable)]
  names(ind)[2] <- "value"

  var_lbl <- .var_labels[variable] %||% variable
  col     <- .lineage_colors[variable] %||% "steelblue"
  if (is.null(title)) title <- paste0("Profil individuel vs population — ", var_lbl)

  ggplot(pop, aes(x = time_d)) +
    geom_ribbon(aes(ymin = p05, ymax = p95), fill = col, alpha = 0.15) +
    geom_ribbon(aes(ymin = p25, ymax = p75), fill = col, alpha = 0.30) +
    geom_line(aes(y = p50), colour = col, linewidth = 0.8, linetype = "dashed") +
    geom_line(data = ind, aes(x = time_d, y = value),
              colour = "black", linewidth = 1.1) +
    labs(
      title    = title,
      subtitle = "Noir = individu | Pointillé = médiane population | Ombres = IC25-75% et IC5-95%",
      x        = "Temps (jours)",
      y        = var_lbl
    ) +
    .pkpd_theme()
}

# =============================================================================
# SECTION 7 — SAUVEGARDE GROUPÉE
# =============================================================================

#' Sauvegarde tous les plots principaux en PNG et PDF
#'
#' @param pop_result   Résultat de simulate_population()
#' @param scores_df    Résultat de score_population()
#' @param grade_dist   Résultat de grade_distribution()
#' @param pop_summary  Résultat de summarize_population()
#' @param out_dir      Répertoire de sortie
#' @param prefix       Préfixe des fichiers (défaut "pkpd")
#' @param width        Largeur (pouces)
#' @param height       Hauteur (pouces)
#'
#' @return Invisible : liste des chemins générés
save_all_plots <- function(pop_result,
                           scores_df,
                           grade_dist,
                           pop_summary,
                           out_dir = "results_pkpd",
                           prefix  = "pkpd",
                           width   = 12,
                           height  = 9) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- list()

  message("Génération des graphiques...")

  # 1. Dynamique populationnelle — lignées hématologiques
  p1 <- plot_population_dynamics(
    pop_summary[pop_summary$variable %in% c("Neut", "Plt", "Ret", "RBC"), ]
  )
  paths$pop_dynamics <- file.path(out_dir, paste0(prefix, "_01_population_dynamics.png"))
  ggsave(paths$pop_dynamics, p1, width = width, height = height, dpi = 150)

  # 2. Dynamique populationnelle — progéniteurs
  if (any(c("MPP","CMP","MEP") %in% pop_summary$variable)) {
    p2 <- plot_population_dynamics(
      pop_summary[pop_summary$variable %in% c("MPP","CMP","MEP","Damage"), ],
      title = "Dynamique des progéniteurs et dommage ADN",
      ncol  = 2
    )
    paths$progenitors <- file.path(out_dir, paste0(prefix, "_02_progenitors.png"))
    ggsave(paths$progenitors, p2, width = width, height = height, dpi = 150)
  }

  # 3. Distribution des grades
  p3 <- plot_grade_distribution(grade_dist)
  paths$grades <- file.path(out_dir, paste0(prefix, "_03_grade_distribution.png"))
  ggsave(paths$grades, p3, width = 10, height = 6, dpi = 150)

  # 4. Nadirs
  p4 <- plot_nadir_boxplot(scores_df)
  paths$nadirs <- file.path(out_dir, paste0(prefix, "_04_nadir_boxplot.png"))
  ggsave(paths$nadirs, p4, width = 10, height = 7, dpi = 150)

  message(sprintf("Fichiers sauvegardés dans '%s/'", out_dir))
  invisible(paths)
}
