############################################################
# vpc.R
# Visual Predictive Check — Fornari 2019
#
# Améliorations vs version originale :
#   1. Zones de grade NCI-CTCAE superposées sur les panneaux VPC
#   2. Grades calculés dans la même boucle que le VPC (pas de double simulation)
#   3. IIV étendue : Slope + baselines + MTT_Neut/Plt + kcirc_Neut/Plt
#   4. Figure combinée 3 panneaux : Neut VPC | Plt VPC | barplot grades
############################################################
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# ── Valeurs sigma (Table S4 — Fornari 2019) ─────────────
SIGMA <- c(
  MPP  = 0.33,
  CMP  = 0.19,
  MEP  = 0.33,
  Neut = 0.17,
  Mono = 0.17,
  Plt  = 0.17,
  Ret  = 0.36,
  RBC  = 0.06
)

# ── Seuils NCI-CTCAE v5.0 ───────────────────────────────
NEUT_THRESHOLDS <- c(G1 = 2.0, G2 = 1.5, G3 = 1.0, G4 = 0.5)
PLT_THRESHOLDS  <- c(G1 = 150, G2 = 75,  G3 = 50,  G4 = 25)
GRADE_COLORS    <- c(G1 = "#fee08b", G2 = "#fc8d59",
                     G3 = "#d73027", G4 = "#7b0404")

# ── Attribution de grade ─────────────────────────────────
assign_grade_neut <- function(x) {
  ifelse(x <  NEUT_THRESHOLDS["G4"], 4L,
  ifelse(x <  NEUT_THRESHOLDS["G3"], 3L,
  ifelse(x <  NEUT_THRESHOLDS["G2"], 2L,
  ifelse(x <  NEUT_THRESHOLDS["G1"], 1L, 0L))))
}
assign_grade_plt <- function(x) {
  ifelse(x <  PLT_THRESHOLDS["G4"], 4L,
  ifelse(x <  PLT_THRESHOLDS["G3"], 3L,
  ifelse(x <  PLT_THRESHOLDS["G2"], 2L,
  ifelse(x <  PLT_THRESHOLDS["G1"], 1L, 0L))))
}

# ════════════════════════════════════════════════════════
# 1. GÉNÉRER LES SIMULATIONS VPC (erreur résiduelle seule)
#    Utilisé pour le VPC rat (Figure 3)
# ════════════════════════════════════════════════════════
run_vpc <- function(sim, n_sim = 1000, seed = 42) {

  set.seed(seed)
  cells <- names(SIGMA)
  n_t   <- nrow(sim)

  vpc_stats <- lapply(cells, function(cell) {
    sigma  <- SIGMA[cell]
    Y_pred <- sim[[cell]]
    eps    <- matrix(rnorm(n_sim * n_t, 0, sigma), nrow = n_sim, ncol = n_t)
    Y_mat  <- sweep(exp(eps), 2, Y_pred, FUN = "*")
    data.frame(
      days = sim$days,
      p05  = apply(Y_mat, 2, quantile, 0.05, na.rm = TRUE),
      p50  = apply(Y_mat, 2, quantile, 0.50, na.rm = TRUE),
      p95  = apply(Y_mat, 2, quantile, 0.95, na.rm = TRUE),
      pred = Y_pred,
      cell = cell
    )
  })
  names(vpc_stats) <- cells
  vpc_stats
}

# ════════════════════════════════════════════════════════
# 2. PANNEAU VPC INDIVIDUEL
#    Nouveau : zones de grade NCI-CTCAE en fond coloré
# ════════════════════════════════════════════════════════
plot_vpc_panel <- function(vpc_df, title, color = "steelblue",
                           baseline   = NULL,
                           dose_days  = NULL,
                           obs        = NULL,
                           thresholds = NULL) {  # c(G1, G2, G3, G4) décroissant

  all_vals <- c(vpc_df$p05, vpc_df$p95)
  if (!is.null(obs) && length(obs$value) > 0)
    all_vals <- c(all_vals, obs$value[obs$value > 0])

  all_pos <- all_vals[is.finite(all_vals) & all_vals > 0]
  ymin     <- 10^floor(log10(min(all_pos) * 0.3))
  ymax     <- 10^ceiling(log10(max(all_pos) * 3.0))
  breaks_y <- 10^seq(log10(ymin), log10(ymax), by = 1)

  vpc_df$p05  <- pmax(vpc_df$p05,  ymin * 0.5)
  vpc_df$p95  <- pmax(vpc_df$p95,  ymin * 0.5)
  vpc_df$p50  <- pmax(vpc_df$p50,  ymin * 0.5)
  vpc_df$pred <- pmax(vpc_df$pred, ymin * 0.5)

  p <- ggplot(vpc_df, aes(x = days))

  # ── Zones de grade (fond coloré, avant le ruban) ──────
  if (!is.null(thresholds) && length(thresholds) == 4) {
    thr <- as.numeric(thresholds)
    p <- p +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = ymin * 0.5,  ymax = thr[4],
               fill = GRADE_COLORS["G4"], alpha = 0.10) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = thr[4], ymax = thr[3],
               fill = GRADE_COLORS["G3"], alpha = 0.10) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = thr[3], ymax = thr[2],
               fill = GRADE_COLORS["G2"], alpha = 0.10) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = thr[2], ymax = thr[1],
               fill = GRADE_COLORS["G1"], alpha = 0.08) +
      # Lignes de seuil
      geom_hline(yintercept = thr,
                 linetype   = "dashed",
                 linewidth  = 0.35,
                 colour     = unname(GRADE_COLORS)) +
      # Labels G1–G4
      annotate("text",
               x      = max(vpc_df$days) * 0.97,
               y      = thr * 1.18,
               hjust  = 1, size = 2.5,
               label  = names(thresholds),
               colour = unname(GRADE_COLORS))
  }

  # ── Ruban + courbes ─────────────────────────────────────
  p <- p +
    geom_ribbon(aes(ymin = p05, ymax = p95),
                fill = color, alpha = 0.20) +
    geom_line(aes(y = p50),
              color = "red", linewidth = 0.9, alpha = 0.9) +
    geom_line(aes(y = pred),
              color = "black", linewidth = 0.4,
              linetype = "dashed", alpha = 0.5)

  # ── Points observés ─────────────────────────────────────
  if (!is.null(obs) && length(obs$time) > 0) {
    df_obs <- data.frame(x = obs$time, y = pmax(obs$value, ymin * 0.5))
    p <- p + geom_point(data = df_obs, aes(x = x, y = y),
                        shape = 21, fill = "white", color = "black",
                        size = 1.8, stroke = 0.65, alpha = 0.85)
  }

  # ── Baseline et doses ────────────────────────────────────
  if (!is.null(baseline) && baseline > 0)
    p <- p + geom_hline(yintercept = baseline,
                        linetype = "dashed", color = "grey50",
                        linewidth = 0.45, alpha = 0.6)

  if (!is.null(dose_days)) {
    vd <- dose_days[dose_days <= max(vpc_df$days)]
    if (length(vd) > 0)
      p <- p + geom_vline(xintercept = vd, linetype = "dotted",
                          color = "grey55", linewidth = 0.35, alpha = 0.55)
  }

  p +
    scale_y_log10(limits  = c(ymin, ymax),
                  breaks  = breaks_y,
                  labels  = trans_format("log10", math_format(10^.x))) +
    labs(title = title,
         x     = "Time (d)",
         y     = expression(10^9~cells~L^{-1})) +
    theme_bw(base_size = 9.5) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      plot.title  = element_text(face = "bold", size = 9, hjust = 0.5),
      axis.title  = element_text(size = 7.5),
      axis.text   = element_text(size = 7)
    )
}

# ════════════════════════════════════════════════════════
# 3. GRILLE COMPLÈTE 4×2 VPC (rat / validation)
# ════════════════════════════════════════════════════════
plot_vpc_all <- function(vpc_stats, pars,
                         titre     = "VPC — Haematopoiesis",
                         dose_days = NULL,
                         obs_list  = NULL) {

  get_obs <- function(name) {
    if (!is.null(obs_list) && name %in% names(obs_list)) obs_list[[name]]
    else NULL
  }

  colors    <- c(MPP = "#2166ac", CMP = "#4dac26", MEP = "#7b3294",
                 Neut = "#d7191c", Mono = "#a6611a", Plt = "#1b9e77",
                 Ret = "#e377c2", RBC = "#cc3399")
  baselines <- list(MPP = pars$MPP0, CMP = pars$CMP0, MEP = pars$MEP0,
                    Neut = pars$Neut0, Mono = pars$Mono0, Plt = pars$Plt0,
                    Ret = pars$Ret0, RBC = pars$RBC0)
  titles    <- c(MPP = "Multi-potent progenitors", CMP = "Common myeloid progenitors",
                 MEP = "MEP", Neut = "Neutrophils", Mono = "Monocytes",
                 Plt = "Platelets", Ret = "Reticulocytes", RBC = "Red blood cells")
  order     <- c("MPP", "Neut", "CMP", "Mono", "MEP", "Plt", "Ret", "RBC")

  panels <- lapply(order, function(cell) {
    thr <- switch(cell, Neut = NEUT_THRESHOLDS, Plt = PLT_THRESHOLDS, NULL)
    plot_vpc_panel(vpc_stats[[cell]], titles[cell], colors[cell],
                   baselines[[cell]], dose_days, get_obs(cell),
                   thresholds = thr)
  })

  legend_grob <- textGrob(
    paste0("Ruban = [5e–95e percentile]  |  Rouge = médiane  |  ",
           "Tiret = prédiction déterministe  |  ● = obs  |  ",
           "Zones = grades NCI-CTCAE"),
    gp = gpar(fontsize = 7.5, col = "grey30")
  )
  grid.arrange(grobs  = panels, ncol = 2,
               top    = textGrob(titre, gp = gpar(fontface = "bold", fontsize = 11)),
               bottom = legend_grob)
}

# ════════════════════════════════════════════════════════
# 4. SAUVEGARDE PDF (rat)
# ════════════════════════════════════════════════════════
save_vpc <- function(sim, pars, file, titre, n_sim = 1000,
                     dose_days = NULL, obs_list = NULL,
                     width = 12, height = 11) {
  cat(sprintf("  → VPC (%d simulations)...\n", n_sim))
  vpc_stats <- run_vpc(sim, n_sim = n_sim)
  pdf(file, width = width, height = height)
  plot_vpc_all(vpc_stats, pars, titre = titre,
               dose_days = dose_days, obs_list = obs_list)
  dev.off()
  message("✓ VPC sauvegardé : ", file)
  invisible(vpc_stats)
}

# ════════════════════════════════════════════════════════
# 5. BARPLOT DISTRIBUTION DES GRADES (panneau interne)
# ════════════════════════════════════════════════════════
.plot_grade_bar <- function(grade_neut, grade_plt, n_patients,
                            auc_target, interval_h, n_cycles) {

  pct_neut <- sapply(1:4, function(g) 100 * mean(grade_neut == g, na.rm = TRUE))
  pct_plt  <- sapply(1:4, function(g) 100 * mean(grade_plt  == g, na.rm = TRUE))

  df <- data.frame(
    Grade   = rep(paste0("G", 1:4), 2),
    Pct     = c(pct_neut, pct_plt),
    Lineage = rep(c("Neutropenia", "Thrombocytopenia"), each = 4)
  )
  df$Grade   <- factor(df$Grade,   levels = paste0("G", 1:4))
  df$Lineage <- factor(df$Lineage, levels = c("Neutropenia", "Thrombocytopenia"))

  ggplot(df, aes(x = Grade, y = Pct, fill = Grade)) +
    geom_col(width = 0.65, alpha = 0.88) +
    geom_text(aes(label = sprintf("%.1f%%", Pct)),
              vjust = -0.4, size = 3.0) +
    facet_wrap(~ Lineage) +
    scale_fill_manual(values = unname(GRADE_COLORS), guide = "none") +
    scale_y_continuous(limits = c(0, max(df$Pct, 1) * 1.30),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = "Grades NCI-CTCAE au nadir",
      subtitle = sprintf(
        "AUC=%g  Q%dD\u00d7%d  n=%d  IIV\u00a0: Slope\u00a0+\u00a0MTT\u00a0+\u00a0kcirc\u00a0+\u00a0baselines",
        auc_target, round(interval_h / 24), n_cycles, n_patients),
      x = NULL, y = "% patients"
    ) +
    theme_bw(base_size = 9.5) +
    theme(
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 9),
      plot.subtitle    = element_text(hjust = 0.5, size = 7, color = "grey40"),
      panel.grid.minor = element_blank()
    )
}

# ════════════════════════════════════════════════════════
# 6. VPC HUMAIN — Figure combinée (Neut VPC + Plt VPC + Grades)
#
#  Améliorations vs version originale :
#    • Zones de grade NCI-CTCAE sur les panneaux Neut et Plt
#    • Grades calculés dans la même boucle (pas de double simulation)
#    • IIV étendue :
#        – Slope_MPP/CMP/MEP (Table S4)
#        – Neut0, Plt0       (Table S4)
#        – MTT_Neut, MTT_Plt (CV rat Table 1)
#        – kcirc_Neut, kcirc_Plt (CV rat Table 1)
#    • Figure 3 panneaux : [Neut | Plt] en haut, [Grades] en bas
#    • Retourne invisiblement les matrices et vecteurs de grade
# ════════════════════════════════════════════════════════
save_vpc_human <- function(sim, pars, file, titre,
                           n_sim      = 1000,
                           dose_days  = NULL,
                           obs_list   = NULL,
                           auc_target = 5,
                           gfr_fixed  = 125,
                           times      = NULL,
                           interval_h = 21 * 24,
                           n_cycles   = 2,
                           width      = 10,
                           height     = 11) {

  cat(sprintf(
    "  → VPC humain + grades : %d patients (GFR=%g fixe, IIV étendue)...\n",
    n_sim, gfr_fixed))
  set.seed(42)

  if (is.null(times)) times <- seq(0, 63 * 24, by = 1)
  n_t       <- length(times)
  dose_fixe <- auc_target * (gfr_fixed + 25)

  # ── Paramètres IIV ──────────────────────────────────────
  # Slopes & baselines : sigma log-additif (Table S4 Fornari 2019)
  omega_MPP   <- SIGMA["MPP"]
  omega_CMP   <- SIGMA["CMP"]
  omega_MEP   <- SIGMA["MEP"]
  omega_Neut0 <- SIGMA["Neut"]
  omega_Plt0  <- SIGMA["Plt"]

  # MTT : CV rat Table 1 → omega = sqrt(log(1 + CV²))
  omega_MTTNeut <- sqrt(log(1 + 0.56^2))   # CV = 56 %
  omega_MTTPlt  <- sqrt(log(1 + 0.056^2))  # CV =  6 %

  # kcirc : CV rat Table 1
  omega_kcircNeut <- sqrt(log(1 + 0.45^2)) # CV = 45 %
  omega_kcircPlt  <- sqrt(log(1 + 0.46^2)) # CV = 46 %

  # ── Matrices de résultats ────────────────────────────────
  mat_Neut   <- matrix(NA_real_, nrow = n_sim, ncol = n_t)
  mat_Plt    <- matrix(NA_real_, nrow = n_sim, ncol = n_t)
  grade_neut <- integer(n_sim)
  grade_plt  <- integer(n_sim)

  for (i in seq_len(n_sim)) {

    pars_i <- pars

    # Slopes (sensibilité médicament)
    pars_i$Slope_MPP <- pars$Slope_MPP * exp(rnorm(1, 0, omega_MPP))
    pars_i$Slope_CMP <- pars$Slope_CMP * exp(rnorm(1, 0, omega_CMP))
    pars_i$Slope_MEP <- pars$Slope_MEP * exp(rnorm(1, 0, omega_MEP))

    # Baselines individuelles
    neut0_i        <- pars$Neut0 * exp(rnorm(1, 0, omega_Neut0))
    plt0_i         <- pars$Plt0  * exp(rnorm(1, 0, omega_Plt0))
    pars_i$Neut0   <- neut0_i
    pars_i$Plt0    <- plt0_i

    # MTT individuels
    pars_i$MTT_Neut <- pars$MTT_Neut * exp(rnorm(1, 0, omega_MTTNeut))
    pars_i$MTT_Plt  <- pars$MTT_Plt  * exp(rnorm(1, 0, omega_MTTPlt))

    # kcirc individuels
    pars_i$k_circ_Neut <- pars$k_circ_Neut * exp(rnorm(1, 0, omega_kcircNeut))
    pars_i$k_circ_Plt  <- pars$k_circ_Plt  * exp(rnorm(1, 0, omega_kcircPlt))

    # Schéma posologique (PK identique pour tous — Supp. S11/S12)
    pars_i$rate_fun <- make_repeated_infusion(
      dose_mg    = dose_fixe,
      Tinfu_h    = 1,
      interval_h = interval_h,
      n_cycles   = n_cycles
    )

    # État initial rééquilibré aux baselines individuelles (Eq. S3)
    a_Neut  <- 3 / pars_i$MTT_Neut
    a_Plt   <- 3 / pars_i$MTT_Plt
    T_Neut  <- pars_i$k_circ_Neut * neut0_i / a_Neut
    T_Plt   <- pars_i$k_circ_Plt  * plt0_i  / a_Plt
    T1_Plt  <- T_Plt / pars$lambda2

    state_i <- init_state
    state_i["Neut"]    <- neut0_i
    state_i["Plt"]     <- plt0_i
    state_i["T1_Neut"] <- T_Neut
    state_i["T2_Neut"] <- T_Neut
    state_i["T3_Neut"] <- T_Neut
    state_i["T1_Plt"]  <- T1_Plt
    state_i["T2_Plt"]  <- T_Plt
    state_i["T3_Plt"]  <- T_Plt

    tryCatch({
      out_i <- as.data.frame(lsoda(
        y        = state_i,
        times    = times,
        func     = pkpd_fornari,
        parms    = pars_i,
        rtol     = 1e-4,
        atol     = 1e-6,
        maxsteps = 10000
      ))
      mat_Neut[i, ]  <- out_i$Neut
      mat_Plt[i, ]   <- out_i$Plt
      grade_neut[i]  <- assign_grade_neut(min(out_i$Neut, na.rm = TRUE))
      grade_plt[i]   <- assign_grade_plt( min(out_i$Plt,  na.rm = TRUE))
    }, error = function(e) NULL)

    if (i %% 200 == 0)
      cat(sprintf("    %d/%d patients simulés\n", i, n_sim))
  }

  # ── Percentiles VPC ──────────────────────────────────────
  days_vec   <- times / 24
  make_stats <- function(mat, Yref) data.frame(
    days = days_vec,
    p05  = apply(mat, 2, quantile, 0.05, na.rm = TRUE),
    p50  = apply(mat, 2, quantile, 0.50, na.rm = TRUE),
    p95  = apply(mat, 2, quantile, 0.95, na.rm = TRUE),
    pred = Yref
  )
  vpc_neut <- make_stats(mat_Neut, sim$Neut)
  vpc_plt  <- make_stats(mat_Plt,  sim$Plt)

  # ── Résumé console ───────────────────────────────────────
  pct_n <- sapply(1:4, function(g) 100 * mean(grade_neut == g, na.rm = TRUE))
  pct_p <- sapply(1:4, function(g) 100 * mean(grade_plt  == g, na.rm = TRUE))
  cat(sprintf("  Neutropénie  : G1=%.1f%%  G2=%.1f%%  G3=%.1f%%  G4=%.1f%%\n",
              pct_n[1], pct_n[2], pct_n[3], pct_n[4]))
  cat(sprintf("  Thrombopénie : G1=%.1f%%  G2=%.1f%%  G3=%.1f%%  G4=%.1f%%\n",
              pct_p[1], pct_p[2], pct_p[3], pct_p[4]))

  # ── Panneaux VPC ─────────────────────────────────────────
  get_obs <- function(nm) {
    if (!is.null(obs_list) && nm %in% names(obs_list)) obs_list[[nm]]
    else NULL
  }

  p_neut <- plot_vpc_panel(
    vpc_neut, "Neutrophils", "#d7191c",
    baseline   = pars$Neut0,
    dose_days  = dose_days,
    obs        = get_obs("Neut"),
    thresholds = NEUT_THRESHOLDS
  )
  p_plt <- plot_vpc_panel(
    vpc_plt, "Platelets", "#1b9e77",
    baseline   = pars$Plt0,
    dose_days  = dose_days,
    obs        = get_obs("Plt"),
    thresholds = PLT_THRESHOLDS
  )

  # ── Panneau grades ───────────────────────────────────────
  p_grades <- .plot_grade_bar(grade_neut, grade_plt, n_sim,
                               auc_target, interval_h, n_cycles)

  # ── Assemblage figure ────────────────────────────────────
  legend_grob <- textGrob(
    paste0(
      "Ruban = [5e\u201395e percentile]  |  Rouge = m\u00e9diane simul\u00e9e  |  ",
      "Tiret = pr\u00e9diction d\u00e9terministe  |  \u25cf = donn\u00e9es observ\u00e9es  |  ",
      "Zones color\u00e9es = grades NCI-CTCAE"
    ),
    gp = gpar(fontsize = 7, col = "grey30")
  )

  pdf(file, width = width, height = height)
  grid.arrange(
    p_neut, p_plt, p_grades,
    layout_matrix = matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
    top    = textGrob(titre, gp = gpar(fontface = "bold", fontsize = 11)),
    bottom = legend_grob
  )
  dev.off()

  message("✓ VPC humain + grades sauvegardé : ", file)
  invisible(list(
    vpc        = list(Neut = vpc_neut, Plt = vpc_plt),
    grade_neut = grade_neut,
    grade_plt  = grade_plt
  ))
}
