############################################################
# vpc.R
# Visual Predictive Check — Fornari 2019
# Erreur log-additive (Table S4) sur 1000 simulations
# sigma : MPP=0.33, CMP=0.19, MEP=0.33, Neut=0.17,
#         Mono=0.17, Plt=0.17, Ret=0.36, RBC=0.06
############################################################
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# ── Valeurs sigma (Table S4) ────────────────────────────
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

# ════════════════════════════════════════════════════════
# 1. GÉNÉRER LES SIMULATIONS VPC
#    Entrée  : sim       = data.frame issu de simulate_all()
#              n_sim     = nombre de simulations (défaut 1000)
#    Sortie  : liste de data.frames de percentiles par cellule
# ════════════════════════════════════════════════════════
run_vpc <- function(sim, n_sim = 1000, seed = 42) {

  set.seed(seed)
  cells <- names(SIGMA)
  n_t   <- nrow(sim)

  # Pour chaque type cellulaire, générer n_sim profils bruités
  vpc_stats <- lapply(cells, function(cell) {

    sigma   <- SIGMA[cell]
    Y_pred  <- sim[[cell]]           # prédiction déterministe (vecteur)

    # Matrice n_sim × n_t de simulations
    # Y_sim[i,j] = Y_pred[j] × exp( ε[i,j] ),  ε ~ N(0, σ²)
    eps    <- matrix(rnorm(n_sim * n_t, mean = 0, sd = sigma),
                     nrow = n_sim, ncol = n_t)
    Y_mat  <- sweep(exp(eps), 2, Y_pred, FUN = "*")

    # Percentiles colonne par colonne
    data.frame(
      days = sim$days,
      p05  = apply(Y_mat, 2, quantile, probs = 0.05, na.rm = TRUE),
      p50  = apply(Y_mat, 2, quantile, probs = 0.50, na.rm = TRUE),
      p95  = apply(Y_mat, 2, quantile, probs = 0.95, na.rm = TRUE),
      pred = Y_pred,
      cell = cell
    )
  })

  names(vpc_stats) <- cells
  vpc_stats
}

# ════════════════════════════════════════════════════════
# 2. PANNEAU VPC INDIVIDUEL
#    Ruban bleu   = [5e, 95e] percentile simulé
#    Ligne rouge  = médiane simulée (p50)
#    Ligne noire  = prédiction déterministe
#    Points noirs = données observées (obs)
#    Tiret gris   = baseline
# ════════════════════════════════════════════════════════
plot_vpc_panel <- function(vpc_df, title, color = "steelblue",
                           baseline  = NULL,
                           dose_days = NULL,
                           obs       = NULL) {

  # Calcul des limites dynamiques (simulations + obs)
  all_vals <- c(vpc_df$p05, vpc_df$p95)
  if (!is.null(obs) && length(obs$value) > 0)
    all_vals <- c(all_vals, obs$value[obs$value > 0])

  all_pos <- all_vals[is.finite(all_vals) & all_vals > 0]
  ymin    <- 10^floor(log10(min(all_pos) * 0.3))
  ymax    <- 10^ceiling(log10(max(all_pos) * 3.0))
  breaks_y <- 10^seq(log10(ymin), log10(ymax), by = 1)

  # Clamp les valeurs à ymin
  vpc_df$p05  <- pmax(vpc_df$p05,  ymin * 0.5)
  vpc_df$p95  <- pmax(vpc_df$p95,  ymin * 0.5)
  vpc_df$p50  <- pmax(vpc_df$p50,  ymin * 0.5)
  vpc_df$pred <- pmax(vpc_df$pred, ymin * 0.5)

  p <- ggplot(vpc_df, aes(x = days)) +

    # Ruban bleu — intervalle de prédiction [5e, 95e]
    geom_ribbon(aes(ymin = p05, ymax = p95),
                fill = color, alpha = 0.20) +

    # Médiane simulée (rouge, comme Figure 3 Fornari)
    geom_line(aes(y = p50),
              color = "red", linewidth = 0.9, alpha = 0.9) +

    # Prédiction déterministe (noire fine, optionnel)
    geom_line(aes(y = pred),
              color = "black", linewidth = 0.4,
              linetype = "dashed", alpha = 0.5)

  # Points observés
  if (!is.null(obs) && length(obs$time) > 0) {
    df_obs <- data.frame(
      x = obs$time,
      y = pmax(obs$value, ymin * 0.5)
    )
    p <- p + geom_point(data  = df_obs, aes(x = x, y = y),
                        shape = 21, fill = "white",
                        color = "black", size = 1.8,
                        stroke = 0.65, alpha = 0.85)
  }

  # Baseline (tiretée grise)
  if (!is.null(baseline) && baseline > 0)
    p <- p + geom_hline(yintercept = baseline,
                        linetype = "dashed",
                        color = "grey50",
                        linewidth = 0.45, alpha = 0.6)

  # Lignes de dose (pointillées)
  if (!is.null(dose_days)) {
    vd <- dose_days[dose_days <= max(vpc_df$days)]
    if (length(vd) > 0)
      p <- p + geom_vline(xintercept = vd, linetype = "dotted",
                          color = "grey55",
                          linewidth = 0.35, alpha = 0.55)
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
# 3. GRILLE COMPLÈTE 4×2 VPC
# ════════════════════════════════════════════════════════
plot_vpc_all <- function(vpc_stats, pars,
                         titre     = "VPC — Haematopoiesis",
                         dose_days = NULL,
                         obs_list  = NULL) {

  get_obs <- function(name) {
    if (!is.null(obs_list) && name %in% names(obs_list))
      obs_list[[name]]
    else NULL
  }

  # Couleurs identiques à plots.R pour cohérence visuelle
  colors <- c(MPP  = "#2166ac", CMP  = "#4dac26", MEP  = "#7b3294",
              Neut = "#d7191c", Mono = "#a6611a", Plt  = "#1b9e77",
              Ret  = "#e377c2", RBC  = "#cc3399")

  baselines <- list(MPP = pars$MPP0, CMP = pars$CMP0, MEP = pars$MEP0,
                    Neut = pars$Neut0, Mono = pars$Mono0, Plt = pars$Plt0,
                    Ret = pars$Ret0, RBC = pars$RBC0)

  titles <- c(MPP  = "Multi-potent progenitors",
              CMP  = "Common myeloid progenitors",
              MEP  = "MEP",
              Neut = "Neutrophils",
              Mono = "Monocytes",
              Plt  = "Platelets",
              Ret  = "Reticulocytes",
              RBC  = "Red blood cells")

  # Ordre des panneaux identique à Figure 3 Fornari
  order <- c("MPP", "Neut", "CMP", "Mono", "MEP", "Plt", "Ret", "RBC")

  panels <- lapply(order, function(cell) {
    plot_vpc_panel(
      vpc_df    = vpc_stats[[cell]],
      title     = titles[cell],
      color     = colors[cell],
      baseline  = baselines[[cell]],
      dose_days = dose_days,
      obs       = get_obs(cell)
    )
  })

  # Légende commune sous forme de texte
  legend_grob <- textGrob(
    paste0("Ruban bleu = [5e–95e percentile]  |  ",
           "Ligne rouge = médiane simulée  |  ",
           "Tiret noir = prédiction déterministe  |  ",
           "● = données observées"),
    gp = gpar(fontsize = 7.5, col = "grey30")
  )

  grid.arrange(
    grobs = panels,
    ncol  = 2,
    top   = textGrob(titre,
                     gp = gpar(fontface = "bold", fontsize = 11)),
    bottom = legend_grob
  )
}

# ════════════════════════════════════════════════════════
# 4. SAUVEGARDE PDF
# ════════════════════════════════════════════════════════
save_vpc <- function(sim, pars, file, titre,
                     n_sim     = 1000,
                     dose_days = NULL,
                     obs_list  = NULL,
                     width     = 12,
                     height    = 11) {

  cat(sprintf("  → Génération VPC (%d simulations)...\n", n_sim))
  vpc_stats <- run_vpc(sim, n_sim = n_sim)

  pdf(file, width = width, height = height)
  plot_vpc_all(vpc_stats, pars,
               titre     = titre,
               dose_days = dose_days,
               obs_list  = obs_list)
  dev.off()

  message("✓ VPC sauvegardé : ", file)
  invisible(vpc_stats)
}

# ════════════════════════════════════════════════════════
# 5. VPC HUMAIN — Neutrophiles + Plaquettes (Figure 4)
#    Deux panneaux empilés verticalement (1 colonne × 2 lignes)
#    comme dans la Figure 4 de Fornari 2019
# ════════════════════════════════════════════════════════
plot_vpc_human <- function(vpc_stats, pars,
                           titre     = "VPC — Human (Figure 4)",
                           dose_days = NULL,
                           obs_list  = NULL) {

  get_obs <- function(name) {
    if (!is.null(obs_list) && name %in% names(obs_list))
      obs_list[[name]]
    else NULL
  }

  p_neut <- plot_vpc_panel(
    vpc_df    = vpc_stats[["Neut"]],
    title     = "Neutrophils",
    color     = "#d7191c",
    baseline  = pars$Neut0,
    dose_days = dose_days,
    obs       = get_obs("Neut")
  )

  p_plt <- plot_vpc_panel(
    vpc_df    = vpc_stats[["Plt"]],
    title     = "Platelets",
    color     = "#1b9e77",
    baseline  = pars$Plt0,
    dose_days = dose_days,
    obs       = get_obs("Plt")
  )

  legend_grob <- textGrob(
    paste0("Ruban coloré = [5e–95e percentile]  |  ",
           "Ligne rouge = médiane simulée  |  ",
           "Tiret noir = prédiction déterministe  |  ",
           "● = données observées"),
    gp = gpar(fontsize = 7.5, col = "grey30")
  )

  grid.arrange(
    p_neut, p_plt,
    ncol   = 1,
    top    = textGrob(titre,
                      gp = gpar(fontface = "bold", fontsize = 11)),
    bottom = legend_grob
  )
}

save_vpc_human <- function(sim, pars, file, titre,
                           n_sim      = 1000,
                           dose_days  = NULL,
                           obs_list   = NULL,
                           auc_target = 5,
                           gfr_fixed  = 125,
                           times      = NULL,
                           interval_h = 21 * 24,
                           n_cycles   = 2,
                           width      = 7,
                           height     = 9) {

  # Supp. S11 : GFR=125 fixe, même PK pour tous les patients
  # Optimisation : run_vpc() applique le bruit résiduel de façon vectorisée
  # sur la simulation déterministe (sim) — 1 seul appel ODE au lieu de n_sim.
  cat(sprintf("  → VPC humain : %d simulations (bruit résiduel vectorisé)...\n", n_sim))

  set.seed(42)
  vpc_all <- run_vpc(sim, n_sim = n_sim)

  vpc_stats <- list(
    Neut = vpc_all[["Neut"]][, c("days","p05","p50","p95","pred")],
    Plt  = vpc_all[["Plt"]] [, c("days","p05","p50","p95","pred")]
  )

  pdf(file, width = width, height = height)
  plot_vpc_human(vpc_stats, pars,
                 titre     = titre,
                 dose_days = dose_days,
                 obs_list  = obs_list)
  dev.off()

  message("✓ VPC humain sauvegardé : ", file)
  invisible(vpc_stats)
}

# ════════════════════════════════════════════════════════
# 6. VPC ÉLARGI — variabilité des baselines Neut0 + Plt0
#    Neut0 ~ log-uniforme [2, 7]   (range adulte, Fornari Table 1)
#    Plt0  ~ log-uniforme [150, 400] (range adulte, littérature)
#    + erreur résiduelle log-additive (SIGMA, Table S4)
# ════════════════════════════════════════════════════════
save_vpc_human_extended <- function(sim, pars, file, titre,
                                    n_sim      = 1000,
                                    dose_days  = NULL,
                                    obs_list   = NULL,
                                    auc_target = 5,
                                    gfr_fixed  = 78,
                                    times      = NULL,
                                    interval_h = 21 * 24,
                                    n_cycles   = 2,
                                    neut0_range = c(2.0, 7.0),
                                    plt0_range  = c(150.0, 400.0),
                                    width      = 7,
                                    height     = 9) {

  cat(sprintf(
    "  → VPC élargi : %d patients | Neut0 ~ LogU[%.1f,%.1f] | Plt0 ~ LogU[%.1f,%.1f]\n",
    n_sim, neut0_range[1], neut0_range[2], plt0_range[1], plt0_range[2]))

  set.seed(42)
  if (is.null(times)) times <- seq(0, 63 * 24, by = 1)
  n_t <- length(times)

  dose_fixe <- auc_target * (gfr_fixed + 25)

  # Tirages log-uniformes sur les ranges physiologiques
  neut0_vec <- exp(runif(n_sim,
                         log(neut0_range[1]), log(neut0_range[2])))
  plt0_vec  <- exp(runif(n_sim,
                         log(plt0_range[1]),  log(plt0_range[2])))

  mat_Neut <- matrix(NA, nrow = n_sim, ncol = n_t)
  mat_Plt  <- matrix(NA, nrow = n_sim, ncol = n_t)

  for (i in seq_len(n_sim)) {
    pars_i        <- pars
    pars_i$Neut0  <- neut0_vec[i]
    pars_i$Plt0   <- plt0_vec[i]
    pars_i$rate_fun <- make_repeated_infusion(
      dose_mg = dose_fixe, Tinfu_h = 1,
      interval_h = interval_h, n_cycles = n_cycles
    )

    # État initial rééquilibré pour ce patient
    a_Neut <- 3 / pars_i$MTT_Neut
    T_Neut <- pars_i$k_circ_Neut * neut0_vec[i] / a_Neut
    a_Plt  <- 3 / pars_i$MTT_Plt
    T_Plt  <- pars_i$k_circ_Plt  * plt0_vec[i]  / a_Plt
    T1_Plt <- T_Plt / pars_i$lambda2

    ss_T <- function(k, c0, mtt) k * c0 / (3 / mtt)
    state_i <- c(
      C1=0, C2=0, Damage=0,
      MPP=pars$MPP0, CMP=pars$CMP0, MEP=pars$MEP0,
      T1_Neut=T_Neut, T2_Neut=T_Neut, T3_Neut=T_Neut,
      Neut=neut0_vec[i],
      T1_Mono=ss_T(pars$k_circ_Mono,pars$Mono0,pars$MTT_Mono),
      T2_Mono=ss_T(pars$k_circ_Mono,pars$Mono0,pars$MTT_Mono),
      T3_Mono=ss_T(pars$k_circ_Mono,pars$Mono0,pars$MTT_Mono),
      Mono=pars$Mono0,
      T1_Ret=ss_T(pars$k_circ_RBC,pars$Ret0,pars$MTT_Ret),
      T2_Ret=ss_T(pars$k_circ_RBC,pars$Ret0,pars$MTT_Ret),
      T3_Ret=ss_T(pars$k_circ_RBC,pars$Ret0,pars$MTT_Ret),
      Ret=pars$Ret0, RBC=pars$RBC0,
      T1_Plt=T1_Plt, T2_Plt=T_Plt, T3_Plt=T_Plt,
      Plt=plt0_vec[i]
    )

    tryCatch({
      out_i <- as.data.frame(lsoda(
        y=state_i, times=times, func=pkpd_fornari, parms=pars_i,
        rtol=1e-3, atol=1e-5, maxsteps=50000
      ))
      mat_Neut[i,] <- out_i$Neut * exp(rnorm(n_t, 0, SIGMA["Neut"]))
      mat_Plt[i,]  <- out_i$Plt  * exp(rnorm(n_t, 0, SIGMA["Plt"]))
    }, error = function(e) NULL)

    if (i %% 50 == 0)
      cat(sprintf("    %d/%d patients simulés\n", i, n_sim))
  }

  days_vec   <- times / 24
  make_stats <- function(mat, Yref) {
    # Interpoler Yref si grille temporelle différente de mat
    pred_interp <- if (length(Yref) == length(days_vec)) {
      Yref
    } else {
      approx(sim$time / 24, Yref, days_vec, rule = 2)$y
    }
    data.frame(
      days = days_vec,
      p05  = apply(mat, 2, quantile, probs=0.05, na.rm=TRUE),
      p50  = apply(mat, 2, quantile, probs=0.50, na.rm=TRUE),
      p95  = apply(mat, 2, quantile, probs=0.95, na.rm=TRUE),
      pred = pred_interp
    )
  }
  vpc_stats <- list(
    Neut = make_stats(mat_Neut, sim$Neut),
    Plt  = make_stats(mat_Plt,  sim$Plt)
  )

  pdf(file, width = width, height = height)
  plot_vpc_human(vpc_stats, pars,
                 titre     = titre,
                 dose_days = dose_days,
                 obs_list  = obs_list)
  dev.off()

  message("✓ VPC élargi sauvegardé : ", file)
  invisible(vpc_stats)
}

# ════════════════════════════════════════════════════════
# 7. VPC POPULATION avec IIV PD complète (Figure 4 style)
#    IIV : omega_CL=0.35, omega_Slope_MEP=0.547, etc.
#    Bandes [5e–95e] + [25e–75e] + médiane + seuils grade
# ════════════════════════════════════════════════════════
plot_vpc_iiv_panel <- function(vpc_df, title, color,
                               thresholds = NULL,
                               baseline   = NULL,
                               dose_days  = NULL,
                               obs        = NULL) {

  all_vals <- c(vpc_df$p05, vpc_df$p95)
  if (!is.null(obs) && length(obs$value) > 0)
    all_vals <- c(all_vals, obs$value[obs$value > 0])
  all_pos  <- all_vals[is.finite(all_vals) & all_vals > 0]
  ymin     <- 10^floor(log10(min(all_pos) * 0.3))
  ymax     <- 10^ceiling(log10(max(all_pos) * 3.0))
  breaks_y <- 10^seq(log10(ymin), log10(ymax), by = 1)

  vpc_df <- as.data.frame(lapply(vpc_df, function(x)
    if (is.numeric(x)) pmax(x, ymin * 0.5) else x))

  p <- ggplot(vpc_df, aes(x = days)) +
    geom_ribbon(aes(ymin = p05, ymax = p95),
                fill = color, alpha = 0.15) +
    geom_ribbon(aes(ymin = p25, ymax = p75),
                fill = color, alpha = 0.25) +
    geom_line(aes(y = p50), color = color, linewidth = 1.0)

  if (!is.null(thresholds)) {
    grade_cols <- c("#fee08b", "#fc8d59", "#d73027", "#7b0404")
    grade_lbls <- names(thresholds)
    for (g in seq_along(thresholds)) {
      p <- p + geom_hline(yintercept = thresholds[g],
                          linetype = "dashed", linewidth = 0.4,
                          color = grade_cols[g], alpha = 0.8) +
               annotate("text", x = max(vpc_df$days) * 0.02,
                        y = thresholds[g] * 1.18,
                        label = grade_lbls[g], size = 2.5,
                        color = grade_cols[g], hjust = 0)
    }
  }

  if (!is.null(obs) && length(obs$time) > 0) {
    df_obs <- data.frame(x = obs$time,
                         y = pmax(obs$value, ymin * 0.5))
    p <- p + geom_point(data = df_obs, aes(x = x, y = y),
                        shape = 21, fill = "white", color = "black",
                        size = 1.8, stroke = 0.65, alpha = 0.85)
  }

  if (!is.null(baseline) && baseline > 0)
    p <- p + geom_hline(yintercept = baseline, linetype = "dashed",
                        color = "grey50", linewidth = 0.45, alpha = 0.6)

  if (!is.null(dose_days)) {
    vd <- dose_days[dose_days <= max(vpc_df$days)]
    if (length(vd) > 0)
      p <- p + geom_vline(xintercept = vd, linetype = "dotted",
                          color = "grey55", linewidth = 0.35, alpha = 0.55)
  }

  p + scale_y_log10(limits = c(ymin, ymax), breaks = breaks_y,
                    labels = trans_format("log10", math_format(10^.x))) +
    labs(title = title, x = "Time (d)",
         y = expression(10^9~cells~L^{-1})) +
    theme_bw(base_size = 9.5) +
    theme(panel.grid.minor  = element_blank(),
          panel.grid.major  = element_line(color = "grey92"),
          plot.title        = element_text(face = "bold", size = 9, hjust = 0.5),
          axis.title        = element_text(size = 7.5),
          axis.text         = element_text(size = 7))
}

save_vpc_human_pd_iiv <- function(pars, init_state_arg,
                                  file,
                                  titre      = "VPC — IIV PD complète",
                                  n_sim      = 500,
                                  auc_target = 5,
                                  gfr_fixed  = 125,
                                  n_cycles   = 2,
                                  interval_h = 21 * 24,
                                  dose_days  = NULL,
                                  obs_list   = NULL,
                                  seed       = 42,
                                  width      = 7,
                                  height     = 9) {

  set.seed(seed)
  times <- seq(0, n_cycles * interval_h + 21 * 24, by = 1)
  n_t   <- length(times)

  omega_CL        <- 0.35
  omega_Slope_CMP <- 0.624
  omega_Slope_MEP <- 0.547
  omega_Slope_MPP <- 0.624
  omega_Neut0     <- 0.326
  omega_Plt0      <- 0.268

  CL_typical <- (gfr_fixed + 25) * 60 / 1000
  dose_fixe  <- auc_target * (gfr_fixed + 25)

  eta_CL    <- rnorm(n_sim, 0, omega_CL)
  eta_CMP   <- rnorm(n_sim, 0, omega_Slope_CMP)
  eta_MEP   <- rnorm(n_sim, 0, omega_Slope_MEP)
  eta_MPP   <- rnorm(n_sim, 0, omega_Slope_MPP)
  eta_Neut0 <- rnorm(n_sim, 0, omega_Neut0)
  eta_Plt0  <- rnorm(n_sim, 0, omega_Plt0)

  neut0_i <- pars$Neut0 * exp(eta_Neut0)
  plt0_i  <- pars$Plt0  * exp(eta_Plt0)

  mat_Neut <- matrix(NA, nrow = n_sim, ncol = n_t)
  mat_Plt  <- matrix(NA, nrow = n_sim, ncol = n_t)

  cat(sprintf("  → VPC IIV PD : %d patients (GFR=%g, AUC=%g)...\n",
              n_sim, gfr_fixed, auc_target))

  for (i in seq_len(n_sim)) {
    pars_i           <- pars
    pars_i$CL        <- CL_typical * exp(eta_CL[i])
    pars_i$Slope_MPP <- pars$Slope_MPP * exp(eta_MPP[i])
    pars_i$Slope_CMP <- pars$Slope_CMP * exp(eta_CMP[i])
    pars_i$Slope_MEP <- pars$Slope_MEP * exp(eta_MEP[i])
    pars_i$Neut0     <- neut0_i[i]
    pars_i$Plt0      <- plt0_i[i]
    pars_i$rate_fun  <- make_repeated_infusion(
      dose_mg = dose_fixe, Tinfu_h = 1,
      interval_h = interval_h, n_cycles = n_cycles
    )

    a_Neut  <- 3 / pars_i$MTT_Neut
    a_Plt   <- 3 / pars_i$MTT_Plt
    T_Neut  <- pars_i$k_circ_Neut * neut0_i[i] / a_Neut
    T_Plt   <- pars_i$k_circ_Plt  * plt0_i[i]  / a_Plt
    T1_Plt  <- T_Plt / pars_i$lambda2

    state_i              <- init_state_arg
    state_i["Neut"]      <- neut0_i[i]
    state_i["Plt"]       <- plt0_i[i]
    state_i["T1_Neut"]   <- T_Neut
    state_i["T2_Neut"]   <- T_Neut
    state_i["T3_Neut"]   <- T_Neut
    state_i["T1_Plt"]    <- T1_Plt
    state_i["T2_Plt"]    <- T_Plt
    state_i["T3_Plt"]    <- T_Plt

    tryCatch({
      out_i <- as.data.frame(lsoda(
        y = state_i, times = times, func = pkpd_fornari, parms = pars_i,
        rtol = 1e-3, atol = 1e-5, maxsteps = 50000
      ))
      mat_Neut[i, ] <- out_i$Neut
      mat_Plt[i, ]  <- out_i$Plt
    }, error = function(e) NULL)

    if (i %% 100 == 0)
      cat(sprintf("    %d/%d\n", i, n_sim))
  }

  days_vec <- times / 24
  make_stats <- function(mat) data.frame(
    days = days_vec,
    p05  = apply(mat, 2, quantile, probs = 0.05, na.rm = TRUE),
    p25  = apply(mat, 2, quantile, probs = 0.25, na.rm = TRUE),
    p50  = apply(mat, 2, quantile, probs = 0.50, na.rm = TRUE),
    p75  = apply(mat, 2, quantile, probs = 0.75, na.rm = TRUE),
    p95  = apply(mat, 2, quantile, probs = 0.95, na.rm = TRUE)
  )

  vpc_stats <- list(
    Neut = make_stats(mat_Neut),
    Plt  = make_stats(mat_Plt)
  )

  neut_thr <- c(G1 = 2.0, G2 = 1.5, G3 = 1.0, G4 = 0.5)
  plt_thr  <- c(G1 = 150, G2 = 75,  G3 = 50,  G4 = 25)

  get_obs <- function(name)
    if (!is.null(obs_list) && name %in% names(obs_list)) obs_list[[name]] else NULL

  p_neut <- plot_vpc_iiv_panel(
    vpc_df    = vpc_stats$Neut,
    title     = "Neutrophils — population IIV",
    color     = "#d7191c",
    thresholds = neut_thr,
    baseline  = pars$Neut0,
    dose_days = dose_days,
    obs       = get_obs("Neut")
  )
  p_plt <- plot_vpc_iiv_panel(
    vpc_df    = vpc_stats$Plt,
    title     = "Platelets — population IIV",
    color     = "#1b9e77",
    thresholds = plt_thr,
    baseline  = pars$Plt0,
    dose_days = dose_days,
    obs       = get_obs("Plt")
  )

  legend_grob <- textGrob(
    paste0("Ruban clair = [5e–95e percentile]  |  ",
           "Ruban foncé = [25e–75e]  |  ",
           "Ligne = médiane  |  ● = données"),
    gp = gpar(fontsize = 7.5, col = "grey30")
  )

  pdf(file, width = width, height = height)
  grid.arrange(p_neut, p_plt, ncol = 1,
               top    = textGrob(titre,
                                 gp = gpar(fontface = "bold", fontsize = 11)),
               bottom = legend_grob)
  dev.off()

  # Résumé des percentiles au nadir
  plt_min_p05 <- min(vpc_stats$Plt$p05, na.rm = TRUE)
  plt_min_p50 <- min(vpc_stats$Plt$p50, na.rm = TRUE)
  plt_min_p95 <- min(vpc_stats$Plt$p95, na.rm = TRUE)
  cat(sprintf("  Plt nadir : p05=%.1f  p50=%.1f  p95=%.1f\n",
              plt_min_p05, plt_min_p50, plt_min_p95))
  cat(sprintf("  → p05 Plt < 150 (G1)? %s  < 75 (G2)? %s  < 50 (G3)? %s\n",
              ifelse(plt_min_p05 < 150, "OUI", "non"),
              ifelse(plt_min_p05 < 75,  "OUI", "non"),
              ifelse(plt_min_p05 < 50,  "OUI", "non")))

  message("✓ VPC IIV PD sauvegardé : ", file)
  invisible(vpc_stats)
}
