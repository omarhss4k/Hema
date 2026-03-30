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
  cat(sprintf("  → VPC humain : %d patients (GFR=%g fixe — Supp. S11) + erreur résiduelle...\n",
              n_sim, gfr_fixed))

  set.seed(42)

  # Grille temporelle identique à run_human.R si non fournie
  if (is.null(times)) times <- seq(0, 63 * 24, by = 1)
  n_t <- length(times)

  # ── Même PK pour tous (Supp. S11/S12) ──────────────────────
  dose_fixe <- auc_target * (gfr_fixed + 25)   # Calvert dose fixe
  mat_Neut <- matrix(NA, nrow = n_sim, ncol = n_t)
  mat_Plt  <- matrix(NA, nrow = n_sim, ncol = n_t)

  # ── IIV (proxy ω = σ Fornari Table S4, même convention que plots_grades.R) ──
  omega_Slope_MPP <- SIGMA["MPP"]
  omega_Slope_CMP <- SIGMA["CMP"]
  omega_Slope_MEP <- SIGMA["MEP"]
  omega_Neut0     <- SIGMA["Neut"]
  omega_Plt0      <- SIGMA["Plt"]

  # Tirages IIV (log-normaux, indépendants) — 1 par patient
  eta_MPP   <- rnorm(n_sim, 0, omega_Slope_MPP)
  eta_CMP   <- rnorm(n_sim, 0, omega_Slope_CMP)
  eta_MEP   <- rnorm(n_sim, 0, omega_Slope_MEP)
  eta_Neut0 <- rnorm(n_sim, 0, omega_Neut0)
  eta_Plt0  <- rnorm(n_sim, 0, omega_Plt0)

  ss_T <- function(k, c0, mtt) k * c0 / (3 / mtt)

  for (i in seq_len(n_sim)) {
    dose_i <- dose_fixe

    # IIV PD : sensibilité médicament (log-normale, Fornari Table S4 comme proxy ω)
    pars_i             <- pars
    pars_i$Slope_MPP   <- pars$Slope_MPP * exp(eta_MPP[i])
    pars_i$Slope_CMP   <- pars$Slope_CMP * exp(eta_CMP[i])
    pars_i$Slope_MEP   <- pars$Slope_MEP * exp(eta_MEP[i])

    # IIV PD : baselines individuelles (log-normales)
    neut0_i            <- pars$Neut0 * exp(eta_Neut0[i])
    plt0_i             <- pars$Plt0  * exp(eta_Plt0[i])
    pars_i$Neut0       <- neut0_i
    pars_i$Plt0        <- plt0_i

    pars_i$rate_fun <- make_repeated_infusion(
      dose_mg    = dose_i,
      Tinfu_h    = 1,
      interval_h = interval_h,
      n_cycles   = n_cycles
    )

    state0 <- c(
      C1=0, C2=0, Damage=0,
      MPP=pars$MPP0, CMP=pars$CMP0, MEP=pars$MEP0,
      T1_Neut=ss_T(pars$k_circ_Neut,neut0_i,pars$MTT_Neut),
      T2_Neut=ss_T(pars$k_circ_Neut,neut0_i,pars$MTT_Neut),
      T3_Neut=ss_T(pars$k_circ_Neut,neut0_i,pars$MTT_Neut),
      Neut=neut0_i,
      T1_Mono=ss_T(pars$k_circ_Mono,pars$Mono0,pars$MTT_Mono),
      T2_Mono=ss_T(pars$k_circ_Mono,pars$Mono0,pars$MTT_Mono),
      T3_Mono=ss_T(pars$k_circ_Mono,pars$Mono0,pars$MTT_Mono),
      Mono=pars$Mono0,
      T1_Ret=ss_T(pars$k_circ_RBC,pars$Ret0,pars$MTT_Ret),
      T2_Ret=ss_T(pars$k_circ_RBC,pars$Ret0,pars$MTT_Ret),
      T3_Ret=ss_T(pars$k_circ_RBC,pars$Ret0,pars$MTT_Ret),
      Ret=pars$Ret0, RBC=pars$RBC0,
      T1_Plt=ss_T(pars$k_circ_Plt,plt0_i,pars$MTT_Plt),
      T2_Plt=ss_T(pars$k_circ_Plt,plt0_i,pars$MTT_Plt),
      T3_Plt=ss_T(pars$k_circ_Plt,plt0_i,pars$MTT_Plt),
      Plt=plt0_i
    )

    tryCatch({
      out_i <- as.data.frame(lsoda(
        y=state0, times=times, func=pkpd_fornari, parms=pars_i,
        rtol=1e-4, atol=1e-6, maxsteps=10000
      ))
      mat_Neut[i,] <- out_i$Neut * exp(rnorm(n_t, 0, SIGMA["Neut"]))
      mat_Plt[i,]  <- out_i$Plt  * exp(rnorm(n_t, 0, SIGMA["Plt"]))
    }, error = function(e) NULL)

    if (i %% 200 == 0)
      cat(sprintf("    %d/%d patients simulés\n", i, n_sim))
  }

  # ── Percentiles ──────────────────────────────────────────────
  days_vec  <- times / 24
  make_stats <- function(mat, Yref) data.frame(
    days = days_vec,
    p05  = apply(mat, 2, quantile, probs=0.05, na.rm=TRUE),
    p50  = apply(mat, 2, quantile, probs=0.50, na.rm=TRUE),
    p95  = apply(mat, 2, quantile, probs=0.95, na.rm=TRUE),
    pred = Yref
  )
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

  message("✓ VPC humain sauvegardé : ", file)
  invisible(vpc_stats)
}
