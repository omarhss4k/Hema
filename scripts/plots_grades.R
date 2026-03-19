############################################################
# plots_grades.R
# Grades NCI-CTCAE v5.0 — Neutropénie + Thrombocytopénie
# Fonctions :
#   save_grade_plots()    → rat  (Figure 3 style)
#   save_grade_figure4c() → humain population (Figure 4c)
############################################################
library(ggplot2)
library(gridExtra)
library(grid)
library(deSolve)

# ════════════════════════════════════════════════════════
# SEUILS NCI-CTCAE v5.0
# ════════════════════════════════════════════════════════
# Neutropénie (10^9 cells/L)
NEUT_THRESHOLDS <- c(G1 = 2.0, G2 = 1.5, G3 = 1.0, G4 = 0.5)

# Thrombocytopénie (10^9 cells/L)
PLT_THRESHOLDS  <- c(G1 = 150, G2 = 75, G3 = 50, G4 = 25)

# Couleurs des grades
GRADE_COLORS <- c(G1 = "#fee08b", G2 = "#fc8d59",
                  G3 = "#d73027", G4 = "#7b0404")

# ════════════════════════════════════════════════════════
# 1. FONCTION UTILITAIRE — grade d'une valeur scalaire
# ════════════════════════════════════════════════════════
assign_grade_neut <- function(x) {
  ifelse(x <  NEUT_THRESHOLDS["G4"], 4,
  ifelse(x <  NEUT_THRESHOLDS["G3"], 3,
  ifelse(x <  NEUT_THRESHOLDS["G2"], 2,
  ifelse(x <  NEUT_THRESHOLDS["G1"], 1, 0))))
}

assign_grade_plt <- function(x) {
  ifelse(x <  PLT_THRESHOLDS["G4"], 4,
  ifelse(x <  PLT_THRESHOLDS["G3"], 3,
  ifelse(x <  PLT_THRESHOLDS["G2"], 2,
  ifelse(x <  PLT_THRESHOLDS["G1"], 1, 0))))
}

# Nadir grade = grade de la valeur minimale
nadir_grade_neut <- function(vec) assign_grade_neut(min(vec, na.rm = TRUE))
nadir_grade_plt  <- function(vec) assign_grade_plt( min(vec, na.rm = TRUE))

# ════════════════════════════════════════════════════════
# 2. PANNEAU GRADE — courbe + bandes de grade colorées
# ════════════════════════════════════════════════════════
plot_grade_panel <- function(sim, yvar, title, color,
                             thresholds,
                             baseline  = NULL,
                             dose_days = NULL) {

  y     <- sim[[yvar]]
  days  <- sim$days
  ymin  <- min(thresholds) * 0.1
  ymax  <- max(y, na.rm = TRUE) * 3
  ymin  <- max(ymin, 1e-3)

  breaks_y <- 10^seq(floor(log10(ymin)), ceiling(log10(ymax)), by = 1)

  df <- data.frame(x = days, y = pmax(y, ymin * 0.5))

  p <- ggplot() +

    # Bandes de grade (du plus sévère au moins sévère)
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=ymin,              ymax=thresholds["G4"],
             fill=GRADE_COLORS["G4"], alpha=0.15) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=thresholds["G4"],  ymax=thresholds["G3"],
             fill=GRADE_COLORS["G3"], alpha=0.13) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=thresholds["G3"],  ymax=thresholds["G2"],
             fill=GRADE_COLORS["G2"], alpha=0.11) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=thresholds["G2"],  ymax=thresholds["G1"],
             fill=GRADE_COLORS["G1"], alpha=0.09) +

    # Lignes de seuil
    geom_hline(yintercept = as.numeric(thresholds),
               linetype = "dashed", linewidth = 0.35,
               color = c(GRADE_COLORS["G1"], GRADE_COLORS["G2"],
                         GRADE_COLORS["G3"], GRADE_COLORS["G4"]),
               alpha = 0.8) +

    # Labels grades (à droite)
    annotate("text", x = max(days)*0.97, y = thresholds["G1"]*1.15,
             label="G1", size=2.5, color=GRADE_COLORS["G1"], hjust=1) +
    annotate("text", x = max(days)*0.97, y = thresholds["G2"]*1.15,
             label="G2", size=2.5, color=GRADE_COLORS["G2"], hjust=1) +
    annotate("text", x = max(days)*0.97, y = thresholds["G3"]*1.15,
             label="G3", size=2.5, color=GRADE_COLORS["G3"], hjust=1) +
    annotate("text", x = max(days)*0.97, y = thresholds["G4"]*1.15,
             label="G4", size=2.5, color=GRADE_COLORS["G4"], hjust=1) +

    # Courbe simulation
    geom_line(data = df, aes(x = x, y = y),
              color = color, linewidth = 1.0)

  # Baseline
  if (!is.null(baseline) && baseline > 0)
    p <- p + geom_hline(yintercept = baseline, linetype = "dashed",
                        color = "grey50", linewidth = 0.4, alpha = 0.5)

  # Doses
  if (!is.null(dose_days)) {
    vd <- dose_days[dose_days <= max(days)]
    if (length(vd) > 0)
      p <- p + geom_vline(xintercept = vd, linetype = "dotted",
                          color = "grey55", linewidth = 0.35, alpha = 0.55)
  }

  p +
    scale_y_log10(limits = c(ymin, ymax), breaks = breaks_y,
                  labels = scales::trans_format("log10",
                                                scales::math_format(10^.x))) +
    labs(title = title, x = "Time (d)",
         y = expression(10^9~cells~L^{-1})) +
    theme_bw(base_size = 9.5) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92"),
          plot.title  = element_text(face = "bold", size = 9, hjust = 0.5),
          axis.title  = element_text(size = 7.5),
          axis.text   = element_text(size = 7))
}

# ════════════════════════════════════════════════════════
# 3. SAVE_GRADE_PLOTS — rat (Figure 3 style)
#    Retourne liste avec nadir/grade/pct pour console
# ════════════════════════════════════════════════════════
save_grade_plots <- function(sim, pars, file, titre_base,
                             dose_days = NULL,
                             width = 10, height = 5) {

  p_neut <- plot_grade_panel(
    sim, "Neut", "Neutrophils (NCI-CTCAE)", "#d7191c",
    NEUT_THRESHOLDS, pars$Neut0, dose_days
  )
  p_plt <- plot_grade_panel(
    sim, "Plt", "Platelets (NCI-CTCAE)", "#1b9e77",
    PLT_THRESHOLDS, pars$Plt0, dose_days
  )

  pdf(file, width = width, height = height)
  grid.arrange(p_neut, p_plt, ncol = 2,
               top = grid::textGrob(titre_base,
                                    gp = grid::gpar(fontface="bold",
                                                    fontsize=11)))
  dev.off()
  message("✓ Grades sauvegardés : ", file)

  # Statistiques résumées
  nadir_neut <- min(sim$Neut, na.rm = TRUE)
  nadir_plt  <- min(sim$Plt,  na.rm = TRUE)
  g_neut     <- nadir_grade_neut(sim$Neut)
  g_plt      <- nadir_grade_plt(sim$Plt)

  # % temps passé dans chaque grade (proxy de % patients)
  pct_neut <- sapply(1:4, function(g) {
    thr_lo <- c(0, NEUT_THRESHOLDS["G4"], NEUT_THRESHOLDS["G3"],
                NEUT_THRESHOLDS["G2"])[g]
    thr_hi <- c(NEUT_THRESHOLDS["G4"], NEUT_THRESHOLDS["G3"],
                NEUT_THRESHOLDS["G2"], NEUT_THRESHOLDS["G1"])[g]
    100 * mean(sim$Neut >= thr_lo & sim$Neut < thr_hi, na.rm=TRUE)
  })
  pct_plt <- sapply(1:4, function(g) {
    thr_lo <- c(0, PLT_THRESHOLDS["G4"], PLT_THRESHOLDS["G3"],
                PLT_THRESHOLDS["G2"])[g]
    thr_hi <- c(PLT_THRESHOLDS["G4"], PLT_THRESHOLDS["G3"],
                PLT_THRESHOLDS["G2"], PLT_THRESHOLDS["G1"])[g]
    100 * mean(sim$Plt >= thr_lo & sim$Plt < thr_hi, na.rm=TRUE)
  })

  invisible(list(
    nadir_neut       = nadir_neut,
    grade_nadir_neut = g_neut,
    pct_neut         = pct_neut,
    nadir_plt        = nadir_plt,
    grade_nadir_plt  = g_plt,
    pct_plt          = pct_plt
  ))
}

# ════════════════════════════════════════════════════════
# 4. SAVE_GRADE_FIGURE4C — humain population (Figure 4c)
#    Simule n_patients avec GFR variable → % par grade
# ════════════════════════════════════════════════════════
save_grade_figure4c <- function(base_pars, init_state,
                                auc_target = 5,
                                n_cycles   = 2,
                                interval_h = 21 * 24,
                                n_patients = 500,
                                gfr_mean   = 78,
                                gfr_sd     = 20,
                                file       = "Figure4c_grades.pdf",
                                titre      = "% patients par grade",
                                seed       = 42,
                                width      = 8,
                                height     = 5) {

  set.seed(seed)
  cat(sprintf("  → Figure 4c : %d patients (GFR~N(%g,%g) + IIV CL/Neut0/Plt0)...\n",
              n_patients, gfr_mean, gfr_sd))

  # ── Sources de variabilité inter-individuelle ──────────────
  # 1. IIV log-normale sur CL résiduel (au-delà de la prédiction GFR)
  #    omega_CL = 0.35 (Zandvliet 2008 Table 3 : ~30-35% CV)
  omega_CL   <- 0.35
  # 2. IIV log-normale sur les Slopes (sensibilité médicament)
  #    omega_Slope = 0.40 (variabilité PD interindividuelle, cohérent avec
  #    les CV% de Table 2 : 11-21% rat → ~30-40% en clinique)
  omega_Slope <- 0.40
  # 3. Variabilité des baselines (Table 1 : Neut0 2-7, Plt0 150-500)
  neut0_mean <- base_pars$Neut0;  neut0_sd <- 1.5
  plt0_mean  <- base_pars$Plt0;   plt0_sd  <- 100

  times    <- seq(0, n_cycles * interval_h + 21*24, by = 1)
  gfr_vals <- pmax(20, pmin(150,
               rnorm(n_patients, mean = gfr_mean, sd = gfr_sd)))
  eta_CL    <- rnorm(n_patients, 0, omega_CL)
  eta_Slope <- rnorm(n_patients, 0, omega_Slope)
  neut0_i   <- pmax(1.5, pmin(10, rnorm(n_patients, neut0_mean, neut0_sd)))
  plt0_i    <- pmax(75,  pmin(700, rnorm(n_patients, plt0_mean, plt0_sd)))

  grade_neut <- integer(n_patients)
  grade_plt  <- integer(n_patients)

  for (i in seq_len(n_patients)) {
    # Dose Calvert sur GFR individuel
    dose_i <- auc_target * (gfr_vals[i] + 25)

    # PK avec IIV sur CL (variabilité résiduelle post-Calvert)
    pars_i     <- base_pars
    pars_i$CL  <- base_pars$CL * exp(eta_CL[i])

    # PD : IIV sur sensibilité médicament (Slope)
    pars_i$Slope_MPP <- base_pars$Slope_MPP * exp(eta_Slope[i])
    pars_i$Slope_CMP <- base_pars$Slope_CMP * exp(eta_Slope[i])
    pars_i$Slope_MEP <- base_pars$Slope_MEP * exp(eta_Slope[i])

    # PD : baselines individuelles
    pars_i$Neut0 <- neut0_i[i]
    pars_i$Plt0  <- plt0_i[i]

    pars_i$rate_fun <- make_repeated_infusion(
      dose_mg    = dose_i,
      Tinfu_h    = 1,
      interval_h = interval_h,
      n_cycles   = n_cycles
    )

    # État initial adapté aux baselines individuelles
    state_i        <- init_state
    state_i["Neut"] <- neut0_i[i]
    state_i["Plt"]  <- plt0_i[i]
    # Rééquilibrer les compartiments transit Neut et Plt
    a_Neut <- 3 / pars_i$MTT_Neut
    a_Plt  <- 3 / pars_i$MTT_Plt
    T_Neut <- pars_i$k_circ_Neut * neut0_i[i] / a_Neut
    T_Plt  <- pars_i$k_circ_Plt  * plt0_i[i]  / a_Plt
    state_i["T1_Neut"] <- T_Neut
    state_i["T2_Neut"] <- T_Neut
    state_i["T3_Neut"] <- T_Neut
    T1_Plt <- T_Plt / pars_i$lambda2
    state_i["T1_Plt"] <- T1_Plt
    state_i["T2_Plt"] <- T_Plt
    state_i["T3_Plt"] <- T_Plt

    tryCatch({
      out_i <- as.data.frame(lsoda(
        y        = state_i,
        times    = times,
        func     = pkpd_fornari,
        parms    = pars_i,
        rtol     = 1e-6, atol = 1e-8,
        maxsteps = 100000
      ))
      grade_neut[i] <- nadir_grade_neut(out_i$Neut)
      grade_plt[i]  <- nadir_grade_plt(out_i$Plt)
    }, error = function(e) {
      grade_neut[i] <<- NA
      grade_plt[i]  <<- NA
    })

    if (i %% 100 == 0)
      cat(sprintf("    %d/%d patients\n", i, n_patients))
  }

  # % par grade (1 à 4)
  pct_neut <- sapply(1:4, function(g)
    100 * mean(grade_neut == g, na.rm = TRUE))
  pct_plt  <- sapply(1:4, function(g)
    100 * mean(grade_plt  == g, na.rm = TRUE))

  cat(sprintf("  Neutropénie  : G1=%.1f%%  G2=%.1f%%  G3=%.1f%%  G4=%.1f%%\n",
              pct_neut[1], pct_neut[2], pct_neut[3], pct_neut[4]))
  cat(sprintf("  Thrombopénie : G1=%.1f%%  G2=%.1f%%  G3=%.1f%%  G4=%.1f%%\n",
              pct_plt[1],  pct_plt[2],  pct_plt[3],  pct_plt[4]))

  # ── Barplot Figure 4c ─────────────────────────────────
  df_plot <- data.frame(
    Grade    = rep(paste0("G", 1:4), 2),
    Pct      = c(pct_neut, pct_plt),
    Lineage  = rep(c("Neutropenia", "Thrombocytopenia"), each = 4),
    stringsAsFactors = FALSE
  )
  df_plot$Grade   <- factor(df_plot$Grade,   levels = paste0("G", 1:4))
  df_plot$Lineage <- factor(df_plot$Lineage,
                            levels = c("Neutropenia", "Thrombocytopenia"))

  p <- ggplot(df_plot, aes(x = Grade, y = Pct, fill = Grade)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", Pct)),
              vjust = -0.4, size = 3) +
    facet_wrap(~ Lineage) +
    scale_fill_manual(values = unname(GRADE_COLORS), guide = "none") +
    scale_y_continuous(limits = c(0, max(df_plot$Pct) * 1.25),
                       labels = function(x) paste0(x, "%")) +
    labs(title = titre,
         subtitle = sprintf("AUC=%.0f  Q%dD×%d  n=%d patients  GFR~N(%g,%g²)",
                            auc_target, interval_h/24, n_cycles,
                            n_patients, gfr_mean, gfr_sd),
         x = "NCI-CTCAE v5.0 Grade",
         y = "% patients") +
    theme_bw(base_size = 10) +
    theme(strip.text      = element_text(face = "bold"),
          plot.title      = element_text(face = "bold", hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5, size = 8,
                                         color = "grey40"),
          panel.grid.minor = element_blank())

  pdf(file, width = width, height = height)
  print(p)
  dev.off()
  message("✓ Figure 4c sauvegardée : ", file)

  invisible(list(pct_neut = pct_neut, pct_plt = pct_plt,
                 grade_neut = grade_neut, grade_plt = grade_plt))
}
