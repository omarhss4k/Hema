############################################################
# compare_options.R
# Comparaison 3 options pour reproduire Figure 4c (Fornari 2019)
#
# Option 1 : Dose fixe (GFR_ref=125 → 750 mg), IIV sur CL + Slopes
# Option 2 : GFR ~ N(125, 25²) → dose moyenne = 750 mg, IIV sur CL + Slopes
# Option 3 : GFR ~ N(78, 20²) + omega_Slope très forte (1.20)
#
# Cible Fornari (barres jaunes IC50-scaled, Figure 4c) :
#   Neutropénie  : G1≈12%  G2≈15%  G3≈15%  G4≈1%
#   Thrombopénie : G1≈31%  G2≈11%  G3≈10%  G4≈4%
############################################################
setwd("/home/user/Hema/scripts")

source("pkpd_model_FORNARI.R")
source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")
source("plots_grades.R")

library(ggplot2)
library(dplyr)
library(deSolve)

N_PAT      <- 500    # patients par option (réduit pour la vitesse)
AUC_TARGET <- 5
INTERVAL_H <- 21 * 24
N_CYCLES   <- 2
GFR_REF    <- 125    # GFR du patient de référence (Supplementary S11)
TIMES      <- seq(0, N_CYCLES * INTERVAL_H + 21 * 24, by = 1)

# ── Cibles Fornari Figure 4c (barres jaunes = IC50-scaled) ─
FORNARI_NEUT <- c(G1 = 12, G2 = 15, G3 = 15, G4 = 1)
FORNARI_PLT  <- c(G1 = 31, G2 = 11, G3 = 10, G4 = 4)

# ════════════════════════════════════════════════════════
# FONCTION COMMUNE : simuler n patients, retourne grades
# ════════════════════════════════════════════════════════
simulate_grades <- function(option_name,
                            gfr_mean   = 78,
                            gfr_sd     = 20,
                            dose_fixe  = NULL,   # si non NULL : dose identique pour tous
                            omega_CL   = 0.35,
                            omega_Slope = 0.55,
                            seed       = 42) {

  cat(sprintf("\n══ Option : %s ══\n", option_name))
  set.seed(seed)

  gfr_vals  <- pmax(20, pmin(200, rnorm(N_PAT, mean = gfr_mean, sd = gfr_sd)))
  eta_CL    <- rnorm(N_PAT, 0, omega_CL)
  eta_Slope <- rnorm(N_PAT, 0, omega_Slope)
  neut0_i   <- pmax(1.0, pmin(10,  rnorm(N_PAT, init_pars$Neut0, 1.5)))
  plt0_i    <- pmax(50,  pmin(700, rnorm(N_PAT, init_pars$Plt0,  100)))

  grade_neut <- integer(N_PAT)
  grade_plt  <- integer(N_PAT)

  a_Plt  <- 3 / init_pars$MTT_Plt
  a_Neut <- 3 / init_pars$MTT_Neut

  for (i in seq_len(N_PAT)) {

    # ── Dose ───────────────────────────────────────────────
    dose_i <- if (!is.null(dose_fixe)) {
      dose_fixe
    } else {
      AUC_TARGET * (gfr_vals[i] + 25)
    }

    # ── Paramètres individuels ─────────────────────────────
    pars_i              <- init_pars
    pars_i$CL           <- init_pars$CL * exp(eta_CL[i])
    pars_i$Slope_MPP    <- init_pars$Slope_MPP * exp(eta_Slope[i])
    pars_i$Slope_CMP    <- init_pars$Slope_CMP * exp(eta_Slope[i])
    pars_i$Slope_MEP    <- init_pars$Slope_MEP * exp(eta_Slope[i])
    pars_i$Neut0        <- neut0_i[i]
    pars_i$Plt0         <- plt0_i[i]
    pars_i$rate_fun     <- make_repeated_infusion(
      dose_mg    = dose_i,
      Tinfu_h    = 1,
      interval_h = INTERVAL_H,
      n_cycles   = N_CYCLES
    )

    # ── État initial adapté aux baselines individuelles ────
    state_i           <- init_state
    state_i["Neut"]   <- neut0_i[i]
    state_i["Plt"]    <- plt0_i[i]
    T_Neut_i          <- pars_i$k_circ_Neut * neut0_i[i] / (3 / pars_i$MTT_Neut)
    T_Plt_i           <- pars_i$k_circ_Plt  * plt0_i[i]  / a_Plt
    state_i["T1_Neut"] <- T_Neut_i
    state_i["T2_Neut"] <- T_Neut_i
    state_i["T3_Neut"] <- T_Neut_i
    state_i["T1_Plt"]  <- T_Plt_i / pars_i$lambda2
    state_i["T2_Plt"]  <- T_Plt_i
    state_i["T3_Plt"]  <- T_Plt_i

    tryCatch({
      out_i <- as.data.frame(lsoda(
        y        = state_i,
        times    = TIMES,
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
      cat(sprintf("  %d/%d patients\n", i, N_PAT))
  }

  pct_neut <- sapply(1:4, function(g) 100 * mean(grade_neut == g, na.rm = TRUE))
  pct_plt  <- sapply(1:4, function(g) 100 * mean(grade_plt  == g, na.rm = TRUE))

  cat(sprintf("  Neut : G1=%.1f%% G2=%.1f%% G3=%.1f%% G4=%.1f%%\n",
              pct_neut[1], pct_neut[2], pct_neut[3], pct_neut[4]))
  cat(sprintf("  Plt  : G1=%.1f%% G2=%.1f%% G3=%.1f%% G4=%.1f%%\n",
              pct_plt[1],  pct_plt[2],  pct_plt[3],  pct_plt[4]))

  list(option = option_name,
       pct_neut = pct_neut,
       pct_plt  = pct_plt,
       grade_neut = grade_neut,
       grade_plt  = grade_plt)
}

# ════════════════════════════════════════════════════════
# OPTION 1 : Dose fixe (GFR=125 → 750 mg), IIV CL + Slopes
# ════════════════════════════════════════════════════════
dose_ref <- AUC_TARGET * (GFR_REF + 25)   # 750 mg
res1 <- simulate_grades(
  option_name = "Option 1\n(Dose fixe 750mg)",
  dose_fixe   = dose_ref,
  omega_CL    = 0.35,
  omega_Slope = 0.55,
  seed        = 42
)

# ════════════════════════════════════════════════════════
# OPTION 2 : GFR ~ N(125, 25²) → dose moyenne ≈ 750 mg
# ════════════════════════════════════════════════════════
res2 <- simulate_grades(
  option_name = "Option 2\n(GFR~N(125,25))",
  gfr_mean    = 125,
  gfr_sd      = 25,
  omega_CL    = 0.35,
  omega_Slope = 0.55,
  seed        = 42
)

# ════════════════════════════════════════════════════════
# OPTION 3 : GFR ~ N(78, 20²) + omega_Slope = 1.20
# ════════════════════════════════════════════════════════
res3 <- simulate_grades(
  option_name  = "Option 3\n(omega_Slope=1.20)",
  gfr_mean     = 78,
  gfr_sd       = 20,
  omega_CL     = 0.35,
  omega_Slope  = 1.20,
  seed         = 42
)

# ════════════════════════════════════════════════════════
# FIGURE COMPARATIVE — 3 options + cible Fornari
# ════════════════════════════════════════════════════════
make_comparison_figure <- function(res_list, out_file) {

  grades  <- paste0("G", 1:4)
  options <- sapply(res_list, function(r) r$option)

  # ── Data.frame long ──────────────────────────────────
  rows <- list()
  for (r in res_list) {
    for (g in 1:4) {
      rows[[length(rows) + 1]] <- data.frame(
        Grade    = grades[g],
        Option   = r$option,
        Neut_pct = r$pct_neut[g],
        Plt_pct  = r$pct_plt[g],
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, rows)

  # Ajouter la cible Fornari
  df_fornari <- data.frame(
    Grade    = grades,
    Option   = "Cible Fornari\n(IC50-scaled)",
    Neut_pct = as.numeric(FORNARI_NEUT),
    Plt_pct  = as.numeric(FORNARI_PLT),
    stringsAsFactors = FALSE
  )
  df_all <- rbind(df_fornari, df)

  df_all$Grade  <- factor(df_all$Grade, levels = grades)
  df_all$Option <- factor(df_all$Option,
                           levels = c("Cible Fornari\n(IC50-scaled)",
                                      unique(sapply(res_list, `[[`, "option"))))

  # Couleurs
  option_colors <- c(
    "Cible Fornari\n(IC50-scaled)"      = "#2E75B6",
    "Option 1\n(Dose fixe 750mg)"       = "#ED7D31",
    "Option 2\n(GFR~N(125,25))"         = "#70AD47",
    "Option 3\n(omega_Slope=1.20)"      = "#9B59B6"
  )

  # ── Neut ──────────────────────────────────────────────
  p_neut <- ggplot(df_all, aes(x = Grade, y = Neut_pct, fill = Option)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.88) +
    geom_text(aes(label = sprintf("%.0f", Neut_pct)),
              position = position_dodge(width = 0.75),
              vjust = -0.4, size = 2.8) +
    scale_fill_manual(values = option_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 40), labels = function(x) paste0(x, "%")) +
    labs(title = "Neutropénie — % patients par grade nadir",
         x = "Grade NCI-CTCAE", y = "% patients") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom",
          legend.text   = element_text(size = 8),
          panel.grid.minor = element_blank())

  # ── Plt ───────────────────────────────────────────────
  p_plt <- ggplot(df_all, aes(x = Grade, y = Plt_pct, fill = Option)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.88) +
    geom_text(aes(label = sprintf("%.0f", Plt_pct)),
              position = position_dodge(width = 0.75),
              vjust = -0.4, size = 2.8) +
    scale_fill_manual(values = option_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 50), labels = function(x) paste0(x, "%")) +
    labs(title = "Thrombocytopénie — % patients par grade nadir",
         x = "Grade NCI-CTCAE", y = "% patients") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom",
          legend.text   = element_text(size = 8),
          panel.grid.minor = element_blank())

  pdf(out_file, width = 13, height = 7)
  gridExtra::grid.arrange(
    p_neut, p_plt,
    ncol = 2,
    top  = grid::textGrob(
      sprintf("Comparaison des 3 options — Carboplatine AUC=%d Q21D×%d (%d patients/option)",
              AUC_TARGET, N_CYCLES, N_PAT),
      gp = grid::gpar(fontface = "bold", fontsize = 13)
    )
  )
  dev.off()
  message("✓ Figure comparative sauvegardée : ", out_file)
}

# ════════════════════════════════════════════════════════
# FIGURES INDIVIDUELLES PAR OPTION (style Figure 4c)
# ════════════════════════════════════════════════════════
make_individual_figure <- function(res, out_file, subtitle = "") {
  grades <- paste0("G", 1:4)
  df <- data.frame(
    Grade   = rep(grades, 2),
    Lineage = rep(c("Neutropenia", "Thrombocytopenia"), each = 4),
    Pct     = c(res$pct_neut, res$pct_plt),
    stringsAsFactors = FALSE
  )
  df$Grade   <- factor(df$Grade, levels = grades)
  df$Lineage <- factor(df$Lineage, levels = c("Neutropenia", "Thrombocytopenia"))

  GRADE_COLORS <- c(G1 = "#fee08b", G2 = "#fc8d59", G3 = "#d73027", G4 = "#7b0404")

  p <- ggplot(df, aes(x = Grade, y = Pct, fill = Grade)) +
    geom_col(width = 0.65, alpha = 0.88) +
    geom_text(aes(label = sprintf("%.1f%%", Pct)), vjust = -0.4, size = 3) +
    facet_wrap(~ Lineage) +
    scale_fill_manual(values = unname(GRADE_COLORS), guide = "none") +
    scale_y_continuous(limits = c(0, max(df$Pct, 5) * 1.3),
                       labels = function(x) paste0(x, "%")) +
    labs(title    = gsub("\n", " — ", res$option),
         subtitle = subtitle,
         x = "NCI-CTCAE v5.0 Grade", y = "% patients") +
    theme_bw(base_size = 10) +
    theme(strip.text       = element_text(face = "bold"),
          plot.title        = element_text(face = "bold", hjust = 0.5),
          plot.subtitle     = element_text(hjust = 0.5, size = 8, color = "grey40"),
          panel.grid.minor  = element_blank())

  pdf(out_file, width = 8, height = 5)
  print(p)
  dev.off()
  message("✓ Sauvegardé : ", out_file)
}

if (!dir.exists("results_HUMAN")) dir.create("results_HUMAN")

make_individual_figure(res1, "results_HUMAN/Compare_Option1_doseFixe.pdf",
  subtitle = sprintf("Dose fixe = %.0f mg (AUC=5, GFR_ref=125), omega_CL=0.35, omega_Slope=0.55, n=%d",
                     dose_ref, N_PAT))

make_individual_figure(res2, "results_HUMAN/Compare_Option2_GFR125.pdf",
  subtitle = sprintf("GFR~N(125,25), omega_CL=0.35, omega_Slope=0.55, n=%d", N_PAT))

make_individual_figure(res3, "results_HUMAN/Compare_Option3_slopeLarge.pdf",
  subtitle = sprintf("GFR~N(78,20), omega_CL=0.35, omega_Slope=1.20, n=%d", N_PAT))

make_comparison_figure(
  list(res1, res2, res3),
  "results_HUMAN/Compare_3options_vs_Fornari.pdf"
)

# ════════════════════════════════════════════════════════
# TABLEAU DE SYNTHÈSE
# ════════════════════════════════════════════════════════
cat("\n\n══════════════════════════════════════════════════════\n")
cat("TABLEAU DE SYNTHÈSE — % patients par grade (nadir)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("%-32s  %5s  %5s  %5s  %5s\n", "", "G1", "G2", "G3", "G4"))
cat("── Neutropénie ────────────────────────────────────\n")
cat(sprintf("%-32s  %5.1f  %5.1f  %5.1f  %5.1f\n",
            "Cible Fornari (IC50-scaled)",
            FORNARI_NEUT["G1"], FORNARI_NEUT["G2"],
            FORNARI_NEUT["G3"], FORNARI_NEUT["G4"]))
for (r in list(res1, res2, res3))
  cat(sprintf("%-32s  %5.1f  %5.1f  %5.1f  %5.1f\n",
              gsub("\n", " ", r$option),
              r$pct_neut[1], r$pct_neut[2], r$pct_neut[3], r$pct_neut[4]))
cat("── Thrombocytopénie ───────────────────────────────\n")
cat(sprintf("%-32s  %5.1f  %5.1f  %5.1f  %5.1f\n",
            "Cible Fornari (IC50-scaled)",
            FORNARI_PLT["G1"], FORNARI_PLT["G2"],
            FORNARI_PLT["G3"], FORNARI_PLT["G4"]))
for (r in list(res1, res2, res3))
  cat(sprintf("%-32s  %5.1f  %5.1f  %5.1f  %5.1f\n",
              gsub("\n", " ", r$option),
              r$pct_plt[1], r$pct_plt[2], r$pct_plt[3], r$pct_plt[4]))

# Score RMSE vs Fornari (sur les 8 grades G1-G4 x 2 lignées)
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
cat("\n── Score RMSE vs Fornari (plus petit = mieux) ─────\n")
for (r in list(res1, res2, res3)) {
  sc <- rmse(c(FORNARI_NEUT, FORNARI_PLT),
             c(r$pct_neut, r$pct_plt))
  cat(sprintf("%-32s  RMSE = %.2f %%\n", gsub("\n", " ", r$option), sc))
}
cat("══════════════════════════════════════════════════════\n")
