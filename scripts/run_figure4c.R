############################################################
# run_figure4c.R
# Figure 4c de Fornari 2019 — Population simulation
#
# Selon Supp. S12 :
#   Les barres "Friberg" = simulation Friberg (Schmitt 2010)
#   avec la même PK que S11 (AUC=5, GFR=125, dose fixe).
#   Les barres "QSP" = modèle Fornari (IC50-scaled Slopes).
#
# Ce script génère les deux distributions côte à côte.
############################################################
setwd("/home/user/Hema/scripts")

source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")
source("pkpd_model_FORNARI.R")
source("plots_grades.R")
source("friberg_schmitt.R")

library(ggplot2)
library(gridExtra)
library(grid)
library(deSolve)

# ══════════════════════════════════════════════════════════
# PARAMÈTRES COMMUNS
# ══════════════════════════════════════════════════════════
N_PAT      <- 1000
AUC_TARGET <- 5
GFR_FIXED  <- 125
INTERVAL_H <- 21 * 24
N_CYCLES   <- 2   # Fornari 2019 Fig4c : "Two cycles of 21 days" (légende Fig4)
SEED       <- 42

cat("══════════════════════════════════════════════════════\n")
cat("Figure 4c — Population simulation (n=", N_PAT, ")\n")
cat("  AUC=", AUC_TARGET, " GFR=", GFR_FIXED, " 2 cycles (Q21D)\n")
cat("══════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════
# 1. SIMULATION FRIBERG (Schmitt 2010 via Supp. S12)
# ══════════════════════════════════════════════════════════
cat("── 1. Friberg model (Schmitt 2010) ──\n")
res_friberg <- simulate_friberg_population(
  n_patients  = N_PAT,
  auc_target  = AUC_TARGET,
  gfr_fixed   = GFR_FIXED,
  n_cycles    = N_CYCLES,
  interval_h  = INTERVAL_H,
  seed        = SEED,
  verbose     = TRUE
)

# ══════════════════════════════════════════════════════════
# 2. SIMULATION QSP FORNARI (IC50-scaled Slopes)
# ══════════════════════════════════════════════════════════
cat("\n── 2. QSP Fornari (IC50-scaled Slopes, Supp. S11) ──\n")
res_qsp <- save_grade_figure4c(
  base_pars  = init_pars,
  init_state = init_state,
  auc_target = AUC_TARGET,
  n_cycles   = N_CYCLES,
  interval_h = INTERVAL_H,
  n_patients = N_PAT,
  gfr_fixed  = GFR_FIXED,
  file       = "results_HUMAN/Figure4c_QSP.pdf",
  titre      = "Fornari QSP (IC50-scaled) — Figure 4c",
  seed       = SEED,
  width      = 8,
  height     = 5
)

# ══════════════════════════════════════════════════════════
# 3. FIGURE COMPARATIVE (style Figure 4c)
# ══════════════════════════════════════════════════════════
grades <- paste0("G", 1:4)

df <- data.frame(
  Grade   = rep(grades, 4),
  Model   = rep(c("Friberg\n(Schmitt 2010)", "QSP Fornari\n(IC50-scaled)"), each = 8),
  Lineage = rep(c("Neutropenia", "Thrombocytopenia",
                  "Neutropenia", "Thrombocytopenia"), each = 4),
  Pct     = c(
    res_friberg$pct_neut,
    res_friberg$pct_plt,
    res_qsp$pct_neut,
    res_qsp$pct_plt
  ),
  stringsAsFactors = FALSE
)
df$Grade   <- factor(df$Grade, levels = grades)
df$Model   <- factor(df$Model, levels = c("Friberg\n(Schmitt 2010)",
                                          "QSP Fornari\n(IC50-scaled)"))
df$Lineage <- factor(df$Lineage, levels = c("Neutropenia", "Thrombocytopenia"))

MODEL_COLORS <- c(
  "Friberg\n(Schmitt 2010)" = "#2E75B6",
  "QSP Fornari\n(IC50-scaled)" = "#E8A838"
)

p <- ggplot(df, aes(x = Grade, y = Pct, fill = Model)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, alpha = 0.90) +
  geom_text(aes(label = sprintf("%.0f%%", Pct)),
            position = position_dodge(width = 0.72),
            vjust = -0.4, size = 3.0) +
  facet_wrap(~ Lineage) +
  scale_fill_manual(values = MODEL_COLORS, name = NULL) +
  scale_y_continuous(limits = c(0, max(df$Pct, na.rm = TRUE) * 1.25),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Figure 4c — Carboplatin AUC=5 Q21D×2 (n=1000)",
    subtitle = "Friberg/Schmitt 2010 vs QSP Fornari 2019  |  GFR=125 | 2 cycles (Fig4 Fornari 2019)",
    x        = "NCI-CTCAE v5.0 Grade (nadir)",
    y        = "% patients"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold", size = 11),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle     = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position   = "bottom",
    legend.text       = element_text(size = 9),
    panel.grid.minor  = element_blank()
  )

if (!dir.exists("results_HUMAN")) dir.create("results_HUMAN")
pdf("results_HUMAN/Figure4c_2cycles.pdf", width = 10, height = 6)
print(p)
dev.off()
message("✓ Figure4c_Friberg_vs_QSP.pdf")

# ══════════════════════════════════════════════════════════
# TABLEAU DE SYNTHÈSE
# ══════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════════════════════\n")
cat("TABLEAU DE SYNTHÈSE — % patients par grade nadir\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("%-28s  %5s %5s %5s %5s\n", "", "G1", "G2", "G3", "G4"))
cat("── Neutropénie ────────────────────────────────────\n")
cat(sprintf("%-28s  %5.1f %5.1f %5.1f %5.1f\n", "Friberg (Schmitt 2010)",
            res_friberg$pct_neut[1], res_friberg$pct_neut[2],
            res_friberg$pct_neut[3], res_friberg$pct_neut[4]))
cat(sprintf("%-28s  %5.1f %5.1f %5.1f %5.1f\n", "QSP Fornari",
            res_qsp$pct_neut[1], res_qsp$pct_neut[2],
            res_qsp$pct_neut[3], res_qsp$pct_neut[4]))
cat("── Thrombocytopénie ───────────────────────────────\n")
cat(sprintf("%-28s  %5.1f %5.1f %5.1f %5.1f\n", "Friberg (Schmitt 2010)",
            res_friberg$pct_plt[1], res_friberg$pct_plt[2],
            res_friberg$pct_plt[3], res_friberg$pct_plt[4]))
cat(sprintf("%-28s  %5.1f %5.1f %5.1f %5.1f\n", "QSP Fornari",
            res_qsp$pct_plt[1], res_qsp$pct_plt[2],
            res_qsp$pct_plt[3], res_qsp$pct_plt[4]))
cat("══════════════════════════════════════════════════════\n")
