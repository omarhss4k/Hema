############################################################
# run_single_dose_human.R
# S11 — Carboplatin PK in human plasma (Supplementary S11)
#
# 2-compartment PK model (Zandvliet et al. 2008)
# Single dose, Calvert formula (Eq. S12):
#   Dose (mg) = AUC_target × (GFR + 25)
# Target AUC(0-24h) = 5 (mg·min)/mL, GFR = 125 mL/min
############################################################
source("pkpd_model_FORNARI.R")
source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")

library(ggplot2)
library(gridExtra)
library(grid)

if (!dir.exists("results_HUMAN")) dir.create("results_HUMAN")

# ── Calvert formula (Eq. S12) ────────────────────────────
AUC_target_mgmin_mL <- 5     # (mg·min)/mL  ← target AUC(0-24h)
GFR_mL_min          <- 125   # mL/min — normal renal function (S11)

dose_mg <- AUC_target_mgmin_mL * (GFR_mL_min + 25)   # 750 mg

cat(sprintf(
  "=== S11 — Single dose Calvert ===\n  AUC_target = %g (mg·min)/mL\n  GFR = %g mL/min\n  Dose = %.0f mg\n\n",
  AUC_target_mgmin_mL, GFR_mL_min, dose_mg
))

# ── PK parameters adjusted for GFR = 125 mL/min ─────────
# Zandvliet 2008: CL = GFR + 25  (Calvert relationship)
#   CL = (125 + 25) mL/min = 150 mL/min = 9.0 L/h
CL_mL_min <- GFR_mL_min + 25                 # 150 mL/min
CL_L_h    <- CL_mL_min * 60 / 1000           # 9.0 L/h

pars_s11      <- init_pars
pars_s11$CL   <- CL_L_h   # override default CL (which uses median GFR=78)

cat(sprintf("  CL adjusted for GFR=%g: %.0f mL/min = %.4f L/h\n",
            GFR_mL_min, CL_mL_min, CL_L_h))
cat(sprintf("  Expected AUC(0-inf) = Dose/CL = %.0f mg / %.4f L/h = %.2f mg·h/L\n",
            dose_mg, CL_L_h, dose_mg / CL_L_h))
cat(sprintf("  = %.2f (mg·min)/mL ← target = %g (mg·min)/mL ✓\n\n",
            dose_mg / CL_L_h * 60 / 1000, AUC_target_mgmin_mL))

# ── Single dose infusion (1-hour, day 0) ─────────────────
Tinfu_h           <- 1
pars_s11$rate_fun <- function(t) {
  if (t >= 0 && t < Tinfu_h) dose_mg / Tinfu_h else 0
}

# ── Time grid ────────────────────────────────────────────
# Fine grid for PK (0–48 h), extended for PD effects (0–42 days)
times_pk_h  <- c(seq(0, 48, by = 0.5),   seq(49, 42*24, by = 1))
times_all   <- sort(unique(times_pk_h))

# ── Simulation ───────────────────────────────────────────
cat("=== Running single-dose PKPD simulation ===\n")
sim_s11 <- simulate_all(times_all, pars_s11, init_state)

# ── Compute numerical AUC(0-24h) ────────────────────────
mask24   <- sim_s11$time <= 24
t24      <- sim_s11$time[mask24]
C24      <- sim_s11$C1[mask24]
auc24_mgL_h    <- sum(diff(t24) * (head(C24, -1) + tail(C24, -1)) / 2)  # mg·h/L
auc24_mgmin_mL <- auc24_mgL_h * 60 / 1000   # (mg·min)/mL  [mg·h/L × 60min/h ÷ 1000mL/L]

cat(sprintf("\n  Numerical AUC(0-24h) = %.2f mg·h/L = %.2f (mg·min)/mL\n",
            auc24_mgL_h, auc24_mgmin_mL))

# ── Plot 1: PK concentration-time profile (0–48 h) ──────
df_pk <- data.frame(
  time_h  = sim_s11$time[sim_s11$time <= 48],
  C1_mgL  = sim_s11$C1[sim_s11$time <= 48],
  C1_uM   = sim_s11$C1[sim_s11$time <= 48] * 1000 / pars_s11$MW_carboplatin
)

p_pk <- ggplot(df_pk, aes(x = time_h)) +
  geom_line(aes(y = C1_mgL), color = "#2166ac", linewidth = 1.1) +
  annotate("rect", xmin = 0, xmax = Tinfu_h,
           ymin = -Inf, ymax = Inf,
           fill = "#fdae61", alpha = 0.25) +
  annotate("text", x = Tinfu_h / 2, y = max(df_pk$C1_mgL) * 0.92,
           label = "Infusion", size = 3, color = "#e08214", hjust = 0.5) +
  annotate("text",
           x = 24, y = max(df_pk$C1_mgL) * 0.55,
           label = sprintf("AUC(0-24h) = %.2f\n(mg·min)/mL", auc24_mgmin_mL),
           size = 3.2, color = "#555555", hjust = 0.5) +
  labs(
    title    = sprintf(
      "Carboplatin plasma PK — Single dose %.0f mg (AUC_target=%g, GFR=%g mL/min)",
      dose_mg, AUC_target_mgmin_mL, GFR_mL_min),
    subtitle = sprintf(
      "2-compartment model (Zandvliet et al. 2008)  |  CL=%.1f L/h  V1=%.2f L  Q=%.2f L/h  V2=%.2f L",
      CL_L_h, pars_s11$V1, pars_s11$Q, pars_s11$V2),
    x = "Time (h)",
    y = "Carboplatin plasma concentration (mg/L)"
  ) +
  scale_x_continuous(breaks = seq(0, 48, by = 6)) +
  theme_bw(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size  = 8, color = "grey40"),
    panel.grid.minor = element_blank()
  )

# ── Plot 2: Damage kinetics ──────────────────────────────
df_dmg <- data.frame(
  time_h  = sim_s11$time[sim_s11$time <= 48],
  Damage  = sim_s11$Damage[sim_s11$time <= 48]
)

p_dmg <- ggplot(df_dmg, aes(x = time_h, y = Damage)) +
  geom_line(color = "#d73027", linewidth = 1.0) +
  labs(
    title = "Platinum-DNA damage kinetics",
    x     = "Time (h)",
    y     = "Damage (a.u.)"
  ) +
  scale_x_continuous(breaks = seq(0, 48, by = 6)) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())

# ── Plot 3: PD — Neutrophils & Platelets (42 days) ──────
df_pd <- data.frame(
  days = sim_s11$days,
  Neut = sim_s11$Neut,
  Plt  = sim_s11$Plt
)

p_neut <- ggplot(df_pd, aes(x = days, y = Neut)) +
  geom_hline(yintercept = pars_s11$Neut0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 1.0, linetype = "dotted",
             color = "#d73027", linewidth = 0.6) +
  geom_line(color = "#d7191c", linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "dotted",
             color = "grey55", linewidth = 0.4) +
  annotate("text", x = 1, y = 1.0 * 1.15, label = "G3 (<1.0)",
           size = 2.8, color = "#d73027", hjust = 0) +
  scale_y_log10(limits = c(0.05, 10),
                breaks  = c(0.1, 0.5, 1, 2, 5, 10),
                labels  = c("0.1", "0.5", "1", "2", "5", "10")) +
  labs(title = "Neutrophils", x = "Time (days)",
       y = expression(10^9 ~ cells ~ L^{-1})) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        panel.grid.minor = element_blank())

p_plt <- ggplot(df_pd, aes(x = days, y = Plt)) +
  geom_hline(yintercept = pars_s11$Plt0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 150, linetype = "dotted",
             color = "#d73027", linewidth = 0.6) +
  geom_line(color = "#1b9e77", linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "dotted",
             color = "grey55", linewidth = 0.4) +
  annotate("text", x = 1, y = 150 * 1.08, label = "G1 (<150)",
           size = 2.8, color = "#d73027", hjust = 0) +
  scale_y_log10(limits = c(10, 600),
                breaks  = c(10, 25, 50, 100, 200, 345, 600),
                labels  = c("10", "25", "50", "100", "200", "345", "600")) +
  labs(title = "Platelets", x = "Time (days)",
       y = expression(10^9 ~ cells ~ L^{-1})) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        panel.grid.minor = element_blank())

# ── Save PDF ──────────────────────────────────────────────
pdf("results_HUMAN/FigureS11_single_dose_PK.pdf", width = 10, height = 5)
grid.arrange(p_pk, p_dmg, ncol = 2,
             top = textGrob(
               "S11 – Carboplatin PK in Human Plasma (Single dose, Calvert formula)",
               gp = gpar(fontface = "bold", fontsize = 12)))
dev.off()
message("✓ Saved: results_HUMAN/FigureS11_single_dose_PK.pdf")

pdf("results_HUMAN/FigureS11_single_dose_PD.pdf", width = 10, height = 5)
grid.arrange(p_neut, p_plt, ncol = 2,
             top = textGrob(
               sprintf("S11 – Hematopoiesis after Single Dose Carboplatin %.0f mg (AUC=%g, GFR=%g mL/min)",
                       dose_mg, AUC_target_mgmin_mL, GFR_mL_min),
               gp = gpar(fontface = "bold", fontsize = 11)))
dev.off()
message("✓ Saved: results_HUMAN/FigureS11_single_dose_PD.pdf")

# ── Summary ───────────────────────────────────────────────
cat("\n=======================================================\n")
cat("S11 — Single dose simulation summary\n")
cat(sprintf("  Dose       : %.0f mg\n",         dose_mg))
cat(sprintf("  GFR        : %g mL/min\n",        GFR_mL_min))
cat(sprintf("  CL         : %.1f L/h\n",         CL_L_h))
cat(sprintf("  AUC(0-24h) : %.2f (mg·min)/mL\n", auc24_mgmin_mL))
cat(sprintf("  Neut nadir : %.3f × 10⁹/L\n",     min(sim_s11$Neut)))
cat(sprintf("  Plt  nadir : %.1f × 10⁹/L\n",     min(sim_s11$Plt)))
cat("Output: results_HUMAN/FigureS11_single_dose_PK.pdf\n")
cat("Output: results_HUMAN/FigureS11_single_dose_PD.pdf\n")
cat("=======================================================\n")
