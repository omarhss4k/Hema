############################################################
# plot_grades_neutropenie.R
# Graphique distribution des grades de neutropénie
# Modèle T-DXd (N=200) vs FDA DESTINY-Breast01 (n=184)
############################################################
library(deSolve)

setwd("/home/user/Hema/scripts_tdxd")
source("parameters_tdxd_rat.R")
source("parameters_tdxd_human.R")
source("pkpd_tdxd_rat.R")

dir.create("results_PKPD_human", showWarnings = FALSE)

set.seed(42)
N_patients <- 200

# ── Paramètres typiques + Slope_CMP calibré ──────────────
pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL
pars_typ   <- c(pars_pd_hu, tdxd_pars_hu)
pars_typ$Slope_CMP <- Slope_CMP_tdxd_human   # 12.0

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)
times_hu    <- seq(0, 126 * 24, by = 12)

omega_CL    <- 0.35
omega_V1    <- 0.20
omega_SCMP  <- 0.33
omega_SMEP  <- 0.33

ctcae_neut <- function(neut) {
  if      (neut < 0.5) "G4"
  else if (neut < 1.0) "G3"
  else if (neut < 1.5) "G2"
  else if (neut < 2.0) "G1"
  else                  "G0"
}

# RBC0 = 5000 × 10⁹/L (Fornari humain) → proxy Hgb
# Seuils Hgb → fraction RBC : G1<10g/dL→<83%, G2<9→<75%, G3<8→<67%, G4<6.5→<54%
ctcae_anemia <- function(rbc, rbc0) {
  frac <- rbc / rbc0
  if      (frac < 0.54) "G4"
  else if (frac < 0.67) "G3"
  else if (frac < 0.75) "G2"
  else if (frac < 0.83) "G1"
  else                   "G0"
}

# ── Simulation population ─────────────────────────────────
cat(sprintf("Simulation N=%d patients...\n", N_patients))

neut_nadirs <- numeric(N_patients)
rbc_nadirs  <- numeric(N_patients)
grades_neut <- character(N_patients)
grades_ane  <- character(N_patients)
n_ok        <- 0

for (i in 1:N_patients) {
  pars_i <- pars_typ
  pars_i$CL_ADC    <- pars_typ$CL_ADC  * exp(rnorm(1,0,omega_CL))
  pars_i$V1_ADC    <- pars_typ$V1_ADC  * exp(rnorm(1,0,omega_V1))
  pars_i$Slope_CMP <- pars_typ$Slope_CMP * exp(rnorm(1,0,omega_SCMP))
  pars_i$Slope_MEP <- pars_typ$Slope_MEP * exp(rnorm(1,0,omega_SMEP))
  pars_i$rate_fun  <- make_tdxd_infusion(5.4,70,1.5,21*24,6)

  out <- tryCatch(
    suppressMessages(suppressWarnings(as.data.frame(lsoda(
      y=state0_hu, times=times_hu, func=pkpd_tdxd_fornari, parms=pars_i,
      rtol=1e-4, atol=1e-6, maxsteps=500000)))),
    error=function(e) NULL)

  if (!is.null(out) && nrow(out) > 10) {
    n_ok <- n_ok + 1
    neut_nadirs[i] <- min(out$Neut, na.rm=TRUE)
    rbc_nadirs[i]  <- min(out$RBC,  na.rm=TRUE)
    grades_neut[i] <- ctcae_neut(neut_nadirs[i])
    grades_ane[i]  <- ctcae_anemia(rbc_nadirs[i], pars_typ$RBC0)
  }
}

# Filtrer les NA
ok_idx      <- grades_neut != ""
grades_ok   <- grades_neut[ok_idx]
grades_ane_ok <- grades_ane[ok_idx]
nadirs_ok   <- neut_nadirs[ok_idx]
rbc_ok      <- rbc_nadirs[ok_idx]

cat(sprintf("Patients simulés avec succès : %d/%d\n\n", n_ok, N_patients))

# ── Comptages ─────────────────────────────────────────────
grade_order <- c("G0","G1","G2","G3","G4")
tab  <- table(factor(grades_ok, levels=grade_order))
pct  <- as.numeric(100 * tab / n_ok)

tab_ane  <- table(factor(grades_ane_ok, levels=grade_order))
pct_ane  <- as.numeric(100 * tab_ane / n_ok)

# Données FDA DESTINY-Breast01 (BLA 761139, n=184)
# Source : FDA Clinical Review — Grades CTCAE v5.0
pct_fda <- c(71.2, 6.5, 7.1, 12.5, 2.7)   # G0, G1, G2, G3, G4

cat("═══════════════════════════════════════════════════════\n")
cat(sprintf("  Distribution grades neutropénie — N=%d patients\n", n_ok))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("  %-6s  %-12s  %-12s  %s\n",
            "Grade", "Modèle (%)", "FDA (%)", "Critère (×10⁹/L)"))
cat(sprintf("  %-6s  %-12s  %-12s  %s\n",
            "─────","──────────","──────────","────────────────"))
criteres <- c("≥ 2.0 (normal)", "1.5–2.0", "1.0–1.5", "0.5–1.0", "< 0.5 (sévère)")
for (j in seq_along(grade_order)) {
  cat(sprintf("  %-6s  %-12.1f  %-12.1f  %s\n",
              grade_order[j], pct[j], pct_fda[j], criteres[j]))
}
cat(sprintf("\n  G3-4 total : Modèle=%.1f%%  FDA=%.1f%%\n",
            pct[4]+pct[5], pct_fda[4]+pct_fda[5]))
cat("═══════════════════════════════════════════════════════\n")

# ══════════════════════════════════════════════════════════
# FIGURE : Barplot double (modèle vs FDA)
# ══════════════════════════════════════════════════════════
pdf("results_PKPD_human/grades_neutropenie_modele_vs_FDA.pdf",
    width = 8, height = 6)

par(mar = c(5, 5, 4, 2), bg = "white")

# Couleurs par grade (bleu=léger → rouge=sévère)
grade_cols <- c(
  G0 = "#4393c3",   # bleu (normal)
  G1 = "#92c5de",   # bleu clair
  G2 = "#f4a582",   # orange
  G3 = "#d6604d",   # rouge-orange
  G4 = "#b2182b"    # rouge foncé
)

mat   <- rbind(pct, pct_fda)
cols  <- rep(grade_cols, each=2)

bp <- barplot(
  mat,
  beside    = TRUE,
  col       = c("#2c7bb6","#d7191c"),   # bleu=modèle, rouge=FDA
  names.arg = grade_order,
  ylim      = c(0, 85),
  main      = "Distribution des grades de neutropénie\nT-DXd 5.4 mg/kg Q3W × 6 cycles",
  xlab      = "Grade CTCAE v5.0  (Neutrophiles × 10⁹/L)",
  ylab      = "Patients (%)",
  cex.main  = 1.15,
  cex.lab   = 1.0,
  cex.axis  = 0.95,
  cex.names = 1.05,
  border    = "white",
  space     = c(0.1, 0.6)
)

# Seuils en µg/mL sur l'axe x (annotations)
seuils <- c("≥2.0\n(normal)", "1.5–2.0", "1.0–1.5", "0.5–1.0", "<0.5\n(sévère)")
mtext(seuils, side=1, at=colMeans(bp), line=2.8, cex=0.72, col="grey40")

# Valeurs au-dessus des barres
for (j in 1:5) {
  text(bp[1,j], pct[j]     + 1.5, sprintf("%.0f%%", pct[j]),     cex=0.85, col="#2c7bb6", font=2)
  text(bp[2,j], pct_fda[j] + 1.5, sprintf("%.0f%%", pct_fda[j]), cex=0.85, col="#d7191c", font=2)
}

# Légende
legend("topright",
       legend = c(sprintf("Modèle T-DXd (N=%d)", n_ok),
                  "FDA DESTINY-Breast01 (n=184)"),
       fill   = c("#2c7bb6","#d7191c"),
       border = "white",
       bty    = "n",
       cex    = 0.95)

# Ligne G3-4 total
g34_mod <- pct[4]+pct[5]
g34_fda <- pct_fda[4]+pct_fda[5]
mtext(sprintf("G3-4 total — Modèle : %.1f%%   FDA : %.1f%%",
              g34_mod, g34_fda),
      side=3, line=-0.2, cex=0.85, col="grey30")

# Boîte annotation Slope calibré
legend("topleft",
       legend = c(sprintf("Slope_CMP = %.1f (calibré)", Slope_CMP_tdxd_human),
                  sprintf("IC50_ADC = %.1f µg/mL (PFB-10)", pars_typ$IC50_ADC_ugmL)),
       bty = "n", cex = 0.82, text.col = "grey40")

dev.off()
cat("  -> results_PKPD_human/grades_neutropenie_modele_vs_FDA.pdf\n")

# ── Figure 2 : Distribution continue du nadir ─────────────
pdf("results_PKPD_human/nadir_neutrophiles_distribution.pdf",
    width = 7, height = 5)

par(mar = c(5, 5, 3.5, 2), bg = "white")

hist(nadirs_ok,
     breaks  = 25,
     col     = "#92c5de",
     border  = "white",
     main    = "Distribution du nadir de neutrophiles\nT-DXd 5.4 mg/kg Q3W × 6 cycles",
     xlab    = "Neut nadir [× 10⁹/L]",
     ylab    = "Fréquence",
     cex.main = 1.1,
     xlim    = c(0, max(nadirs_ok, na.rm=TRUE) * 1.05))

# Seuils CTCAE
abline(v = c(2.0, 1.5, 1.0, 0.5),
       lty = 2, lwd = 1.5,
       col = c("grey50","#f4a582","#d6604d","#b2182b"))

usr <- par("usr")
text(c(1.75, 1.25, 0.75, 0.25), usr[4]*0.92,
     c("G1","G2","G3","G4"),
     col = c("grey50","#f4a582","#d6604d","#b2182b"),
     cex = 0.9, font = 2)
text(2.2, usr[4]*0.92, "G0", col="grey50", cex=0.9, font=2)

# Médiane
med_nadir <- median(nadirs_ok, na.rm=TRUE)
abline(v = med_nadir, col="#2c7bb6", lwd=2.5)
text(med_nadir + 0.05, usr[4]*0.80,
     sprintf("Médiane\n%.2f×10⁹/L", med_nadir),
     col="#2c7bb6", cex=0.85, adj=0)

# Baseline
abline(v = pars_typ$Neut0, col="darkgreen", lwd=2, lty=3)
text(pars_typ$Neut0 - 0.05, usr[4]*0.70,
     sprintf("Baseline\n%.1f", pars_typ$Neut0),
     col="darkgreen", cex=0.85, adj=1)

dev.off()
cat("  -> results_PKPD_human/nadir_neutrophiles_distribution.pdf\n")

# ══════════════════════════════════════════════════════════
# FIGURE 3 : Grades anémie — modèle vs FDA
# ══════════════════════════════════════════════════════════

# Données FDA DESTINY-Breast01 — anémie (proxy Hgb)
# Source : FDA BLA 761139 Safety Review (CTCAE v5.0)
pct_fda_ane <- c(30.4, 37.0, 22.8, 8.2, 1.6)   # G0, G1, G2, G3, G4

cat("\n═══════════════════════════════════════════════════════\n")
cat(sprintf("  Distribution grades anémie — N=%d patients\n", n_ok))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("  %-6s  %-12s  %-12s  %s\n",
            "Grade", "Modèle (%)", "FDA (%)", "Critère (Hgb proxy)"))
cat(sprintf("  %-6s  %-12s  %-12s  %s\n",
            "─────","──────────","──────────","────────────────"))
crit_ane <- c("≥83% RBC0 (normal)", "75–83% RBC0", "67–75% RBC0",
              "54–67% RBC0", "<54% RBC0 (sévère)")
for (j in seq_along(grade_order)) {
  cat(sprintf("  %-6s  %-12.1f  %-12.1f  %s\n",
              grade_order[j], pct_ane[j], pct_fda_ane[j], crit_ane[j]))
}
cat(sprintf("\n  G3-4 total : Modèle=%.1f%%  FDA=%.1f%%\n",
            pct_ane[4]+pct_ane[5], pct_fda_ane[4]+pct_fda_ane[5]))
cat("═══════════════════════════════════════════════════════\n")

pdf("results_PKPD_human/grades_anemie_modele_vs_FDA.pdf",
    width = 8, height = 6)

par(mar = c(5, 5, 4, 2), bg = "white")

mat_ane <- rbind(pct_ane, pct_fda_ane)

bp2 <- barplot(
  mat_ane,
  beside    = TRUE,
  col       = c("#d95f02","#1b9e77"),   # orange=modèle, vert=FDA
  names.arg = grade_order,
  ylim      = c(0, 85),
  main      = "Distribution des grades d'anémie\nT-DXd 5.4 mg/kg Q3W × 6 cycles",
  xlab      = "Grade CTCAE v5.0  (proxy Hémoglobine via RBC)",
  ylab      = "Patients (%)",
  cex.main  = 1.15,
  cex.lab   = 1.0,
  cex.axis  = 0.95,
  cex.names = 1.05,
  border    = "white",
  space     = c(0.1, 0.6)
)

seuils_ane <- c("≥83%\nRBC0\n(normal)", "75–83%\nRBC0", "67–75%\nRBC0",
                "54–67%\nRBC0", "<54%\nRBC0\n(sévère)")
mtext(seuils_ane, side=1, at=colMeans(bp2), line=3.2, cex=0.68, col="grey40")

for (j in 1:5) {
  text(bp2[1,j], pct_ane[j]     + 1.5, sprintf("%.0f%%", pct_ane[j]),
       cex=0.85, col="#d95f02", font=2)
  text(bp2[2,j], pct_fda_ane[j] + 1.5, sprintf("%.0f%%", pct_fda_ane[j]),
       cex=0.85, col="#1b9e77", font=2)
}

legend("topright",
       legend = c(sprintf("Modèle T-DXd (N=%d)", n_ok),
                  "FDA DESTINY-Breast01 (n=184)"),
       fill   = c("#d95f02","#1b9e77"),
       border = "white",
       bty    = "n",
       cex    = 0.95)

g34_mod_ane <- pct_ane[4]+pct_ane[5]
g34_fda_ane <- pct_fda_ane[4]+pct_fda_ane[5]
mtext(sprintf("G3-4 total — Modèle : %.1f%%   FDA : %.1f%%",
              g34_mod_ane, g34_fda_ane),
      side=3, line=-0.2, cex=0.85, col="grey30")

legend("topleft",
       legend = c(sprintf("Slope_MEP = %.2f (Fornari carboplatin)", pars_typ$Slope_MEP),
                  sprintf("IC50_ADC ery = %.1f µg/mL (PFB-10)", pars_typ$IC50_ADC_ery_ugmL)),
       bty = "n", cex = 0.82, text.col = "grey40")

dev.off()
cat("  -> results_PKPD_human/grades_anemie_modele_vs_FDA.pdf\n")

# ── Figure 4 : Distribution nadir RBC ────────────────────
rbc_pct_drop <- 100 * (rbc_ok - pars_typ$RBC0) / pars_typ$RBC0

pdf("results_PKPD_human/nadir_rbc_distribution.pdf",
    width = 7, height = 5)

par(mar = c(5, 5, 3.5, 2), bg = "white")

hist(rbc_pct_drop,
     breaks  = 25,
     col     = "#fdae6b",
     border  = "white",
     main    = "Distribution du nadir de GR (%Δ baseline)\nT-DXd 5.4 mg/kg Q3W × 6 cycles",
     xlab    = "RBC nadir [%Δ vs baseline]",
     ylab    = "Fréquence",
     cex.main = 1.1)

# Seuils anémie (% chute)
seuils_drop <- -100*c(1-0.83, 1-0.75, 1-0.67, 1-0.54)  # -17, -25, -33, -46 %
abline(v = seuils_drop, lty=2, lwd=1.5,
       col = c("grey50","#f4a582","#d6604d","#b2182b"))

usr2 <- par("usr")
text(seuils_drop + 0.5, usr2[4]*0.9,
     c("G1","G2","G3","G4"),
     col = c("grey50","#f4a582","#d6604d","#b2182b"),
     cex = 0.9, font=2, adj=0)
text(-2, usr2[4]*0.9, "G0", col="grey50", cex=0.9, font=2)

med_rbc <- median(rbc_pct_drop, na.rm=TRUE)
abline(v = med_rbc, col="#d95f02", lwd=2.5)
text(med_rbc - 0.5, usr2[4]*0.78,
     sprintf("Médiane\n%.1f%%", med_rbc),
     col="#d95f02", cex=0.85, adj=1)

dev.off()
cat("  -> results_PKPD_human/nadir_rbc_distribution.pdf\n")
