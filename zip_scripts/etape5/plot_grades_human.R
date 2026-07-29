############################################################
# plot_grades_human.R
# Graphique grades CTCAE -- Neutropénie + Anémie
# T-DXd 5.4 mg/kg Q3W × 6 -- N=300 patients
# Charge depuis RDS si disponible, sinon re-simule
############################################################

rds_path <- "results/population_results.rds"

if (!exists("results") || nrow(results) == 0) {
  if (file.exists(rds_path)) {
    results <- readRDS(rds_path)
    n_ok    <- nrow(results)
  } else {
    stop("Lancez d'abord run_pkpd_tdxd_human_population.R")
  }
}

grade_order <- c("G0", "G1", "G2", "G3", "G4")
n_ok <- nrow(results)

tab_neut   <- table(factor(results$Grade_Neut,   levels = grade_order))
tab_anemia <- table(factor(results$Grade_Anemia, levels = grade_order))
pct_n <- as.numeric(round(100 * tab_neut   / n_ok, 1))
pct_a <- as.numeric(round(100 * tab_anemia / n_ok, 1))
names(pct_n) <- grade_order
names(pct_a) <- grade_order

# FDA BLA 761139 -- DESTINY-Breast01 (U201, n=184)
fda_neut  <- c(G0=71, G1=7, G2=7, G3=13, G4=3)
# FDA BLA 761139 -- anémie (lab grading, n=184)
fda_anemia <- c(G0=30, G1=37, G2=24, G3=8, G4=1)

# Palette par grade
grade_cols <- c(G0="#4575b4", G1="#91bfdb", G2="#fee090", G3="#fc8d59", G4="#d73027")

pdf("results_PKPD_human/grades_neut_anemia.pdf", width = 12, height = 6)

layout(matrix(1:2, nrow = 1))
par(mar = c(5, 5, 4, 1.5), mgp = c(3.2, 0.8, 0))

plot_grade_panel <- function(pct_mod, pct_fda, title, fda_label = "FDA (DESTINY-B01, n=184)") {
  xs   <- seq_along(grade_order)
  xw   <- 0.35
  xgap <- 0.12

  plot(NA, xlim = c(0.5, 5.5), ylim = c(0, 100),
       xaxt = "n", xlab = "", ylab = "Patients (%)",
       main = title, cex.main = 1.3, cex.axis = 1.1, cex.lab = 1.15,
       las = 1)
  axis(1, at = xs, labels = grade_order, cex.axis = 1.2, font = 2)
  abline(h = seq(0, 100, 20), col = "grey90", lwd = 0.8)

  for (gi in seq_along(grade_order)) {
    g   <- grade_order[gi]
    col <- grade_cols[g]

    # Barre modèle (pleine)
    xl <- xs[gi] - xgap/2 - xw
    xr <- xs[gi] - xgap/2
    rect(xl, 0, xr, pct_mod[g], col = col, border = "white", lwd = 0.5)
    if (pct_mod[g] >= 3)
      text((xl+xr)/2, pct_mod[g] + 2.5,
           sprintf("%.0f%%", pct_mod[g]), cex = 0.95, font = 2, col = col)

    # Barre FDA (hachurée)
    xl2 <- xs[gi] + xgap/2
    xr2 <- xs[gi] + xgap/2 + xw
    rect(xl2, 0, xr2, pct_fda[g], col = NA, border = col, lwd = 2.5, lty = 1)
    rect(xl2, 0, xr2, pct_fda[g],
         density = 18, angle = 45, col = col, border = col, lwd = 0.5)
    if (pct_fda[g] >= 3)
      text((xl2+xr2)/2, pct_fda[g] + 2.5,
           sprintf("%.0f%%", pct_fda[g]), cex = 0.95, col = col)
  }

  legend("topright",
         legend = c(sprintf("Modèle (N=%d)", n_ok), fda_label),
         fill   = c("grey60", NA),
         border = c(NA, "grey40"),
         density = c(NA, 18),
         angle   = c(NA, 45),
         bty = "n", cex = 1.05)

  # Annotation G3-4
  g34_mod <- pct_mod["G3"] + pct_mod["G4"]
  g34_fda <- pct_fda["G3"] + pct_fda["G4"]
  mtext(sprintf("G3-4 : modèle %.1f%%  |  FDA %.1f%%",
                g34_mod, g34_fda),
        side = 1, line = 3.8, cex = 1.0, col = "grey30")
}

plot_grade_panel(pct_n, fda_neut,
                 "Neutropénie -- T-DXd 5.4 mg/kg Q3W × 6")
plot_grade_panel(pct_a, fda_anemia,
                 "Anémie (proxy RBC) -- T-DXd 5.4 mg/kg Q3W × 6")

mtext("T-DXd 5.4 mg/kg Q3W × 6 cycles  |  Grades CTCAE v5  |  Modèle vs FDA BLA 761139",
      outer = TRUE, line = -1.2, cex = 0.95, col = "grey40")

dev.off()
cat("  -> results_PKPD_human/grades_neut_anemia.pdf\n")
