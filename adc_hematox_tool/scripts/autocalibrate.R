############################################################
# scripts/autocalibrate.R
# Calibration AUTOMATIQUE : ajuste les parametres PD cles de CALIB pour
# minimiser l'ecart moindres-carres entre le modele et les profils observes
# de la REFERENCE (au lieu de tuner a la main).
#
# Parametres optimises : IC50ADC_scale, Emax_CMP, k_depl_direct, k_rep.
# (ED50_kill, Emax_MEP restent aux valeurs de CALIB.)
#
# Sorties : results/calibration.rds (lue par predict.R) + valeurs a copier
#           dans config/compounds.R + results/calibration_fit.pdf.
#
#     Rscript scripts/autocalibrate.R
############################################################
source("config/compounds.R"); source("model/engine.R"); load_model("model")
ref <- REFERENCE
dir.create("results", showWarnings = FALSE)

if (!file.exists(ref$profiles)) stop("Profils absents : ", ref$profiles)
obs <- read.csv(ref$profiles, stringsAsFactors = FALSE)

# fold observe par animal (baseline = 1re valeur), lignees Neut + Mono
baseline <- function(aid, L) { s <- obs[obs$Animal_Id == aid, ]; s <- s[order(s$time_day), ]
  v <- s[[L]][!is.na(s[[L]])]; if (length(v)) v[1] else NA }
mk <- function(L, b0) {
  df <- data.frame(dose = obs$Dose_mgkg, day = obs$time_day,
                   fo = mapply(function(a, x) x / baseline(a, L), obs$Animal_Id, obs[[L]]),
                   lign = L, b0 = b0, stringsAsFactors = FALSE)
  # signal de depletion, MAIS on exclut les points terminaux (animaux sacrifies :
  # zeros persistants tardifs) qui fausseraient la recuperation (k_rep -> 0).
  term <- df$day > 22 & df$fo < 0.05
  df[is.finite(df$fo) & df$day > 0 & df$fo <= 1 & !term, ]
}
OBS <- rbind(mk("Neut", init_pars$Neut0), mk("Mono", init_pars$Mono0))
n_term <- sum(obs$time_day > 22 & (obs$Neut / mapply(function(a,x) baseline(a,"Neut"), obs$Animal_Id, obs$Neut) < 0.05), na.rm = TRUE)
cat(sprintf("Points utilises pour le fit : %d (Neut+Mono, <= baseline, hors terminaux)\n", nrow(OBS)))

# -- sim allegee (tol relachee, pas 12h, hmax libre) : rapide pour le fit --
sim_fast <- function(p, dose) {
  st <- c(adc_hu_state0, init_state[!names(init_state) %in% c("C1", "C2", "Damage")])
  p$rate_fun <- make_adc_infusion(dose_mgkg = dose, BW_kg = ref$BW_kg,
                  Tinfu_h = ref$Tinfu_h, interval_h = 21 * 24, n_cycles = 1)
  o <- as.data.frame(lsoda(y = st, times = seq(0, 32 * 24, by = 12),
         func = pkpd_adc_fornari, parms = p, rtol = 1e-3, atol = 1e-5, maxsteps = 1e5))
  o$td <- o$time / 24; o
}

# -- SSE modele vs observe pour un jeu CALIB donne --
sse <- function(calib) {
  tot <- 0
  for (d in unique(OBS$dose)) {
    p <- tryCatch(build_pars(ref, ref$IC50_myelo, calib, ref$IC50_myelo), error = function(e) NULL)
    if (is.null(p)) return(1e6)
    o <- tryCatch(sim_fast(p, d), error = function(e) NULL)
    if (is.null(o) || any(!is.finite(o$Neut))) return(1e6)
    for (L in c("Neut", "Mono")) {
      b0 <- if (L == "Neut") init_pars$Neut0 else init_pars$Mono0
      sub <- OBS[OBS$dose == d & OBS$lign == L, ]
      if (!nrow(sub)) next
      pred <- approx(o$td, o[[L]] / b0, xout = sub$day, rule = 2)$y
      tot <- tot + sum((pred - sub$fo)^2)
    }
  }
  tot
}

# -- transformations BORNEES : chaque parametre reste dans une plage
#    physiquement plausible -> evite les solutions degenerees (identifiabilite).
BND <- list(IC50ADC_scale = c(1, 300),      # ug/mL (potency)
            Emax_CMP      = c(0.50, 0.999),  # blocage proliferation (reel pour un cytotoxique)
            k_depl_direct = c(0.02, 2.0),    # deplation directe
            k_rep         = c(5e-4, 0.05))   # reparation / recuperation
enc <- function(v, b) qlogis(min(max((v - b[1]) / (b[2] - b[1]), 1e-4), 1 - 1e-4))
dec <- function(th, b) b[1] + (b[2] - b[1]) * plogis(th)
to_calib <- function(th) { c <- CALIB
  c$IC50ADC_scale <- dec(th[1], BND$IC50ADC_scale); c$Emax_CMP <- dec(th[2], BND$Emax_CMP)
  c$k_depl_direct <- dec(th[3], BND$k_depl_direct); c$k_rep <- dec(th[4], BND$k_rep); c }
th0 <- c(enc(CALIB$IC50ADC_scale, BND$IC50ADC_scale), enc(CALIB$Emax_CMP, BND$Emax_CMP),
         enc(CALIB$k_depl_direct, BND$k_depl_direct), enc(CALIB$k_rep, BND$k_rep))

sse0 <- sse(CALIB)
cat(sprintf("SSE initiale (CALIB actuel) = %.4f\n", sse0))
cat("Optimisation en cours (Nelder-Mead)...\n")
fit <- optim(th0, function(th) sse(to_calib(th)), method = "Nelder-Mead",
             control = list(maxit = 120, reltol = 1e-5))
best <- to_calib(fit$par)
cat(sprintf("SSE finale = %.4f  (amelioration : -%.1f%%)\n",
            fit$value, 100 * (sse0 - fit$value) / sse0))

cat("\n--- CALIB : avant -> apres (a copier dans config/compounds.R) ---\n")
show <- function(n, fmt) cat(sprintf(paste0("  %-14s ", fmt, "  ->  ", fmt, "\n"), n, CALIB[[n]], best[[n]]))
show("IC50ADC_scale", "%8.3f"); show("Emax_CMP", "%8.4f")
show("k_depl_direct", "%8.4f"); show("k_rep", "%8.5f")

# -- alerte : parametre bute sur une borne -> fit mal contraint --
hit <- c()
for (n in names(BND)) { b <- BND[[n]]; v <- best[[n]]; tol <- 0.02 * (b[2] - b[1])
  if (v <= b[1] + tol) hit <- c(hit, sprintf("%s (borne BASSE %.4g)", n, b[1]))
  if (v >= b[2] - tol) hit <- c(hit, sprintf("%s (borne HAUTE %.4g)", n, b[2])) }
if (length(hit)) {
  cat("\n[!] Parametres qui BUTENT sur une borne -> fit mal contraint / degenere :\n")
  for (h in hit) cat("    -", h, "\n")
  cat("    Interpretation prudente : ces valeurs ne sont pas identifiables par tes\n",
      "    donnees. Envisage de FIXER l'ancrage de puissance ou d'ajouter des points.\n", sep = "")
} else cat("\n[ok] Aucun parametre ne bute sur une borne.\n")

saveRDS(best, "results/calibration.rds")
cat("\n-> results/calibration.rds  (lue par predict.R)\n")

# -- figure de controle (courbes optimisees vs points) --
pal <- grDevices::colorRampPalette(c("#2166ac", "#1a9641", "#f0a000", "#b2182b"))
cols <- setNames(pal(length(ref$doses)), as.character(ref$doses))
bN <- if (!is.null(ref$baseline_neut)) ref$baseline_neut else 2.4
bM <- if (!is.null(ref$baseline_mono)) ref$baseline_mono else 0.8
pdf("results/calibration_fit.pdf", width = 13, height = 5.4)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3, 1))
for (cc in list(c("Neut", init_pars$Neut0, "Neutrophiles"),
                c("Mono", init_pars$Mono0, "Monocytes"))) {
  L <- cc[1]; b <- as.numeric(cc[2])
  plot(NA, xlim = c(0, 32), ylim = c(0, 3.5), xlab = "Jour",
       ylab = "fold vs baseline", main = paste(cc[3], "—", ref$name, "(auto-calage)"))
  abline(h = 1, lty = 2, col = "grey60")
  draw_grade_lines(if (L == "Neut") bN else bM, ymax = 3.5, xmax = 32)
  for (d in ref$doses) {
    p <- build_pars(ref, ref$IC50_myelo, best, ref$IC50_myelo)
    o <- simulate(p, d, ref$BW_kg, ref$Tinfu_h, tmax_day = 35)
    lines(o$td, o[[L]] / b, col = cols[as.character(d)], lwd = 2.5)
    od <- obs[obs$Dose_mgkg == d, ]
    for (k in 1:nrow(od)) { bb <- baseline(od$Animal_Id[k], L); fo <- od[[L]][k] / bb
      if (!isTRUE(fo <= 1)) next
      points(od$time_day[k], fo, col = cols[as.character(d)], pch = 19, cex = 1.2) }
  }
  legend("topright", paste0(names(cols), " mg/kg"), col = cols, lwd = 2.5, pch = 19, bty = "n", cex = .8)
}
dev.off()
cat("-> results/calibration_fit.pdf\n")
