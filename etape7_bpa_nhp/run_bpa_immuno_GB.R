############################################################
# run_bpa_immuno_GB.R
# Simulation du modele IMMUNO-PD (BPA, ADC-STING) sur les
# GLOBULES BLANCS (neutrophiles + monocytes) vs donnees NHP.
#
# Reproduit : rebond dose-dependant a ~j7 (granulopoiese
# d'urgence, pilotee par le mediateur aigu Sa) PUIS effondrement
# myeloide retarde (pilote par le mediateur chronique Sc).
#
# PK ADC 2-cmt fixee aux valeurs Shiny (singe, BW 2.5 kg).
#
# Lancer depuis etape7_bpa_nhp/ :   Rscript run_bpa_immuno_GB.R
#
# CONFIDENTIEL -- donnees BPA, usage interne (data_obs/ gitignore).
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("../etape5_tdxd_humain/parameters_tdxd_human.R")
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")   # make_tdxd_infusion
  source("pkpd_bpa_immuno.R")
}))

## ═══════ [1] PK SHINY (singe, 2-cmt) + POIDS ═══════
BW_KG   <- 2.5
CL_ADC  <- 0.00484307629862215   # L/h
V1_ADC  <- 0.04907800106132178   # L
Q_ADC   <- 0.01627593959372783   # L/h
V2_ADC  <- 0.1170226727748038    # L
TINFU_H <- 1.5                    # duree perfusion (h)
DOSES   <- c(0.3, 1, 3)           # mg/kg

## ═══════ [2] PARAMETRES IMMUNO-PD CALIBRES ═══════
# Sa = mediateur AIGU (pulse, pic ~j4-6)  -> pilote le REBOND
# Sc = mediateur CHRONIQUE (retarde, pic ~j11-15) -> pilote l'EFFONDREMENT
mkp <- function() {
  p <- init_pars
  p$CL_ADC <- CL_ADC; p$V1_ADC <- V1_ADC; p$Q_ADC <- Q_ADC; p$V2_ADC <- V2_ADC
  p$EC50_ADC <- 30;   p$Smax <- 1
  # -- mediateurs --
  p$k_a    <- 0.05;   p$koff_a <- 0.012    # Sa : pulse aigu (T1/2 ~2.4 j)
  p$k_c    <- 0.030;  p$koff_c <- 0.004    # Sc : chronique auto-amplifie, lent
  # -- rebond (granulopoiese d'urgence, regime lineaire) --
  p$Emax_stim <- 300; p$EC50_stim <- 3
  # -- deplation myeloide selective (erythroide epargnee) --
  p$k_depl_CMP <- 0.60; p$k_depl_MPP <- 0.12; p$k_depl_MEP <- 1e-4
  p
}

## ═══════ [3] SIMULATION ═══════
sim <- function(dose) {
  p <- mkp()
  pd0 <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
  s0  <- c(C_ADC1 = 0, C_ADC2 = 0, Sa = 0, Sc = 0, pd0)
  p$rate_fun <- make_tdxd_infusion(dose_mgkg = dose, BW_kg = BW_KG,
                  Tinfu_h = TINFU_H, interval_h = 21 * 24, n_cycles = 1)
  o <- as.data.frame(lsoda(y = s0, times = seq(0, 30 * 24, by = 3),
         func = pkpd_bpa_immuno, parms = p,
         rtol = 1e-6, atol = 1e-8, maxsteps = 1e6, hmax = 0.75))
  o$td <- o$time / 24
  o
}
p0   <- mkp()
sims <- setNames(lapply(DOSES, sim), as.character(DOSES))

## ═══════ [4] DONNEES OBSERVEES (moyenne des 2 singes/dose) ═══════
obs <- read.csv("data_obs/bpa_hema.csv", stringsAsFactors = FALSE)
obs$jr <- round(obs$time_day)
doseanim <- list("0.3" = c(3601, 3602), "1" = c(1501, 1502), "3" = c(2501, 2502))
meanfold <- function(dose, col) {
  aids <- doseanim[[dose]]
  fl <- lapply(aids, function(a) {
    s <- obs[obs$Animal_Id == a, ]; s <- s[order(s$time_day), ]
    b <- s[[col]][!is.na(s[[col]])][1]
    data.frame(jr = s$jr, f = s[[col]] / b)
  })
  d <- do.call(rbind, fl); ag <- aggregate(f ~ jr, d, mean, na.rm = TRUE)
  ag[order(ag$jr), ]
}

## ═══════ [5] FIGURE ═══════
cols <- c("0.3" = "#2166ac", "1" = "#1a9641", "3" = "#b2182b")
dir.create("results", showWarnings = FALSE)
pdf("results/BPA_immuno_GB.pdf", width = 13, height = 5.2)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
for (cc in list(c("Neut", p0$Neut0, "Neutrophiles"),
                c("Mono", p0$Mono0, "Monocytes"))) {
  col <- cc[1]; b <- as.numeric(cc[2]); lab <- cc[3]
  ym <- max(3.5, max(unlist(lapply(names(doseanim), function(d) meanfold(d, col)$f)), na.rm = TRUE) * 1.05)
  plot(NA, xlim = c(0, 30), ylim = c(0, ym), xlab = "Jour",
       ylab = "Fold-change (moy. 2 singes)", main = paste(lab, "— globules blancs"))
  abline(h = 1, lty = 2, col = "grey60")
  for (d in names(sims)) {
    o <- sims[[d]]; lines(o$td, o[[col]] / b, col = cols[d], lwd = 2.4)
    mf <- meanfold(d, col); points(mf$jr, mf$f, pch = 19, col = cols[d], cex = 1.4)
  }
  legend("topright", paste0(names(cols), " mg/kg"), col = cols, lwd = 2.4, pch = 19, bty = "n", cex = 0.85)
}
dev.off()

## ═══════ [6] RESUME CONSOLE ═══════
cat("\n=== BPA immuno-PD | globules blancs ===\n")
for (d in names(sims)) {
  o <- sims[[d]]
  cat(sprintf("dose %-4s mg/kg : Neut rebond x%.1f@j%.0f  nadir=%.2f@j%.0f  |  Mono nadir=%.2f@j%.0f\n",
    d, max(o$Neut) / p0$Neut0, o$td[which.max(o$Neut)],
    min(o$Neut) / p0$Neut0, o$td[which.min(o$Neut)],
    min(o$Mono) / p0$Mono0, o$td[which.min(o$Mono)]))
}
cat("-> results/BPA_immuno_GB.pdf\n")
