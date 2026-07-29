############################################################
# run_ADCcyto_G4.R
# ADC CYTOTOXIQUE (donnees NHP completes : PK + IC50 + profils)
# Modele Fornari IC50-driven + deplation directe des progeniteurs.
#
# Reproduit la myelosuppression G4 sur les 2 lignees myeloides
# (neutrophiles + monocytes), avec le bon timing (nadir tardif
# j16-22, via reparation lente k_rep), et epargne erythroide.
#
# Choix de modelisation :
#  - le blocage de proliferation seul plafonne a G3 ; on ajoute
#    une deplation directe (k_depl_direct, PK-couplee via Damage)
#    pour atteindre le G4 MONOCYTAIRE observe (points reels,
#    NON terminaux : chute des j8-15 avec animal vivant).
#  - la neutrophilie de STRESS precoce (j6-15) n'est PAS modelisee
#    (un cytotoxique ne fait que deprimer) -> attendu.
#
# Lancer depuis etape7_bpa_nhp/ :  Rscript run_ADCcyto_G4.R
#
# CONFIDENTIEL -- donnees internes (data_obs/ gitignore).
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("../etape5_tdxd_humain/parameters_tdxd_human.R")
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")   # make_tdxd_infusion
  source("../etape3_tdxd_rat/pkpd_tdxd_rat.R")          # pkpd_tdxd_fornari
}))

## ═══════ [1] PK ADC cytotoxique (singe, 2-cmt) ═══════
BW_KG   <- 2.5
CL_ADC  <- 0.00484307629862215   # L/h
V1_ADC  <- 0.04907800106132178   # L
Q_ADC   <- 0.01627593959372783   # L/h
V2_ADC  <- 0.1170226727748038    # L
TINFU_H <- 1.5
DOSES   <- c(0.3, 1, 3)          # mg/kg

## ═══════ [2] PARAMETRES PD CALIBRES ═══════
K_REP        <- 0.002    # reparation lente -> nadir tardif (j16-22)
ED50_KILL    <- 0.1
EMAX_CMP     <- 0.982    # blocage proliferation myeloide (quasi-max)
K_DEPL_DIR   <- 0.25     # deplation directe progeniteurs -> atteint le G4
EMAX_MEP     <- 1e-4     # erythroide EPARGNE (IC50 erythro >> myelo)
# Ancres de POTENCY de la reference (pour la translation vers d'autres ADC) :
IC50_REF     <- 0.008615 # nM     IC50 myelo in vitro de CET ADC cytotoxique
IC50ADC_REF  <- 27.70    # ug/mL  echelle IC50_ADC du modele a cette calibration

mkp <- function() {
  p <- c(init_pars, tdxd_pars_hu)
  p$CL_ADC <- CL_ADC; p$V1_ADC <- V1_ADC; p$Q_ADC <- Q_ADC; p$V2_ADC <- V2_ADC
  p$k_int <- 0; p$krel_power <- 0; p$krel_factor <- 1     # pas de cible / linker T-DXd
  p$use_ADC_driver <- TRUE; p$Damage_threshold <- 0
  p$k_rep <- K_REP; p$ED50_kill <- ED50_KILL
  p$Emax_CMP_kill <- EMAX_CMP; p$Emax_MPP_kill <- 1; p$Emax_MEP_kill <- EMAX_MEP
  p$Slope_CMP <- 1; p$Slope_MPP <- 1; p$Slope_MEP <- 1e-6
  p$k_depl_direct <- K_DEPL_DIR; p$k_depl_direct_MPP <- K_DEPL_DIR * 0.3
  p
}

## ═══════ [3] SIMULATION ═══════
sim <- function(dose) {
  p  <- mkp()
  st <- c(tdxd_hu_state0, init_state[!names(init_state) %in% c("C1", "C2", "Damage")])
  p$rate_fun <- make_tdxd_infusion(dose_mgkg = dose, BW_kg = BW_KG,
                  Tinfu_h = TINFU_H, interval_h = 21 * 24, n_cycles = 1)
  o <- as.data.frame(lsoda(y = st, times = seq(0, 35 * 24, by = 3),
         func = pkpd_tdxd_fornari, parms = p,
         rtol = 1e-5, atol = 1e-7, maxsteps = 5e5, hmax = 0.75))
  o$td <- o$time / 24
  o
}
sims <- setNames(lapply(DOSES, sim), as.character(DOSES))

## ═══════ [4] DONNEES OBSERVEES ═══════
obs <- read.csv("data_obs/bpa_hema.csv", stringsAsFactors = FALSE)
obs$jr <- round(obs$time_day)
bi <- function(aid, L) {
  s <- obs[obs$Animal_Id == aid, ]; s <- s[order(s$time_day), ]
  v <- s[[L]][!is.na(s[[L]])]; if (!length(v)) NA else v[1]
}

## ═══════ [5] FIGURE (neutrophiles + monocytes) ═══════
cols <- c("0.3" = "#2166ac", "1" = "#1a9641", "3" = "#b2182b")
dir.create("results", showWarnings = FALSE)
pdf("results/ADCcyto_G4_neut_mono.pdf", width = 13, height = 5.4)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3, 1))
for (cc in list(c("Neut", init_pars$Neut0, "Neutrophiles"),
                c("Mono", init_pars$Mono0, "Monocytes"))) {
  L <- cc[1]; b <- as.numeric(cc[2])
  plot(NA, xlim = c(0, 32), ylim = c(0, 3.5), xlab = "Jour",
       ylab = "fold vs baseline", main = paste(cc[3], "— ADC cytotoxique"))
  abline(h = 1, lty = 2, col = "grey60")
  abline(h = 0.21, lty = 3, col = "red"); text(31, 0.28, "G4", col = "red", cex = .7, adj = 1)
  for (d in names(sims)) {
    o <- sims[[d]]; lines(o$td, o[[L]] / b, col = cols[d], lwd = 2.5)
    od <- obs[obs$Dose_mgkg == as.numeric(d), ]
    for (k in 1:nrow(od)) {
      bb <- bi(od$Animal_Id[k], L); fo <- od[[L]][k] / bb
      term <- (od$jr[k] > 22 & fo < 0.05)   # cercle = terminal/sacrifice
      points(od$jr[k], fo, col = cols[d], pch = ifelse(term, 1, 19), cex = 1.2, lwd = 1.5)
    }
  }
  legend("topright", paste0(names(cols), " mg/kg"), col = cols, lwd = 2.5, pch = 19, bty = "n", cex = .8)
}
dev.off()

## ═══════ [6] RESUME CONSOLE (grades CTCAE) ═══════
cn <- function(x) if (x < 0.5) "G4" else if (x < 1) "G3" else if (x < 1.5) "G2" else if (x < 2) "G1" else "G0"
BASE_NEUT <- 2.4; BASE_MONO <- 0.8   # baselines moyennes singe (x10^9/L) pour le grade absolu
cat("\n=== ADC cytotoxique | myelosuppression (Fornari + deplation directe) ===\n")
for (d in names(sims)) {
  o <- sims[[d]]
  nN <- min(o$Neut) / init_pars$Neut0; nM <- min(o$Mono) / init_pars$Mono0
  cat(sprintf("dose %-4s mg/kg : Neut nadir=%.2f@j%.0f (%s) | Mono nadir=%.2f@j%.0f (%s)\n",
    d, nN, o$td[which.min(o$Neut)], cn(nN * BASE_NEUT),
    nM, o$td[which.min(o$Mono)], cn(nM * BASE_MONO)))
}
cat("-> results/ADCcyto_G4_neut_mono.pdf\n")

## ═══════ [7] SAUVEGARDE DE LA CALIBRATION (alimente run_BPA_prediction.R) ═══════
calib <- list(K_REP = K_REP, ED50_KILL = ED50_KILL, EMAX_CMP = EMAX_CMP,
              K_DEPL_DIR = K_DEPL_DIR, EMAX_MEP = EMAX_MEP,
              IC50_REF = IC50_REF, IC50ADC_REF = IC50ADC_REF)
saveRDS(calib, "results/calib_cyto.rds")
cat("-> results/calib_cyto.rds  (calibration -> lue par run_BPA_prediction.R)\n")
