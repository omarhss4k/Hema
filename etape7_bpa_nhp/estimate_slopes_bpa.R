############################################################
# estimate_slopes_bpa.R
# Estimation Slope_CMP / Slope_MEP / D0 du BPA sur les comptages
# hemato NHP reels (fold-change vs baseline individuelle).
# PK 2-cmt FIXEE aux valeurs Shiny (a remplir ci-dessous).
#
# CONFIDENTIEL -- donnees composé BPA, usage interne.
############################################################
library(deSolve)
source("../etape2_carboplatin_humain/parameters_human.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("../etape5_tdxd_humain/parameters_tdxd_human.R")
source("../etape3_tdxd_rat/parameters_tdxd_rat.R")   # make_tdxd_infusion + pkpd_tdxd_fornari (2-cmt)

## ════════════════════════════════════════════════════════
## [1] PARAMS PK SHINY (singe) -- A REMPLIR
## ════════════════════════════════════════════════════════
BW_KG   <- 4.0        # poids moyen des singes (kg)  <-- verifier
CL_ADC  <- 0.000      # L/h        <-- REMPLIR (Shiny). Si Shiny en L/kg/h : x BW_KG
V1_ADC  <- 0.000      # L          <-- REMPLIR. Si Shiny en L/kg : x BW_KG
Q_ADC   <- 0.000      # L/h        <-- REMPLIR (2-cmt)
V2_ADC  <- 0.000      # L          <-- REMPLIR (2-cmt)
TINFU_H <- 1.0        # duree d'infusion (h)         <-- verifier (bolus court ?)

## ════════════════════════════════════════════════════════
## [2] DONNEES OBSERVEES
## ════════════════════════════════════════════════════════
obs <- read.csv("data_obs/bpa_hema.csv", stringsAsFactors = FALSE)
# lignees a ajuster (Neut via CMP ; Ret/RBC via MEP)
lineages <- c("Neut","RBC","Ret")     # Plt/Mono optionnels
# Exclure la phase de NEUTROPHILIE precoce (effet STING, hors modele Fornari) :
FIT_FROM_DAY <- 8      # ne fitte que jour >= 8 (ajuster/mettre 0 pour tout inclure)

## ════════════════════════════════════════════════════════
## [3] MODELE (PK 2-cmt fixe + Fornari, kill lineaire Slope)
## ════════════════════════════════════════════════════════
build_pars <- function(SlopeCMP, SlopeMEP, D0) {
  p <- c(init_pars, tdxd_pars_hu)
  p$CL_ADC<-CL_ADC; p$V1_ADC<-V1_ADC; p$Q_ADC<-Q_ADC; p$V2_ADC<-V2_ADC
  p$k_int <- 0                      # pas de cible HER2
  p$krel_power<-0; p$krel_factor<-1 # pas de linker T-DXd
  p$use_ADC_driver <- TRUE          # driver = C_ADC1 ; IC50_ADC herite (echelle Damage)
  p$Damage_threshold <- D0
  p$Slope_CMP<-SlopeCMP; p$Slope_MEP<-SlopeMEP; p$Slope_MPP<-SlopeCMP*1.39
  p[["k_dam"]]<-NULL; p[["k_rep"]]<-NULL
  p
}
state0 <- c(tdxd_hu_state0, init_state[!names(init_state) %in% c("C1","C2","Damage")])
base_mod <- c(Neut=init_pars$Neut0, RBC=init_pars$RBC0, Ret=init_pars$Ret0,
              Plt=init_pars$Plt0, Mono=init_pars$Mono0)

obs_days <- sort(unique(obs$time_day))
times <- seq(0, max(obs_days)+2, by=0.25)*24   # en heures
sim_fold <- function(p, dose) {
  p$rate_fun <- make_tdxd_infusion(dose_mgkg=dose, BW_kg=BW_KG, Tinfu_h=TINFU_H,
                  interval_h=21*24, n_cycles=1)
  sol <- tryCatch(as.data.frame(lsoda(y=state0, times=times, func=pkpd_tdxd_fornari,
           parms=p, rtol=1e-5, atol=1e-7, maxsteps=5e5, hmax=TINFU_H/2)), error=function(e) NULL)
  if (is.null(sol)) return(NULL)
  sol$td <- sol$time/24; sol
}

# baseline individuelle = 1re mesure de l'animal (j~1, proxy pre-dose)
base_ind <- function(aid, L) {
  s <- obs[obs$Animal_Id==aid, ]; s <- s[order(s$time_day), ]
  v <- s[[L]][!is.na(s[[L]])]; if(!length(v)) return(NA_real_); v[1]
}

## ════════════════════════════════════════════════════════
## [4] OBJECTIF (SSR fold-change) + ESTIMATION
## ════════════════════════════════════════════════════════
objective <- function(theta) {
  SlopeCMP<-exp(theta[1]); SlopeMEP<-exp(theta[2]); D0<-exp(theta[3])
  p <- build_pars(SlopeCMP, SlopeMEP, D0)
  ss<-0
  for (d in unique(obs$Dose_mgkg)) {
    sf <- sim_fold(p, d); if(is.null(sf)) return(1e6)
    od <- obs[obs$Dose_mgkg==d & obs$time_day>=FIT_FROM_DAY, ]
    for (k in 1:nrow(od)) {
      j<-od$time_day[k]; aid<-od$Animal_Id[k]; ip<-which.min(abs(sf$td-j))
      for (L in lineages) {
        bi<-base_ind(aid,L); ov<-od[[L]][k]; if(is.na(bi)||is.na(ov)||bi==0) next
        pf <- sf[[L]][ip]/base_mod[L]
        ss <- ss + (ov/bi - pf)^2
      }
    }
  }
  ss
}

cat("Estimation (Slope_CMP, Slope_MEP, D0) ...\n")
fit <- optim(par=log(c(1, 1, 0.05)), fn=objective, method="Nelder-Mead",
             control=list(maxit=500, reltol=1e-7))
S_CMP<-exp(fit$par[1]); S_MEP<-exp(fit$par[2]); D0<-exp(fit$par[3])
cat(sprintf("\n=== PARAMETRES ESTIMES (BPA, donnees NHP reelles) ===\n"))
cat(sprintf("  Slope_CMP = %.3f\n", S_CMP))
cat(sprintf("  Slope_MEP = %.3f\n", S_MEP))
cat(sprintf("  D0 (Damage_threshold) = %.4f\n", D0))
cat(sprintf("  SSR finale = %.4f  (convergence=%d)\n", fit$value, fit$convergence))

## ════════════════════════════════════════════════════════
## [5] QUALITE D'AJUSTEMENT + FIGURE obs vs modele
## ════════════════════════════════════════════════════════
p <- build_pars(S_CMP, S_MEP, D0)
cat("\nRMSE fold-change par lignee :\n")
for (L in lineages) {
  num<-0;n<-0
  for (d in unique(obs$Dose_mgkg)) {
    sf<-sim_fold(p,d); od<-obs[obs$Dose_mgkg==d & obs$time_day>=FIT_FROM_DAY,]
    for (k in 1:nrow(od)) {
      bi<-base_ind(od$Animal_Id[k],L); ov<-od[[L]][k]; if(is.na(bi)||is.na(ov)||bi==0) next
      pf<-sf[[L]][which.min(abs(sf$td-od$time_day[k]))]/base_mod[L]
      num<-num+(ov/bi-pf)^2; n<-n+1
    }
  }
  cat(sprintf("  %-5s : %.1f%%  (n=%d)\n", L, 100*sqrt(num/max(n,1)), n))
}

dir.create("results", showWarnings=FALSE)
pdf("results/BPA_fit_obs_vs_model.pdf", width=13, height=4.2)
par(mfrow=c(1,length(lineages)), mar=c(4,4.2,3,1))
cols<-c("0.3"="#2166ac","1"="#1a9641","3"="#b2182b")
for (L in lineages) {
  plot(NA, xlim=c(0,max(obs_days)), ylim=c(0, 2.5),
       xlab="Jour", ylab="Fold-change vs baseline", main=paste(L,"(pts=obs, lignes=modele)"))
  abline(h=1,lty=2,col="grey60"); abline(v=FIT_FROM_DAY,lty=3,col="grey70")
  for (d in unique(obs$Dose_mgkg)) {
    sf<-sim_fold(p,d); lines(sf$td, sf[[L]]/base_mod[L], col=cols[as.character(d)], lwd=2)
    od<-obs[obs$Dose_mgkg==d,]
    for (k in 1:nrow(od)){ bi<-base_ind(od$Animal_Id[k],L); ov<-od[[L]][k]
      if(!is.na(bi)&&!is.na(ov)&&bi>0) points(od$time_day[k], ov/bi, col=cols[as.character(d)], pch=19) }
  }
  legend("topright", legend=paste0(names(cols)," mg/kg"), col=cols, pch=19, lwd=2, bty="n", cex=0.8)
}
dev.off()
cat("\n-> results/BPA_fit_obs_vs_model.pdf\n")
