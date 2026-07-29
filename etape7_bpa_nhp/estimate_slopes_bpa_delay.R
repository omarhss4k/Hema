############################################################
# estimate_slopes_bpa_delay.R
# Estimation Slope_CMP / Slope_MEP / k_delay du BPA -- modele a
# SIGNAL TOXIQUE RETARDE (capture la toxicite decouplee de la PK).
# PK 2-cmt FIXEE aux valeurs Shiny (a remplir).
#
# CONFIDENTIEL -- donnees BPA, usage interne.
############################################################
library(deSolve)
source("../etape2_carboplatin_humain/parameters_human.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("../etape5_tdxd_humain/parameters_tdxd_human.R")
source("../etape3_tdxd_rat/parameters_tdxd_rat.R")   # make_tdxd_infusion
source("pkpd_bpa_delay.R")

## ═══════ [1] PARAMS PK SHINY (singe) -- A REMPLIR ═══════
BW_KG   <- 4.0
CL_ADC  <- 0.000      # L/h   <-- REMPLIR (x BW_KG si /kg)
V1_ADC  <- 0.000      # L     <-- REMPLIR
Q_ADC   <- 0.000      # L/h   <-- REMPLIR
V2_ADC  <- 0.000      # L     <-- REMPLIR
TINFU_H <- 1.0

## ═══════ [2] DONNEES + REGLAGES ═══════
obs <- read.csv("data_obs/bpa_hema.csv", stringsAsFactors = FALSE)
lineages <- c("RBC","Ret")     # ajouter "Neut" ensuite
FIT_FROM_DAY <- 10             # exclut la hausse STING (pic j8) ; le delai gere le reste
D0_FIX <- 0.0                  # seuil sur le SIGNAL (0 = pas de seuil, le delai suffit)

## ═══════ [3] MODELE ═══════
build_pars <- function(SlopeCMP, SlopeMEP, k_delay) {
  p <- c(init_pars, tdxd_pars_hu)
  p$CL_ADC<-CL_ADC; p$V1_ADC<-V1_ADC; p$Q_ADC<-Q_ADC; p$V2_ADC<-V2_ADC
  p$k_int<-0; p$krel_power<-0; p$krel_factor<-1
  p$use_ADC_driver<-TRUE
  p$Damage_threshold<-D0_FIX
  p$k_delay <- k_delay
  p$Slope_CMP<-SlopeCMP; p$Slope_MEP<-SlopeMEP; p$Slope_MPP<-SlopeCMP*1.39
  p[["k_dam"]]<-NULL; p[["k_rep"]]<-NULL
  p
}
# etat initial : + Sig1, Sig2 (=0) apres Damage
pd0 <- init_state[!names(init_state) %in% c("C1","C2","Damage")]
state0 <- c(C_ADC1=0,C_ADC2=0,C_DXd=0,C_DXd_ic=0,Damage=0, Sig1=0,Sig2=0, pd0)
base_mod <- c(Neut=init_pars$Neut0, RBC=init_pars$RBC0, Ret=init_pars$Ret0,
              Plt=init_pars$Plt0, Mono=init_pars$Mono0)

obs_days <- sort(unique(obs$time_day))
times <- seq(0, max(obs_days)+2, by=0.25)*24
sim_fold <- function(p, dose) {
  p$rate_fun <- make_tdxd_infusion(dose_mgkg=dose, BW_kg=BW_KG, Tinfu_h=TINFU_H, interval_h=21*24, n_cycles=1)
  sol <- tryCatch(as.data.frame(lsoda(y=state0, times=times, func=pkpd_bpa_delay, parms=p,
           rtol=1e-5, atol=1e-7, maxsteps=5e5, hmax=TINFU_H/2)), error=function(e) NULL)
  if (is.null(sol)) return(NULL); sol$td<-sol$time/24; sol
}
base_ind <- function(aid,L){s<-obs[obs$Animal_Id==aid,];s<-s[order(s$time_day),];v<-s[[L]][!is.na(s[[L]])];if(!length(v))return(NA_real_);v[1]}

## ═══════ [4] OBJECTIF + ESTIMATION (Slope_CMP, Slope_MEP, k_delay) ═══════
objective <- function(theta) {
  SlopeCMP<-exp(theta[1]); SlopeMEP<-exp(theta[2]); kdel<-exp(theta[3])
  p <- build_pars(SlopeCMP, SlopeMEP, kdel); ss<-0
  for (d in unique(obs$Dose_mgkg)) {
    sf<-sim_fold(p,d); if(is.null(sf)) return(1e6)
    od<-obs[obs$Dose_mgkg==d & obs$time_day>=FIT_FROM_DAY,]
    for (k in 1:nrow(od)) {
      j<-od$time_day[k]; aid<-od$Animal_Id[k]; ip<-which.min(abs(sf$td-j))
      for (L in lineages) {
        bi<-base_ind(aid,L); ov<-od[[L]][k]; if(is.na(bi)||is.na(ov)||bi==0) next
        ss<-ss+(ov/bi - sf[[L]][ip]/base_mod[L])^2
      }
    }
  }
  ss
}
cat(sprintf("Estimation (Slope_CMP, Slope_MEP, k_delay) | D0=%.3f | fit j>=%d ...\n", D0_FIX, FIT_FROM_DAY))
# init : k_delay ~ 0.008 /h -> delai ~ 2/0.008/24 ~ 10 j
fit <- optim(par=log(c(1,1,0.008)), fn=objective, method="Nelder-Mead", control=list(maxit=600, reltol=1e-7))
S_CMP<-exp(fit$par[1]); S_MEP<-exp(fit$par[2]); KDEL<-exp(fit$par[3])
delai_j <- 2/KDEL/24
cat(sprintf("\n=== PARAMETRES ESTIMES (BPA, modele a delai) ===\n"))
cat(sprintf("  Slope_CMP = %.3f\n", S_CMP))
cat(sprintf("  Slope_MEP = %.3f\n", S_MEP))
cat(sprintf("  k_delay   = %.5f /h  -> delai signal ~ %.1f jours\n", KDEL, delai_j))
cat(sprintf("  SSR finale = %.4f  (conv=%d)\n", fit$value, fit$convergence))

## ═══════ [5] RMSE + FIGURE ═══════
p<-build_pars(S_CMP,S_MEP,KDEL)
cat("\nRMSE fold-change par lignee :\n")
for (L in lineages) { num<-0;n<-0
  for (d in unique(obs$Dose_mgkg)) { sf<-sim_fold(p,d); od<-obs[obs$Dose_mgkg==d & obs$time_day>=FIT_FROM_DAY,]
    for (k in 1:nrow(od)){bi<-base_ind(od$Animal_Id[k],L);ov<-od[[L]][k];if(is.na(bi)||is.na(ov)||bi==0)next
      num<-num+(ov/bi-sf[[L]][which.min(abs(sf$td-od$time_day[k]))]/base_mod[L])^2;n<-n+1}}
  cat(sprintf("  %-5s : %.1f%%  (n=%d)\n",L,100*sqrt(num/max(n,1)),n)) }
dir.create("results",showWarnings=FALSE)
pdf("results/BPA_fit_delay.pdf",width=13,height=4.2); par(mfrow=c(1,length(lineages)),mar=c(4,4.2,3,1))
cols<-c("0.3"="#2166ac","1"="#1a9641","3"="#b2182b")
for (L in lineages) {
  plot(NA,xlim=c(0,max(obs_days)),ylim=c(0,2.5),xlab="Jour",ylab="Fold-change",main=paste(L,"(delai)"))
  abline(h=1,lty=2,col="grey60"); abline(v=FIT_FROM_DAY,lty=3,col="grey70")
  for (d in unique(obs$Dose_mgkg)) { sf<-sim_fold(p,d); lines(sf$td,sf[[L]]/base_mod[L],col=cols[as.character(d)],lwd=2)
    od<-obs[obs$Dose_mgkg==d,]; for(k in 1:nrow(od)){bi<-base_ind(od$Animal_Id[k],L);ov<-od[[L]][k]
      if(!is.na(bi)&&!is.na(ov)&&bi>0)points(od$time_day[k],ov/bi,col=cols[as.character(d)],pch=19)} }
  legend("topright",legend=paste0(names(cols)," mg/kg"),col=cols,pch=19,lwd=2,bty="n",cex=0.8)
}
dev.off(); cat("\n-> results/BPA_fit_delay.pdf\n")
