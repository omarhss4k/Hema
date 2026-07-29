############################################################
# run_bpa_v2_profiles.R
# Profils temporels hemato (8 lignees, P10/P50/P90) -- BPA NHP
# version v2 (architecture etape5 : mixture + kill lineaire, ratio IC50).
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("../etape5_tdxd_humain/parameters_tdxd_human.R")
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")
  source("pkpd_bpa_1cmt.R")
  source("../shared/bpa_build_params.R")
}))

DOSES<-c(0.1,0.3,0.5); SC_LAB<-"S2 (IC50_CMP=0.015 nM, ratio-Slope)"
BW_KG<-4; TINFU_H<-1.5; INTERVAL_H<-21*24; N_CYCLES<-4
FREQ<-sprintf("Q%dW",round(INTERVAL_H/(7*24)))
omega_CL<-0.35; omega_V1<-0.20; omega_Slope_CMP<-0.33; omega_Slope_MEP<-0.33

IC50_MEP_tdxd<-27.3/148000*1e6; IC50_CMP_tdxd<-28.1/148000*1e6
ratio_MEP<-IC50_MEP_tdxd/155; ratio_CMP<-IC50_CMP_tdxd/0.015

p<-build_bpa_pars(BW_target_kg=BW_KG,IC50_CMP_nM=7.8,ED50_kill=NULL,verbose=FALSE)
pars_pd<-init_pars; pars_pd[["k_dam"]]<-NULL; pars_pd[["k_rep"]]<-NULL
pt0<-modifyList(pars_pd,p)
state0<-c(c(C_ADC1=0,C_DXd=0,C_DXd_ic=0,Damage=0),init_state[!names(init_state)%in%c("C1","C2","Damage")])
times<-seq(0,N_CYCLES*INTERVAL_H+7*24,by=6); n_times<-length(times); times_j<-times/24

dir.create("results",showWarnings=FALSE)
pdf("results/BPA_v2_profiles_doses_S2.pdf",width=16,height=9)
for (DOSE in DOSES) {
set.seed(42); N<-100
mk<-function() matrix(NA_real_,n_times,N)
mMPP<-mk();mCMP<-mk();mMEP<-mk();mNeut<-mk();mMono<-mk();mRet<-mk();mRBC<-mk();mPlt<-mk()
for (i in 1:N) {
  pi<-pt0
  pi$CL_ADC<-pt0$CL_ADC*exp(rnorm(1,0,omega_CL)); pi$V1_ADC<-pt0$V1_ADC*exp(rnorm(1,0,omega_V1))
  rc<-runif(1); pr<-1-p_sensitive_tdxd-p_moderate_tdxd
  Sb<-if(rc<pr)Slope_resist_tdxd else if(rc<pr+p_moderate_tdxd)Slope_moderate_tdxd else Slope_sensitive_tdxd
  pi$Slope_CMP<-max(0,Sb*ratio_CMP*exp(rnorm(1,0,omega_Slope_CMP)))
  rm<-runif(1); pmr<-1-p_MEP_light_tdxd-p_sensitive_mep_tdxd
  Smb<-if(rm<pmr)Slope_MEP_resist_tdxd else if(rm<pmr+p_MEP_light_tdxd)Slope_MEP_light_tdxd else Slope_MEP_sensitive_tdxd
  pi$Slope_MEP<-max(0,Smb*ratio_MEP*exp(rnorm(1,0,omega_Slope_MEP)))
  pi$rate_fun<-make_tdxd_infusion(dose_mgkg=DOSE,BW_kg=BW_KG,Tinfu_h=TINFU_H,interval_h=INTERVAL_H,n_cycles=N_CYCLES)
  out<-tryCatch(as.data.frame(lsoda(y=state0,times=times,func=pkpd_bpa_fornari_1cmt,parms=pi,rtol=1e-5,atol=1e-7,maxsteps=5e5,hmax=TINFU_H/2)),error=function(e)NULL)
  if(!is.null(out)&&nrow(out)>=n_times){n<-n_times
    mMPP[,i]<-out$MPP[1:n];mCMP[,i]<-out$CMP[1:n];mMEP[,i]<-out$MEP[1:n];mNeut[,i]<-out$Neut[1:n]
    mMono[,i]<-out$Mono[1:n];mRet[,i]<-out$Ret[1:n];mRBC[,i]<-out$RBC[1:n];mPlt[,i]<-out$Plt[1:n]}
  if(i%%25==0) cat(sprintf("  %d/%d\n",i,N))
}
.pp<-function(mat,baseline,col,main,ylab,hl=NULL,hc=NULL){
  p10<-apply(mat,1,quantile,.10,na.rm=TRUE);p50<-apply(mat,1,quantile,.50,na.rm=TRUE);p90<-apply(mat,1,quantile,.90,na.rm=TRUE)
  ymin<-max(1e-3,min(p10[p10>0],na.rm=TRUE)*.5);ymax<-max(p90,na.rm=TRUE)*2
  plot(times_j,p50,type="n",log="y",ylim=c(ymin,ymax),xlim=range(times_j),xlab="Time (d)",ylab=ylab,main=main,
    panel.first={abline(v=seq(0,max(times_j),25),col="grey90");abline(h=c(1e-3,1e-2,.1,1,10,1e2,1e3,1e4),col="grey90")})
  polygon(c(times_j,rev(times_j)),c(pmax(p10,ymin*.9),rev(pmax(p90,ymin*.9))),col=adjustcolor(col,.25),border=NA)
  lines(times_j,pmax(p50,ymin*.9),col=col,lwd=2)
  if(!is.null(baseline)&&baseline>0)abline(h=baseline,lty=2,col="grey40",lwd=1.2)
  if(!is.null(hl))abline(h=hl,lty=2,col=hc,lwd=1.2)
  abline(v=(0:(N_CYCLES-1))*(INTERVAL_H/24),lty=3,col="grey60")
}
par(mfrow=c(2,4),mar=c(3.5,3.8,2.8,.5),oma=c(0,0,2.5,0))
.pp(mMPP,pt0$MPP0,"#4575b4","Multi-potent progenitors","10^9/L")
.pp(mNeut,pt0$Neut0,"#d73027","Neutrophils","10^9/L",c(2,1,.5),c("grey50","red","darkred"))
.pp(mCMP,pt0$CMP0,"#1a9641","Common myeloid progenitors","10^9/L")
.pp(mMono,pt0$Mono0,"#8c510a","Monocytes","10^9/L")
.pp(mMEP,pt0$MEP0,"#762a83","MEP","10^9/L")
.pp(mPlt,pt0$Plt0,"#35978f","Platelets","10^9/L",c(150,75,50,25),c("grey50","orange","red","darkred"))
.pp(mRet,pt0$Ret0,"#f1b6da","Reticulocytes","10^9/L")
.pp(mRBC,pt0$RBC0,"#c51b7d","Red blood cells","10^9/L",c(pt0$RBC0*.67,pt0$RBC0*.8),c("darkred","orange"))
mtext(sprintf("BPA NHP v2 (ratio-Slope) -- %.1f mg/kg %s x%d  |  %s  (N=%d)",DOSE,FREQ,N_CYCLES,SC_LAB,N),outer=TRUE,cex=1.1,font=2)
}
dev.off()
cat("-> results/BPA_v2_profiles_doses_S2.pdf\n")
