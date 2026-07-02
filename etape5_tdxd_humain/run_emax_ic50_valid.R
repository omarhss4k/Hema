############################################################
# run_emax_ic50_valid.R
# Validation 1ere passe : mecanisme Emax-IC50 (deplation directe)
# vs grades FDA (etape5 T-DXd 5.4 mg/kg Q3W x6).
# Mixture kmax de 1ere passe (NON finement calibree).
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("parameters_tdxd_human.R")
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")
  source("pkpd_emax_ic50_test.R")
}))

pars_pd <- init_pars; pars_pd[["k_dam"]]<-NULL; pars_pd[["k_rep"]]<-NULL
pt0 <- c(pars_pd, tdxd_pars_hu)
pt0$IC50_CMP_nM<-189.9; pt0$IC50_MEP_nM<-184.5; pt0$IC50_MPP_nM<-189.9
pt0$kmax_MPP<-0; pt0$kmax_MEP<-0; pt0$kmax_plt_prol<-0
state0 <- c(tdxd_hu_state0, init_state[!names(init_state) %in% c("C1","C2","Damage")])
times <- seq(0, 126*24, by=6)

cn<-function(x)if(x<0.5)"G4"else if(x<1)"G3"else if(x<1.5)"G2"else if(x<2)"G1"else"G0"
ca<-function(x,x0){f<-x/x0;if(f<0.54)"G4"else if(f<0.67)"G3"else if(f<0.8)"G2"else if(f<0.9)"G1"else"G0"}
go<-c("G0","G1","G2","G3","G4")

set.seed(42); N<-100
gN<-gA<-character(N)
for (i in 1:N) {
  pt <- pt0
  pt$CL_ADC <- pt0$CL_ADC*exp(rnorm(1,0,0.35)); pt$V1_ADC <- pt0$V1_ADC*exp(rnorm(1,0,0.20))
  # Mixture kmax_CMP (neutropenie) : proportions FDA (66/14/20)
  rc<-runif(1)
  pt$kmax_CMP <- if(rc<0.66) 0.02*exp(rnorm(1,0,0.3))
                 else if(rc<0.80) 0.15*exp(rnorm(1,0,0.3))
                 else 0.48*exp(rnorm(1,0,0.3))
  # Mixture kmax_ret_prol (anemie) : proportions FDA (34/33/33)
  rm<-runif(1)
  pt$kmax_ret_prol <- if(rm<0.34) 0.005*exp(rnorm(1,0,0.3))
                      else if(rm<0.67) 0.04*exp(rnorm(1,0,0.3))
                      else 0.18*exp(rnorm(1,0,0.5))
  pt$rate_fun <- make_tdxd_infusion(dose_mgkg=5.4, BW_kg=70, Tinfu_h=1.5, interval_h=21*24, n_cycles=6)
  out<-tryCatch(as.data.frame(lsoda(y=state0,times=times,func=pkpd_emax_ic50,parms=pt,rtol=1e-5,atol=1e-7,maxsteps=5e5,hmax=0.75)),error=function(e)NULL)
  if(!is.null(out)){gN[i]<-cn(min(out$Neut,na.rm=T));gA[i]<-ca(min(out$RBC,na.rm=T),pt0$RBC0)}else gN[i]<-NA
  if(i%%20==0) cat(sprintf("  %d/%d\n",i,N))
}
ok<-!is.na(gN)
tn<-round(100*table(factor(gN[ok],levels=go))/sum(ok))
tan<-round(100*table(factor(gA[ok],levels=go))/sum(ok))
cat("\n=== VALIDATION Emax-IC50 vs FDA (T-DXd 5.4 mg/kg, N=",sum(ok),") ===\n")
cat("NEUTROPENIE :", paste(go,tn,sep="=",collapse="  "), sprintf("| G3-4=%d%% (FDA 16%%)\n", tn["G3"]+tn["G4"]))
cat("ANEMIE      :", paste(go,tan,sep="=",collapse="  "), sprintf("| G3-4=%d%% (FDA 9%%)\n", tan["G3"]+tan["G4"]))
