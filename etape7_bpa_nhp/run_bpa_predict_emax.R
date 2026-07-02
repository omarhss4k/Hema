############################################################
# run_bpa_predict_emax.R
# PREDICTION HEMATOTOXICITE BPA -- approche Emax-IC50 (IC50 = EC50 direct)
#
# kmax CALIBRES sur T-DXd (etape5, reproduit grades FDA) -- CONSERVES ici
#   (amplitude = propriete payload/mecanisme cytotoxique).
# On change SEULEMENT les IC50 (potence BPA entre via l'EC50) :
#   IC50_CMP_nM = 7.8 (S1) / 0.015 (S2)   [myeloide]
#   IC50_MEP_nM = 155                      [erythroide, commun]
# + PK 1-cmt allometrique souris->singe (build_bpa_pars).
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("../etape5_tdxd_humain/parameters_tdxd_human.R")
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")
  source("pkpd_emax_ic50_1cmt.R")
  source("../shared/bpa_build_params.R")
}))

# -- PK BPA NHP (allometrie) : on ne garde que CL/V1/k_int etc. ----------------
BW_KG<-4; TINFU_H<-1.5; INTERVAL_H<-21*24; N_CYCLES<-4
FREQ <- sprintf("Q%dW", round(INTERVAL_H/(7*24)))
p_pk <- build_bpa_pars(BW_target_kg=BW_KG, IC50_CMP_nM=7.8, ED50_kill=NULL, verbose=FALSE)
pars_pd <- init_pars; pars_pd[["k_dam"]]<-NULL; pars_pd[["k_rep"]]<-NULL
pt0 <- modifyList(pars_pd, p_pk)
pt0$IC50_MEP_nM <- 155         # erythroide BPA (commun S1/S2)
pt0$IC50_MPP_nM <- 7.8         # myeloide amont (= CMP)
pt0$kmax_MPP<-0; pt0$kmax_MEP<-0; pt0$kmax_plt_prol<-0

state0 <- c(c(C_ADC1=0,C_DXd=0,C_DXd_ic=0,Damage=0),
            init_state[!names(init_state) %in% c("C1","C2","Damage")])
times  <- seq(0, N_CYCLES*INTERVAL_H + 7*24, by=6)

cn<-function(x)if(x<0.5)"G4"else if(x<1)"G3"else if(x<1.5)"G2"else if(x<2)"G1"else"G0"
ca<-function(x,x0){f<-x/x0;if(f<0.54)"G4"else if(f<0.67)"G3"else if(f<0.8)"G2"else if(f<0.9)"G1"else"G0"}
cp<-function(x)if(x<25)"G4"else if(x<50)"G3"else if(x<75)"G2"else if(x<150)"G1"else"G0"
go<-c("G0","G1","G2","G3","G4")

# kmax CALIBRES (etape5/FDA) : mixture par sous-groupe
run_dose <- function(dose, IC50_CMP, N=100, seed=42) {
  set.seed(seed); gN<-gA<-gP<-character(N)
  for (i in 1:N) {
    pt <- pt0
    pt$IC50_CMP_nM <- IC50_CMP
    pt$CL_ADC <- pt0$CL_ADC*exp(rnorm(1,0,0.35)); pt$V1_ADC <- pt0$V1_ADC*exp(rnorm(1,0,0.20))
    rc<-runif(1)   # mixture kmax_CMP (neutro) -- valeurs calibrees T-DXd
    pt$kmax_CMP <- if(rc<0.66) 0.02*exp(rnorm(1,0,0.3))
                   else if(rc<0.80) 0.13*exp(rnorm(1,0,0.3))
                   else 0.80*exp(rnorm(1,0,0.35))
    rm<-runif(1)   # mixture kmax_ret_prol (anemie)
    pt$kmax_ret_prol <- if(rm<0.34) 0.005*exp(rnorm(1,0,0.3))
                        else if(rm<0.67) 0.035*exp(rnorm(1,0,0.3))
                        else 0.125*exp(rnorm(1,0,0.5))
    pt$rate_fun <- make_tdxd_infusion(dose_mgkg=dose, BW_kg=BW_KG,
                    Tinfu_h=TINFU_H, interval_h=INTERVAL_H, n_cycles=N_CYCLES)
    out<-tryCatch(as.data.frame(lsoda(y=state0,times=times,func=pkpd_emax_ic50_1cmt,
          parms=pt,rtol=1e-5,atol=1e-7,maxsteps=5e5,hmax=TINFU_H/2)),error=function(e)NULL)
    if(!is.null(out)){gN[i]<-cn(min(out$Neut,na.rm=T));gA[i]<-ca(min(out$RBC,na.rm=T),pt0$RBC0)
                      gP[i]<-cp(min(out$Plt,na.rm=T))}else gN[i]<-NA
  }
  ok<-!is.na(gN)
  list(neut=round(100*table(factor(gN[ok],levels=go))/sum(ok)),
       anem=round(100*table(factor(gA[ok],levels=go))/sum(ok)),
       plt =round(100*table(factor(gP[ok],levels=go))/sum(ok)))
}

doses <- c(0.1,0.25,0.5,1.0,2.0,3.6,5.4)
print_sc <- function(lab, IC50_CMP) {
  cat(sprintf("\n\n=== BPA NHP Emax-IC50 -- %s ===\n", lab))
  cat(sprintf("kmax calibres FDA | IC50_CMP=%.4g nM (EC50) | %s x%d | N=100/dose\n", IC50_CMP, FREQ, N_CYCLES))
  cat(sprintf("%-8s | %-24s | %-24s | %-24s\n","Dose","NEUTRO (G0..G4)","ANEMIE (G0..G4)","THROMBO (G0..G4)"))
  cat(paste(rep("-",92),collapse=""),"\n")
  for (d in doses) {
    r<-run_dose(d,IC50_CMP); f<-function(t)paste(sprintf("%2.0f",t),collapse="/"); g<-function(t)t["G3"]+t["G4"]
    cat(sprintf("%-8s | %-14s G34=%3.0f%% | %-14s G34=%3.0f%% | %-14s G34=%3.0f%%\n",
        d,f(r$neut),g(r$neut),f(r$anem),g(r$anem),f(r$plt),g(r$plt)))
  }
}
print_sc("SCENARIO 1 (IC50_CMP=7.8 nM)",  7.8)
print_sc("SCENARIO 2 (IC50_CMP=0.015 nM)", 0.015)
cat("\n\nNOTE : kmax figes depuis la calibration FDA T-DXd ; seules les IC50 changent.\n")
cat("Potence BPA portee par l'EC50 (IC50 basse -> effet sature a basse dose).\n")
