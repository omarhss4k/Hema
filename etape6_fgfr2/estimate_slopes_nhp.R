############################################################
# estimate_slopes_nhp.R  --  TEST de la technique d'estimation
# Estime Slope_CMP / Slope_MEP du composé interne par ajustement
# du modèle Fornari NHP aux COMPTAGES hemato observes (4 doses,
# 4 lignees). Meme principe que nlmixr2 (minimisation modele-donnees),
# via optim() sur le modele deSolve existant.
#
# NORMALISATION FOLD-CHANGE : chaque animal est normalise par sa
# baseline pre-dose (jour<0) ; les predictions par la baseline du
# modele. -> on compare des variations relatives, pas des valeurs
# absolues (sinon les ecarts de baseline inter-animaux dominent).
#
# CONFIDENTIEL : donnees composé interne -- usage interne uniquement.
############################################################
library(deSolve)
source("parameters_fgfr2_nhp.R")
source("parameters_nhp.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("pkpd_fgfr2_nhp.R")
source("data_nhp.R")     # -> obs_data

build_pars <- function() {
  p <- init_pars
  for (nm in names(fgfr2_nhp)) p[[nm]] <- fgfr2_nhp[[nm]]
  p$k_dam_DXd <- 0.017; p$rate_fun <- NULL; p
}
pd_state <- init_state[setdiff(names(init_state), c("C1","C2","Damage"))]
nhp_state0 <- c(C_ADC1=0, C_ADC2=0, Damage=0, pd_state)

BW<-4; TINFU<-0.5
doses <- sort(unique(obs_data$dose_mgkg))
obs_days <- sort(unique(obs_data$jour[obs_data$jour>=0]))
times <- seq(0, (max(obs_days)+3)*24, by=3)
lineages <- c("Neut","Plt","RBC","Ret")
base_mod <- c(Neut=init_pars$Neut0, Plt=init_pars$Plt0, RBC=init_pars$RBC0, Ret=init_pars$Ret0)

# -- baselines INDIVIDUELLES (moyenne des mesures pre-dose jour<0) ---------
pre <- obs_data[obs_data$jour < 0, ]
base_ind <- function(animal, L) {
  v <- pre[[L]][pre$Animal_Id==animal]; if(!length(v)||all(is.na(v))) return(NA_real_); mean(v,na.rm=TRUE)
}

# -- prediction fold-change par dose (pred/baseline_modele) ---------------
sim_fold <- function(p, dose) {
  p$rate_fun <- make_nhp_infusion(dose_mgkg=dose, BW_kg=BW, Tinfu_h=TINFU,
                  interval_h=fgfr2_nhp$interval_h, n_cycles=1)
  sol <- tryCatch(as.data.frame(ode(y=nhp_state0, times=times, func=pkpd_nhp_ode,
           parms=p, method="lsoda", hmax=TINFU/2)), error=function(e) NULL)
  if (is.null(sol)) return(NULL)
  sol$time_d <- sol$time/24
  idx <- sapply(obs_days, function(d) which.min(abs(sol$time_d-d)))
  fold <- sapply(lineages, function(L) sol[[L]][idx]/base_mod[L])
  data.frame(jour=obs_days, fold, check.names=FALSE)
}

# -- objectif : SSR sur fold-change (obs/baseline_ind vs pred/baseline_mod) --
objective <- function(theta) {
  p <- build_pars(); p$Slope_CMP<-exp(theta[1]); p$Slope_MEP<-exp(theta[2]); p$Slope_MPP<-exp(theta[1])*1.39
  ss<-0
  for (d in doses) {
    pf <- sim_fold(p,d); if(is.null(pf)) return(1e6)
    od <- obs_data[obs_data$dose_mgkg==d & obs_data$jour>=0, ]
    for (k in 1:nrow(od)) {
      j <- od$jour[k]; a <- od$Animal_Id[k]; pj <- pf[pf$jour==j, ]
      for (L in lineages) {
        bi <- base_ind(a,L); if(is.na(bi)||is.na(od[[L]][k])) next
        r <- od[[L]][k]/bi - pj[[L]]; ss <- ss + r^2
      }
    }
  }
  ss
}

theta0 <- log(c(init_pars$Slope_CMP, init_pars$Slope_MEP))
cat(sprintf("Slopes initiaux (main) : CMP=%.2f MEP=%.2f MPP=%.2f\n",
    init_pars$Slope_CMP, init_pars$Slope_MEP, init_pars$Slope_MPP))
cat(sprintf("SSR initiale (fold-change) = %.4f\n\n", objective(theta0)))

cat("Estimation (optim Nelder-Mead)...\n")
fit <- optim(par=log(c(0.5,0.6)), fn=objective, method="Nelder-Mead",
             control=list(maxit=400, reltol=1e-7))
S_CMP<-exp(fit$par[1]); S_MEP<-exp(fit$par[2])
cat(sprintf("\n=== SLOPES ESTIMES (fold-change) ===\n"))
cat(sprintf("  Slope_CMP = %.3f   (initial 0.57)\n", S_CMP))
cat(sprintf("  Slope_MEP = %.3f   (initial 0.66)\n", S_MEP))
cat(sprintf("  Slope_MPP = %.3f   (lie 1.39 x CMP)\n", S_CMP*1.39))
cat(sprintf("  SSR finale = %.4f  (vs init %.4f)  conv=%d\n",
    fit$value, objective(theta0), fit$convergence))

# -- RMSE fold-change par lignee ------------------------------------------
p <- build_pars(); p$Slope_CMP<-S_CMP; p$Slope_MEP<-S_MEP; p$Slope_MPP<-S_CMP*1.39
cat("\nRMSE fold-change par lignee (modele estime vs obs) :\n")
for (L in lineages) {
  num<-0;n<-0
  for (d in doses) {
    pf<-sim_fold(p,d); od<-obs_data[obs_data$dose_mgkg==d & obs_data$jour>=0,]
    for (k in 1:nrow(od)) {
      bi<-base_ind(od$Animal_Id[k],L); if(is.na(bi)||is.na(od[[L]][k])) next
      pj<-pf[pf$jour==od$jour[k],]
      num<-num+(od[[L]][k]/bi - pj[[L]])^2; n<-n+1
    }
  }
  cat(sprintf("  %-5s : %.1f%% fold  (n=%d)\n", L, 100*sqrt(num/n), n))
}
