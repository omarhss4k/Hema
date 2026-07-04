############################################################
# estimate_slopes_nhp.R  --  TEST de la technique d'estimation
# Estime Slope_CMP / Slope_MEP du composé interne par ajustement
# du modèle Fornari NHP aux COMPTAGES hemato observes (4 doses,
# 4 lignees). Meme principe que nlmixr2 (minimisation modele-donnees),
# via optim() sur le modele deSolve existant.
#
# CONFIDENTIEL : donnees composé interne -- usage interne uniquement.
############################################################
library(deSolve)
source("parameters_fgfr2_nhp.R")
source("parameters_nhp.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("pkpd_fgfr2_nhp.R")
source("data_nhp.R")     # -> obs_data (Neut, Plt, RBC, Ret par animal/dose/jour)

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
tmax_h <- (max(obs_days)+3)*24
times <- seq(0, tmax_h, by=3)

base <- c(Neut=init_pars$Neut0, Plt=init_pars$Plt0, RBC=init_pars$RBC0, Ret=init_pars$Ret0)

# -- simulate une dose, retourne predictions aux jours observes -----------
sim_pred <- function(p, dose) {
  p$rate_fun <- make_nhp_infusion(dose_mgkg=dose, BW_kg=BW, Tinfu_h=TINFU,
                  interval_h=fgfr2_nhp$interval_h, n_cycles=1)
  sol <- tryCatch(as.data.frame(ode(y=nhp_state0, times=times, func=pkpd_nhp_ode,
           parms=p, method="lsoda", hmax=TINFU/2)), error=function(e) NULL)
  if (is.null(sol)) return(NULL)
  sol$time_d <- sol$time/24
  idx <- sapply(obs_days, function(d) which.min(abs(sol$time_d - d)))
  data.frame(jour=obs_days, Neut=sol$Neut[idx], Plt=sol$Plt[idx],
             RBC=sol$RBC[idx], Ret=sol$Ret[idx])
}

# -- objectif : SSR relative (normalisee par baseline lignee) -------------
lineages <- c("Neut","Plt","RBC","Ret")
objective <- function(theta) {
  p <- build_pars()
  p$Slope_CMP <- exp(theta[1]); p$Slope_MEP <- exp(theta[2])
  p$Slope_MPP <- exp(theta[1]) * 1.39   # lie a CMP (ratio actuel), reduit dimension
  ss <- 0
  for (d in doses) {
    pr <- sim_pred(p, d); if (is.null(pr)) return(1e6)
    od <- obs_data[obs_data$dose_mgkg==d & obs_data$jour>=0, ]
    for (L in lineages) {
      m <- merge(od[,c("jour",L)], pr[,c("jour",L)], by="jour", suffixes=c(".o",".p"))
      r <- (m[[paste0(L,".o")]] - m[[paste0(L,".p")]]) / base[L]
      ss <- ss + sum(r^2, na.rm=TRUE)
    }
  }
  ss
}

cat("Slopes initiaux (calibres main) : CMP=%.2f MEP=%.2f MPP=%.2f\n")
cat(sprintf("  CMP=%.2f  MEP=%.2f  MPP=%.2f\n", init_pars$Slope_CMP, init_pars$Slope_MEP, init_pars$Slope_MPP))
cat(sprintf("SSR initiale = %.4f\n\n", objective(log(c(init_pars$Slope_CMP, init_pars$Slope_MEP)))))

cat("Estimation en cours (optim Nelder-Mead)...\n")
fit <- optim(par=log(c(0.5,0.6)), fn=objective, method="Nelder-Mead",
             control=list(maxit=300, reltol=1e-6))

S_CMP<-exp(fit$par[1]); S_MEP<-exp(fit$par[2])
cat(sprintf("\n=== SLOPES ESTIMES sur comptages NHP ===\n"))
cat(sprintf("  Slope_CMP = %.3f  (initial 0.57)\n", S_CMP))
cat(sprintf("  Slope_MEP = %.3f  (initial 0.66)\n", S_MEP))
cat(sprintf("  Slope_MPP = %.3f  (lie = 1.39 x CMP)\n", S_CMP*1.39))
cat(sprintf("  SSR finale = %.4f  (vs initiale %.4f)  convergence=%d\n",
    fit$value, objective(log(c(init_pars$Slope_CMP,init_pars$Slope_MEP))), fit$convergence))

# -- qualite d'ajustement par lignee (RMSE relatif) -----------------------
p <- build_pars(); p$Slope_CMP<-S_CMP; p$Slope_MEP<-S_MEP; p$Slope_MPP<-S_CMP*1.39
cat("\nRMSE relatif par lignee (modele estime vs obs) :\n")
for (L in lineages) {
  num<-0; n<-0
  for (d in doses) {
    pr<-sim_pred(p,d); od<-obs_data[obs_data$dose_mgkg==d & obs_data$jour>=0,]
    m<-merge(od[,c("jour",L)], pr[,c("jour",L)], by="jour", suffixes=c(".o",".p"))
    num<-num+sum(((m[[paste0(L,".o")]]-m[[paste0(L,".p")]])/base[L])^2,na.rm=TRUE)
    n<-n+sum(!is.na(m[[paste0(L,".o")]]))
  }
  cat(sprintf("  %-5s : %.1f%%  (n=%d)\n", L, 100*sqrt(num/n), n))
}
