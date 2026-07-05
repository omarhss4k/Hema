############################################################
# plot_obs_profiles.R  -- DIAGNOSTIC (usage interne, CONFIDENTIEL)
# Profils hemato OBSERVES (fold-change vs baseline individuelle)
# superposes a la prediction du modele Fornari (Slopes actuels).
############################################################
library(deSolve)
source("parameters_fgfr2_nhp.R")
source("parameters_nhp.R")
source("../shared/parameters_FORNARI_CORRECT.R")
source("pkpd_fgfr2_nhp.R")
source("data_nhp.R")

build_pars <- function() { p<-init_pars; for(nm in names(fgfr2_nhp)) p[[nm]]<-fgfr2_nhp[[nm]]
  p$k_dam_DXd<-0.017; p$rate_fun<-NULL; p }
pd_state <- init_state[setdiff(names(init_state),c("C1","C2","Damage"))]
nhp_state0 <- c(C_ADC1=0,C_ADC2=0,Damage=0,pd_state)
BW<-4; TINFU<-0.5
doses <- sort(unique(obs_data$dose_mgkg))
lineages <- c("Neut","Plt","RBC","Ret")
base_mod <- c(Neut=init_pars$Neut0,Plt=init_pars$Plt0,RBC=init_pars$RBC0,Ret=init_pars$Ret0)
pre <- obs_data[obs_data$jour<0,]
base_ind <- function(a,L){v<-pre[[L]][pre$Animal_Id==a];if(!length(v)||all(is.na(v)))return(NA_real_);mean(v,na.rm=TRUE)}

tmax <- max(obs_data$jour)+3
times <- seq(0, tmax*24, by=3)
sim_fold <- function(p,dose){
  p$rate_fun<-make_nhp_infusion(dose_mgkg=dose,BW_kg=BW,Tinfu_h=TINFU,interval_h=fgfr2_nhp$interval_h,n_cycles=1)
  sol<-as.data.frame(ode(y=nhp_state0,times=times,func=pkpd_nhp_ode,parms=p,method="lsoda",hmax=TINFU/2))
  sol$td<-sol$time/24; sol
}
p0 <- build_pars()   # Slopes actuels 0.57/0.66/0.79
sims <- lapply(doses, function(d) sim_fold(p0,d)); names(sims)<-as.character(doses)

cols <- c("#2166ac","#1a9641","#d6604d","#b2182b"); names(cols)<-as.character(doses)
dir.create("results",showWarnings=FALSE)
pdf("results/OBS_profiles_diagnostic.pdf", width=12, height=9)
par(mfrow=c(2,2), mar=c(4,4.5,3,1))
for (L in lineages) {
  # points observes fold-change
  fx<-c(); fy<-c(); fd<-c()
  for (k in 1:nrow(obs_data)) {
    if (obs_data$jour[k]<0) next
    bi<-base_ind(obs_data$Animal_Id[k],L); v<-obs_data[[L]][k]
    if(is.na(bi)||is.na(v)) next
    fx<-c(fx,obs_data$jour[k]); fy<-c(fy,v/bi); fd<-c(fd,as.character(obs_data$dose_mgkg[k]))
  }
  plot(fx,fy,pch=19,col=cols[fd],cex=1.2,xlab="Jour",ylab="Fold-change vs baseline",
       main=sprintf("%s  (obs points + modele lignes)",L),
       ylim=range(c(0.3,fy,1.3),na.rm=TRUE))
  abline(h=1,lty=2,col="grey50")
  for (d in as.character(doses)) {
    s<-sims[[d]]; lines(s$td, s[[L]]/base_mod[L], col=cols[d], lwd=2)
  }
  legend("topright", legend=paste0(doses," mg/kg"), col=cols, pch=19, lwd=2, bty="n", cex=0.8)
}
mtext("Profils hemato NHP -- OBSERVES (points) vs MODELE Fornari Slopes actuels (lignes)",
      outer=TRUE, line=-1.5, font=2, cex=1.05)
dev.off()
cat("-> results/OBS_profiles_diagnostic.pdf\n")
