############################################################
# run_bpa_predict_v2.R
# PREDICTION HEMATOTOXICITE BPA -- architecture repliquee de l'etape 5
# (T-DXd humain, VALIDEE vs FDA : Neut 15.7% vs 16.2%, Anemie 8.7% vs 9%).
#
# PRINCIPE (identique etape 5, seul le composant BPA change) :
#   - kill LINEAIRE  : kill_X = Slope_X * D_kill   (pas d'ED50_kill/Emax)
#   - mixture de sensibilite par lignee (proportions heritees T-DXd)
#   - driver Damage  : C_ADC1 (use_ADC_driver=TRUE), IC50_ADC=27.7 ug/mL
#
# SPECIFIQUE BPA :
#   - PK 1-cmt allometrique souris->cible (build_bpa_pars)
#   - Slopes = Slopes_tdxd(etape5) * ratio IC50   (potence relative FLT3)
#       ratio_CMP = IC50_CMP_tdxd/IC50_CMP_bpa   (myeloide, fort)
#       ratio_MEP = IC50_MEP_tdxd/IC50_MEP_bpa   (erythroide, faible)
#     -> selectivite myeloide (FLT3) portee par les Slopes, comme etape 5.
#
# HYPOTHESES (pas de donnees cliniques BPA) :
#   - proportions de mixture identiques au T-DXd (distribution de sensibilite
#     supposee propriete de la population hematopoietique, non du compose)
#   - IC50_ADC=27.7 ug/mL (echelle Damage, heritee T-DXd)
#   - thrombopenie directe : mecanisme T-DXd conserve (FLT3 faible sur
#     megacaryocytes -> probablement surestime, borne haute)
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("../etape5_tdxd_humain/parameters_tdxd_human.R")   # mixtures T-DXd + Slopes
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")        # make_tdxd_infusion
  source("pkpd_bpa_1cmt.R")
  source("../shared/bpa_build_params.R")
}))

# -- IIV (log-normal) : memes valeurs que l'etape 5 --------------------------
omega_CL        <- 0.35   # Yin 2020
omega_V1        <- 0.20   # Yin 2020
omega_Slope_CMP <- 0.33   # Fornari Table S4
omega_Slope_MEP <- 0.33   # Fornari Table S4

# -- Ratios de puissance IC50 (BPA vs T-DXd) ----------------------------------
IC50_MEP_tdxd <- 27.3 / 148000 * 1e6   # 184.5 nM
IC50_CMP_tdxd <- 28.1 / 148000 * 1e6   # 189.9 nM
ratio_MEP     <- IC50_MEP_tdxd / 155   # ~1.19  (erythroide)
ratio_CMP_s1  <- IC50_CMP_tdxd / 7.8   # ~24.3  (myeloide, sous-traitant)
ratio_CMP_s2  <- IC50_CMP_tdxd / 0.015 # ~12660 (myeloide, interne, extreme)

# -- Protocole NHP ------------------------------------------------------------
BW_KG <- 4; TINFU_H <- 1.5; INTERVAL_H <- 21*24; N_CYCLES <- 4
FREQ  <- sprintf("Q%dW", round(INTERVAL_H/(7*24)))

# -- Etat / grille ------------------------------------------------------------
state_pd <- init_state[!names(init_state) %in% c("C1","C2","Damage")]
state0   <- c(c(C_ADC1=0, C_DXd=0, C_DXd_ic=0, Damage=0), state_pd)
times    <- seq(0, N_CYCLES*INTERVAL_H + 7*24, by=6)

ctcae_neut   <- function(x) if(x<0.5)"G4" else if(x<1.0)"G3" else if(x<1.5)"G2" else if(x<2.0)"G1" else "G0"
ctcae_anemia <- function(x,x0){f<-x/x0; if(f<0.54)"G4" else if(f<0.67)"G3" else if(f<0.80)"G2" else if(f<0.90)"G1" else "G0"}
ctcae_plt    <- function(x) if(x<25)"G4" else if(x<50)"G3" else if(x<75)"G2" else if(x<150)"G1" else "G0"
grade_order  <- c("G0","G1","G2","G3","G4")

# -- Un scenario x une dose : N patients avec mixture (style etape 5) ---------
run_dose <- function(dose_mgkg, ratio_CMP, N=100, seed=42) {
  set.seed(seed)
  # PK BPA (mode LINEAIRE : ED50_kill=NULL -> supprime la branche Emax)
  p <- build_bpa_pars(BW_target_kg=BW_KG, IC50_CMP_nM=7.8, ED50_kill=NULL, verbose=FALSE)
  pars_pd <- init_pars; pars_pd[["k_dam"]]<-NULL; pars_pd[["k_rep"]]<-NULL
  pars_typ <- modifyList(pars_pd, p)

  gN<-gA<-gP<-character(N); ndN<-ndR<-ndP<-numeric(N)
  for (i in 1:N) {
    pi <- pars_typ
    pi$CL_ADC <- pars_typ$CL_ADC * exp(rnorm(1,0,omega_CL))
    pi$V1_ADC <- pars_typ$V1_ADC * exp(rnorm(1,0,omega_V1))

    # Mixture CMP (neutropenie) : proportions T-DXd, Slopes x ratio_CMP
    rc <- runif(1); pr <- 1 - p_sensitive_tdxd - p_moderate_tdxd
    S_CMP_base <- if (rc < pr) Slope_resist_tdxd
                  else if (rc < pr + p_moderate_tdxd) Slope_moderate_tdxd
                  else Slope_sensitive_tdxd
    pi$Slope_CMP <- max(0, S_CMP_base * ratio_CMP * exp(rnorm(1,0,omega_Slope_CMP)))

    # Mixture MEP (anemie) : proportions T-DXd, Slopes x ratio_MEP
    rm <- runif(1); pmr <- 1 - p_MEP_light_tdxd - p_sensitive_mep_tdxd
    S_MEP_base <- if (rm < pmr) Slope_MEP_resist_tdxd
                  else if (rm < pmr + p_MEP_light_tdxd) Slope_MEP_light_tdxd
                  else Slope_MEP_sensitive_tdxd
    pi$Slope_MEP <- max(0, S_MEP_base * ratio_MEP * exp(rnorm(1,0,omega_Slope_MEP)))

    # Mixture thrombopenie directe (mecanisme T-DXd conserve -- hypothese)
    pi$Slope_Plt_direct <- if (runif(1) < p_plt_susceptible)
                             max(0, Slope_Plt_direct_sus * exp(rnorm(1,0,omega_Slope_Plt_direct)))
                           else 0

    pi$rate_fun <- make_tdxd_infusion(dose_mgkg=dose_mgkg, BW_kg=BW_KG,
                    Tinfu_h=TINFU_H, interval_h=INTERVAL_H, n_cycles=N_CYCLES)
    out <- tryCatch(suppressMessages(suppressWarnings(as.data.frame(lsoda(
      y=state0, times=times, func=pkpd_bpa_fornari_1cmt, parms=pi,
      rtol=1e-5, atol=1e-7, maxsteps=5e5, hmax=TINFU_H/2)))), error=function(e) NULL)
    if (!is.null(out) && nrow(out)>10) {
      nn<-min(out$Neut,na.rm=TRUE); nr<-min(out$RBC,na.rm=TRUE); np<-min(out$Plt,na.rm=TRUE)
      gN[i]<-ctcae_neut(nn); gA[i]<-ctcae_anemia(nr,pars_typ$RBC0); gP[i]<-ctcae_plt(np)
      ndN[i]<-nn; ndR[i]<-nr; ndP[i]<-np
    } else {gN[i]<-NA}
  }
  ok<-!is.na(gN)
  list(neut=round(100*table(factor(gN[ok],levels=grade_order))/sum(ok)),
       anem=round(100*table(factor(gA[ok],levels=grade_order))/sum(ok)),
       plt =round(100*table(factor(gP[ok],levels=grade_order))/sum(ok)),
       neut_nadir=median(ndN,na.rm=TRUE))
}

# -- Balayage doses -----------------------------------------------------------
doses <- c(0.1, 0.25, 0.5, 1.0, 2.0, 3.6, 5.4)

print_scenario <- function(label, ratio_CMP) {
  cat(sprintf("\n\n=== PREDICTION BPA NHP -- %s ===\n", label))
  cat(sprintf("Architecture etape5 (mixture + kill lineaire) | %s x%d | N=100/dose\n", FREQ, N_CYCLES))
  cat(sprintf("%-8s | %-26s | %-26s | %-26s\n","Dose","NEUTROPENIE (G0..G4)","ANEMIE (G0..G4)","THROMBO (G0..G4)"))
  cat(paste(rep("-",96),collapse=""),"\n")
  for (d in doses) {
    r <- run_dose(d, ratio_CMP)
    f <- function(t) paste(sprintf("%2.0f",t),collapse="/")
    g34 <- function(t) t["G3"]+t["G4"]
    cat(sprintf("%-8s | %-16s G34=%3.0f%% | %-16s G34=%3.0f%% | %-16s G34=%3.0f%%\n",
        d, f(r$neut), g34(r$neut), f(r$anem), g34(r$anem), f(r$plt), g34(r$plt)))
  }
}

print_scenario("SCENARIO 1 (IC50_CMP=7.8 nM, sous-traitant)", ratio_CMP_s1)
print_scenario("SCENARIO 2 (IC50_CMP=0.015 nM, interne)",    ratio_CMP_s2)

cat("\n\nNOTE : architecture repliquee de l'etape 5 (validee FDA). Selectivite\n")
cat("myeloide (FLT3) via ratio IC50 CMP(x24) >> MEP(x1.2) -> neutropenie attendue dominante.\n")
