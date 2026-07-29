############################################################
# predict_hematotox_nhp.R
# LIVRABLE MISSION : prediction d'hematotoxicite NHP a partir
#   des donnees PK (Shiny souris -> allometrie singe) + IC50 (CFU)
#
# VERSION CORRIGEE (prediction quantitative defendable) :
#   [1] IC50 nM : conversion ug/mL correcte (pas de x1000)
#   [2] k_int   : = 0 (pas de cible HER2 pour BPA/FLT3)
#   [3] kill    : modele Emax borne (ED50_kill) au lieu du lineaire
#   [4] hmax    : borne le pas solveur (sinon cycles 2-4 sautes)
#
# Sortie : tableau dose x grade CTCAE pour les 2 scenarios IC50
############################################################
suppressMessages(suppressWarnings({
  library(deSolve)
  source("../etape2_carboplatin_humain/parameters_human.R")
  source("../shared/parameters_FORNARI_CORRECT.R")
  source("../etape5_tdxd_humain/parameters_tdxd_human.R")
  source("../etape3_tdxd_rat/parameters_tdxd_rat.R")
  source("pkpd_bpa_1cmt.R")
  source("parameters_bpa_nhp.R")   # base (valeurs buggees ecrasees ci-dessous)
}))

# ══════════════════════════════════════════════════════════
# [1] IC50 T-DXd correctes (nM) : C[nM] = C[ug/mL]*1e6/MW
# ══════════════════════════════════════════════════════════
IC50_MEP_tdxd_nM <- 27.3 / 148000 * 1e6   # 184.5 nM
IC50_CMP_tdxd_nM <- 28.1 / 148000 * 1e6   # 189.9 nM

# Slopes recalcules (ratio IC50 vs T-DXd)
Slope_MEP_BPA    <- Slope_MEP_sensitive_tdxd * (IC50_MEP_tdxd_nM / 155)     # ~1.19
Slope_CMP_BPA_s1 <- Slope_sensitive_tdxd     * (IC50_CMP_tdxd_nM / 7.8)     # ~1439
Slope_CMP_BPA_s2 <- Slope_sensitive_tdxd     * (IC50_CMP_tdxd_nM / 0.015)   # ~748068

# ══════════════════════════════════════════════════════════
# [2]+[3] PK/PD corrige
# ══════════════════════════════════════════════════════════
pars_pd <- init_pars; pars_pd[["k_dam"]] <- NULL; pars_pd[["k_rep"]] <- NULL
pars_base <- c(pars_pd, bpa_pars_nhp)
pars_base$k_int      <- 0        # [2] pas de cible HER2
pars_base$krel_power <- 0        #     relargage non decroissant (BPA != linker T-DXd)
pars_base$krel_factor<- 1.0
pars_base$ED50_kill  <- 0.10     # [3] kill Emax borne

bpa_state0 <- c(C_ADC1 = 0, C_DXd = 0, C_DXd_ic = 0, Damage = 0)
state_pd   <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0     <- c(bpa_state0, state_pd)

times   <- seq(0, (N_CYCLES * INTERVAL_H + 7 * 24), by = 6)
n_times <- length(times)

ctcae_neut   <- function(x) if (x<0.5)"G4" else if (x<1.0)"G3" else if (x<1.5)"G2" else if (x<2.0)"G1" else "G0"
ctcae_anemia <- function(x,x0){f<-x/x0; if(f<0.54)"G4" else if(f<0.67)"G3" else if(f<0.80)"G2" else if(f<0.90)"G1" else "G0"}
ctcae_plt    <- function(x) if (x<25)"G4" else if (x<50)"G3" else if (x<75)"G2" else if (x<150)"G1" else "G0"
grade_order  <- c("G0","G1","G2","G3","G4")

# ══════════════════════════════════════════════════════════
# Fonction : un scenario x une dose -> distribution de grades
# ══════════════════════════════════════════════════════════
run_dose <- function(dose_mgkg, Slope_CMP_typ, N = 100, seed = 42) {
  set.seed(seed)
  gN <- gA <- gP <- character(N)
  nadN <- nadR <- numeric(N)
  for (i in 1:N) {
    pi <- pars_base
    pi$CL_ADC    <- pars_base$CL_ADC * exp(rnorm(1,0,omega_CL))
    pi$V1_ADC    <- pars_base$V1_ADC * exp(rnorm(1,0,omega_V1))
    pi$Slope_CMP <- max(0, Slope_CMP_typ * exp(rnorm(1,0,omega_Slope_CMP)))
    pi$Slope_MEP <- max(0, Slope_MEP_BPA * exp(rnorm(1,0,omega_Slope_MEP)))
    pi$rate_fun  <- make_tdxd_infusion(dose_mgkg=dose_mgkg, BW_kg=BW_KG,
                     Tinfu_h=TINFU_H, interval_h=INTERVAL_H, n_cycles=N_CYCLES)
    out <- tryCatch(suppressMessages(suppressWarnings(as.data.frame(lsoda(
      y=state0, times=times, func=pkpd_bpa_fornari_1cmt, parms=pi,
      rtol=1e-5, atol=1e-7, maxsteps=5e5, hmax=TINFU_H/2)))), error=function(e) NULL)  # [4]
    if (!is.null(out) && nrow(out) > 10) {
      nn <- min(out$Neut,na.rm=TRUE); nr <- min(out$RBC,na.rm=TRUE); np <- min(out$Plt,na.rm=TRUE)
      gN[i]<-ctcae_neut(nn); gA[i]<-ctcae_anemia(nr,pars_base$RBC0); gP[i]<-ctcae_plt(np)
      nadN[i]<-nn; nadR[i]<-nr
    } else { gN[i]<-NA; gA[i]<-NA; gP[i]<-NA; nadN[i]<-NA; nadR[i]<-NA }
  }
  ok <- !is.na(gN)
  list(
    n_ok = sum(ok),
    neut = round(100*table(factor(gN[ok],levels=grade_order))/sum(ok)),
    anem = round(100*table(factor(gA[ok],levels=grade_order))/sum(ok)),
    plt  = round(100*table(factor(gP[ok],levels=grade_order))/sum(ok)),
    neut_nadir_med = median(nadN,na.rm=TRUE)
  )
}

# ══════════════════════════════════════════════════════════
# Balayage de doses x 2 scenarios
# ══════════════════════════════════════════════════════════
doses <- c(0.1, 0.25, 0.5, 1.0, 2.0, 3.6, 5.4)

print_scenario <- function(label, Slope_CMP_typ) {
  cat("\n\n╔══════════════════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  PREDICTION HEMATOTOXICITE NHP -- %-36s║\n", label))
  cat(sprintf("║  Slope_CMP=%.1f  Slope_MEP=%.2f  |  %s x%d  N=100/dose%s║\n",
      Slope_CMP_typ, Slope_MEP_BPA,
      sprintf("Q%dW", round(INTERVAL_H/(7*24))), N_CYCLES, strrep(" ", 10)))
  cat("╚══════════════════════════════════════════════════════════════════════╝\n")
  cat(sprintf("%-9s | %-22s | %-22s | %-22s\n",
      "Dose", "NEUTROPENIE G3-4", "ANEMIE G3-4", "THROMBO G3-4"))
  cat(sprintf("%-9s | %-22s | %-22s | %-22s\n",
      "(mg/kg)", "(G0/G1/G2/G3/G4)", "(G0/G1/G2/G3/G4)", "(G0/G1/G2/G3/G4)"))
  cat(paste(rep("-",85), collapse=""), "\n")
  for (d in doses) {
    r <- run_dose(d, Slope_CMP_typ)
    fmt <- function(t) paste(sprintf("%2.0f", t), collapse="/")
    gr <- function(t) t["G3"]+t["G4"]
    cat(sprintf("%-9s | %-14s G3-4=%3.0f%% | %-14s G3-4=%3.0f%% | %-14s G3-4=%3.0f%%\n",
        d, fmt(r$neut), gr(r$neut), fmt(r$anem), gr(r$anem), fmt(r$plt), gr(r$plt)))
  }
}

print_scenario("SCENARIO 1 (IC50_CMP=7.8 nM, sous-traitant)", Slope_CMP_BPA_s1)
print_scenario("SCENARIO 2 (IC50_CMP=0.015 nM, interne)",    Slope_CMP_BPA_s2)

cat("\n\nNOTE : prediction quantitative -- version corrigee (IC50, k_int=0, kill Emax, hmax)\n")
cat("Toxicite limitante attendue : ANEMIE (RBC vie longue -> cumul inter-cycles)\n")
