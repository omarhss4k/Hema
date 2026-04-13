############################################################
# quick_scan_n300.R — Scan rapide N=100/slope pour trouver
# la valeur exacte donnant G3-4 ~20% en population N=300
############################################################
library(deSolve)
setwd("/home/user/Hema/scripts_tdxd")
source("parameters_tdxd_rat.R")
source("parameters_tdxd_human.R")
source("pkpd_tdxd_rat.R")

set.seed(42)

pars_pd_hu <- init_pars; pars_pd_hu[["k_dam"]] <- NULL; pars_pd_hu[["k_rep"]] <- NULL
pars_base <- c(pars_pd_hu, tdxd_pars_hu)
state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)
times_hu    <- seq(0, 126 * 24, by = 12)

omega_CL <- 0.35; omega_V1 <- 0.20; omega_SCMP <- 0.33
N <- 100

run_pop <- function(slope_typ) {
  set.seed(123)
  g34 <- 0; n_ok <- 0
  for (i in 1:N) {
    pars_i <- pars_base
    pars_i$CL_ADC    <- pars_base$CL_ADC  * exp(rnorm(1,0,omega_CL))
    pars_i$V1_ADC    <- pars_base$V1_ADC  * exp(rnorm(1,0,omega_V1))
    pars_i$Slope_CMP <- slope_typ         * exp(rnorm(1,0,omega_SCMP))
    pars_i$rate_fun  <- make_tdxd_infusion(5.4,70,1.5,21*24,6)
    out <- tryCatch(suppressMessages(suppressWarnings(as.data.frame(lsoda(
      y=state0_hu, times=times_hu, func=pkpd_tdxd_fornari, parms=pars_i,
      rtol=1e-4, atol=1e-6, maxsteps=500000)))), error=function(e) NULL)
    if (!is.null(out)) { n_ok <- n_ok+1; if(min(out$Neut,na.rm=T)<1.0) g34 <- g34+1 }
  }
  list(g34_pct = g34/n_ok*100, n_ok = n_ok)
}

slopes <- c(8, 10, 12)
cat(sprintf("%-10s  %-10s\n", "Slope_CMP", "G3-4 %"))
cat(sprintf("%-10s  %-10s\n", "─────────", "──────"))
for (s in slopes) {
  cat(sprintf("Slope=%2d ...", s)); flush.console()
  r <- run_pop(s)
  cat(sprintf("\r%-10d  %-10.1f  (n_ok=%d)\n", s, r$g34_pct, r$n_ok))
}
