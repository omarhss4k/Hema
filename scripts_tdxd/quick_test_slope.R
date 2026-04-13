############################################################
# quick_test_slope.R — Test rapide : quel Slope_CMP donne G3-4 ?
# Un seul patient (pas d'IIV), pas de 6h → 12h
############################################################
library(deSolve)
setwd("/home/user/Hema/scripts_tdxd")
source("parameters_tdxd_rat.R")
source("parameters_tdxd_human.R")
source("pkpd_tdxd_rat.R")

pars_pd_hu <- init_pars
pars_pd_hu[["k_dam"]] <- NULL
pars_pd_hu[["k_rep"]] <- NULL
pars_base <- c(pars_pd_hu, tdxd_pars_hu)

state_pd_hu <- init_state[!names(init_state) %in% c("C1", "C2", "Damage")]
state0_hu   <- c(tdxd_hu_state0, state_pd_hu)
times_hu    <- seq(0, 126 * 24, by = 12)  # 12h pour rapidité

# Test des slopes élevés — patient typique (pas d'IIV)
slopes_test <- c(2, 3, 5, 8, 12, 20)

cat("Test patient typique (pas d'IIV) — Neut nadir par Slope_CMP :\n")
cat(sprintf("  %-12s  %-12s  %-10s  %s\n", "Slope_CMP", "Neut nadir", "% baisse", "Grade"))
cat(sprintf("  %-12s  %-12s  %-10s  %s\n", "─────────", "──────────","────────","─────"))

for (s in slopes_test) {
  pars_i <- pars_base
  pars_i$Slope_CMP <- s
  pars_i$rate_fun  <- make_tdxd_infusion(
    dose_mgkg = 5.4, BW_kg = 70, Tinfu_h = 1.5,
    interval_h = 21 * 24, n_cycles = 6
  )
  out <- tryCatch(
    suppressMessages(suppressWarnings(
      as.data.frame(lsoda(
        y = state0_hu, times = times_hu,
        func = pkpd_tdxd_fornari, parms = pars_i,
        rtol = 1e-4, atol = 1e-6, maxsteps = 500000
      ))
    )),
    error = function(e) { cat("  ERROR at Slope=", s, "\n"); NULL }
  )
  if (!is.null(out)) {
    n_nadir <- min(out$Neut, na.rm=TRUE)
    pct <- 100*(n_nadir - pars_base$Neut0)/pars_base$Neut0
    grade <- if(n_nadir < 0.5) "G4" else if(n_nadir < 1.0) "G3" else if(n_nadir < 1.5) "G2" else if(n_nadir < 2.0) "G1" else "G0"
    cat(sprintf("  %-12.1f  %-12.3f  %-10.1f  %s\n", s, n_nadir, pct, grade))
  }
}
