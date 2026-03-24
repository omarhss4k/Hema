############################################################
# friberg_schmitt.R
# Modèle de Friberg (Clin Pharmacol Ther 2002) pour la
# simulation de la Figure 4c de Fornari 2019 (Supp. S12)
#
# S12 : "we compared our model results with simulated
# clinical data generated with Friberg model and parameter
# values reported in [Schmitt 2010, réf. 5]"
#
# PK  : identique à Fornari — 2 compartiments, carboplatin total mg/L
#       CL ajustée au GFR=125 (= 150 mL/min → 9 L/h)
# PD  : Friberg standard — Prol + 3 Transit + Circ pour Neut et Plt
#       Drug effect : E(t) = Slope × Damage(t)   [Damage = proxy adducts ADN]
#       Damage partagé avec le modèle QSP Fornari (k_dam = k_rep = 0.017 /h)
#
# NOTE sur le driver Damage (pas C_free instantané) :
#   Avec la PK rapide de Fornari (t1/2_β ≈ 3.6h), le carboplatin libre
#   est éliminé en ~24h → E(t) = Slope × C_free(t) donne nadir Grade 0 car
#   l'exposition est trop brève vs MTT_N = 152.9h.
#   En utilisant Damage (t1/2 = 41h), l'effet persiste ~5 jours,
#   produisant un nadir G1-G3 selon le Slope individuel — cohérent avec
#   les observations cliniques de Schmitt 2010.
#
# Slopes calibrés (non directement de Schmitt 2010) :
#   Slope_N = 0.20 /µM  →  nadir déterministe G1 (1.96)  →  G3/G4 ≈ 17% pop.
#   Slope_P = 0.40 /µM  →  nadir déterministe G2 (66)    →  Plt distribution
############################################################

# ══════════════════════════════════════════════════════════
# PARAMÈTRES FRIBERG / SCHMITT 2010
# ══════════════════════════════════════════════════════════
schmitt_pars <- list(

  # ── PK (Zandvliet 2008, identique à parameters_human.R) ──
  V1          = 7.87,             # L
  V2          = 8.06,             # L
  Q           = 1.98,             # L/h
  mgL_to_uM   = 1000 / 371.25,   # 2.694 µM/mg/L

  # ── Damage (même que Fornari QSP, Table 1) ──
  k_dam = 0.017,   # /h (formation adducts)
  k_rep = 0.017,   # /h (réparation)

  # ── PD neutrophiles (Schmitt 2010 Table III) ──
  MTT_N   = 152.9,   # h — mean transit time
  gamma_N = 0.209,   # feedback exponent
  Circ0_N = 4.57,    # ×10⁹/L — ANC baseline moyen
  Slope_N = 0.20,    # /µM — calibré sur Damage (nadir G1 → G3/G4≈17%)

  # ── PD plaquettes (Schmitt 2010 Table III) ──
  MTT_P   = 175.3,   # h
  gamma_P = 0.686,   # feedback exponent
  Circ0_P = 248.0,   # ×10⁹/L — PLT baseline moyen (Schmitt 2010)
  Slope_P = 0.40,    # /µM — calibré sur Damage

  # ── IIV (Schmitt 2010 / De Carlo 2025) — SD log-normale ──
  omega_Circ0_N = 0.326,
  omega_Slope_N = 0.624,
  omega_Circ0_P = 0.268,
  omega_Slope_P = 0.547,

  # Placeholder
  rate_fun = NULL
)

# ══════════════════════════════════════════════════════════
# ODE FRIBERG — Damage-driven drug effect
# État : (C1, C2, Damage,
#         Prol_N, T1_N, T2_N, T3_N, Circ_N,
#         Prol_P, T1_P, T2_P, T3_P, Circ_P)
# ══════════════════════════════════════════════════════════
friberg_ode <- function(time, state, pars) {
  with(as.list(c(state, pars)), {

    # ── PK 2-compartiments ──
    rate_in <- if (!is.null(pars$rate_fun)) pars$rate_fun(time) else 0
    dC1 <- rate_in / V1 - (CL / V1) * C1 - (Q / V1) * C1 + (Q / V2) * C2
    dC2 <- (Q / V1) * C1 - (Q / V2) * C2

    # ── Damage (adducts ADN, t1/2 = 41h) ──
    dDamage <- k_dam * C1 * mgL_to_uM - k_rep * Damage

    # ── Drug effect via Damage ──
    E_N <- min(1.0, Slope_N * Damage)
    E_P <- min(1.0, Slope_P * Damage)

    # ── Neutrophiles — Friberg 4 compartiments ──
    k_tr_N  <- 4.0 / MTT_N
    FB_N    <- (Circ0_N / max(Circ_N, 0.001))^gamma_N
    dProl_N <- k_tr_N * Prol_N * (1 - E_N) * FB_N - k_tr_N * Prol_N
    dT1_N   <- k_tr_N * (Prol_N - T1_N)
    dT2_N   <- k_tr_N * (T1_N   - T2_N)
    dT3_N   <- k_tr_N * (T2_N   - T3_N)
    dCirc_N <- k_tr_N * T3_N - k_tr_N * Circ_N

    # ── Plaquettes — Friberg 4 compartiments ──
    k_tr_P  <- 4.0 / MTT_P
    FB_P    <- (Circ0_P / max(Circ_P, 0.001))^gamma_P
    dProl_P <- k_tr_P * Prol_P * (1 - E_P) * FB_P - k_tr_P * Prol_P
    dT1_P   <- k_tr_P * (Prol_P - T1_P)
    dT2_P   <- k_tr_P * (T1_P   - T2_P)
    dT3_P   <- k_tr_P * (T2_P   - T3_P)
    dCirc_P <- k_tr_P * T3_P - k_tr_P * Circ_P

    list(c(dC1, dC2, dDamage,
           dProl_N, dT1_N, dT2_N, dT3_N, dCirc_N,
           dProl_P, dT1_P, dT2_P, dT3_P, dCirc_P))
  })
}

# ══════════════════════════════════════════════════════════
# ÉTAT INITIAL — steady-state avant traitement
# ══════════════════════════════════════════════════════════
make_friberg_init_state <- function(pars) {
  c(
    C1      = 0,
    C2      = 0,
    Damage  = 0,
    Prol_N  = pars$Circ0_N,
    T1_N    = pars$Circ0_N,
    T2_N    = pars$Circ0_N,
    T3_N    = pars$Circ0_N,
    Circ_N  = pars$Circ0_N,
    Prol_P  = pars$Circ0_P,
    T1_P    = pars$Circ0_P,
    T2_P    = pars$Circ0_P,
    T3_P    = pars$Circ0_P,
    Circ_P  = pars$Circ0_P
  )
}

# ══════════════════════════════════════════════════════════
# SIMULATION POPULATION — Figure 4c (Supp. S12)
# Même protocole que S11 : AUC=5, GFR=125 fixe pour tous
# ══════════════════════════════════════════════════════════
simulate_friberg_population <- function(
    n_patients  = 1000,
    auc_target  = 5,
    gfr_fixed   = 125,
    n_cycles    = 1,
    interval_h  = 21 * 24,
    seed        = 42,
    verbose     = TRUE
) {
  library(deSolve)
  set.seed(seed)

  CL_fixed <- (gfr_fixed + 25) * 60 / 1000   # mL/min → L/h
  dose_mg  <- auc_target * (gfr_fixed + 25)   # Calvert
  if (verbose) cat(sprintf(
    "  Friberg S12 : n=%d, AUC=%g, GFR=%g, dose=%.0f mg, CL=%.2f L/h\n",
    n_patients, auc_target, gfr_fixed, dose_mg, CL_fixed))

  times <- seq(0, n_cycles * interval_h + 4 * 7 * 24, by = 1)

  eta_Circ0_N <- rnorm(n_patients, 0, schmitt_pars$omega_Circ0_N)
  eta_Slope_N <- rnorm(n_patients, 0, schmitt_pars$omega_Slope_N)
  eta_Circ0_P <- rnorm(n_patients, 0, schmitt_pars$omega_Circ0_P)
  eta_Slope_P <- rnorm(n_patients, 0, schmitt_pars$omega_Slope_P)

  grade_neut <- integer(n_patients)
  grade_plt  <- integer(n_patients)

  for (i in seq_len(n_patients)) {
    pars_i <- schmitt_pars
    pars_i$CL      <- CL_fixed
    pars_i$Circ0_N <- schmitt_pars$Circ0_N * exp(eta_Circ0_N[i])
    pars_i$Slope_N <- schmitt_pars$Slope_N  * exp(eta_Slope_N[i])
    pars_i$Circ0_P <- schmitt_pars$Circ0_P * exp(eta_Circ0_P[i])
    pars_i$Slope_P <- schmitt_pars$Slope_P  * exp(eta_Slope_P[i])
    pars_i$rate_fun <- make_repeated_infusion(
      dose_mg    = dose_mg,
      Tinfu_h    = 1,
      interval_h = interval_h,
      n_cycles   = n_cycles
    )
    state_i <- make_friberg_init_state(pars_i)

    tryCatch({
      out_i <- as.data.frame(lsoda(
        y        = state_i,
        times    = times,
        func     = friberg_ode,
        parms    = pars_i,
        rtol     = 1e-6,
        atol     = 1e-8,
        maxsteps = 100000
      ))
      nadir_N <- min(out_i$Circ_N, na.rm = TRUE)
      nadir_P <- min(out_i$Circ_P, na.rm = TRUE)
      grade_neut[i] <- ifelse(nadir_N < 0.5, 4,
                       ifelse(nadir_N < 1.0, 3,
                       ifelse(nadir_N < 1.5, 2,
                       ifelse(nadir_N < 2.0, 1, 0))))
      grade_plt[i]  <- ifelse(nadir_P <  25, 4,
                       ifelse(nadir_P <  50, 3,
                       ifelse(nadir_P <  75, 2,
                       ifelse(nadir_P < 150, 1, 0))))
    }, error = function(e) {
      grade_neut[i] <<- NA
      grade_plt[i]  <<- NA
    })

    if (verbose && i %% 200 == 0)
      cat(sprintf("    %d/%d\n", i, n_patients))
  }

  pct_neut <- sapply(1:4, function(g) 100 * mean(grade_neut == g, na.rm = TRUE))
  pct_plt  <- sapply(1:4, function(g) 100 * mean(grade_plt  == g, na.rm = TRUE))

  if (verbose) {
    cat(sprintf("  Neut: G1=%.1f%% G2=%.1f%% G3=%.1f%% G4=%.1f%%\n",
                pct_neut[1], pct_neut[2], pct_neut[3], pct_neut[4]))
    cat(sprintf("  Plt:  G1=%.1f%% G2=%.1f%% G3=%.1f%% G4=%.1f%%\n",
                pct_plt[1],  pct_plt[2],  pct_plt[3],  pct_plt[4]))
  }

  list(pct_neut   = pct_neut,
       pct_plt    = pct_plt,
       grade_neut = grade_neut,
       grade_plt  = grade_plt)
}
