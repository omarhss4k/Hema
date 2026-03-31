# =============================================================================
# 03_simulation_core.R
# Moteur de simulation populationnelle — Modèle PKPD Carboplatine
#
# Fonctionnalités :
#   1. Définition du schéma posologique (dose unique ou cycles multiples)
#   2. Construction des conditions initiales à l'état stationnaire
#   3. Simulation d'un patient individuel (ODE solver deSolve/lsoda)
#   4. Simulation d'une population (N patients avec IIV)
#
# Dépendances :
#   - deSolve   (CRAN)
#   - 01_model_ode.R
#   - 02_parameters_library.R
# =============================================================================

# Opérateur null-coalesce
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Détection du répertoire courant ---
.module_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "pkpd_model")

# --- Chargement des modules ---
source(file.path(.module_dir, "01_model_ode.R"))
source(file.path(.module_dir, "02_parameters_library.R"))

if (!requireNamespace("deSolve", quietly = TRUE))
  stop("Package 'deSolve' requis. Installez-le avec : install.packages('deSolve')")
library(deSolve)

# =============================================================================
# SECTION 1 — SCHÉMA POSOLOGIQUE
# =============================================================================

#' Crée un schéma posologique
#'
#' @param dose_mg    Dose en mg (calculée via Calvert si NULL)
#' @param AUC_target AUC cible (mg.min/mL) pour calcul Calvert
#' @param GFR        GFR (mL/min) pour calcul Calvert
#' @param schedule   Vecteur de temps d'administration (h), défaut c(0)
#' @param infusion_h Durée de la perfusion (h), défaut 0.5
#' @param species    "human" ou "rat"
#'
#' @return Liste : dose_umol, schedule_h, infusion_h, rate_umol_h
make_dosing <- function(dose_mg    = NULL,
                        AUC_target = 5,
                        GFR        = 125,
                        schedule   = c(0),
                        infusion_h = 0.5,
                        species    = "human") {

  if (is.null(dose_mg)) {
    if (species == "human") {
      dose_mg <- calc_calvert_dose(AUC_target, GFR)
    } else {
      stop("Pour le rat, spécifiez dose_mg explicitement (ex: 30 mg/kg * poids)")
    }
  }

  dose_umol     <- mg_to_umol(dose_mg)
  rate_umol_h   <- dose_umol / infusion_h

  list(
    dose_mg    = dose_mg,
    dose_umol  = dose_umol,
    schedule_h = schedule,
    infusion_h = infusion_h,
    rate_umol_h = rate_umol_h
  )
}

# =============================================================================
# SECTION 2 — CONDITIONS INITIALES À L'ÉTAT STATIONNAIRE
# =============================================================================

#' Construit le vecteur d'état initial Y0 à partir des paramètres
#'
#' @param p Liste de paramètres complète (issue de build_parameters())
#' @return Vecteur nommé de 23 valeurs
build_initial_state <- function(p) {

  c(
    Cen    = 0,
    Per    = 0,
    Damage = 0,

    MPP    = p$MPP0,
    CMP    = p$CMP0,
    MEP    = p$MEP0,

    T1Neut = p$T1N0,
    T2Neut = p$T1N0,   # T1=T2=T3 à l'état stationnaire
    T3Neut = p$T1N0,

    T1Mono = p$T1Mo0,
    T2Mono = p$T1Mo0,
    T3Mono = p$T1Mo0,

    T1Ret  = p$T1R0,
    T2Ret  = p$T2R0,
    T3Ret  = p$T3R0,

    T1Plt  = p$T1P0,
    T2Plt  = p$T2P0,
    T3Plt  = p$T3P0,

    Neut   = p$Neut0,
    Mono   = p$Mono0,
    Ret    = p$Ret0,
    Plt    = p$Plt0,
    RBC    = p$RBC0
  )
}

# =============================================================================
# SECTION 3 — SIMULATION INDIVIDUELLE
# =============================================================================

#' Simule la dynamique PKPD pour un patient
#'
#' @param p          Liste de paramètres (build_parameters())
#' @param dosing     Liste de schéma posologique (make_dosing())
#' @param t_end      Durée de simulation (h), défaut 30*24 = 720 h
#' @param t_res      Résolution temporelle (h), défaut 1 h
#' @param solver_opt Liste d'options deSolve (rtol, atol, maxsteps)
#'
#' @return data.frame avec colonnes : time + toutes les variables d'état
simulate_patient <- function(p,
                             dosing,
                             t_end      = 30 * 24,
                             t_res      = 1,
                             solver_opt = list(rtol    = 1e-6,
                                               atol    = 1e-8,
                                               maxsteps = 50000)) {

  # Vecteur de temps de sortie (heures)
  times <- seq(0, t_end, by = t_res)

  # Conditions initiales
  y0 <- build_initial_state(p)

  # --- Intégration de dose_fun dans parms (deSolve passe uniquement parms) ---
  p$dose_fun <- .make_dose_function(dosing)

  # --- Résolution ODE avec lsoda ---
  out <- tryCatch({
    ode(
      y        = y0,
      times    = times,
      func     = pkpd_ode,
      parms    = p,
      method   = "lsoda",
      rtol     = solver_opt$rtol,
      atol     = solver_opt$atol,
      maxsteps = solver_opt$maxsteps
    )
  }, error = function(e) {
    warning("ODE solver failed: ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(out)) return(NULL)

  df <- as.data.frame(out)
  names(df)[1] <- "time_h"
  df$time_d <- df$time_h / 24   # temps en jours pour les plots
  df
}


#' Crée la fonction de débit de perfusion à partir du schéma posologique
#'
#' @param dosing  Liste make_dosing()
#' @return Fonction f(t) → débit en µmol/h
.make_dose_function <- function(dosing) {

  t_starts <- dosing$schedule_h
  t_ends   <- t_starts + dosing$infusion_h
  rate     <- dosing$rate_umol_h

  function(t) {
    active <- any(t >= t_starts & t < t_ends)
    if (active) rate else 0
  }
}

# =============================================================================
# SECTION 4 — SIMULATION POPULATIONNELLE (IIV)
# =============================================================================

#' Génère les paramètres individuels d'un patient en appliquant l'IIV
#'
#' @param p_pop   Paramètres population (build_parameters())
#' @param seed_i  Seed aléatoire individuelle (pour reproductibilité)
#' @return Liste de paramètres individuels avec variabilité appliquée
.sample_individual_params <- function(p_pop, seed_i = NULL) {

  if (!is.null(seed_i)) set.seed(seed_i)

  cv <- IIV_params_CV
  p  <- p_pop

  # Variabilité log-normale : X_i = X_pop * exp(eta),  eta ~ N(0, omega^2)
  # omega ≈ sqrt(log(1 + CV^2))  (approximation log-normale)
  .lnorm_sample <- function(mu, cv_val) {
    omega <- sqrt(log(1 + cv_val^2))
    mu * exp(rnorm(1, 0, omega))
  }

  # Paramètres physiologiques avec IIV
  p$Neut0    <- .lnorm_sample(p_pop$Neut0,    cv$cv_Neut0)
  p$Mono0    <- .lnorm_sample(p_pop$Mono0,    cv$cv_Mono0)
  p$Plt0     <- .lnorm_sample(p_pop$Plt0,     cv$cv_Plt0)
  p$Ret0     <- .lnorm_sample(p_pop$Ret0,     cv$cv_Ret0)
  p$RBC0     <- .lnorm_sample(p_pop$RBC0,     cv$cv_RBC0)
  p$MPP0     <- .lnorm_sample(p_pop$MPP0,     cv$cv_MPP0)
  p$CMP0     <- .lnorm_sample(p_pop$CMP0,     cv$cv_CMP0)
  p$MEP0     <- .lnorm_sample(p_pop$MEP0,     cv$cv_MEP0)

  p$kcircNeut <- .lnorm_sample(p_pop$kcircNeut, cv$cv_kcircNeut)
  p$kcircMono <- .lnorm_sample(p_pop$kcircMono, cv$cv_kcircMono)
  p$kcircPlt  <- .lnorm_sample(p_pop$kcircPlt,  cv$cv_kcircPlt)
  p$kcircRBC  <- .lnorm_sample(p_pop$kcircRBC,  cv$cv_kcircRBC)

  p$MTTNeut  <- .lnorm_sample(p_pop$MTTNeut,  cv$cv_MTTNeut)
  p$MTTMono  <- .lnorm_sample(p_pop$MTTMono,  cv$cv_MTTMono)
  p$MTTPlt   <- .lnorm_sample(p_pop$MTTPlt,   cv$cv_MTTPlt)
  p$MTTRet   <- .lnorm_sample(p_pop$MTTRet,   cv$cv_MTTRet)

  # Paramètres drug avec IIV
  p$SlopeMPP  <- .lnorm_sample(p_pop$SlopeMPP,  cv$cv_SlopeMPP)
  p$SlopeCMP  <- .lnorm_sample(p_pop$SlopeCMP,  cv$cv_SlopeCMP)
  p$SlopeMEP  <- .lnorm_sample(p_pop$SlopeMEP,  cv$cv_SlopeMEP)
  p$delta_Ret <- .lnorm_sample(p_pop$delta_Ret, cv$cv_delta_Ret)
  p$delta_Plt <- .lnorm_sample(p_pop$delta_Plt, cv$cv_delta_Plt)

  # IIV sur PK (log-normale, omega de Zandvliet 2008)
  pk_iiv <- IIV_pk_omega
  p$CL   <- p_pop$CL   * exp(rnorm(1, 0, pk_iiv$omega_CL))
  p$VCen <- p_pop$VCen * exp(rnorm(1, 0, pk_iiv$omega_VCen))
  p$VPer <- p_pop$VPer * exp(rnorm(1, 0, pk_iiv$omega_VPer))
  p$QC   <- p_pop$QC   * exp(rnorm(1, 0, pk_iiv$omega_QC))

  # Recalcul des paramètres à l'état stationnaire avec les nouvelles valeurs
  p <- .compute_steady_state(p)

  p
}


#' Simule une population de N patients
#'
#' @param n_patients   Nombre de patients (défaut 1000)
#' @param species      "human" ou "rat"
#' @param dosing       Liste make_dosing()
#' @param t_end        Durée de simulation (h)
#' @param t_res        Résolution temporelle (h)
#' @param seed         Seed global pour reproductibilité
#' @param patient_cov  data.frame de covariables (colonnes: age,weight,creat,sex)
#'                     Si NULL, valeurs population médiane utilisées
#' @param verbose      Afficher la progression (défaut TRUE)
#'
#' @return Liste :
#'   $simulations : liste de data.frames (un par patient)
#'   $params      : liste de listes de paramètres individuels
simulate_population <- function(n_patients  = 1000,
                                species     = "human",
                                dosing      = make_dosing(),
                                t_end       = 60 * 24,
                                t_res       = 2,
                                seed        = 42,
                                patient_cov = NULL,
                                verbose     = TRUE) {

  set.seed(seed)

  # Paramètres population de base
  p_pop <- build_parameters(species = species)

  sims   <- vector("list", n_patients)
  params <- vector("list", n_patients)

  if (verbose) {
    message(sprintf("Simulation de %d patients (%s, dose=%.1f mg, t_end=%d j)...",
                    n_patients, species, dosing$dose_mg, round(t_end / 24)))
    pb <- txtProgressBar(min = 0, max = n_patients, style = 3)
  }

  for (i in seq_len(n_patients)) {

    # Covariables patient individuelles (si fournies)
    if (!is.null(patient_cov) && i <= nrow(patient_cov)) {
      pat <- as.list(patient_cov[i, ])
    } else {
      pat <- list(age = 60, weight = 70, creat = 80, sex = "M")
    }

    # Paramètres population pour ce patient
    p_base <- build_parameters(species = species, patient = pat)

    # Application de l'IIV
    p_i <- .sample_individual_params(p_base, seed_i = seed * 1000 + i)

    params[[i]] <- p_i

    # Simulation individuelle
    sim <- simulate_patient(
      p      = p_i,
      dosing = dosing,
      t_end  = t_end,
      t_res  = t_res
    )

    sims[[i]] <- sim

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    n_failed <- sum(sapply(sims, is.null))
    message(sprintf("Terminé. %d/%d simulations réussies.", n_patients - n_failed, n_patients))
  }

  list(simulations = sims, params = params)
}


# =============================================================================
# SECTION 5 — AGRÉGATION DES RÉSULTATS POPULATION
# =============================================================================

#' Calcule les percentiles population (5%, 25%, 50%, 75%, 95%) par variable
#'
#' @param pop_result  Résultat de simulate_population()
#' @param variables   Noms des colonnes à résumer
#'                    (défaut: c("Neut","Mono","Ret","Plt","RBC","Damage"))
#'
#' @return data.frame long avec colonnes : time_h, time_d, variable, p05, p25, p50, p75, p95
summarize_population <- function(pop_result,
                                 variables = c("Neut", "Mono", "Ret", "Plt", "RBC", "Damage")) {

  sims <- pop_result$simulations
  # Filtrer les simulations nulles
  sims <- Filter(Negate(is.null), sims)

  if (length(sims) == 0) stop("Aucune simulation disponible.")

  # Grille de temps commune
  time_grid <- sims[[1]]$time_h

  results <- lapply(variables, function(var) {

    # Matrice : lignes = temps, colonnes = patients
    mat <- do.call(cbind, lapply(sims, function(df) {
      if (is.null(df) || !(var %in% names(df))) return(rep(NA, length(time_grid)))
      df[[var]]
    }))

    # Percentiles par ligne (chaque instant)
    p_mat <- apply(mat, 1, function(row) {
      quantile(row, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    })

    data.frame(
      time_h   = time_grid,
      time_d   = time_grid / 24,
      variable = var,
      p05      = p_mat["5%",  ],
      p25      = p_mat["25%", ],
      p50      = p_mat["50%", ],
      p75      = p_mat["75%", ],
      p95      = p_mat["95%", ]
    )
  })

  do.call(rbind, results)
}
