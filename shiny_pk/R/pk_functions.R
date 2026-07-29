library(PKNCA)
library(dplyr)
library(tidyr)

# ── NCA ───────────────────────────────────────────────────────────────────────

#' Run non-compartmental analysis via PKNCA
#' @param df         data.frame with PK observations
#' @param dose_col   column name for dose amount
#' @param time_col   column name for time
#' @param conc_col   column name for concentration
#' @param animal_col column name for subject / animal ID
#' @param route      "iv" or "ev" (extravascular)
#' @return data.frame: one row per animal with AUClast, AUCinf, Cmax, Tmax,
#'         t½, CL, Vz, MRT
run_nca <- function(df, dose_col, time_col, conc_col, animal_col, route = "iv") {

  route_pknca <- if (route == "iv") "intravascular" else "extravascular"

  conc_df <- data.frame(
    subject = df[[animal_col]],
    time    = df[[time_col]],
    conc    = df[[conc_col]],
    stringsAsFactors = FALSE
  )

  # One dosing row per subject at time 0
  dose_df <- df %>%
    group_by(subject_tmp = .data[[animal_col]]) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      subject = subject_tmp,
      dose    = .data[[dose_col]],
      time    = 0
    )

  conc_obj <- PKNCAconc(conc_df,  conc ~ time | subject)
  dose_obj <- PKNCAdose(dose_df,  dose ~ time | subject, route = route_pknca)
  pk_data  <- PKNCAdata(conc_obj, dose_obj)

  suppressMessages(results <- pk.nca(pk_data))

  params_wanted <- c(
    "auclast", "aucinf.obs", "cmax", "tmax",
    "half.life", "cl.obs", "vz.obs",
    "mrt.iv.last", "mrt.last"
  )

  res_df <- as.data.frame(results) %>%
    filter(PPTESTCD %in% params_wanted) %>%
    mutate(PPORRES = suppressWarnings(as.numeric(as.character(PPORRES)))) %>%
    select(subject, parameter = PPTESTCD, value = PPORRES) %>%
    # Harmonise MRT regardless of route-specific name
    mutate(parameter = dplyr::recode(parameter,
      "mrt.iv.last" = "MRT", "mrt.last" = "MRT")) %>%
    group_by(subject, parameter) %>%
    summarise(value = dplyr::first(value[!is.na(value)]), .groups = "drop") %>%
    pivot_wider(names_from = parameter, values_from = value)

  names(res_df)[names(res_df) == "subject"] <- animal_col
  res_df
}

# ── ODE systems (deSolve-compatible signatures) ───────────────────────────────

#' 1-compartment IV ODE  —  state: A1 (amount) ; parms: CL, V1
#' IV bolus initial condition: A1(0) = dose, observe C1 = A1/V1
pk_ode_1comp <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dA1 <- -(CL / V1) * A1
    list(c(dA1))
  })
}

#' 2-compartment IV ODE  —  state: A1, A2 ; parms: CL, V1, V2, Q
#' IV bolus initial condition: A1(0) = dose, A2(0) = 0, observe C1 = A1/V1
pk_ode_2comp <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    k10 <- CL / V1
    k12 <- Q  / V1
    k21 <- Q  / V2
    dA1 <- -(k10 + k12) * A1 + k21 * A2
    dA2 <-   k12        * A1 - k21 * A2
    list(c(dA1, dA2))
  })
}

# ── Analytical C(t) predictors (internal, used by fit + simulate) ─────────────

.pred_1comp <- function(dose, CL, V1, times) {
  (dose / V1) * exp(-(CL / V1) * times)
}

.pred_2comp <- function(dose, CL, V1, V2, Q, times) {
  k10  <- CL / V1
  k12  <- Q  / V1
  k21  <- Q  / V2
  disc <- (k10 + k12 + k21)^2 - 4 * k10 * k21
  if (disc <= 0 || is.na(disc)) return(rep(NA_real_, length(times)))
  s     <- sqrt(disc)
  alpha <- (k10 + k12 + k21 + s) / 2
  beta  <- (k10 + k12 + k21 - s) / 2
  if (abs(alpha - beta) < 1e-12)  return(rep(NA_real_, length(times)))
  A <- (dose / V1) * (alpha - k21) / (alpha - beta)
  B <- (dose / V1) * (k21  - beta) / (alpha - beta)
  A * exp(-alpha * times) + B * exp(-beta * times)
}

# ── Model fitting ─────────────────────────────────────────────────────────────

#' Fit 1- or 2-compartment IV PK model via nlminb (log-space optimisation)
#'
#' @param df         data.frame
#' @param dose_col   dose column name (can differ per animal)
#' @param time_col   time column name
#' @param conc_col   concentration column name
#' @param animal_col subject / animal ID column name
#' @param n_comp     1 or 2
#' @return list with:
#'   \item{params}   named vector of primary PK parameters
#'   \item{rse}      percent relative standard errors (%RSE) from Hessian
#'   \item{derived}  derived macro-constants (k10, k12, k21, alpha, beta, t½, Vss)
#'   \item{AIC, BIC} information criteria (log-normal residual model)
#'   \item{RSS}      residual sum of squared log-residuals at optimum
#'   \item{convergence, message} nlminb status
fit_pk_model <- function(df, dose_col, time_col, conc_col, animal_col,
                         n_comp = 2, init_params = NULL) {

  stopifnot(n_comp %in% c(1L, 2L))

  animals  <- unique(df[[animal_col]])
  doses    <- tapply(df[[dose_col]], df[[animal_col]], mean, na.rm = TRUE)

  df_clean <- df[df[[conc_col]] > 0 & is.finite(df[[conc_col]]), ]
  n_obs    <- nrow(df_clean)
  if (n_obs < n_comp * 2)
    stop("Too few valid observations for a ", n_comp, "-compartment model.")

  # ── Objective: sum of squared log-residuals (log-normal error, all animals) ──
  make_obj <- function(nc) {
    function(log_theta) {
      theta <- exp(log_theta)
      rss   <- 0
      for (id in animals) {
        sub    <- df_clean[df_clean[[animal_col]] == id, ]
        if (nrow(sub) == 0) next
        t_obs  <- sub[[time_col]]
        C_obs  <- sub[[conc_col]]
        di     <- unname(doses[as.character(id)])
        C_pred <- if (nc == 1L)
          .pred_1comp(di, theta[1], theta[2], t_obs)
        else
          .pred_2comp(di, theta[1], theta[2], theta[3], theta[4], t_obs)
        ok <- is.finite(C_pred) & C_pred > 0
        if (!any(ok)) return(1e10)
        rss <- rss + mean((log(C_obs[ok]) - log(C_pred[ok]))^2)
      }
      rss
    }
  }

  # ── Auto-initialise from data ──
  cmax_est  <- max(df_clean[[conc_col]], na.rm = TRUE)
  dose_mean <- mean(unlist(doses))
  # Scale up: for IV bolus, observed Cmax < true C0 due to early distribution
  V1_init   <- dose_mean / cmax_est * 2

  # Terminal β slope: use only the later 50% of the time course to avoid
  # contamination from the rapid α (distribution) phase
  t_all    <- df_clean[[time_col]]
  t_cutoff <- quantile(t_all, 0.5)
  df_term  <- df_clean[t_all >= t_cutoff, , drop = FALSE]
  if (nrow(df_term) < 3) df_term <- df_clean

  beta_est <- tryCatch({
    lm_df  <- data.frame(t    = df_term[[time_col]],
                         logc = log(df_term[[conc_col]]))
    lm_fit <- lm(logc ~ t, data = lm_df)
    slope  <- coef(lm_fit)[["t"]]
    if (is.finite(slope) && slope < 0) abs(slope) else NULL
  }, error = function(e) NULL)

  CL_init <- if (!is.null(beta_est)) beta_est * V1_init else V1_init * 0.05

  if (n_comp == 1L) {
    param_names <- c("CL", "V1")
    n_p         <- 2L
    log_starts  <- list(log(c(CL_init, V1_init)))
  } else {
    param_names <- c("CL", "V1", "V2", "Q")
    n_p         <- 4L
    # Grid over k12 = Q/V1 spanning 3 orders of magnitude.
    # Antibody PK can have fast distribution (t½α ~ 1h → k12 ~ 0.5 h⁻¹)
    # so we must try large k12 values, not just fractions of k10.
    log_starts <- lapply(c(0.005, 0.05, 0.5, 5.0), function(k12) {
      log(c(CL_init, V1_init, V1_init * 2, k12 * V1_init))
    })
  }

  obj_fn <- make_obj(n_comp)

  # Add user-supplied starting point if provided
  if (!is.null(init_params)) {
    log_init_user <- if (n_comp == 1L)
      log(pmax(c(init_params$CL, init_params$V1), 1e-12))
    else
      log(pmax(c(init_params$CL, init_params$V1,
                 init_params$V2, init_params$Q), 1e-12))
    log_starts <- c(log_starts, list(log_init_user))
  }

  .run_nlminb <- function(log_start) {
    tryCatch(
      nlminb(start     = log_start,
             objective = obj_fn,
             control   = list(eval.max = 3000, iter.max = 1500,
                              rel.tol = 1e-12, x.tol = 1e-12)),
      error = function(e)
        list(objective = Inf, convergence = 99L,
             par = log_start, message = conditionMessage(e))
    )
  }

  opts <- lapply(log_starts, .run_nlminb)
  opt  <- opts[[which.min(sapply(opts, `[[`, "objective"))]]

  params <- setNames(exp(opt$par), param_names)

  # ── Hessian → %RSE ──
  # mean()-based objective deflates the Hessian by n_obs/n_groups.
  # Rescale to the equivalent sum-based Hessian before inverting.
  n_groups <- length(animals)
  n_eff    <- n_obs / n_groups          # avg observations per group
  H <- tryCatch(
    optimHess(opt$par, obj_fn),
    error = function(e) matrix(NA_real_, n_p, n_p)
  )
  rse <- tryCatch({
    H_scaled <- H * n_eff               # rescale to sum-based curvature
    se_log   <- sqrt(abs(diag(solve(H_scaled))))
    setNames(se_log * 100, paste0("%RSE_", param_names))
  }, error = function(e)
    setNames(rep(NA_real_, n_p), paste0("%RSE_", param_names))
  )

  # ── AIC / BIC (log-normal residual MLE) ──
  rss    <- opt$objective
  sigma2 <- rss / n_obs
  logL   <- -n_obs / 2 * (1 + log(2 * pi * sigma2))
  aic    <- -2 * logL + 2  * n_p
  bic    <- -2 * logL + n_p * log(n_obs)

  list(
    params      = params,
    rse         = rse,
    derived     = .compute_derived(params, n_comp),
    AIC         = aic,
    BIC         = bic,
    RSS         = rss,
    n_obs       = n_obs,
    convergence = opt$convergence,
    message     = opt$message
  )
}

# ── Derived macro-constants ───────────────────────────────────────────────────

.compute_derived <- function(params, n_comp) {
  CL  <- unname(params["CL"])
  V1  <- unname(params["V1"])
  k10 <- CL / V1

  if (n_comp == 1L) {
    return(list(k10 = k10, t_half = log(2) / k10))
  }

  V2  <- unname(params["V2"])
  Q   <- unname(params["Q"])
  k12 <- Q / V1
  k21 <- Q / V2
  s   <- sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21)

  alpha <- (k10 + k12 + k21 + s) / 2
  beta  <- (k10 + k12 + k21 - s) / 2

  list(
    k10          = k10,
    k12          = k12,
    k21          = k21,
    alpha        = alpha,
    beta         = beta,
    t_half_alpha = log(2) / alpha,
    t_half_beta  = log(2) / beta,
    Vss          = V1 + V2
  )
}

# ── Individual fitting ────────────────────────────────────────────────────────

#' Fit PK model individually for each animal
#' @return data.frame : one row per animal with params, %RSE, AIC, BIC
fit_pk_individual <- function(df, dose_col, time_col, conc_col, animal_col,
                              n_comp = 2, init_params = NULL) {
  animals <- unique(df[[animal_col]])

  rows <- lapply(animals, function(id) {
    sub <- df[df[[animal_col]] == id, ]
    fit <- tryCatch(
      fit_pk_model(sub, dose_col, time_col, conc_col, animal_col,
                   n_comp = n_comp, init_params = init_params),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      row <- data.frame(animal = as.character(id), stringsAsFactors = FALSE)
    } else {
      rse_clean <- setNames(unname(fit$rse), gsub("^%RSE_", "", names(fit$rse)))
      row <- data.frame(
        animal      = as.character(id),
        CL          = unname(fit$params["CL"]),
        V1          = unname(fit$params["V1"]),
        stringsAsFactors = FALSE
      )
      if (n_comp == 2L) {
        row$V2 <- unname(fit$params["V2"])
        row$Q  <- unname(fit$params["Q"])
      }
      derived <- unlist(fit$derived)
      for (nm in names(derived)) row[[nm]] <- unname(derived[nm])
      row$AIC  <- fit$AIC
      row$BIC  <- fit$BIC
      row$`%RSE_CL` <- unname(rse_clean["CL"])
      row$`%RSE_V1` <- unname(rse_clean["V1"])
      if (n_comp == 2L) {
        row$`%RSE_V2` <- unname(rse_clean["V2"])
        row$`%RSE_Q`  <- unname(rse_clean["Q"])
      }
    }
    row
  })

  do.call(dplyr::bind_rows, rows)
}

# ── Simulation ────────────────────────────────────────────────────────────────

#' Simulate IV bolus PK profile via analytical solution (exact for 1- and 2-cmt)
#' @param params named vector: CL, V1 [, V2, Q]
#' @param dose   IV bolus dose
#' @param times  numeric vector of sampling times
#' @param n_comp 1 or 2
#' @return data.frame(time, conc)
simulate_pk <- function(params, dose, times, n_comp = 2) {
  stopifnot(n_comp %in% c(1L, 2L))
  conc <- if (n_comp == 1L)
    .pred_1comp(dose, params["CL"], params["V1"], times)
  else
    .pred_2comp(dose, params["CL"], params["V1"],
                params["V2"], params["Q"], times)
  data.frame(time = times, conc = conc)
}
