# =============================================================================
# 04_scoring_toxicity.R
# Calcul des nadirs et attribution des grades de toxicité NCI-CTCAE v5
#
# Ce module calcule :
#   1. Le nadir (valeur minimale post-traitement) pour chaque lignée
#   2. Le grade de toxicité NCI-CTCAE (0 à 4)
#   3. Le temps au nadir et le temps de récupération
#   4. La distribution des grades sur une population
#
# Seuils NCI-CTCAE v5.0 :
#   Neutrophiles (Neut, 10^9/L) :
#     Grade 0 : ≥ 2.0
#     Grade 1 : 1.5 – 1.99
#     Grade 2 : 1.0 – 1.49
#     Grade 3 : 0.5 – 0.99
#     Grade 4 : < 0.5
#
#   Plaquettes (Plt, 10^9/L) :
#     Grade 0 : ≥ 100
#     Grade 1 : 75 – 99
#     Grade 2 : 50 – 74
#     Grade 3 : 25 – 49
#     Grade 4 : < 25
#
#   Hémoglobine / RBC — approx par Ret (10^9/L) :
#     Grade 0 : ≥ 100  (RBC; on utilise RBC en 10^9/L → adaptez si Hgb)
#     Grade 1 : 80 – 99
#     Grade 2 : 65 – 79  (indication transfusionnelle)
#     Grade 3 : < 65     (transfusion urgente)
#     Grade 4 : conséquences engageant le pronostic vital (rare avec carbo seul)
#
# Références : NCI-CTCAE v5.0 (2017), De Carlo 2025
# =============================================================================

# =============================================================================
# SECTION 1 — SEUILS NCI-CTCAE
# =============================================================================

#' Seuils de toxicité NCI-CTCAE (10^9 cellules/L)
CTCAE_THRESHOLDS <- list(

  Neut = list(
    lineage = "Neutrophiles",
    unit    = "10^9/L",
    cuts    = c(2.0, 1.5, 1.0, 0.5),   # G0>=2.0, G1:1.5-1.99, G2:1.0-1.49, G3:0.5-0.99, G4:<0.5
    labels  = c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4")
  ),

  Plt = list(
    lineage = "Plaquettes",
    unit    = "10^9/L",
    cuts    = c(100, 75, 50, 25),       # G0>=100, G1:75-99, G2:50-74, G3:25-49, G4:<25
    labels  = c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4")
  ),

  Mono = list(
    lineage = "Monocytes",
    unit    = "10^9/L",
    cuts    = c(0.20, 0.10, 0.05, 0.01),
    labels  = c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4")
  ),

  Ret = list(
    lineage = "Réticulocytes",
    unit    = "10^9/L",
    cuts    = c(40, 30, 20, 10),        # seuils illustratifs
    labels  = c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4")
  ),

  RBC = list(
    lineage = "GR / Hémoglobine (approx.)",
    unit    = "10^9/L",
    cuts    = c(3800, 3200, 2700, 2000), # seuils en 10^9/L (Hgb équivalent)
    labels  = c("Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4")
  )
)

# =============================================================================
# SECTION 2 — ATTRIBUTION DU GRADE POUR UNE VALEUR
# =============================================================================

#' Attribue un grade NCI-CTCAE à une valeur cellulaire
#'
#' @param value    Valeur numérique (10^9/L)
#' @param lineage  Nom de la lignée : "Neut", "Plt", "Mono", "Ret", "RBC"
#'
#' @return Entier 0-4
assign_grade <- function(value, lineage = "Neut") {

  thresh <- CTCAE_THRESHOLDS[[lineage]]
  if (is.null(thresh)) stop("Lignée inconnue : ", lineage)

  cuts <- thresh$cuts   # vecteur décroissant de 4 seuils

  if      (value >= cuts[1]) return(0L)
  else if (value >= cuts[2]) return(1L)
  else if (value >= cuts[3]) return(2L)
  else if (value >= cuts[4]) return(3L)
  else                       return(4L)
}

#' Vectorise assign_grade sur un vecteur de valeurs
#'
#' @param values   Vecteur numérique
#' @param lineage  Nom de la lignée
#'
#' @return Vecteur d'entiers 0-4
assign_grade_vec <- function(values, lineage = "Neut") {
  vapply(values, assign_grade, integer(1), lineage = lineage)
}

# =============================================================================
# SECTION 3 — ANALYSE D'UN PROFIL INDIVIDUEL
# =============================================================================

#' Calcule les métriques de toxicité pour un patient
#'
#' @param sim_df     data.frame issu de simulate_patient()
#' @param lineages   Vecteur de lignées à analyser
#' @param t_start_h  Début de la fenêtre d'analyse (h après dose), défaut 0
#' @param t_end_h    Fin de la fenêtre d'analyse (h), NULL = toute la simulation
#' @param baseline_frac Fraction de la valeur basale définissant la récupération (0.8 = 80%)
#'
#' @return data.frame avec colonnes :
#'   lineage, nadir_value, nadir_time_h, nadir_time_d,
#'   grade_nadir, time_to_recovery_h, time_to_recovery_d, AUC_below_normal
score_patient <- function(sim_df,
                          lineages       = c("Neut", "Plt", "Ret", "RBC"),
                          t_start_h      = 0,
                          t_end_h        = NULL,
                          baseline_frac  = 0.80) {

  if (is.null(sim_df)) return(NULL)
  if (is.null(t_end_h)) t_end_h <- max(sim_df$time_h)

  # Fenêtre temporelle
  window <- sim_df[sim_df$time_h >= t_start_h & sim_df$time_h <= t_end_h, ]

  results <- lapply(lineages, function(lin) {

    if (!(lin %in% names(window))) {
      warning("Variable '", lin, "' absente du data.frame.")
      return(NULL)
    }

    vals  <- window[[lin]]
    times <- window$time_h

    # Valeur basale = valeur au t_start
    baseline <- vals[which.min(abs(times - t_start_h))]

    # Nadir : valeur minimale post-dose
    idx_nadir    <- which.min(vals)
    nadir_val    <- vals[idx_nadir]
    nadir_time_h <- times[idx_nadir]

    # Grade au nadir
    grade <- assign_grade(nadir_val, lineage = lin)

    # Temps de récupération : premier temps > baseline_frac * baseline APRÈS le nadir
    t_recovery_h <- NA_real_
    post_nadir   <- window[times >= nadir_time_h, ]
    recov_thresh <- baseline_frac * baseline
    recov_idx    <- which(post_nadir[[lin]] >= recov_thresh)
    if (length(recov_idx) > 0) {
      t_recovery_h <- post_nadir$time_h[recov_idx[1]]
    }

    # AUC sous le seuil normal (grade ≥ 1) — intégration trapézoïdale
    normal_thresh <- CTCAE_THRESHOLDS[[lin]]$cuts[1]
    auc_below     <- .trapz_below_threshold(times, vals, threshold = normal_thresh)

    data.frame(
      lineage          = lin,
      baseline         = round(baseline,  2),
      nadir_value      = round(nadir_val,  2),
      nadir_time_h     = round(nadir_time_h, 1),
      nadir_time_d     = round(nadir_time_h / 24, 1),
      grade_nadir      = grade,
      recovery_thresh  = round(recov_thresh, 2),
      time_to_recovery_h = round(t_recovery_h, 1),
      time_to_recovery_d = round(t_recovery_h / 24, 1),
      AUC_below_normal = round(auc_below, 2),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), results))
}


#' Intégration trapézoïdale de (threshold - x) pour x < threshold
.trapz_below_threshold <- function(times, values, threshold) {
  deficit <- pmax(0, threshold - values)
  # Règle des trapèzes
  n   <- length(times)
  auc <- 0
  for (i in seq_len(n - 1)) {
    dt   <- times[i + 1] - times[i]
    auc  <- auc + (deficit[i] + deficit[i + 1]) / 2 * dt
  }
  auc
}

# =============================================================================
# SECTION 4 — ANALYSE POPULATIONNELLE
# =============================================================================

#' Calcule les scores de toxicité pour toute une population
#'
#' @param pop_result  Résultat de simulate_population()
#' @param lineages    Lignées à analyser
#' @param ...         Autres arguments passés à score_patient()
#'
#' @return data.frame long avec une ligne par (patient, lignée)
score_population <- function(pop_result,
                             lineages = c("Neut", "Plt", "Ret", "RBC"),
                             ...) {

  sims <- pop_result$simulations
  n    <- length(sims)

  all_scores <- lapply(seq_len(n), function(i) {
    sc <- score_patient(sims[[i]], lineages = lineages, ...)
    if (!is.null(sc)) sc$patient_id <- i
    sc
  })

  df <- do.call(rbind, Filter(Negate(is.null), all_scores))
  df
}


#' Distribue les grades sur la population (fréquence par grade et par lignée)
#'
#' @param scores_df  data.frame issu de score_population()
#'
#' @return data.frame : lineage, grade, n, pct
grade_distribution <- function(scores_df) {

  lineages <- unique(scores_df$lineage)

  dist <- lapply(lineages, function(lin) {
    sub   <- scores_df[scores_df$lineage == lin, ]
    total <- nrow(sub)
    tab   <- table(sub$grade_nadir)
    grades <- 0:4
    counts <- sapply(as.character(grades), function(g) {
      ifelse(g %in% names(tab), tab[[g]], 0L)
    })
    data.frame(
      lineage = lin,
      grade   = grades,
      n       = as.integer(counts),
      pct     = round(as.numeric(counts) / total * 100, 1),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, dist)
}


#' Résumé statistique des nadirs pour une population
#'
#' @param scores_df  data.frame issu de score_population()
#'
#' @return data.frame : lineage + statistiques descriptives du nadir
nadir_summary <- function(scores_df) {

  lineages <- unique(scores_df$lineage)

  summ <- lapply(lineages, function(lin) {
    sub <- scores_df[scores_df$lineage == lin, ]
    q   <- quantile(sub$nadir_value, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    data.frame(
      lineage     = lin,
      n           = nrow(sub),
      mean_nadir  = round(mean(sub$nadir_value, na.rm = TRUE), 2),
      sd_nadir    = round(sd(sub$nadir_value,   na.rm = TRUE), 2),
      p05         = round(q["5%"],   2),
      p25         = round(q["25%"],  2),
      p50         = round(q["50%"],  2),
      p75         = round(q["75%"],  2),
      p95         = round(q["95%"],  2),
      pct_grade3_4 = round(mean(sub$grade_nadir >= 3, na.rm = TRUE) * 100, 1),
      pct_grade4   = round(mean(sub$grade_nadir == 4, na.rm = TRUE) * 100, 1),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summ)
}
