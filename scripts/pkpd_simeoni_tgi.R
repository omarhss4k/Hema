# =============================================================================
# Modele PK/PD TGI — Simeoni (2004) — Ajustement sequentiel
# Script generique — xenogreffe souris, dose unique IV bolus
#
# PK  : 2 compartiments IV bolus lineaire (parametres fixes issus du CSV)
#   dA1/dt = -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
#   dA2/dt =  (Q/V1)*A1 - (Q/V2)*A2
#   C1 = A1/V1
#
# PD  : Modele de Simeoni (compartiments de transit)
#   dx1/dt = [g(w)/w - k2*C1] * x1
#   dx2/dt =  k2*C1*x1 - k1*x2
#   dx3/dt =  k1*x2 - k1*x3
#   dx4/dt =  k1*x3 - k1*x4
#   g(w) = L0*w / (1 + (L0*w/L1)^psi)^(1/psi)   psi = 20
#
#   k2 [L/ug/j] : taux de cytotoxicite — estime directement en Etape B
#   TSC = L0/k2 [ug/L] : concentration seuil (derive, pour interpretation)
#   MTT = 4/k1  (jours)
#
# Ajustement SEQUENTIEL :
#   Etape A — Croissance  : L0 et L1 sur le groupe vehicule uniquement
#   Etape B — Efficacite  : k1 et k2 sur les groupes traites (L0/L1 fixes)
# =============================================================================

library(deSolve)
library(DEoptim)
library(ggplot2)
library(readxl)

# =============================================================================
# CONFIG — A MODIFIER POUR CHAQUE NOUVEAU PRODUIT / NOUVELLE COHORTE
# =============================================================================

# --- Fichiers de donnees ---
FICHIER_PK    <- "pk_resultats_20260720_150418.csv"   # CSV issu de l'analyse PK
FICHIER_TUMOR <- "Tumor_volum_SNU.xlsx"               # Excel volume tumoral

# --- Identite du produit ---
NOM_PRODUIT <- "Fc-silent FGFR2-huBPA-LP1"
NOM_COHORTE <- "SNU"

# --- Noms des groupes dans le CSV PK (colonne Animal) ---
#   Ces noms servent a extraire les Cmax par dose depuis le fichier PK.
NOM_DOSE_1 <- "groupe 4"   # dose la plus faible
NOM_DOSE_2 <- "grp 5"      # dose intermediaire
NOM_DOSE_3 <- "grp 7"      # dose la plus elevee

# --- Doses en ug/kg (meme ordre que NOM_DOSE_1/2/3) ---
DOSES_UGKG <- c(1200, 5600, 11800)   # [ug/kg]

# --- Bornes MTT pour l'Etape B ---
MTT_MIN_J <- 5    # [jours] — MTT minimum biologique plausible
MTT_MAX_J <- 25   # [jours] — MTT maximum biologique plausible

# --- Bornes k2 pour l'Etape B ---
#   k2 est estime directement (TSC = L0/k2 est derive pour interpretation).
#
#   k2_lower : drug tres faible — TSC jusqu'a K2_TSC_MAX_X * Cmax_dose3
#     Augmenter K2_TSC_MAX_X si le drug est tres peu potent.
#   k2_upper : drug tres puissant — TSC aussi bas que Cmax_dose1 / K2_TSC_MIN_DIV
#     Augmenter K2_TSC_MIN_DIV si le modele sur-supprime a faible dose.
K2_TSC_MAX_X   <- 1    # TSC_max = K2_TSC_MAX_X   * Cmax_dose3  (borne basse de k2)
K2_TSC_MIN_DIV <- 50   # TSC_min = Cmax_dose1 / K2_TSC_MIN_DIV (borne haute de k2)

# --- Poids de groupe pour la dose la plus faible (Etape B) ---
#   W_DOSE_FAIBLE = 1.0 : poids normal
#   W_DOSE_FAIBLE < 1.0 : reduit l'influence de ce groupe sur l'estimation de TSC.
#   Utile quand ce groupe a un comportement atypique (croissance initiale + stabilisation)
#   non reproductible avec un seul k2 = L0/TSC commun.
W_DOSE_FAIBLE <- 0.1

# --- Position des groupes dans le fichier Excel ---
#   Deux formats sont detectes automatiquement :
#
#   Format TALL (temps en col 1, groupes en colonnes) :
#     Col 1=Temps | Col 2=Ctrl | Col 3=Dose1 | Col 4=Dose2 | Col 5=Dose3
#     LIGNE_* = index de colonne de groupe (1=col2, 2=col3, 3=col4, 4=col5)
#
#   Format WIDE (groupes en lignes, temps en colonnes) :
#     Row 1 = header "Group" | t0 | t1 | ...
#     Row 2+ = groupe | val | val | ...
#     LIGNE_* = offset de ligne depuis l'en-tete (1=premiere ligne de donnees)
#
#   Dans les deux cas : 1=vehicule/ctrl, 2=dose1, 3=dose2, 4=dose3 par defaut.
LIGNE_CTRL  <- 1   # vehicule
LIGNE_DOSE1 <- 2   # dose la plus faible
LIGNE_DOSE2 <- 3   # dose intermediaire
LIGNE_DOSE3 <- 4   # dose la plus elevee

# =============================================================================
# FIN CONFIG — NE PAS MODIFIER CE QUI SUIT
# =============================================================================

# Titre du graphique (construit automatiquement depuis CONFIG)
TITRE_GRAPHE <- paste0("PK/PD TGI - Simeoni (2004) - ", NOM_PRODUIT, " - ", NOM_COHORTE)

# Etiquettes des groupes pour le graphe (construites depuis CONFIG)
LABEL_DOSE_1 <- paste0(DOSES_UGKG[1] / 1000, " mg/kg")
LABEL_DOSE_2 <- paste0(DOSES_UGKG[2] / 1000, " mg/kg")
LABEL_DOSE_3 <- paste0(DOSES_UGKG[3] / 1000, " mg/kg")

# =============================================================================
# 1. CHARGEMENT DES PARAMETRES PK
# =============================================================================

pk_data <- read.csv(FICHIER_PK, stringsAsFactors = FALSE)

# x1000 : mL/kg -> L/kg    x24 : /h -> /j
pk_fixed <- c(
  CL = (pk_data$CL[1] * 1000) * 24,
  V1 =  pk_data$V1[1] * 1000,
  V2 =  pk_data$V2[1] * 1000,
  Q  = (pk_data$Q[1]  * 1000) * 24
)

cat("=== Parametres PK fixes (souris, unites journalieres) ===\n")
cat(sprintf("CL = %.5e L/j/kg\n", pk_fixed["CL"]))
cat(sprintf("V1 = %.5e L/kg\n",   pk_fixed["V1"]))
cat(sprintf("V2 = %.5e L/kg\n",   pk_fixed["V2"]))
cat(sprintf("Q  = %.5e L/j/kg\n", pk_fixed["Q"]))

# Cmax (ng/mL = ug/L) issues du CSV
Cmax_dose1 <- pk_data$cmax[pk_data$Animal == NOM_DOSE_1]
Cmax_dose2 <- pk_data$cmax[pk_data$Animal == NOM_DOSE_2]
Cmax_dose3 <- pk_data$cmax[pk_data$Animal == NOM_DOSE_3]

cat(sprintf("\nCmax (ug/L) : %s = %.0f | %s = %.0f | %s = %.0f\n",
            LABEL_DOSE_1, Cmax_dose1,
            LABEL_DOSE_2, Cmax_dose2,
            LABEL_DOSE_3, Cmax_dose3))

# =============================================================================
# 2. DONNEES TUMORALES
#
# Deux formats detectes automatiquement (voir CONFIG LIGNE_*) :
#   TALL : temps en col 1, groupes en colonnes (format habituel)
#   WIDE : temps en ligne d'en-tete, groupes en lignes
# =============================================================================

raw_tumor <- suppressMessages(
  read_xlsx(FICHIER_TUMOR, col_names = FALSE, .name_repair = "minimal")
)

# Localiser les deux lignes d'en-tete (contiennent "Time" ou "Group")
header_rows <- which(grepl("Time|Group",
                           trimws(as.character(raw_tumor[[1]])),
                           ignore.case = TRUE) &
                     !is.na(raw_tumor[[1]]))
stopifnot("2 blocs attendus (moyennes + SEM)" = length(header_rows) == 2)

fix_sem <- function(x) ifelse(is.na(x) | x <= 0, 1, x)

# Detecter le format : col2 de la ligne d'en-tete est numerique → WIDE
is_wide <- !is.na(suppressWarnings(as.numeric(as.character(
  raw_tumor[[2]][header_rows[1]]
))))

if (is_wide) {
  # Format WIDE — temps dans la ligne d'en-tete, groupes dans les lignes suivantes
  times_d <- suppressWarnings(as.numeric(as.character(raw_tumor[header_rows[1], ])))
  times_d <- times_d[!is.na(times_d)]
  n_t     <- length(times_d)

  read_row_w <- function(h_row, offset) {
    suppressWarnings(as.numeric(as.character(raw_tumor[h_row + offset, 2:(n_t + 1)])))
  }
  tv_ctrl  <- read_row_w(header_rows[1], LIGNE_CTRL)
  tv_dose1 <- read_row_w(header_rows[1], LIGNE_DOSE1)
  tv_dose2 <- read_row_w(header_rows[1], LIGNE_DOSE2)
  tv_dose3 <- read_row_w(header_rows[1], LIGNE_DOSE3)
  sem_ctrl  <- fix_sem(read_row_w(header_rows[2], LIGNE_CTRL))
  sem_dose1 <- fix_sem(read_row_w(header_rows[2], LIGNE_DOSE1))
  sem_dose2 <- fix_sem(read_row_w(header_rows[2], LIGNE_DOSE2))
  sem_dose3 <- fix_sem(read_row_w(header_rows[2], LIGNE_DOSE3))

} else {
  # Format TALL — temps en col 1, groupes en colonnes
  extract_block <- function(raw, r1, r2) {
    blk <- raw[seq(r1, r2), ]
    blk[!is.na(suppressWarnings(as.numeric(as.character(blk[[1]])))), ]
  }
  b1 <- extract_block(raw_tumor, header_rows[1] + 1, header_rows[2] - 1)
  b2 <- extract_block(raw_tumor, header_rows[2] + 1, nrow(raw_tumor))

  num_col <- function(df, j) suppressWarnings(as.numeric(as.character(df[[j]])))

  times_d  <- num_col(b1, 1)
  tv_ctrl  <- num_col(b1, LIGNE_CTRL  + 1)
  tv_dose1 <- num_col(b1, LIGNE_DOSE1 + 1)
  tv_dose2 <- num_col(b1, LIGNE_DOSE2 + 1)
  tv_dose3 <- num_col(b1, LIGNE_DOSE3 + 1)
  sem_ctrl  <- fix_sem(num_col(b2, LIGNE_CTRL  + 1))
  sem_dose1 <- fix_sem(num_col(b2, LIGNE_DOSE1 + 1))
  sem_dose2 <- fix_sem(num_col(b2, LIGNE_DOSE2 + 1))
  sem_dose3 <- fix_sem(num_col(b2, LIGNE_DOSE3 + 1))
}

get_tv0 <- function(tv) {
  v <- tv[times_d == 0]
  if (length(v) == 0 || is.na(v[1])) tv[!is.na(tv)][1] else v[1]
}

tv0_ctrl  <- get_tv0(tv_ctrl)
tv0_dose1 <- get_tv0(tv_dose1)
tv0_dose2 <- get_tv0(tv_dose2)
tv0_dose3 <- get_tv0(tv_dose3)

ok_ctrl  <- !is.na(tv_ctrl)
ok_dose1 <- !is.na(tv_dose1)
ok_dose2 <- !is.na(tv_dose2)
ok_dose3 <- !is.na(tv_dose3)

cat(sprintf("\nTV0 (j0) : Ctrl=%.0f | %s=%.0f | %s=%.0f | %s=%.0f mm3\n",
            tv0_ctrl,
            LABEL_DOSE_1, tv0_dose1,
            LABEL_DOSE_2, tv0_dose2,
            LABEL_DOSE_3, tv0_dose3))

# =============================================================================
# 3. EQUATIONS DU MODELE DE SIMEONI  <- NE PAS MODIFIER
# =============================================================================

PSI <- 20

simeoni_rhs <- function(t, state, parms) {
  A1 <- max(state["A1"], 0); A2 <- max(state["A2"], 0)
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)

  CL <- parms["CL"]; V1 <- parms["V1"]
  V2 <- parms["V2"]; Q  <- parms["Q"]
  L0 <- parms["L0"]; L1 <- parms["L1"]
  k1 <- parms["k1"]; k2 <- parms["k2"]

  C1  <- A1 / V1
  dA1 <- -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
  dA2 <-  (Q/V1)*A1 - (Q/V2)*A2

  w  <- x1 + x2 + x3 + x4
  gw <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  dx1 <- (growth_rate - k2 * C1) * x1
  dx2 <- k2 * C1 * x1 - k1 * x2
  dx3 <- k1 * x2 - k1 * x3
  dx4 <- k1 * x3 - k1 * x4

  list(c(A1 = dA1, A2 = dA2, x1 = dx1, x2 = dx2, x3 = dx3, x4 = dx4))
}

simeoni_ctrl_rhs <- function(t, state, parms) {
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)
  L0 <- parms["L0"]; L1 <- parms["L1"]

  w  <- x1 + x2 + x3 + x4
  gw <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  list(c(x1 = growth_rate * x1, x2 = 0, x3 = 0, x4 = 0))
}

# =============================================================================
# 4. FONCTIONS DE SIMULATION  <- NE PAS MODIFIER
# =============================================================================

sim_ctrl_fn <- function(L0, L1, tv0, times_out) {
  t_all <- sort(unique(c(0, times_out)))
  out <- tryCatch(
    as.data.frame(lsoda(
      y     = c(x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
      times = t_all,
      func  = simeoni_ctrl_rhs,
      parms = c(L0 = L0, L1 = L1)
    )),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(times_out)))
  w_tot <- out$x1 + out$x2 + out$x3 + out$x4
  approx(out$time, w_tot, xout = times_out, rule = 2)$y
}

sim_treated <- function(dose_ugkg, tv0, params, times_out) {
  t_all  <- sort(unique(c(0, times_out)))
  params <- setNames(as.numeric(params),
                     c("CL", "V1", "V2", "Q", "L0", "L1", "k1", "k2"))
  out <- tryCatch(
    as.data.frame(lsoda(
      y        = c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
      times    = t_all,
      func     = simeoni_rhs,
      parms    = params,
      atol     = 1e-6, rtol = 1e-6,
      maxsteps = 50000
    )),
    error   = function(e) { message("lsoda error: ", conditionMessage(e)); NULL },
    warning = function(w) {
      suppressWarnings(
        as.data.frame(lsoda(
          y        = c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
          times    = t_all,
          func     = simeoni_rhs,
          parms    = params,
          atol     = 1e-4, rtol = 1e-4,
          maxsteps = 100000
        ))
      )
    }
  )
  if (is.null(out) || any(is.na(out$x1))) return(rep(NA_real_, length(times_out)))
  w_tot <- pmax(out$x1, 0) + pmax(out$x2, 0) + pmax(out$x3, 0) + pmax(out$x4, 0)
  approx(out$time, w_tot, xout = times_out, rule = 2)$y
}

# =============================================================================
# 5A. ETAPE A — Croissance : L0 et L1 sur le groupe vehicule uniquement
#     Bornes L0 : doublement dans [3, 14 j]
#     Bornes L1 : [TV_max observe, 200 x TV_max]
# =============================================================================

cat("\n===========================================================\n")
cat("ETAPE A — Ajustement vehicule : L0 et L1\n")
cat("===========================================================\n")

lm_ctrl <- lm(log(tv_ctrl[ok_ctrl]) ~ times_d[ok_ctrl])
L0_init <- max(coef(lm_ctrl)[2], 0.005)

lower_A <- c(log(0.05),
             log(max(tv_ctrl[ok_ctrl], na.rm = TRUE)))
upper_A <- c(log(0.25),
             log(max(tv_ctrl[ok_ctrl], na.rm = TRUE) * 200))

objective_vehicule <- function(logpar) {
  L0 <- exp(logpar[1])
  L1 <- exp(logpar[2])
  if (L0 <= 0 || L1 <= 0) return(1e12)

  pred <- sim_ctrl_fn(L0, L1, tv0_ctrl, times_d[ok_ctrl])
  if (any(is.na(pred))) return(1e12)
  pred <- pmax(pred, 0.1)

  sum(1 / sem_ctrl[ok_ctrl]^2 * (log(tv_ctrl[ok_ctrl]) - log(pred))^2,
      na.rm = TRUE)
}

set.seed(42)
fit_de_A <- DEoptim(
  fn      = objective_vehicule,
  lower   = lower_A,
  upper   = upper_A,
  control = DEoptim.control(
    NP      = 80,
    itermax = 600,
    F       = 0.8,
    CR      = 0.9,
    trace   = 100,
    reltol  = 1e-8,
    steptol = 150
  )
)

fit_local_A <- nlminb(
  start     = fit_de_A$optim$bestmem,
  objective = objective_vehicule,
  lower     = lower_A,
  upper     = upper_A,
  control   = list(eval.max = 2000, iter.max = 800,
                   rel.tol  = 1e-12, x.tol = 1e-12)
)

if (fit_local_A$objective < fit_de_A$optim$bestval) {
  best_par_A <- fit_local_A$par
  best_val_A <- fit_local_A$objective
  src_A      <- "nlminb"
} else {
  best_par_A <- as.numeric(fit_de_A$optim$bestmem)
  best_val_A <- fit_de_A$optim$bestval
  src_A      <- "DEoptim"
}

L0_est <- exp(best_par_A[1])
L1_est <- exp(best_par_A[2])

cat(sprintf("\n-> L0 = %.6f /j   [doublement = %.1f j]  (%s)\n",
            L0_est, log(2) / L0_est, src_A))
cat(sprintf("-> L1 = %.2f mm3/j  (%s)\n", L1_est, src_A))
cat(sprintf("   Objectif DEoptim : %.8f\n", fit_de_A$optim$bestval))
cat(sprintf("   Objectif retenu  : %.8f\n", best_val_A))

# =============================================================================
# 5B. ETAPE B — Efficacite : k1 et TSC (L0/L1 fixes depuis Etape A)
#
#   k2 derive : k2 = L0_est / TSC
#   Bornes TSC : [Cmax_dose1 / FACTEUR_TSC_LOWER, Cmax_dose2 / FACTEUR_TSC_UPPER]
#   Bornes k1  : MTT dans [MTT_MIN_J, MTT_MAX_J]
#   Poids dose la plus faible : W_DOSE_FAIBLE
# =============================================================================

cat("\n===========================================================\n")
cat("ETAPE B — Ajustement groupes traites : k1 et k2\n")
cat("===========================================================\n")
cat(sprintf("L0 fixe = %.6f /j  |  L1 fixe = %.2f mm3/j\n", L0_est, L1_est))

# Bornes k2 calculees depuis CONFIG et Cmax observes
TSC_min   <- Cmax_dose1 / K2_TSC_MIN_DIV          # drug le plus puissant
TSC_max   <- Cmax_dose3 * K2_TSC_MAX_X             # drug le plus faible
k2_lower  <- L0_est / TSC_max
k2_upper  <- L0_est / TSC_min
k1_lower  <- 4 / MTT_MAX_J
k1_upper  <- 4 / MTT_MIN_J

cat(sprintf("Poids %s : %.2f\n", LABEL_DOSE_1, W_DOSE_FAIBLE))
cat(sprintf("Bornes TSC equivalentes : [%.0f, %.0f] ug/L\n", TSC_min, TSC_max))
cat(sprintf("Bornes k2 : [%.3e, %.3e] L/ug/j\n", k2_lower, k2_upper))
cat(sprintf("Bornes k1 : [%.4f, %.4f] /j    (MTT dans [%d, %d] j)\n",
            k1_lower, k1_upper, MTT_MIN_J, MTT_MAX_J))

lower_B <- c(log(k1_lower), log(k2_lower))
upper_B <- c(log(k1_upper), log(k2_upper))

objective_traites <- function(logpar) {
  k1 <- exp(logpar[1])
  k2 <- exp(logpar[2])
  if (k1 <= 0 || k2 <= 0) return(1e12)

  params_all <- c(
    CL = unname(pk_fixed["CL"]), V1 = unname(pk_fixed["V1"]),
    V2 = unname(pk_fixed["V2"]), Q  = unname(pk_fixed["Q"]),
    L0 = L0_est, L1 = L1_est, k1 = k1, k2 = k2
  )

  pred1 <- sim_treated(DOSES_UGKG[1], tv0_dose1, params_all, times_d[ok_dose1])
  if (any(is.na(pred1))) return(1e12)
  pred1 <- pmax(pred1, 0.1)

  pred2 <- sim_treated(DOSES_UGKG[2], tv0_dose2, params_all, times_d[ok_dose2])
  if (any(is.na(pred2))) return(1e12)
  pred2 <- pmax(pred2, 0.1)

  pred3 <- sim_treated(DOSES_UGKG[3], tv0_dose3, params_all, times_d[ok_dose3])
  if (any(is.na(pred3))) return(1e12)
  pred3 <- pmax(pred3, 0.1)

  W_DOSE_FAIBLE * sum(1 / sem_dose1[ok_dose1]^2 * (log(tv_dose1[ok_dose1]) - log(pred1))^2, na.rm = TRUE) +
                  sum(1 / sem_dose2[ok_dose2]^2 * (log(tv_dose2[ok_dose2]) - log(pred2))^2, na.rm = TRUE) +
                  sum(1 / sem_dose3[ok_dose3]^2 * (log(tv_dose3[ok_dose3]) - log(pred3))^2, na.rm = TRUE)
}

# Verification avant DEoptim
params_sanity <- c(
  CL = unname(pk_fixed["CL"]), V1 = unname(pk_fixed["V1"]),
  V2 = unname(pk_fixed["V2"]), Q  = unname(pk_fixed["Q"]),
  L0 = L0_est, L1 = L1_est,
  k1 = 4 / 15,
  k2 = L0_est / 3000
)
test_ok <- sim_treated(DOSES_UGKG[1], tv0_dose1, params_sanity, times_d[ok_dose1])
if (all(is.na(test_ok))) {
  stop("ERREUR : sim_treated retourne NA avec des parametres raisonnables.")
}
cat(sprintf("  [Sanity check] sim_treated %s : OK  (pred[1] = %.1f mm3)\n",
            LABEL_DOSE_1, test_ok[1]))

set.seed(42)
fit_de_B <- DEoptim(
  fn      = objective_traites,
  lower   = lower_B,
  upper   = upper_B,
  control = DEoptim.control(
    NP      = 80,
    itermax = 800,
    F       = 0.8,
    CR      = 0.9,
    trace   = 100,
    reltol  = 1e-8,
    steptol = 200
  )
)

fit_local_B <- nlminb(
  start     = fit_de_B$optim$bestmem,
  objective = objective_traites,
  lower     = lower_B,
  upper     = upper_B,
  control   = list(eval.max = 3000, iter.max = 1000,
                   rel.tol  = 1e-12, x.tol = 1e-12)
)

if (fit_local_B$objective < fit_de_B$optim$bestval) {
  best_par_B <- fit_local_B$par
  best_val_B <- fit_local_B$objective
  src_B      <- "nlminb"
} else {
  best_par_B <- as.numeric(fit_de_B$optim$bestmem)
  best_val_B <- fit_de_B$optim$bestval
  src_B      <- "DEoptim"
}

k1_est  <- exp(best_par_B[1])
k2_est  <- exp(best_par_B[2])
TSC_est <- L0_est / k2_est

cat(sprintf("\n-> k1  = %.5f /j   [MTT = %.1f j]  (%s)\n", k1_est, 4 / k1_est, src_B))
cat(sprintf("-> TSC = %.0f ug/L  (%s)\n", TSC_est, src_B))
cat(sprintf("   k2 derive = %.3e L/ug/j\n", k2_est))
cat(sprintf("   Objectif DEoptim : %.8f\n", fit_de_B$optim$bestval))
cat(sprintf("   Objectif retenu  : %.8f\n", best_val_B))

# Avertissements si parametre sur borne
tol <- 0.02
rng_k2 <- log(k2_upper) - log(k2_lower)
rng_k1 <- log(k1_upper) - log(k1_lower)
if (abs(log(k2_est) - log(k2_lower)) < tol * rng_k2)
  cat(sprintf("  [AVERT] k2 sur borne inferieure (TSC=%.0f ug/L — drug tres faible)\n", TSC_est))
if (abs(log(k2_est) - log(k2_upper))  < tol * rng_k2)
  cat(sprintf("  [AVERT] k2 sur borne superieure (TSC=%.0f ug/L — drug tres puissant)\n", TSC_est))
if (abs(log(k1_est) - log(k1_upper))  < tol * rng_k1)
  cat(sprintf("  [AVERT] k1 sur borne superieure (MTT=%.1fj)\n", 4/k1_est))
if (abs(log(k1_est) - log(k1_lower))  < tol * rng_k1)
  cat(sprintf("  [AVERT] k1 sur borne inferieure (MTT=%.1fj)\n", 4/k1_est))

# =============================================================================
# 6. RESUME FINAL
# =============================================================================

mtt <- 4 / k1_est

cat("\n===========================================================\n")
cat(sprintf("PARAMETRES FINAUX — %s  (ajustement sequentiel)\n", NOM_COHORTE))
cat("===========================================================\n")
cat(sprintf("L0  = %.6f /j         [Etape A — %s]\n", L0_est,  src_A))
cat(sprintf("L1  = %.2f mm3/j      [Etape A — %s]\n", L1_est,  src_A))
cat(sprintf("k1  = %.5f /j         [Etape B — %s]\n", k1_est,  src_B))
cat(sprintf("TSC = %.0f ug/L       [Etape B — %s]\n", TSC_est, src_B))
cat("-----------------------------------------------------------\n")
cat(sprintf("MTT = %.1f j\n", mtt))
cat(sprintf("k2  = %.3e L/ug/j  (= L0/TSC)\n", k2_est))
cat(sprintf("Plausibilite : Cmax(%s)=%.0f | TSC=%.0f | Cmax(%s)=%.0f ug/L\n",
            LABEL_DOSE_1, Cmax_dose1, TSC_est, LABEL_DOSE_2, Cmax_dose2))
cat(sprintf("  Cmax/TSC @ %s : %.1f\n", LABEL_DOSE_1, Cmax_dose1 / TSC_est))
cat(sprintf("  Cmax/TSC @ %s : %.1f\n", LABEL_DOSE_2, Cmax_dose2 / TSC_est))
cat(sprintf("  Cmax/TSC @ %s : %.1f\n", LABEL_DOSE_3, Cmax_dose3 / TSC_est))

# =============================================================================
# 7. SIMULATION FINALE
# =============================================================================

params_best <- c(
  CL = unname(pk_fixed["CL"]), V1 = unname(pk_fixed["V1"]),
  V2 = unname(pk_fixed["V2"]), Q  = unname(pk_fixed["Q"]),
  L0 = L0_est, L1 = L1_est, k1 = k1_est, k2 = k2_est
)

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

pred_ctrl_sim  <- sim_ctrl_fn(L0_est, L1_est, tv0_ctrl,  times_sim)
pred_dose1_sim <- sim_treated(DOSES_UGKG[1], tv0_dose1, params_best, times_sim)
pred_dose2_sim <- sim_treated(DOSES_UGKG[2], tv0_dose2, params_best, times_sim)
pred_dose3_sim <- sim_treated(DOSES_UGKG[3], tv0_dose3, params_best, times_sim)

niv <- c("Vehicule", LABEL_DOSE_1, LABEL_DOSE_2, LABEL_DOSE_3)

df_sim <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_sim,  Groupe = "Vehicule"),
  data.frame(jour = times_sim, TV = pred_dose1_sim, Groupe = LABEL_DOSE_1),
  data.frame(jour = times_sim, TV = pred_dose2_sim, Groupe = LABEL_DOSE_2),
  data.frame(jour = times_sim, TV = pred_dose3_sim, Groupe = LABEL_DOSE_3)
)
df_sim$Groupe <- factor(df_sim$Groupe, levels = niv)

df_obs <- rbind(
  data.frame(jour = times_d[ok_ctrl],  TV = tv_ctrl[ok_ctrl],    sem = sem_ctrl[ok_ctrl],   Groupe = "Vehicule"),
  data.frame(jour = times_d[ok_dose1], TV = tv_dose1[ok_dose1],  sem = sem_dose1[ok_dose1], Groupe = LABEL_DOSE_1),
  data.frame(jour = times_d[ok_dose2], TV = tv_dose2[ok_dose2],  sem = sem_dose2[ok_dose2], Groupe = LABEL_DOSE_2),
  data.frame(jour = times_d[ok_dose3], TV = tv_dose3[ok_dose3],  sem = sem_dose3[ok_dose3], Groupe = LABEL_DOSE_3)
)
df_obs$Groupe <- factor(df_obs$Groupe, levels = niv)

# =============================================================================
# 8. GRAPHIQUE
# =============================================================================

cols <- c("Vehicule" = "#888888")
palette_doses <- c("#74ADD1", "#4393C3", "#2166AC")
cols[LABEL_DOSE_1] <- palette_doses[1]
cols[LABEL_DOSE_2] <- palette_doses[2]
cols[LABEL_DOSE_3] <- palette_doses[3]

ymax <- max(df_obs$TV + df_obs$sem, na.rm = TRUE)

subtitle_txt <- paste0(
  "A: L0=", round(L0_est, 4), "/j | L1=", format(round(L1_est), big.mark = " "), "mm3/j",
  "  ||  ",
  "B: k1=", round(k1_est, 3), "/j | TSC=", round(TSC_est, 0), "ug/L | MTT=", round(mtt, 1), "j"
)

p_simeoni <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(data    = df_sim[df_sim$Groupe != "Vehicule", ],
            mapping = aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1) +
  geom_line(data    = df_sim[df_sim$Groupe == "Vehicule", ],
            mapping = aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1, linetype = "dashed") +
  geom_errorbar(data    = df_obs,
                mapping = aes(x = jour, ymin = TV - sem, ymax = TV + sem, color = Groupe),
                width = 0.8, linewidth = 0.5) +
  geom_point(data    = df_obs,
             mapping = aes(x = jour, y = TV, color = Groupe),
             size = 2.5) +
  annotate("point", x = 0,   y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text",  x = 1.5, y = -ymax * 0.06,
           label = "= Dose unique (j0)", hjust = 0, size = 3.2, color = "#CC0000") +
  scale_x_continuous(breaks = seq(0, max(times_d, na.rm = TRUE), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = TITRE_GRAPHE,
    subtitle = subtitle_txt,
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm3)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 9, color = "grey50"),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

print(p_simeoni)

nom_base <- paste0("plot_PKPD_simeoni_", NOM_COHORTE)
ggsave(paste0("scripts/", nom_base, ".png"), p_simeoni,
       width = 9, height = 5.5, dpi = 150)
cat(sprintf("\nGraphique sauvegarde : scripts/%s.png\n", nom_base))

# =============================================================================
# 9. SAUVEGARDE
# =============================================================================

simeoni_results <- list(
  produit       = NOM_PRODUIT,
  cohorte       = NOM_COHORTE,
  pk_fixed      = pk_fixed,
  L0            = L0_est,
  L1            = L1_est,
  k1            = k1_est,
  TSC_ugL       = TSC_est,
  k2            = k2_est,
  MTT_days      = mtt,
  obj_A_DEoptim = fit_de_A$optim$bestval,
  obj_A_final   = best_val_A,
  obj_B_DEoptim = fit_de_B$optim$bestval,
  obj_B_final   = best_val_B
)

nom_rdata <- paste0("scripts/resultats_PKPD_simeoni_", NOM_COHORTE, ".RData")
save(simeoni_results, file = nom_rdata)
cat(sprintf("Resultats sauvegardes : %s\n", nom_rdata))
