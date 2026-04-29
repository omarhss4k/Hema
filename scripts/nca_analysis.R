# =============================================================================
# Analyse PK — modèle 2 compartiments IV bolus avec rxode2
# =============================================================================
# Modèle :   d/dt(A1) = -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
#            d/dt(A2) =   (Q/V1)*A1 - (Q/V2)*A2
#            C = A1 / V1
# Unités  :  A1, A2 [ng/kg]  |  V1, V2 [mL/kg]  |  CL, Q [mL/h/kg]  |  C [ng/mL]
# =============================================================================

# install.packages(c("rxode2", "dplyr", "tidyr", "ggplot2"))
library(rxode2)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. Données ----------------------------------------------------------------
# BLQ terminaux (336h, 504h) → NA  |  96h et 168h sont quantifiables

pk_raw <- bind_rows(

  data.frame(
    subject    = "Animal_01",
    dose_mg_kg = 3,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     69800, 48500, 17700, 6730, 3550, 2200, 1060, NA,  NA )
  ),

  data.frame(
    subject    = "Animal_02",
    dose_mg_kg = 3,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     73200, 50200, 16600, 6370, 3680, 2180, 1020, NA,  NA )
  )
)

# ── 2. Définition du modèle rxode2 --------------------------------------------
pk2cmt <- rxode2({
  d/dt(A1) <- -(CL / V1 + Q / V1) * A1 + (Q / V2) * A2
  d/dt(A2) <-  (Q / V1) * A1 - (Q / V2) * A2
  C        <- A1 / V1
})

# ── 3. Ajustement par animal (moindres carrés sur log-concentrations) ---------
# La fonction objectif minimise Σ [ln(Cobs) - ln(Cpred)]²
# Les paramètres sont optimisés sur l'échelle log (→ contrainte de positivité).

fit_animal <- function(animal_data) {

  dose  <- animal_data$dose_mg_kg[1]
  A1_0  <- dose * 1e6        # mg/kg → ng/kg  (1 mg = 1e6 ng)
  obs   <- animal_data %>% filter(time > 0, !is.na(conc))

  objective <- function(log_p) {
    p <- setNames(exp(log_p), c("CL", "V1", "Q", "V2"))
    tryCatch({
      sol  <- rxSolve(pk2cmt,
                      params = p,
                      inits  = c(A1 = A1_0, A2 = 0),
                      events = et(obs$time))
      pred <- sol$C
      if (any(!is.finite(pred) | pred <= 0)) return(1e10)
      sum((log(obs$conc) - log(pred))^2)
    }, error = function(e) 1e10)
  }

  # Valeurs initiales — estimées à partir des données :
  #   V1 ≈ dose/Cmax ≈ 3e6/70000 ≈ 43 mL/kg
  #   CL ≈ 1.73 mL/h/kg (NCA préliminaire)
  #   Q  ≈ 5 mL/h/kg, V2 ≈ 130 mL/kg
  log_p0 <- log(c(CL = 1.73, V1 = 43, Q = 5, V2 = 130))

  fit <- optim(log_p0, objective,
               method  = "Nelder-Mead",
               control = list(maxit = 20000, reltol = 1e-12))

  exp(fit$par) |> setNames(c("CL", "V1", "Q", "V2"))
}

# ── 4. Paramètres dérivés (analytiques, modèle 2-cmt) -------------------------
# Les constantes microscopiques donnent les valeurs propres alpha (rapide)
# et beta (terminale) du système biexponentiel :
#   C(t) = A·exp(-alpha·t) + B·exp(-beta·t)
# La demi-vie terminale est t1/2 = ln(2) / beta.
# Le volume de distribution terminal : Vz = CL / beta.

derived_params <- function(p) {
  k10  <- p["CL"] / p["V1"]
  k12  <- p["Q"]  / p["V1"]
  k21  <- p["Q"]  / p["V2"]
  S    <- k10 + k12 + k21
  disc <- sqrt(S^2 - 4 * k10 * k21)
  alpha <- (S + disc) / 2
  beta  <- (S - disc) / 2
  list(
    half_life_alpha = log(2) / alpha,
    half_life_beta  = log(2) / beta,
    Vz              = p["CL"] / beta
  )
}

# ── 5. AUClast (trapèzes linéaires) ------------------------------------------
auclast_animal <- function(animal_data) {
  d <- animal_data %>% filter(!is.na(conc)) %>% arrange(time)
  n <- nrow(d)
  sum(diff(d$time) * (d$conc[-n] + d$conc[-1]) / 2)
}

# ── 6. Boucle sur les animaux -------------------------------------------------
animals  <- unique(pk_raw$subject)
results  <- list()

for (anim in animals) {
  cat("\n--- Ajustement :", anim, "---\n")
  dat   <- pk_raw %>% filter(subject == anim)
  dose  <- dat$dose_mg_kg[1]

  p     <- fit_animal(dat)
  drv   <- derived_params(p)
  aucl  <- auclast_animal(dat)

  # AUCinf = AUClast + C_last / beta  (C_last = dernière conc quantifiable)
  c_last <- dat %>% filter(!is.na(conc)) %>% slice_max(time, n=1) %>% pull(conc)
  beta   <- log(2) / drv$half_life_beta
  aucinf <- aucl + c_last / beta

  results[[anim]] <- data.frame(
    Dose_mg_kg      = dose,
    Animal_Id       = anim,
    Half_life_h     = round(drv$half_life_beta, 1),
    Cmax_ng_mL      = round(max(dat$conc, na.rm=TRUE), 0),
    Cmax_D          = round(max(dat$conc, na.rm=TRUE) / dose, 0),
    AUClast_h_ng_mL = round(aucl,   0),
    AUCinf_h_ng_mL  = round(aucinf, 0),
    Vz_mL_kg        = round(drv$Vz, 0),
    CL_mL_h_kg      = round(p["CL"], 2)
  )
}

nca_summary <- bind_rows(results)
row.names(nca_summary) <- NULL

units_row <- data.frame(
  Dose_mg_kg = "mg/kg", Animal_Id = "", Half_life_h = "h",
  Cmax_ng_mL = "ng/mL", Cmax_D = "ng/mL/mg/kg",
  AUClast_h_ng_mL = "h*ng/mL", AUCinf_h_ng_mL = "h*ng/mL",
  Vz_mL_kg = "mL/kg", CL_mL_h_kg = "mL/h/kg",
  stringsAsFactors = FALSE
)

cat("\n=== Paramètres PK — tableau de sortie ===\n")
print(rbind(units_row, nca_summary), row.names = FALSE)
write.csv(nca_summary, "nca_results.csv", row.names = FALSE)
cat("\nTableau exporté : nca_results.csv\n")

# ── 7. Courbes ajustées + données observées -----------------------------------
# Générer les prédictions du modèle ajusté pour chaque animal

pred_list <- list()
for (anim in animals) {
  dat  <- pk_raw %>% filter(subject == anim)
  dose <- dat$dose_mg_kg[1]
  p    <- fit_animal(dat)

  times_pred <- seq(0.083, 168, length.out = 300)
  sol <- rxSolve(pk2cmt,
                 params = p,
                 inits  = c(A1 = dose * 1e6, A2 = 0),
                 events = et(times_pred))

  pred_list[[anim]] <- data.frame(
    subject = anim, time = sol$time, conc = sol$C
  )
}
pred_df <- bind_rows(pred_list)

obs_df  <- pk_raw %>% filter(!is.na(conc), conc > 0)

theme_pk <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# 7a. Linéaire
p_linear <- ggplot() +
  geom_line(data = pred_df, aes(x=time, y=conc, color=subject), linewidth=0.9) +
  geom_point(data = obs_df, aes(x=time, y=conc, color=subject), size=3) +
  labs(title = "Profil PK — modèle 2 compartiments (rxode2)",
       subtitle = "Points = observations ; lignes = modèle ajusté",
       x = "Temps (h)", y = "Concentration (ng/mL)", color = "Animal") +
  theme_pk
print(p_linear)
ggsave("nca_linear.png", plot = p_linear, width = 8, height = 5, dpi = 300)

# 7b. Semi-logarithmique
p_semilog <- ggplot() +
  geom_line(data = pred_df, aes(x=time, y=conc, color=subject), linewidth=0.9) +
  geom_point(data = obs_df, aes(x=time, y=conc, color=subject), size=3) +
  scale_y_log10() +
  labs(title = "Profil PK — modèle 2 compartiments (rxode2) — échelle semi-log",
       subtitle = "Points = observations ; lignes = modèle ajusté",
       x = "Temps (h)", y = "Concentration (ng/mL) — log", color = "Animal") +
  theme_pk
print(p_semilog)
ggsave("nca_semilog.png", plot = p_semilog, width = 8, height = 5, dpi = 300)

cat("\nGraphiques exportés : nca_linear.png  nca_semilog.png\n")
cat("Script terminé.\n")
