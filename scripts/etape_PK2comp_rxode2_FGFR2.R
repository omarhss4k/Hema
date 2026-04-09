# =============================================================================
# PK 2-compartiments — rxode2 + nlminb
# Fc-silent B/C huBPA-LP1 (FGFR2)
# IV bolus, 3 doses : 1, 5, 10 mg/kg
# Valeurs initiales : résultats de la méthode des résidus
# =============================================================================

library(rxode2)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. DONNÉES
# =============================================================================

raw    <- read_xlsx("data/PK souris FGFR2.xlsx")
time_h <- as.numeric(raw[[1]])
conc10 <- as.numeric(raw[[2]])
conc5  <- as.numeric(raw[[3]])
conc1  <- as.numeric(raw[[4]])

dose10 <- 10 * 1000   # µg/kg
dose5  <-  5 * 1000
dose1  <-  1 * 1000

mk <- function(t, c, d, label) {
  ok <- !is.na(c) & c > 0
  data.frame(t = t[ok], C = c[ok], dose = d, Dose = label)
}
df <- rbind(
  mk(time_h, conc10, dose10, "10 mg/kg"),
  mk(time_h, conc5,  dose5,  "5 mg/kg"),
  mk(time_h, conc1,  dose1,  "1 mg/kg")
)
df$Dose <- factor(df$Dose, levels = c("10 mg/kg", "5 mg/kg", "1 mg/kg"))

d10 <- df[df$Dose == "10 mg/kg", ]
d5  <- df[df$Dose == "5 mg/kg",  ]
d1  <- df[df$Dose == "1 mg/kg",  ]

# =============================================================================
# 2. MODÈLE rxode2 — 2 compartiments, IV bolus
#
#   A1 = quantité dans le compartiment central  (µg/kg)
#   A2 = quantité dans le compartiment périphérique
#   C1 = A1 / V1  (concentration observée)
#
#   dA1/dt = -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
#   dA2/dt =  (Q/V1)*A1 - (Q/V2)*A2
# =============================================================================

mod2comp <- rxode2({
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
  d/dt(A2) <-  (Q/V1) * A1 - (Q/V2) * A2
})

# Simulation d'un groupe (dose unique IV bolus à t=0)
sim_dose <- function(dose, params, times) {
  ev <- eventTable()
  ev$add.dosing(dose = dose, nbr.doses = 1, dosing.to = 1)
  ev$add.sampling(sort(unique(c(0, times))))

  out <- tryCatch(
    rxSolve(mod2comp, params, ev),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(times)))
  approx(out$time, out$C1, xout = times, rule = 2)$y
}

# =============================================================================
# 3. VALEURS INITIALES — chargées depuis la méthode des résidus
# =============================================================================

load("scripts/resultats_PK2comp_FGFR2.RData")   # → pk2comp

init_params <- c(
  CL = pk2comp$CL,
  V1 = pk2comp$V1,
  V2 = pk2comp$V2,
  Q  = pk2comp$Q
)

cat("=== Valeurs initiales (méthode des résidus) ===\n")
cat("CL =", round(init_params["CL"], 6), "\n")
cat("V1 =", round(init_params["V1"], 5), "\n")
cat("V2 =", round(init_params["V2"], 5), "\n")
cat("Q  =", round(init_params["Q"],  6), "\n")

# =============================================================================
# 4. FONCTION OBJECTIVE — résidus log, poids égaux par groupe
# =============================================================================

objective <- function(logpar) {
  par        <- exp(logpar)
  names(par) <- c("CL", "V1", "V2", "Q")

  p10 <- tryCatch(sim_dose(dose10, par, d10$t), error = function(e) NULL)
  p5  <- tryCatch(sim_dose(dose5,  par, d5$t),  error = function(e) NULL)
  p1  <- tryCatch(sim_dose(dose1,  par, d1$t),  error = function(e) NULL)

  if (is.null(p10) || any(is.na(p10) | p10 <= 0)) return(1e10)
  if (is.null(p5)  || any(is.na(p5)  | p5  <= 0)) return(1e10)
  if (is.null(p1)  || any(is.na(p1)  | p1  <= 0)) return(1e10)

  mean((log(d10$C) - log(p10))^2) +
  mean((log(d5$C)  - log(p5))^2)  +
  mean((log(d1$C)  - log(p1))^2)
}

# =============================================================================
# 5. OPTIMISATION — nlminb en log-espace
# =============================================================================

cat("\nOptimisation nlminb en cours...\n")

fit <- nlminb(
  start     = log(init_params),
  objective = objective,
  control   = list(eval.max = 3000, iter.max = 1500,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

best_par        <- exp(fit$par)
names(best_par) <- c("CL", "V1", "V2", "Q")

# =============================================================================
# 6. PARAMÈTRES DÉRIVÉS
# =============================================================================

k10 <- best_par["CL"] / best_par["V1"]
k12 <- best_par["Q"]  / best_par["V1"]
k21 <- best_par["Q"]  / best_par["V2"]
Vss <- best_par["V1"] + best_par["V2"]

# Valeurs propres → α (rapide) et β (lent)
sum_k  <- k10 + k12 + k21
disc   <- sqrt((k10 + k12 - k21)^2 + 4 * k12 * k21)
alpha  <- (sum_k + disc) / 2
beta   <- (sum_k - disc) / 2

t_half_alpha <- log(2) / alpha
t_half_beta  <- log(2) / beta

cat("\n=== Paramètres PK 2-compartiments (rxode2 + nlminb) ===\n")
cat("CL  =", round(best_par["CL"], 6), "L/h/kg  →",
             round(best_par["CL"] * 24, 4), "L/j/kg\n")
cat("V1  =", round(best_par["V1"], 5), "L/kg\n")
cat("V2  =", round(best_par["V2"], 5), "L/kg\n")
cat("Q   =", round(best_par["Q"],  6), "L/h/kg\n")
cat("Vss =", round(Vss, 5), "L/kg\n")
cat("k10 =", round(k10, 6), "/h\n")
cat("k12 =", round(k12, 6), "/h\n")
cat("k21 =", round(k21, 6), "/h\n")
cat("t½α =", round(t_half_alpha, 2), "h\n")
cat("t½β =", round(t_half_beta,  1), "h  =",
             round(t_half_beta / 24, 2), "jours\n")
cat("Objectif final :", round(fit$objective, 6), "\n")

# =============================================================================
# 7. SIMULATION FINALE ET GRAPHIQUE
# =============================================================================

times_full <- seq(0, max(df$t) * 1.05, by = 1)

df_sim <- do.call(rbind, lapply(
  list(list(dose10, "10 mg/kg"),
       list(dose5,  "5 mg/kg"),
       list(dose1,  "1 mg/kg")),
  function(x) data.frame(
    t    = times_full,
    C    = sim_dose(x[[1]], best_par, times_full),
    Dose = x[[2]]
  )
))
df_sim$Dose <- factor(df_sim$Dose, levels = c("10 mg/kg", "5 mg/kg", "1 mg/kg"))

ggplot() +
  geom_line(data = df_sim, aes(x = t, y = C, color = Dose), linewidth = 1) +
  geom_point(data = df,    aes(x = t, y = C, color = Dose), size = 2.5) +
  scale_y_log10() +
  labs(
    title    = "PK 2-compartiments (rxode2) — Fc-silent B/C huBPA-LP1 (FGFR2)",
    subtitle = paste0(
      "V1 = ", round(best_par["V1"], 4),
      "  V2 = ", round(best_par["V2"], 4),
      "  CL = ", round(best_par["CL"], 6),
      "  t½α = ", round(t_half_alpha, 1), "h",
      "  t½β = ", round(t_half_beta,  1), "h"
    ),
    x = "Temps (heures)",
    y = "Concentration (échelle log)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PK2comp_rxode2_FGFR2.png", width = 8, height = 5, dpi = 150)
cat("\nGraphique → scripts/plot_PK2comp_rxode2_FGFR2.png\n")

# =============================================================================
# 8. SAUVEGARDE
# =============================================================================

pk2comp_rxode2 <- list(
  CL = best_par["CL"], V1 = best_par["V1"],
  V2 = best_par["V2"], Q  = best_par["Q"],
  k10 = k10, k12 = k12, k21 = k21, Vss = Vss,
  alpha = alpha, beta = beta,
  t_half_alpha = t_half_alpha, t_half_beta = t_half_beta
)

save(pk2comp_rxode2, df,
     file = "scripts/resultats_PK2comp_rxode2_FGFR2.RData")
cat("Résultats → scripts/resultats_PK2comp_rxode2_FGFR2.RData\n")
