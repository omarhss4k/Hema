# =============================================================================
# PKPD — Fc-silent B/C huBPA-LP1 (FGFR2) — Modèle Simeoni 2004
# 3 doses : 1, 5, 10 mg/kg
# Paramètres PD de croissance fixés (stagiaire) : l0, l1, p
# Paramètres à estimer : CL, V1, V2, Q, k1, k2
# =============================================================================

library(deSolve)
library(ggplot2)
library(readxl)
library(DEoptim)

# =============================================================================
# 1. LECTURE ET PRÉPARATION DES DONNÉES
# =============================================================================

raw <- read_xlsx("data/PK souris FGFR2.xlsx")

# Colonnes : A=temps(h), B=10mg/kg, C=5mg/kg, D=1mg/kg  (colonne E ignorée)
data_raw <- data.frame(
  time_h  = as.numeric(raw[[1]]),
  d10     = as.numeric(raw[[2]]),   # 10 mg/kg
  d5      = as.numeric(raw[[3]]),   # 5  mg/kg
  d1      = as.numeric(raw[[4]])    # 1  mg/kg
)

# Conversion heures → jours
data_raw$time <- data_raw$time_h / 24

# Retirer les lignes entièrement NA sur les colonnes de mesure
data_raw <- data_raw[rowSums(!is.na(data_raw[, c("d10","d5","d1")])) > 0, ]
data_raw <- data_raw[order(data_raw$time), ]

cat("Temps (jours) :", round(data_raw$time, 3), "\n")

# Séparation par groupe (retirer NA)
data10 <- data_raw[!is.na(data_raw$d10), c("time","d10")]; names(data10) <- c("x","y")
data5  <- data_raw[!is.na(data_raw$d5),  c("time","d5")];  names(data5)  <- c("x","y")
data1  <- data_raw[!is.na(data_raw$d1),  c("time","d1")];  names(data1)  <- c("x","y")

# =============================================================================
# 2. ESTIMATION DE w0 (masse tumorale initiale)
# Régression log-linéaire sur les points pré-traitement (time < t_admin)
# =============================================================================

t_admin <- 13   # jour d'administration (jours)

pre10 <- data10[data10$x < t_admin, ]
pre5  <- data5[data5$x  < t_admin, ]
pre1  <- data1[data1$x  < t_admin, ]

estimate_w0 <- function(df) {
  if (nrow(df) < 2) return(NA)
  fit <- lm(log(y) ~ x, data = df)
  exp(predict(fit, newdata = data.frame(x = 0)))
}

w0_10 <- estimate_w0(pre10)
w0_5  <- estimate_w0(pre5)
w0_1  <- estimate_w0(pre1)

# Moyenne des w0 disponibles
w0_vals <- c(w0_10, w0_5, w0_1)
w0 <- mean(w0_vals, na.rm = TRUE)
cat("\nw0 estimés :", round(w0_vals, 4), "\n")
cat("w0 commun  :", round(w0, 4), "\n")

# =============================================================================
# 3. PARAMÈTRES PD DE CROISSANCE (fixés — stagiaire)
# =============================================================================

l0 <- 0.146
l1 <- 0.334
p  <- 20

# Doses (converties en µg/kg depuis mg/kg)
dose10 <- 10 * 1000   # µg/kg
dose5  <-  5 * 1000
dose1  <-  1 * 1000

# =============================================================================
# 4. MODÈLE ODE PKPD — 3 groupes, 21 équations
# États : C1_i, C2_i (PK), x1_i..x4_i (PD), w_i (masse totale)
# i = 1 (10mg/kg), 2 (5mg/kg), 3 (1mg/kg)
# =============================================================================

combined_model <- function(time, state, params) {
  with(as.list(c(state, params)), {

    # ---- Groupe 1 : 10 mg/kg ----
    dC1_1 <- -(CL/V1)*C1_1 - (Q/V1)*C1_1 + (Q/V2)*C2_1
    dC2_1 <-  (Q/V1)*C1_1  - (Q/V2)*C2_1
    gr1   <- l0*x1_1 / ((1 + (l0/l1 * w_1)^p)^(1/p))
    dx1_1 <- gr1 - k2*C1_1*x1_1
    dx2_1 <- k2*C1_1*x1_1 - k1*x2_1
    dx3_1 <- k1*(x2_1 - x3_1)
    dx4_1 <- k1*(x3_1 - x4_1)
    dw_1  <- dx1_1 + dx2_1 + dx3_1 + dx4_1

    # ---- Groupe 2 : 5 mg/kg ----
    dC1_2 <- -(CL/V1)*C1_2 - (Q/V1)*C1_2 + (Q/V2)*C2_2
    dC2_2 <-  (Q/V1)*C1_2  - (Q/V2)*C2_2
    gr2   <- l0*x1_2 / ((1 + (l0/l1 * w_2)^p)^(1/p))
    dx1_2 <- gr2 - k2*C1_2*x1_2
    dx2_2 <- k2*C1_2*x1_2 - k1*x2_2
    dx3_2 <- k1*(x2_2 - x3_2)
    dx4_2 <- k1*(x3_2 - x4_2)
    dw_2  <- dx1_2 + dx2_2 + dx3_2 + dx4_2

    # ---- Groupe 3 : 1 mg/kg ----
    dC1_3 <- -(CL/V1)*C1_3 - (Q/V1)*C1_3 + (Q/V2)*C2_3
    dC2_3 <-  (Q/V1)*C1_3  - (Q/V2)*C2_3
    gr3   <- l0*x1_3 / ((1 + (l0/l1 * w_3)^p)^(1/p))
    dx1_3 <- gr3 - k2*C1_3*x1_3
    dx2_3 <- k2*C1_3*x1_3 - k1*x2_3
    dx3_3 <- k1*(x2_3 - x3_3)
    dx4_3 <- k1*(x3_3 - x4_3)
    dw_3  <- dx1_3 + dx2_3 + dx3_3 + dx4_3

    list(c(
      dC1_1, dC2_1, dx1_1, dx2_1, dx3_1, dx4_1, dw_1,
      dC1_2, dC2_2, dx1_2, dx2_2, dx3_2, dx4_2, dw_2,
      dC1_3, dC2_3, dx1_3, dx2_3, dx3_3, dx4_3, dw_3
    ))
  })
}

# =============================================================================
# 5. SIMULATION EN 2 SEGMENTS
# =============================================================================

simulate_model <- function(times, params) {

  state0 <- c(
    C1_1=0, C2_1=0, x1_1=w0, x2_1=0, x3_1=0, x4_1=0, w_1=w0,
    C1_2=0, C2_2=0, x1_2=w0, x2_2=0, x3_2=0, x4_2=0, w_2=w0,
    C1_3=0, C2_3=0, x1_3=w0, x2_3=0, x3_3=0, x4_3=0, w_3=w0
  )

  # Segment 1 : 0 → t_admin (croissance libre)
  t_before <- sort(unique(c(0, times[times <= t_admin], t_admin)))
  out1 <- lsoda(y=state0, times=t_before,
                func=combined_model, parms=params,
                rtol=1e-8, atol=1e-10)
  df1 <- as.data.frame(out1)

  # Injection à t_admin
  state_inj        <- as.numeric(tail(df1, 1)[, -1])
  names(state_inj) <- names(state0)
  state_inj["C1_1"] <- dose10 / params["V1"]
  state_inj["C1_2"] <- dose5  / params["V1"]
  state_inj["C1_3"] <- dose1  / params["V1"]

  # Segment 2 : t_admin → fin
  t_after <- sort(unique(c(t_admin, times[times >= t_admin])))
  out2 <- lsoda(y=state_inj, times=t_after,
                func=combined_model, parms=params,
                rtol=1e-8, atol=1e-10)
  df2 <- as.data.frame(out2)

  # Fusion sans doublon à t_admin
  rbind(df1[df1$time < t_admin, ], df2)
}

# =============================================================================
# 6. FONCTION OBJECTIVE (résidus log)
# =============================================================================

objective_function <- function(par) {
  params <- c(CL=par[1], V1=par[2], V2=par[3], Q=par[4],
              k1=par[5], k2=par[6])

  all_times <- sort(unique(c(data10$x, data5$x, data1$x)))
  sim <- tryCatch(simulate_model(all_times, params), error = function(e) NULL)
  if (is.null(sim) || nrow(sim) < 2) return(1e10)

  pred10 <- approx(sim$time, sim$w_1, xout = data10$x)$y
  pred5  <- approx(sim$time, sim$w_2, xout = data5$x)$y
  pred1  <- approx(sim$time, sim$w_3, xout = data1$x)$y

  if (any(is.na(pred10) | pred10 <= 0)) return(1e10)
  if (any(is.na(pred5)  | pred5  <= 0)) return(1e10)
  if (any(is.na(pred1)  | pred1  <= 0)) return(1e10)

  val <- sum((log(data10$y) - log(pred10))^2) +
         sum((log(data5$y)  - log(pred5))^2)  +
         sum((log(data1$y)  - log(pred1))^2)

  if (!is.finite(val)) return(1e10)
  val
}

# =============================================================================
# 7. OPTIMISATION GLOBALE — DEoptim
# Paramètres : CL, V1, V2, Q, k1, k2
# Unités PK en jours (cohérent avec t_admin en jours)
# =============================================================================

set.seed(42)
res_de <- DEoptim(
  fn    = objective_function,
  lower = c(CL=0.01, V1=0.01, V2=0.005, Q=0.001, k1=0.1,  k2=1e-7),
  upper = c(CL=50,   V1=5,    V2=5,     Q=10,    k1=5.0,  k2=0.01),
  control = DEoptim.control(
    NP        = 100,
    itermax   = 600,
    F         = 0.8,
    CR        = 0.9,
    trace     = 100,
    parallelType = 0
  )
)

best_de <- res_de$optim$bestmem
cat("\n--- DEoptim ---\n")
cat("CL =", round(best_de["CL"], 5), "\n")
cat("V1 =", round(best_de["V1"], 5), "\n")
cat("V2 =", round(best_de["V2"], 5), "\n")
cat("Q  =", round(best_de["Q"],  5), "\n")
cat("k1 =", round(best_de["k1"], 5), "\n")
cat("k2 =", round(best_de["k2"], 8), "\n")

# =============================================================================
# 8. AFFINAGE LOCAL — nlminb
# =============================================================================

fit <- nlminb(
  start   = best_de,
  objective = objective_function,
  lower   = c(0.01, 0.01, 0.005, 0.001, 0.1,  1e-7),
  upper   = c(50,   5,    5,     10,    5.0,  0.01),
  control = list(eval.max=2000, iter.max=1000,
                 rel.tol=1e-12, x.tol=1e-12)
)

cat("\n--- nlminb (affiné) ---\n")
cat("CL =", round(fit$par["CL"], 5), "\n")
cat("V1 =", round(fit$par["V1"], 5), "\n")
cat("V2 =", round(fit$par["V2"], 5), "\n")
cat("Q  =", round(fit$par["Q"],  5), "\n")
cat("k1 =", round(fit$par["k1"], 5), "\n")
cat("k2 =", round(fit$par["k2"], 8), "\n")
cat("Objectif final :", fit$objective, "\n")

best_params <- fit$par

# =============================================================================
# 9. SIMULATION FINALE
# =============================================================================

times_full <- seq(0, 30, by = 0.1)
sim_final  <- simulate_model(times_full, best_params)

# =============================================================================
# 10. GRAPHIQUE
# =============================================================================

df_sim <- data.frame(
  time = rep(sim_final$time, 3),
  w    = c(sim_final$w_1, sim_final$w_2, sim_final$w_3),
  Dose = rep(c("10 mg/kg", "5 mg/kg", "1 mg/kg"), each = nrow(sim_final))
)

df_obs <- data.frame(
  time = c(data10$x, data5$x, data1$x),
  w    = c(data10$y, data5$y, data1$y),
  Dose = c(rep("10 mg/kg", nrow(data10)),
           rep("5 mg/kg",  nrow(data5)),
           rep("1 mg/kg",  nrow(data1)))
)

df_sim$Dose <- factor(df_sim$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))
df_obs$Dose <- factor(df_obs$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))

ggplot() +
  geom_line(data = df_sim, aes(x=time, y=w, color=Dose), linewidth=1) +
  geom_point(data = df_obs, aes(x=time, y=w, color=Dose), size=2.5) +
  geom_vline(xintercept = t_admin, linetype="dashed", color="grey50") +
  annotate("text", x=t_admin+0.3, y=max(df_obs$w, na.rm=TRUE)*0.95,
           label="Administration", hjust=0, size=3.5, color="grey40") +
  labs(
    title    = "PKPD — Fc-silent B/C huBPA-LP1 (FGFR2)",
    subtitle = paste0(
      "CL=", round(best_params["CL"],3),
      "  V1=", round(best_params["V1"],3),
      "  k1=", round(best_params["k1"],4),
      "  k2=", round(best_params["k2"],7)
    ),
    x = "Temps (jours)", y = "Masse tumorale (g)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PKPD_FGFR2.png", width=8, height=5, dpi=150)

# =============================================================================
# 11. SAUVEGARDE
# =============================================================================

save(best_params, w0, times_full, sim_final, data10, data5, data1,
     file = "scripts/resultats_PKPD_FGFR2.RData")

cat("\nTerminé. Résultats sauvegardés dans scripts/resultats_PKPD_FGFR2.RData\n")
