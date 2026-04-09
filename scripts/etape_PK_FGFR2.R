# =============================================================================
# ÉTAPE PK — Estimation (CL, V1, V2, Q)
# Médicament : Fc-silent B/C huBPA-LP1 (FGFR2)
# Modèle : 2 compartiments, administration IV bolus
# Doses   : 1, 5, 10 mg/kg
# Temps   : heures
# Conc.   : unité inconnue (cohérente dans tout le fichier)
# =============================================================================

library(deSolve)
library(ggplot2)
library(readxl)
library(DEoptim)

# =============================================================================
# 1. LECTURE ET PRÉPARATION DES DONNÉES
# =============================================================================

raw <- read_xlsx("data/PK souris FGFR2.xlsx")

# Colonnes : A=temps(h), B=10mg/kg, C=5mg/kg, D=1mg/kg
time_h <- as.numeric(raw[[1]])
conc10 <- as.numeric(raw[[2]])   # 10 mg/kg
conc5  <- as.numeric(raw[[3]])   #  5 mg/kg
conc1  <- as.numeric(raw[[4]])   #  1 mg/kg

# Séparation par groupe (retirer NA) et tri
mk_df <- function(t, c) {
  df <- data.frame(t = t, c = c)
  df <- df[!is.na(df$c) & df$c > 0, ]
  df[order(df$t), ]
}

d10 <- mk_df(time_h, conc10)
d5  <- mk_df(time_h, conc5)
d1  <- mk_df(time_h, conc1)

cat("Groupe 10 mg/kg :", nrow(d10), "points\n")
cat("Groupe  5 mg/kg :", nrow(d5),  "points\n")
cat("Groupe  1 mg/kg :", nrow(d1),  "points\n")

# Doses (mg/kg)
dose10 <- 10
dose5  <-  5
dose1  <-  1

# =============================================================================
# 2. MODÈLE ODE — 2 compartiments, 3 groupes (6 équations)
# dC1/dt = -(CL/V1)*C1 - (Q/V1)*C1 + (Q/V2)*C2
# dC2/dt =  (Q/V1)*C1  - (Q/V2)*C2
# =============================================================================

pk_model <- function(time, state, params) {
  with(as.list(c(state, params)), {

    dC1_10 <- -(CL/V1)*C1_10 - (Q/V1)*C1_10 + (Q/V2)*C2_10
    dC2_10 <-  (Q/V1)*C1_10  - (Q/V2)*C2_10

    dC1_5  <- -(CL/V1)*C1_5  - (Q/V1)*C1_5  + (Q/V2)*C2_5
    dC2_5  <-  (Q/V1)*C1_5   - (Q/V2)*C2_5

    dC1_1  <- -(CL/V1)*C1_1  - (Q/V1)*C1_1  + (Q/V2)*C2_1
    dC2_1  <-  (Q/V1)*C1_1   - (Q/V2)*C2_1

    list(c(dC1_10, dC2_10,
           dC1_5,  dC2_5,
           dC1_1,  dC2_1))
  })
}

# =============================================================================
# 3. SIMULATION
# Administration IV bolus à t=0 : C1(0) = dose / V1
# =============================================================================

simulate_pk <- function(times, params) {
  state0 <- c(
    C1_10 = dose10 / params["V1"],
    C2_10 = 0,
    C1_5  = dose5  / params["V1"],
    C2_5  = 0,
    C1_1  = dose1  / params["V1"],
    C2_1  = 0
  )

  out <- lsoda(y=state0, times=times,
               func=pk_model, parms=params,
               rtol=1e-8, atol=1e-10)
  as.data.frame(out)
}

# =============================================================================
# 4. FONCTION OBJECTIVE (résidus log)
# =============================================================================

objective_pk <- function(par) {
  params <- c(CL=par[1], V1=par[2], V2=par[3], Q=par[4])

  all_times <- sort(unique(c(0, d10$t, d5$t, d1$t)))
  sim <- tryCatch(simulate_pk(all_times, params), error = function(e) NULL)
  if (is.null(sim) || nrow(sim) < 2) return(1e10)

  p10 <- approx(sim$time, sim$C1_10, xout = d10$t)$y
  p5  <- approx(sim$time, sim$C1_5,  xout = d5$t)$y
  p1  <- approx(sim$time, sim$C1_1,  xout = d1$t)$y

  if (any(is.na(p10) | p10 <= 0)) return(1e10)
  if (any(is.na(p5)  | p5  <= 0)) return(1e10)
  if (any(is.na(p1)  | p1  <= 0)) return(1e10)

  val <- sum((log(d10$c) - log(p10))^2) +
         sum((log(d5$c)  - log(p5))^2)  +
         sum((log(d1$c)  - log(p1))^2)

  if (!is.finite(val)) return(1e10)
  val
}

# =============================================================================
# 5. OPTIMISATION GLOBALE — DEoptim
# =============================================================================

set.seed(42)
res_de <- DEoptim(
  fn    = objective_pk,
  lower = c(CL=1e-4, V1=1e-4, V2=1e-4, Q=1e-5),
  upper = c(CL=100,  V1=50,   V2=50,   Q=50),
  control = DEoptim.control(
    NP        = 80,
    itermax   = 500,
    F         = 0.8,
    CR        = 0.9,
    trace     = 100,
    parallelType = 0
  )
)

best_de <- res_de$optim$bestmem
cat("\n--- DEoptim ---\n")
cat("CL =", round(best_de["CL"], 5), "(mg/kg)/(conc·h)\n")
cat("V1 =", round(best_de["V1"], 5), "(mg/kg)/conc\n")
cat("V2 =", round(best_de["V2"], 5), "(mg/kg)/conc\n")
cat("Q  =", round(best_de["Q"],  5), "(mg/kg)/(conc·h)\n")

# =============================================================================
# 6. AFFINAGE LOCAL — nlminb
# =============================================================================

fit <- nlminb(
  start     = best_de,
  objective = objective_pk,
  lower     = c(1e-4, 1e-4, 1e-4, 1e-5),
  upper     = c(100,  50,   50,   50),
  control   = list(eval.max=2000, iter.max=1000,
                   rel.tol=1e-12, x.tol=1e-12)
)

cat("\n--- nlminb (affiné) ---\n")
cat("CL =", round(fit$par["CL"], 5), "\n")
cat("V1 =", round(fit$par["V1"], 5), "\n")
cat("V2 =", round(fit$par["V2"], 5), "\n")
cat("Q  =", round(fit$par["Q"],  5), "\n")
cat("Objectif final :", fit$objective, "\n")

# Demi-vie terminale approximative (heures)
pk_params <- fit$par
t_half <- log(2) * pk_params["V1"] / pk_params["CL"]
cat("\nDemi-vie approx. :", round(t_half, 1), "h  =",
    round(t_half/24, 2), "jours\n")

best_params_pk <- fit$par

# =============================================================================
# 7. SIMULATION FINALE ET GRAPHIQUE
# =============================================================================

times_full <- seq(0, max(c(d10$t, d5$t, d1$t)) * 1.05, by = 1)
sim_final  <- simulate_pk(times_full, best_params_pk)

df_sim <- data.frame(
  time = rep(sim_final$time, 3),
  conc = c(sim_final$C1_10, sim_final$C1_5, sim_final$C1_1),
  Dose = rep(c("10 mg/kg","5 mg/kg","1 mg/kg"), each = nrow(sim_final))
)

df_obs <- data.frame(
  time = c(d10$t, d5$t, d1$t),
  conc = c(d10$c, d5$c, d1$c),
  Dose = c(rep("10 mg/kg", nrow(d10)),
           rep("5 mg/kg",  nrow(d5)),
           rep("1 mg/kg",  nrow(d1)))
)

df_sim$Dose <- factor(df_sim$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))
df_obs$Dose <- factor(df_obs$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))

ggplot() +
  geom_line(data = df_sim, aes(x=time, y=conc, color=Dose), linewidth=1) +
  geom_point(data = df_obs, aes(x=time, y=conc, color=Dose), size=2.5) +
  scale_y_log10() +
  labs(
    title    = "PK — Fc-silent B/C huBPA-LP1 (FGFR2)",
    subtitle = paste0(
      "CL=", round(best_params_pk["CL"], 4),
      "  V1=", round(best_params_pk["V1"], 4),
      "  V2=", round(best_params_pk["V2"], 4),
      "  Q=",  round(best_params_pk["Q"],  4),
      "  t½≈", round(t_half, 1), "h"
    ),
    x = "Temps (heures)",
    y = "Concentration (unité inconnue, échelle log)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PK_FGFR2.png", width=8, height=5, dpi=150)

# =============================================================================
# 8. SAUVEGARDE — paramètres PK à réutiliser pour l'étape PKPD
# =============================================================================

save(best_params_pk, t_half, sim_final, d10, d5, d1,
     file = "scripts/resultats_PK_FGFR2.RData")

cat("\nParamètres PK sauvegardés → scripts/resultats_PK_FGFR2.RData\n")
cat("Ces paramètres seront fixés lors de l'estimation PKPD (k1, k2)\n")
cat("quand les données de volumes tumoraux seront disponibles.\n")
