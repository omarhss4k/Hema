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

# Doses converties en µg/kg pour cohérence avec les concentrations observées
# (ng/mL ou unité équivalente) → V1 se retrouve dans une plage raisonnable (~0.1-10)
dose10 <- 10 * 1000   # µg/kg
dose5  <-  5 * 1000
dose1  <-  1 * 1000

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
  V1_val <- as.numeric(params["V1"])   # strip name to éviter "C1_10.V1"
  state0 <- c(
    C1_10 = dose10 / V1_val,
    C2_10 = 0,
    C1_5  = dose5  / V1_val,
    C2_5  = 0,
    C1_1  = dose1  / V1_val,
    C2_1  = 0
  )

  out <- tryCatch(
    lsoda(y=state0, times=times,
          func=pk_model, parms=params,
          rtol=1e-6, atol=1e-8),
    error   = function(e) NULL,
    warning = function(w) {
      # lsoda peut émettre un warning et retourner un résultat partiel
      suppressWarnings(
        lsoda(y=state0, times=times,
              func=pk_model, parms=params,
              rtol=1e-6, atol=1e-8)
      )
    }
  )
  if (is.null(out)) return(NULL)
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
  # Vérifier que tous les temps demandés sont dans la sortie
  if (max(sim$time) < max(all_times) * 0.99) return(1e10)

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
# 5. VALEURS INITIALES ESTIMÉES DEPUIS LES DONNÉES
# =============================================================================

# V1 : dose / concentration au premier point (t=4h ≈ t=0 pour t_half long)
V1_init <- dose10 / d10$c[1]

# k_elim : pente terminale log-linéaire sur les 3 derniers points de la dose max
n_last  <- min(3, nrow(d10))
lm_term <- lm(log(c) ~ t, data = tail(d10, n_last))
k_elim  <- abs(coef(lm_term)[2])   # /h
CL_init <- k_elim * V1_init

# V2, Q : plusieurs hypothèses pour briser la symétrie V1=V2
# Pour un anticorps : V2 > V1 (distribution tissulaire), Q >> CL (échange rapide)
V2_init <- 2 * V1_init
Q_init  <- 20 * CL_init

cat("Valeurs initiales estimées depuis les données :\n")
cat("  V1 =", round(V1_init,  5), "\n")
cat("  CL =", round(CL_init,  8), "\n")
cat("  V2 =", round(V2_init,  5), "(= 2 × V1)\n")
cat("  Q  =", round(Q_init,   6), "(= 20 × CL)\n")

# =============================================================================
# 6. OPTIMISATION — nlminb sur paramètres log-transformés (multi-départ)
# Travailler en log-espace : pas de bornes, gère tous les ordres de grandeur
# =============================================================================

objective_pk_log <- function(logpar) {
  par <- exp(logpar)
  names(par) <- c("CL", "V1", "V2", "Q")
  objective_pk(par)
}

# 5 points de départ avec des asymétries V1/V2 différentes
start_pts <- list(
  # 1. Typique anticorps : V2>V1, Q>>CL
  c(CL=CL_init,        V1=V1_init,     V2=2*V1_init,   Q=20*CL_init),
  # 2. V2 encore plus grand
  c(CL=CL_init,        V1=V1_init,     V2=5*V1_init,   Q=50*CL_init),
  # 3. V2 petit, Q grand (distribution rapide)
  c(CL=CL_init,        V1=V1_init,     V2=0.5*V1_init, Q=30*CL_init),
  # 4. CL plus grand, V2 moyen
  c(CL=CL_init*3,      V1=V1_init,     V2=3*V1_init,   Q=10*CL_init),
  # 5. Départ "1-comp" comme avant — gardé pour comparaison
  c(CL=CL_init,        V1=V1_init,     V2=V1_init,     Q=CL_init)
)

best_obj <- Inf
best_fit <- NULL
for (i in seq_along(start_pts)) {
  s <- pmax(start_pts[[i]], 1e-12)   # éviter log(0)
  cat("  Départ", i, "...\n")
  fit_try <- tryCatch(
    nlminb(log(s), objective_pk_log,
           control = list(eval.max=3000, iter.max=1500,
                          rel.tol=1e-12, x.tol=1e-12)),
    error = function(e) NULL
  )
  if (!is.null(fit_try) && is.finite(fit_try$objective) &&
      fit_try$objective < best_obj) {
    best_obj <- fit_try$objective
    best_fit <- fit_try
  }
}

if (is.null(best_fit)) stop("Aucun point de départ n'a convergé.")

best_params_pk        <- exp(best_fit$par)
names(best_params_pk) <- c("CL", "V1", "V2", "Q")

cat("\n--- nlminb (multi-départ, log-espace) ---\n")
cat("CL =", round(best_params_pk["CL"], 8), "\n")
cat("V1 =", round(best_params_pk["V1"], 5), "\n")
cat("V2 =", round(best_params_pk["V2"], 5), "\n")
cat("Q  =", round(best_params_pk["Q"],  5), "\n")
cat("Objectif final :", best_fit$objective, "\n")

# Demi-vie terminale approximative
t_half <- log(2) * best_params_pk["V1"] / best_params_pk["CL"]
cat("\nDemi-vie approx. :", round(t_half, 1), "h  =",
    round(t_half/24, 2), "jours\n")

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
