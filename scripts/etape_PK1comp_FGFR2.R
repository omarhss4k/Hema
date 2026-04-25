# =============================================================================
# PK 1-compartiment — régression log-linéaire
# C(t) = (dose/V1) * exp(-k * t)
# log(C/dose) = -log(V1) - k*t  → régression linéaire simple
# Paramètres : V1, k_elim (→ CL = k * V1)
# =============================================================================

library(ggplot2)
library(readxl)

# --- Données ---
raw <- read_xlsx("data/PK souris FGFR2.xlsx")

time_h <- as.numeric(raw[[1]])
conc10 <- as.numeric(raw[[2]])
conc5  <- as.numeric(raw[[3]])
conc1  <- as.numeric(raw[[4]])

dose10 <- 10 * 1000   # µg/kg
dose5  <-  5 * 1000
dose1  <-  1 * 1000

# Assembler les 3 groupes (retirer NA)
mk <- function(t, c, d, label) {
  ok <- !is.na(c) & c > 0
  data.frame(t=t[ok], C=c[ok], dose=d, Dose=label)
}
df <- rbind(
  mk(time_h, conc10, dose10, "10 mg/kg"),
  mk(time_h, conc5,  dose5,  "5 mg/kg"),
  mk(time_h, conc1,  dose1,  "1 mg/kg")
)
df$Dose <- factor(df$Dose, levels=c("10 mg/kg","5 mg/kg","1 mg/kg"))

cat("Points totaux :", nrow(df), "\n")

# Seuil : on ne garde que la phase terminale (après la distribution)
t_terminal <- 168   # heures — premier point hors phase de distribution
df_term <- df[df$t >= t_terminal, ]
cat("Points phase terminale (t ≥", t_terminal, "h) :", nrow(df_term), "\n")

# =============================================================================
# ESTIMATION — régression sur la phase terminale uniquement
# Hypothèse : PK linéaire, 1-comp, phase de distribution terminée à t≥168h
# =============================================================================

fit <- lm(log(C / dose) ~ t, data = df_term)

# Valeurs initiales depuis la régression
k_init  <- -coef(fit)[["t"]]
V1_init <-  exp(-coef(fit)[["(Intercept)"]])

cat("\n--- Régression log-linéaire phase terminale (départ) ---\n")
cat("V1 =", round(V1_init, 5), "  k =", round(k_init, 6), "/h\n")
cat("R² =", round(summary(fit)$r.squared, 4), "\n")

# =============================================================================
# OPTIMISATION — nlminb sur paramètres log-transformés
# Objectif : résidus log pondérés équitablement par groupe
# =============================================================================

d10 <- df_term[df_term$Dose == "10 mg/kg", ]
d5  <- df_term[df_term$Dose == "5 mg/kg",  ]
d1  <- df_term[df_term$Dose == "1 mg/kg",  ]

objective_1comp <- function(logpar) {
  k  <- exp(logpar[1])
  V1 <- exp(logpar[2])

  pred <- function(dose, t) (dose / V1) * exp(-k * t)

  r10 <- log(d10$C) - log(pred(dose10, d10$t))
  r5  <- log(d5$C)  - log(pred(dose5,  d5$t))
  r1  <- log(d1$C)  - log(pred(dose1,  d1$t))

  # Poids égaux par groupe (indépendant du nombre de points par groupe)
  mean(r10^2) + mean(r5^2) + mean(r1^2)
}

fit_opt <- nlminb(
  start   = log(c(k_init, V1_init)),
  objective = objective_1comp,
  control = list(eval.max = 2000, iter.max = 1000,
                 rel.tol = 1e-12, x.tol = 1e-12)
)

k_elim <- exp(fit_opt$par[1])
V1     <- exp(fit_opt$par[2])
CL     <- k_elim * V1
t_half <- log(2) / k_elim

cat("\n=== nlminb (affiné, poids égaux par groupe) ===\n")
cat("V1     =", round(V1,     5), "\n")
cat("k_elim =", round(k_elim, 6), "/h\n")
cat("CL     =", round(CL,     7), "\n")
cat("t½     =", round(t_half, 1), "h  =", round(t_half/24, 2), "jours\n")
cat("CL (L/j/kg) =", round(CL * 24, 5), "\n")
cat("Objectif    =", round(fit_opt$objective, 6), "\n")

# =============================================================================
# GRAPHIQUE
# =============================================================================

times_full <- seq(0, max(df$t) * 1.05, by = 1)

df_sim <- do.call(rbind, lapply(
  list(c(dose10,"10 mg/kg"), c(dose5,"5 mg/kg"), c(dose1,"1 mg/kg")),
  function(x) data.frame(
    t    = times_full,
    C    = (as.numeric(x[1]) / V1) * exp(-k_elim * times_full),
    Dose = x[2]
  )
))
df_sim$Dose <- factor(df_sim$Dose, levels=c("10 mg/kg","5 mg/kg","1 mg/kg"))

ggplot() +
  geom_line(data=df_sim, aes(x=t, y=C, color=Dose), linewidth=1) +
  # Tous les points observés (phase distribution en transparent)
  geom_point(data=df[df$t < t_terminal, ],
             aes(x=t, y=C, color=Dose), size=2.5, alpha=0.3, shape=1) +
  # Points phase terminale (pleins)
  geom_point(data=df[df$t >= t_terminal, ],
             aes(x=t, y=C, color=Dose), size=2.5) +
  geom_vline(xintercept=t_terminal, linetype="dashed", color="grey50") +
  annotate("text", x=t_terminal+10, y=max(df$C)*0.7,
           label=paste0("t = ", t_terminal, "h\n(phase terminale)"),
           hjust=0, size=3, color="grey40") +
  scale_y_log10() +
  labs(
    title    = "PK 1-comp phase terminale — Fc-silent B/C huBPA-LP1 (FGFR2)",
    subtitle = paste0("V1=", round(V1,4),
                      "  CL=", round(CL,6),
                      "  t½=", round(t_half,1), "h  (ajusté sur t≥",
                      t_terminal, "h)"),
    x = "Temps (heures)",
    y = "Concentration (échelle log)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PK1comp_FGFR2.png", width=8, height=5, dpi=150)

# --- Sauvegarde ---
save(V1, CL, k_elim, t_half, df,
     file = "scripts/resultats_PK1comp_FGFR2.RData")

cat("\nSauvegardé → scripts/resultats_PK1comp_FGFR2.RData\n")
