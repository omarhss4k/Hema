# =============================================================================
# PK 2-compartiments — Méthode des résidus (feathering)
# C(t) = A·exp(-α·t) + B·exp(-β·t)
# Pas d'ODE, pas d'optimisation — 2 régressions linéaires
# =============================================================================

library(ggplot2)
library(readxl)

# --- Données ---
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
  data.frame(t=t[ok], C=c[ok], dose=d, Dose=label)
}
df <- rbind(
  mk(time_h, conc10, dose10, "10 mg/kg"),
  mk(time_h, conc5,  dose5,  "5 mg/kg"),
  mk(time_h, conc1,  dose1,  "1 mg/kg")
)
df$Dose <- factor(df$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))

t_terminal <- 168   # h — séparation distribution / élimination

# =============================================================================
# ÉTAPE 1 — Phase terminale (t ≥ 168h)
# log(C/dose) = log(b) - β·t   →   β, b = B/dose
# =============================================================================

df_term  <- df[df$t >= t_terminal, ]
fit_beta <- lm(log(C / dose) ~ t, data = df_term)

beta         <- -coef(fit_beta)[["t"]]
b            <-  exp(coef(fit_beta)[["(Intercept)"]])
t_half_beta  <-  log(2) / beta

cat("=== Étape 1 — Phase terminale ===\n")
cat("β      =", round(beta, 7), "/h\n")
cat("B/dose =", round(b,    7), "\n")
cat("t½β    =", round(t_half_beta, 1), "h  =", round(t_half_beta/24, 2), "jours\n")
cat("R²     =", round(summary(fit_beta)$r.squared, 4), "\n")

# =============================================================================
# ÉTAPE 2 — Résidus sur les points précoces (t < 168h)
# C_resid = C_obs - b·dose·exp(-β·t)
# log(C_resid/dose) = log(a) - α·t   →   α, a = A/dose
# =============================================================================

df_early          <- df[df$t < t_terminal, ]
df_early$C_term   <- b * df_early$dose * exp(-beta * df_early$t)
df_early$C_resid  <- df_early$C - df_early$C_term

df_resid <- df_early[df_early$C_resid > 0, ]
cat("\nPoints résidus positifs :", nrow(df_resid), "/", nrow(df_early), "\n")

fit_alpha    <- lm(log(C_resid / dose) ~ t, data = df_resid)
alpha        <- -coef(fit_alpha)[["t"]]
a            <-  exp(coef(fit_alpha)[["(Intercept)"]])
t_half_alpha <-  log(2) / alpha

cat("\n=== Étape 2 — Phase de distribution ===\n")
cat("α      =", round(alpha, 6), "/h\n")
cat("A/dose =", round(a,     6), "\n")
cat("t½α    =", round(t_half_alpha, 2), "h\n")
cat("R²     =", round(summary(fit_alpha)$r.squared, 4), "\n")

# =============================================================================
# ÉTAPE 3 — Paramètres PK micro et macro
# =============================================================================

V1  <- 1 / (a + b)
k21 <- (a * beta + b * alpha) / (a + b)
k10 <- alpha * beta / k21
k12 <- alpha + beta - k10 - k21
CL  <- k10 * V1
Q   <- k12 * V1
V2  <- Q   / k21
Vss <- V1  + V2

cat("\n=== Paramètres PK 2-compartiments ===\n")
cat("V1  =", round(V1,  5), "L/kg\n")
cat("V2  =", round(V2,  5), "L/kg\n")
cat("Vss =", round(Vss, 5), "L/kg\n")
cat("CL  =", round(CL,  7), "L/h/kg  →", round(CL * 24, 5), "L/j/kg\n")
cat("Q   =", round(Q,   5), "L/h/kg  →", round(Q  * 24, 5), "L/j/kg\n")
cat("k10 =", round(k10, 6), "/h\n")
cat("k12 =", round(k12, 6), "/h\n")
cat("k21 =", round(k21, 6), "/h\n")
cat("t½α =", round(t_half_alpha, 2), "h\n")
cat("t½β =", round(t_half_beta,  1), "h  =", round(t_half_beta/24, 2), "jours\n")

# =============================================================================
# GRAPHIQUE
# =============================================================================

times_full <- seq(0, max(df$t) * 1.05, by = 1)

pred_2comp <- function(dose, t) {
  dose * (a * exp(-alpha * t) + b * exp(-beta * t))
}

df_sim <- do.call(rbind, lapply(
  list(c(dose10,"10 mg/kg"), c(dose5,"5 mg/kg"), c(dose1,"1 mg/kg")),
  function(x) data.frame(
    t    = times_full,
    C    = pred_2comp(as.numeric(x[1]), times_full),
    Dose = x[2]
  )
))
df_sim$Dose <- factor(df_sim$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))

ggplot() +
  geom_line(data = df_sim, aes(x=t, y=C, color=Dose), linewidth=1) +
  geom_point(data = df,    aes(x=t, y=C, color=Dose), size=2.5) +
  geom_vline(xintercept = t_terminal, linetype="dashed", color="grey50") +
  annotate("text", x=t_terminal+15, y=max(df$C)*0.6,
           label="t = 168h", hjust=0, size=3.5, color="grey40") +
  scale_y_log10() +
  labs(
    title    = "PK 2-comp (méthode des résidus) — Fc-silent B/C huBPA-LP1",
    subtitle = paste0(
      "V1=", round(V1,4), "  V2=", round(V2,4),
      "  CL=", round(CL,6),
      "  t½α=", round(t_half_alpha,1), "h",
      "  t½β=", round(t_half_beta,1), "h"
    ),
    x = "Temps (heures)",
    y = "Concentration (échelle log)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PK2comp_residus_FGFR2.png", width=8, height=5, dpi=150)

# --- Sauvegarde ---
pk2comp <- list(V1=V1, V2=V2, Vss=Vss, CL=CL, Q=Q,
                k10=k10, k12=k12, k21=k21,
                alpha=alpha, beta=beta, a=a, b=b,
                t_half_alpha=t_half_alpha, t_half_beta=t_half_beta)
save(pk2comp, df, file = "scripts/resultats_PK2comp_FGFR2.RData")
cat("\nSauvegardé → scripts/resultats_PK2comp_FGFR2.RData\n")
