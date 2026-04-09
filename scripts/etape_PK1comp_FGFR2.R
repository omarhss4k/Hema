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

# =============================================================================
# ESTIMATION — une seule régression linéaire sur log(C/dose) ~ t
# Hypothèse : PK linéaire (mêmes V1 et k pour les 3 doses)
# =============================================================================

fit <- lm(log(C / dose) ~ t, data = df)

k_elim <- -coef(fit)[["t"]]              # /h
V1     <-  exp(-coef(fit)[["(Intercept)"]])
CL     <-  k_elim * V1                   # même unités/h
t_half <-  log(2) / k_elim               # heures

cat("\n=== PK 1-compartiment ===\n")
cat("V1     =", round(V1,     5), "\n")
cat("k_elim =", round(k_elim, 6), "/h\n")
cat("CL     =", round(CL,     7), "\n")
cat("t½     =", round(t_half, 1), "h  =", round(t_half/24, 2), "jours\n")
cat("CL (L/j/kg) =", round(CL * 24, 5), "\n")
cat("R²     =", round(summary(fit)$r.squared, 4), "\n")

# =============================================================================
# GRAPHIQUE
# =============================================================================

times_full <- seq(0, max(df$t) * 1.05, by = 1)

df_sim <- do.call(rbind, lapply(
  list(c(dose10,"10 mg/kg"), c(dose5,"5 mg/kg"), c(dose1,"1 mg/kg")),
  function(x) data.frame(
    t    = times_full,
    C    = as.numeric(x[1]) * exp(coef(fit)[1]) * exp(-k_elim * times_full),
    Dose = x[2]
  )
))
df_sim$Dose <- factor(df_sim$Dose, levels=c("10 mg/kg","5 mg/kg","1 mg/kg"))

ggplot() +
  geom_line(data=df_sim, aes(x=t, y=C, color=Dose), linewidth=1) +
  geom_point(data=df,    aes(x=t, y=C, color=Dose), size=2.5) +
  scale_y_log10() +
  labs(
    title    = "PK 1-comp — Fc-silent B/C huBPA-LP1 (FGFR2)",
    subtitle = paste0("V1=", round(V1,4),
                      "  CL=", round(CL,6),
                      "  t½=", round(t_half,1), "h"),
    x = "Temps (heures)",
    y = "Concentration (échelle log)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PK1comp_FGFR2.png", width=8, height=5, dpi=150)

# --- Sauvegarde ---
save(V1, CL, k_elim, t_half, df,
     file = "scripts/resultats_PK1comp_FGFR2.RData")

cat("\nSauvegardé → scripts/resultats_PK1comp_FGFR2.RData\n")
