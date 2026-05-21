# =============================================================================
# Schéma pédagogique — Modèle PD Emax (TGI)
# =============================================================================

library(ggplot2)
library(grid)
library(gridExtra)

# =============================================================================
# PANNEAU 1 : Schéma compartimentiel
# =============================================================================

p_schema <- ggplot() +
  xlim(0, 10) + ylim(0, 10) +
  coord_fixed() +

  # ── Boîte TV ──────────────────────────────────────────────────────────────
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 3.8, ymax = 6.2,
           fill = "#D6EAF8", color = "#2980B9", linewidth = 1.2) +
  annotate("text", x = 5, y = 5.5, label = "TV",
           size = 9, fontface = "bold", color = "#1A5276") +
  annotate("text", x = 5, y = 4.5, label = "Volume tumoral",
           size = 3.5, color = "#1A5276") +

  # ── Flèche croissance (gauche → boîte) ────────────────────────────────────
  annotate("segment", x = 1, xend = 3.4, y = 5, yend = 5,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "#27AE60", linewidth = 1.2) +
  annotate("text", x = 2.2, y = 5.5, label = "kg · TV",
           size = 4, color = "#27AE60", fontface = "bold") +
  annotate("text", x = 2.2, y = 4.6, label = "Croissance",
           size = 3.2, color = "#27AE60") +

  # ── Flèche mort (boîte → droite) ──────────────────────────────────────────
  annotate("segment", x = 6.6, xend = 9, y = 5, yend = 5,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "#E74C3C", linewidth = 1.2) +
  annotate("text", x = 7.8, y = 5.5,
           label = "E(C) · TV",
           size = 4, color = "#E74C3C", fontface = "bold") +
  annotate("text", x = 7.8, y = 4.6, label = "Mort cellulaire",
           size = 3.2, color = "#E74C3C") +

  # ── Boîte concentration C ─────────────────────────────────────────────────
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 0.5, ymax = 2.0,
           fill = "#FDEBD0", color = "#E67E22", linewidth = 1) +
  annotate("text", x = 5, y = 1.5, label = "C(t)",
           size = 7, fontface = "bold", color = "#784212") +
  annotate("text", x = 5, y = 0.9,
           label = "Concentration plasmatique\n(modèle PK)",
           size = 3, color = "#784212") +

  # ── Flèche C → flèche mort ────────────────────────────────────────────────
  annotate("segment", x = 7.8, xend = 7.8, y = 4.5, yend = 2.5,
           linetype = "dashed", color = "#E67E22", linewidth = 0.9) +
  annotate("segment", x = 6.6, xend = 7.8, y = 1.25, yend = 2.5,
           arrow = arrow(length = unit(0.25, "cm"), type = "open"),
           color = "#E67E22", linewidth = 0.9) +
  annotate("text", x = 8.5, y = 3.5, label = "module\nl'effet",
           size = 3, color = "#E67E22", hjust = 0.5) +

  # ── Équation différentielle ───────────────────────────────────────────────
  annotate("rect", xmin = 0.3, xmax = 9.7, ymin = 7.2, ymax = 9.5,
           fill = "#F8F9FA", color = "grey60", linewidth = 0.8,
           linetype = "dashed") +
  annotate("text", x = 5, y = 9.0,
           label = "Équation du modèle",
           size = 4, fontface = "bold", color = "grey30") +
  annotate("text", x = 5, y = 8.2,
           label = "dTV/dt  =  kg · TV  –  E(C) · TV",
           size = 5, color = "black", fontface = "italic") +
  annotate("text", x = 5, y = 7.5,
           label = "=  [ kg  –  E(C) ]  · TV",
           size = 4.5, color = "grey40", fontface = "italic") +

  theme_void() +
  labs(title = "Modèle PD Emax — compartiment unique")

# =============================================================================
# PANNEAU 2 : Courbe Emax
# =============================================================================

C_seq  <- seq(0, 5, length.out = 300)
Emax   <- 1          # effet maximum normalisé
EC50   <- 1          # concentration pour 50% d'effet

df_emax <- data.frame(
  C = C_seq,
  E = Emax * C_seq / (EC50 + C_seq)
)

p_emax <- ggplot(df_emax, aes(x = C, y = E)) +
  geom_line(color = "#E74C3C", linewidth = 1.5) +

  # Ligne asymptote Emax
  geom_hline(yintercept = Emax, linetype = "dashed",
             color = "grey50", linewidth = 0.8) +

  # Point EC50
  geom_segment(x = EC50, xend = EC50, y = 0, yend = 0.5,
               linetype = "dotted", color = "#E67E22", linewidth = 1) +
  geom_segment(x = 0, xend = EC50, y = 0.5, yend = 0.5,
               linetype = "dotted", color = "#E67E22", linewidth = 1) +
  geom_point(x = EC50, y = 0.5, size = 4, color = "#E67E22") +

  # Annotations
  annotate("text", x = 4.5, y = 1.05, label = "Emax",
           size = 4, color = "grey40", fontface = "bold") +
  annotate("text", x = EC50, y = -0.07, label = "EC50",
           size = 4, color = "#E67E22", fontface = "bold") +
  annotate("text", x = -0.15, y = 0.5, label = "Emax/2",
           size = 3.5, color = "#E67E22", hjust = 1) +

  # Zone faible dose (linéaire)
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = Emax * 0.5 / (EC50 + 0.5),
           fill = "#ABEBC6", alpha = 0.3) +
  annotate("text", x = 0.25, y = 0.08, label = "≈ linéaire\n(C << EC50)",
           size = 2.8, color = "#1E8449", hjust = 0.5) +

  # Zone saturation
  annotate("rect", xmin = 3, xmax = 5, ymin = Emax * 3/(EC50+3), ymax = Emax,
           fill = "#F1948A", alpha = 0.3) +
  annotate("text", x = 4, y = 0.88, label = "Saturation\n(C >> EC50)",
           size = 2.8, color = "#922B21", hjust = 0.5) +

  scale_x_continuous(expand = c(0.02, 0),
                     labels = function(x) ifelse(x == 0, "0", paste0(x, "·EC50"))) +
  scale_y_continuous(expand = c(0.02, 0), limits = c(-0.12, 1.15),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c("0", "0.25·Emax", "0.5·Emax", "0.75·Emax", "Emax")) +
  labs(
    title    = "Fonction E(C) — Effet Emax (Hill n=1)",
    subtitle = expression(E(C) == frac(E[max] %.% C, EC[50] + C)),
    x        = "Concentration C",
    y        = "Effet E(C)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.subtitle = element_text(size = 13, color = "grey30"),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# PANNEAU 3 : Simulation TV selon différents régimes
# =============================================================================

sim_emax <- function(kg, Emax_val, EC50_val, C_const, TV0 = 150, t_end = 50) {
  dt <- 0.1
  t  <- seq(0, t_end, by = dt)
  TV <- numeric(length(t))
  TV[1] <- TV0
  for (i in seq_along(t)[-1]) {
    E  <- Emax_val * C_const / (EC50_val + C_const)
    TV[i] <- TV[i-1] + dt * (kg - E) * TV[i-1]
    if (TV[i] < 0) TV[i] <- 0
  }
  data.frame(t = t, TV = TV)
}

kg_val    <- 0.05
Emax_val  <- 0.15
EC50_val  <- 50

df_ctrl    <- sim_emax(kg_val, Emax_val, EC50_val, C_const = 0)
df_tsc     <- sim_emax(kg_val, Emax_val, EC50_val, C_const = kg_val * EC50_val / (Emax_val - kg_val))
df_inhib   <- sim_emax(kg_val, Emax_val, EC50_val, C_const = 200)
df_regress <- sim_emax(kg_val, Emax_val, EC50_val, C_const = 1000)

df_sim3 <- rbind(
  cbind(df_ctrl,    Groupe = "Contrôle (C = 0)"),
  cbind(df_tsc,     Groupe = "Stase (C = TSC)"),
  cbind(df_inhib,   Groupe = "Inhibition partielle"),
  cbind(df_regress, Groupe = "Régression (C >> TSC)")
)

df_sim3$Groupe <- factor(df_sim3$Groupe,
  levels = c("Contrôle (C = 0)", "Inhibition partielle",
             "Stase (C = TSC)", "Régression (C >> TSC)"))

cols3 <- c("Contrôle (C = 0)"      = "#888888",
           "Inhibition partielle"   = "#4393C3",
           "Stase (C = TSC)"        = "#F39C12",
           "Régression (C >> TSC)"  = "#E74C3C")

p_sim <- ggplot(df_sim3, aes(x = t, y = TV, color = Groupe)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 150, linetype = "dotted", color = "grey60") +
  annotate("text", x = 48, y = 160, label = "TV₀", size = 3.5, color = "grey50") +
  scale_color_manual(values = cols3) +
  labs(
    title    = "Dynamique tumorale selon la concentration",
    subtitle = "kg fixé, Emax et EC50 fixés — C constante (simulation illustrative)",
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm³)",
    color    = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 9, color = "grey50")
  )

# =============================================================================
# ASSEMBLAGE
# =============================================================================

layout <- rbind(c(1, 2),
                c(3, 3))

g <- arrangeGrob(p_schema, p_emax, p_sim, layout_matrix = layout,
                 top = textGrob(
                   "Modèle PD Emax — fonctionnement",
                   gp = gpar(fontsize = 15, fontface = "bold")
                 ))

ggsave("scripts/schema_PD_Emax.png", g,
       width = 13, height = 10, dpi = 150)
cat("Schéma → scripts/schema_PD_Emax.png\n")
