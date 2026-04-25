# =============================================================================
# Schéma du modèle PK 2-compartiments — Fc-silent B/C huBPA-LP1
# =============================================================================

library(ggplot2)

# Charger les paramètres estimés
load("scripts/resultats_PK2comp_FGFR2.RData")

# Arrondir les paramètres pour affichage
V1  <- round(pk2comp$V1,  4)
V2  <- round(pk2comp$V2,  4)
CL  <- round(pk2comp$CL * 24, 4)   # L/j/kg
Q   <- round(pk2comp$Q  * 24, 4)   # L/j/kg
k10 <- round(pk2comp$k10, 5)
k12 <- round(pk2comp$k12, 5)
k21 <- round(pk2comp$k21, 5)

# ── Coordonnées des éléments ──────────────────────────────────────────────────
# Dose (flèche entrante)
# Compartiment central (C1)
# Compartiment périphérique (C2)
# Élimination (flèche sortante)

box_h  <- 1.4    # hauteur des boîtes
box_w  <- 2.8    # largeur des boîtes
gap    <- 2.2    # espace entre les boîtes

# Centres
x_c1 <- 3.0;   y_c1 <- 3.0
x_c2 <- 3.0 + box_w + gap;  y_c2 <- 3.0

# Flèche dose
x_dose <- x_c1 - 2.5
y_dose <- y_c1 + box_h/2 + 0.5

# Couleurs
col_box   <- "#2E86AB"
col_elim  <- "#E84855"
col_dose  <- "#3BB273"
col_arrow <- "#444444"

ggplot() +

  # ── Compartiment central ──────────────────────────────────────────────────
  annotate("rect",
           xmin = x_c1 - box_w/2, xmax = x_c1 + box_w/2,
           ymin = y_c1 - box_h/2, ymax = y_c1 + box_h/2,
           fill = col_box, color = "white", linewidth = 1.2, alpha = 0.9) +
  annotate("text", x = x_c1, y = y_c1 + 0.25,
           label = "Compartiment central",
           color = "white", fontface = "bold", size = 4.2) +
  annotate("text", x = x_c1, y = y_c1 - 0.28,
           label = paste0("V1 = ", V1, " L/kg"),
           color = "white", size = 3.8) +

  # ── Compartiment périphérique ─────────────────────────────────────────────
  annotate("rect",
           xmin = x_c2 - box_w/2, xmax = x_c2 + box_w/2,
           ymin = y_c2 - box_h/2, ymax = y_c2 + box_h/2,
           fill = "#5C6BC0", color = "white", linewidth = 1.2, alpha = 0.9) +
  annotate("text", x = x_c2, y = y_c2 + 0.25,
           label = "Compartiment\npériphérique",
           color = "white", fontface = "bold", size = 4.2, lineheight = 0.9) +
  annotate("text", x = x_c2, y = y_c2 - 0.28,
           label = paste0("V2 = ", V2, " L/kg"),
           color = "white", size = 3.8) +

  # ── Flèche k12 : C1 → C2 ─────────────────────────────────────────────────
  annotate("segment",
           x    = x_c1 + box_w/2,
           xend = x_c2 - box_w/2,
           y    = y_c1 + 0.2,
           yend = y_c2 + 0.2,
           arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
           color = col_arrow, linewidth = 1) +
  annotate("text",
           x = (x_c1 + box_w/2 + x_c2 - box_w/2) / 2,
           y = y_c1 + 0.55,
           label = paste0("k12 = ", k12, " /h\n(Q = ", Q, " L/j/kg)"),
           size = 3.3, color = col_arrow, fontface = "italic") +

  # ── Flèche k21 : C2 → C1 ─────────────────────────────────────────────────
  annotate("segment",
           x    = x_c2 - box_w/2,
           xend = x_c1 + box_w/2,
           y    = y_c1 - 0.2,
           yend = y_c2 - 0.2,
           arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
           color = col_arrow, linewidth = 1) +
  annotate("text",
           x = (x_c1 + box_w/2 + x_c2 - box_w/2) / 2,
           y = y_c1 - 0.6,
           label = paste0("k21 = ", k21, " /h"),
           size = 3.3, color = col_arrow, fontface = "italic") +

  # ── Flèche dose IV ────────────────────────────────────────────────────────
  annotate("segment",
           x = x_c1 - box_w/2 - 1.5, xend = x_c1 - box_w/2,
           y = y_c1, yend = y_c1,
           arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
           color = col_dose, linewidth = 1.3) +
  annotate("label",
           x = x_c1 - box_w/2 - 0.8, y = y_c1 + 0.45,
           label = "Dose IV\n(bolus t=0)",
           fill = col_dose, color = "white",
           size = 3.5, fontface = "bold", label.size = 0) +

  # ── Flèche élimination k10 ────────────────────────────────────────────────
  annotate("segment",
           x = x_c1, xend = x_c1,
           y = y_c1 - box_h/2, yend = y_c1 - box_h/2 - 1.4,
           arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
           color = col_elim, linewidth = 1.3) +
  annotate("label",
           x = x_c1 + 1.1, y = y_c1 - box_h/2 - 0.75,
           label = paste0("k10 = ", k10, " /h\nCL = ", CL, " L/j/kg"),
           fill = col_elim, color = "white",
           size = 3.5, fontface = "bold", label.size = 0) +

  # ── Équation du modèle ────────────────────────────────────────────────────
  annotate("text",
           x = (x_c1 + x_c2) / 2, y = 0.8,
           label = "C(t) = A·exp(-α·t) + B·exp(-β·t)",
           size = 4, color = "grey30", fontface = "italic") +
  annotate("text",
           x = (x_c1 + x_c2) / 2, y = 0.3,
           label = paste0("t½α = ", round(pk2comp$t_half_alpha, 1),
                          " h     t½β = ", round(pk2comp$t_half_beta, 1),
                          " h (≈ ", round(pk2comp$t_half_beta/24, 1), " j)"),
           size = 3.8, color = "grey30") +

  # ── Mise en page ──────────────────────────────────────────────────────────
  xlim(0.2, 10.8) + ylim(-0.2, 5.2) +
  labs(title = "Modèle PK 2-compartiments — Fc-silent B/C huBPA-LP1 (FGFR2)") +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(hjust=0.5, face="bold", size=13, margin=margin(b=10)),
    plot.margin   = margin(15, 15, 15, 15),
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("scripts/schema_PK2comp_FGFR2.png", width=10, height=5, dpi=150)
cat("Schéma sauvegardé → scripts/schema_PK2comp_FGFR2.png\n")
