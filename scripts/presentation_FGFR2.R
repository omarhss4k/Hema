# =============================================================================
# PRÉSENTATION TUTRICE — Fc-silent FGFR2-huBPA-LP1
# Génère 9 figures PNG prêtes à intégrer dans une présentation :
#   01_donnees_PK_brutes.png
#   02_PK_predit_vs_observe.png
#   03_parametres_PK.png
#   04_schema_modele_PK.png
#   05_croissance_vs_inhibition.png
#   06_schema_modele_Simeoni.png
#   07_PKPD_predit_vs_observe.png
#   08_tableau_TGI.png
#   09_parametres_PKPD.png
# =============================================================================

library(rxode2)
library(deSolve)
library(DEoptim)
library(ggplot2)
library(readxl)
library(patchwork)

OUT <- "scripts/presentation"
dir.create(OUT, showWarnings = FALSE)

SAVE <- function(name, p, w = 10, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h,
         dpi = 180, bg = "white")
  cat("  →", file.path(OUT, paste0(name, ".png")), "\n")
}

COLS <- c("10 mg/kg" = "#1B4F9E",
          "5 mg/kg"  = "#2E86C1",
          "1 mg/kg"  = "#85C1E9",
          "Contrôle" = "#888888",
          "3 mg/kg"  = "#27AE60")

THEME <- theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        legend.background = element_rect(fill = "white", color = NA))

# =============================================================================
# CHARGEMENT DES RÉSULTATS PRÉ-CALCULÉS
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # pk2comp_rxode2, df_pk
load("scripts/resultats_PKPD_rxode2_FGFR2.RData")       # k1_fixed, k2_global,
                                                          # params_fixed, w0,
                                                          # dat_ctrl, dat_d3, dat_d10

PK <- c(CL = unname(pk2comp_rxode2$CL),
        V1 = unname(pk2comp_rxode2$V1),
        V2 = unname(pk2comp_rxode2$V2),
        Q  = unname(pk2comp_rxode2$Q))

# Doses PK
dose10_pk <- 10 * 1000
dose5_pk  <-  5 * 1000
dose1_pk  <-  1 * 1000

# Doses PD
dose10_pd <- 10 * 1000
dose3_pd  <-  3 * 1000
dose0_pd  <-  0

T_CENSOR_CTRL <- 35 * 24
DOSE_DAYS     <- c(0, 14, 28, 42)
PD_P          <- 1

# =============================================================================
# MODÈLES (redéfinition locale, pas de dépendance globale)
# =============================================================================

mod2comp <- rxode2({
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
  d/dt(A2) <-  (Q/V1) * A1 - (Q/V2) * A2
})

sim_pk <- function(dose, params, times) {
  ev <- eventTable()
  ev$add.dosing(dose = dose, nbr.doses = 1, dosing.to = 1)
  ev$add.sampling(sort(unique(c(0, times))))
  out <- tryCatch(rxSolve(mod2comp, params, ev), error = function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, length(times)))
  approx(out$time, out$C1, xout = times, rule = 2)$y
}

pkpd_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    C1     <- A1 / V1
    w      <- x1 + x2 + x3 + x4
    growth <- l0 * x1 / (1 + (l0/l1 * w)^p)^(1/p)
    list(c(
      -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2,
       (Q/V1)*A1 - (Q/V2)*A2,
      growth - k2*C1*x1,
      k2*C1*x1  - k1*x2,
      k1*(x2 - x3),
      k1*(x3 - x4)
    ))
  })
}

sim_multi <- function(dose, params, obs_times, n_doses = 4, interval_h = 14*24) {
  state0    <- c(A1 = dose, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)
  events_df <- if (n_doses > 1)
    data.frame(var    = "A1",
               time   = seq(interval_h, (n_doses-1)*interval_h, by = interval_h),
               value  = dose, method = "add")
  else NULL
  times <- sort(unique(c(0, obs_times)))
  out <- tryCatch(
    as.data.frame(lsoda(state0, times, pkpd_ode, as.list(params),
                        events = if (!is.null(events_df)) list(data = events_df) else NULL,
                        rtol = 1e-6, atol = 1e-8)),
    error = function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, length(obs_times)))
  w_vec <- with(out, x1 + x2 + x3 + x4)
  approx(out$time, pmax(w_vec, 1e-9), xout = obs_times, rule = 2)$y
}

par_pkpd <- c(params_fixed, k1 = k1_fixed, k2 = k2_global)

# =============================================================================
# FIG 1 — DONNÉES PK BRUTES
# =============================================================================
cat("\n[1/9] Données PK brutes\n")

raw_pk <- read_xlsx("PK souris FGFR2.xlsx")
time_h <- as.numeric(raw_pk[[1]])

.mk <- function(t, c, lbl) {
  ok <- !is.na(c) & c > 0
  data.frame(Temps_h = t[ok], Concentration = c[ok], Dose = lbl)
}
df_pk_raw <- rbind(
  .mk(time_h, as.numeric(raw_pk[[2]]), "10 mg/kg"),
  .mk(time_h, as.numeric(raw_pk[[3]]), "5 mg/kg"),
  .mk(time_h, as.numeric(raw_pk[[4]]), "1 mg/kg")
)
df_pk_raw$Dose <- factor(df_pk_raw$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))

p1a <- ggplot(df_pk_raw, aes(x = Temps_h, y = Concentration, color = Dose)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(aes(group = Dose), linewidth = 0.6, linetype = "dotted") +
  scale_color_manual(values = COLS) +
  labs(title = "Données PK observées — Échelle linéaire",
       x = "Temps (heures)", y = "Concentration (µg/L)", color = NULL) +
  THEME

p1b <- ggplot(df_pk_raw, aes(x = Temps_h, y = Concentration, color = Dose)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(aes(group = Dose), linewidth = 0.6, linetype = "dotted") +
  scale_y_log10() +
  scale_color_manual(values = COLS) +
  labs(title = "Données PK observées — Échelle log",
       x = "Temps (heures)", y = "Concentration log(µg/L)", color = NULL) +
  THEME

SAVE("01_donnees_PK_brutes", p1a + p1b + plot_layout(guides = "collect") &
       theme(legend.position = "bottom"), w = 13, h = 5)

# =============================================================================
# FIG 2 — AJUSTEMENT PK PRÉDIT VS OBSERVÉ
# =============================================================================
cat("[2/9] Ajustement PK prédit vs observé\n")

t_full_pk <- seq(0, max(df_pk_raw$Temps_h) * 1.05, by = 1)
df_pk_sim <- do.call(rbind, list(
  data.frame(Temps_h = t_full_pk,
             Concentration = sim_pk(dose10_pk, PK, t_full_pk), Dose = "10 mg/kg"),
  data.frame(Temps_h = t_full_pk,
             Concentration = sim_pk(dose5_pk,  PK, t_full_pk), Dose = "5 mg/kg"),
  data.frame(Temps_h = t_full_pk,
             Concentration = sim_pk(dose1_pk,  PK, t_full_pk), Dose = "1 mg/kg")
))
df_pk_sim$Dose <- factor(df_pk_sim$Dose, levels = c("10 mg/kg","5 mg/kg","1 mg/kg"))

p2 <- ggplot() +
  geom_line(data = df_pk_sim, aes(x = Temps_h, y = Concentration, color = Dose),
            linewidth = 1.2) +
  geom_point(data = df_pk_raw, aes(x = Temps_h, y = Concentration, color = Dose),
             size = 3) +
  scale_y_log10() +
  scale_color_manual(values = COLS) +
  labs(
    title    = "Modèle PK 2-compartiments — Prédit vs Observé",
    subtitle = paste0(
      "V1 = ", round(PK["V1"], 3), " L/kg   ",
      "V2 = ", round(PK["V2"], 3), " L/kg   ",
      "CL = ", round(PK["CL"]*24, 3), " L/j/kg   ",
      "t½β = ", round(pk2comp_rxode2$t_half_beta/24, 1), " j"
    ),
    x = "Temps (heures)", y = "Concentration (µg/L — échelle log)",
    color = NULL,
    caption = "Lignes = modèle   ●  Points = données observées"
  ) +
  THEME

SAVE("02_PK_predit_vs_observe", p2, w = 9, h = 5.5)

# =============================================================================
# FIG 3 — TABLEAU PARAMÈTRES PK
# =============================================================================
cat("[3/9] Tableau paramètres PK\n")

pk_tbl <- data.frame(
  Paramètre   = c("CL",  "V1",  "V2",  "Vss", "Q",   "k10", "k12", "k21",
                  "α",   "β",   "t½α", "t½β"),
  Description = c("Clairance systémique",
                  "Volume central",
                  "Volume périphérique",
                  "Volume distribution ss",
                  "Clairance intercompart.",
                  "Constante élimination",
                  "Transfert C1 → C2",
                  "Transfert C2 → C1",
                  "Constante phase rapide",
                  "Constante phase lente",
                  "Demi-vie distribution",
                  "Demi-vie élimination"),
  Valeur      = c(
    sprintf("%.4f L/h/kg  (%.3f L/j/kg)", PK["CL"], PK["CL"]*24),
    sprintf("%.4f L/kg",  PK["V1"]),
    sprintf("%.4f L/kg",  PK["V2"]),
    sprintf("%.4f L/kg",  pk2comp_rxode2$Vss),
    sprintf("%.4f L/h/kg (%.3f L/j/kg)", PK["Q"], PK["Q"]*24),
    sprintf("%.5f /h",    pk2comp_rxode2$k10),
    sprintf("%.5f /h",    pk2comp_rxode2$k12),
    sprintf("%.5f /h",    pk2comp_rxode2$k21),
    sprintf("%.5f /h",    pk2comp_rxode2$alpha),
    sprintf("%.6f /h",    pk2comp_rxode2$beta),
    sprintf("%.1f h",     pk2comp_rxode2$t_half_alpha),
    sprintf("%.1f h  (%.1f jours)", pk2comp_rxode2$t_half_beta,
            pk2comp_rxode2$t_half_beta/24)
  )
)

# Construction grille ggplot
nr <- nrow(pk_tbl)
nc <- 3
cells3 <- do.call(rbind, lapply(seq_len(nr), function(i) {
  bg <- if (i %% 2 == 0) "#EBF5FB" else "white"
  hl <- i > 9   # demi-vies en surbrillance
  data.frame(row  = i, col = 1:nc,
             lbl  = c(pk_tbl$Paramètre[i], pk_tbl$Description[i], pk_tbl$Valeur[i]),
             fill = c(if (hl) "#D6EAF8" else bg),
             bold = hl, stringsAsFactors = FALSE)
}))
hdr3 <- data.frame(row = 0, col = 1:nc,
                   lbl  = c("Paramètre", "Description", "Valeur"),
                   fill = "#1A5276", bold = TRUE)
cells3 <- rbind(hdr3, cells3)
cells3$y    <- max(cells3$row) - cells3$row
cells3$tcol <- ifelse(cells3$row == 0, "white", "black")
cells3$face <- ifelse(cells3$bold, "bold", "plain")
xw <- c(1, 1, 2.2)
cells3$xpos <- xw[cells3$col]

p3 <- ggplot(cells3, aes(x = xpos, y = y)) +
  geom_tile(aes(fill = fill, width = c(0.9,0.9,1.9)[cells3$col]),
            color = "grey80", linewidth = 0.3) +
  geom_text(aes(label = lbl, color = tcol, fontface = face),
            size = 3.3, hjust = 0.5) +
  scale_fill_identity() + scale_color_identity() +
  labs(title    = "Paramètres PK — Modèle 2-compartiments IV bolus",
       subtitle = "Fc-silent FGFR2-huBPA-LP1 | Souris | Optimisation nlminb (log-espace)") +
  theme_void(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13, margin = margin(b=4)),
        plot.subtitle = element_text(size = 9, color = "grey40", margin = margin(b=8)),
        plot.margin   = margin(12, 20, 12, 20))

SAVE("03_parametres_PK", p3, w = 11, h = 5.5)

# =============================================================================
# FIG 4 — SCHÉMA MODÈLE PK 2-COMPARTIMENTS
# =============================================================================
cat("[4/9] Schéma modèle PK 2-compartiments\n")

arrow <- function(x1,y1,x2,y2)
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2),
               arrow = arrow(length = unit(0.3,"cm"), type="closed"),
               linewidth = 1, color = "#2C3E50")

box <- function(xmin,xmax,ymin,ymax,fill="#D6EAF8",col="#2980B9")
  annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
           fill=fill, color=col, linewidth=1.2)

lbl <- function(x,y,txt,sz=4.5,bold=FALSE,col="black")
  annotate("text", x=x, y=y, label=txt, size=sz,
           fontface=if(bold)"bold"else"plain", color=col, hjust=0.5)

p4 <- ggplot() + xlim(0,10) + ylim(0,6) +
  # Dose IV bolus
  box(0.2, 1.8, 3.8, 5.2, fill="#FDEBD0", col="#E67E22") +
  lbl(1, 4.8, "Dose IV bolus", sz=4, bold=TRUE, col="#E67E22") +
  lbl(1, 4.3, "D = dose × kg", sz=3.5, col="#7F8C8D") +

  # Compartiment central
  box(2.8, 6.2, 3.2, 5.8, fill="#D6EAF8", col="#2980B9") +
  lbl(4.5, 5.1, "Compartiment central", sz=4.5, bold=TRUE, col="#1A5276") +
  lbl(4.5, 4.6, "A1  (quantité, µg/kg)", sz=3.8, col="#2C3E50") +
  lbl(4.5, 4.0, "C1 = A1 / V1", sz=3.8, col="#2C3E50") +

  # Compartiment périphérique
  box(2.8, 6.2, 0.5, 2.8, fill="#D5F5E3", col="#27AE60") +
  lbl(4.5, 2.1, "Compartiment périphérique", sz=4.5, bold=TRUE, col="#1E8449") +
  lbl(4.5, 1.6, "A2  (quantité, µg/kg)", sz=3.8, col="#2C3E50") +

  # Élimination
  box(7.2, 9.8, 3.8, 5.2, fill="#FADBD8", col="#E74C3C") +
  lbl(8.5, 4.8, "Élimination", sz=4, bold=TRUE, col="#C0392B") +
  lbl(8.5, 4.3, "CL · C1", sz=3.5, col="#7F8C8D") +

  # Flèches
  arrow(1.8, 4.5, 2.8, 4.5) +   # dose → C1
  arrow(6.2, 4.5, 7.2, 4.5) +   # C1 → élim

  # Q : C1 ↔ C2
  arrow(4.5, 3.2, 4.5, 2.8) +   # C1 → C2
  arrow(4.9, 2.8, 4.9, 3.2) +   # C2 → C1

  # Labels flèches Q
  lbl(3.9, 3.0, "Q/V1", sz=3.5, col="#27AE60") +
  lbl(5.5, 3.0, "Q/V2", sz=3.5, col="#27AE60") +

  # Équations
  annotate("label", x=4.5, y=0.25,
           label = "dA1/dt = -(CL/V1 + Q/V1)·A1 + (Q/V2)·A2\ndA2/dt =  (Q/V1)·A1 - (Q/V2)·A2",
           size = 3.5, hjust = 0.5, fill = "#F8F9FA", color = "#2C3E50",
           label.padding = unit(0.4,"lines")) +

  labs(title    = "Structure du modèle PK — 2 compartiments (IV bolus)",
       subtitle = "Paramètres estimés : CL, V1, V2, Q  |  Optimisation nlminb en log-espace") +
  theme_void(base_size = 12) +
  theme(plot.title    = element_text(face="bold", size=14, margin=margin(b=4)),
        plot.subtitle = element_text(size=9, color="grey40", margin=margin(b=6)),
        plot.margin   = margin(10,10,10,10))

SAVE("04_schema_modele_PK", p4, w = 11, h = 6)

# =============================================================================
# FIG 5 — CROISSANCE TUMORALE VS INHIBITION
# =============================================================================
cat("[5/9] Croissance tumorale vs inhibition\n")

read_tumor_excel <- function(path) {
  raw <- suppressMessages(read_xlsx(path, col_names=FALSE, .name_repair="minimal"))
  raw <- raw[, colSums(!is.na(raw)) > 0]
  hidx <- which(raw[[1]] == "Group")
  times <- suppressWarnings(as.numeric(as.character(unlist(raw[hidx[1], -1]))))
  valid <- !is.na(times);  times <- times[valid]
  extr  <- function(s, e) {
    b <- raw[s:e, ];  b <- b[!is.na(b[[1]]) & b[[1]] != "", ]
    dat <- b[, c(TRUE, valid)]
    for (j in 2:ncol(dat)) dat[[j]] <- as.numeric(dat[[j]])
    dat
  }
  list(times=times,
       blk_mean=extr(hidx[1]+1, hidx[2]-1),
       blk_sem =extr(hidx[2]+1, nrow(raw)))
}
.gr <- function(blk, pat)
  as.numeric(unlist(blk[grep(pat, blk[[1]], ignore.case=TRUE)[1], -1]))

tv <- read_tumor_excel("TumorVolume_FGFR2.xlsx")
times_d <- tv$times

w_ctrl   <- .gr(tv$blk_mean, "Group 01")
w_d10_tv <- .gr(tv$blk_mean, "Group 03")
w_d3_tv  <- .gr(tv$blk_mean, "Group 04")
sem_ctrl <- .gr(tv$blk_sem,  "Group 01")
sem_d10  <- .gr(tv$blk_sem,  "Group 03")
sem_d3   <- .gr(tv$blk_sem,  "Group 04")

df_tv <- rbind(
  data.frame(jour=times_d, mean=w_ctrl,   sem=sem_ctrl, Groupe="Contrôle"),
  data.frame(jour=times_d, mean=w_d10_tv, sem=sem_d10,  Groupe="10 mg/kg"),
  data.frame(jour=times_d, mean=w_d3_tv,  sem=sem_d3,   Groupe="3 mg/kg")
)
df_tv <- df_tv[!is.na(df_tv$mean), ]
df_tv$Groupe <- factor(df_tv$Groupe, levels=c("Contrôle","3 mg/kg","10 mg/kg"))

ymax <- max(df_tv$mean + df_tv$sem, na.rm=TRUE)

p5 <- ggplot(df_tv, aes(x=jour, y=mean, color=Groupe, group=Groupe)) +
  geom_vline(xintercept=DOSE_DAYS, linetype="dashed", color="grey70", linewidth=0.4) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem, fill=Groupe),
              alpha=0.15, color=NA) +
  geom_line(linewidth=1.2) +
  geom_point(size=2.8) +
  annotate("point", x=DOSE_DAYS, y=-ymax*0.07,
           shape=17, size=4, color="#C0392B") +
  annotate("text", x=max(DOSE_DAYS)+1.5, y=-ymax*0.07,
           label="= Traitement Q2W", hjust=0, size=3.2, color="#C0392B") +
  scale_color_manual(values=COLS) +
  scale_fill_manual(values=COLS) +
  scale_x_continuous(breaks=seq(0, max(df_tv$jour), by=7)) +
  coord_cartesian(ylim=c(-ymax*0.14, ymax*1.1), clip="off") +
  labs(
    title    = "Croissance tumorale vs Inhibition par le traitement",
    subtitle = "Fc-silent FGFR2-huBPA-LP1 | mean ± SEM | IV Q2W × 4",
    x = "Temps (jours)", y = "Volume tumoral (mm³)", color=NULL, fill=NULL
  ) + THEME +
  theme(legend.position="bottom", plot.margin=margin(5,10,30,5))

SAVE("05_croissance_vs_inhibition", p5, w=10, h=6)

# =============================================================================
# FIG 6 — SCHÉMA MODÈLE SIMEONI PKPD
# =============================================================================
cat("[6/9] Schéma modèle Simeoni PKPD\n")

p6 <- ggplot() + xlim(0,14) + ylim(0,8) +

  # ---- Bloc PK ----
  annotate("rect", xmin=0.2, xmax=3.8, ymin=5.5, ymax=7.5,
           fill="#D6EAF8", color="#2980B9", linewidth=1.2) +
  lbl(2, 6.9, "PK 2-compartiments", sz=4.2, bold=TRUE, col="#1A5276") +
  lbl(2, 6.4, "C1(t) = A1(t) / V1",  sz=3.6, col="#2C3E50") +
  lbl(2, 5.9, "k12, k21, CL, V1, V2", sz=3.4, col="#7F8C8D") +

  # flèche PK → PD
  arrow(3.8, 6.5, 5.2, 4.7) +
  lbl(4.8, 5.8, "C1(t)", sz=3.5, col="#E74C3C") +

  # ---- Cellules proliférantes x1 ----
  annotate("rect", xmin=5.0, xmax=8.0, ymin=3.5, ymax=5.5,
           fill="#D5F5E3", color="#27AE60", linewidth=1.2) +
  lbl(6.5, 5.0, "x1", sz=5, bold=TRUE, col="#1E8449") +
  lbl(6.5, 4.5, "Prolifération", sz=3.8, col="#2C3E50") +
  lbl(6.5, 4.0, "l0·x1 / (1+(l0/l1·w)^p)^(1/p)", sz=2.9, col="#7F8C8D") +
  lbl(6.5, 3.7, "−k2·C1·x1", sz=3.1, col="#E74C3C") +

  # flèche x1 → x2
  arrow(8.0, 4.5, 9.0, 4.5) +
  lbl(8.5, 4.8, "k2·C1", sz=3.2, col="#E74C3C") +

  # ---- Transit x2 ----
  annotate("rect", xmin=9.0, xmax=10.5, ymin=3.5, ymax=5.5,
           fill="#FDEBD0", color="#E67E22", linewidth=1.0) +
  lbl(9.75, 4.5, "x2", sz=4.5, bold=TRUE, col="#784212") +
  lbl(9.75, 3.8, "k1·x2", sz=3, col="#7F8C8D") +

  arrow(10.5, 4.5, 11.0, 4.5) +

  # ---- Transit x3 ----
  annotate("rect", xmin=11.0, xmax=12.0, ymin=3.5, ymax=5.5,
           fill="#FDEBD0", color="#E67E22", linewidth=1.0) +
  lbl(11.5, 4.5, "x3", sz=4.5, bold=TRUE, col="#784212") +

  arrow(12.0, 4.5, 12.5, 4.5) +

  # ---- Transit x4 ----
  annotate("rect", xmin=12.5, xmax=13.5, ymin=3.5, ymax=5.5,
           fill="#FADBD8", color="#E74C3C", linewidth=1.0) +
  lbl(13.0, 4.5, "x4", sz=4.5, bold=TRUE, col="#922B21") +

  # ---- Masse tumorale totale ----
  annotate("rect", xmin=5.5, xmax=12.8, ymin=0.3, ymax=1.8,
           fill="#F9F9F9", color="#7F8C8D", linewidth=0.8, linetype="dashed") +
  lbl(9.2, 1.4, "w(t) = x1 + x2 + x3 + x4", sz=4, bold=TRUE, col="#2C3E50") +
  lbl(9.2, 0.8, "Volume tumoral total simulé (comparé aux données)", sz=3.3, col="#7F8C8D") +

  # flèche w → x1 (rétroaction croissance)
  geom_curve(aes(x=9.0, y=1.8, xend=6.5, yend=3.5),
             arrow=arrow(length=unit(0.28,"cm"), type="closed"),
             curvature=-0.3, linewidth=0.8, color="#27AE60") +
  lbl(7.3, 2.8, "rétroaction\ncroissance", sz=3, col="#1E8449") +

  # ---- Légende équations ----
  annotate("label", x=2, y=2.5,
           label=paste0(
             "Paramètres PD :\n",
             "  l0  = taux croissance exponentielle (/h)\n",
             "  l1  = taux croissance linéaire (g/h)\n",
             "  k1  = transit cellules endommagées (/h)\n",
             "  k2  = activité cytotoxique [L/(µg·h)]\n",
             "  p   = coeff. Hill (p=1, Simeoni 2004)"
           ),
           size=3.2, hjust=0, fill="#F4F6F7", color="#2C3E50",
           label.padding=unit(0.5,"lines")) +

  labs(title    = "Structure du modèle PKPD — Simeoni 2004",
       subtitle = "Croissance tumorale + compartiments de transit pour les cellules endommagées") +
  theme_void(base_size=12) +
  theme(plot.title    = element_text(face="bold", size=14, margin=margin(b=4)),
        plot.subtitle = element_text(size=9, color="grey40", margin=margin(b=6)),
        plot.margin   = margin(10,10,10,10))

SAVE("06_schema_modele_Simeoni", p6, w=13, h=7.5)

# =============================================================================
# FIG 7 — PKPD PRÉDIT VS OBSERVÉ
# =============================================================================
cat("[7/9] PKPD prédit vs observé\n")

times_full_pd <- seq(0, 49*24, by=4)
lev_pd <- c("Contrôle","3 mg/kg","10 mg/kg")

df_sim_pd <- do.call(rbind, list(
  data.frame(t=times_full_pd/24,
             w=sim_multi(dose0_pd,  par_pkpd, times_full_pd), Groupe="Contrôle"),
  data.frame(t=times_full_pd/24,
             w=sim_multi(dose3_pd,  par_pkpd, times_full_pd), Groupe="3 mg/kg"),
  data.frame(t=times_full_pd/24,
             w=sim_multi(dose10_pd, par_pkpd, times_full_pd), Groupe="10 mg/kg")
))
df_sim_pd$Groupe <- factor(df_sim_pd$Groupe, levels=lev_pd)

df_obs_pd <- rbind(
  data.frame(t=dat_ctrl$t/24, w=dat_ctrl$w, Groupe="Contrôle"),
  data.frame(t=dat_d3$t/24,   w=dat_d3$w,   Groupe="3 mg/kg"),
  data.frame(t=dat_d10$t/24,  w=dat_d10$w,  Groupe="10 mg/kg")
)
df_obs_pd$Groupe <- factor(df_obs_pd$Groupe, levels=lev_pd)

# Résidus relatifs
df_res <- do.call(rbind, lapply(lev_pd, function(g) {
  obs <- df_obs_pd[df_obs_pd$Groupe == g, ]
  t_obs <- obs$t
  w_sim <- sim_multi(
    switch(g, "Contrôle"=dose0_pd, "3 mg/kg"=dose3_pd, "10 mg/kg"=dose10_pd),
    par_pkpd, t_obs*24)
  data.frame(t=t_obs, resid=(obs$w - w_sim)/w_sim*100, Groupe=g)
}))
df_res$Groupe <- factor(df_res$Groupe, levels=lev_pd)

p7a <- ggplot() +
  geom_vline(xintercept=DOSE_DAYS, linetype="dashed", color="grey70", linewidth=0.4) +
  geom_line(data=df_sim_pd, aes(x=t, y=w, color=Groupe), linewidth=1.2) +
  geom_point(data=df_obs_pd, aes(x=t, y=w, color=Groupe, shape=Groupe), size=3) +
  scale_color_manual(values=COLS) +
  scale_shape_manual(values=c(16,17,15)) +
  labs(title="PKPD Simeoni 2004 — Prédit vs Observé",
       subtitle=paste0("k1=", round(k1_fixed,4), " /h  |  k2=",
                       formatC(k2_global, format="e", digits=2),
                       "  |  Protocole Q2W × 4 doses"),
       x="Temps (jours)", y="Volume tumoral (g)",
       color=NULL, shape=NULL,
       caption="Lignes = modèle   ●▲■ = données observées") +
  THEME + theme(legend.position="bottom")

p7b <- ggplot(df_res, aes(x=t, y=resid, color=Groupe)) +
  geom_hline(yintercept=0, linewidth=0.8, color="grey40") +
  geom_hline(yintercept=c(-30,30), linetype="dashed", color="grey70") +
  geom_point(size=2.5, alpha=0.9) +
  scale_color_manual(values=COLS) +
  labs(title="Résidus relatifs",
       subtitle="(Observé − Prédit) / Prédit × 100",
       x="Temps (jours)", y="Résidu relatif (%)", color=NULL) +
  THEME + theme(legend.position="none")

SAVE("07_PKPD_predit_vs_observe",
     p7a + p7b + plot_layout(widths=c(2,1)), w=14, h=6)

# =============================================================================
# FIG 8 — TABLEAU TGI OBSERVÉ + PRÉDIT
# =============================================================================
cat("[8/9] Tableau TGI\n")

w0_ctrl <- w_ctrl[times_d == 0]
w0_d10  <- w_d10_tv[times_d == 0]
w0_d3   <- w_d3_tv[times_d == 0]
ok_tgi  <- !is.na(w_ctrl) & !is.na(w_d10_tv) & !is.na(w_d3_tv) & times_d > 0

tgi_f <- function(wt, wt0, wc, wc0) {
  dc <- wc - wc0
  if (is.na(dc) || dc <= 0) return(NA_real_)
  round((1 - (wt - wt0)/dc)*100, 1)
}

# TGI prédit
t_seq_d <- seq(0, 49, by=1) * 24
wc_sim  <- sim_multi(dose0_pd,  par_pkpd, t_seq_d)
w3_sim  <- sim_multi(dose3_pd,  par_pkpd, t_seq_d)
w10_sim <- sim_multi(dose10_pd, par_pkpd, t_seq_d)
t_match <- times_d[ok_tgi]
wc_s  <- approx(t_seq_d/24, wc_sim,  xout=t_match)$y
w3_s  <- approx(t_seq_d/24, w3_sim,  xout=t_match)$y
w10_s <- approx(t_seq_d/24, w10_sim, xout=t_match)$y
w0_g  <- w0 * 1000   # g → mm³

df_tgi_full <- data.frame(
  Temps      = paste0("j", times_d[ok_tgi]),
  TGI_3_obs  = mapply(tgi_f, w_d3_tv[ok_tgi],  w0_d3,  w_ctrl[ok_tgi], w0_ctrl),
  TGI_10_obs = mapply(tgi_f, w_d10_tv[ok_tgi], w0_d10, w_ctrl[ok_tgi], w0_ctrl),
  TGI_3_pred = mapply(function(w3, wc)
    tgi_f(w3*1000, w0_g, wc*1000, w0_g), w3_s, wc_s),
  TGI_10_pred= mapply(function(w10, wc)
    tgi_f(w10*1000, w0_g, wc*1000, w0_g), w10_s, wc_s)
)

tgi_hdr <- c("Temps",
             "TGI 3mg obs (%)", "TGI 3mg pred (%)",
             "TGI 10mg obs (%)", "TGI 10mg pred (%)")
nc8 <- length(tgi_hdr)
nr8 <- nrow(df_tgi_full)

.fill_tgi <- function(v, bg) {
  if (is.na(v))   return(bg)
  if (v >= 60)   return("#A9DFBF")
  if (v >= 30)   return("#FAD7A0")
  bg
}

rows8 <- lapply(seq_len(nr8), function(i) {
  bg  <- if (i %% 2 == 0) "#F2F3F4" else "white"
  r   <- df_tgi_full[i, ]
  lbs <- c(r$Temps,
            sprintf("%.1f %%", r$TGI_3_obs),
            sprintf("%.1f %%", r$TGI_3_pred),
            sprintf("%.1f %%", r$TGI_10_obs),
            sprintf("%.1f %%", r$TGI_10_pred))
  fls <- c(bg,
            .fill_tgi(r$TGI_3_obs, bg), .fill_tgi(r$TGI_3_pred, bg),
            .fill_tgi(r$TGI_10_obs, bg), .fill_tgi(r$TGI_10_pred, bg))
  data.frame(row=i, col=1:nc8, lbl=lbs, fill=fls, bold=(i==nr8))
})

hdr8 <- data.frame(row=0, col=1:nc8, lbl=tgi_hdr,
                   fill="#1A5276", bold=TRUE)
cells8 <- rbind(hdr8, do.call(rbind, rows8))
cells8$y    <- max(cells8$row) - cells8$row
cells8$tcol <- ifelse(cells8$row==0, "white", "black")
cells8$face <- ifelse(cells8$bold, "bold", "plain")

p8 <- ggplot(cells8, aes(x=col, y=y)) +
  geom_tile(aes(fill=fill), color="grey70", linewidth=0.3) +
  geom_text(aes(label=lbl, color=tcol, fontface=face),
            size=3.2, lineheight=0.95) +
  scale_fill_identity() + scale_color_identity() +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.05,0)) +
  labs(
    title    = "Tumour Growth Inhibition — Observé vs Prédit",
    subtitle = paste0("Vert : TGI ≥ 60 %  ·  Orange : 30–60 %  ·  Dernier point en gras\n",
                      "TGI = (1 − ΔW_traité / ΔW_contrôle) × 100")
  ) +
  theme_void(base_size=11) +
  theme(plot.title    = element_text(face="bold", size=13, margin=margin(b=4)),
        plot.subtitle = element_text(size=9, color="grey40", margin=margin(b=6)),
        plot.margin   = margin(12,12,12,12))

SAVE("08_tableau_TGI", p8,
     w=11, h=0.5*(nr8+2)+2)

# =============================================================================
# FIG 9 — TABLEAU RÉCAPITULATIF PARAMÈTRES PKPD
# =============================================================================
cat("[9/9] Tableau paramètres PKPD\n")

pkpd_par_tbl <- data.frame(
  Paramètre   = c("l0",  "l1",  "p",  "k1",   "k2",
                  "V1",  "V2",  "CL", "Q",    "w0"),
  Type        = c("PD","PD","PD","PD","PD",
                  "PK (fixé)","PK (fixé)","PK (fixé)","PK (fixé)","Initial"),
  Description = c("Taux croissance exponentielle",
                  "Taux croissance linéaire",
                  "Coeff. Hill (Simeoni 2004)",
                  "Transition cellules endommagées",
                  "Activité cytotoxique",
                  "Volume compartiment central",
                  "Volume périphérique",
                  "Clairance systémique",
                  "Clairance intercompart.",
                  "Masse tumorale initiale"),
  Valeur      = c(
    sprintf("%.4f /j  (%.6f /h)", params_fixed["l0"]*24, params_fixed["l0"]),
    sprintf("%.4f g/j (%.6f g/h)", params_fixed["l1"]*24, params_fixed["l1"]),
    sprintf("%.0f", params_fixed["p"]),
    sprintf("%.4f /j  (%.6f /h)", k1_fixed*24, k1_fixed),
    formatC(k2_global, format="e", digits=3),
    sprintf("%.4f L/kg", PK["V1"]),
    sprintf("%.4f L/kg", PK["V2"]),
    sprintf("%.4f L/h/kg", PK["CL"]),
    sprintf("%.4f L/h/kg", PK["Q"]),
    sprintf("%.4f g", w0)
  ),
  Estimation  = c("DEoptim","DEoptim","Fixé (littérature)","Fixé (Simeoni 2004)",
                  "DEoptim","nlminb (PK)","nlminb (PK)","nlminb (PK)","nlminb (PK)",
                  "Moyenne groupes j0")
)

type_fill <- c("PD"="#D5F5E3", "PK (fixé)"="#D6EAF8", "Initial"="#FEF9E7")

nr9 <- nrow(pkpd_par_tbl)
nc9 <- 5
rows9 <- lapply(seq_len(nr9), function(i) {
  bg <- type_fill[pkpd_par_tbl$Type[i]]
  data.frame(row=i, col=1:nc9,
             lbl=c(pkpd_par_tbl$Paramètre[i], pkpd_par_tbl$Type[i],
                   pkpd_par_tbl$Description[i], pkpd_par_tbl$Valeur[i],
                   pkpd_par_tbl$Estimation[i]),
             fill=bg, bold=FALSE)
})
hdr9 <- data.frame(row=0, col=1:nc9,
                   lbl=c("Paramètre","Type","Description","Valeur","Méthode d'estimation"),
                   fill="#1A5276", bold=TRUE)
cells9 <- rbind(hdr9, do.call(rbind, rows9))
cells9$y    <- max(cells9$row) - cells9$row
cells9$tcol <- ifelse(cells9$row==0, "white", "black")
cells9$face <- ifelse(cells9$bold, "bold", "plain")

p9 <- ggplot(cells9, aes(x=col, y=y)) +
  geom_tile(aes(fill=fill), color="grey70", linewidth=0.3) +
  geom_text(aes(label=lbl, color=tcol, fontface=face),
            size=3.2, lineheight=0.95) +
  scale_fill_identity() + scale_color_identity() +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.05,0)) +
  labs(
    title    = "Récapitulatif des paramètres du modèle PKPD",
    subtitle = "Vert = paramètres PD estimés  ·  Bleu = paramètres PK fixés  ·  Jaune = conditions initiales"
  ) +
  theme_void(base_size=11) +
  theme(plot.title    = element_text(face="bold", size=13, margin=margin(b=4)),
        plot.subtitle = element_text(size=9, color="grey40", margin=margin(b=6)),
        plot.margin   = margin(12,12,12,12))

SAVE("09_parametres_PKPD", p9, w=13, h=5.5)

# =============================================================================
cat("\n=== Toutes les figures générées dans", OUT, "===\n")
list.files(OUT, pattern="\\.png$")
