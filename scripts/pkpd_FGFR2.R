# =============================================================================
# PK/PD — Fc-silent B/C huBPA-LP1 (FGFR2)
# 1. PK 2-compartiments  : rxode2 + nlminb
# 2. Courbe tumorale     : lecture Excel + TGI
# 3. PKPD Simeoni 2004   : deSolve/lsoda + DEoptim
# =============================================================================

library(rxode2)
library(deSolve)
library(DEoptim)
library(ggplot2)
library(readxl)

# =============================================================================
# CONSTANTES GLOBALES
# =============================================================================

T_CENSOR_CTRL_J <- 35        # jours au-delà desquels le contrôle est censuré
                              # (souris euthanasées → moyenne chute)
T_CENSOR_CTRL   <- T_CENSOR_CTRL_J * 24   # en heures

PD_P            <- 1         # coefficient de Hill Simeoni 2004 (original = 1)
                              # p=20 provoque des gradients extrêmes → instabilité

DOSE_DAYS       <- c(0, 14, 28, 42)   # protocole Q2W × 4

# =============================================================================
# 1. DONNÉES PK
# =============================================================================

raw    <- read_xlsx("PK souris FGFR2.xlsx")
time_h <- as.numeric(raw[[1]])
conc10 <- as.numeric(raw[[2]])
conc5  <- as.numeric(raw[[3]])
conc1  <- as.numeric(raw[[4]])

dose10_pk <- 10 * 1000   # µg/kg
dose5_pk  <-  5 * 1000
dose1_pk  <-  1 * 1000

.mk_pk <- function(t, c, d, label) {
  ok <- !is.na(c) & c > 0
  data.frame(t = t[ok], C = c[ok], dose = d, Dose = label)
}
df_pk <- rbind(
  .mk_pk(time_h, conc10, dose10_pk, "10 mg/kg"),
  .mk_pk(time_h, conc5,  dose5_pk,  "5 mg/kg"),
  .mk_pk(time_h, conc1,  dose1_pk,  "1 mg/kg")
)
df_pk$Dose <- factor(df_pk$Dose, levels = c("10 mg/kg", "5 mg/kg", "1 mg/kg"))

d10 <- df_pk[df_pk$Dose == "10 mg/kg", ]
d5  <- df_pk[df_pk$Dose == "5 mg/kg",  ]
d1  <- df_pk[df_pk$Dose == "1 mg/kg",  ]

# =============================================================================
# 2. MODÈLE PK rxode2 — 2 compartiments, IV bolus
# =============================================================================

mod2comp <- rxode2({
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
  d/dt(A2) <-  (Q/V1) * A1 - (Q/V2) * A2
})

sim_dose_pk <- function(dose, params, times) {
  ev <- eventTable()
  ev$add.dosing(dose = dose, nbr.doses = 1, dosing.to = 1)
  ev$add.sampling(sort(unique(c(0, times))))
  out <- tryCatch(rxSolve(mod2comp, params, ev), error = function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, length(times)))
  approx(out$time, out$C1, xout = times, rule = 2)$y
}

# =============================================================================
# 3. VALEURS INITIALES PK — méthode des résidus
# =============================================================================

load("scripts/resultats_PK2comp_FGFR2.RData")   # → pk2comp

init_params <- c(CL = pk2comp$CL, V1 = pk2comp$V1,
                 V2 = pk2comp$V2, Q  = pk2comp$Q)

cat("=== Valeurs initiales (méthode des résidus) ===\n")
cat("CL =", round(init_params["CL"], 6), "L/h/kg\n")
cat("V1 =", round(init_params["V1"], 5), "L/kg\n")
cat("V2 =", round(init_params["V2"], 5), "L/kg\n")
cat("Q  =", round(init_params["Q"],  6), "L/h/kg\n")

# =============================================================================
# 4. FONCTION OBJECTIVE PK — résidus log, poids égaux par groupe
# =============================================================================

objective_pk <- function(logpar) {
  par        <- exp(logpar)
  names(par) <- c("CL", "V1", "V2", "Q")

  p10 <- tryCatch(sim_dose_pk(dose10_pk, par, d10$t), error = function(e) NULL)
  p5  <- tryCatch(sim_dose_pk(dose5_pk,  par, d5$t),  error = function(e) NULL)
  p1  <- tryCatch(sim_dose_pk(dose1_pk,  par, d1$t),  error = function(e) NULL)

  if (is.null(p10) || any(is.na(p10) | p10 <= 0)) return(1e10)
  if (is.null(p5)  || any(is.na(p5)  | p5  <= 0)) return(1e10)
  if (is.null(p1)  || any(is.na(p1)  | p1  <= 0)) return(1e10)

  mean((log(d10$C) - log(p10))^2) +
  mean((log(d5$C)  - log(p5))^2)  +
  mean((log(d1$C)  - log(p1))^2)
}

# =============================================================================
# 5. OPTIMISATION PK — nlminb en log-espace
# =============================================================================

cat("\nOptimisation nlminb en cours...\n")

fit_pk <- nlminb(
  start     = log(init_params),
  objective = objective_pk,
  control   = list(eval.max = 3000, iter.max = 1500,
                   rel.tol = 1e-12, x.tol = 1e-12)
)

best_par        <- exp(fit_pk$par)
names(best_par) <- c("CL", "V1", "V2", "Q")

# =============================================================================
# 6. PARAMÈTRES DÉRIVÉS PK
# =============================================================================

k10 <- best_par["CL"] / best_par["V1"]
k12 <- best_par["Q"]  / best_par["V1"]
k21 <- best_par["Q"]  / best_par["V2"]
Vss <- best_par["V1"] + best_par["V2"]

# Valeurs propres — forme équivalente à sqrt((Σk)²−4·k10·k21)
disc         <- sqrt((k10 + k12 - k21)^2 + 4 * k12 * k21)
sum_k        <- k10 + k12 + k21
alpha        <- (sum_k + disc) / 2
beta         <- (sum_k - disc) / 2
t_half_alpha <- log(2) / alpha
t_half_beta  <- log(2) / beta

cat("\n=== Paramètres PK 2-compartiments ===\n")
cat("CL  =", round(best_par["CL"], 6), "L/h/kg  =", round(best_par["CL"]*24, 4), "L/j/kg\n")
cat("V1  =", round(best_par["V1"], 5), "L/kg\n")
cat("V2  =", round(best_par["V2"], 5), "L/kg\n")
cat("Vss =", round(Vss,            5), "L/kg\n")
cat("Q   =", round(best_par["Q"],  6), "L/h/kg  =", round(best_par["Q"]*24, 4), "L/j/kg\n")
cat("k10 =", round(k10, 6), "/h\n")
cat("k12 =", round(k12, 6), "/h\n")
cat("k21 =", round(k21, 6), "/h\n")
cat("t½α =", round(t_half_alpha, 2), "h\n")
cat("t½β =", round(t_half_beta,  1), "h =", round(t_half_beta/24, 2), "jours\n")
cat("Objectif :", round(fit_pk$objective, 6), "\n")

# =============================================================================
# 7. GRAPHIQUE PK
# =============================================================================

times_pk_full <- seq(0, max(df_pk$t) * 1.05, by = 1)

df_sim_pk <- do.call(rbind, lapply(
  list(list(dose10_pk, "10 mg/kg"),
       list(dose5_pk,  "5 mg/kg"),
       list(dose1_pk,  "1 mg/kg")),
  function(x) data.frame(
    t    = times_pk_full,
    C    = sim_dose_pk(x[[1]], best_par, times_pk_full),
    Dose = x[[2]]
  )
))
df_sim_pk$Dose <- factor(df_sim_pk$Dose, levels = c("10 mg/kg", "5 mg/kg", "1 mg/kg"))

ggplot() +
  geom_line(data = df_sim_pk, aes(x = t, y = C, color = Dose), linewidth = 1) +
  geom_point(data = df_pk,    aes(x = t, y = C, color = Dose), size = 2.5) +
  scale_y_log10() +
  labs(
    title    = "PK 2-compartiments (rxode2) — Fc-silent B/C huBPA-LP1 (FGFR2)",
    subtitle = paste0(
      "V1 = ", round(best_par["V1"], 4), " L/kg",
      "   V2 = ", round(best_par["V2"], 4), " L/kg",
      "   CL = ", round(best_par["CL"]*24, 4), " L/j/kg",
      "   t½α = ", round(t_half_alpha, 1), " h",
      "   t½β = ", round(t_half_beta/24, 1), " j"
    ),
    x = "Temps (heures)", y = "Concentration (échelle log)"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_PK2comp_rxode2_FGFR2.png", width = 8, height = 5, dpi = 150)
cat("\nGraphique → scripts/plot_PK2comp_rxode2_FGFR2.png\n")

# =============================================================================
# 8. SAUVEGARDE PK
# =============================================================================

pk2comp_rxode2 <- list(
  CL = best_par["CL"], V1 = best_par["V1"],
  V2 = best_par["V2"], Q  = best_par["Q"],
  k10 = k10, k12 = k12, k21 = k21, Vss = Vss,
  alpha = alpha, beta = beta,
  t_half_alpha = t_half_alpha, t_half_beta = t_half_beta
)
save(pk2comp_rxode2, df_pk,
     file = "scripts/resultats_PK2comp_rxode2_FGFR2.RData")
cat("Résultats → scripts/resultats_PK2comp_rxode2_FGFR2.RData\n")

# =============================================================================
# HELPER : lecture robuste du fichier TumorVolume (deux blocs means + SEM)
# Utilisé par les sections courbe tumorale, TGI et PKPD.
# =============================================================================

read_tumor_excel <- function(path) {
  raw <- suppressMessages(
    read_xlsx(path, col_names = FALSE, .name_repair = "minimal")
  )
  raw <- raw[, colSums(!is.na(raw)) > 0]

  header_idx <- which(raw[[1]] == "Group")
  stopifnot("Deux blocs 'Group' attendus (means + SEM)" = length(header_idx) == 2)

  times <- suppressWarnings(
    as.numeric(as.character(unlist(raw[header_idx[1], -1])))
  )
  valid <- !is.na(times)
  times <- times[valid]

  extract_block <- function(start, end) {
    blk <- raw[start:end, ]
    blk <- blk[!is.na(blk[[1]]) & blk[[1]] != "", ]
    dat <- blk[, c(TRUE, valid)]
    for (j in seq(2, ncol(dat))) dat[[j]] <- as.numeric(dat[[j]])
    dat
  }

  list(
    times    = times,
    blk_mean = extract_block(header_idx[1] + 1, header_idx[2] - 1),
    blk_sem  = extract_block(header_idx[2] + 1, nrow(raw))
  )
}

.get_row <- function(blk, pattern)
  as.numeric(unlist(blk[grep(pattern, blk[[1]], ignore.case = TRUE)[1], -1]))

# =============================================================================
# 9. COURBE DE CROISSANCE TUMORALE
# =============================================================================

tv <- read_tumor_excel("TumorVolume_FGFR2.xlsx")

cat("\nGroupes (moyennes) :", tv$blk_mean[[1]], "\n")
cat("Temps (jours)      :", tv$times, "\n")

patterns_tv <- c("Group 01", "Group 03", "Group 04")
labels_tv   <- c("Contrôle", "10 mg/kg", "3 mg/kg")

df_tv <- do.call(rbind, lapply(seq_along(patterns_tv), function(i) {
  mn  <- .get_row(tv$blk_mean, patterns_tv[i])
  sem <- .get_row(tv$blk_sem,  patterns_tv[i])
  ok  <- !is.na(mn)
  data.frame(jour = tv$times[ok], mean = mn[ok], sem = sem[ok],
             Groupe = labels_tv[i])
}))
df_tv$Groupe <- factor(df_tv$Groupe, levels = c("Contrôle", "3 mg/kg", "10 mg/kg"))

ymax <- max(df_tv$mean + df_tv$sem, na.rm = TRUE)
cols_tv <- c("Contrôle" = "#888888", "3 mg/kg" = "#4393C3", "10 mg/kg" = "#2166AC")

p_tv <- ggplot(df_tv, aes(x = jour, y = mean, color = Groupe, group = Groupe)) +
  geom_vline(xintercept = DOSE_DAYS, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.8, linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  annotate("point", x = DOSE_DAYS, y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text", x = max(DOSE_DAYS) + 1.5, y = -ymax * 0.06,
           label = "= Traitement", hjust = 0, size = 3.2, color = "#CC0000") +
  scale_x_continuous(breaks = seq(0, max(df_tv$jour), by = 7)) +
  scale_color_manual(values = cols_tv) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "Fc-silent FGFR2-huBPA-LP1 — Courbe de croissance tumorale",
    subtitle = "mean ± SEM",
    x = "Temps (jours)", y = "Volume tumoral (mm³)", color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 10, color = "grey50"),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

ggsave("scripts/plot_courbe_tumorale_FGFR2.png", p_tv,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_courbe_tumorale_FGFR2.png\n")

# SEM pour pondération PKPD ultérieure
dat_sem <- list(
  time = tv$times,
  ctrl = .get_row(tv$blk_sem, "Group 01"),
  d3   = .get_row(tv$blk_sem, "Group 04"),
  d10  = .get_row(tv$blk_sem, "Group 03")
)
save(dat_sem, file = "scripts/sem_tumorale_FGFR2.RData")
cat("SEM sauvegardées → scripts/sem_tumorale_FGFR2.RData\n")

# =============================================================================
# 10. TABLEAU TGI OBSERVÉ
# =============================================================================

times_d  <- tv$times
w_ctrl   <- .get_row(tv$blk_mean, "Group 01")
w_d10_tv <- .get_row(tv$blk_mean, "Group 03")
w_d3_tv  <- .get_row(tv$blk_mean, "Group 04")
sem_ctrl <- .get_row(tv$blk_sem,  "Group 01")
sem_d10  <- .get_row(tv$blk_sem,  "Group 03")
sem_d3   <- .get_row(tv$blk_sem,  "Group 04")

w0_ctrl_tv <- w_ctrl[times_d == 0]
w0_d10_tv  <- w_d10_tv[times_d == 0]
w0_d3_tv   <- w_d3_tv[times_d == 0]
ok_tgi     <- !is.na(w_ctrl) & !is.na(w_d10_tv) & !is.na(w_d3_tv) & times_d > 0

tgi_pct <- function(wt, wt0, wc, wc0)
  round((1 - (wt - wt0) / (wc - wc0)) * 100, 1)

df_tgi <- data.frame(
  jour      = times_d[ok_tgi],
  tgi_3     = tgi_pct(w_d3_tv[ok_tgi],  w0_d3_tv,  w_ctrl[ok_tgi], w0_ctrl_tv),
  tgi_10    = tgi_pct(w_d10_tv[ok_tgi], w0_d10_tv, w_ctrl[ok_tgi], w0_ctrl_tv),
  mean_ctrl = round(w_ctrl[ok_tgi],    0),
  sem_ctrl  = round(sem_ctrl[ok_tgi],  0),
  mean_3    = round(w_d3_tv[ok_tgi],   0),
  sem_3     = round(sem_d3[ok_tgi],    0),
  mean_10   = round(w_d10_tv[ok_tgi],  0),
  sem_10    = round(sem_d10[ok_tgi],   0)
)

cat(sprintf("\n%-6s  %12s  %14s  %12s  %14s\n",
            "Temps", "Ctrl (mm³)", "3 mg/kg (mm³)", "TGI 3mg%", "TGI 10mg%"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(df_tgi))) {
  cat(sprintf("j%-5d  %5.0f±%-5.0f  %5.0f±%-6.0f  %8.1f %%  %10.1f %%\n",
              df_tgi$jour[i],
              df_tgi$mean_ctrl[i], df_tgi$sem_ctrl[i],
              df_tgi$mean_3[i],    df_tgi$sem_3[i],
              df_tgi$tgi_3[i],     df_tgi$tgi_10[i]))
}
cat(strrep("-", 65), "\n")
cat(sprintf("Dernier point (j%d)  →  3 mg/kg : %.1f %%  |  10 mg/kg : %.1f %%\n\n",
            tail(df_tgi$jour, 1), tail(df_tgi$tgi_3, 1), tail(df_tgi$tgi_10, 1)))

# Tableau visuel ggplot2
col_names_tbl <- c("Temps",
                   "Contrôle\nmoy (mm³)", "Contrôle\nSEM",
                   "3 mg/kg\nmoy (mm³)",  "3 mg/kg\nSEM",  "TGI 3mg\n(%)",
                   "10 mg/kg\nmoy (mm³)", "10 mg/kg\nSEM", "TGI 10mg\n(%)")
n_col <- length(col_names_tbl)
n_row <- nrow(df_tgi)

header_tbl <- data.frame(
  row = 0, col = seq_len(n_col), label = col_names_tbl,
  fill = "#2c3e50", tcol = "white", bold = TRUE
)

rows_tbl <- lapply(seq_len(n_row), function(i) {
  tgi3  <- df_tgi$tgi_3[i];  tgi10 <- df_tgi$tgi_10[i]
  bg    <- if (i %% 2 == 0) "#f2f2f2" else "white"
  fill_tgi3  <- if (tgi3  >= 60) "#c3e6cb" else if (tgi3  >= 30) "#ffeeba" else bg
  fill_tgi10 <- if (tgi10 >= 60) "#c3e6cb" else if (tgi10 >= 30) "#ffeeba" else bg
  data.frame(
    row   = i, col = seq_len(n_col),
    label = c(paste0("j", df_tgi$jour[i]),
              df_tgi$mean_ctrl[i], df_tgi$sem_ctrl[i],
              df_tgi$mean_3[i],    df_tgi$sem_3[i],   sprintf("%.1f %%", tgi3),
              df_tgi$mean_10[i],   df_tgi$sem_10[i],  sprintf("%.1f %%", tgi10)),
    fill  = c(rep(bg, 5), fill_tgi3, rep(bg, 2), fill_tgi10),
    tcol  = "black",
    bold  = (i == n_row)
  )
})

cells <- do.call(rbind, c(list(header_tbl), rows_tbl))
cells$y    <- max(cells$row) - cells$row
cells$face <- ifelse(cells$bold, "bold", "plain")

p_tbl <- ggplot(cells, aes(x = col, y = y)) +
  geom_tile(aes(fill = fill), color = "grey70", linewidth = 0.3) +
  geom_text(aes(label = label, color = tcol, fontface = face),
            size = 3.0, lineheight = 0.95) +
  scale_fill_identity() + scale_color_identity() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.05, 0)) +
  labs(
    title    = "TGI observé — Fc-silent FGFR2-huBPA-LP1",
    subtitle = "Q2W (j0/j14/j28/j42) | Vert : TGI ≥ 60 %  ·  Jaune : 30–60 %  ·  Dernier point en gras"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey40", margin = margin(b = 6)),
    plot.margin   = margin(12, 12, 12, 12)
  )

ggsave("scripts/tableau_TGI_FGFR2.png", p_tbl,
       width = 11, height = 0.55 * (n_row + 2) + 1.5, dpi = 150)
cat("Tableau → scripts/tableau_TGI_FGFR2.png\n")

write.csv(
  data.frame(
    Temps        = paste0("j", df_tgi$jour),
    Ctrl_moy_mm3 = df_tgi$mean_ctrl, Ctrl_SEM     = df_tgi$sem_ctrl,
    D3_moy_mm3   = df_tgi$mean_3,    D3_SEM       = df_tgi$sem_3,
    TGI_3mg_pct  = df_tgi$tgi_3,
    D10_moy_mm3  = df_tgi$mean_10,   D10_SEM      = df_tgi$sem_10,
    TGI_10mg_pct = df_tgi$tgi_10
  ),
  "scripts/tableau_TGI_FGFR2.csv", row.names = FALSE
)
cat("CSV → scripts/tableau_TGI_FGFR2.csv\n")

# =============================================================================
# 11. PARAMÈTRES PK FIXÉS POUR LE PKPD
# =============================================================================

load("resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2

PK <- c(
  CL = unname(pk2comp_rxode2$CL),
  V1 = unname(pk2comp_rxode2$V1),
  V2 = unname(pk2comp_rxode2$V2),
  Q  = unname(pk2comp_rxode2$Q)
)

cat("\n=== Paramètres PK fixés ===\n")
cat("CL =", round(PK["CL"], 6), "L/h/kg\n")
cat("V1 =", round(PK["V1"], 5), "L/kg\n")
cat("V2 =", round(PK["V2"], 5), "L/kg\n")
cat("Q  =", round(PK["Q"],  6), "L/h/kg\n")

C1_max_10 <- 10000 / PK["V1"]   # µg/L — repère pour les bornes k2

# =============================================================================
# 12. DONNÉES TUMORALES PKPD (mm³ → g)
# =============================================================================

# Réutilisation de la fonction read_tumor_excel définie en section helper
tv_pd <- read_tumor_excel("TumorVolume_FGFR2.xlsx")

time_d_pd <- tv_pd$times
time_h_pd <- time_d_pd * 24

.mk_pd <- function(t, w) {
  ok <- !is.na(w) & w > 0
  data.frame(t = t[ok], w = w[ok] / 1000)   # mm³ → g
}

dat_ctrl <- .mk_pd(time_h_pd, .get_row(tv_pd$blk_mean, "Group 01"))
dat_d10  <- .mk_pd(time_h_pd, .get_row(tv_pd$blk_mean, "Group 03"))
dat_d3   <- .mk_pd(time_h_pd, .get_row(tv_pd$blk_mean, "Group 04"))

dose10_pd <- 10 * 1000   # µg/kg
dose3_pd  <-  3 * 1000
dose0_pd  <-  0

w0 <- mean(c(
  dat_ctrl$w[which.min(dat_ctrl$t)],
  dat_d10$w[which.min(dat_d10$t)],
  dat_d3$w[which.min(dat_d3$t)]
), na.rm = TRUE)
cat("\nw0 :", round(w0, 4), "g\n")

# =============================================================================
# 13. MODÈLES PKPD — Simeoni 2004
#
#  Modèle A (linéaire)  : effet = k2 * C1               — utilisé pour l0/l1
#  Modèle B (Emax)      : effet = Emax * C1 / (EC50+C1) — modèle principal
#
#  Le modèle Emax capture la saturation de l'effet à forte dose :
#   - À faible C1 : effet ≈ (Emax/EC50) * C1  (linéaire, comme k2)
#   - À forte C1 : effet → Emax               (plafond biologique)
#  Cela explique pourquoi un k2 unique ne peut pas ajuster 3 et 10 mg/kg.
# =============================================================================

# Modèle A — linéaire (k2), utilisé uniquement pour la calibration l0/l1
pkpd_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    C1     <- A1 / V1
    w      <- x1 + x2 + x3 + x4
    growth <- l0 * x1 / (1 + (l0/l1 * w)^p)^(1/p)

    dA1 <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
    dA2 <-  (Q/V1) * A1 - (Q/V2) * A2
    dx1 <- growth - k2 * C1 * x1
    dx2 <- k2 * C1 * x1 - k1 * x2
    dx3 <- k1 * (x2 - x3)
    dx4 <- k1 * (x3 - x4)

    list(c(dA1, dA2, dx1, dx2, dx3, dx4))
  })
}

# Modèle B — Emax (modèle principal)
pkpd_ode_emax <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    C1     <- A1 / V1
    w      <- x1 + x2 + x3 + x4
    growth <- l0 * x1 / (1 + (l0/l1 * w)^p)^(1/p)
    kill   <- Emax * C1 / (EC50 + C1)   # effet saturant [/h]

    dA1 <- -(CL/V1 + Q/V1) * A1 + (Q/V2) * A2
    dA2 <-  (Q/V1) * A1 - (Q/V2) * A2
    dx1 <- growth - kill * x1
    dx2 <- kill * x1 - k1 * x2
    dx3 <- k1 * (x2 - x3)
    dx4 <- k1 * (x3 - x4)

    list(c(dA1, dA2, dx1, dx2, dx3, dx4))
  })
}

# =============================================================================
# 14. PRÉ-CALIBRATION l0, l1 DEPUIS LE GROUPE CONTRÔLE
# =============================================================================

obj_growth <- function(logpar) {
  l0_t <- unname(exp(logpar[1]));  l1_t <- unname(exp(logpar[2]))
  ok   <- dat_ctrl$t <= T_CENSOR_CTRL
  pars   <- c(as.list(PK), l0 = l0_t, l1 = l1_t, p = PD_P, k1 = 1, k2 = 0)
  state0 <- c(A1 = 0, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)
  out <- tryCatch(
    as.data.frame(lsoda(state0, sort(unique(c(0, dat_ctrl$t[ok]))),
                        pkpd_ode, pars, rtol = 1e-6, atol = 1e-8)),
    error = function(e) NULL
  )
  if (is.null(out)) return(1e10)
  w_vec <- with(out, x1 + x2 + x3 + x4)
  if (any(!is.finite(w_vec))) return(1e10)
  pc  <- approx(out$time, pmax(w_vec, 1e-9), xout = dat_ctrl$t[ok], rule = 2)$y
  val <- mean((log(dat_ctrl$w[ok]) - log(pc))^2)
  if (!is.finite(val)) 1e10 else val
}

cat("\nPré-calibration l0, l1 (DEoptim)...\n")
fit_growth <- DEoptim(
  obj_growth,
  lower   = c(log(0.001/24), log(0.005/24)),
  upper   = c(log(0.5/24),   log(5/24)),
  control = DEoptim.control(NP = 20L, itermax = 500L,
                            F = 0.8, CR = 0.9, strategy = 2L, trace = FALSE)
)

l0 <- unname(exp(fit_growth$optim$bestmem[1]))
l1 <- unname(exp(fit_growth$optim$bestmem[2]))

cat("  l0 =", round(l0*24, 4), "/j  (", round(l0, 7), "/h)\n")
cat("  l1 =", round(l1*24, 4), "g/j (", round(l1, 6), "g/h)\n")
cat("  Objectif contrôle :", round(fit_growth$optim$bestval, 5), "\n")

k2_ref <- l0 / C1_max_10
cat("  Repère k2 : k2_ref ≈", formatC(k2_ref, format = "e", digits = 2), "\n")

# =============================================================================
# 15. FONCTIONS DE SIMULATION PKPD — Modèle Emax
# =============================================================================

.PENALTY <- 1e6

.sim_emax <- function(state0, pars_list, obs_times, events_df = NULL) {
  times <- sort(unique(c(0, obs_times)))
  out <- tryCatch(
    as.data.frame(lsoda(state0, times, pkpd_ode_emax, pars_list,
                        events = if (!is.null(events_df)) list(data = events_df) else NULL,
                        rtol = 1e-6, atol = 1e-8)),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(.PENALTY, length(obs_times)))
  w_vec <- with(out, x1 + x2 + x3 + x4)
  if (any(!is.finite(w_vec)) || max(out$time) < max(obs_times) - 1)
    return(rep(.PENALTY, length(obs_times)))
  approx(out$time, pmax(w_vec, 1e-9), xout = obs_times, rule = 2)$y
}

simulate_single <- function(dose, params, obs_times) {
  state0 <- c(A1 = dose, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)
  .sim_emax(state0, as.list(params), obs_times)
}

simulate_multi <- function(dose, params, obs_times,
                           n_doses = 4, interval_h = 14 * 24) {
  state0    <- c(A1 = dose, A2 = 0, x1 = w0, x2 = 0, x3 = 0, x4 = 0)
  events_df <- if (n_doses > 1)
    data.frame(var    = "A1",
               time   = seq(interval_h, (n_doses - 1) * interval_h, by = interval_h),
               value  = dose, method = "add")
  else NULL
  .sim_emax(state0, as.list(params), obs_times, events_df)
}

# =============================================================================
# 16. FONCTION OBJECTIVE — Modèle Emax
#     Paramètres estimés : Emax (/h), EC50 (µg/L)
#     k1 fixé (Simeoni 2004), l0/l1/p fixés (calibration contrôle)
# =============================================================================

params_fixed <- c(PK, l0 = l0, l1 = l1, p = PD_P)

# k1 fixé — transit moyen = 4/k1 = 8 jours (Simeoni 2004)
k1_fixed <- 0.5 / 24   # /h

.residuals_emax <- function(par) {
  ok_c <- dat_ctrl$t <= T_CENSOR_CTRL
  pc  <- tryCatch(simulate_multi(dose0_pd,  par, dat_ctrl$t), error = function(e) NULL)
  p3  <- tryCatch(simulate_multi(dose3_pd,  par, dat_d3$t),   error = function(e) NULL)
  p10 <- tryCatch(simulate_multi(dose10_pd, par, dat_d10$t),  error = function(e) NULL)

  .bad <- function(p) is.null(p) || !length(p) || any(!is.finite(p) | p <= 0)
  if (.bad(pc) || .bad(p3) || .bad(p10)) return(1e10)

  val <- mean((log(dat_ctrl$w[ok_c]) - log(pc[ok_c]))^2) +
         mean((log(dat_d3$w)         - log(p3))^2)        +
         mean((log(dat_d10$w)        - log(p10))^2)

  if (!is.finite(val)) 1e10 else val
}

obj_emax <- function(logpar) {
  Emax_t <- unname(exp(logpar[1]))
  EC50_t <- unname(exp(logpar[2]))
  # Garde-fou numérique : Emax ne doit pas dépasser 1000·l0
  if (Emax_t > 1000 * l0) return(1e10)
  par <- c(params_fixed, k1 = k1_fixed, Emax = Emax_t, EC50 = EC50_t)
  .residuals_emax(par)
}

# =============================================================================
# 17. OPTIMISATION EMAX — DEoptim (log-espace, 2 paramètres)
#
#   Emax  : taux de destruction maximal [/h]
#           bornes : [l0 * 0.01,  l0 * 500]  — de très faible à très fort effet
#   EC50  : concentration à effet semi-maximal [µg/L]
#           bornes : [1,  C1_max_10 * 100]   — très sensible à très résistant
# =============================================================================

cat("\n=== Repères pour les bornes Emax ===\n")
cat("  l0         =", formatC(l0,         format="e", digits=2), "/h\n")
cat("  C1_max_10  =", formatC(C1_max_10,  format="e", digits=2), "µg/L\n")
cat("  Emax_lower =", formatC(l0*0.01,    format="e", digits=2), "/h\n")
cat("  Emax_upper =", formatC(l0*500,     format="e", digits=2), "/h\n")
cat("  EC50_lower = 1 µg/L\n")
cat("  EC50_upper =", formatC(C1_max_10*100, format="e", digits=2), "µg/L\n")

run_deoptim <- function(obj_fn, lower, upper, label, NP = NULL) {
  np  <- length(lower)
  NP  <- if (is.null(NP)) as.integer(10 * np) else as.integer(NP)
  cat("\nDEoptim —", label, "(NP =", NP, ")...\n")
  out <- tryCatch(
    DEoptim(obj_fn, lower = lower, upper = upper,
            control = DEoptim.control(NP = NP, itermax = 2000L,
                                      F = 0.8, CR = 0.9,
                                      strategy = 2L, trace = FALSE)),
    error = function(e) { message("DEoptim erreur : ", e$message); NULL }
  )
  if (is.null(out)) stop(paste("DEoptim échoué pour", label))
  out
}

fit_emax <- run_deoptim(
  obj_emax,
  lower = c(log(l0 * 0.01),  log(1)),
  upper = c(log(l0 * 500),   log(C1_max_10 * 100)),
  label = "Emax + EC50 (k1 fixé)"
)

Emax_fit <- unname(exp(fit_emax$optim$bestmem[1]))
EC50_fit <- unname(exp(fit_emax$optim$bestmem[2]))

cat("\n=== Modèle Emax — Paramètres estimés ===\n")
cat("  Emax =", formatC(Emax_fit, format="e", digits=3), "/h",
    " =", round(Emax_fit*24, 4), "/j\n")
cat("  EC50 =", round(EC50_fit, 1), "µg/L\n")
cat("  k1   =", round(k1_fixed*24, 4), "/j  (fixé)\n")
cat("  Objectif :", round(fit_emax$optim$bestval, 5), "\n")

# Rapport Emax/l0 : indique si l'effet peut dépasser la croissance
cat("  Emax/l0 =", round(Emax_fit/l0, 2),
    if (Emax_fit > l0) "→ régression possible (Emax > l0)"
    else "→ inhibition partielle seulement (Emax < l0)", "\n")

# Fraction de C1_max utilisée à EC50
cat("  EC50 / C1_max_10 =", round(EC50_fit / C1_max_10 * 100, 1),
    "% → effet à 10 mg/kg =",
    round(C1_max_10 / (EC50_fit + C1_max_10) * 100, 1), "% Emax\n")

# =============================================================================
# 18. GRAPHIQUES PKPD — Modèle Emax
# =============================================================================

times_pd_full <- seq(0, 49 * 24, by = 4)
lev_pd   <- c("Contrôle", "3 mg/kg", "10 mg/kg")
doses_pd <- list(dose0_pd, dose3_pd, dose10_pd)

par_emax <- c(params_fixed, k1 = k1_fixed, Emax = Emax_fit, EC50 = EC50_fit)

df_sim_emax <- do.call(rbind, mapply(function(d, g)
  data.frame(t      = times_pd_full / 24,
             w      = simulate_multi(d, par_emax, times_pd_full),
             Groupe = g),
  doses_pd, lev_pd, SIMPLIFY = FALSE))
df_sim_emax$Groupe <- factor(df_sim_emax$Groupe, levels = lev_pd)

df_obs_pd <- rbind(
  data.frame(t = dat_ctrl$t/24, w = dat_ctrl$w, Groupe = "Contrôle"),
  data.frame(t = dat_d3$t/24,   w = dat_d3$w,   Groupe = "3 mg/kg"),
  data.frame(t = dat_d10$t/24,  w = dat_d10$w,  Groupe = "10 mg/kg")
)
df_obs_pd$Groupe <- factor(df_obs_pd$Groupe, levels = lev_pd)

COLS_PD <- c("Contrôle" = "#888888", "3 mg/kg" = "#27AE60", "10 mg/kg" = "#1B4F9E")

# Courbe Emax prédit vs observé
p_emax <- ggplot() +
  geom_vline(xintercept = DOSE_DAYS, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(data  = df_sim_emax, aes(x = t, y = w * 1000, color = Groupe),
            linewidth = 1.2) +
  geom_point(data = df_obs_pd,   aes(x = t, y = w * 1000, color = Groupe,
                                     shape = Groupe), size = 3) +
  scale_color_manual(values = COLS_PD) +
  scale_shape_manual(values = c(16, 17, 15)) +
  annotate("text", x = DOSE_DAYS,
           y = max(df_obs_pd$w * 1000, na.rm = TRUE) * 1.08,
           label = paste0("j", DOSE_DAYS), size = 2.8, color = "grey40") +
  labs(
    title    = "PKPD Simeoni-Emax — Prédit vs Observé",
    subtitle = paste0(
      "Fc-silent FGFR2-huBPA-LP1  |  Protocole Q2W x4  |  ",
      "Emax = ", formatC(Emax_fit, format="e", digits=2), " /h",
      "   EC50 = ", round(EC50_fit, 0), " µg/L",
      "   k1 = ", round(k1_fixed*24, 3), " /j"
    ),
    x = "Temps (jours)", y = "Volume tumoral (mm³)",
    color = NULL, shape = NULL,
    caption = "Lignes = modèle Emax   ●▲■ = données observées"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("scripts/plot_PKPD_Emax_FGFR2.png", p_emax,
       width = 9, height = 5.5, dpi = 150)
cat("\nGraphique → scripts/plot_PKPD_Emax_FGFR2.png\n")

# Résidus relatifs
df_res_emax <- do.call(rbind, lapply(lev_pd, function(g) {
  obs <- df_obs_pd[df_obs_pd$Groupe == g, ]
  dose_g <- switch(g, "Contrôle" = dose0_pd, "3 mg/kg" = dose3_pd,
                   "10 mg/kg" = dose10_pd)
  w_sim <- simulate_multi(dose_g, par_emax, obs$t * 24)
  data.frame(t = obs$t,
             resid = (obs$w - w_sim) / w_sim * 100,
             Groupe = g)
}))
df_res_emax$Groupe <- factor(df_res_emax$Groupe, levels = lev_pd)

p_res <- ggplot(df_res_emax, aes(x = t, y = resid, color = Groupe)) +
  geom_hline(yintercept = 0,      linewidth = 0.8, color = "grey40") +
  geom_hline(yintercept = c(-30, 30), linetype = "dashed", color = "grey70") +
  geom_point(size = 2.5) +
  scale_color_manual(values = COLS_PD) +
  labs(title    = "Résidus relatifs — Modèle Emax",
       subtitle = "(Observé − Prédit) / Prédit × 100",
       x = "Temps (jours)", y = "Résidu (%)", color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("scripts/plot_residus_Emax_FGFR2.png", p_res,
       width = 6, height = 4, dpi = 150)

# Courbe Emax — fonction de C1 (pour comprendre l'effet à chaque dose)
C1_seq <- seq(0, C1_max_10 * 1.1, length.out = 300)
df_emax_curve <- data.frame(
  C1     = C1_seq,
  effet  = Emax_fit * C1_seq / (EC50_fit + C1_seq),
  pct    = Emax_fit * C1_seq / (EC50_fit + C1_seq) / Emax_fit * 100
)
C1_d3  <- (dose3_pd  / PK["V1"])
C1_d10 <- (dose10_pd / PK["V1"])

p_ec <- ggplot(df_emax_curve, aes(x = C1, y = pct)) +
  geom_line(linewidth = 1.2, color = "#1A5276") +
  geom_vline(xintercept = C1_d3,  linetype = "dashed", color = "#27AE60") +
  geom_vline(xintercept = C1_d10, linetype = "dashed", color = "#1B4F9E") +
  annotate("text", x = C1_d3,  y = 5, label = "C1 max\n3 mg/kg",
           hjust = -0.1, size = 3.5, color = "#27AE60") +
  annotate("text", x = C1_d10, y = 5, label = "C1 max\n10 mg/kg",
           hjust = -0.1, size = 3.5, color = "#1B4F9E") +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey50") +
  annotate("text", x = EC50_fit, y = 53, label = paste0("EC50 = ", round(EC50_fit, 0), " µg/L"),
           hjust = 0, size = 3.5, color = "grey40") +
  labs(
    title    = "Courbe effet-concentration — Modèle Emax",
    subtitle = paste0("Emax = ", formatC(Emax_fit, format="e", digits=2),
                      " /h   EC50 = ", round(EC50_fit, 0), " µg/L"),
    x = "C1 (µg/L)", y = "% Emax atteint"
  ) +
  theme_bw(base_size = 13)

ggsave("scripts/plot_Emax_curve_FGFR2.png", p_ec,
       width = 7, height = 4.5, dpi = 150)
cat("Graphique → scripts/plot_Emax_curve_FGFR2.png\n")

# =============================================================================
# 19. TGI PRÉDIT — Modèle Emax
# =============================================================================

tgi_model <- function(w_ctrl_end, w_treat_end, w0_val) {
  delta_ctrl <- w_ctrl_end - w0_val
  if (is.na(delta_ctrl) || delta_ctrl <= 0) return(NA_real_)
  round((1 - (w_treat_end - w0_val) / delta_ctrl) * 100, 1)
}

t_seq_tgi <- seq(0, 49 * 24, by = 24)
wc_e  <- tail(simulate_multi(dose0_pd,  par_emax, t_seq_tgi), 1)
w3_e  <- tail(simulate_multi(dose3_pd,  par_emax, t_seq_tgi), 1)
w10_e <- tail(simulate_multi(dose10_pd, par_emax, t_seq_tgi), 1)

tgi3_emax  <- tgi_model(wc_e, w3_e,  w0)
tgi10_emax <- tgi_model(wc_e, w10_e, w0)

cat("\n=== TGI prédit à j49 — Modèle Emax ===\n")
cat("   3 mg/kg  :", ifelse(is.na(tgi3_emax),  "NA",
    ifelse(tgi3_emax > 100,  paste0(tgi3_emax,  " % (régression)"),
                             paste0(tgi3_emax,  " %"))), "\n")
cat("  10 mg/kg  :", ifelse(is.na(tgi10_emax), "NA",
    ifelse(tgi10_emax > 100, paste0(tgi10_emax, " % (régression)"),
                             paste0(tgi10_emax, " %"))), "\n")

# =============================================================================
# 20. SAUVEGARDE PKPD
# =============================================================================

save(k1_fixed, Emax_fit, EC50_fit, par_emax,
     params_fixed, w0,
     dat_ctrl, dat_d3, dat_d10,
     file = "scripts/resultats_PKPD_Emax_FGFR2.RData")
cat("\nRésultats → scripts/resultats_PKPD_Emax_FGFR2.RData\n")
