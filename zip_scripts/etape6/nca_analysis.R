# =============================================================================
# Analyse PK — modèle 2 compartiments IV bolus avec rxode2
# =============================================================================
# Modèle :   d/dt(A1) = -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
#            d/dt(A2) =   (Q/V1)*A1 - (Q/V2)*A2
#            C = A1 / V1
# Unités  :  A1, A2 [ng/kg]  |  V1, V2 [mL/kg]  |  CL, Q [mL/h/kg]  |  C [ng/mL]
# =============================================================================

# install.packages(c("rxode2", "dplyr", "tidyr", "ggplot2"))
library(rxode2)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. Données ----------------------------------------------------------------
# BLQ terminaux (336h, 504h) → NA  |  96h et 168h sont quantifiables

# Facteur de correction dose : 1.297 (dose réelle = dose nominale × 1.297)
# Doses corrigées : 3→4 | 10→13 | 20→26 | 30→39 mg/kg
# Concentrations corrigées : × 1.297

pk_raw <- bind_rows(

  # ── 4 mg/kg (nominale 3 mg/kg × 1.297) ──────────────
  data.frame(
    subject    = "Animal_01",
    dose_mg_kg = 3 * 1.297,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     69800, 48500, 17700, 6730, 3550, 2200, 1060, NA,  NA ) * 1.297
  ),

  data.frame(
    subject    = "Animal_02",
    dose_mg_kg = 3 * 1.297,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     73200, 50200, 16600, 6370, 3680, 2180, 1020, NA,  NA ) * 1.297
  ),

  # ── 13 mg/kg (nominale 10 mg/kg × 1.297) — à remplir ─
  # data.frame(
  #   subject    = "Animal_03",
  #   dose_mg_kg = 10 * 1.297,
  #   time = c(),   # heures depuis la dose
  #   conc = c() * 1.297
  # ),
  # data.frame(
  #   subject    = "Animal_04",
  #   dose_mg_kg = 10 * 1.297,
  #   time = c(),
  #   conc = c() * 1.297
  # ),

  # ── 26 mg/kg (nominale 20 mg/kg × 1.297) — à remplir ─
  # data.frame(
  #   subject    = "Animal_05",
  #   dose_mg_kg = 20 * 1.297,
  #   time = c(),
  #   conc = c() * 1.297
  # ),
  # data.frame(
  #   subject    = "Animal_06",
  #   dose_mg_kg = 20 * 1.297,
  #   time = c(),
  #   conc = c() * 1.297
  # ),

  # ── 39 mg/kg (nominale 30 mg/kg × 1.297) — à remplir ─
  # data.frame(
  #   subject    = "Animal_07",
  #   dose_mg_kg = 30 * 1.297,
  #   time = c(),
  #   conc = c() * 1.297
  # ),
  # data.frame(
  #   subject    = "Animal_08",
  #   dose_mg_kg = 30 * 1.297,
  #   time = c(),
  #   conc = c() * 1.297
  # )
)

# ── 2. Définition du modèle rxode2 --------------------------------------------
pk2cmt <- rxode2({
  d/dt(A1) <- -(CL / V1 + Q / V1) * A1 + (Q / V2) * A2
  d/dt(A2) <-  (Q / V1) * A1 - (Q / V2) * A2
  C        <- A1 / V1
})

# ── 3. Ajustement par animal (moindres carrés sur log-concentrations) ---------
# La fonction objectif minimise Σ [ln(Cobs) - ln(Cpred)]²
# Les paramètres sont optimisés sur l'échelle log (→ contrainte de positivité).

fit_animal <- function(animal_data) {

  dose  <- animal_data$dose_mg_kg[1]
  A1_0  <- dose * 1e6        # mg/kg → ng/kg  (1 mg = 1e6 ng)
  obs   <- animal_data %>% filter(time > 0, !is.na(conc))

  objective <- function(log_p) {
    p <- setNames(exp(log_p), c("CL", "V1", "Q", "V2"))
    tryCatch({
      sol  <- rxSolve(pk2cmt,
                      params = p,
                      inits  = c(A1 = A1_0, A2 = 0),
                      events = et(obs$time))
      pred <- sol$C
      if (any(!is.finite(pred) | pred <= 0)) return(1e10)
      sum((log(obs$conc) - log(pred))^2)
    }, error = function(e) 1e10)
  }

  # Valeurs initiales — estimées à partir des données :
  #   V1 ≈ dose/Cmax ≈ 3e6/70000 ≈ 43 mL/kg
  #   CL ≈ 1.73 mL/h/kg (NCA préliminaire)
  #   Q  ≈ 5 mL/h/kg, V2 ≈ 130 mL/kg
  log_p0 <- log(c(CL = 1.73, V1 = 43, Q = 5, V2 = 130))

  fit <- optim(log_p0, objective,
               method  = "Nelder-Mead",
               control = list(maxit = 20000, reltol = 1e-12))

  exp(fit$par) |> setNames(c("CL", "V1", "Q", "V2"))
}

# ── 4. Paramètres dérivés (analytiques, modèle 2-cmt) -------------------------
# Les constantes microscopiques donnent les valeurs propres alpha (rapide)
# et beta (terminale) du système biexponentiel :
#   C(t) = A·exp(-alpha·t) + B·exp(-beta·t)
# La demi-vie terminale est t1/2 = ln(2) / beta.
# Le volume de distribution terminal : Vz = CL / beta.

derived_params <- function(p) {
  k10  <- p["CL"] / p["V1"]
  k12  <- p["Q"]  / p["V1"]
  k21  <- p["Q"]  / p["V2"]
  S    <- k10 + k12 + k21
  disc <- sqrt(S^2 - 4 * k10 * k21)
  alpha <- (S + disc) / 2
  beta  <- (S - disc) / 2
  list(
    half_life_alpha = log(2) / alpha,
    half_life_beta  = log(2) / beta,
    Vz              = p["CL"] / beta
  )
}

# ── 5. AUClast (trapèzes linéaires) ------------------------------------------
auclast_animal <- function(animal_data) {
  d <- animal_data %>% filter(!is.na(conc)) %>% arrange(time)
  n <- nrow(d)
  sum(diff(d$time) * (d$conc[-n] + d$conc[-1]) / 2)
}

# ── 6. Boucle sur les animaux -------------------------------------------------
animals    <- unique(pk_raw$subject)
results    <- list()
pk2cmt_raw <- list()   # paramètres bruts du fit 2-cmt
fit_params <- list()   # paramètres stockés pour les graphiques (évite double fit)

for (anim in animals) {
  cat("\n--- Ajustement :", anim, "---\n")
  dat   <- pk_raw %>% filter(subject == anim)
  dose  <- dat$dose_mg_kg[1]

  p     <- fit_animal(dat)
  fit_params[[anim]] <- p     # stocker pour réutilisation dans les graphiques
  drv   <- derived_params(p)
  aucl  <- auclast_animal(dat)

  # AUCinf = AUClast + C_last / beta  (C_last = dernière conc quantifiable)
  c_last <- dat %>% filter(!is.na(conc)) %>% slice_max(time, n=1) %>% pull(conc)
  beta   <- log(2) / drv$half_life_beta
  aucinf <- aucl + c_last / beta

  results[[anim]] <- data.frame(
    Dose_mg_kg      = dose,
    Animal_Id       = anim,
    Half_life_h     = round(drv$half_life_beta, 1),
    Cmax_ng_mL      = round(max(dat$conc, na.rm=TRUE), 0),
    Cmax_D          = round(max(dat$conc, na.rm=TRUE) / dose, 0),
    AUClast_h_ng_mL = round(aucl,   0),
    AUCinf_h_ng_mL  = round(aucinf, 0),
    Vz_mL_kg        = round(drv$Vz, 0),
    CL_mL_h_kg      = round(p["CL"], 2)
  )

  # Stocker les 4 paramètres du modèle 2-cmt (unités : mL/h/kg et mL/kg)
  pk2cmt_raw[[anim]] <- data.frame(
    Animal_Id  = anim,
    dose_mg_kg = dose,
    CL_mL_h_kg = round(p["CL"], 4),
    V1_mL_kg   = round(p["V1"], 4),
    Q_mL_h_kg  = round(p["Q"],  4),
    V2_mL_kg   = round(p["V2"], 4),
    t12_beta_h = round(drv$half_life_beta, 2),
    Vz_mL_kg   = round(drv$Vz, 2)
  )
}

nca_summary <- bind_rows(results)
row.names(nca_summary) <- NULL

units_row <- data.frame(
  Dose_mg_kg = "mg/kg", Animal_Id = "", Half_life_h = "h",
  Cmax_ng_mL = "ng/mL", Cmax_D = "ng/mL/mg/kg",
  AUClast_h_ng_mL = "h*ng/mL", AUCinf_h_ng_mL = "h*ng/mL",
  Vz_mL_kg = "mL/kg", CL_mL_h_kg = "mL/h/kg",
  stringsAsFactors = FALSE
)

cat("\n=== Paramètres PK — tableau de sortie ===\n")
print(rbind(units_row, nca_summary), row.names = FALSE)
write.csv(nca_summary, "nca_results.csv", row.names = FALSE)

# ── Export paramètres 2-cmt + ligne moyenne ───────────────────────────────
pk2cmt_df <- bind_rows(pk2cmt_raw)
mean_row   <- data.frame(
  Animal_Id  = "Mean",
  dose_mg_kg = NA,
  CL_mL_h_kg = mean(pk2cmt_df$CL_mL_h_kg),
  V1_mL_kg   = mean(pk2cmt_df$V1_mL_kg),
  Q_mL_h_kg  = mean(pk2cmt_df$Q_mL_h_kg),
  V2_mL_kg   = mean(pk2cmt_df$V2_mL_kg),
  t12_beta_h  = mean(pk2cmt_df$t12_beta_h),
  Vz_mL_kg   = mean(pk2cmt_df$Vz_mL_kg)
)
pk2cmt_df <- bind_rows(pk2cmt_df, mean_row)
row.names(pk2cmt_df) <- NULL

cat("\n=== Paramètres modèle 2-cmt (rxode2 fit) ===\n")
print(pk2cmt_df, row.names = FALSE)
write.csv(pk2cmt_df, "pk2cmt_params.csv", row.names = FALSE)
cat("Tableau exporté : pk2cmt_params.csv\n")
cat("\nTableau exporté : nca_results.csv\n")

# ── 7. Courbes ajustées + données observées -----------------------------------
# Générer les prédictions du modèle ajusté pour chaque animal

t_max_obs  <- max(pk_raw$time, na.rm = TRUE)
times_pred <- seq(0.083, t_max_obs * 1.05, length.out = 400)

pred_list <- list()
for (anim in animals) {
  dat  <- pk_raw %>% filter(subject == anim)
  dose <- dat$dose_mg_kg[1]
  p    <- fit_params[[anim]]   # réutilise les paramètres du premier fit

  sol <- rxSolve(pk2cmt,
                 params = p,
                 inits  = c(A1 = dose * 1e6, A2 = 0),
                 events = et(times_pred))

  pred_list[[anim]] <- data.frame(
    subject = anim, time = sol$time, conc = sol$C
  )
}
pred_df <- bind_rows(pred_list)

obs_df  <- pk_raw %>% filter(!is.na(conc), conc > 0) %>%
  mutate(Dose = paste0(round(dose_mg_kg), " mg/kg"))
pred_df <- pred_df %>%
  left_join(pk_raw %>% select(subject, dose_mg_kg) %>% distinct(), by = "subject") %>%
  mutate(Dose = paste0(round(dose_mg_kg), " mg/kg"))

dose_cols_nca <- c("4 mg/kg"  = "#2166ac",
                   "13 mg/kg" = "#4dac26",
                   "26 mg/kg" = "#f4a582",
                   "39 mg/kg" = "#d6604d")

# Formes manuelles : même forme par dose (2 animaux/dose), plein vs creux
animal_shapes <- c(
  "Animal_01" = 16, "Animal_02" = 1,   # 4  mg/kg  : cercle plein / creux
  "Animal_03" = 17, "Animal_04" = 2,   # 13 mg/kg  : triangle plein / creux
  "Animal_05" = 15, "Animal_06" = 0,   # 26 mg/kg  : carré plein / creux
  "Animal_07" = 18, "Animal_08" = 5    # 39 mg/kg  : losange plein / creux
)

theme_pk <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# 7a. Linéaire
p_linear <- ggplot() +
  geom_line(data  = pred_df, aes(x=time, y=conc, color=Dose, group=subject), linewidth=0.9) +
  geom_point(data = obs_df,  aes(x=time, y=conc, color=Dose, shape=subject), size=3, na.rm=TRUE) +
  scale_color_manual(values = dose_cols_nca) +
  scale_shape_manual(values = animal_shapes) +
  labs(title    = "PK Profile — 2-Compartment Model (rxode2 fit)",
       subtitle = "Points = observations ; lines = fitted model | ●○ = animal 1/2 per dose",
       x = "Time (h)", y = "Concentration (ng/mL)",
       color = "Dose", shape = "Animal") +
  theme_pk
print(p_linear)
ggsave("nca_linear.png", plot = p_linear, width = 8, height = 5, dpi = 300)

# 7b. Semi-logarithmique
p_semilog <- ggplot() +
  geom_line(data  = pred_df, aes(x=time, y=conc, color=Dose, group=subject), linewidth=0.9) +
  geom_point(data = obs_df,  aes(x=time, y=conc, color=Dose, shape=subject), size=3, na.rm=TRUE) +
  scale_color_manual(values = dose_cols_nca) +
  scale_shape_manual(values = animal_shapes) +
  scale_y_log10() +
  labs(title    = "PK Profile — 2-Compartment Model (rxode2 fit) — Semi-log scale",
       subtitle = "Points = observations ; lines = fitted model | ●○ = animal 1/2 per dose",
       x = "Time (h)", y = "Concentration (ng/mL) — log",
       color = "Dose", shape = "Animal") +
  theme_pk
print(p_semilog)
ggsave("nca_semilog.png", plot = p_semilog, width = 8, height = 5, dpi = 300)

cat("\nGraphiques exportés : nca_linear.png  nca_semilog.png\n")

# ── 8. Tableau récapitulatif paramètres PK ────────────────────────────────────
# Colonnes : Animal | Dose | Cmax | t½β | CL | V1 | Q | V2
# Ligne finale : moyenne des 8 animaux

recap_df <- merge(
  nca_summary[, c("Animal_Id", "Dose_mg_kg", "Cmax_ng_mL", "Half_life_h")],
  pk2cmt_df[pk2cmt_df$Animal_Id != "Mean",
            c("Animal_Id", "CL_mL_h_kg", "V1_mL_kg", "Q_mL_h_kg", "V2_mL_kg")],
  by = "Animal_Id"
)
recap_df <- recap_df[order(recap_df$Dose_mg_kg), ]

mean_row2 <- pk2cmt_df[pk2cmt_df$Animal_Id == "Mean", ]

# Format numérique → caractère pour affichage
fmt <- function(x, digits = 2, big = FALSE) {
  if (big) formatC(round(x), format = "d", big.mark = " ")
  else     sprintf(paste0("%.", digits, "f"), x)
}

build_display <- function(df) {
  data.frame(
    "Animal"        = df$Animal_Id,
    "Dose\n(mg/kg)" = fmt(df$Dose_mg_kg, 1),
    "Cmax\n(ng/mL)" = fmt(df$Cmax_ng_mL, 0, big = TRUE),
    "t½β\n(h)"      = fmt(df$Half_life_h, 1),
    "CL\n(mL/h/kg)" = fmt(df$CL_mL_h_kg, 2),
    "V1\n(mL/kg)"   = fmt(df$V1_mL_kg,   1),
    "Q\n(mL/h/kg)"  = fmt(df$Q_mL_h_kg,  2),
    "V2\n(mL/kg)"   = fmt(df$V2_mL_kg,   1),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

tbl_body <- build_display(recap_df)

tbl_mean <- data.frame(
  "Animal"        = "Mean",
  "Dose\n(mg/kg)" = "—",
  "Cmax\n(ng/mL)" = "—",
  "t½β\n(h)"      = fmt(mean_row2$t12_beta_h, 1),
  "CL\n(mL/h/kg)" = fmt(mean_row2$CL_mL_h_kg, 2),
  "V1\n(mL/kg)"   = fmt(mean_row2$V1_mL_kg,   1),
  "Q\n(mL/h/kg)"  = fmt(mean_row2$Q_mL_h_kg,  2),
  "V2\n(mL/kg)"   = fmt(mean_row2$V2_mL_kg,   1),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

display_tbl <- rbind(tbl_body, tbl_mean)

# ── Rendu ggplot2 (sans dépendance externe) ───────────────────────────────────
n_col  <- ncol(display_tbl)
n_body <- nrow(display_tbl)
n_row  <- n_body + 1   # + ligne d'en-tête

col_names <- colnames(display_tbl)
col_w     <- c(2.2, 1.3, 1.8, 1.2, 1.4, 1.3, 1.4, 1.3)   # largeur relative par col
col_x     <- cumsum(c(0, col_w[-n_col])) + col_w / 2       # centre de chaque col

# Long format — ordre ligne par ligne (t + as.vector corrige le column-major d'unlist)
row_major <- function(df) as.vector(t(as.matrix(df)))

cells <- data.frame(
  ri  = c(rep(1, n_col),
          rep(2:(n_row - 1), each = n_col),
          rep(n_row, n_col)),
  ci  = rep(1:n_col, n_row),
  val = c(col_names,
          row_major(tbl_body),
          row_major(tbl_mean)),
  stringsAsFactors = FALSE
)

# Géométrie
cells$cx  <- col_x[cells$ci]
cells$cy  <- (n_row + 1) - cells$ri   # y inverse pour lire de haut en bas

# Style cellules
cells$bg   <- ifelse(cells$ri == 1,     "#2c3e50",
              ifelse(cells$ri == n_row, "#dfe6e9",
              ifelse(cells$ri %% 2 == 0, "#f8f9fa", "white")))
cells$fcol <- ifelse(cells$ri == 1, "white", "#2d3436")
cells$face <- ifelse(cells$ri == 1, "bold",
              ifelse(cells$ri == n_row, "bold.italic", "plain"))
cells$sz   <- ifelse(cells$ri == 1, 3.6, 3.3)

p_tbl <- ggplot(cells, aes(x = cx, y = cy)) +
  geom_tile(aes(fill = I(bg), width = col_w[ci]),
            height = 0.85, color = "white", linewidth = 0.6) +
  geom_text(aes(label = val, color = I(fcol), fontface = face, size = I(sz)),
            lineheight = 0.85) +
  scale_x_continuous(limits = c(0, sum(col_w)), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0.5, n_row + 0.5), expand = c(0, 0)) +
  labs(
    title    = "PK Parameters — 2-Compartment Model (rxode2 fit)",
    subtitle = "IV bolus · NHP (cynomolgus macaque) · FGFR2 inhibitor  |  CL, V weight-normalized (mL/h/kg, mL/kg)"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13,
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey45", size = 10,
                                 margin = margin(b = 6)),
    plot.margin   = margin(10, 10, 10, 10)
  )

ggsave("pk_recap_table.pdf", plot = p_tbl,
       width = 11, height = 4.2, device = cairo_pdf)
ggsave("pk_recap_table.png", plot = p_tbl,
       width = 11, height = 4.2, dpi = 300)

cat("Tableau exporté : pk_recap_table.pdf  pk_recap_table.png\n")
cat("Script terminé.\n")
