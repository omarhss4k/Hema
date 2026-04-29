# =============================================================================
# Analyse Non-Compartimentale (NCA) — multi-animaux
# =============================================================================
# Structure : un bloc de données par animal → PKNCA calcule les paramètres
# indépendamment pour chaque sujet → une ligne par animal dans le tableau.
# =============================================================================

# ── 1. Packages ---------------------------------------------------------------
# install.packages(c("PKNCA", "tidyverse"))
library(PKNCA)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 2. Données brutes — AJOUTER / MODIFIER ICI --------------------------------
# Ajoutez autant de blocs que d'animaux.
# BLQ terminal → NA   |   BLQ à T=0 → 0   (convention NCA)

pk_raw <- bind_rows(

  # ---- Animal 1 ----
  data.frame(
    subject = "Animal_01",
    dose_mg_kg = 3,          # dose de cet animal (mg/kg)
    time = c(0,      0.083, 4,  24, 48, 72, 96, 168),
    conc = c(0,      NA,    NA, NA, NA, NA, NA, NA )
    #                ↑ remplacez les NA par vos concentrations (ng/mL)
    #                                             ↑ BLQ → laisser NA
  ),

  # ---- Animal 2 ----
  data.frame(
    subject = "Animal_02",
    dose_mg_kg = 3,          # dose de cet animal (mg/kg)
    time = c(0,      0.083, 4,  24, 48, 72, 96, 168),
    conc = c(0,      NA,    NA, NA, NA, NA, NA, NA )
    #                ↑ remplacez les NA par vos concentrations (ng/mL)
  )

  # Pour un 3e animal, copiez-collez un bloc supplémentaire ici
)

# ── 3. Vérification des données -----------------------------------------------
cat("=== Données PK ===\n")
print(pk_raw)

# ── 4. Objets PKNCA -----------------------------------------------------------
pk_conc <- PKNCAconc(
  data    = pk_raw,
  formula = conc ~ time | subject
)

pk_dose_data <- pk_raw %>%
  distinct(subject, dose_mg_kg) %>%
  mutate(time = 0) %>%
  rename(dose = dose_mg_kg)

pk_dose <- PKNCAdose(
  data    = pk_dose_data,
  formula = dose ~ time | subject
)

# ── 5. Intervalles et paramètres ----------------------------------------------
intervals <- data.frame(
  start      = 0,
  end        = 72,   # dernier point quantifiable (avant les BLQ)
  cmax       = TRUE,
  tmax       = TRUE,
  auclast    = TRUE,
  aucinf.obs = TRUE,
  half.life  = TRUE,
  cl.obs     = TRUE,
  vz.obs     = TRUE
)

pk_data_obj <- PKNCAdata(
  data.conc = pk_conc,
  data.dose = pk_dose,
  intervals = intervals,
  options   = list(
    auc.method = "linear",
    # Restreindre lambda_z aux points 24-72h (phase d'élimination terminale).
    # Sans cette contrainte, PKNCA peut inclure le point 4h (distribution)
    # et sur-estimer lambda_z → sous-estimer t1/2.
    lambda.z.time.range = c(24, 72)
  )
)

# ── 6. Calcul NCA -------------------------------------------------------------
pk_results  <- pk.nca(pk_data_obj)
results_df  <- as.data.frame(pk_results$result)

# ── 7. Tableau de sortie multi-animaux ----------------------------------------
# Pivot : une ligne par sujet, une colonne par paramètre
wide <- results_df %>%
  filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "aucinf.obs",
                         "half.life", "cl.obs", "vz.obs")) %>%
  select(subject, PPTESTCD, PPORRES) %>%
  mutate(PPORRES = as.numeric(PPORRES)) %>%
  pivot_wider(names_from = PPTESTCD, values_from = PPORRES)

# Rattacher la dose par sujet
dose_map <- pk_raw %>% distinct(subject, dose_mg_kg)
wide <- left_join(wide, dose_map, by = "subject")

dose_ok_vec <- !is.na(wide$dose_mg_kg) & wide$dose_mg_kg > 0

# Conversions d'unités (dose mg/kg, concentrations ng/mL) :
#   CL  [mg/kg / h*ng/mL] × 1e6  → mL/h/kg  (1 mg = 1e6 ng)
#   Vz  [mg*mL / kg*ng]   × 1000 → mL/g      (1 kg = 1000 g)
wide <- wide %>%
  mutate(
    cl.obs  = ifelse(dose_ok_vec, cl.obs  * 1e6,  NA_real_),
    vz.obs  = ifelse(dose_ok_vec, vz.obs  * 1000, NA_real_),
    Cmax_D  = ifelse(dose_ok_vec, cmax / dose_mg_kg, NA_real_)
  )

nca_summary <- wide %>%
  transmute(
    Dose_mg_kg      = round(dose_mg_kg, 3),
    Animal_Id       = subject,
    Half_life_h     = round(half.life,  3),
    Cmax_ng_mL      = round(cmax,       2),
    Cmax_D          = round(Cmax_D,     4),
    AUClast_h_ng_mL = round(auclast,    2),
    AUCinf_h_ng_mL  = round(aucinf.obs, 2),
    Vz_mL_g         = round(vz.obs,     4),
    CL_mL_h_kg      = round(cl.obs,     4)
  )

units_row <- data.frame(
  Dose_mg_kg      = "mg/kg",
  Animal_Id       = "",
  Half_life_h     = "h",
  Cmax_ng_mL      = "ng/mL",
  Cmax_D          = "kg*ng/mL / mg/kg",
  AUClast_h_ng_mL = "h*ng/mL",
  AUCinf_h_ng_mL  = "h*ng/mL",
  Vz_mL_g         = "mL/g",
  CL_mL_h_kg      = "mL/h/kg",
  stringsAsFactors = FALSE
)

cat("\n=== Paramètres PK — tableau de sortie ===\n")
print(rbind(units_row, nca_summary), row.names = FALSE)

write.csv(nca_summary, "nca_results.csv", row.names = FALSE)
cat("\nTableau exporté : nca_results.csv\n")

# ── 8. Graphiques multi-animaux -----------------------------------------------
pk_plot_data <- pk_raw %>% filter(!is.na(conc))

theme_pk <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# 8a. Échelle linéaire
p_linear <- ggplot(pk_plot_data, aes(x = time, y = conc,
                                     color = subject, group = subject)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  labs(
    title    = "Profil concentration-temps — Échelle linéaire",
    subtitle = "BLQ (96h, 168h) exclus",
    x        = "Temps (h)", y = "Concentration (ng/mL)", color = "Animal"
  ) +
  theme_pk
print(p_linear)
ggsave("nca_linear.png", plot = p_linear, width = 8, height = 5, dpi = 300)

# 8b. Échelle semi-logarithmique
p_semilog <- ggplot(pk_plot_data %>% filter(conc > 0),
                    aes(x = time, y = conc, color = subject, group = subject)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_y_log10() +
  labs(
    title    = "Profil concentration-temps — Échelle semi-logarithmique",
    subtitle = "T=0 et BLQ exclus ; pente terminale = lambda_z (24–72h)",
    x        = "Temps (h)", y = "Concentration (ng/mL) — log", color = "Animal"
  ) +
  theme_pk
print(p_semilog)
ggsave("nca_semilog.png", plot = p_semilog, width = 8, height = 5, dpi = 300)

cat("\nGraphiques exportés : nca_linear.png  nca_semilog.png\n")
cat("Script terminé.\n")
