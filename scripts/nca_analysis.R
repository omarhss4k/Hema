# =============================================================================
# Analyse Non-Compartimentale (NCA) — multi-animaux
# =============================================================================
library(PKNCA)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 2. Données brutes ---------------------------------------------------------
# BLQ terminal (336h, 504h) → NA
# 96h et 168h sont QUANTIFIABLES → inclus dans l'AUC et la régression terminale

pk_raw <- bind_rows(

  # ---- Animal 1 ----
  data.frame(
    subject    = "Animal_01",
    dose_mg_kg = 3,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     69800, 48500, 17700, 6730, 3550, 2200, 1060, NA,  NA )
  ),

  # ---- Animal 2 ----
  data.frame(
    subject    = "Animal_02",
    dose_mg_kg = 3,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     73200, 50200, 16600, 6370, 3680, 2180, 1020, NA,  NA )
  )
  # Pour ajouter un animal, copiez-collez un bloc supplémentaire ici
)

cat("=== Données PK ===\n")
print(pk_raw)

# ── 4. Objets PKNCA -----------------------------------------------------------
pk_conc <- PKNCAconc(pk_raw, conc ~ time | subject)

pk_dose_data <- pk_raw %>%
  distinct(subject, dose_mg_kg) %>%
  mutate(time = 0) %>%
  rename(dose = dose_mg_kg)

pk_dose <- PKNCAdose(pk_dose_data, dose ~ time | subject)

# ── 5. Intervalles et paramètres ----------------------------------------------
# end = 168 : dernier point quantifiable (BLQ à 336h et 504h exclus)
# lambda.z.time.range = c(96, 168) : restreint la régression log-linéaire aux
#   deux points terminaux de la phase d'élimination. Les points 24–72h sont
#   encore en phase de distribution (pente trop raide) et biaisent lambda_z.
#   Vérification manuelle :
#     Animal 1 : λz = ln(2200/1060)/(168-96) = 0.00987 h⁻¹ → t½ = 70.2h
#     Animal 2 : λz = ln(2180/1020)/(168-96) = 0.01028 h⁻¹ → t½ = 67.4h

intervals <- data.frame(
  start      = 0,
  end        = 168,
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
    auc.method          = "linear",
    lambda.z.time.range = c(96, 168),  # phase d'élimination terminale uniquement
    min.hl.points       = 2            # autoriser 2 points (96h et 168h)
  )
)

# ── 6. Calcul NCA -------------------------------------------------------------
pk_results <- pk.nca(pk_data_obj)
results_df <- as.data.frame(pk_results$result)

# ── 7. Tableau de sortie ------------------------------------------------------
wide <- results_df %>%
  filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "aucinf.obs",
                         "half.life", "cl.obs", "vz.obs")) %>%
  select(subject, PPTESTCD, PPORRES) %>%
  mutate(PPORRES = as.numeric(PPORRES)) %>%
  pivot_wider(names_from = PPTESTCD, values_from = PPORRES)

dose_map    <- pk_raw %>% distinct(subject, dose_mg_kg)
wide        <- left_join(wide, dose_map, by = "subject")
dose_ok_vec <- !is.na(wide$dose_mg_kg) & wide$dose_mg_kg > 0

# Conversions d'unités (dose mg/kg, concentrations ng/mL) :
#   CL [mg/kg / h*ng/mL] × 1e6 → mL/h/kg   (1 mg = 1e6 ng)
#   Vz [mg*mL / kg*ng]   × 1e6 → mL/kg     (1 mg = 1e6 ng)
wide <- wide %>%
  mutate(
    cl.obs = ifelse(dose_ok_vec, cl.obs * 1e6, NA_real_),
    vz.obs = ifelse(dose_ok_vec, vz.obs * 1e6, NA_real_),
    Cmax_D = ifelse(dose_ok_vec, cmax / dose_mg_kg, NA_real_)
  )

nca_summary <- wide %>%
  transmute(
    Dose_mg_kg      = round(dose_mg_kg, 3),
    Animal_Id       = subject,
    Half_life_h     = round(half.life,  1),
    Cmax_ng_mL      = round(cmax,       0),
    Cmax_D          = round(Cmax_D,     0),
    AUClast_h_ng_mL = round(auclast,    0),
    AUCinf_h_ng_mL  = round(aucinf.obs, 0),
    Vz_mL_kg        = round(vz.obs,     0),
    CL_mL_h_kg      = round(cl.obs,     2)
  )

units_row <- data.frame(
  Dose_mg_kg      = "mg/kg", Animal_Id = "",
  Half_life_h     = "h",     Cmax_ng_mL = "ng/mL",
  Cmax_D          = "ng/mL/mg/kg",
  AUClast_h_ng_mL = "h*ng/mL", AUCinf_h_ng_mL = "h*ng/mL",
  Vz_mL_kg        = "mL/kg",   CL_mL_h_kg = "mL/h/kg",
  stringsAsFactors = FALSE
)

cat("\n=== Paramètres PK — tableau de sortie ===\n")
print(rbind(units_row, nca_summary), row.names = FALSE)

write.csv(nca_summary, "nca_results.csv", row.names = FALSE)
cat("\nTableau exporté : nca_results.csv\n")

# ── 8. Graphiques -------------------------------------------------------------
pk_plot_data <- pk_raw %>% filter(!is.na(conc))

theme_pk <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

p_linear <- ggplot(pk_plot_data, aes(x = time, y = conc,
                                     color = subject, group = subject)) +
  geom_line(linewidth = 0.9) + geom_point(size = 3) +
  labs(title = "Profil concentration-temps — Échelle linéaire",
       subtitle = "BLQ (336h, 504h) exclus",
       x = "Temps (h)", y = "Concentration (ng/mL)", color = "Animal") +
  theme_pk
print(p_linear)
ggsave("nca_linear.png", plot = p_linear, width = 8, height = 5, dpi = 300)

p_semilog <- ggplot(pk_plot_data %>% filter(conc > 0),
                    aes(x = time, y = conc, color = subject, group = subject)) +
  geom_line(linewidth = 0.9) + geom_point(size = 3) +
  scale_y_log10() +
  labs(title = "Profil concentration-temps — Échelle semi-logarithmique",
       subtitle = "T=0 et BLQ exclus ; lambda_z estimé sur 96–168h",
       x = "Temps (h)", y = "Concentration (ng/mL) — log", color = "Animal") +
  theme_pk
print(p_semilog)
ggsave("nca_semilog.png", plot = p_semilog, width = 8, height = 5, dpi = 300)

cat("\nGraphiques exportés : nca_linear.png  nca_semilog.png\n")
cat("Script terminé.\n")
