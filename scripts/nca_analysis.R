# =============================================================================
# Analyse Non-Compartimentale (NCA) — multi-animaux
# =============================================================================
library(PKNCA)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 2. Données brutes ---------------------------------------------------------
# BLQ terminal (336h, 504h) → NA
# 96h et 168h sont QUANTIFIABLES

pk_raw <- bind_rows(

  data.frame(
    subject    = "Animal_01",
    dose_mg_kg = 3,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     69800, 48500, 17700, 6730, 3550, 2200, 1060, NA,  NA )
  ),

  data.frame(
    subject    = "Animal_02",
    dose_mg_kg = 3,
    time = c(0,     0.083, 4,     24,    48,   72,   96,   168,  336, 504),
    conc = c(0,     73200, 50200, 16600, 6370, 3680, 2180, 1020, NA,  NA )
  )
)

cat("=== Données PK ===\n")
print(pk_raw)

# ── 3. Lambda_z manuel par animal ---------------------------------------------
# On utilise les deux derniers points quantifiables (96h et 168h) qui
# appartiennent à la phase d'élimination terminale mono-exponentielle.
# λz = -pente de ln(C) vs t  =>  λz = ln(C1/C2) / (t2 - t1)
# t1/2 = ln(2) / λz

terminal_times <- c(96, 168)   # points de la phase terminale

lambda_z_df <- pk_raw %>%
  filter(time %in% terminal_times, !is.na(conc)) %>%
  group_by(subject) %>%
  arrange(time) %>%
  summarise(
    t1      = time[1],
    t2      = time[2],
    c1      = conc[1],
    c2      = conc[2],
    lambda_z = log(c1 / c2) / (t2 - t1),
    half_life = log(2) / lambda_z,
    .groups = "drop"
  )

cat("\n=== Lambda_z et demi-vie (régression terminale manuelle) ===\n")
print(lambda_z_df)

# ── 4. PKNCA pour AUClast, Cmax, Tmax ----------------------------------------
pk_conc <- PKNCAconc(pk_raw, conc ~ time | subject)

pk_dose_data <- pk_raw %>%
  distinct(subject, dose_mg_kg) %>%
  mutate(time = 0) %>%
  rename(dose = dose_mg_kg)
pk_dose <- PKNCAdose(pk_dose_data, dose ~ time | subject)

intervals <- data.frame(
  start   = 0,
  end     = 168,
  cmax    = TRUE,
  tmax    = TRUE,
  auclast = TRUE
)

pk_data_obj <- PKNCAdata(
  data.conc = pk_conc,
  data.dose = pk_dose,
  intervals = intervals,
  options   = list(auc.method = "linear")
)

pk_results <- pk.nca(pk_data_obj)
results_df <- as.data.frame(pk_results$result)

# ── 5. Calculs manuels : AUCinf, CL, Vz ---------------------------------------
# AUCinf = AUClast + C_last / λz
# CL     = Dose / AUCinf            [mg/kg / h*ng/mL × 1e6 → mL/h/kg]
# Vz     = CL / λz                  [mL/h/kg / h⁻¹ → mL/kg]

wide <- results_df %>%
  filter(PPTESTCD %in% c("cmax", "tmax", "auclast")) %>%
  select(subject, PPTESTCD, PPORRES) %>%
  mutate(PPORRES = as.numeric(PPORRES)) %>%
  pivot_wider(names_from = PPTESTCD, values_from = PPORRES)

# C_last : dernière concentration quantifiable par animal
c_last_df <- pk_raw %>%
  filter(!is.na(conc)) %>%
  group_by(subject) %>%
  slice_max(time, n = 1) %>%
  select(subject, c_last = conc)

wide <- wide %>%
  left_join(lambda_z_df %>% select(subject, lambda_z, half_life), by = "subject") %>%
  left_join(c_last_df, by = "subject") %>%
  left_join(pk_raw %>% distinct(subject, dose_mg_kg), by = "subject") %>%
  mutate(
    aucinf   = auclast + c_last / lambda_z,
    cl_raw   = dose_mg_kg / aucinf,
    CL       = cl_raw * 1e6,        # mL/h/kg
    Vz       = CL / lambda_z,       # mL/kg
    Cmax_D   = cmax / dose_mg_kg
  )

# ── 6. Tableau de sortie -------------------------------------------------------
nca_summary <- wide %>%
  transmute(
    Dose_mg_kg      = round(dose_mg_kg, 3),
    Animal_Id       = subject,
    Half_life_h     = round(half_life, 1),
    Cmax_ng_mL      = round(cmax,      0),
    Cmax_D          = round(Cmax_D,    0),
    AUClast_h_ng_mL = round(auclast,   0),
    AUCinf_h_ng_mL  = round(aucinf,    0),
    Vz_mL_kg        = round(Vz,        0),
    CL_mL_h_kg      = round(CL,        2)
  )

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
cat("\nTableau exporté : nca_results.csv\n")

# ── 7. Graphiques -------------------------------------------------------------
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
       subtitle = "lambda_z estimé sur 96–168h (phase terminale)",
       x = "Temps (h)", y = "Concentration (ng/mL) — log", color = "Animal") +
  theme_pk
print(p_semilog)
ggsave("nca_semilog.png", plot = p_semilog, width = 8, height = 5, dpi = 300)

cat("\nGraphiques exportés : nca_linear.png  nca_semilog.png\n")
cat("Script terminé.\n")
