# =============================================================================
# Analyse Non-Compartimentale (NCA) - Script R
# =============================================================================
# Ce script réalise une NCA complète avec :
#   - Gestion des BLQ (Below Limit of Quantification)
#   - Calcul des paramètres PK via PKNCA
#   - Estimation de lambda_z sur les points terminaux en phase log-linéaire
#   - Graphiques en échelle linéaire et semi-logarithmique
# =============================================================================

# ── 1. Packages ---------------------------------------------------------------
# Installer PKNCA si nécessaire :
# install.packages("PKNCA")
# install.packages("tidyverse")

library(PKNCA)
library(dplyr)
library(ggplot2)

# ── 2. Données brutes ---------------------------------------------------------
# IMPORTANT : Remplacez les valeurs NA par vos concentrations réelles
# pour les points 0.083h, 4h, 24h, 48h et 72h.
# Les deux derniers points (96h et 168h) sont BLQ → traités comme NA
# afin de ne pas biaiser l'estimation de la pente d'élimination (lambda_z).

time_h <- c(0, 0.083, 4, 24, 48, 72, 96, 168)

conc_raw <- c(
  0,      # T = 0    : valeur BLQ remplacée par 0 (convention NCA)
  NA,     # T = 0.083 h → REMPLACEZ par votre valeur (ng/mL)
  NA,     # T = 4    h → REMPLACEZ par votre valeur (ng/mL)
  NA,     # T = 24   h → REMPLACEZ par votre valeur (ng/mL)
  NA,     # T = 48   h → REMPLACEZ par votre valeur (ng/mL)
  NA,     # T = 72   h → REMPLACEZ par votre valeur (ng/mL)
  NA,     # T = 96   h : BLQ → NA (exclure de la régression terminale)
  NA      # T = 168  h : BLQ → NA (exclure de la régression terminale)
)

# ── 3. Construction du dataframe ----------------------------------------------
pk_data <- data.frame(
  time = time_h,
  conc = conc_raw,
  subject = "S01"        # identifiant sujet requis par PKNCA
)

# Afficher un résumé des données (utile pour vérification)
cat("=== Données PK ===\n")
print(pk_data)

# ── 4. Objets PKNCA -----------------------------------------------------------
# PKNCAconc  : objet contenant les concentrations et les temps
# PKNCAdose  : objet contenant la dose administrée
#              → remplacez 1 par la dose réelle (même unité que pour CL/Vz)

pk_conc <- PKNCAconc(
  data    = pk_data,
  formula = conc ~ time | subject
)

# Dose administrée (à renseigner) : si inconnu, laisser NA et CL/Vz seront NA
dose_value <- NA_real_   # Ex. : 100  pour 100 mg (ou ng selon unité souhaitée)

pk_dose_data <- data.frame(
  time    = 0,
  dose    = dose_value,
  subject = "S01"
)

pk_dose <- PKNCAdose(
  data    = pk_dose_data,
  formula = dose ~ time | subject
)

# ── 5. Intervalles d'intérêt et paramètres à calculer ------------------------
# AUClast  : AUC du premier au dernier point quantifiable (trapèzes linéaires)
# AUCinf   : AUC extrapolée jusqu'à l'infini via lambda_z
# Cmax     : concentration maximale observée
# Tmax     : temps correspondant à Cmax
# half.life: demi-vie terminale (t1/2 = ln2 / lambda_z)
# cl.obs   : clairance = dose / AUCinf  (NA si dose non renseignée)
# vz.obs   : volume de distribution terminal = dose / (lambda_z × AUCinf)

intervals <- data.frame(
  start    = 0,
  end      = 72,       # dernier point quantifiable (avant les BLQ)
  cmax     = TRUE,
  tmax     = TRUE,
  auclast  = TRUE,
  aucinf.obs = TRUE,
  half.life  = TRUE,
  cl.obs     = TRUE,
  vz.obs     = TRUE
)

pk_data_obj <- PKNCAdata(
  data.conc  = pk_conc,
  data.dose  = pk_dose,
  intervals  = intervals,
  # Méthode des trapèzes linéaires pour le calcul de l'AUC
  options    = list(auc.method = "linear")
)

# ── 6. Calcul des paramètres NCA ----------------------------------------------
# PKNCA estime automatiquement lambda_z (pente d'élimination) par régression
# log-linéaire sur les points terminaux en phase monoexponentielle :
#   ln(C) = ln(C0_terminal) - lambda_z × t
# L'algorithme sélectionne le meilleur sous-ensemble de points terminaux
# (r² ajusté maximal) en excluant Cmax et les BLQ/NA.
# lambda_z est ensuite utilisé pour calculer t1/2, AUCinf, CL et Vz.

pk_results <- pk.nca(pk_data_obj)

# Afficher les résultats
cat("\n=== Résultats NCA ===\n")
results_df <- as.data.frame(pk_results$result)
print(results_df)

# Extraire les paramètres clés dans un tableau propre
params_of_interest <- c("cmax", "tmax", "auclast", "aucinf.obs",
                        "half.life", "cl.obs", "vz.obs")
summary_table <- results_df %>%
  filter(PPTESTCD %in% params_of_interest) %>%
  select(PPTESTCD, PPORRES)

cat("\n=== Paramètres PK clés ===\n")
print(summary_table)

# ── 7. Graphiques -------------------------------------------------------------
# Palette et thème communs
pk_plot_data <- pk_data %>%
  filter(!is.na(conc))   # exclure les NA pour le tracé

theme_pk <- theme_bw(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.position = "none"
  )

# 7a. Echelle linéaire ---------------------------------------------------------
p_linear <- ggplot(pk_plot_data, aes(x = time, y = conc)) +
  geom_line(color = "#2c7bb6", linewidth = 0.9) +
  geom_point(color = "#2c7bb6", size = 3) +
  labs(
    title    = "Profil concentration-temps — Échelle linéaire",
    subtitle = "Les points BLQ (96h, 168h) sont exclus",
    x        = "Temps (h)",
    y        = "Concentration (ng/mL)"
  ) +
  theme_pk

print(p_linear)
ggsave("nca_linear.png", plot = p_linear, width = 8, height = 5, dpi = 300)

# 7b. Echelle semi-logarithmique -----------------------------------------------
# On retire les concentrations nulles (T=0) pour éviter log(0) = -Inf
pk_semilog_data <- pk_plot_data %>%
  filter(conc > 0)

p_semilog <- ggplot(pk_semilog_data, aes(x = time, y = conc)) +
  geom_line(color = "#d7191c", linewidth = 0.9) +
  geom_point(color = "#d7191c", size = 3) +
  scale_y_log10() +
  labs(
    title    = "Profil concentration-temps — Échelle semi-logarithmique",
    subtitle = "T=0 (C=0) et BLQ exclus ; la phase terminale linéaire correspond à lambda_z",
    x        = "Temps (h)",
    y        = "Concentration (ng/mL) — échelle log"
  ) +
  theme_pk

print(p_semilog)
ggsave("nca_semilog.png", plot = p_semilog, width = 8, height = 5, dpi = 300)

cat("\nGraphiques exportés : nca_linear.png et nca_semilog.png\n")
cat("Script terminé.\n")
