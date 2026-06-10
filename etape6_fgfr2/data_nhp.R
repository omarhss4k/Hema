############################################################
# data_nhp.R
# Données hématologiques NHP observées -- FGFR2 inhibiteur
#
# Remplir nhp_hema_data.csv avec les valeurs réelles.
#
# Groupes / doses :
#   1001, 1002 →  4 mg/kg  (t : -3, 2, 8, 15, 22 j)
#   2001, 2002 → 13 mg/kg  (t : -3, 2, 8, 15, 22 j)
#   4001, 4002 → 26 mg/kg  (t : -3, 2, 8, 15, 22 j)
#   3101, 3002 → 39 mg/kg  (t : -3, 2, 8, 12 j)
#
# Conversions unités → modèle (10⁹/L) :
#   Neut  : 10³/µL  × 1      = 10⁹/L
#   Plt   : 10³/µL  × 1      = 10⁹/L
#   RBC   : 10⁶/µL  × 1000   = 10⁹/L
#   Ret   : 10⁹/L   × 1      = 10⁹/L  (déjà en unités modèle)
############################################################

nhp_raw <- read.csv("nhp_hema_data.csv", stringsAsFactors = FALSE)

obs_data <- data.frame(
  Animal_Id = nhp_raw$Animal_Id,
  dose_mgkg = nhp_raw$Dose_mgkg,
  jour      = nhp_raw$time_day,
  Neut      = nhp_raw$Neut_1e3_uL,
  Plt       = nhp_raw$Plt_1e3_uL,
  RBC       = nhp_raw$RBC_1e6_uL * 1000,
  Ret       = nhp_raw$Ret_1e9_L,
  stringsAsFactors = FALSE
)

# Résumé console
n_total <- nrow(nhp_raw)
cat(sprintf("Données NHP chargées : %d lignes (animaux × temps)\n", n_total))
for (v in c("Neut", "Plt", "RBC", "Ret")) {
  n_ok <- sum(!is.na(obs_data[[v]]))
  cat(sprintf("  %-5s : %d valeurs non-NA\n", v, n_ok))
}
