# install_packages.R
# Installer toutes les dependances du Shiny PK avant le premier lancement
# Lancer UNE SEULE FOIS : source("install_packages.R")

pkgs <- c(
  "shiny", "bslib", "shinyjs",
  "plotly", "ggplot2", "gridExtra",
  "DT", "rhandsontable",
  "readxl", "writexl",
  "dplyr", "tidyr",
  "PKNCA"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  message("Installation de : ", paste(to_install, collapse = ", "))
  install.packages(to_install)
} else {
  message("Tous les packages sont deja installes.")
}

message("\nPour lancer l'application : shiny::runApp('app.R')")
