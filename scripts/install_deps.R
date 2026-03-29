# install_deps.R
# Installation des packages R nécessaires au projet Hema
# Lancer une seule fois avant la première exécution des scripts

# Packages actuellement utilisés
pkgs_current <- c(
  "deSolve",    # Solveur d'équations différentielles (lsoda)
  "ggplot2",    # Visualisation
  "gridExtra",  # Mise en page des plots (grid.arrange)
  "grid",       # Annotations et layouts
  "scales"      # Transformations d'axes (log10, etc.)
)

# Packages à venir (Phase 2 : estimation des paramètres)
pkgs_future <- c(
  "nlmixr2",    # Estimation NLME en R (alternative à NONMEM)
  "saemix",     # SAEM pour modèles non-linéaires à effets mixtes
  "rxode2"      # Moteur ODE pour nlmixr2
)

# Installation des packages actuels
message("Installation des packages courants...")
install.packages(pkgs_current, repos = "https://cloud.r-project.org")

# Décommenter pour installer les packages futurs (Phase 2)
# message("Installation des packages pour l'estimation (Phase 2)...")
# install.packages(pkgs_future, repos = "https://cloud.r-project.org")

message("Installation terminée.")
message("Vérification :")
for (pkg in pkgs_current) {
  status <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-12s %s\n", pkg, ifelse(status, "OK", "ECHEC")))
}
