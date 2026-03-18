############################################################
# data_fornari.R
# Données digitalisées depuis Figure 3 — Fornari 2019
# Carboplatin 40 mg/kg Q14D × 8 cycles
# x = jours, y = 10⁹ cells/L
############################################################

read_dig <- function(path) {
  df <- read.csv(path, header = TRUE, strip.white = TRUE)
  colnames(df) <- c("time", "value")
  df <- df[order(df$time), ]
  df
}

MPP_obs  <- read_dig("MPP.csv")
CMP_obs  <- read_dig("CMP.csv")

# MEP : deux séries digitalisées (symboles carrés + losanges Fig. 3)
MEP_a    <- read_dig("MEP.csv")
MEP_b    <- read_dig("MEPF.csv")
MEP_obs  <- rbind(MEP_a, MEP_b)
MEP_obs  <- MEP_obs[order(MEP_obs$time), ]

Neut_obs <- read_dig("Neut.csv")
Mono_obs <- read_dig("MONO.csv")
Plt_obs  <- read_dig("Plt.csv")
Ret_obs  <- read_dig("Reut.csv")
RBC_obs  <- read_dig("RBC.csv")

obs_fornari <- list(
  MPP  = MPP_obs,
  CMP  = CMP_obs,
  MEP  = MEP_obs,
  Neut = Neut_obs,
  Mono = Mono_obs,
  Plt  = Plt_obs,
  Ret  = Ret_obs,
  RBC  = RBC_obs
)

cat(sprintf("  Donnees Fornari chargees : %d types cellulaires\n",
            length(obs_fornari)))
for (nm in names(obs_fornari))
  cat(sprintf("    %-5s : %d points  [%.1f – %.1f j]\n",
              nm,
              nrow(obs_fornari[[nm]]),
              min(obs_fornari[[nm]]$time),
              max(obs_fornari[[nm]]$time)))
