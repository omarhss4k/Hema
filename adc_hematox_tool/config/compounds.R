############################################################
# config/compounds.R  --  LE SEUL FICHIER A EDITER pour un usage courant.
# Definis ici : l'ADC de reference, la calibration PD, et le(s)
# nouveau(x) composE(s) a predire.
############################################################

## ═══════ ADC DE REFERENCE (donnees completes : PK + IC50 + profils) ═══════
## Sert a CALIBRER le modele. Ses profils hemato vont dans data/.
REFERENCE <- list(
  name       = "ADC cytotoxique",
  BW_kg      = 2.5,
  # -- PK 2-cmt (unites : L, L/h) --
  CL = 0.00484307629862215, V1 = 0.04907800106132178,
  Q  = 0.01627593959372783, V2 = 0.1170226727748038,
  Tinfu_h    = 1.5,
  IC50_myelo = 0.008615,          # nM  (IC50 myeloide in vitro)
  doses      = c(0.3, 1, 3),      # mg/kg
  profiles   = "data/reference_profiles.csv"   # CSV confidentiel (voir README)
)

## ═══════ CALIBRATION PD (calee sur REFERENCE) ═══════
## Ces valeurs sont le resultat du calage. Pour re-caler : ajuste-les et
## relance scripts/calibrate.R (qui les re-sauvegarde dans results/).
CALIB <- list(
  k_rep         = 0.002,   # reparation lente -> nadir tardif (j16-22)
  ED50_kill     = 0.1,
  Emax_CMP      = 0.982,   # blocage proliferation myeloide (quasi-max)
  k_depl_direct = 0.25,    # deplation directe progeniteurs -> atteint le G4
  Emax_MEP      = 1e-4,    # erythroide EPARGNE
  IC50ADC_scale = 27.70    # ug/mL : echelle de potency du modele a cette calibration
)

## ═══════ NOUVEAU COMPOSE a predire ═══════
## Remplis PK + IC50. Le modele predit sa toxicite SANS donnees in vivo.
NEW_COMPOUND <- list(
  name       = "BPA",
  BW_kg      = 2.5,
  # -- PK 2-cmt du composE (depuis ton Shiny) --
  CL = 0.000467305265769336, V1 = 0.0586062942310418,
  Q  = 0.00102316504406689,  V2 = 0.0424186702071202,
  Tinfu_h    = 1.5,
  IC50_myelo = c(0.015, 7.8),     # nM -- une ou plusieurs valeurs (fourchette)
  doses      = c(0.3, 1, 3),      # mg/kg
  # baselines absolues (x10^9/L) pour convertir en grade CTCAE :
  baseline_neut = 2.4, baseline_mono = 0.8,
  # note mecanistique (affichee dans la sortie) :
  note = "Si mecanisme immuno/STING : prediction cytotoxique = PLANCHER."
)
