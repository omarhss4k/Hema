############################################################
# make_pptx_8slides.R
# PowerPoint 8 slides - T-DXd PK/PD Results
# Encoding: UTF-8 (ASCII-only strings to avoid corruption)
# Generated with officer package
############################################################

options(encoding = "UTF-8")
Sys.setlocale("LC_ALL", "en_US.UTF-8")

library(officer)

setwd("/home/user/Hema/scripts_tdxd")

# ── Couleurs ──────────────────────────────────────────────
col_bleu   <- "#2c7bb6"
col_rouge  <- "#d7191c"
col_vert   <- "#1b9e77"
col_orange <- "#d95f02"
col_gris   <- "#555555"

# ── Repertoire images ─────────────────────────────────────
img_dir <- "results_PKPD_human/pptx_png"

imgs <- list(
  schema       = file.path(img_dir, "schema_modele.png"),
  pk_valid     = file.path(img_dir, "pk_validation.png"),
  pd_typique   = file.path(img_dir, "pd_neut_typique.png"),
  pop_grades   = file.path(img_dir, "population_grades.png"),
  nadir_dist   = file.path(img_dir, "nadir_distribution.png"),
  tableau      = file.path(img_dir, "tableau_params.png")
)

n_imgs <- sum(sapply(imgs, file.exists))
cat(sprintf("Images disponibles : %d / %d\n", n_imgs, length(imgs)))

# ── Presentation ──────────────────────────────────────────
pres <- read_pptx()

# ══════════════════════════════════════════════════════════
# Slide 1 : Titre
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title Slide", master = "Office Theme")
pres <- ph_with(pres,
  value = "Modele PK/PD T-DXd\nSimulation population - Toxicite hematologique",
  location = ph_location_type(type = "ctrTitle"))
pres <- ph_with(pres,
  value = paste0(
    "Validation FDA BLA 761139 | DESTINY-Breast01\n",
    "N = 200 patients simules | Neutropenie & Anemie CTCAE"
  ),
  location = ph_location_type(type = "subTitle"))

# ══════════════════════════════════════════════════════════
# Slide 2 : Schema du modele
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Schema du modele PK/PD T-DXd",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$schema)) {
  pres <- ph_with(pres,
    value = external_img(imgs$schema, width = 9.0, height = 5.5),
    location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
} else {
  schema_txt <- block_list(
    fpar(ftext("PK -- T-DXd (ADC 2 compartiments, Yin 2020)",
               fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
    fpar(ftext("  ADC serum -> DXd plasma -> DXd intracellulaire -> Dommages ADN (Damage)",
               fp_text(font.size = 14))),
    fpar(ftext(" ")),
    fpar(ftext("PD -- Hematopoiese (Fornari 2019)",
               fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
    fpar(ftext("  MPP -> CMP -> Neutrophiles",
               fp_text(font.size = 14))),
    fpar(ftext("  MPP -> MEP -> Erythroblastes",
               fp_text(font.size = 14))),
    fpar(ftext(" ")),
    fpar(ftext("Lien PK->PD : Damage inhibe proliferation des progeniteurs",
               fp_text(font.size = 14, italic = TRUE, color = col_gris)))
  )
  pres <- ph_with(pres, value = schema_txt,
                  location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
}

# ══════════════════════════════════════════════════════════
# Slide 3 : Validation PK
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Validation PK -- ADC & DXd vs FDA BLA 761139",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$pk_valid)) {
  pres <- ph_with(pres,
    value = external_img(imgs$pk_valid, width = 9.0, height = 5.5),
    location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
} else {
  pk_txt <- block_list(
    fpar(ftext("ADC Cmax : 121 ug/mL  (FDA = 122 ug/mL)  [OK]",
               fp_text(font.size = 16, color = col_vert))),
    fpar(ftext("DXd Cmax :   3.8 ng/mL  (FDA =   4.4 ng/mL)  [OK]",
               fp_text(font.size = 16, color = col_vert))),
    fpar(ftext("ADC T1/2 : ~23 jours  (FDA : 5-6 jours mesure, modele pop)",
               fp_text(font.size = 14, color = col_gris)))
  )
  pres <- ph_with(pres, value = pk_txt,
                  location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
}

# ══════════════════════════════════════════════════════════
# Slide 4 : PD patient typique
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Profil PD -- Patient typique (parametres medianes)",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$pd_typique)) {
  pres <- ph_with(pres,
    value = external_img(imgs$pd_typique, width = 9.0, height = 5.5),
    location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
}

# ══════════════════════════════════════════════════════════
# Slide 5 : Simulation population N=100 - Grades
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Simulation population N=200 -- Distribution grades CTCAE",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$pop_grades)) {
  pres <- ph_with(pres,
    value = external_img(imgs$pop_grades, width = 9.0, height = 5.5),
    location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
}

# ══════════════════════════════════════════════════════════
# Slide 6 : Distribution des nadirs
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Distribution des nadirs -- Neutrophiles & Hemoglobine",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$nadir_dist)) {
  pres <- ph_with(pres,
    value = external_img(imgs$nadir_dist, width = 9.0, height = 5.5),
    location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
}

# ══════════════════════════════════════════════════════════
# Slide 7 : Tableau parametres
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Parametres du modele PK/PD -- T-DXd humain",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$tableau)) {
  pres <- ph_with(pres,
    value = external_img(imgs$tableau, width = 9.0, height = 5.5),
    location = ph_location(left = 0.5, top = 1.4, width = 9.0, height = 5.5))
}

# ══════════════════════════════════════════════════════════
# Slide 8 : Synthese et prochaines etapes
# ══════════════════════════════════════════════════════════
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres,
  value = "Synthese et prochaines etapes",
  location = ph_location_type(type = "title"))

synth_txt <- block_list(
  fpar(ftext("[OK]  Validation PK : ADC Cmax = 121 ug/mL (FDA = 122)",
             fp_text(font.size = 16, bold = TRUE, color = col_vert))),
  fpar(ftext("[OK]  Neutropenie G3-4 : 18%   (FDA : 16-20%)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("[OK]  Anemie G3-4      :  7%   (FDA : ~9%)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("[OK]  Nadir Neut median : 1.22 x 10^9/L",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext(" ")),
  fpar(ftext("[!]  Limitation : 0% G0 modele vs 71% clinique",
             fp_text(font.size = 15, bold = TRUE, color = col_orange))),
  fpar(ftext("     Cause : T1/2 ADC (23j) > Q3W (21j) -> accumulation entre cycles",
             fp_text(font.size = 14, color = col_orange))),
  fpar(ftext(" ")),
  fpar(ftext("Prochaines etapes :",
             fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
  fpar(ftext("  1. Augmenter IC50_ADC pour corriger distribution G0",
             fp_text(font.size = 14))),
  fpar(ftext("  2. Simuler dose 6.4 mg/kg Q3W",
             fp_text(font.size = 14))),
  fpar(ftext("  3. Comparaison T-DXd vs T-DM1 (meme modele)",
             fp_text(font.size = 14))),
  fpar(ftext("  4. Adapter structure aux ADC longue demi-vie (Ait-Oudhia 2017)",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Sources : Yin 2020 (PopPK) | Fornari 2019 (PD) | FDA BLA 761139",
             fp_text(font.size = 12, italic = TRUE, color = col_gris)))
)

pres <- ph_with(pres, value = synth_txt,
                location = ph_location_type(type = "body"))

# ── Sauvegarde ────────────────────────────────────────────
out_path <- "results_PKPD_human/TDXd_PKPD_Results.pptx"
print(pres, target = out_path)

n_slides <- length(pres)
cat(sprintf("\nPowerPoint genere : %s\n", out_path))
cat(sprintf("  %d diapositives\n", n_slides))
cat("  Encodage : UTF-8 / ASCII-only strings\n")
cat("  Images incluses :\n")
for (nm in names(imgs)) {
  status <- if (file.exists(imgs[[nm]])) "[OK]" else "[ABSENT]"
  cat(sprintf("    %s  %s\n", status, basename(imgs[[nm]])))
}
