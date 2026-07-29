############################################################
# make_presentation.R
# PowerPoint - T-DXd PK/PD : Resultats rat + humain + suite
# NOTE: ASCII-only strings to avoid UTF-8 encoding issues
# Lancer depuis le dossier presentation/
############################################################
library(officer)

# Couleurs
col_bleu   <- "#2c7bb6"
col_rouge  <- "#d7191c"
col_vert   <- "#1b9e77"
col_orange <- "#d95f02"
col_gris   <- "#555555"

# Convertir PDFs en PNG pour insertion
convert_pdf <- function(pdf_path, out_png) {
  if (!file.exists(out_png)) {
    cmd <- sprintf("pdftoppm -r 150 -png -l 1 '%s' '%s_tmp' && mv '%s_tmp'-1.png '%s'",
                   pdf_path, out_png, out_png, out_png)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (!file.exists(out_png)) {
      system(sprintf("convert -density 150 '%s'[0] -quality 95 '%s'",
                     pdf_path, out_png),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
  }
  out_png
}

img_dir <- "tmp_pptx_imgs"
dir.create(img_dir, showWarnings = FALSE)

imgs <- list(
  rat_summary  = convert_pdf("../etape3_tdxd_rat/results/rat_validation_FDA_summary.pdf",
                              file.path(img_dir, "rat_summary.png")),
  rat_profils  = convert_pdf("../etape3_tdxd_rat/results/rat_profils_temporels.pdf",
                              file.path(img_dir, "rat_profils.png")),
  hu_neut      = convert_pdf("../etape5_tdxd_humain/results/grades_neutropenie_modele_vs_FDA.pdf",
                              file.path(img_dir, "hu_neut.png")),
  hu_anemie    = convert_pdf("../etape5_tdxd_humain/results/grades_anemie_modele_vs_FDA.pdf",
                              file.path(img_dir, "hu_anemie.png"))
)

n_imgs <- sum(sapply(imgs, file.exists))
cat("Images converties :", n_imgs, "/ 4\n")

# ======================================================
# Construction du PowerPoint
# ======================================================
pres <- read_pptx()

# -- Slide 1 : Titre -----------------------------------
pres <- add_slide(pres, layout = "Title Slide", master = "Office Theme")
pres <- ph_with(pres,
  value = "Modele PK/PD T-DXd\nToxicite hematologique",
  location = ph_location_type(type = "ctrTitle"))
pres <- ph_with(pres,
  value = "Rat (preclinique) -> Humain (clinique)\nDESTINY-Breast01 | FDA BLA 761139",
  location = ph_location_type(type = "subTitle"))

# -- Slide 2 : Plan ------------------------------------
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Plan",
                location = ph_location_type(type = "title"))
plan_txt <- block_list(
  fpar(ftext("1.  Modele rat -- validation vs FDA BLA 761139",
             fp_text(font.size = 20, bold = TRUE, color = col_bleu))),
  fpar(ftext("     Calibration PK + seuils histopathologiques (20 / 60 / 197 mg/kg)",
             fp_text(font.size = 16, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("2.  Modele humain -- validation vs DESTINY-Breast01",
             fp_text(font.size = 20, bold = TRUE, color = col_bleu))),
  fpar(ftext("     Grade 3-4 neutropenie & anemie  |  N = 200 patients simules",
             fp_text(font.size = 16, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("3.  Limitation principale & perspectives",
             fp_text(font.size = 20, bold = TRUE, color = col_bleu))),
  fpar(ftext("     Accumulation ADC entre cycles -> surrepresentation G2",
             fp_text(font.size = 16, color = col_gris)))
)
pres <- ph_with(pres, value = plan_txt,
                location = ph_location_type(type = "body"))

# -- Slide 3 : Structure du modele ---------------------
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Structure du modele PK/PD",
                location = ph_location_type(type = "title"))
struct_txt <- block_list(
  fpar(ftext("PK -- T-DXd (ADC 2 compartiments, Yin 2020)",
             fp_text(font.size = 18, bold = TRUE, color = col_bleu))),
  fpar(ftext("  ADC serum -> DXd plasma -> DXd intracellulaire -> Dommages ADN (Damage)",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("PD -- Hematopoiese (Fornari 2019)",
             fp_text(font.size = 18, bold = TRUE, color = col_bleu))),
  fpar(ftext("  MPP -> CMP -> Neutrophiles / Monocytes",
             fp_text(font.size = 15))),
  fpar(ftext("  MPP -> MEP -> Erythroblastes / Plaquettes",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("Lien PK->PD : Damage inhibe la proliferation des progeniteurs",
             fp_text(font.size = 15, italic = TRUE, color = col_gris))),
  fpar(ftext("  dCMP/dt ~ (1 - Slope_CMP x Damage) x CMP",
             fp_text(font.size = 14, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("Driver toxicite : C_ADC1 [ug/mL]  (assay PFB-10 sur HPC, IC50 = 27.7 ug/mL)",
             fp_text(font.size = 14, italic = TRUE, color = col_orange)))
)
pres <- ph_with(pres, value = struct_txt,
                location = ph_location_type(type = "body"))

# -- Slide 4 : RAT -- validation seuils FDA ------------
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Rat -- Validation vs seuils FDA BLA 761139 (Table 6)",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$rat_summary)) {
  pres <- ph_with(pres,
    value = external_img(imgs$rat_summary, width = 9.5, height = 5.6),
    location = ph_location(left = 0.2, top = 1.4, width = 9.5, height = 5.6))
}

# -- Slide 5 : RAT -- profils temporels ---------------
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Rat -- Profils temporels sur 3 cycles (63 jours)",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$rat_profils)) {
  pres <- ph_with(pres,
    value = external_img(imgs$rat_profils, width = 9.5, height = 5.2),
    location = ph_location(left = 0.2, top = 1.5, width = 9.5, height = 5.2))
}

# -- Slide 6 : RAT -- resume chiffre ------------------
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Rat -- Resultats cles",
                location = ph_location_type(type = "title"))
rat_txt <- block_list(
  fpar(ftext("Calibration Slope_MEP = 1.00  (carboplatin Fornari = 2.19)",
             fp_text(font.size = 16, bold = TRUE))),
  fpar(ftext(" ")),
  fpar(ftext("Dose        Ret (sang)    MEP (moelle)   Neut      FDA histopathologie",
             fp_text(font.size = 14, bold = TRUE, color = col_bleu))),
  fpar(ftext("20 mg/kg    -14.0%        -9.8%          -4.3%     Ret chute >= 20 mg/kg  [OK]",
             fp_text(font.size = 14))),
  fpar(ftext("60 mg/kg    -32.5%        -23.7%         -10.9%    MEP chute >= 60 mg/kg  [OK]",
             fp_text(font.size = 14))),
  fpar(ftext("197 mg/kg   -59.6%        -46.9%         -24.3%    Neut chute >= 197 mg/kg  [OK]",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Les 3 seuils de dose FDA sont reproduits par le modele.",
             fp_text(font.size = 16, bold = TRUE, color = col_vert)))
)
pres <- ph_with(pres, value = rat_txt,
                location = ph_location_type(type = "body"))

# -- Slide 7 : HUMAIN -- neutropenie grades ------------
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Humain -- Grades neutropenie : modele vs DESTINY-Breast01",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$hu_neut)) {
  pres <- ph_with(pres,
    value = external_img(imgs$hu_neut, width = 8.5, height = 5.4),
    location = ph_location(left = 0.7, top = 1.4, width = 8.5, height = 5.4))
}

# -- Slide 8 : HUMAIN -- anemie grades ----------------
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Humain -- Grades anemie : modele vs DESTINY-Breast01",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$hu_anemie)) {
  pres <- ph_with(pres,
    value = external_img(imgs$hu_anemie, width = 8.5, height = 5.4),
    location = ph_location(left = 0.7, top = 1.4, width = 8.5, height = 5.4))
}

# -- Slide 9 : HUMAIN -- resume chiffre ---------------
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Humain -- Resultats cles (N = 200 patients)",
                location = ph_location_type(type = "title"))
hu_txt <- block_list(
  fpar(ftext("PK -- Validation vs FDA BLA 761139 :",
             fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
  fpar(ftext("  ADC Cmax = 121 ug/mL  (FDA = 122)  [OK]",
             fp_text(font.size = 15))),
  fpar(ftext("  DXd Cmax = 3.8 ng/mL  (FDA = 4.4)  [OK]",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("PD -- Toxicite hematologique (Slope_CMP = 12, calibre) :",
             fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Neutropenie G3-4 : 18%   (FDA : 16-20%)  [OK]",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("  Anemie G3-4      :  7%   (FDA :  ~9%)    [OK]",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("  Nadir Neut median : 1.22 x 10^9/L",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("[!]  Tout grade : 100% (modele) vs 29% (FDA)  -> limitation",
             fp_text(font.size = 15, bold = TRUE, color = col_orange)))
)
pres <- ph_with(pres, value = hu_txt,
                location = ph_location_type(type = "body"))

# -- Slide 10 : Limitation -----------------------------
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Limitation : accumulation de l'ADC entre cycles",
                location = ph_location_type(type = "title"))
lim_txt <- block_list(
  fpar(ftext("Cause :",
             fp_text(font.size = 16, bold = TRUE, color = col_rouge))),
  fpar(ftext("  T1/2 terminale ADC ~ 23 jours  >  intervalle Q3W = 21 jours",
             fp_text(font.size = 15))),
  fpar(ftext("  -> Le medicament s'accumule entre les cycles",
             fp_text(font.size = 15))),
  fpar(ftext("  -> Signal Damage jamais nul -> tous les patients ont une chute de Neut",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("Consequence :",
             fp_text(font.size = 16, bold = TRUE, color = col_rouge))),
  fpar(ftext("  0% G0 dans le modele vs 71% en clinique",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("Reference :",
             fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Ait-Oudhia et al., AAPS Journal 2017 (PMID: 28646408)",
             fp_text(font.size = 14, italic = TRUE))),
  fpar(ftext("  -> Modele Fornari/Friberg non adapte aux ADC a longue demi-vie",
             fp_text(font.size = 14))),
  fpar(ftext("  Fornari 2019 : valide uniquement sur le carboplatine (T1/2 ~ 2h)",
             fp_text(font.size = 14)))
)
pres <- ph_with(pres, value = lim_txt,
                location = ph_location_type(type = "body"))

# -- Slide 11 : Perspectives ---------------------------
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Perspectives",
                location = ph_location_type(type = "title"))
suite_txt <- block_list(
  fpar(ftext("Option 1 -- Corriger la distribution des grades",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Augmenter IC50_ADC in vivo (les progeniteurs HPC expriment peu HER2)",
             fp_text(font.size = 14))),
  fpar(ftext("  -> Recalibrer Slope_CMP sur les nouvelles donnees de Damage",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Option 2 -- Explorer d'autres doses",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Simuler 6.4 mg/kg Q3W (dose superieure exploree en essai)",
             fp_text(font.size = 14))),
  fpar(ftext("  Comparer les profils de toxicite entre doses",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Option 3 -- Comparaison T-DXd vs T-DM1",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  T-DM1 (Kadcyla) : meme cible HER2, payload different (DM1 vs DXd)",
             fp_text(font.size = 14))),
  fpar(ftext("  Comparer les profils de myelosuppression avec le meme modele",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Option 4 -- Adapter le modele aux ADC (structure Ait-Oudhia 2017)",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Integrer un compartiment de relargage intratumoral",
             fp_text(font.size = 14)))
)
pres <- ph_with(pres, value = suite_txt,
                location = ph_location_type(type = "body"))

# -- Slide 12 : Conclusion -----------------------------
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Conclusion",
                location = ph_location_type(type = "title"))
concl_txt <- block_list(
  fpar(ftext("[OK]  Modele PK/PD T-DXd construit et calibre sur 2 especes",
             fp_text(font.size = 16, bold = TRUE, color = col_vert))),
  fpar(ftext(" ")),
  fpar(ftext("[OK]  Rat : 3 seuils FDA reproduits (20 / 60 / 197 mg/kg)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("[OK]  Humain : G3-4 neutropenie 18% (FDA 16-20%),  anemie 7% (FDA 9%)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("[OK]  PK validee : ADC Cmax = 121 ug/mL (FDA = 122)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext(" ")),
  fpar(ftext("[!]  Limitation : distribution des grades (G0 sous-estime)",
             fp_text(font.size = 15, bold = TRUE, color = col_orange))),
  fpar(ftext("     Cause : T1/2 ADC > Q3W -> accumulation entre cycles",
             fp_text(font.size = 14, color = col_orange))),
  fpar(ftext("     Reference : Ait-Oudhia et al. AAPS J 2017",
             fp_text(font.size = 14, italic = TRUE, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("Sources : Yin 2020 (PopPK) | Fornari 2019 (PD) | FDA BLA 761139",
             fp_text(font.size = 12, italic = TRUE, color = col_gris)))
)
pres <- ph_with(pres, value = concl_txt,
                location = ph_location_type(type = "body"))

# -- Sauvegarde ----------------------------------------
out_path <- "TDXD_presentation.pptx"
print(pres, target = out_path)
cat(sprintf("\nPowerPoint genere : %s\n", out_path))
cat(sprintf("  %d slides\n", length(pres)))
