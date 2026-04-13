############################################################
# make_presentation.R
# PowerPoint — T-DXd PK/PD : Résultats rat + humain + suite
############################################################
library(officer)

setwd("/home/user/Hema/scripts_tdxd")

# Couleurs
col_bleu   <- "#2c7bb6"
col_rouge  <- "#d7191c"
col_vert   <- "#1b9e77"
col_orange <- "#d95f02"
col_gris   <- "#555555"
col_fond   <- "#f7f9fc"

# Helper : ajouter un slide avec titre + image PDF convertie en PNG
add_slide_img <- function(pres, layout, title, img_path,
                          subtitle = NULL,
                          img_left = 0.4, img_top = 1.5,
                          img_w = 9.2, img_h = 5.4) {

  pres <- add_slide(pres, layout = layout, master = "Office Theme")

  # Titre
  pres <- ph_with(pres, value = title,
                  location = ph_location_type(type = "title"))

  # Sous-titre / corps si fourni
  if (!is.null(subtitle)) {
    pres <- ph_with(pres, value = subtitle,
                    location = ph_location_type(type = "body"))
  }

  # Image
  if (!is.null(img_path) && file.exists(img_path)) {
    pres <- ph_with(pres,
                    value = external_img(img_path, width = img_w, height = img_h),
                    location = ph_location(left = img_left, top = img_top,
                                           width = img_w, height = img_h))
  }
  pres
}

# Convertir PDFs en PNG pour insertion
convert_pdf <- function(pdf_path, out_png, width = 1400) {
  if (!file.exists(out_png)) {
    cmd <- sprintf("pdftoppm -r 150 -png -l 1 '%s' '%s_tmp' && mv '%s_tmp'-1.png '%s'",
                   pdf_path, out_png, out_png, out_png)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    # Fallback avec convert (ImageMagick)
    if (!file.exists(out_png)) {
      system(sprintf("convert -density 150 '%s'[0] -quality 95 '%s'",
                     pdf_path, out_png),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
  }
  out_png
}

# Préparer les images PNG
img_dir <- "tmp_pptx_imgs"
dir.create(img_dir, showWarnings = FALSE)

imgs <- list(
  rat_summary  = convert_pdf("results_PKPD/rat_validation_FDA_summary.pdf",
                              file.path(img_dir, "rat_summary.png")),
  rat_profils  = convert_pdf("results_PKPD/rat_profils_temporels.pdf",
                              file.path(img_dir, "rat_profils.png")),
  hu_neut      = convert_pdf("results_PKPD_human/grades_neutropenie_modele_vs_FDA.pdf",
                              file.path(img_dir, "hu_neut.png")),
  hu_anemie    = convert_pdf("results_PKPD_human/grades_anemie_modele_vs_FDA.pdf",
                              file.path(img_dir, "hu_anemie.png"))
)

cat("Images converties :", sum(sapply(imgs, file.exists)), "/ 4\n")

# ══════════════════════════════════════════════════════════
# Construction du PowerPoint
# ══════════════════════════════════════════════════════════
pres <- read_pptx()

# ── Slide 1 : Titre ──────────────────────────────────────
pres <- add_slide(pres, layout = "Title Slide", master = "Office Theme")
pres <- ph_with(pres,
  value = "Modèle PK/PD T-DXd\nToxicité hématologique",
  location = ph_location_type(type = "ctrTitle"))
pres <- ph_with(pres,
  value = "Rat (préclinique) → Humain (clinique)\nDESTINY-Breast01 | FDA BLA 761139",
  location = ph_location_type(type = "subTitle"))

# ── Slide 2 : Plan ────────────────────────────────────────
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Plan",
                location = ph_location_type(type = "title"))
plan_txt <- block_list(
  fpar(ftext("1.  Modèle rat — validation vs FDA BLA 761139",
             fp_text(font.size = 20, bold = TRUE, color = col_bleu))),
  fpar(ftext("     Calibration PK + seuils histopathologiques (20 / 60 / 197 mg/kg)",
             fp_text(font.size = 16, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("2.  Modèle humain — validation vs DESTINY-Breast01",
             fp_text(font.size = 20, bold = TRUE, color = col_bleu))),
  fpar(ftext("     Grade 3-4 neutropénie & anémie  |  N = 200 patients simulés",
             fp_text(font.size = 16, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("3.  Limitation principale & perspectives",
             fp_text(font.size = 20, bold = TRUE, color = col_bleu))),
  fpar(ftext("     Accumulation ADC entre cycles → surreprésentation G2",
             fp_text(font.size = 16, color = col_gris)))
)
pres <- ph_with(pres, value = plan_txt,
                location = ph_location_type(type = "body"))

# ── Slide 3 : Modèle (structure) ─────────────────────────
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Structure du modèle PK/PD",
                location = ph_location_type(type = "title"))
struct_txt <- block_list(
  fpar(ftext("PK — T-DXd (ADC 2 compartiments, Yin 2020)",
             fp_text(font.size = 18, bold = TRUE, color = col_bleu))),
  fpar(ftext("  ADC sérum  →  DXd plasma  →  DXd intracellulaire  →  Dommages ADN (Damage)",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("PD — Hématopoïèse (Fornari 2019)",
             fp_text(font.size = 18, bold = TRUE, color = col_bleu))),
  fpar(ftext("  MPP  →  CMP  →  Neutrophiles / Monocytes",
             fp_text(font.size = 15))),
  fpar(ftext("  MPP  →  MEP  →  Érythroblastes / Plaquettes",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("Lien PK→PD : Damage inhibe la prolifération des progéniteurs",
             fp_text(font.size = 15, italic = TRUE, color = col_gris))),
  fpar(ftext("  dCMP/dt ∝ (1 − Slope_CMP × Damage) × CMP",
             fp_text(font.size = 14, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("Driver toxicité : C_ADC1 [µg/mL]  (assay PFB-10 sur HPC, IC50 = 27.7 µg/mL)",
             fp_text(font.size = 14, italic = TRUE, color = col_orange)))
)
pres <- ph_with(pres, value = struct_txt,
                location = ph_location_type(type = "body"))

# ── Slide 4 : RAT — validation seuils FDA ────────────────
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Rat — Validation vs seuils FDA BLA 761139 (Table 6)",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$rat_summary)) {
  pres <- ph_with(pres,
    value = external_img(imgs$rat_summary, width = 9.5, height = 5.6),
    location = ph_location(left = 0.2, top = 1.4, width = 9.5, height = 5.6))
}

# ── Slide 5 : RAT — profils temporels ────────────────────
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Rat — Profils temporels sur 3 cycles (63 jours)",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$rat_profils)) {
  pres <- ph_with(pres,
    value = external_img(imgs$rat_profils, width = 9.5, height = 5.2),
    location = ph_location(left = 0.2, top = 1.5, width = 9.5, height = 5.2))
}

# ── Slide 6 : RAT — résumé chiffré ───────────────────────
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Rat — Résultats clés",
                location = ph_location_type(type = "title"))
rat_txt <- block_list(
  fpar(ftext("Calibration Slope_MEP = 1.00  (carboplatin Fornari = 2.19)",
             fp_text(font.size = 16, bold = TRUE))),
  fpar(ftext(" ")),
  fpar(ftext("Dose        Ret (sang)    MEP (moelle)   Neut      FDA histopathologie",
             fp_text(font.size = 14, bold = TRUE, color = col_bleu))),
  fpar(ftext("20 mg/kg    −14.0%        −9.8%          −4.3%     ↓ Ret ≥ 20 mg/kg  ✓",
             fp_text(font.size = 14))),
  fpar(ftext("60 mg/kg    −32.5%        −23.7%         −10.9%    ↓ MEP ≥ 60 mg/kg  ✓",
             fp_text(font.size = 14))),
  fpar(ftext("197 mg/kg   −59.6%        −46.9%         −24.3%    ↓ Neut ≥ 197 mg/kg ✓",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Les 3 seuils de dose FDA sont reproduits par le modèle.",
             fp_text(font.size = 16, bold = TRUE, color = col_vert)))
)
pres <- ph_with(pres, value = rat_txt,
                location = ph_location_type(type = "body"))

# ── Slide 7 : HUMAIN — neutropénie grades ────────────────
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Humain — Grades neutropénie : modèle vs DESTINY-Breast01",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$hu_neut)) {
  pres <- ph_with(pres,
    value = external_img(imgs$hu_neut, width = 8.5, height = 5.4),
    location = ph_location(left = 0.7, top = 1.4, width = 8.5, height = 5.4))
}

# ── Slide 8 : HUMAIN — anémie grades ─────────────────────
pres <- add_slide(pres, layout = "Title Only", master = "Office Theme")
pres <- ph_with(pres,
  value = "Humain — Grades anémie : modèle vs DESTINY-Breast01",
  location = ph_location_type(type = "title"))
if (file.exists(imgs$hu_anemie)) {
  pres <- ph_with(pres,
    value = external_img(imgs$hu_anemie, width = 8.5, height = 5.4),
    location = ph_location(left = 0.7, top = 1.4, width = 8.5, height = 5.4))
}

# ── Slide 9 : HUMAIN — résumé chiffré ────────────────────
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Humain — Résultats clés (N = 200 patients)",
                location = ph_location_type(type = "title"))
hu_txt <- block_list(
  fpar(ftext("PK — Validation vs FDA BLA 761139 :",
             fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
  fpar(ftext("  ADC Cmax = 121 µg/mL  (FDA = 122)  ✓",
             fp_text(font.size = 15))),
  fpar(ftext("  DXd Cmax = 3.8 ng/mL  (FDA = 4.4)  ✓",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("PD — Toxicité hématologique (Slope_CMP = 12, calibré) :",
             fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Neutropénie G3-4 : 18%   (FDA : 16-20%)  ✓",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("  Anémie G3-4      :  7%   (FDA :  ~9%)    ✓",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("  Nadir Neut médian : 1.22 × 10⁹/L",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("⚠  Tout grade : 100% (modèle) vs 29% (FDA)  → limitation",
             fp_text(font.size = 15, bold = TRUE, color = col_orange)))
)
pres <- ph_with(pres, value = hu_txt,
                location = ph_location_type(type = "body"))

# ── Slide 10 : Limitation ────────────────────────────────
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Limitation principale : accumulation de l'ADC",
                location = ph_location_type(type = "title"))
lim_txt <- block_list(
  fpar(ftext("Cause :",
             fp_text(font.size = 16, bold = TRUE, color = col_rouge))),
  fpar(ftext("  T½ terminale ADC ≈ 23 jours  >  intervalle Q3W = 21 jours",
             fp_text(font.size = 15))),
  fpar(ftext("  → Le médicament s'accumule entre les cycles",
             fp_text(font.size = 15))),
  fpar(ftext("  → Signal Damage jamais nul → tous les patients ont une chute de Neut",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("Conséquence :",
             fp_text(font.size = 16, bold = TRUE, color = col_rouge))),
  fpar(ftext("  0% G0 dans le modèle vs 71% en clinique",
             fp_text(font.size = 15))),
  fpar(ftext(" ")),
  fpar(ftext("Référence :",
             fp_text(font.size = 16, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Ait-Oudhia et al., AAPS Journal 2017 (PMID: 28646408)",
             fp_text(font.size = 14, italic = TRUE))),
  fpar(ftext("  → Modèle Fornari/Friberg non adapté aux ADC à longue demi-vie",
             fp_text(font.size = 14))),
  fpar(ftext("  Fornari 2019 : validé uniquement sur le carboplatine (T½ ≈ 2h)",
             fp_text(font.size = 14)))
)
pres <- ph_with(pres, value = lim_txt,
                location = ph_location_type(type = "body"))

# ── Slide 11 : Suite / Perspectives ──────────────────────
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Perspectives",
                location = ph_location_type(type = "title"))
suite_txt <- block_list(
  fpar(ftext("Option 1 — Corriger la distribution des grades",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Augmenter IC50_ADC in vivo (les progéniteurs HPC expriment peu HER2)",
             fp_text(font.size = 14))),
  fpar(ftext("  → Recalibrer Slope_CMP sur les nouvelles données de Damage",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Option 2 — Explorer d'autres doses",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Simuler 6.4 mg/kg Q3W (dose supérieure explorée en essai)",
             fp_text(font.size = 14))),
  fpar(ftext("  Comparer les profils de toxicité entre doses",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Option 3 — Comparaison T-DXd vs T-DM1",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  T-DM1 (Kadcyla) : même cible HER2, payload différent (DM1 vs DXd)",
             fp_text(font.size = 14))),
  fpar(ftext("  Comparer les profils de myélosuppression avec le même modèle",
             fp_text(font.size = 14))),
  fpar(ftext(" ")),
  fpar(ftext("Option 4 — Adapter le modèle aux ADC (structure Ait-Oudhia 2017)",
             fp_text(font.size = 17, bold = TRUE, color = col_bleu))),
  fpar(ftext("  Intégrer un compartiment de relargage intratumoral",
             fp_text(font.size = 14)))
)
pres <- ph_with(pres, value = suite_txt,
                location = ph_location_type(type = "body"))

# ── Slide 12 : Conclusion ─────────────────────────────────
pres <- add_slide(pres, layout = "Title and Content", master = "Office Theme")
pres <- ph_with(pres, value = "Conclusion",
                location = ph_location_type(type = "title"))
concl_txt <- block_list(
  fpar(ftext("✓  Modèle PK/PD T-DXd construit et calibré sur 2 espèces",
             fp_text(font.size = 16, bold = TRUE, color = col_vert))),
  fpar(ftext(" ")),
  fpar(ftext("✓  Rat : 3 seuils FDA reproduits (20 / 60 / 197 mg/kg)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("✓  Humain : G3-4 neutropénie 18% (FDA 16-20%),  anémie 7% (FDA 9%)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext("✓  PK validée : ADC Cmax = 121 µg/mL (FDA = 122)",
             fp_text(font.size = 15, color = col_vert))),
  fpar(ftext(" ")),
  fpar(ftext("⚠  Limitation : distribution des grades (G0 sous-estimé)",
             fp_text(font.size = 15, bold = TRUE, color = col_orange))),
  fpar(ftext("     Cause : T½ ADC > Q3W → accumulation entre cycles",
             fp_text(font.size = 14, color = col_orange))),
  fpar(ftext("     Documenté : Ait-Oudhia et al. AAPS J 2017",
             fp_text(font.size = 14, italic = TRUE, color = col_gris))),
  fpar(ftext(" ")),
  fpar(ftext("Sources : Yin 2020 (PopPK) | Fornari 2019 (PD) | FDA BLA 761139",
             fp_text(font.size = 12, italic = TRUE, color = col_gris)))
)
pres <- ph_with(pres, value = concl_txt,
                location = ph_location_type(type = "body"))

# ── Sauvegarde ────────────────────────────────────────────
out_path <- "TDXD_presentation.pptx"
print(pres, target = out_path)
cat(sprintf("\nPowerPoint généré : %s\n", out_path))
cat(sprintf("  %d slides\n", length(pres)))
