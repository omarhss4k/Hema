#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_pptx_8slides.py
PowerPoint 8 slides - T-DXd PK/PD Results
Encodage : UTF-8 natif Python (aucun probleme d'accents)

8 diapositives :
  1. Titre
  2. Schema du modele
  3. Validation PK
  4. PD patient typique
  5. Simulation population (grades CTCAE)
  6. Distribution des nadirs
  7. Tableau parametres
  8. Synthese et prochaines etapes
"""

import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

# ── Couleurs ──────────────────────────────────────────────
BLEU   = RGBColor(0x2c, 0x7b, 0xb6)
ROUGE  = RGBColor(0xd7, 0x19, 0x1c)
VERT   = RGBColor(0x1b, 0x9e, 0x77)
ORANGE = RGBColor(0xd9, 0x5f, 0x02)
GRIS   = RGBColor(0x55, 0x55, 0x55)
NOIR   = RGBColor(0x00, 0x00, 0x00)
BLANC  = RGBColor(0xFF, 0xFF, 0xFF)

# ── Chemins ───────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
IMG_DIR    = os.path.join(SCRIPT_DIR, "results_PKPD_human", "pptx_png")
OUT_PATH   = os.path.join(SCRIPT_DIR, "results_PKPD_human", "TDXd_PKPD_Results.pptx")

IMGS = {
    "schema":     os.path.join(IMG_DIR, "schema_modele.png"),
    "pk_valid":   os.path.join(IMG_DIR, "pk_validation.png"),
    "pd_typique": os.path.join(IMG_DIR, "pd_neut_typique.png"),
    "pop_grades": os.path.join(IMG_DIR, "population_grades.png"),
    "nadir_dist": os.path.join(IMG_DIR, "nadir_distribution.png"),
    "tableau":    os.path.join(IMG_DIR, "tableau_params.png"),
}

# ── Dimensions (16:9) ─────────────────────────────────────
SLIDE_W = Inches(10)
SLIDE_H = Inches(7.5)


def new_prs():
    """Cree une presentation 16:9 vierge."""
    prs = Presentation()
    prs.slide_width  = SLIDE_W
    prs.slide_height = SLIDE_H
    return prs


def blank_layout(prs):
    return prs.slide_layouts[6]   # Blank


def add_title_box(slide, text, top=Inches(0.3), color=BLEU, size=Pt(28), bold=True):
    """Ajoute un titre en haut du slide."""
    txBox = slide.shapes.add_textbox(Inches(0.3), top, Inches(9.4), Inches(0.9))
    tf = txBox.text_frame
    tf.word_wrap = True
    p  = tf.paragraphs[0]
    run = p.add_run()
    run.text = text
    run.font.size  = size
    run.font.bold  = bold
    run.font.color.rgb = color
    return txBox


def add_image(slide, img_path, left=Inches(0.5), top=Inches(1.3),
              width=Inches(9.0), height=Inches(5.8)):
    """Ajoute une image si elle existe."""
    if os.path.exists(img_path):
        slide.shapes.add_picture(img_path, left, top, width, height)
        return True
    return False


def add_bullet_box(slide, lines, top=Inches(1.3), height=Inches(5.8)):
    """Ajoute une zone de texte multi-lignes.

    lines : liste de (texte, taille, gras, couleur)
    """
    txBox = slide.shapes.add_textbox(Inches(0.4), top, Inches(9.2), height)
    tf = txBox.text_frame
    tf.word_wrap = True
    first = True
    for (text, size, bold, color) in lines:
        if first:
            p = tf.paragraphs[0]
            first = False
        else:
            p = tf.add_paragraph()
        p.space_before = Pt(2)
        run = p.add_run()
        run.text = text
        run.font.size  = Pt(size)
        run.font.bold  = bold
        run.font.color.rgb = color


# ══════════════════════════════════════════════════════════
# Slide 1 : Titre
# ══════════════════════════════════════════════════════════
def slide_titre(prs):
    slide = prs.slides.add_slide(blank_layout(prs))

    # Fond bleu fonce
    from pptx.util import Emu
    from pptx.dml.color import RGBColor
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = RGBColor(0x1a, 0x3a, 0x5c)

    # Titre principal
    tb1 = slide.shapes.add_textbox(Inches(0.6), Inches(2.0), Inches(8.8), Inches(1.8))
    tf1 = tb1.text_frame
    tf1.word_wrap = True
    p1 = tf1.paragraphs[0]
    p1.alignment = PP_ALIGN.CENTER
    r1 = p1.add_run()
    r1.text = "Modele PK/PD T-DXd"
    r1.font.size  = Pt(36)
    r1.font.bold  = True
    r1.font.color.rgb = BLANC

    p1b = tf1.add_paragraph()
    p1b.alignment = PP_ALIGN.CENTER
    r1b = p1b.add_run()
    r1b.text = "Simulation population - Toxicite hematologique"
    r1b.font.size  = Pt(24)
    r1b.font.bold  = False
    r1b.font.color.rgb = RGBColor(0xaa, 0xcc, 0xff)

    # Sous-titre
    tb2 = slide.shapes.add_textbox(Inches(0.6), Inches(4.2), Inches(8.8), Inches(1.2))
    tf2 = tb2.text_frame
    tf2.word_wrap = True
    p2 = tf2.paragraphs[0]
    p2.alignment = PP_ALIGN.CENTER
    r2 = p2.add_run()
    r2.text = "Validation FDA BLA 761139 | DESTINY-Breast01"
    r2.font.size  = Pt(18)
    r2.font.color.rgb = RGBColor(0xdd, 0xee, 0xff)

    p2b = tf2.add_paragraph()
    p2b.alignment = PP_ALIGN.CENTER
    r2b = p2b.add_run()
    r2b.text = "N = 200 patients simules | Neutropenie & Anemie CTCAE"
    r2b.font.size  = Pt(16)
    r2b.font.color.rgb = RGBColor(0xdd, 0xee, 0xff)


# ══════════════════════════════════════════════════════════
# Slide 2 : Schema du modele
# ══════════════════════════════════════════════════════════
def slide_schema(prs):
    slide = prs.slides.add_slide(blank_layout(prs))
    add_title_box(slide, "Schema du modele PK/PD T-DXd")
    if not add_image(slide, IMGS["schema"]):
        # Fallback texte
        lines = [
            ("PK -- T-DXd (ADC 2 compartiments, Yin 2020)", 17, True, BLEU),
            ("  ADC serum -> DXd plasma -> DXd intracellulaire -> Dommages ADN", 14, False, NOIR),
            ("", 8, False, NOIR),
            ("PD -- Hematopoiese (Fornari 2019)", 17, True, BLEU),
            ("  MPP -> CMP -> Neutrophiles / Monocytes", 14, False, NOIR),
            ("  MPP -> MEP -> Erythroblastes / Plaquettes", 14, False, NOIR),
            ("", 8, False, NOIR),
            ("Lien PK->PD : Damage inhibe proliferation des progeniteurs", 14, True, GRIS),
            ("  dCMP/dt ~ (1 - Slope_CMP x Damage) x CMP", 13, False, GRIS),
            ("", 8, False, NOIR),
            ("Driver toxicite : C_ADC1 [ug/mL]  (IC50 = 27.7 ug/mL)", 13, True, ORANGE),
        ]
        add_bullet_box(slide, lines)


# ══════════════════════════════════════════════════════════
# Slide 3 : Validation PK
# ══════════════════════════════════════════════════════════
def slide_pk(prs):
    slide = prs.slides.add_slide(blank_layout(prs))
    add_title_box(slide, "Validation PK -- ADC & DXd vs FDA BLA 761139")
    if not add_image(slide, IMGS["pk_valid"]):
        lines = [
            ("ADC Cmax : 121 ug/mL   (FDA = 122 ug/mL)  [OK]", 18, True, VERT),
            ("DXd Cmax :   3.8 ng/mL  (FDA =   4.4 ng/mL)  [OK]", 18, True, VERT),
            ("", 8, False, NOIR),
            ("ADC T1/2 ~ 23 jours  (Yin 2020 PopPK)", 15, False, GRIS),
            ("DXd T1/2 ~  5.8 jours", 15, False, GRIS),
            ("", 8, False, NOIR),
            ("Protocole : dose 5.4 mg/kg IV Q3W x 6 cycles", 14, False, NOIR),
            ("Modele PopPK : 2 compartiments ADC + 1 compartiment DXd", 14, False, NOIR),
        ]
        add_bullet_box(slide, lines)


# ══════════════════════════════════════════════════════════
# Slide 4 : PD patient typique
# ══════════════════════════════════════════════════════════
def slide_pd_typique(prs):
    slide = prs.slides.add_slide(blank_layout(prs))
    add_title_box(slide, "Profil PD -- Patient typique (parametres medianes)")
    add_image(slide, IMGS["pd_typique"])


# ══════════════════════════════════════════════════════════
# Slide 5 : Simulation population - Grades CTCAE
# ══════════════════════════════════════════════════════════
def slide_pop_grades(prs):
    slide = prs.slides.add_slide(blank_layout(prs))
    add_title_box(slide, "Simulation population N=200 -- Grades CTCAE neutropenie & anemie")
    add_image(slide, IMGS["pop_grades"])


# ══════════════════════════════════════════════════════════
# Slide 6 : Distribution des nadirs
# ══════════════════════════════════════════════════════════
def slide_nadirs(prs):
    slide = prs.slides.add_slide(blank_layout(prs))
    add_title_box(slide, "Distribution des nadirs -- Neutrophiles & Hemoglobine")
    add_image(slide, IMGS["nadir_dist"])


# ══════════════════════════════════════════════════════════
# Slide 7 : Tableau parametres
# ══════════════════════════════════════════════════════════
def slide_tableau(prs):
    slide = prs.slides.add_slide(blank_layout(prs))
    add_title_box(slide, "Parametres du modele PK/PD T-DXd -- humain (Yin 2020 + calibration)")
    if not add_image(slide, IMGS["tableau"]):
        lines = [
            ("Parametres PK (Yin 2020 population typique)", 16, True, BLEU),
            ("  CL     = 0.338 L/h       V1  = 3.13 L", 14, False, NOIR),
            ("  Q      = 0.198 L/h       V2  = 2.78 L", 14, False, NOIR),
            ("  Krel   = 0.105 h-1  (ADC -> DXd)", 14, False, NOIR),
            ("  CL_DXd = 17.8  L/h       Vd_DXd = 172 L", 14, False, NOIR),
            ("", 8, False, NOIR),
            ("Parametres PD (Fornari 2019, calibres T-DXd)", 16, True, BLEU),
            ("  Slope_CMP = 12.0   (calibre sur DESTINY-Breast01 Neut G3-4 = 16-20%)", 14, False, NOIR),
            ("  Slope_MEP =  1.0   (calibre sur FDA BLA 761139 rat, seuil Ret -20%)", 14, False, NOIR),
            ("  IC50_ADC  = 27.7 ug/mL  (assay PFB-10 HPC)", 14, False, NOIR),
        ]
        add_bullet_box(slide, lines)


# ══════════════════════════════════════════════════════════
# Slide 8 : Synthese et prochaines etapes
# ══════════════════════════════════════════════════════════
def slide_synthese(prs):
    slide = prs.slides.add_slide(blank_layout(prs))
    add_title_box(slide, "Synthese et prochaines etapes")

    lines = [
        ("[OK]  Validation PK : ADC Cmax = 121 ug/mL  (FDA = 122)", 15, True,  VERT),
        ("[OK]  Neutropenie G3-4 : 18%   (FDA : 16-20%)", 15, False, VERT),
        ("[OK]  Anemie G3-4      :  7%   (FDA : ~9%)", 15, False, VERT),
        ("[OK]  Nadir Neut median : 1.22 x 10^9/L", 15, False, VERT),
        ("", 6, False, NOIR),
        ("[!]  Limitation : 0% G0 modele vs 71% clinique", 15, True,  ORANGE),
        ("     T1/2 ADC (23j) > Q3W (21j) -> accumulation entre cycles", 14, False, ORANGE),
        ("", 6, False, NOIR),
        ("Prochaines etapes :", 15, True, BLEU),
        ("  1.  Augmenter IC50_ADC pour corriger la distribution G0", 14, False, NOIR),
        ("  2.  Simuler 6.4 mg/kg Q3W (dose superieure exploree en essai)", 14, False, NOIR),
        ("  3.  Comparaison T-DXd vs T-DM1 (meme modele, payloads differents)", 14, False, NOIR),
        ("  4.  Adapter structure aux ADC longue demi-vie (Ait-Oudhia 2017)", 14, False, NOIR),
        ("", 6, False, NOIR),
        ("Sources : Yin 2020 (PopPK) | Fornari 2019 (PD) | FDA BLA 761139", 11, True, GRIS),
    ]
    add_bullet_box(slide, lines, top=Inches(1.2), height=Inches(6.0))


# ══════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════
def main():
    # Verification images
    n_ok = sum(os.path.exists(p) for p in IMGS.values())
    print(f"Images disponibles : {n_ok} / {len(IMGS)}")
    for name, path in IMGS.items():
        status = "[OK]" if os.path.exists(path) else "[ABSENT]"
        print(f"  {status}  {os.path.basename(path)}")

    prs = new_prs()

    slide_titre(prs)
    slide_schema(prs)
    slide_pk(prs)
    slide_pd_typique(prs)
    slide_pop_grades(prs)
    slide_nadirs(prs)
    slide_tableau(prs)
    slide_synthese(prs)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    prs.save(OUT_PATH)

    print(f"\nPowerPoint genere : {OUT_PATH}")
    print(f"  {len(prs.slides)} diapositives")
    print("  Encodage : UTF-8 Python natif (pas de corruption d'accents)")


if __name__ == "__main__":
    main()
