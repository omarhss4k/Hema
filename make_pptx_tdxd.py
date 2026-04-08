"""
make_pptx_tdxd.py
Generates Point_Stage_TDXd.pptx — internship meeting presentation on T-DXd PK/PD modeling.
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import os

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DARK_BLUE   = RGBColor(0x1F, 0x38, 0x64)
MID_BLUE    = RGBColor(0x2E, 0x54, 0x9A)
ACCENT_BLUE = RGBColor(0x4A, 0x86, 0xC8)
WHITE       = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT_GREY  = RGBColor(0xF2, 0xF2, 0xF2)
DARK_GREY   = RGBColor(0x40, 0x40, 0x40)
BLACK       = RGBColor(0x00, 0x00, 0x00)
TABLE_HDR   = RGBColor(0x1F, 0x38, 0x64)
TABLE_ALT   = RGBColor(0xD9, 0xE2, 0xF3)

IMG_DIR = "/home/user/Hema/scripts_tdxd/results_PKPD"
OUT_PATH = "/home/user/Hema/Point_Stage_TDXd.pptx"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def new_prs():
    prs = Presentation()
    prs.slide_width  = Inches(13.33)
    prs.slide_height = Inches(7.5)
    return prs


def blank_slide(prs):
    blank_layout = prs.slide_layouts[6]  # completely blank
    return prs.slides.add_slide(blank_layout)


def fill_shape(shape, color):
    shape.fill.solid()
    shape.fill.fore_color.rgb = color


def add_rect(slide, left, top, width, height, color):
    shape = slide.shapes.add_shape(
        1,  # MSO_SHAPE_TYPE.RECTANGLE
        Inches(left), Inches(top), Inches(width), Inches(height)
    )
    fill_shape(shape, color)
    shape.line.fill.background()
    return shape


def add_textbox(slide, text, left, top, width, height,
                font_name="Calibri", font_size=18, bold=False, italic=False,
                color=BLACK, align=PP_ALIGN.LEFT, word_wrap=True):
    txBox = slide.shapes.add_textbox(
        Inches(left), Inches(top), Inches(width), Inches(height)
    )
    tf = txBox.text_frame
    tf.word_wrap = word_wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.name = font_name
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return txBox


def add_bullet_para(tf, text, level=0,
                    font_name="Calibri", font_size=16,
                    bold=False, color=DARK_GREY, space_before=6):
    from pptx.util import Pt
    from pptx.oxml.ns import qn
    from lxml import etree

    p = tf.add_paragraph()
    p.level = level
    p.space_before = Pt(space_before)

    # bullet character
    pPr = p._p.get_or_add_pPr()
    buChar = etree.SubElement(pPr, qn('a:buChar'))
    buChar.set('char', '•')

    run = p.add_run()
    run.text = text
    run.font.name = font_name
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.color.rgb = color
    return p


def add_title_bar(slide, title_text, bar_height=1.1):
    """Dark-blue title bar across the top."""
    bar = add_rect(slide, 0, 0, 13.33, bar_height, DARK_BLUE)
    # accent stripe
    add_rect(slide, 0, bar_height - 0.07, 13.33, 0.07, ACCENT_BLUE)

    txBox = slide.shapes.add_textbox(
        Inches(0.35), Inches(0.12), Inches(12.5), Inches(bar_height - 0.2)
    )
    tf = txBox.text_frame
    tf.word_wrap = False
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = title_text
    run.font.name = "Calibri"
    run.font.size = Pt(28)
    run.font.bold = True
    run.font.color.rgb = WHITE
    return bar


def set_slide_background(slide, color):
    from pptx.oxml.ns import qn
    from lxml import etree
    background = slide.background
    fill = background.fill
    fill.solid()
    fill.fore_color.rgb = color


# ---------------------------------------------------------------------------
# Slide 1 — Title
# ---------------------------------------------------------------------------

def slide1(prs):
    slide = blank_slide(prs)
    set_slide_background(slide, DARK_BLUE)

    # decorative accent rectangle bottom
    add_rect(slide, 0, 6.8, 13.33, 0.7, MID_BLUE)
    add_rect(slide, 0, 7.3, 13.33, 0.2, ACCENT_BLUE)

    # horizontal rule
    add_rect(slide, 0.5, 3.55, 12.33, 0.06, ACCENT_BLUE)

    # Main title
    tb = slide.shapes.add_textbox(Inches(0.5), Inches(1.5), Inches(12.33), Inches(1.6))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    run = p.add_run()
    run.text = "Point de Stage"
    run.font.name = "Calibri"
    run.font.size = Pt(44)
    run.font.bold = True
    run.font.color.rgb = WHITE

    p2 = tf.add_paragraph()
    p2.alignment = PP_ALIGN.CENTER
    run2 = p2.add_run()
    run2.text = "Modélisation PK/PD T-DXd"
    run2.font.name = "Calibri"
    run2.font.size = Pt(38)
    run2.font.bold = True
    run2.font.color.rgb = ACCENT_BLUE

    # Subtitle
    tb2 = slide.shapes.add_textbox(Inches(0.5), Inches(3.7), Inches(12.33), Inches(1.2))
    tf2 = tb2.text_frame
    tf2.word_wrap = True
    p3 = tf2.paragraphs[0]
    p3.alignment = PP_ALIGN.CENTER
    run3 = p3.add_run()
    run3.text = "T-DXd (Trastuzumab Deruxtecan) — Toxicité Hématologique Rat"
    run3.font.name = "Calibri"
    run3.font.size = Pt(22)
    run3.font.bold = False
    run3.font.color.rgb = RGBColor(0xC5, 0xD5, 0xEE)

    # Date
    tb3 = slide.shapes.add_textbox(Inches(0.5), Inches(6.85), Inches(12.33), Inches(0.55))
    tf3 = tb3.text_frame
    p4 = tf3.paragraphs[0]
    p4.alignment = PP_ALIGN.CENTER
    run4 = p4.add_run()
    run4.text = "Avril 2026"
    run4.font.name = "Calibri"
    run4.font.size = Pt(16)
    run4.font.color.rgb = WHITE


# ---------------------------------------------------------------------------
# Slide 2 — Contexte & Objectif
# ---------------------------------------------------------------------------

def slide2(prs):
    slide = blank_slide(prs)
    set_slide_background(slide, LIGHT_GREY)
    add_title_bar(slide, "Contexte & Objectif")

    bullets = [
        "T-DXd : ADC (DAR=8), topoisomérase I inhibiteur (DXd)",
        "Indication : cancer du sein HER2+ (FDA approuvé)",
        "Problématique : toxicité hématologique dose-limitante (neutropénie, thrombopénie)",
        "Objectif : modèle PK/PD mécanistique RAT pour prédire les nadirs hématologiques",
        "Base : modèle de Fornari 2019 (carboplatin/rat) → étendu à T-DXd",
    ]

    tb = slide.shapes.add_textbox(Inches(0.55), Inches(1.3), Inches(12.2), Inches(5.8))
    tf = tb.text_frame
    tf.word_wrap = True

    first = True
    for b in bullets:
        if first:
            p = tf.paragraphs[0]
            first = False
        else:
            p = tf.add_paragraph()
        p.space_before = Pt(10)
        p.level = 0

        from pptx.oxml.ns import qn
        from lxml import etree
        pPr = p._p.get_or_add_pPr()
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', '▶')

        run = p.add_run()
        run.text = b
        run.font.name = "Calibri"
        run.font.size = Pt(19)
        run.font.color.rgb = DARK_GREY


# ---------------------------------------------------------------------------
# Slide 3 — Architecture du modèle
# ---------------------------------------------------------------------------

def slide3(prs):
    from lxml import etree
    from pptx.oxml.ns import qn

    slide = blank_slide(prs)
    set_slide_background(slide, WHITE)
    add_title_bar(slide, "Architecture du Modèle PK/PD (25 états)")

    # Left column box
    left_box = add_rect(slide, 0.4, 1.25, 5.9, 5.85, RGBColor(0xD9, 0xE2, 0xF3))
    left_box.line.color.rgb = MID_BLUE
    left_box.line.width = Pt(1)

    # Right column box
    right_box = add_rect(slide, 6.9, 1.25, 5.9, 5.85, RGBColor(0xFD, 0xF0, 0xD5))
    right_box.line.color.rgb = RGBColor(0xC9, 0x9A, 0x00)
    right_box.line.width = Pt(1)

    # Left header
    lh = add_rect(slide, 0.4, 1.25, 5.9, 0.5, MID_BLUE)
    lh.line.fill.background()
    add_textbox(slide, "Chaîne PK T-DXd (5 états)", 0.42, 1.27, 5.86, 0.45,
                font_size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Right header
    rh = add_rect(slide, 6.9, 1.25, 5.9, 0.5, RGBColor(0xBF, 0x8F, 0x00))
    rh.line.fill.background()
    add_textbox(slide, "Cascade PD Fornari (20 états)", 6.92, 1.27, 5.86, 0.45,
                font_size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Left bullets
    left_items = [
        "ADC sérum (2 compartiments : Yin 2020)",
        "DXd plasma (libération via Krel temps-dép.)",
        "DXd intracellulaire (moelle, Vasalou 2024)",
        "Damage ADN (γH2AX, modèle Emax)",
    ]
    tb_l = slide.shapes.add_textbox(Inches(0.65), Inches(1.85), Inches(5.4), Inches(5.0))
    tf_l = tb_l.text_frame
    tf_l.word_wrap = True
    for i, item in enumerate(left_items):
        p = tf_l.paragraphs[0] if i == 0 else tf_l.add_paragraph()
        p.space_before = Pt(10)
        pPr = p._p.get_or_add_pPr()
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', '•')
        run = p.add_run()
        run.text = item
        run.font.name = "Calibri"
        run.font.size = Pt(17)
        run.font.color.rgb = DARK_GREY

    # Right bullets
    right_items = [
        "MPP → CMP / MEP",
        "CMP → Neutrophiles, Monocytes",
        "MEP → Réticulocytes → GR, Plaquettes",
        "Connexion : Damage → Slope_MPP/CMP/MEP",
    ]
    tb_r = slide.shapes.add_textbox(Inches(7.15), Inches(1.85), Inches(5.4), Inches(5.0))
    tf_r = tb_r.text_frame
    tf_r.word_wrap = True
    for i, item in enumerate(right_items):
        p = tf_r.paragraphs[0] if i == 0 else tf_r.add_paragraph()
        p.space_before = Pt(10)
        pPr = p._p.get_or_add_pPr()
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', '•')
        run = p.add_run()
        run.text = item
        run.font.name = "Calibri"
        run.font.size = Pt(17)
        run.font.color.rgb = DARK_GREY

    # Arrow connecting the two boxes
    add_textbox(slide, "⟹", 6.25, 3.5, 0.6, 0.6,
                font_size=28, bold=True, color=DARK_BLUE, align=PP_ALIGN.CENTER)


# ---------------------------------------------------------------------------
# Slide 4 — Paramètres PK clés (table)
# ---------------------------------------------------------------------------

def slide4(prs):
    from pptx.util import Inches, Pt
    from pptx.oxml.ns import qn
    from lxml import etree

    slide = blank_slide(prs)
    set_slide_background(slide, LIGHT_GREY)
    add_title_bar(slide, "Paramètres PK — Allométrie Yin 2020 → Rat (250 g)")

    headers = ["Paramètre", "Humain", "Rat (allométrie)"]
    rows = [
        ["CL_ADC",        "0.421 L/j",    "2.56×10⁻⁴ L/h"],
        ["V1_ADC",        "2.77 L",        "9.89 mL"],
        ["DAR",           "8",             "8"],
        ["Krel (cycle 1)","0.0159 h⁻¹",   "0.0159 h⁻¹"],
        ["IC50_DXd",      "0.31 µM",       "0.31 µM"],
        ["k_dam = k_rep", "—",             "0.017 h⁻¹"],
    ]

    col_widths = [3.5, 3.5, 3.5]
    row_height = 0.62
    table_left = 1.4
    table_top  = 1.4
    n_rows = len(rows) + 1
    n_cols = 3

    tbl = slide.shapes.add_table(
        n_rows, n_cols,
        Inches(table_left), Inches(table_top),
        Inches(sum(col_widths)), Inches(row_height * n_rows)
    ).table

    # Column widths
    for ci, w in enumerate(col_widths):
        tbl.columns[ci].width = Inches(w)

    def set_cell(cell, text, bg_color, txt_color, bold=False, font_size=16, align=PP_ALIGN.CENTER):
        cell.fill.solid()
        cell.fill.fore_color.rgb = bg_color
        tf = cell.text_frame
        tf.word_wrap = False
        p = tf.paragraphs[0]
        p.alignment = align
        run = p.add_run()
        run.text = text
        run.font.name = "Calibri"
        run.font.size = Pt(font_size)
        run.font.bold = bold
        run.font.color.rgb = txt_color
        # cell margins
        cell._tc.get_or_add_tcPr()

    # Header row
    for ci, h in enumerate(headers):
        set_cell(tbl.cell(0, ci), h, TABLE_HDR, WHITE, bold=True, font_size=17)

    # Data rows
    for ri, row in enumerate(rows):
        bg = TABLE_ALT if ri % 2 == 0 else WHITE
        for ci, val in enumerate(row):
            align = PP_ALIGN.LEFT if ci == 0 else PP_ALIGN.CENTER
            set_cell(tbl.cell(ri + 1, ci), val, bg, DARK_GREY,
                     bold=(ci == 0), font_size=16, align=align)

    # Note below table
    add_textbox(slide, "Allométrie basée sur Yin et al. 2020 — poids rat = 250 g",
                1.4, table_top + row_height * n_rows + 0.15, 10.5, 0.45,
                font_size=13, italic=True, color=RGBColor(0x70, 0x70, 0x70))


# ---------------------------------------------------------------------------
# Slide 5 — Résultats PK 5 mg/kg
# ---------------------------------------------------------------------------

def slide5(prs):
    slide = blank_slide(prs)
    set_slide_background(slide, WHITE)
    add_title_bar(slide, "Résultats PK — T-DXd 5 mg/kg Q3W × 3 cycles (Rat)")

    img_path = os.path.join(IMG_DIR, "PK_5mgkg.png")
    if os.path.exists(img_path):
        slide.shapes.add_picture(
            img_path,
            Inches(0.4), Inches(1.25),
            Inches(12.53), Inches(5.55)
        )

    add_textbox(slide,
                "ADC sérum  |  DXd intracellulaire  |  Damage ADN — Krel temps-dépendant (Yin 2020)",
                0.4, 6.88, 12.53, 0.45,
                font_size=13, italic=True, color=RGBColor(0x50, 0x50, 0x50),
                align=PP_ALIGN.CENTER)


# ---------------------------------------------------------------------------
# Slide 6 — Résultats PD — comparaison des doses
# ---------------------------------------------------------------------------

def slide6(prs):
    slide = blank_slide(prs)
    set_slide_background(slide, WHITE)
    add_title_bar(slide, "Résultats PD — Comparaison 5 mg/kg vs 10 mg/kg")

    img5  = os.path.join(IMG_DIR, "PD_5mgkg.png")
    img10 = os.path.join(IMG_DIR, "PD_10mgkg.png")

    if os.path.exists(img5):
        slide.shapes.add_picture(img5,  Inches(0.25), Inches(1.25), Inches(6.35), Inches(5.35))
    if os.path.exists(img10):
        slide.shapes.add_picture(img10, Inches(6.73), Inches(1.25), Inches(6.35), Inches(5.35))

    # Labels above images
    add_textbox(slide, "5 mg/kg", 0.25, 1.25, 6.35, 0.0,
                font_size=14, bold=True, color=MID_BLUE, align=PP_ALIGN.CENTER)
    add_textbox(slide, "10 mg/kg", 6.73, 1.25, 6.35, 0.0,
                font_size=14, bold=True, color=MID_BLUE, align=PP_ALIGN.CENTER)

    add_textbox(slide,
                "Damage max : 0.069 (5 mg/kg) → 0.14 (10 mg/kg) — "
                "effet PD modéré attendu chez le rat (absence d'effet bystander HER2+)",
                0.25, 6.88, 12.83, 0.45,
                font_size=13, italic=True, color=RGBColor(0x50, 0x50, 0x50),
                align=PP_ALIGN.CENTER)


# ---------------------------------------------------------------------------
# Slide 7 — Translatabilité clinique
# ---------------------------------------------------------------------------

def slide7(prs):
    from lxml import etree
    from pptx.oxml.ns import qn

    slide = blank_slide(prs)
    set_slide_background(slide, LIGHT_GREY)
    add_title_bar(slide, "Translatabilité — Cadre Clinique")

    bullets = [
        ("Cavg clinique T-DXd (6.4 mg/kg Q3W) = 33.3 µg/mL", False),
        ("IC50 ADC ex vivo (CFU assay PFB-10) : 27.3 µg/mL (érythroïde), 28.1 µg/mL (CFU-GM)", False),
        ("Ratio Cavg/IC50 ≈ 1.2 → zone d'effet pharmacologique", True),
        ("Prochaine étape : driver PD = C_ADC1 avec IC50_ADC = 27.3 µg/mL (au lieu de DXd_ic)", False),
        ("Calibration sur données précliniques rat à venir", False),
    ]

    tb = slide.shapes.add_textbox(Inches(0.55), Inches(1.35), Inches(12.2), Inches(5.8))
    tf = tb.text_frame
    tf.word_wrap = True

    for i, (text, highlight) in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_before = Pt(12)
        pPr = p._p.get_or_add_pPr()
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', '▶')
        run = p.add_run()
        run.text = text
        run.font.name = "Calibri"
        run.font.size = Pt(18)
        run.font.bold = highlight
        run.font.color.rgb = DARK_BLUE if highlight else DARK_GREY


# ---------------------------------------------------------------------------
# Slide 8 — Structure du code & Prochaines étapes
# ---------------------------------------------------------------------------

def slide8(prs):
    from lxml import etree
    from pptx.oxml.ns import qn

    slide = blank_slide(prs)
    set_slide_background(slide, WHITE)
    add_title_bar(slide, "Structure du Code & Perspectives")

    # Left column box — Structure R
    left_box = add_rect(slide, 0.4, 1.25, 5.9, 5.85, RGBColor(0xD9, 0xE2, 0xF3))
    left_box.line.color.rgb = MID_BLUE
    left_box.line.width = Pt(1)

    lh = add_rect(slide, 0.4, 1.25, 5.9, 0.5, MID_BLUE)
    lh.line.fill.background()
    add_textbox(slide, "Structure R", 0.42, 1.27, 5.86, 0.45,
                font_size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Right column box — Prochaines étapes
    right_box = add_rect(slide, 6.9, 1.25, 5.9, 5.85, RGBColor(0xE2, 0xF0, 0xD9))
    right_box.line.color.rgb = RGBColor(0x37, 0x86, 0x10)
    right_box.line.width = Pt(1)

    rh = add_rect(slide, 6.9, 1.25, 5.9, 0.5, RGBColor(0x37, 0x86, 0x10))
    rh.line.fill.background()
    add_textbox(slide, "Prochaines étapes", 6.92, 1.27, 5.86, 0.45,
                font_size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Left content
    left_lines = [
        ("scripts/ : modèle Fornari (carboplatin)", 0),
        ("scripts_tdxd/ : extension T-DXd", 0),
        ("parameters_tdxd_rat.R (allométrie)", 1),
        ("pkpd_tdxd_rat.R (ODE fusionné 25 états)", 1),
        ("run_pkpd_tdxd_rat.R (scénarios)", 1),
    ]

    tb_l = slide.shapes.add_textbox(Inches(0.65), Inches(1.85), Inches(5.4), Inches(5.1))
    tf_l = tb_l.text_frame
    tf_l.word_wrap = True

    for i, (text, level) in enumerate(left_lines):
        p = tf_l.paragraphs[0] if i == 0 else tf_l.add_paragraph()
        p.space_before = Pt(8)
        p.level = level
        pPr = p._p.get_or_add_pPr()
        if level == 0:
            buChar = etree.SubElement(pPr, qn('a:buChar'))
            buChar.set('char', '•')
        else:
            buChar = etree.SubElement(pPr, qn('a:buChar'))
            buChar.set('char', '◦')
        run = p.add_run()
        run.text = ("    " if level == 1 else "") + text
        run.font.name = "Calibri"
        run.font.size = Pt(15 if level == 0 else 14)
        run.font.bold = (level == 0)
        run.font.color.rgb = DARK_GREY

    # Right content
    right_lines = [
        "Remplacer driver DXd_ic → C_ADC1 (IC50=27.3 µg/mL)",
        "Calibrer Slope_MPP/CMP/MEP sur données rat",
        "Validation sur données précliniques observées",
        "Extension : simulations sensibilité (dose, schéma)",
    ]

    tb_r = slide.shapes.add_textbox(Inches(7.15), Inches(1.85), Inches(5.4), Inches(5.1))
    tf_r = tb_r.text_frame
    tf_r.word_wrap = True

    for i, text in enumerate(right_lines):
        p = tf_r.paragraphs[0] if i == 0 else tf_r.add_paragraph()
        p.space_before = Pt(10)
        pPr = p._p.get_or_add_pPr()
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', '✓')
        run = p.add_run()
        run.text = text
        run.font.name = "Calibri"
        run.font.size = Pt(16)
        run.font.color.rgb = DARK_GREY


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    prs = new_prs()
    slide1(prs)
    slide3(prs)
    slide4(prs)
    slide5(prs)
    slide6(prs)
    slide7(prs)
    slide8(prs)
    prs.save(OUT_PATH)
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
