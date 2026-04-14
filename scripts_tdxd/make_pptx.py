#!/usr/bin/env python3
"""
make_pptx.py — Presentation PK/PD T-DXd : resultats et prochaines etapes
"""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import os

# ── Couleurs ─────────────────────────────────────────────
BLUE_DARK  = RGBColor(0x1a, 0x3a, 0x5c)   # titre
BLUE_MED   = RGBColor(0x21, 0x66, 0xac)   # accent
RED        = RGBColor(0xb2, 0x18, 0x2b)   # alerte
GREEN      = RGBColor(0x1a, 0x7a, 0x3c)   # OK
ORANGE     = RGBColor(0xe6, 0x7e, 0x00)   # warning
GREY_LIGHT = RGBColor(0xf0, 0xf4, 0xf8)   # fond slide
GREY_TEXT  = RGBColor(0x44, 0x44, 0x44)

def new_prs():
    prs = Presentation()
    prs.slide_width  = Inches(13.33)
    prs.slide_height = Inches(7.5)
    return prs

def blank_slide(prs):
    layout = prs.slide_layouts[6]   # totalement vide
    return prs.slides.add_slide(layout)

def bg(slide, color=GREY_LIGHT):
    fill = slide.background.fill
    fill.solid()
    fill.fore_color.rgb = color

def box(slide, l, t, w, h, text, font_size=14, bold=False,
        color=GREY_TEXT, bg_color=None, align=PP_ALIGN.LEFT,
        italic=False, wrap=True):
    txBox = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf = txBox.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    if bg_color:
        fill = txBox.fill
        fill.solid()
        fill.fore_color.rgb = bg_color
    return txBox

def rect(slide, l, t, w, h, fill_color, line_color=None):
    from pptx.util import Inches
    shape = slide.shapes.add_shape(
        1,  # MSO_SHAPE_TYPE.RECTANGLE
        Inches(l), Inches(t), Inches(w), Inches(h)
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    if line_color:
        shape.line.color.rgb = line_color
        shape.line.width = Pt(1)
    else:
        shape.line.fill.background()
    return shape

def hline(slide, l, t, w, color=BLUE_MED, width_pt=1.5):
    line = slide.shapes.add_shape(2, Inches(l), Inches(t), Inches(w), Pt(2))
    line.fill.background()
    line.line.color.rgb = color
    line.line.width = Pt(width_pt)

# ══════════════════════════════════════════════════════════════
# SLIDE 1 — TITRE
# ══════════════════════════════════════════════════════════════
def slide_titre(prs):
    sl = blank_slide(prs)
    bg(sl, RGBColor(0x1a, 0x3a, 0x5c))

    box(sl, 0.5, 0.8, 12.3, 1.5,
        "Modele PK/PD T-DXd — Resultats et Prochaines Etapes",
        font_size=32, bold=True, color=RGBColor(0xff,0xff,0xff),
        align=PP_ALIGN.CENTER)

    box(sl, 0.5, 2.5, 12.3, 0.6,
        "Trastuzumab deruxtecan (ENHERTU) — Toxicite hematologique",
        font_size=18, color=RGBColor(0xaa, 0xcc, 0xee),
        align=PP_ALIGN.CENTER)

    hline(sl, 2.0, 3.3, 9.3, color=RGBColor(0x74,0xad,0xd1), width_pt=1)

    items = [
        "Modele : PK 2-compartiments ADC + Fornari 2019 (hematopoiese 20 etats)",
        "Population : N=300 patients, IIV log-normale (Yin 2020)",
        "Validation : FDA BLA 761139 — DESTINY-Breast01 (n=184)",
        "Dose : 5.4 mg/kg Q3W x6 cycles",
    ]
    for i, item in enumerate(items):
        box(sl, 2.5, 3.6 + i*0.62, 9.0, 0.55, "  " + item,
            font_size=14, color=RGBColor(0xdd, 0xee, 0xff))

# ══════════════════════════════════════════════════════════════
# SLIDE 2 — SCHEMA DU MODELE
# ══════════════════════════════════════════════════════════════
def slide_schema(prs):
    sl = blank_slide(prs)
    bg(sl)

    box(sl, 0.4, 0.2, 12.5, 0.6, "Structure du modele PK/PD",
        font_size=22, bold=True, color=BLUE_DARK)
    hline(sl, 0.4, 0.85, 12.5)

    # Chaine PK
    blocks_pk = [
        ("ADC serum\n2 compartiments\n(Yin 2020)", 0.5, 1.1, 2.2, 1.4, RGBColor(0x21,0x66,0xac)),
        ("DXd plasma\n1 compartiment\nCmax=3.8 ng/mL", 3.2, 1.1, 2.2, 1.4, RGBColor(0x44,0x88,0xbb)),
        ("DXd intra-\ncellulaire\n(x35)", 5.9, 1.1, 2.0, 1.4, RGBColor(0x66,0x99,0xcc)),
        ("Dommages ADN\n(Damage)", 8.4, 1.1, 2.0, 1.4, RGBColor(0x99,0x55,0x88)),
    ]
    for txt, l, t, w, h, col in blocks_pk:
        rect(sl, l, t, w, h, col)
        box(sl, l+0.05, t+0.1, w-0.1, h-0.15, txt,
            font_size=11, bold=True, color=RGBColor(0xff,0xff,0xff),
            align=PP_ALIGN.CENTER)

    # fleches PK
    for x in [2.72, 5.42, 8.12]:
        box(sl, x, 1.6, 0.5, 0.4, "->", font_size=16, bold=True,
            color=GREY_TEXT, align=PP_ALIGN.CENTER)

    # Driver
    box(sl, 0.5, 2.75, 10.0, 0.4,
        "Driver toxicite : C_ADC1 [ug/mL]  (surrogate — justification FDA BLA 761139 p.95)",
        font_size=12, italic=True, color=ORANGE)

    # Fornari
    rect(sl, 0.5, 3.3, 12.0, 3.4, RGBColor(0xf5,0xf5,0xf5),
         line_color=BLUE_MED)
    box(sl, 0.7, 3.35, 5.0, 0.4, "Modele hematopoiese Fornari 2019 (20 etats)",
        font_size=13, bold=True, color=BLUE_DARK)

    cell_types = [
        ("MPP", 0.8, 3.9, 1.5, 0.7, RGBColor(0x6b,0xae,0xd6)),
        ("CMP", 2.6, 3.9, 1.5, 0.7, RGBColor(0x74,0xad,0xd1)),
        ("MEP", 4.4, 3.9, 1.5, 0.7, RGBColor(0x74,0xad,0xd1)),
        ("Neut", 2.4, 5.0, 1.5, 0.7, RGBColor(0xd6,0x60,0x4d)),
        ("RBC/Ret", 4.2, 5.0, 1.6, 0.7, RGBColor(0xf4,0xa5,0x82)),
        ("Plt", 6.0, 5.0, 1.5, 0.7, RGBColor(0xfd,0xae,0x61)),
    ]
    for txt, l, t, w, h, col in cell_types:
        rect(sl, l, t, w, h, col)
        box(sl, l, t+0.1, w, h-0.15, txt,
            font_size=13, bold=True, color=RGBColor(0xff,0xff,0xff),
            align=PP_ALIGN.CENTER)

    box(sl, 8.0, 3.7, 4.2, 2.8,
        "Kill terms :\n\nSlope_CMP x Damage\n-> CMP (neutrophiles)\n\nSlope_MEP x Damage\n-> MEP (erythrocytes)",
        font_size=12, color=GREY_TEXT)

# ══════════════════════════════════════════════════════════════
# SLIDE 3 — RESULTATS GRADES
# ══════════════════════════════════════════════════════════════
def slide_grades(prs):
    sl = blank_slide(prs)
    bg(sl)

    box(sl, 0.4, 0.2, 12.5, 0.6,
        "Resultats : Grades CTCAE vs FDA DESTINY-Breast01",
        font_size=22, bold=True, color=BLUE_DARK)
    hline(sl, 0.4, 0.85, 12.5)

    # Neutropenie
    rect(sl, 0.4, 1.0, 5.9, 5.8, RGBColor(0xeb,0xf3,0xfb), line_color=BLUE_MED)
    box(sl, 0.6, 1.05, 5.5, 0.5, "NEUTROPENIE",
        font_size=16, bold=True, color=BLUE_DARK)

    neut_data = [
        ("G0  (>=2.0)", "0%",    "71%",  False),
        ("G1  (1.5-2.0)", "14%", "7%",   False),
        ("G2  (1.0-1.5)", "70%", "7%",   False),
        ("G3  (0.5-1.0)", "16%", "13%",  False),
        ("G4  (<0.5)",    "0%",  "3%",   False),
        ("G3-4 TOTAL",   "15.7%","16-20%", True),
        ("Tout grade",   "100%", "29%",  True),
    ]
    headers = ["Grade", "Modele", "FDA", ""]
    for j, h in enumerate(headers[:3]):
        box(sl, 0.6+j*1.7, 1.6, 1.6, 0.35, h,
            font_size=11, bold=True, color=BLUE_DARK)
    hline(sl, 0.6, 1.95, 5.5, color=BLUE_MED, width_pt=0.8)

    for i, (grade, mod, fda, is_total) in enumerate(neut_data):
        y = 2.05 + i*0.47
        if is_total:
            rect(sl, 0.5, y-0.03, 5.8, 0.43,
                 RGBColor(0xd0,0xe8,0xf5), line_color=None)
        box(sl, 0.6, y, 1.65, 0.4, grade,
            font_size=11, bold=is_total, color=GREY_TEXT)
        # couleur modele
        mc = GREEN if (grade=="G3-4 TOTAL") else (RED if mod=="100%" else GREY_TEXT)
        box(sl, 2.3, y, 1.2, 0.4, mod,
            font_size=11, bold=is_total, color=mc, align=PP_ALIGN.CENTER)
        box(sl, 4.0, y, 1.2, 0.4, fda,
            font_size=11, bold=is_total, color=GREY_TEXT, align=PP_ALIGN.CENTER)
        # tick/cross
        if grade == "G3-4 TOTAL":
            box(sl, 5.3, y, 0.6, 0.4, "OK", font_size=11, bold=True, color=GREEN)
        elif grade == "Tout grade":
            box(sl, 5.3, y, 0.6, 0.4, "(!)", font_size=11, bold=True, color=RED)

    # Anemie
    rect(sl, 6.9, 1.0, 5.9, 5.8, RGBColor(0xfe, 0xf0, 0xeb), line_color=RED)
    box(sl, 7.1, 1.05, 5.5, 0.5, "ANEMIE",
        font_size=16, bold=True, color=RED)

    anem_data = [
        ("G0  (RBC >90%)", "0%",  "30%",  False),
        ("G1  (RBC 83-90%)", "0%","37%",  False),
        ("G2  (RBC 67-83%)", "94%","24%", False),
        ("G3  (RBC 54-67%)", "6%", "8%",  False),
        ("G4  (RBC <54%)", "0%",  "1%",   False),
        ("G3-4 TOTAL",     "6.0%","~9%",  True),
        ("Tout grade",     "100%","~70%", True),
    ]
    for j, h in enumerate(headers[:3]):
        box(sl, 7.1+j*1.7, 1.6, 1.6, 0.35, h,
            font_size=11, bold=True, color=GREY_TEXT)
    hline(sl, 7.1, 1.95, 5.5, color=RED, width_pt=0.8)

    for i, (grade, mod, fda, is_total) in enumerate(anem_data):
        y = 2.05 + i*0.47
        if is_total:
            rect(sl, 6.9, y-0.03, 5.8, 0.43,
                 RGBColor(0xf9,0xdd,0xd5), line_color=None)
        box(sl, 7.1, y, 1.65, 0.4, grade,
            font_size=11, bold=is_total, color=GREY_TEXT)
        mc = GREEN if grade=="G3-4 TOTAL" else (ORANGE if mod=="100%" else GREY_TEXT)
        box(sl, 8.8, y, 1.2, 0.4, mod,
            font_size=11, bold=is_total, color=mc, align=PP_ALIGN.CENTER)
        box(sl, 10.5, y, 1.2, 0.4, fda,
            font_size=11, bold=is_total, color=GREY_TEXT, align=PP_ALIGN.CENTER)
        if grade == "G3-4 TOTAL":
            box(sl, 11.8, y, 0.7, 0.4, "~OK", font_size=11, bold=True, color=GREEN)
        elif grade == "Tout grade":
            box(sl, 11.8, y, 0.7, 0.4, "(!)", font_size=11, bold=True, color=ORANGE)

# ══════════════════════════════════════════════════════════════
# SLIDE 4 — CE QUI MARCHE / CE QUI NE MARCHE PAS
# ══════════════════════════════════════════════════════════════
def slide_bilan(prs):
    sl = blank_slide(prs)
    bg(sl)

    box(sl, 0.4, 0.2, 12.5, 0.6, "Bilan : Forces et Limitations du Modele",
        font_size=22, bold=True, color=BLUE_DARK)
    hline(sl, 0.4, 0.85, 12.5)

    # OK
    rect(sl, 0.4, 1.0, 5.9, 5.8, RGBColor(0xeb,0xf7,0xed), line_color=GREEN)
    box(sl, 0.6, 1.05, 5.5, 0.5, "Ce qui fonctionne",
        font_size=16, bold=True, color=GREEN)
    ok_items = [
        "PK ADC validee : Cmax 121 ug/mL (FDA 122) +0.8%",
        "G3-4 neutropenie : 15.7%  (FDA 16-20%)",
        "G3-4 anemie : 6.0%  (FDA ~9%) proche",
        "Chaine mecanistique complete :\n  ADC -> DXd -> Damage -> Kill",
        "IIV log-normale calibree\n  (Yin 2020 / Fornari 2019)",
        "Scaling rat -> humain coherent",
    ]
    for i, item in enumerate(ok_items):
        box(sl, 0.8, 1.65 + i*0.78, 5.3, 0.7,
            "  " + item, font_size=12, color=GREY_TEXT)
        box(sl, 0.55, 1.65 + i*0.78, 0.3, 0.4,
            "v", font_size=14, bold=True, color=GREEN)

    # Problemes
    rect(sl, 6.9, 1.0, 5.9, 5.8, RGBColor(0xfe, 0xf0, 0xeb), line_color=RED)
    box(sl, 7.1, 1.05, 5.5, 0.5, "Limitations identifiees",
        font_size=16, bold=True, color=RED)
    pb_items = [
        ("100% tout grade neutropenie\n  (FDA ~29%) — surestimation",   RED),
        ("Driver toxicite : C_ADC1\n  (FDA utilise DXd Cavg — p.95)", ORANGE),
        ("DXd plasma T½ modele = 1h\n  vs T½ apparent FDA = 5.8j",     ORANGE),
        ("G0 neutropenie : 0%\n  (FDA ~71%) — tous les patients\n  ont une chute de neutrophiles", RED),
        ("Anémie tout grade : 100%\n  (FDA ~70%) — proche mais\n  distribution G0/G1/G2 decalee", ORANGE),
    ]
    for i, (item, col) in enumerate(pb_items):
        box(sl, 7.3, 1.65 + i*0.92, 5.2, 0.85,
            "  " + item, font_size=12, color=GREY_TEXT)
        box(sl, 7.05, 1.65 + i*0.92, 0.3, 0.4,
            "x", font_size=14, bold=True, color=col)

# ══════════════════════════════════════════════════════════════
# SLIDE 5 — CAUSE RACINE : TOUT GRADE
# ══════════════════════════════════════════════════════════════
def slide_cause(prs):
    sl = blank_slide(prs)
    bg(sl)

    box(sl, 0.4, 0.2, 12.5, 0.6,
        "Cause Racine : Pourquoi 100% Tout Grade ?",
        font_size=22, bold=True, color=BLUE_DARK)
    hline(sl, 0.4, 0.85, 12.5)

    box(sl, 0.5, 1.0, 12.0, 0.5,
        "Probleme fondamental : T½ ADC ≈ 23 jours  >  Intervalle Q3W = 21 jours",
        font_size=15, bold=True, color=RED, align=PP_ALIGN.CENTER)

    # Schema PK
    rect(sl, 0.5, 1.65, 12.0, 2.2, RGBColor(0xf0,0xf4,0xf8), line_color=BLUE_MED)
    box(sl, 0.7, 1.7, 8.0, 0.4, "Profil C_ADC1 sur 6 cycles Q3W :", font_size=12, bold=True, color=BLUE_DARK)
    box(sl, 0.7, 2.1, 11.5, 1.5,
        "  Cycle 1       Cycle 2       Cycle 3       Cycle 4       Cycle 5       Cycle 6\n"
        "   Cmax=127     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "   --------    Trough ≈ 60 ug/mL (E_drug ≈ 68% IC50)    -> jamais nulle",
        font_size=11, color=GREY_TEXT)

    box(sl, 0.5, 4.0, 12.0, 0.45,
        "=> E_drug ne redescend jamais a 0 entre cycles => Damage residuel => tous les patients affectes",
        font_size=13, italic=True, color=RED, align=PP_ALIGN.CENTER)

    rect(sl, 0.5, 4.55, 5.7, 2.5, RGBColor(0xeb,0xf7,0xed), line_color=GREEN)
    box(sl, 0.7, 4.6, 5.3, 0.45, "Solutions explorees :", font_size=13, bold=True, color=GREEN)
    sols = [
        "n_hill=2 (Hill sigmoide) : attenue\n  E_drug sous IC50",
        "Revenir a n_hill=1 + accepter la\n  limitation sur tout grade",
        "Reviser DXd PK (T½ court) pour\n  driver DXd_ic coherent",
    ]
    for i, s in enumerate(sols):
        box(sl, 0.9, 5.1+i*0.6, 5.1, 0.55, s, font_size=11, color=GREY_TEXT)
        box(sl, 0.65, 5.1+i*0.6, 0.3, 0.4, "-", font_size=14, bold=True, color=GREEN)

    rect(sl, 6.7, 4.55, 5.9, 2.5, RGBColor(0xfe, 0xf9, 0xec), line_color=ORANGE)
    box(sl, 6.9, 4.6, 5.5, 0.45, "Choix actuel :", font_size=13, bold=True, color=ORANGE)
    box(sl, 6.9, 5.05, 5.5, 1.9,
        "n_hill=1 + use_ADC_driver=TRUE\n\n"
        "Priorite : reproduire G3-4 correctement\n"
        "Tout grade non reproduit = limitation\n"
        "documentee (acceptable pour calibration\n"
        "dose-finding)",
        font_size=11, color=GREY_TEXT)

# ══════════════════════════════════════════════════════════════
# SLIDE 6 — PROCHAINES ETAPES
# ══════════════════════════════════════════════════════════════
def slide_next(prs):
    sl = blank_slide(prs)
    bg(sl)

    box(sl, 0.4, 0.2, 12.5, 0.6, "Prochaines Etapes",
        font_size=22, bold=True, color=BLUE_DARK)
    hline(sl, 0.4, 0.85, 12.5)

    steps = [
        (
            "1. Corriger le DXd PK (Priorite haute)",
            RED,
            [
                "T½ DXd modele = 1.06h   vs   FDA apparent = 5.8j",
                "Ajouter un 2e compartiment DXd (lié/libre) ou ralentir k_effD",
                "Objectif : DXd_ic Cavg suffisant pour driver toxicite",
                "-> Permettra d'utiliser C_DXd_ic comme driver (aligne FDA p.95)",
            ]
        ),
        (
            "2. Recalibrer Slope_CMP apres correction PK (Priorite haute)",
            RED,
            [
                "Quand DXd_ic Cavg sera realiste -> re-scanner Slope_CMP",
                "Cible inchangee : G3-4 neutropenie = 16-20%",
            ]
        ),
        (
            "3. Reproduire tout grade (Priorite moyenne)",
            ORANGE,
            [
                "Tester n_hill=2 avec DXd_ic driver (une fois PK corrigee)",
                "n_hill=2 reduit E_drug sub-IC50 -> attenue nadir entre cycles",
                "Cible : tout grade ~29% (FDA DESTINY-Breast01)",
            ]
        ),
        (
            "4. Valider PD rat (Priorite basse)",
            BLUE_MED,
            [
                "Aucune donnee experimentale rat utilisee pour PD (nadir Neut/Plt)",
                "Ajouter comparaison avec donnees precliniques T-DXd si disponibles",
            ]
        ),
    ]

    y = 1.05
    for title, col, bullets in steps:
        rect(sl, 0.4, y, 12.5, 0.45, col)
        box(sl, 0.55, y+0.03, 12.0, 0.4, title,
            font_size=13, bold=True, color=RGBColor(0xff,0xff,0xff))
        y += 0.48
        for b in bullets:
            box(sl, 0.8, y, 12.0, 0.38, "   " + b, font_size=11, color=GREY_TEXT)
            box(sl, 0.55, y, 0.3, 0.35, "->", font_size=11, color=col)
            y += 0.37
        y += 0.1

# ══════════════════════════════════════════════════════════════
# SLIDE 7 — RESUME EXECUTIF
# ══════════════════════════════════════════════════════════════
def slide_resume(prs):
    sl = blank_slide(prs)
    bg(sl, RGBColor(0x1a, 0x3a, 0x5c))

    box(sl, 0.5, 0.5, 12.3, 0.65, "Resume Executif",
        font_size=26, bold=True, color=RGBColor(0xff,0xff,0xff),
        align=PP_ALIGN.CENTER)
    hline(sl, 1.5, 1.2, 10.3, color=RGBColor(0x74,0xad,0xd1))

    cards = [
        ("PK ADC", "Cmax = 121 ug/mL\nFDA = 122 ug/mL\n+0.8%", GREEN),
        ("G3-4 Neut", "Modele : 15.7%\nFDA : 16-20%\nOK", GREEN),
        ("G3-4 Anemie", "Modele : 6.0%\nFDA : ~9%\nProche", GREEN),
        ("Tout grade", "Modele : 100%\nFDA : 29%\nA corriger", RED),
        ("Driver PD", "C_ADC1\nFDA : DXd Cavg\nLimitation", ORANGE),
    ]

    for i, (title, body, col) in enumerate(cards):
        x = 0.5 + i * 2.55
        rect(sl, x, 1.5, 2.3, 2.8, RGBColor(0x25,0x4a,0x70), line_color=col)
        rect(sl, x, 1.5, 2.3, 0.5, col)
        box(sl, x, 1.52, 2.3, 0.46, title,
            font_size=13, bold=True, color=RGBColor(0xff,0xff,0xff),
            align=PP_ALIGN.CENTER)
        box(sl, x+0.1, 2.1, 2.1, 2.1, body,
            font_size=12, color=RGBColor(0xdd,0xee,0xff),
            align=PP_ALIGN.CENTER)

    box(sl, 0.5, 4.55, 12.3, 0.5,
        "Modele fonctionnel pour G3-4. Priorite : corriger T½ DXd pour driver mecanistique aligne FDA.",
        font_size=14, bold=True, color=RGBColor(0xff,0xdd,0x77),
        align=PP_ALIGN.CENTER)

    box(sl, 0.5, 5.2, 12.3, 1.8,
        "Branche git : claude/analyze-script-tdxd-Ox5X4\n"
        "Fichiers cles : parameters_tdxd_human.R | pkpd_tdxd_rat.R | run_pkpd_tdxd_human_population.R\n"
        "Resultats : scripts_tdxd/results_PKPD_human/",
        font_size=11, color=RGBColor(0x99,0xbb,0xdd),
        align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    os.chdir("/home/user/Hema/scripts_tdxd")
    os.makedirs("results_PKPD_human", exist_ok=True)

    prs = new_prs()
    slide_titre(prs)
    slide_schema(prs)
    slide_grades(prs)
    slide_bilan(prs)
    slide_cause(prs)
    slide_next(prs)
    slide_resume(prs)

    out = "results_PKPD_human/TDXd_PKPD_Bilan.pptx"
    prs.save(out)
    print(f"-> {out}  ({len(prs.slides)} slides)")
