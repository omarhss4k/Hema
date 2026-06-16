"""
Génération de la présentation Soutenance_Hassani_Omar_final.pptx
Pipeline PK/PD hématotoxicité ADCs — Omar Hassani
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.oxml.ns import qn
from lxml import etree
import os

# ── Palette ──────────────────────────────────────────────────────────────────
NAVY      = RGBColor(0x1B, 0x4F, 0x72)
TEAL      = RGBColor(0x2E, 0x86, 0xAB)
WHITE     = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT_BG  = RGBColor(0xF7, 0xFB, 0xFD)
GRAY      = RGBColor(0x5D, 0x6D, 0x7E)
YELLOW    = RGBColor(0xF3, 0x9C, 0x12)
RED_LIGHT = RGBColor(0xCB, 0x4C, 0x35)
GREEN     = RGBColor(0x1E, 0x8B, 0x4C)

prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]

IMG_DIR = "/home/user/Hema"

# ── Helpers ───────────────────────────────────────────────────────────────────

def set_bg(slide, color):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_rect(slide, l, t, w, h, fill_color=None, line_color=None, line_width_pt=0):
    shape = slide.shapes.add_shape(1, Inches(l), Inches(t), Inches(w), Inches(h))
    if fill_color:
        shape.fill.solid()
        shape.fill.fore_color.rgb = fill_color
    else:
        shape.fill.background()
    if line_color:
        shape.line.color.rgb = line_color
        shape.line.width = Pt(line_width_pt)
    else:
        shape.line.fill.background()
    return shape


def add_tb(slide, text, l, t, w, h,
           font="Calibri", size=18, bold=False, italic=False,
           color=NAVY, align=PP_ALIGN.LEFT,
           bg=None, wrap=True):
    txBox = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf = txBox.text_frame
    tf.word_wrap = wrap
    tf.margin_left  = Inches(0.06)
    tf.margin_top   = Inches(0.04)
    tf.margin_right = Inches(0.04)
    tf.margin_bottom= Inches(0.04)
    for i, line in enumerate(text.split('\n')):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.alignment = align
        run = p.add_run()
        run.text = line
        run.font.name  = font
        run.font.size  = Pt(size)
        run.font.bold  = bold
        run.font.italic= italic
        run.font.color.rgb = color
    if bg:
        txBox.fill.solid()
        txBox.fill.fore_color.rgb = bg
    return txBox


def add_title_bar(slide, text, dark=False):
    add_rect(slide, 0, 0, 13.33, 1.3, fill_color=NAVY)
    tb = add_tb(slide, text, 0.35, 0.12, 12.6, 1.05,
                font="Cambria", size=28, bold=True, color=WHITE, align=PP_ALIGN.LEFT)
    add_rect(slide, 0, 1.3, 13.33, 0.07, fill_color=TEAL)
    return tb


def add_slide_transition(slide):
    xml = '<p:transition xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" spd="med" dur="500"><p:fade/></p:transition>'
    slide._element.append(etree.fromstring(xml))


_id = [200]

def add_animations(slide, shapes):
    if not shapes:
        return
    base = _id[0]
    _id[0] += len(shapes) * 6 + 10

    parts = []
    for idx, sh in enumerate(shapes):
        cid  = base + idx * 6
        sid  = base + idx * 6 + 1
        aid  = base + idx * 6 + 2
        delay = idx * 350
        ntype = "clickEffect" if idx == 0 else "afterEffect"
        parts.append(f'''<p:par xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"
             xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
  <p:cTn id="{cid}" presetID="10" presetClass="entr" presetSubtype="0"
         fill="hold" grpId="0" nodeType="{ntype}">
    <p:stCondLst><p:cond delay="{delay}"/></p:stCondLst>
    <p:childTnLst>
      <p:set>
        <p:cBhvr><p:cTn id="{sid}" dur="1" fill="hold"/>
          <p:tgtEl><p:spTgt spid="{sh.shape_id}"/></p:tgtEl>
          <p:attrNameLst><p:attrName>style.visibility</p:attrName></p:attrNameLst>
        </p:cBhvr>
        <p:to><p:strVal val="visible"/></p:to>
      </p:set>
      <p:animEffect transition="in" filter="fade">
        <p:cBhvr><p:cTn id="{aid}" dur="500"/>
          <p:tgtEl><p:spTgt spid="{sh.shape_id}"/></p:tgtEl>
        </p:cBhvr>
      </p:animEffect>
    </p:childTnLst>
  </p:cTn>
</p:par>''')

    seq_id  = base + len(shapes) * 6
    root_id = base + len(shapes) * 6 + 1
    child_xml = "\n".join(parts)

    timing_xml = f'''<p:timing xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"
    xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
  <p:tnLst><p:par>
    <p:cTn id="{root_id}" dur="indefinite" restart="whenNotActive" nodeType="tmRoot">
      <p:childTnLst>
        <p:seq concurrent="1" nextAc="seek">
          <p:cTn id="{seq_id}" dur="indefinite" nodeType="mainSeq">
            <p:childTnLst>{child_xml}</p:childTnLst>
          </p:cTn>
          <p:prevCondLst><p:cond evt="onPrevClick" delay="0"><p:tn/></p:cond></p:prevCondLst>
          <p:nextCondLst><p:cond evt="onNextClick" delay="0"><p:tn/></p:cond></p:nextCondLst>
        </p:seq>
      </p:childTnLst>
    </p:cTn>
  </p:par></p:tnLst>
</p:timing>'''
    slide._element.append(etree.fromstring(timing_xml))


def add_table(slide, headers, rows, l, t, w, h,
              hdr_bg=NAVY, hdr_fg=WHITE, row_bgs=None, col_widths=None):
    ncols = len(headers)
    nrows = len(rows) + 1
    tbl_shape = slide.shapes.add_table(nrows, ncols,
                                        Inches(l), Inches(t), Inches(w), Inches(h))
    tbl = tbl_shape.table
    if col_widths:
        for i, cw in enumerate(col_widths):
            tbl.columns[i].width = Inches(cw)

    def style_cell(cell, text, bg, fg, bold=False, size=12, align=PP_ALIGN.CENTER):
        cell.fill.solid()
        cell.fill.fore_color.rgb = bg
        tf = cell.text_frame
        tf.word_wrap = True
        p = tf.paragraphs[0]
        p.alignment = align
        run = p.add_run()
        run.text = text
        run.font.name = "Calibri"
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = fg

    for ci, h_txt in enumerate(headers):
        style_cell(tbl.cell(0, ci), h_txt, hdr_bg, hdr_fg, bold=True, size=13)

    for ri, row in enumerate(rows):
        bg = row_bgs[ri] if row_bgs else (LIGHT_BG if ri % 2 == 0 else WHITE)
        for ci, cell_txt in enumerate(row):
            fg = NAVY
            if cell_txt in ("✓",): fg = GREEN
            elif cell_txt in ("⚠",): fg = YELLOW
            style_cell(tbl.cell(ri+1, ci), str(cell_txt), bg, fg, size=12,
                       align=PP_ALIGN.LEFT if ci == 0 else PP_ALIGN.CENTER)
    return tbl_shape


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 1 — Titre
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, NAVY)
add_slide_transition(s)
add_rect(s, 0, 0, 13.33, 0.18, fill_color=TEAL)
add_rect(s, 0, 7.32, 13.33, 0.18, fill_color=TEAL)
add_rect(s, 0.35, 1.1, 0.14, 4.8, fill_color=TEAL)

a = []
a.append(add_tb(s, "Modélisation de l'hématotoxicité\ninduite par les ADCs",
    0.65, 1.0, 12.0, 2.5,
    font="Cambria", size=40, bold=True, color=WHITE, align=PP_ALIGN.LEFT))
a.append(add_tb(s, "Développement d'un pipeline PK/PD de l'animal au patient",
    0.65, 3.6, 12.0, 0.8, size=22, color=TEAL, align=PP_ALIGN.LEFT))
add_rect(s, 0.65, 4.55, 11.5, 0.05, fill_color=TEAL)
a.append(add_tb(s, "Omar Hassani  |  Hervé Perdry  |  Aurélie Petain  |  Pierre Fabre  |  2025–2026",
    0.65, 4.75, 12.0, 0.6, size=16, color=RGBColor(0xAE, 0xD6, 0xF1)))
a.append(add_tb(s, "M2 Science des Données de Santé — Université Paris-Saclay",
    0.65, 5.42, 12.0, 0.5, size=14, italic=True, color=RGBColor(0x7F, 0xB3, 0xD3)))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 2 — Les ADCs
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Les anticorps-drogue conjugués (ADCs)")
a = []

for i, (title, sub, body, bg) in enumerate([
    ("Anticorps monoclonal", "Guidage tumoral", "Reconnaissance spécifique\nde l'antigène tumoral\n(ex. HER2, TROP2, FGFR2)", NAVY),
    ("Linker chimique", "Liant clivable", "Stable en circulation\nCliváble en intra-cellulaire\n(pH acide / protéases)", TEAL),
    ("Payload cytotoxique", "Agent cytotoxique", "Topoisomérase I inhibiteur\n(DXd) ou agent alkylant\n(MMAE, DM1…)", RGBColor(0x1A, 0x5C, 0x38)),
]):
    l = 0.35 + i * 4.32
    a.append(add_rect(s, l, 1.55, 4.1, 2.2, fill_color=bg))
    a.append(add_tb(s, title, l+0.1, 1.6, 3.9, 0.55, size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
    a.append(add_tb(s, sub, l+0.1, 2.15, 3.9, 0.38, size=13, italic=True, color=TEAL if bg!=TEAL else WHITE, align=PP_ALIGN.CENTER))
    add_rect(s, l+0.1, 2.55, 3.9, 0.04, fill_color=WHITE)
    a.append(add_tb(s, body, l+0.1, 2.62, 3.9, 1.0, size=13, color=WHITE))
    if i < 2:
        a.append(add_tb(s, "+", l+4.14, 2.3, 0.3, 0.7, size=28, bold=True, color=TEAL, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 3.95, 12.6, 0.65, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Mécanisme :  Liaison antigène → Internalisation → Libération payload → Dommages ADN → Apoptose",
    0.45, 3.98, 12.4, 0.58, size=14, bold=True, color=NAVY, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 4.75, 12.6, 0.7, fill_color=NAVY))
a.append(add_tb(s, "Exemple clé : T-DXd (Enhertu®)  —  ORR 61%  —  DESTINY-Breast01  —  FDA approval 2019",
    0.45, 4.78, 12.4, 0.62, size=15, color=WHITE, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 5.6, 12.6, 0.55, fill_color=RGBColor(0xF0, 0xF8, 0xFF)))
a.append(add_tb(s, "→ Environ 100+ ADCs en développement clinique — nouvelle classe thérapeutique majeure en oncologie",
    0.45, 5.62, 12.4, 0.5, size=13, italic=True, color=GRAY))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 3 — Hématotoxicité
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Hématotoxicité : mécanisme et enjeux cliniques")
a = []

a.append(add_rect(s, 0.35, 1.52, 12.6, 0.65, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Payload libre → Moelle osseuse → Progéniteurs hématopoïétiques → Cytopénies",
    0.45, 1.55, 12.4, 0.58, size=15, bold=True, color=NAVY, align=PP_ALIGN.CENTER))

for i, (icon, title, body, col) in enumerate([
    ("🔴", "Neutropénie", "Infections opportunistes\nRisque vital si Grade 4\n16% G≥3 — DESTINY-Breast01", RED_LIGHT),
    ("🟠", "Anémie", "Fatigue sévère\nRéduction qualité de vie\nTransfusion si Grade 4", YELLOW),
    ("🔵", "Thrombocytopénie", "Risque hémorragique\nSaignements spontanés\nPlasmaphérèse si G4", RGBColor(0x21, 0x61, 0x8A)),
]):
    l = 0.35 + i * 4.32
    a.append(add_rect(s, l, 2.35, 4.1, 2.7, fill_color=col))
    a.append(add_tb(s, icon + "  " + title, l+0.1, 2.4, 3.9, 0.62, size=17, bold=True, color=WHITE))
    add_rect(s, l+0.1, 3.05, 3.9, 0.04, fill_color=WHITE)
    a.append(add_tb(s, body, l+0.1, 3.12, 3.9, 1.8, size=13, color=WHITE))

a.append(add_rect(s, 0.35, 5.22, 6.0, 0.72, fill_color=RGBColor(0xFD, 0xF2, 0xDB)))
a.append(add_tb(s, "⚠  Chiffre clé : 16% neutropénie G≥3\ndans DESTINY-Breast01 (FDA 2019)",
    0.45, 5.25, 5.8, 0.65, size=13, bold=True, color=RGBColor(0x7D, 0x60, 0x08)))

a.append(add_rect(s, 6.6, 5.22, 6.1, 0.72, fill_color=NAVY))
a.append(add_tb(s, "Double enjeu :\nSécurité patient FIH + Optimisation posologique",
    6.7, 5.25, 5.9, 0.65, size=13, bold=True, color=WHITE))

a.append(add_tb(s, "Les cytopénies surviennent 7–21 jours post-dose — gestion clinique critique — dose-limiting toxicity dans 30% ADCs",
    0.35, 6.12, 12.6, 0.5, size=12, italic=True, color=GRAY))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 4 — Pourquoi modéliser
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Pourquoi un modèle mécanistique ?")
a = []

for i, (title, detail, col) in enumerate([
    ("Données précliniques rares", "n=8 NHP, seulement 2 animaux par groupe de dose → statistiques classiques inapplicables", NAVY),
    ("Modèles empiriques insuffisants", "Régressions dose-réponse sans mécanisme : pas d'extrapolation inter-espèces possible", TEAL),
    ("Hystérèse PK/PD", "Le nadir cellulaire survient plusieurs jours après le pic plasmatique → seul un modèle mécanistique peut le reproduire", RED_LIGHT),
    ("Machine learning inadapté", "Trop peu de données + architecture boîte noire + pas d'extrapolation causale entre espèces", GRAY),
]):
    t = 1.55 + i * 1.22
    a.append(add_rect(s, 0.35, t, 0.1, 1.1, fill_color=col))
    a.append(add_tb(s, "▶  " + title, 0.55, t, 12.0, 0.45, size=16, bold=True, color=col))
    a.append(add_tb(s, detail, 0.55, t+0.46, 12.3, 0.68, size=14, color=NAVY))

a.append(add_rect(s, 0.35, 6.45, 12.6, 0.75, fill_color=NAVY))
a.append(add_tb(s, "✦  Seul un modèle PK/PD semi-mécanistique peut reproduire l'hystérèse et permettre l'extrapolation inter-espèces causalement interprétable",
    0.45, 6.47, 12.4, 0.7, size=14, bold=True, color=WHITE))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 5 — Objectifs
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Objectifs du stage — 4 étapes progressives")
a = []

for i, (num, title, sub, bg) in enumerate([
    ("1", "Rat + carboplatine", "Validation du cadre\nFornari (littérature)", NAVY),
    ("2", "Humain + T-DXd", "Simulation N=300\nValidation FDA BLA 761139", TEAL),
    ("3", "NHP + composé interne", "Données précliniques\nconfidentielles Pierre Fabre", RGBColor(0x1A, 0x5C, 0x38)),
    ("4", "Transposition clinique", "Prédiction dose FIH\nPK allométrique + PD", RGBColor(0x6C, 0x35, 0x83)),
]):
    l = 0.35 + i * 3.22
    a.append(add_rect(s, l, 1.55, 3.1, 3.0, fill_color=bg))
    a.append(add_tb(s, num, l+1.15, 1.6, 0.8, 0.7,
                    font="Cambria", size=32, bold=True, color=bg, bg=WHITE, align=PP_ALIGN.CENTER))
    a.append(add_tb(s, title, l+0.1, 2.42, 2.9, 0.65, size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
    a.append(add_tb(s, sub, l+0.1, 3.12, 2.9, 0.9, size=13, color=WHITE, align=PP_ALIGN.CENTER))
    if i < 3:
        a.append(add_tb(s, "→", l+3.12, 2.75, 0.18, 0.7, size=22, bold=True, color=TEAL, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 4.75, 12.6, 0.06, fill_color=TEAL))
a.append(add_tb(s, "Progression : du plus simple (rat/carboplatine connu) vers le plus complexe (composé interne/FIH)",
    0.35, 4.9, 12.6, 0.48, size=13, italic=True, color=GRAY, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 5.52, 12.6, 1.05, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Fil conducteur : Un seul cadre PK/PD unifié (Friberg/Fornari) — paramétrisé par espèce et composé\nReproductibilité totale : 28 scripts R modulaires, set.seed(42), versionnés sur Git",
    0.45, 5.55, 12.4, 1.0, size=14, color=NAVY))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 6 — Hématopoïèse
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Cascade hématopoïétique — cible du modèle PD")
a = []

nodes = [
    (0.35, 3.2, 2.0, 0.95, "Cellules Souches\n(HSC)", NAVY),
    (2.8,  3.2, 2.1, 0.95, "MPP\n(Multipotent Progenitors)", RGBColor(0x6C, 0x35, 0x83)),
    (5.6,  1.7, 2.0, 0.95, "CMP (Myéloïde)", TEAL),
    (5.6,  4.1, 2.0, 0.95, "MEP (Érythroïde-Mgk)", TEAL),
    (8.4,  0.9, 2.3, 0.85, "Neutrophiles / Monocytes", RGBColor(0x1A, 0x5C, 0x38)),
    (8.4,  2.35, 2.3, 0.85, "Éosinophiles / Basophiles", GRAY),
    (8.4,  3.8, 2.3, 0.85, "Réticulocytes → GR", RED_LIGHT),
    (8.4,  5.2, 2.3, 0.85, "Mégacaryocytes → PLT", RGBColor(0xD4, 0x6A, 0x00)),
]
for l, t, w, h, lbl, bg in nodes:
    a.append(add_rect(s, l, t, w, h, fill_color=bg))
    a.append(add_tb(s, lbl, l+0.05, t+0.12, w-0.1, h-0.15, size=12, bold=True, color=WHITE, align=PP_ALIGN.CENTER))

# Simple arrow lines
arrow_pairs = [
    (2.35, 3.67, 2.8, 3.67),
    (4.9, 3.67, 5.6, 2.17),
    (4.9, 3.67, 5.6, 4.57),
    (7.6, 2.17, 8.4, 1.32),
    (7.6, 2.17, 8.4, 2.77),
    (7.6, 4.57, 8.4, 4.22),
    (7.6, 4.57, 8.4, 5.62),
]
for x1, y1, x2, y2 in arrow_pairs:
    a.append(add_rect(s, x1, y1-0.03, x2-x1, 0.06, fill_color=TEAL))

# Drug target label
a.append(add_rect(s, 2.6, 2.3, 2.5, 0.6, fill_color=RED_LIGHT))
a.append(add_tb(s, "⚡ Cible du médicament\n(MPP — Multipotent Progenitors)", 2.65, 2.32, 2.4, 0.55,
                size=11, bold=True, color=WHITE, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 6.45, 12.6, 0.72, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Rétrocontrôle homéostatique : si Neut↓ → G-CSF↑ → moelle accélère production   |   Durée de vie : Neut ~6–8h | Plt ~9j (NHP) | GR ~120j",
    0.45, 6.48, 12.4, 0.65, size=13, color=NAVY))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 7 — Architecture modèle
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Architecture du modèle PK/PD semi-mécanistique")
a = []

for i, (title, eq, body, bg) in enumerate([
    ("BLOC 1 — Pharmacocinétique",
     "C(t) = A·e^(−αt) + B·e^(−βt)",
     "Modèle 2 compartiments\n4 paramètres : CL, V1, Q, V2\nUnivers des concentrations plasmatiques", NAVY),
    ("BLOC 2 — Dommage cellulaire",
     "dDmg/dt = k_dam×C(µM) − k_rep×Dmg",
     "Accumulation du dommage\nProportionnel à la concentration libre\nRéparation linéaire k_rep", TEAL),
    ("BLOC 3 — Pharmacodynamique",
     "Effet = 1 − Slope × Damage",
     "Cascade hématopoïétique\nNeutrophiles, GR, Plaquettes\nRetour à l'équilibre homéostatique", RGBColor(0x1A, 0x5C, 0x38)),
]):
    l = 0.35 + i * 4.32
    a.append(add_rect(s, l, 1.52, 4.1, 3.9, fill_color=bg))
    a.append(add_tb(s, title, l+0.1, 1.57, 3.9, 0.65, size=15, bold=True, color=WHITE))
    add_rect(s, l+0.1, 2.28, 3.9, 0.05, fill_color=TEAL if bg!=TEAL else WHITE)
    a.append(add_tb(s, eq, l+0.1, 2.38, 3.9, 0.65, size=13, bold=True, italic=True, color=TEAL if bg!=TEAL else WHITE))
    add_rect(s, l+0.1, 3.08, 3.9, 0.04, fill_color=WHITE)
    a.append(add_tb(s, body, l+0.1, 3.17, 3.9, 2.1, size=13, color=WHITE))
    if i < 2:
        a.append(add_tb(s, "→", l+4.14, 3.1, 0.25, 0.75, size=26, bold=True, color=TEAL, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 5.6, 12.6, 0.75, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Chaîne causale :   Dose  →  PK C(t)  →  Dommage cellulaire  →  PD  →  Cellules circulantes observées",
    0.45, 5.62, 12.4, 0.65, size=15, bold=True, color=NAVY, align=PP_ALIGN.CENTER))

a.append(add_tb(s, "Chaque paramètre a une signification biologique — aucun paramètre purement statistique — interprétabilité totale",
    0.35, 6.5, 12.6, 0.5, size=13, italic=True, color=GRAY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 8 — Méthodes PK
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Méthodes d'estimation PK")
a = []

a.append(add_rect(s, 0.35, 1.52, 5.9, 4.0, fill_color=NAVY))
a.append(add_tb(s, "Étape 1 — NCA\n(Analyse Non Compartimentale)", 0.45, 1.55, 5.7, 0.85, size=16, bold=True, color=WHITE))
a.append(add_tb(s,
    "Estimation directe des paramètres :\n\n"
    "  CL = Dose / AUC∞\n\n"
    "  V1 = Dose / C₀\n\n"
    "→ Fournit les valeurs initiales\n"
    "   pour l'optimisation numérique\n\n"
    "→ Pas d'hypothèse de structure",
    0.45, 2.48, 5.7, 2.9, size=14, color=WHITE))

a.append(add_tb(s, "→", 6.3, 3.2, 0.4, 0.85, size=28, bold=True, color=TEAL, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 6.85, 1.52, 6.1, 4.0, fill_color=TEAL))
a.append(add_tb(s, "Étape 2 — Nelder-Mead\n(Optimisation numérique)", 6.95, 1.55, 5.9, 0.85, size=16, bold=True, color=WHITE))
a.append(add_tb(s,
    "Critère : SSR = Σ[log(Cobs) − log(Cpred)]²\n\n"
    "Optimisation sur log(θ) :\n"
    "  • Garantit positivité paramètres\n"
    "  • Symétrie des résidus relatifs\n\n"
    "Multi-start n=10 points de départ\nCV inter-start < 5% sur CL et V1\n\n"
    "→ Convergence globale robuste",
    6.95, 2.48, 5.9, 2.9, size=13, color=WHITE))

a.append(add_rect(s, 0.35, 5.7, 12.6, 0.65, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "✓  Critère de satisfaction : résidus médians < 20%  |  Log-transformation → normalisation résidus relatifs",
    0.45, 5.72, 12.4, 0.6, size=15, bold=True, color=NAVY, align=PP_ALIGN.CENTER))

a.append(add_tb(s, "Appliqué indépendamment à chaque animal/patient — estimation individuelle (non NLME à ce stade)",
    0.35, 6.5, 12.6, 0.5, size=13, italic=True, color=GRAY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 9 — Monte-Carlo & CTCAE
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Simulation Monte-Carlo & Grading CTCAE v5")
a = []

a.append(add_rect(s, 0.35, 1.52, 5.5, 3.3, fill_color=NAVY))
a.append(add_tb(s, "Simulation Monte-Carlo", 0.45, 1.55, 5.3, 0.6, size=17, bold=True, color=WHITE))
a.append(add_tb(s,
    "N = 300 patients virtuels\n\n"
    "θᵢ = θpop × exp(ηᵢ)\n  Distribution log-normale\n\n"
    "Variabilité PK :\n  CV 25–30% (FDA BLA 761139)\n\n"
    "Variabilité PD :\n  CV 30% sur Slopes\n\n"
    "set.seed(42) → reproductibilité exacte",
    0.45, 2.22, 5.3, 2.5, size=13, color=WHITE))

a.append(add_table(s,
    ["Toxicité", "Grade 1", "Grade 2", "Grade 3", "Grade 4"],
    [
        ["Neutrophiles (×10⁹/L)", "1.5–1.9", "1.0–1.5", "0.5–1.0 ⚠", "<0.5 🔴"],
        ["Plaquettes (×10⁹/L)",   "75–LNI",  "50–75",   "25–50 ⚠",   "<25 🔴"],
        ["Hémoglobine (g/dL)",    "10–LNI",  "8–10",    "<8 ⚠",       "Transf. 🔴"],
        ["Leucocytes (×10⁹/L)",   "3.0–3.9", "2.0–3.0", "1.0–2.0 ⚠", "<1.0 🔴"],
    ],
    l=6.05, t=1.52, w=6.95, h=3.3, col_widths=[2.3, 1.0, 1.0, 1.2, 1.45]))

a.append(add_rect(s, 0.35, 5.0, 12.6, 0.58, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "G0 Normal  |  G1 Léger  |  G2 Modéré  |  G3 Sévère (hospit.)  |  G4 Critique (pronostic vital)",
    0.45, 5.02, 12.4, 0.52, size=13, bold=True, color=NAVY, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 5.7, 12.6, 0.72, fill_color=NAVY))
a.append(add_tb(s, "Proportion G≥3 = métrique principale de validation — comparaison directe aux données cliniques FDA BLA 761139",
    0.45, 5.72, 12.4, 0.65, size=14, bold=True, color=WHITE, align=PP_ALIGN.CENTER))

a.append(add_tb(s, "CTCAE v5 NCI — standard réglementaire FDA/EMA pour essais oncologiques",
    0.35, 6.57, 12.6, 0.45, size=12, italic=True, color=GRAY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 10 — Résultats Étape 1
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Étape 1 : Validation rat/carboplatine (Fornari 2011)")
a = []

a.append(add_rect(s, 0.35, 1.52, 12.6, 0.6, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Rat Sprague-Dawley  |  Carboplatine 40 mg/kg  |  8 cycles  |  Données Fornari et al. 2011 (JPKPD)",
    0.45, 1.54, 12.4, 0.55, size=14, bold=True, color=NAVY, align=PP_ALIGN.CENTER))

a.append(add_table(s,
    ["Lignée cellulaire", "Résidu médian (%)", "Interprétation clinique", "Statut"],
    [
        ["Neutrophiles (Neut)",  "2.67%",  "Excellent ajustement — profil nadir reproduit", "✓"],
        ["CMP (Myéloïde)",       "2.37%",  "Excellent ajustement — dynamique correcte", "✓"],
        ["GR (Érythrocytes)",    "2.79%",  "Excellent ajustement — stabilité bien reproduite", "✓"],
        ["MPP (Multipotent)",    "26.96%", "Digitalisation figure originale — non cliniquement pertinent", "⚠"],
        ["MEP (Érythroïde-Mgk)","50.69%", "Digitalisation figure originale — valeurs proches de zéro", "⚠"],
        ["Plaquettes (Plt)",     "25.67%", "Digitalisation figure originale — cinétique globale reproduite", "⚠"],
    ],
    l=0.35, t=2.28, w=12.6, h=3.55,
    col_widths=[2.3, 1.9, 6.0, 1.1],
    row_bgs=[LIGHT_BG, LIGHT_BG, LIGHT_BG,
             RGBColor(0xFD, 0xF2, 0xDB), RGBColor(0xFD, 0xF2, 0xDB), RGBColor(0xFD, 0xF2, 0xDB)]))

a.append(add_rect(s, 0.35, 6.0, 12.6, 0.75, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s,
    "⚠  Résidus élevés MPP/MEP/Plt attribuables à la digitalisation des figures originales (WebPlotDigitizer)\n"
    "    Données Neut/CMP/GR directement disponibles — celles-ci sont cliniquement pertinentes et validées",
    0.45, 6.02, 12.4, 0.7, size=13, color=NAVY))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 11 — Rat → Humain
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Transposition rat → humain (carboplatine)")
a = []

a.append(add_rect(s, 0.35, 1.52, 12.6, 0.55, fill_color=NAVY))
a.append(add_tb(s, "Adaptation PK — Allométrie inter-espèces", 0.45, 1.54, 12.4, 0.5, size=16, bold=True, color=WHITE))

a.append(add_rect(s, 0.35, 2.12, 6.0, 1.35, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "CL_humain = CL_rat × (70 / BW_rat)^0.75\nV proportionnel au poids corporel\nLoi allométrique inter-espèces (West 2002)",
    0.45, 2.15, 5.8, 1.25, size=14, color=NAVY))

a.append(add_rect(s, 6.55, 2.12, 6.1, 1.35, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Adaptation PD — Recalibration des Slopes\nSlopes recalibrées sur données\ncellulaires humaines disponibles (litterature)",
    6.65, 2.15, 5.9, 1.25, size=14, color=NAVY))

a.append(add_rect(s, 0.35, 3.62, 12.6, 0.55, fill_color=TEAL))
a.append(add_tb(s, "Résultats de la transposition", 0.45, 3.64, 12.4, 0.5, size=16, bold=True, color=WHITE))

for i, (label, detail, col) in enumerate([
    ("✓  Neutrophiles", "Concordance excellente sur profil cinétique complet (nadir, récupération, baseline)", GREEN),
    ("✓  Plaquettes — profil global", "Accord global satisfaisant (Δ = 6.4%) — profil cinétique bien reproduit", GREEN),
    ("⚠  Plaquettes — 2e cycle", "Accumulation 2e cycle non reproduite — limite propre à l'architecture Fornari, non à l'implémentation Python", YELLOW),
]):
    t = 4.28 + i * 0.72
    a.append(add_tb(s, label, 0.5, t, 4.5, 0.58, size=15, bold=True, color=col))
    a.append(add_tb(s, detail, 5.1, t, 7.8, 0.58, size=14, color=NAVY))

a.append(add_rect(s, 0.35, 6.52, 12.6, 0.65, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Validation satisfaisante — ouvre la voie à l'application sur T-DXd (N=300 simulation Monte-Carlo)",
    0.45, 6.54, 12.4, 0.6, size=14, bold=True, color=NAVY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 12 — Résultat central T-DXd (SLIDE LE PLUS IMPORTANT)
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, NAVY)
add_slide_transition(s)
add_rect(s, 0, 0, 13.33, 0.18, fill_color=TEAL)
add_rect(s, 0, 7.32, 13.33, 0.18, fill_color=TEAL)
a = []

a.append(add_tb(s, "⭐  Résultat central : T-DXd / DESTINY-Breast01",
    0.35, 0.22, 12.6, 0.88, font="Cambria", size=30, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
a.append(add_tb(s, "T-DXd 5.4 mg/kg Q3W × 6 cycles  |  N=300 patients virtuels (Monte-Carlo)  vs  FDA DESTINY-Breast01 n=184",
    0.35, 1.12, 12.6, 0.5, size=14, color=TEAL, align=PP_ALIGN.CENTER))

comps = [
    ("Neutropénie\ntout grade", "29.3%", "29%"),
    ("Neutropénie\nG3–4", "14%", "16%"),
    ("Anémie\ntout grade", "67.7%", "70%"),
    ("Anémie\nG3–4", "9.6%", "9%"),
]
for i, (label, mod, fda) in enumerate(comps):
    l = 0.35 + i * 3.24
    a.append(add_rect(s, l, 1.75, 3.1, 3.75, fill_color=RGBColor(0x21, 0x61, 0x8A)))
    a.append(add_tb(s, label, l+0.1, 1.8, 2.9, 0.7, size=14, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
    add_rect(s, l+0.15, 2.55, 2.8, 0.05, fill_color=TEAL)
    a.append(add_tb(s, "Modèle", l+0.1, 2.65, 2.9, 0.38, size=12, color=TEAL, align=PP_ALIGN.CENTER))
    a.append(add_tb(s, mod, l+0.1, 3.05, 2.9, 0.85, size=34, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
    a.append(add_tb(s, "FDA", l+0.1, 4.0, 2.9, 0.38, size=12, color=RGBColor(0xAE, 0xD6, 0xF1), align=PP_ALIGN.CENTER))
    a.append(add_tb(s, fda, l+0.1, 4.38, 2.9, 0.72, size=26, bold=True, color=RGBColor(0xAE, 0xD6, 0xF1), align=PP_ALIGN.CENTER))
    a.append(add_tb(s, "✓", l+1.0, 5.1, 1.1, 0.55, size=30, bold=True, color=TEAL, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 5.8, 12.6, 0.9, fill_color=TEAL))
a.append(add_tb(s,
    "Concordance remarquable — toutes proportions dans IC95% Wilson\n"
    "Le pipeline PK/PD reproduit fidèlement les données du dossier BLA FDA 761139",
    0.45, 5.82, 12.4, 0.85, size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 13 — NHP Design + PK
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "NHP cynomolgus : design expérimental & PK")
a = []

a.append(add_rect(s, 0.35, 1.52, 7.9, 1.75, fill_color=NAVY))
a.append(add_tb(s, "Design expérimental", 0.45, 1.55, 7.7, 0.52, size=16, bold=True, color=WHITE))
a.append(add_tb(s,
    "8 singes cynomolgus (Macaca fascicularis) — 4 groupes de 2\n"
    "Doses D1–D4 : 3 / 13 / 26 / 39 mg/kg\n"
    "Composé FGFR2 interne — données confidentielles Pierre Fabre",
    0.45, 2.14, 7.7, 1.08, size=13, color=WHITE))

a.append(add_rect(s, 0.35, 3.4, 7.9, 0.65, fill_color=RED_LIGHT))
a.append(add_tb(s, "⚠  D4 (39 mg/kg) : sacrifiés prématurément → toxicité sévère → PD non analysable",
    0.45, 3.42, 7.7, 0.6, size=13, bold=True, color=WHITE))

a.append(add_rect(s, 0.35, 4.18, 7.9, 2.4, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Paramètres PK moyens (D1–D3)", 0.45, 4.2, 7.7, 0.52, size=15, bold=True, color=NAVY))
for i, (param, val, note) in enumerate([
    ("CL",   "= 1.83 mL/h/kg", "CV inter-animal : 17%"),
    ("V1",   "= 54.2 mL/kg",   "CV inter-animal : 12%"),
    ("t½β",  "= 54.7 h",       "Demi-vie d'élimination terminale"),
    ("Résidus médians", "< 12%", "Tous animaux D1–D3 ✓"),
]):
    t = 4.82 + i * 0.46
    a.append(add_tb(s, param, 0.45, t, 1.8, 0.42, size=14, bold=True, color=TEAL))
    a.append(add_tb(s, val, 2.3, t, 2.5, 0.42, size=14, color=NAVY))
    a.append(add_tb(s, note, 4.9, t, 3.2, 0.42, size=12, italic=True, color=GRAY))

# Image placeholder
img_path = os.path.join(IMG_DIR, "cynomolgus.jpg")
if os.path.exists(img_path):
    s.shapes.add_picture(img_path, Inches(8.45), Inches(1.52), Inches(4.55), Inches(3.25))
else:
    a.append(add_rect(s, 8.45, 1.52, 4.55, 3.25, fill_color=RGBColor(0xD5, 0xE8, 0xF4),
                      line_color=TEAL, line_width_pt=1.5))
    a.append(add_tb(s, "Macaca fascicularis\n(singe cynomolgus)\n\nModèle préclinique standard\npour les ADCs",
                    8.55, 2.5, 4.35, 1.5, size=14, italic=True, color=GRAY, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 6.75, 12.6, 0.5, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "✓  Résidus médians < 12% pour tous les animaux D1–D3 — ajustement PK robuste et reproductible",
    0.45, 6.77, 12.4, 0.45, size=14, bold=True, color=GREEN))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 14 — NHP PD + Hystérèse
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "NHP : dynamique PD & hystérèse temporelle")
a = []

a.append(add_rect(s, 0.35, 1.52, 12.6, 1.12, fill_color=NAVY))
a.append(add_tb(s,
    "Hystérèse PK/PD : le nadir des neutrophiles précède celui des plaquettes de 1–2 semaines\n"
    "Durée de vie : Neutrophiles ~6–8 heures (NHP)  vs  Plaquettes ~9 jours (NHP)\n"
    "→ Signature temporelle distincte par lignée — mécanisme bien capturé par le modèle",
    0.45, 1.55, 12.4, 1.05, size=14, color=WHITE))

a.append(add_rect(s, 0.35, 2.78, 5.9, 2.85, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Le modèle reproduit :", 0.45, 2.81, 5.7, 0.5, size=15, bold=True, color=NAVY))
for i, txt in enumerate([
    "✓  Chute Neut/Plt cohérente D2–D3",
    "✓  Décalage temporel Neut → Plt (1–2 sem.)",
    "✓  Stabilité relative des GR (durée de vie longue)",
    "✓  Dose-réponse croissante D1 → D3",
]):
    a.append(add_tb(s, txt, 0.45, 3.42 + i * 0.52, 5.7, 0.48, size=13, color=GREEN))

a.append(add_rect(s, 6.55, 2.78, 6.4, 2.85, fill_color=RGBColor(0xFD, 0xF2, 0xDB)))
a.append(add_tb(s, "Limites :", 6.65, 2.81, 6.2, 0.5, size=15, bold=True, color=RGBColor(0x7D, 0x60, 0x08)))
for i, txt in enumerate([
    "⚠  n=2/groupe → variabilité non caractérisée",
    "⚠  Calibration PD visuelle sans IC formels",
    "⚠  Rebonds 400–600% non capturés",
    "⚠  D4 exclu → gamme incomplète",
]):
    a.append(add_tb(s, txt, 6.65, 3.42 + i * 0.52, 6.2, 0.48, size=13, color=RGBColor(0x7D, 0x60, 0x08)))

a.append(add_rect(s, 0.35, 5.78, 12.6, 0.55, fill_color=TEAL))
a.append(add_tb(s, "Calibration PD préliminaire — base suffisante pour la transposition clinique FIH",
    0.45, 5.8, 12.4, 0.5, size=14, bold=True, color=WHITE, align=PP_ALIGN.CENTER))

a.append(add_tb(s, "L'hystérèse PK/PD est la signature distinctive de l'hématotoxicité des ADCs — sans alternative empirique",
    0.35, 6.5, 12.6, 0.52, size=13, italic=True, color=GRAY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 15 — Transposition clinique
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Transposition clinique — Prédiction dose FIH")
a = []

for i, (title, eq, body, bg) in enumerate([
    ("PK — Allométrie", "CL_humain = CL_NHP × (70/BW_NHP)^0.75",
     "Validée sur T-DXd :\nerreur < facteur 2 sur CL\n\nApproche réglementairement acceptée\n(FDA MIDD guidance)", NAVY),
    ("PD — Conservatrice", "Slopes NHP → humain directement",
     "Hypothèse conservatrice :\nhumain plus sensible que NHP\n\nMarge de sécurité supplémentaire\npour les patients", TEAL),
    ("Simulation — FIH", "Dose FIH = argmin_d P(G3+) < 10%",
     "Profil toxicité modéré\nComparable à T-DXd\n\nSélection data-driven\nnon empirique", RGBColor(0x1A, 0x5C, 0x38)),
]):
    l = 0.35 + i * 4.32
    a.append(add_rect(s, l, 1.52, 4.1, 4.3, fill_color=bg))
    a.append(add_tb(s, title, l+0.1, 1.57, 3.9, 0.62, size=16, bold=True, color=WHITE))
    add_rect(s, l+0.1, 2.24, 3.9, 0.05, fill_color=TEAL if bg!=TEAL else WHITE)
    a.append(add_tb(s, eq, l+0.1, 2.33, 3.9, 0.62, size=12, bold=True, italic=True, color=TEAL if bg!=TEAL else WHITE))
    add_rect(s, l+0.1, 3.0, 3.9, 0.04, fill_color=WHITE)
    a.append(add_tb(s, body, l+0.1, 3.08, 3.9, 2.6, size=13, color=WHITE))
    if i < 2:
        a.append(add_tb(s, "→", l+4.14, 3.3, 0.25, 0.75, size=26, bold=True, color=TEAL, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 0.35, 6.0, 12.6, 0.65, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Marge de sécurité : Dose tolérée NHP (sans toxicité sévère) → Dose FIH humain par allométrie",
    0.45, 6.02, 12.4, 0.6, size=15, bold=True, color=NAVY, align=PP_ALIGN.CENTER))

a.append(add_tb(s, "FIH = First-In-Human  |  La simulation guide la dose initiale sans exposer les patients à des doses toxiques inutiles",
    0.35, 6.78, 12.6, 0.48, size=12, italic=True, color=GRAY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 16 — Forces
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Forces du pipeline")
a = []

forces = [
    ("① Généricité", "2 composés, 3 espèces,\n1 cadre unifié\n\nRat  →  Humain  →  NHP\nExtensible à tout ADC", NAVY),
    ("② Interprétabilité biologique", "Chaque paramètre =\nsignification mécanistique\n\nPas de boîte noire\nKnowledge-based", TEAL),
    ("③ Simulation in silico", "Doses non testées\nsans expérimentation\n\nN=300 en quelques secondes\nAnalyses de sensibilité", RGBColor(0x1A, 0x5C, 0x38)),
    ("④ Reproductibilité totale", "28 scripts R modulaires\nset.seed(42)\nVersionné sur Git\n\nAudit complet possible", RGBColor(0x6C, 0x35, 0x83)),
]
for i, (title, body, bg) in enumerate(forces):
    l = 0.35 + (i % 2) * 6.54
    t = 1.52 + (i // 2) * 2.62
    a.append(add_rect(s, l, t, 6.2, 2.42, fill_color=bg))
    a.append(add_tb(s, title, l+0.15, t+0.1, 5.9, 0.62, size=19, bold=True, color=WHITE))
    add_rect(s, l+0.15, t+0.78, 5.9, 0.05, fill_color=RGBColor(0xAE, 0xD6, 0xF1))
    a.append(add_tb(s, body, l+0.15, t+0.88, 5.9, 1.42, size=13, color=WHITE))

a.append(add_tb(s, "Validation en 3 contextes indépendants : rat/carboplatine (litérature), humain/T-DXd (FDA), NHP/FGFR2 (interne)",
    0.35, 6.9, 12.6, 0.48, size=13, italic=True, color=GRAY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 17 — Limites
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Limites & points d'attention")
a = []

limits = [
    ("n=2 par groupe NHP",
     "Variabilité inter-individuelle mal caractérisée → impossible de distinguer vraie variabilité biologique du bruit",
     "Élargir cohorte NHP ou données historiques"),
    ("Calibration PD visuelle",
     "Absence d'IC formels sur les Slopes → incertitude sur les paramètres non quantifiée rigoureusement",
     "NLME (nlmixr2) pour IC95% rigoureux"),
    ("Transposition PD inter-espèces",
     "Hypothèse de conservation des Slopes NHP→humain non démontrée → données CFU nécessaires pour validation",
     "Expériences CFU in vitro sur progéniteurs humains"),
    ("Pas de modèle d'efficacité",
     "Pipeline centré toxicité uniquement → analyse bénéfice-risque incomplète pour décision dose FIH finale",
     "Modélisation PK/PD tumorale complémentaire"),
]
for i, (title, detail, solution) in enumerate(limits):
    t = 1.55 + i * 1.27
    a.append(add_rect(s, 0.35, t, 0.32, 1.15, fill_color=RED_LIGHT))
    a.append(add_tb(s, str(i+1), 0.37, t+0.32, 0.28, 0.52, size=18, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
    a.append(add_tb(s, "⚠  " + title, 0.8, t, 11.6, 0.45, size=15, bold=True, color=RGBColor(0x7D, 0x60, 0x08)))
    a.append(add_tb(s, detail, 0.8, t+0.46, 6.8, 0.68, size=13, color=NAVY))
    a.append(add_rect(s, 7.8, t+0.36, 5.1, 0.72, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
    a.append(add_tb(s, "→ " + solution, 7.9, t+0.38, 4.9, 0.68, size=12, italic=True, color=TEAL))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 18 — Perspectives NLME
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Perspectives : vers le NLME (nlmixr2)")
a = []

a.append(add_rect(s, 0.35, 1.52, 5.9, 3.65, fill_color=RGBColor(0xD6, 0xEA, 0xF8)))
a.append(add_tb(s, "Approche actuelle — Limites", 0.45, 1.55, 5.7, 0.55, size=16, bold=True, color=GRAY))
a.append(add_tb(s,
    "Estimation individuelle Nelder-Mead\n→ Moyenne des paramètres individuels\n\n"
    "❌  Pas d'IC formels sur Slopes PD\n"
    "❌  Variabilité inter-individuelle ignorée\n"
    "❌  Données BLQ non intégrées\n"
    "❌  Corrélations entre paramètres perdues",
    0.45, 2.18, 5.7, 2.85, size=13, color=GRAY))

a.append(add_tb(s, "→", 6.38, 3.05, 0.45, 0.9, size=28, bold=True, color=TEAL, align=PP_ALIGN.CENTER))

a.append(add_rect(s, 7.1, 1.52, 5.9, 3.65, fill_color=NAVY))
a.append(add_tb(s, "NLME (nlmixr2) — Apports", 7.2, 1.55, 5.7, 0.55, size=16, bold=True, color=WHITE))
a.append(add_tb(s,
    "θpop + ω² + ε estimés simultanément\n\n"
    "✓  IC95% rigoureux sur les Slopes PD\n"
    "✓  Meilleure caractérisation variabilité\n"
    "✓  Intégration données BLQ (M3 method)\n"
    "✓  Matrice Ω : corrélations paramètres",
    7.2, 2.18, 5.7, 2.85, size=13, color=WHITE))

a.append(add_rect(s, 0.35, 5.32, 12.6, 0.85, fill_color=TEAL))
a.append(add_tb(s,
    "Approche bayésienne : priors rat/T-DXd → régulariser estimation NLME NHP (n=2/groupe)\n"
    "nlmixr2 est l'implémentation R de référence — intégrable directement dans le pipeline existant",
    0.45, 5.35, 12.4, 0.8, size=13, color=WHITE))

a.append(add_tb(s, "NLME = Non-Linear Mixed Effects  |  nlmixr2 : package R CRAN, intégration ODE, algorithme SAEM/FOCE",
    0.35, 6.35, 12.6, 0.48, size=12, italic=True, color=GRAY, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 19 — Impact opérationnel
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, LIGHT_BG)
add_slide_transition(s)
add_title_bar(s, "Impact opérationnel")
a = []

for i, (title, body, bg) in enumerate([
    ("① Dose FIH",
     "Simulation des grades G3-4\nà différentes doses\n\n→ Sélection dose initiale :\n   < 10% G3-4 prédit\n\nDécision data-driven,\nnon empirique", NAVY),
    ("② Monitoring clinique",
     "Prédiction du timing du nadir\npar lignée cellulaire\n\n→ Fréquence hémogrammes\n   adaptée à la cinétique\n\nSurveillance ciblée\npas systématique", TEAL),
    ("③ Comm. MIDD réglementaire",
     "Model-Informed Drug Development\n\n→ Valorisé FDA/EMA dans\n   dossiers IND/CTA\n\nAccélère dialogue\nréglementaire", RGBColor(0x1A, 0x5C, 0x38)),
]):
    l = 0.35 + i * 4.32
    a.append(add_rect(s, l, 1.52, 4.1, 4.65, fill_color=bg))
    a.append(add_tb(s, title, l+0.1, 1.57, 3.9, 0.65, size=17, bold=True, color=WHITE))
    add_rect(s, l+0.1, 2.28, 3.9, 0.05, fill_color=RGBColor(0xAE, 0xD6, 0xF1))
    a.append(add_tb(s, body, l+0.1, 2.38, 3.9, 3.65, size=13, color=WHITE))

a.append(add_rect(s, 0.35, 6.35, 12.6, 0.62, fill_color=RGBColor(0xFD, 0xF2, 0xDB)))
a.append(add_tb(s,
    "Note honnête : pipeline fonctionnel en scripts R (28 scripts modulaires) — pas encore déployé comme outil clé-en-main",
    0.45, 6.37, 12.4, 0.55, size=13, italic=True, color=RGBColor(0x7D, 0x60, 0x08)))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 20 — Conclusion (fond sombre)
# ═══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
set_bg(s, NAVY)
add_slide_transition(s)
add_rect(s, 0, 0, 13.33, 0.18, fill_color=TEAL)
add_rect(s, 0, 7.32, 13.33, 0.18, fill_color=TEAL)
a = []

a.append(add_tb(s, "Conclusions", 0.35, 0.22, 12.6, 0.88,
                font="Cambria", size=36, bold=True, color=WHITE, align=PP_ALIGN.CENTER))
add_rect(s, 2.5, 1.1, 8.33, 0.07, fill_color=TEAL)

msgs = [
    ("Pipeline PK/PD semi-mécanistique",
     "Complet, reproductible, multi-espèces & multi-composés — cadre Friberg/Fornari unifié"),
    ("Validation clinique T-DXd",
     "29.3% vs 29% neutropénie — 67.7% vs 70% anémie — concordance remarquable dossier FDA BLA 761139"),
    ("Application composé interne NHP",
     "PK robuste (résidus <12%) — calibration PD préliminaire — base solide pour la dose FIH"),
    ("Complémentarité SDS × Pharmacologie quantitative",
     "Optimisation posologique, simulation de cohortes, diagnostics statistiques rigoureux"),
]
for i, (title, detail) in enumerate(msgs):
    t = 1.28 + i * 1.3
    a.append(add_rect(s, 0.35, t, 0.14, 1.1, fill_color=TEAL))
    a.append(add_tb(s, title, 0.62, t, 12.0, 0.5, size=16, bold=True, color=TEAL))
    a.append(add_tb(s, detail, 0.62, t+0.52, 12.0, 0.7, size=14, color=RGBColor(0xAE, 0xD6, 0xF1)))

add_rect(s, 0.35, 6.48, 12.6, 0.07, fill_color=TEAL)
a.append(add_tb(s, "Merci de votre attention — je suis à votre disposition pour les questions",
    0.35, 6.62, 12.6, 0.58,
    font="Cambria", size=18, italic=True, color=WHITE, align=PP_ALIGN.CENTER))
add_animations(s, a)

# ═══════════════════════════════════════════════════════════════════════════════
# Save & verify
# ═══════════════════════════════════════════════════════════════════════════════
output_path = "/home/user/Hema/Soutenance_Hassani_Omar_final.pptx"
prs.save(output_path)
print(f"Saved: {output_path}")

from pptx import Presentation as Prs2
chk = Prs2(output_path)
import os
print(f"Slides: {len(chk.slides)} / 20")
print(f"Size:   {os.path.getsize(output_path):,} bytes")
print(f"Dims:   {chk.slide_width.inches:.2f}\" x {chk.slide_height.inches:.2f}\"")
for i, sl in enumerate(chk.slides):
    print(f"  Slide {i+1:02d}: {len(sl.shapes)} shapes")
print("DONE — presentation ready!")
