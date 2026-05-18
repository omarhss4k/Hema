"""
Slide didactique — Modèle de Simeoni (2004)
Inhibition de la croissance tumorale (TGI) — PK/PD préclinique
Même charte graphique que make_pptx.py (projet carboplatine)
"""
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

# ── Palette (identique à make_pptx.py) ───────────────────────────────────────
BLEU_FONCE = RGBColor(0x1A, 0x3A, 0x5C)
BLEU_CLAIR = RGBColor(0x2E, 0x75, 0xB6)
ORANGE     = RGBColor(0xC5, 0x50, 0x1A)
BLANC      = RGBColor(0xFF, 0xFF, 0xFF)
GRIS_CLAIR = RGBColor(0xF2, 0xF2, 0xF2)
VERT       = RGBColor(0x37, 0x86, 0x44)
ROUGE      = RGBColor(0xC0, 0x00, 0x00)
NOIR       = RGBColor(0x1A, 0x1A, 0x1A)
BLEU_PALE  = RGBColor(0xD6, 0xE4, 0xF5)
ORANGE_PL  = RGBColor(0xFF, 0xEB, 0xD6)
VERT_PL    = RGBColor(0xE2, 0xEF, 0xDA)
ROUGE_PL   = RGBColor(0xFF, 0xE8, 0xE8)
JAUNE_PL   = RGBColor(0xFF, 0xF7, 0xCC)
GRIS_BLU   = RGBColor(0xEC, 0xF4, 0xFB)
GRIS_ORANGE= RGBColor(0xFF, 0xF3, 0xE8)

# ── Helpers ───────────────────────────────────────────────────────────────────
prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)
sl = prs.slides.add_slide(prs.slide_layouts[6])   # layout vide

def R(x, y, w, h, fill=None, line=None, lw=Pt(0)):
    shp = sl.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    else:
        shp.fill.background()
    if line:
        shp.line.color.rgb = line
        shp.line.width = lw if lw.pt else Pt(1)
    else:
        shp.line.fill.background()
    return shp

def T(text, x, y, w, h, sz=11, bold=False, col=NOIR,
      al=PP_ALIGN.LEFT, it=False):
    tb = sl.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    p  = tf.paragraphs[0]
    p.alignment = al
    run = p.add_run()
    run.text = text
    run.font.size    = Pt(sz)
    run.font.bold    = bold
    run.font.italic  = it
    run.font.color.rgb = col

def multiT(lines, x, y, w, h):
    """lines = [(text, sz, bold, col, align)]"""
    tb = sl.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    first = True
    for (txt, sz, bd, col, al) in lines:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.alignment = al
        run = p.add_run()
        run.text = txt
        run.font.size  = Pt(sz)
        run.font.bold  = bd
        run.font.color.rgb = col

# =============================================================================
# HEADER
# =============================================================================
R(0, 0, 13.33, 1.06, fill=BLEU_FONCE)
T("Modèle de Simeoni — Inhibition de la croissance tumorale (TGI)",
  0.3, 0.04, 12.6, 0.62, sz=24, bold=True, col=BLANC)
T("Simeoni et al., Cancer Research, 2004  ·  Modèle PK/PD mécanistique préclinique (in vivo)",
  0.3, 0.67, 12.6, 0.3, sz=10, col=RGBColor(0xBD, 0xD7, 0xEE))
R(0, 1.01, 13.33, 0.055, fill=ORANGE)

# =============================================================================
# COLONNE GAUCHE — Concept  (x : 0.15 → 3.55)
# =============================================================================
R(0.15, 1.15, 3.35, 3.75, fill=BLEU_PALE, line=BLEU_CLAIR, lw=Pt(1))
T("Pourquoi ce modèle ?", 0.3, 1.18, 3.1, 0.38, sz=13, bold=True, col=BLEU_FONCE)

# Bloc rouge — modèle simple
R(0.25, 1.63, 3.15, 0.72, fill=ROUGE_PL, line=ROUGE, lw=Pt(0.75))
T("Modèle Emax simple  ✗", 0.35, 1.65, 2.9, 0.28, sz=11, bold=True, col=ROUGE)
T("Effet instantané : pas de délai entre exposition et réponse tumorale",
  0.35, 1.91, 2.95, 0.42, sz=9, col=RGBColor(0x80, 0x00, 0x00), it=True)

# Bloc vert — Simeoni
R(0.25, 2.44, 3.15, 0.72, fill=VERT_PL, line=VERT, lw=Pt(0.75))
T("Modèle de Simeoni  ✓", 0.35, 2.46, 2.9, 0.28, sz=11, bold=True, col=VERT)
T("Compartiments de transit : délai biologique entre dommage et mort",
  0.35, 2.72, 2.95, 0.42, sz=9, col=RGBColor(0x15, 0x55, 0x20), it=True)

# Trois concepts clés
T("Ce que le modèle capture :", 0.3, 3.27, 3.1, 0.3, sz=11, bold=True, col=BLEU_FONCE)

concepts = [
    ("①  Croissance tumorale",
     "Exponentielle puis linéaire (sans traitement)"),
    ("②  Endommagement par le drug",
     "Les cellules lèsées transitent avant de mourir"),
    ("③  Compétition croissance / mort",
     "Régression si effet drug > seuil de croissance"),
]
for i, (titre, detail) in enumerate(concepts):
    yp = 3.60 + i * 0.68
    bg = BLANC if i % 2 == 0 else GRIS_BLU
    R(0.25, yp, 3.15, 0.62, fill=bg, line=BLEU_CLAIR, lw=Pt(0.5))
    T(titre,  0.33, yp + 0.03, 3.0, 0.27, sz=10, bold=True, col=BLEU_FONCE)
    T(detail, 0.33, yp + 0.30, 3.0, 0.30, sz=9,  col=NOIR,  it=True)

# =============================================================================
# COLONNE CENTRE — Schéma  (x : 3.70 → 9.75)
# =============================================================================
R(3.70, 1.15, 6.00, 3.75, fill=GRIS_CLAIR, line=BLEU_CLAIR, lw=Pt(1))
T("Structure des compartiments", 3.85, 1.18, 5.70, 0.38,
  sz=12, bold=True, col=BLEU_FONCE, al=PP_ALIGN.CENTER)

# ── Médicament (PK) ──
R(4.40, 1.65, 1.55, 0.44, fill=BLEU_CLAIR)
T("C(t)",    4.45, 1.67, 0.75, 0.36, sz=13, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
T("Médicament", 5.22, 1.70, 0.68, 0.28, sz=9, col=BLANC)

# Flèche verticale PK → x₁/x₂ junction (tige + pointe)
R(5.09, 2.09, 0.04, 0.40, fill=BLEU_CLAIR)                  # tige
T("▼", 4.99, 2.40, 0.22, 0.22, sz=9, col=BLEU_CLAIR, bold=True, al=PP_ALIGN.CENTER)
T("k₂·C", 5.14, 2.15, 0.5, 0.22, sz=9, col=BLEU_CLAIR, it=True)

# ── Compartiments  y = 2.65 ──
BY = 2.65      # box top
BH = 0.70      # box height

# x₁ — Proliférant (vert)
R(3.90, BY, 1.15, BH, fill=VERT)
T("x₁",      3.93, BY + 0.06, 1.08, 0.34, sz=17, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
T("Prolif.",  3.93, BY + 0.40, 1.08, 0.24, sz=9,  col=BLANC, al=PP_ALIGN.CENTER)
# Flèche de croissance au-dessus de x₁
T("↺  g(w)·x₁", 3.90, BY - 0.26, 1.15, 0.22, sz=9, col=VERT, bold=True, al=PP_ALIGN.CENTER)

# x₂, x₃, x₄ — Transit endommagé (orange)
for i, xi in enumerate(["x₂", "x₃", "x₄"]):
    xp = 5.25 + i * 1.18
    R(xp, BY, 1.00, BH, fill=ORANGE)
    T(xi,        xp + 0.02, BY + 0.06, 0.96, 0.34, sz=17, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
    T("Transit", xp + 0.02, BY + 0.40, 0.96, 0.24, sz=9,  col=BLANC, al=PP_ALIGN.CENTER)

# ✗ — Mort (rouge)
R(8.85, BY, 0.65, BH, fill=ROUGE)
T("✗",    8.87, BY + 0.06, 0.60, 0.36, sz=18, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
T("Mort", 8.87, BY + 0.40, 0.60, 0.24, sz=9,  col=BLANC, al=PP_ALIGN.CENTER)

# ── Flèches horizontales entre compartiments ──
# Positions : x₁ finit à 5.05 ; x₂ à 5.25 ; gaps = 0.20"
# x₂ finit à 6.25 ; x₃ à 6.43 ; gap = 0.18"
# x₃ finit à 7.43 ; x₄ à 7.61 ; gap = 0.18"
# x₄ finit à 8.61 ; ✗ à 8.85 ; gap = 0.24"
arrows = [(5.06, "k₂·C", BLEU_CLAIR),
          (6.27, "k₁",   RGBColor(0x80, 0x40, 0x00)),
          (7.45, "k₁",   RGBColor(0x80, 0x40, 0x00)),
          (8.63, "k₁",   RGBColor(0x80, 0x40, 0x00))]
for (ax, lbl, lcol) in arrows:
    T("→", ax, BY + 0.22, 0.20, 0.30, sz=14, col=BLEU_FONCE, bold=True, al=PP_ALIGN.CENTER)
    T(lbl, ax, BY - 0.06, 0.30, 0.18, sz=8,  col=lcol, it=True, al=PP_ALIGN.CENTER)

# ── Volume tumoral total ──
R(3.80, 3.47, 5.85, 0.32, fill=JAUNE_PL, line=ORANGE, lw=Pt(0.75))
T("Volume tumoral total :   w(t)  =  x₁ + x₂ + x₃ + x₄",
  3.90, 3.49, 5.65, 0.27, sz=11, bold=True, col=ORANGE, al=PP_ALIGN.CENTER)

# ── Loi de croissance ──
R(3.70, 3.88, 6.00, 0.98, fill=BLANC, line=VERT, lw=Pt(1.2))
T("Croissance sans traitement  g(w) :", 3.85, 3.91, 5.70, 0.30,
  sz=11, bold=True, col=VERT)
T("g(w)  =  λ₀  /  ( 1 + (λ₀·w / λ₁)^ψ )^(1/ψ)       avec  ψ = 20  (fixé)",
  3.85, 4.22, 5.70, 0.30, sz=12, bold=True, col=VERT, al=PP_ALIGN.CENTER, it=True)
T("Phase exponentielle (w petit)  ——→  Phase linéaire (w grand)",
  3.85, 4.54, 5.70, 0.26, sz=9.5, col=RGBColor(0x30, 0x70, 0x30),
  al=PP_ALIGN.CENTER, it=True)

# =============================================================================
# COLONNE DROITE — Équations + TSC  (x : 9.85 → 13.10)
# =============================================================================
R(9.85, 1.15, 3.20, 3.75, fill=BLANC, line=BLEU_CLAIR, lw=Pt(1))
T("Équations différentielles", 10.00, 1.18, 3.00, 0.38,
  sz=12, bold=True, col=BLEU_FONCE)

odes = [
    ("dx₁/dt  =  [ g(w) − k₂·C ] · x₁",  BLEU_PALE,  BLEU_FONCE, True),
    ("dx₂/dt  =  k₂·C·x₁  −  k₁·x₂",     GRIS_CLAIR, NOIR,       False),
    ("dx₃/dt  =  k₁·x₂    −  k₁·x₃",     BLANC,      NOIR,       False),
    ("dx₄/dt  =  k₁·x₃    −  k₁·x₄",     GRIS_CLAIR, NOIR,       False),
]
for i, (eq, bg, fc, bd) in enumerate(odes):
    yp = 1.65 + i * 0.46
    R(9.90, yp, 3.10, 0.42, fill=bg)
    T(eq, 9.95, yp + 0.05, 3.00, 0.34, sz=9.5, bold=bd, col=fc, it=True)

# TSC box
R(9.85, 3.51, 3.20, 1.36, fill=ORANGE_PL, line=ORANGE, lw=Pt(1.2))
T("TSC  (Tumor Static Conc.)", 9.95, 3.54, 3.05, 0.30, sz=11, bold=True, col=ORANGE)
T("Concentration seuil qui maintient\nla tumeur stable :", 9.95, 3.86, 3.05, 0.40, sz=10, col=NOIR)
T("TSC  ≈  k₁ / k₂", 9.95, 4.27, 3.05, 0.28,
  sz=14, bold=True, col=ROUGE, al=PP_ALIGN.CENTER)
T("C > TSC  →  régression tumorale",  9.95, 4.57, 3.05, 0.24, sz=9.5, col=ROUGE, bold=True,  al=PP_ALIGN.CENTER)
T("C < TSC  →  croissance ralentie",  9.95, 4.78, 3.05, 0.24, sz=9.5, col=RGBColor(0x80, 0x40, 0x00), al=PP_ALIGN.CENTER)

# =============================================================================
# BAS — Tableau des paramètres  (y : 5.05 → 6.95)
# =============================================================================
R(0.15, 5.05, 13.00, 1.90, fill=GRIS_CLAIR, line=BLEU_CLAIR, lw=Pt(1))
T("Paramètres clés", 0.30, 5.08, 3.50, 0.33, sz=12, bold=True, col=BLEU_FONCE)

# Colonnes : (x_start, width, header_label)
COLS = [
    (0.15, 1.10, "Param."),
    (1.30, 2.90, "Définition"),
    (4.25, 1.45, "Unité"),
    (5.75, 4.30, "Interprétation biologique"),
    (10.10, 2.90, "Estimé depuis"),
]

# En-tête
for (cx, cw, lbl) in COLS:
    R(cx, 5.44, cw, 0.31, fill=BLEU_FONCE)
    T(lbl, cx + 0.04, 5.46, cw - 0.08, 0.27,
      sz=10, bold=True, col=BLANC, al=PP_ALIGN.CENTER)

PARAMS = [
    ("λ₀",   "Taux de croissance exponentielle",  "j⁻¹",    "Vitesse de doublement initial de la tumeur (contrôle)",     "Données contrôle"),
    ("λ₁",   "Taux de croissance linéaire",        "mm³/j",  "Taille limite à saturation — freine la croissance",          "Données contrôle"),
    ("k₂",   "Constante cytotoxique du drug",      "L/µg/j", "Rate d'entrée des cellules en phase endommagée  [= k₂·C]",  "Données traitées + PK"),
    ("k₁",   "Constante de transit",               "j⁻¹",    "= 3/MTT (mean transit time) — durée de la phase de transit", "Données biologiques"),
]
for row_i, (par, defn, unit, interp, src) in enumerate(PARAMS):
    yp = 5.77 + row_i * 0.285
    bg = BLANC if row_i % 2 == 0 else GRIS_BLU
    vals = [par, defn, unit, interp, src]
    for (cx, cw, _), val in zip(COLS, vals):
        R(cx, yp, cw, 0.27, fill=bg, line=RGBColor(0xCC, 0xCC, 0xCC), lw=Pt(0.25))
        is_param = (val == par)
        T(val, cx + 0.04, yp + 0.02, cw - 0.08, 0.23,
          sz=9, bold=is_param,
          col=BLEU_FONCE if is_param else NOIR,
          al=PP_ALIGN.CENTER if is_param else PP_ALIGN.LEFT)

# =============================================================================
# FOOTER
# =============================================================================
R(0, 7.20, 13.33, 0.30, fill=BLEU_FONCE)
T("Réf. : Simeoni M. et al.  —  A Predictive PK-PD Model of Tumor Growth Inhibition in Mouse Xenograft Experiments.  Cancer Research, 2004, 64 : 1094–1101",
  0.20, 7.22, 10.80, 0.26, sz=8, col=BLANC)
T("Réunion FGFR2", 11.20, 7.22, 1.90, 0.26, sz=9, col=BLANC, al=PP_ALIGN.RIGHT)

# =============================================================================
# SAUVEGARDE
# =============================================================================
out = "scripts/slide_simeoni.pptx"
prs.save(out)
print(f"Slide créée : {out}")
