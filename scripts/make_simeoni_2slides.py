"""
2 slides didactiques — Modèle de Simeoni (2004)
Slide 1 : Concept & Structure
Slide 2 : Équations & Paramètres
Même charte graphique que make_pptx.py
"""
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

# ── Palette ───────────────────────────────────────────────────────────────────
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

# ── Helpers (chaque fonction reçoit le slide en premier argument) ─────────────
prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)

def new_slide():
    return prs.slides.add_slide(prs.slide_layouts[6])

def R(sl, x, y, w, h, fill=None, line=None, lw=Pt(0)):
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

def T(sl, text, x, y, w, h, sz=11, bold=False, col=NOIR,
      al=PP_ALIGN.LEFT, it=False):
    tb = sl.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    p  = tf.paragraphs[0]
    p.alignment = al
    run = p.add_run()
    run.text      = text
    run.font.size = Pt(sz)
    run.font.bold = bold
    run.font.italic = it
    run.font.color.rgb = col

def header(sl, titre, sous_titre, page):
    R(sl, 0, 0, 13.33, 1.06, fill=BLEU_FONCE)
    T(sl, titre,      0.30, 0.04, 11.8, 0.62, sz=22, bold=True, col=BLANC)
    T(sl, sous_titre, 0.30, 0.67, 11.5, 0.30, sz=10, col=RGBColor(0xBD, 0xD7, 0xEE))
    T(sl, page, 12.45, 0.10, 0.75, 0.55, sz=20, bold=True,
      col=ORANGE, al=PP_ALIGN.RIGHT)
    R(sl, 0, 1.01, 13.33, 0.055, fill=ORANGE)

def footer(sl):
    R(sl, 0, 7.20, 13.33, 0.30, fill=BLEU_FONCE)
    T(sl, "Réf. : Simeoni M. et al.  —  Cancer Research, 2004, 64 : 1094–1101",
      0.20, 7.22, 10.5, 0.26, sz=8, col=BLANC)
    T(sl, "Réunion FGFR2", 11.2, 7.22, 1.9, 0.26, sz=9,
      col=BLANC, al=PP_ALIGN.RIGHT)


# =============================================================================
# SLIDE 1  —  Concept & Structure
# =============================================================================
sl1 = new_slide()
header(sl1,
       "Modèle de Simeoni (2004)  —  1/2 : Concept & Structure",
       "Pourquoi modéliser le délai entre l'exposition au médicament et la réponse tumorale ?",
       "1/2")

# ── Colonne gauche : Motivation  (0.15 → 4.30) ───────────────────────────────
R(sl1, 0.15, 1.15, 4.05, 5.95, fill=BLEU_PALE, line=BLEU_CLAIR, lw=Pt(1))
T(sl1, "Pourquoi ce modèle ?", 0.28, 1.18, 3.80, 0.38,
  sz=14, bold=True, col=BLEU_FONCE)

# Emax simple — rouge
R(sl1, 0.28, 1.65, 3.80, 0.90, fill=ROUGE_PL, line=ROUGE, lw=Pt(0.75))
T(sl1, "Modèle Emax simple  ✗", 0.38, 1.67, 3.60, 0.30, sz=12, bold=True, col=ROUGE)
T(sl1, "dV/dt = (kg − ke·C) · V",
  0.38, 1.97, 3.60, 0.26, sz=11, col=RGBColor(0x80,0x00,0x00), it=True, al=PP_ALIGN.CENTER)
T(sl1, "Effet instantané — ne capture pas le délai pharmacologique",
  0.38, 2.22, 3.60, 0.26, sz=9, col=RGBColor(0x80,0x00,0x00), it=True)

# Simeoni — vert
R(sl1, 0.28, 2.66, 3.80, 0.90, fill=VERT_PL, line=VERT, lw=Pt(0.75))
T(sl1, "Modèle de Simeoni  ✓", 0.38, 2.68, 3.60, 0.30, sz=12, bold=True, col=VERT)
T(sl1, "x₁ → x₂ → x₃ → x₄  (compartiments de transit)",
  0.38, 2.98, 3.60, 0.26, sz=11, col=RGBColor(0x15,0x55,0x20), it=True, al=PP_ALIGN.CENTER)
T(sl1, "Délai biologique entre dommage cellulaire et mort",
  0.38, 3.23, 3.60, 0.26, sz=9, col=RGBColor(0x15,0x55,0x20), it=True)

# 3 concepts clés
T(sl1, "Ce que capture le modèle :", 0.28, 3.70, 3.80, 0.30,
  sz=12, bold=True, col=BLEU_FONCE)

CONCEPTS = [
    ("①  Croissance tumorale",
     "Phase exponentielle (petit volume)\npuis linéaire (grand volume)"),
    ("②  Endommagement par le drug",
     "k₂·C = taux d'entrée en phase\nendommagée, proportionnel à C(t)"),
    ("③  Compétition croissance / mort",
     "Régression si l'effet du drug\ndépasse le seuil TSC ≈ k₁/k₂"),
]
for i, (titre, detail) in enumerate(CONCEPTS):
    yp = 4.08 + i * 0.96
    bg = BLANC if i % 2 == 0 else GRIS_BLU
    R(sl1, 0.28, yp, 3.80, 0.88, fill=bg, line=BLEU_CLAIR, lw=Pt(0.5))
    T(sl1, titre,  0.36, yp + 0.04, 3.65, 0.28, sz=11, bold=True, col=BLEU_FONCE)
    T(sl1, detail, 0.36, yp + 0.34, 3.65, 0.50, sz=9.5, col=NOIR, it=True)

# ── Colonne droite : Schéma  (4.45 → 13.10) ──────────────────────────────────
R(sl1, 4.45, 1.15, 8.50, 5.95, fill=GRIS_CLAIR, line=BLEU_CLAIR, lw=Pt(1))
T(sl1, "Structure des compartiments", 4.60, 1.18, 8.20, 0.38,
  sz=13, bold=True, col=BLEU_FONCE, al=PP_ALIGN.CENTER)

# Médicament C(t)
R(sl1, 5.35, 1.68, 2.30, 0.50, fill=BLEU_CLAIR)
T(sl1, "C(t)", 5.40, 1.70, 1.15, 0.40, sz=17, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
T(sl1, "Médicament", 6.58, 1.73, 1.02, 0.32, sz=10, col=BLANC)

# Flèche verticale PK → transition x₁/x₂
R(sl1, 6.44, 2.18, 0.045, 0.44, fill=BLEU_CLAIR)
T(sl1, "▼", 6.33, 2.50, 0.25, 0.22, sz=10, col=BLEU_CLAIR, bold=True, al=PP_ALIGN.CENTER)
T(sl1, "k₂·C", 6.50, 2.22, 0.65, 0.24, sz=10, bold=True, col=BLEU_CLAIR, it=True)

# Boîtes compartiments  (BY = 2.78, BH = 1.05)
BY, BH = 2.78, 1.05

# x₁ Proliférant (vert)
R(sl1, 4.62, BY, 1.60, BH, fill=VERT)
T(sl1, "x₁",      4.65, BY + 0.08, 1.52, 0.42, sz=24, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
T(sl1, "Prolif.",  4.65, BY + 0.58, 1.52, 0.30, sz=11, col=BLANC, al=PP_ALIGN.CENTER)
T(sl1, "↺  λ₀ · g(w)", 4.62, BY - 0.30, 1.60, 0.26, sz=10,
  col=VERT, bold=True, al=PP_ALIGN.CENTER)

# x₂, x₃, x₄ Transit (orange)
# x₁ finit à 6.22 → gap → x₂ débute à 6.47
for i, xi in enumerate(["x₂", "x₃", "x₄"]):
    xp = 6.47 + i * 1.57
    R(sl1, xp, BY, 1.40, BH, fill=ORANGE)
    T(sl1, xi,        xp + 0.03, BY + 0.08, 1.33, 0.42, sz=24, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
    T(sl1, "Transit", xp + 0.03, BY + 0.58, 1.33, 0.30, sz=11, col=BLANC, al=PP_ALIGN.CENTER)

# ✗ Mort (rouge)
# x₄ finit à 6.47 + 2*1.57 + 1.40 = 11.01 → gap → ✗ à 11.25
R(sl1, 11.25, BY, 0.95, BH, fill=ROUGE)
T(sl1, "✗",    11.27, BY + 0.08, 0.88, 0.42, sz=24, bold=True, col=BLANC, al=PP_ALIGN.CENTER)
T(sl1, "Mort", 11.27, BY + 0.58, 0.88, 0.30, sz=11, col=BLANC, al=PP_ALIGN.CENTER)

# Flèches horizontales
# x₁(6.22) → x₂(6.47) | x₂(7.87) → x₃(8.04) | x₃(9.44) → x₄(9.61) | x₄(11.01) → ✗(11.25)
ARROWS = [
    (6.22, "k₂·C", BLEU_CLAIR),
    (7.87, "k₁",   RGBColor(0x80, 0x40, 0x00)),
    (9.44, "k₁",   RGBColor(0x80, 0x40, 0x00)),
    (11.01, "k₁",  RGBColor(0x80, 0x40, 0x00)),
]
for (ax, lbl, lcol) in ARROWS:
    T(sl1, "→", ax, BY + 0.30, 0.25, 0.44, sz=17, col=BLEU_FONCE, bold=True, al=PP_ALIGN.CENTER)
    T(sl1, lbl,  ax, BY - 0.10, 0.38, 0.22, sz=9,  col=lcol, it=True, al=PP_ALIGN.CENTER)

# w total
R(sl1, 4.55, BY + BH + 0.08, 7.65, 0.35, fill=JAUNE_PL, line=ORANGE, lw=Pt(0.75))
T(sl1, "Volume tumoral total :   w(t)  =  x₁ + x₂ + x₃ + x₄",
  4.65, BY + BH + 0.10, 7.45, 0.28, sz=12, bold=True, col=ORANGE, al=PP_ALIGN.CENTER)

# g(w) loi de croissance
GY = BY + BH + 0.55
R(sl1, 4.45, GY, 8.50, 1.42, fill=BLANC, line=VERT, lw=Pt(1.2))
T(sl1, "Loi de croissance sans traitement :",
  4.60, GY + 0.04, 8.20, 0.30, sz=12, bold=True, col=VERT)
T(sl1, "g(w)  =  λ₀  /  ( 1 + (λ₀·w / λ₁)^ψ )^(1/ψ)       [ψ = 20, fixé]",
  4.60, GY + 0.37, 8.20, 0.35, sz=14, bold=True, col=VERT, al=PP_ALIGN.CENTER, it=True)
T(sl1, "w petit  →  g(w) ≈ λ₀  (exponentielle)           w grand  →  g(w) ≈ λ₁/w  (linéaire)",
  4.60, GY + 0.76, 8.20, 0.28, sz=10, col=RGBColor(0x20,0x60,0x20), it=True, al=PP_ALIGN.CENTER)
T(sl1, "Représente la croissance naturelle observée dans le groupe contrôle (sans traitement)",
  4.60, GY + 1.06, 8.20, 0.26, sz=9, col=RGBColor(0x30,0x50,0x30), al=PP_ALIGN.CENTER)

# Note biologique
NY = GY + 1.56
R(sl1, 4.45, NY, 8.50, 0.85, fill=GRIS_BLU, line=BLEU_CLAIR, lw=Pt(0.75))
T(sl1, "Interprétation du délai — MTT (Mean Transit Time)",
  4.60, NY + 0.04, 8.20, 0.28, sz=11, bold=True, col=BLEU_FONCE)
T(sl1, "Le MTT = 3/k₁ représente la durée moyenne entre l'endommagement d'une cellule par le médicament "
       "et sa mort effective. Ce paramètre capture le délai pharmacodynamique observé en préclinique.",
  4.60, NY + 0.34, 8.20, 0.44, sz=9.5, col=NOIR)

footer(sl1)


# =============================================================================
# SLIDE 2  —  Équations & Paramètres
# =============================================================================
sl2 = new_slide()
header(sl2,
       "Modèle de Simeoni (2004)  —  2/2 : Équations & Paramètres",
       "Système différentiel couplé  —  interprétation pharmacologique des paramètres",
       "2/2")

# ── Gauche : ODEs (0.15 → 6.55) ──────────────────────────────────────────────
R(sl2, 0.15, 1.15, 6.30, 3.75, fill=BLANC, line=BLEU_CLAIR, lw=Pt(1))
T(sl2, "Système d'équations différentielles", 0.30, 1.18, 6.00, 0.38,
  sz=13, bold=True, col=BLEU_FONCE)

ODES = [
    ("dx₁/dt  =  [ g(w) − k₂·C ] · x₁",  BLEU_PALE,  BLEU_FONCE, True,  "Prolif. + dommage drug"),
    ("dx₂/dt  =  k₂·C · x₁  −  k₁ · x₂", GRIS_CLAIR, NOIR,       False, "Transit 1"),
    ("dx₃/dt  =  k₁ · x₂    −  k₁ · x₃", BLANC,      NOIR,       False, "Transit 2"),
    ("dx₄/dt  =  k₁ · x₃    −  k₁ · x₄", GRIS_CLAIR, NOIR,       False, "Transit 3  →  mort"),
]
for i, (eq, bg, fc, bd, note) in enumerate(ODES):
    yp = 1.68 + i * 0.58
    R(sl2, 0.20, yp, 6.20, 0.52, fill=bg)
    T(sl2, eq,   0.28, yp + 0.07, 4.30, 0.36, sz=12, bold=bd, col=fc, it=True)
    T(sl2, "← " + note, 4.62, yp + 0.10, 1.70, 0.28,
      sz=9, col=RGBColor(0x50,0x50,0x50), it=True)

# w total
W_Y = 1.68 + 4 * 0.58   # = 4.00
R(sl2, 0.20, W_Y, 6.20, 0.42, fill=JAUNE_PL, line=ORANGE, lw=Pt(0.75))
T(sl2, "w(t)  =  x₁ + x₂ + x₃ + x₄       (volume tumoral total)",
  0.30, W_Y + 0.05, 6.00, 0.30, sz=11, bold=True, col=ORANGE, al=PP_ALIGN.CENTER, it=True)

# g(w)
R(sl2, 0.15, W_Y + 0.52, 6.30, 0.85, fill=VERT_PL, line=VERT, lw=Pt(1))
T(sl2, "g(w)  =  λ₀  /  ( 1 + (λ₀·w / λ₁)^ψ )^(1/ψ)      [ψ = 20]",
  0.28, W_Y + 0.55, 6.00, 0.33, sz=13, bold=True, col=VERT, al=PP_ALIGN.CENTER, it=True)
T(sl2, "Exponentielle (w petit)  →  Linéaire (w grand)",
  0.28, W_Y + 0.90, 6.00, 0.26, sz=9.5, col=RGBColor(0x20,0x60,0x20),
  it=True, al=PP_ALIGN.CENTER)

# ── Droite : TSC + lien modèle (6.75 → 13.10) ────────────────────────────────
R(sl2, 6.75, 1.15, 6.20, 3.75, fill=BLANC, line=BLEU_CLAIR, lw=Pt(1))
T(sl2, "Quantité dérivée clé", 6.90, 1.18, 5.90, 0.38,
  sz=13, bold=True, col=BLEU_FONCE)

# TSC
R(sl2, 6.85, 1.65, 6.00, 2.15, fill=ORANGE_PL, line=ORANGE, lw=Pt(1.5))
T(sl2, "TSC — Tumor Static Concentration",
  6.95, 1.68, 5.80, 0.32, sz=12, bold=True, col=ORANGE)
T(sl2, "Concentration seuil pour laquelle la tumeur est stable : croissance = inhibition",
  6.95, 2.02, 5.80, 0.40, sz=10, col=NOIR)
T(sl2, "TSC  ≈  k₁ / k₂",
  6.95, 2.46, 5.80, 0.40, sz=20, bold=True, col=ROUGE, al=PP_ALIGN.CENTER)

R(sl2, 6.90, 2.93, 2.88, 0.36, fill=ROUGE_PL, line=ROUGE, lw=Pt(0.5))
T(sl2, "C > TSC  →  régression", 6.97, 2.95, 2.72, 0.28,
  sz=10, bold=True, col=ROUGE)

R(sl2, 9.87, 2.93, 2.88, 0.36, fill=VERT_PL, line=VERT, lw=Pt(0.5))
T(sl2, "C < TSC  →  ralentissement", 9.94, 2.95, 2.72, 0.28,
  sz=10, bold=True, col=VERT)

# Lien avec notre modèle
R(sl2, 6.85, 3.40, 6.00, 1.47, fill=GRIS_BLU, line=BLEU_CLAIR, lw=Pt(0.75))
T(sl2, "Lien avec le modèle appliqué (FGFR2)",
  6.95, 3.43, 5.80, 0.30, sz=11, bold=True, col=BLEU_FONCE)
T(sl2, "• kg (taux de croissance) est l'équivalent de λ₀",
  6.95, 3.76, 5.80, 0.27, sz=10, col=NOIR)
T(sl2, "• ke · C/(EC50+C) remplace k₂·C  (modèle Emax, sans transit)",
  6.95, 4.03, 5.80, 0.27, sz=10, col=NOIR)
T(sl2, "• TSC correspond à la concentration statique de notre modèle",
  6.95, 4.30, 5.80, 0.27, sz=10, col=NOIR)
T(sl2, "→ Le modèle Emax est un cas simplifié de Simeoni (absence de délai de transit)",
  6.95, 4.60, 5.80, 0.28, sz=9.5, col=RGBColor(0x30,0x30,0x80), it=True)

# ── Bas : Tableau des paramètres (full width) ─────────────────────────────────
R(sl2, 0.15, 5.06, 12.95, 2.02, fill=GRIS_CLAIR, line=BLEU_CLAIR, lw=Pt(1))
T(sl2, "Paramètres du modèle", 0.30, 5.09, 4.0, 0.32,
  sz=12, bold=True, col=BLEU_FONCE)

COLS = [
    (0.15, 1.05, "Param."),
    (1.25, 2.75, "Définition"),
    (4.05, 1.35, "Unité"),
    (5.45, 4.40, "Interprétation pharmacologique"),
    (9.90, 3.10, "Source d'estimation"),
]
for (cx, cw, lbl) in COLS:
    R(sl2, cx, 5.44, cw, 0.32, fill=BLEU_FONCE)
    T(sl2, lbl, cx+0.04, 5.46, cw-0.08, 0.27,
      sz=10, bold=True, col=BLANC, al=PP_ALIGN.CENTER)

PARAMS = [
    ("λ₀", "Taux de croissance exponentielle", "j⁻¹",    "Vitesse initiale de doublement tumoral",               "Groupe contrôle"),
    ("λ₁", "Taux de croissance linéaire",      "mm³/j",  "Vitesse de croissance à saturation (grand volume)",    "Groupe contrôle"),
    ("k₂", "Paramètre cytotoxique du drug",    "L/µg/j", "Rate d'entrée en phase endommagée = k₂ × C(t)",       "Données traitées + PK"),
    ("k₁", "Constante de transit",             "j⁻¹",    "= 3/MTT  —  délai moyen avant mort cellulaire",       "Données biologiques"),
]
for row_i, (par, defn, unit, interp, src) in enumerate(PARAMS):
    yp = 5.78 + row_i * 0.295
    bg = BLANC if row_i % 2 == 0 else GRIS_BLU
    vals = [par, defn, unit, interp, src]
    for (cx, cw, _), val in zip(COLS, vals):
        R(sl2, cx, yp, cw, 0.27, fill=bg,
          line=RGBColor(0xCC,0xCC,0xCC), lw=Pt(0.25))
        is_p = (val == par)
        T(sl2, val, cx+0.04, yp+0.02, cw-0.08, 0.23,
          sz=9, bold=is_p,
          col=BLEU_FONCE if is_p else NOIR,
          al=PP_ALIGN.CENTER if is_p else PP_ALIGN.LEFT)

footer(sl2)

# =============================================================================
# SAUVEGARDE
# =============================================================================
out = "scripts/slide_simeoni_2slides.pptx"
prs.save(out)
print(f"Fichier créé : {out}")
