from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import datetime

# ── Palette Pierre Fabre-ish ─────────────────────────────────────────────────
BLEU       = RGBColor(0x00, 0x33, 0x6B)   # bleu marine
BLEU_CLAIR = RGBColor(0x00, 0x7A, 0xC2)   # bleu accent
BLANC      = RGBColor(0xFF, 0xFF, 0xFF)
GRIS       = RGBColor(0xF2, 0xF2, 0xF2)
VERT       = RGBColor(0x21, 0x8B, 0x4E)   # succès
ORANGE     = RGBColor(0xE8, 0x7B, 0x10)   # attention

W  = Inches(13.33)   # 16:9 widescreen
H  = Inches(7.5)

prs = Presentation()
prs.slide_width  = W
prs.slide_height = H

BLANK = prs.slide_layouts[6]   # complètement vide


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def bg(slide, color):
    fill = slide.background.fill
    fill.solid()
    fill.fore_color.rgb = color

def box(slide, text, x, y, w, h,
        font_size=18, bold=False, italic=False,
        fg=BLANC, bg_color=None, align=PP_ALIGN.LEFT, wrap=True):
    txb = slide.shapes.add_textbox(x, y, w, h)
    tf  = txb.text_frame
    tf.word_wrap = wrap
    p   = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size  = Pt(font_size)
    run.font.bold  = bold
    run.font.color.rgb = fg
    if italic:
        run.font.italic = italic
    if bg_color:
        txb.fill.solid()
        txb.fill.fore_color.rgb = bg_color
    return txb

def rect(slide, x, y, w, h, color):
    shp = slide.shapes.add_shape(1, x, y, w, h)   # MSO_SHAPE_TYPE.RECTANGLE
    shp.fill.solid()
    shp.fill.fore_color.rgb = color
    shp.line.fill.background()
    return shp

def bullet_box(slide, items, x, y, w, h,
               font_size=16, fg=BLEU, marker="●  "):
    txb = slide.shapes.add_textbox(x, y, w, h)
    tf  = txb.text_frame
    tf.word_wrap = True
    for i, (txt, sub) in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.space_before = Pt(4)
        run = p.add_run()
        run.text = marker + txt
        run.font.size  = Pt(font_size)
        run.font.bold  = True
        run.font.color.rgb = fg
        if sub:
            p2 = tf.add_paragraph()
            p2.space_before = Pt(1)
            r2 = p2.add_run()
            r2.text = "      " + sub
            r2.font.size  = Pt(font_size - 2)
            r2.font.color.rgb = RGBColor(0x44, 0x44, 0x44)

def header_bar(slide, title, subtitle=""):
    rect(slide, 0, 0, W, Inches(1.35), BLEU)
    box(slide, title,
        Inches(0.45), Inches(0.15), Inches(12), Inches(0.75),
        font_size=28, bold=True, fg=BLANC)
    if subtitle:
        box(slide, subtitle,
            Inches(0.45), Inches(0.88), Inches(12), Inches(0.4),
            font_size=14, italic=True, fg=BLEU_CLAIR)

def footer(slide, txt="Stage Pierre Fabre — Suivi PKPD | Mars 2026"):
    rect(slide, 0, Inches(7.15), W, Inches(0.35), BLEU)
    box(slide, txt,
        Inches(0.3), Inches(7.17), Inches(12.5), Inches(0.3),
        font_size=9, fg=BLANC, align=PP_ALIGN.CENTER)

def check_badge(slide, x, y, txt, color=VERT):
    rect(slide, x, y, Inches(0.45), Inches(0.38), color)
    box(slide, "✔", x + Inches(0.07), y, Inches(0.35), Inches(0.38),
        font_size=14, bold=True, fg=BLANC, align=PP_ALIGN.CENTER)
    box(slide, txt, x + Inches(0.55), y, Inches(5.0), Inches(0.38),
        font_size=13, bold=True, fg=BLEU)


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 1 — TITRE
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, BLEU)

# Bande blanche centrale
rect(s, Inches(0.6), Inches(1.8), Inches(12.1), Inches(3.8), BLANC)

# Titre principal
box(s, "Réunion de suivi — Stage Pierre Fabre",
    Inches(0.9), Inches(2.0), Inches(11.5), Inches(1.0),
    font_size=32, bold=True, fg=BLEU, align=PP_ALIGN.CENTER)

box(s, "Modélisation PKPD de l'hématotoxicité — ADC Pierre Fabre",
    Inches(0.9), Inches(3.0), Inches(11.5), Inches(0.7),
    font_size=20, fg=BLEU_CLAIR, align=PP_ALIGN.CENTER)

box(s, "Mars 2026",
    Inches(0.9), Inches(3.85), Inches(11.5), Inches(0.5),
    font_size=16, italic=True, fg=BLEU, align=PP_ALIGN.CENTER)

# Bande bleue inférieure
rect(s, 0, Inches(6.9), W, Inches(0.6), BLEU_CLAIR)
box(s, "Pipeline : Souris  →  Rat  →  Humain  |  Modèle de Friberg et al.",
    Inches(0.5), Inches(6.93), Inches(12.3), Inches(0.45),
    font_size=12, fg=BLANC, align=PP_ALIGN.CENTER)


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 2 — RÉCAPITULATIF DES AVANCÉES
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, GRIS)
header_bar(s, "Récapitulatif des avancées",
           "Toutes les étapes franchies depuis le début du stage")
footer(s)

# 2 colonnes de statut
cols = [
    (Inches(0.4),  "Modèle Friberg — Rat",     VERT,
     ["Carboplatin 40 mg/kg Q14D×8",
      "8 lignées cellulaires modélisées",
      "Ajustement modèle/données ✔",
      "125 jours de cinétique reproduits"]),
    (Inches(6.8),  "Transfert Rat → Humain",   VERT,
     ["Paramètres PK transférés depuis rat",
      "Carboplatin Q21D × 2 (clinique)",
      "Nadirs neutrophiles & plaquettes ✔",
      "Chaîne souris → rat → humain validée"]),
]

for cx, titre, col, items in cols:
    rect(s, cx, Inches(1.5), Inches(6.0), Inches(0.5), col)
    box(s, titre, cx + Inches(0.1), Inches(1.52), Inches(5.8), Inches(0.45),
        font_size=13, bold=True, fg=BLANC)
    for i, item in enumerate(items):
        box(s, "▸  " + item,
            cx + Inches(0.1),
            Inches(2.15) + Inches(0.55) * i,
            Inches(5.8), Inches(0.5),
            font_size=12, fg=BLEU)


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 3 — MODÈLE FRIBERG : RAT
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, GRIS)
header_bar(s, "Modèle Friberg — Hématopoïèse Rat",
           "Carboplatin 40 mg/kg Q14D × 8 | 8 lignées cellulaires")
footer(s)

box(s, "Ce qui a été modélisé",
    Inches(0.4), Inches(1.55), Inches(6.0), Inches(0.4),
    font_size=14, bold=True, fg=BLEU)

lignes = [
    ("Progéniteurs",  "MPP → CMP → MEP (hiérarchie hématopoïétique complète)"),
    ("Myéloïde",      "Neutrophiles  +  Monocytes"),
    ("Érythroïde",    "Réticulocytes  +  Globules rouges"),
    ("Thrombocytaire","Plaquettes"),
]
for i, (cat, detail) in enumerate(lignes):
    y = Inches(2.05) + Inches(0.75) * i
    rect(s, Inches(0.4), y, Inches(2.2), Inches(0.55), BLEU_CLAIR)
    box(s, cat, Inches(0.5), y + Inches(0.1), Inches(2.0), Inches(0.38),
        font_size=12, bold=True, fg=BLANC)
    box(s, detail, Inches(2.75), y + Inches(0.1), Inches(9.5), Inches(0.4),
        font_size=12, fg=BLEU)

# Résultats clés
rect(s, Inches(0.4), Inches(5.2), Inches(12.5), Inches(0.04), BLEU_CLAIR)
box(s, "Résultats clés", Inches(0.4), Inches(5.3), Inches(4.0), Inches(0.4),
    font_size=14, bold=True, fg=BLEU)

resultats = [
    "✔  Oscillations cycliques reproduites après chaque injection",
    "✔  Cinétique sur 125 jours cohérente avec les données expérimentales",
    "✔  Modèle PK → Damage → Hématopoïèse validé sur rat",
    "✔  Base solide pour le transfert rat → humain",
]
for i, r in enumerate(resultats):
    box(s, r, Inches(0.5), Inches(5.8) + Inches(0.28) * i, Inches(12.0), Inches(0.28),
        font_size=12, fg=RGBColor(0x11, 0x55, 0x11) if i < 3 else BLEU)


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 5 — TRANSFERT RAT → HUMAIN
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, BLANC)
header_bar(s, "Transfert inter-espèce : Rat → Humain",
           "Carboplatin Q21D × 2 | Validation clinique")
footer(s)

box(s, "Démarche", Inches(0.4), Inches(1.55), Inches(12.0), Inches(0.4),
    font_size=14, bold=True, fg=BLEU)

steps = [
    ("1", BLEU_CLAIR, "Paramètres PK estimés sur rat",
     "CL, V1, V2, Q transférés directement au modèle humain"),
    ("2", BLEU_CLAIR, "Simulation humaine",
     "Carboplatin AUC=5 (750 mg pour GFR=125) — 2 cycles Q21D"),
    ("3", VERT,       "Comparaison données cliniques",
     "Courbes modèle superposées aux observations (neutrophiles + plaquettes)"),
    ("4", VERT,       "Validation",
     "Nadirs bien reproduits — variabilité inter-patient identifiée"),
]
for i, (num, col, titre, desc) in enumerate(steps):
    y = Inches(2.1) + Inches(1.1) * i
    rect(s, Inches(0.4), y, Inches(0.6), Inches(0.6), col)
    box(s, num, Inches(0.4), y, Inches(0.6), Inches(0.6),
        font_size=16, bold=True, fg=BLANC, align=PP_ALIGN.CENTER)
    box(s, titre, Inches(1.15), y, Inches(11.0), Inches(0.35),
        font_size=13, bold=True, fg=BLEU)
    box(s, desc, Inches(1.15), y + Inches(0.35), Inches(11.0), Inches(0.35),
        font_size=11, fg=RGBColor(0x44, 0x44, 0x44))

# Box résultat
rect(s, Inches(0.4), Inches(6.35), Inches(12.5), Inches(0.7), GRIS)
rect(s, Inches(0.4), Inches(6.35), Inches(0.15), Inches(0.7), VERT)
box(s, "Chaîne complète validée :  Souris  →  Rat  →  Humain  ✔",
    Inches(0.7), Inches(6.45), Inches(12.0), Inches(0.45),
    font_size=14, bold=True, fg=BLEU, align=PP_ALIGN.CENTER)


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 6 — VPC & ANALYSE DE SENSIBILITÉ
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, GRIS)
header_bar(s, "VPC — Modèle Humain",
           "Carboplatin AUC=5, Q21D × 2 | 1000 patients simulés")
footer(s)

box(s, "VPC (Visual Predictive Check) — AUC=5, GFR=125 fixe",
    Inches(0.4), Inches(1.55), Inches(12.0), Inches(0.4),
    font_size=14, bold=True, fg=BLEU)

vpc_items = [
    "1000 patients simulés (Supp. S11, GFR=125 fixe)",
    "Intervalles de prédiction 5–95% calculés pour neutrophiles et plaquettes",
    "Erreur résiduelle intégrée dans la simulation populationnelle",
    "Nadirs cycliques bien capturés par les intervalles de prédiction",
]
for i, it in enumerate(vpc_items):
    box(s, "▸  " + it, Inches(0.6), Inches(2.1) + Inches(0.6) * i,
        Inches(12.0), Inches(0.5), font_size=13, fg=BLEU)

rect(s, Inches(0.4), Inches(4.7), Inches(12.5), Inches(0.04), BLEU_CLAIR)

box(s, "Grades NCI-CTCAE — Figure 4c",
    Inches(0.4), Inches(4.85), Inches(12.0), Inches(0.4),
    font_size=14, bold=True, fg=BLEU)

grades = [
    ("Neutrophiles", "Nadir = 2.30 × 10⁹/L  →  Grade 0"),
    ("Plaquettes",   "Nadir = 190 × 10⁹/L   →  Grade 0"),
]
for i, (cell, grade) in enumerate(grades):
    cx = Inches(0.4) + Inches(6.2) * i
    rect(s, cx, Inches(5.35), Inches(5.8), Inches(0.5), VERT)
    box(s, cell, cx + Inches(0.15), Inches(5.38), Inches(5.5), Inches(0.42),
        font_size=13, bold=True, fg=BLANC)
    box(s, grade, cx + Inches(0.15), Inches(5.95), Inches(5.5), Inches(0.4),
        font_size=12, fg=BLEU)


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 7 — PROCHAINE ÉTAPE
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, BLANC)
header_bar(s, "Prochaine étape — Application à l'ADC Pierre Fabre",
           "Étape 5 du pipeline : données précliniques confidentielles")
footer(s)

box(s, "Script R validé — prêt pour les données réelles",
    Inches(0.4), Inches(1.55), Inches(12.0), Inches(0.45),
    font_size=14, bold=True, fg=VERT)

etapes = [
    ("Étape 5", BLEU,
     "Application aux données précliniques de l'ADC",
     "Chargement des données souris confidentielle\n→ Estimation k₁, k₂ spécifiques à la molécule"),
    ("Validation", BLEU_CLAIR,
     "Vérification de la cohérence biologique",
     "Différenciation des doses, comparaison aux valeurs attendues\n→ Même démarche que la validation sur irinotécan"),
    ("Transfert", VERT,
     "Prédiction de la toxicité humaine",
     "Pipeline souris → rat → humain appliqué à l'ADC\n→ Prédiction des nadirs hématologiques en clinique"),
]

for i, (tag, col, titre, desc) in enumerate(etapes):
    y = Inches(2.2) + Inches(1.5) * i
    rect(s, Inches(0.4), y, Inches(1.4), Inches(0.55), col)
    box(s, tag, Inches(0.4), y, Inches(1.4), Inches(0.55),
        font_size=12, bold=True, fg=BLANC, align=PP_ALIGN.CENTER)
    box(s, titre, Inches(2.0), y, Inches(10.7), Inches(0.35),
        font_size=13, bold=True, fg=BLEU)
    box(s, desc, Inches(2.0), y + Inches(0.35), Inches(10.7), Inches(0.75),
        font_size=11, fg=RGBColor(0x33, 0x33, 0x33))

rect(s, Inches(0.4), Inches(6.8), Inches(12.5), Inches(0.04), BLEU_CLAIR)


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 8 — QUESTIONS / DISCUSSION
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, BLEU)

box(s, "Questions & Discussion",
    Inches(1.0), Inches(2.0), Inches(11.3), Inches(1.2),
    font_size=38, bold=True, fg=BLANC, align=PP_ALIGN.CENTER)

box(s, "Stage Pierre Fabre — Modélisation PKPD de l'hématotoxicité — ADC",
    Inches(1.0), Inches(3.4), Inches(11.3), Inches(0.6),
    font_size=16, italic=True, fg=BLEU_CLAIR, align=PP_ALIGN.CENTER)

# Ligne récap
rect(s, Inches(1.5), Inches(4.3), Inches(10.3), Inches(0.04), BLEU_CLAIR)
recap = [
    ("✔ Rat modélisé",    Inches(2.0)),
    ("✔ Humain transféré", Inches(5.8)),
    ("→ ADC en cours",    Inches(9.6)),
]
for txt, x in recap:
    box(s, txt, x, Inches(4.5), Inches(2.6), Inches(0.5),
        font_size=13, bold=True, fg=BLANC, align=PP_ALIGN.CENTER)

rect(s, 0, Inches(6.9), W, Inches(0.6), BLEU_CLAIR)
box(s, "Mars 2026  |  Pipeline : Souris → Rat → Humain  |  Modèle Friberg et al.",
    Inches(0.5), Inches(6.93), Inches(12.3), Inches(0.45),
    font_size=12, fg=BLANC, align=PP_ALIGN.CENTER)


# ─────────────────────────────────────────────────────────────────────────────
out = "/home/user/Hema/reunion_point/reunion_prochaine.pptx"
prs.save(out)
print(f"Sauvegardé : {out}")
