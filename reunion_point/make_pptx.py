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

# Texte gauche
box(s, "Ce qui a été modélisé",
    Inches(0.4), Inches(1.55), Inches(5.5), Inches(0.4),
    font_size=13, bold=True, fg=BLEU)

lignes = [
    ("Progéniteurs",  "MPP → CMP → MEP"),
    ("Myéloïde",      "Neutrophiles  +  Monocytes"),
    ("Érythroïde",    "Réticulocytes  +  GR"),
    ("Thrombocytaire","Plaquettes"),
]
for i, (cat, detail) in enumerate(lignes):
    y = Inches(2.05) + Inches(0.7) * i
    rect(s, Inches(0.4), y, Inches(2.0), Inches(0.5), BLEU_CLAIR)
    box(s, cat, Inches(0.5), y + Inches(0.08), Inches(1.85), Inches(0.35),
        font_size=11, bold=True, fg=BLANC)
    box(s, detail, Inches(2.5), y + Inches(0.08), Inches(3.2), Inches(0.38),
        font_size=11, fg=BLEU)

resultats = [
    "✔  Oscillations cycliques reproduites",
    "✔  125 jours de cinétique validés",
    "✔  PK → Damage → Hématopoïèse ✓",
]
for i, r in enumerate(resultats):
    box(s, r, Inches(0.4), Inches(5.0) + Inches(0.4) * i, Inches(5.5), Inches(0.38),
        font_size=11, fg=RGBColor(0x11, 0x55, 0x11))

# Graphique droite — Figure3_simulation (1800×1650, ~carré)
s.shapes.add_picture("img/Figure3_simulation.png",
                     Inches(5.9), Inches(1.38), Inches(7.1), Inches(5.72))


# ═════════════════════════════════════════════════════════════════════════════
# SLIDE 5 — TRANSFERT RAT → HUMAIN
# ═════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
bg(s, BLANC)
header_bar(s, "Transfert inter-espèce : Rat → Humain",
           "Carboplatin Q21D × 2 | Validation clinique")
footer(s)

box(s, "Carboplatin AUC=5 — Q21D × 2  |  Nadirs neutrophiles & plaquettes",
    Inches(0.4), Inches(1.55), Inches(12.5), Inches(0.4),
    font_size=13, bold=True, fg=BLEU)

# Graphique pleine largeur — Figure4_Q21D_x2 (1500×750, 2:1)
s.shapes.add_picture("img/Figure4_Q21D_x2.png",
                     Inches(0.4), Inches(2.0), Inches(12.5), Inches(4.55))

# Box résultat
rect(s, Inches(0.4), Inches(6.62), Inches(12.5), Inches(0.45), BLEU)
box(s, "Chaîne complète validée :  Souris  →  Rat  →  Humain  ✔",
    Inches(0.7), Inches(6.68), Inches(12.0), Inches(0.35),
    font_size=13, bold=True, fg=BLANC, align=PP_ALIGN.CENTER)


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

# VPC portrait (1050×1350) à gauche
s.shapes.add_picture("img/Figure4_VPC_AUC5.png",
                     Inches(0.3), Inches(1.42), Inches(4.2), Inches(5.65))

# Texte + grades droite
box(s, "1000 patients simulés (GFR=125 fixe)",
    Inches(4.8), Inches(1.6), Inches(8.2), Inches(0.4),
    font_size=12, bold=True, fg=BLEU)

vpc_items = [
    "Intervalles de prédiction 5–95%",
    "Erreur résiduelle intégrée",
    "Nadirs cycliques bien capturés",
]
for i, it in enumerate(vpc_items):
    box(s, "▸  " + it, Inches(4.8), Inches(2.1) + Inches(0.42) * i,
        Inches(8.2), Inches(0.38), font_size=11, fg=BLEU)

rect(s, Inches(4.8), Inches(3.5), Inches(8.2), Inches(0.04), BLEU_CLAIR)
box(s, "Grades NCI-CTCAE (Figure 4c)",
    Inches(4.8), Inches(3.6), Inches(8.2), Inches(0.38),
    font_size=12, bold=True, fg=BLEU)

# Figure4c_grades (1200×750, paysage)
s.shapes.add_picture("img/Figure4c_grades_AUC5.png",
                     Inches(4.8), Inches(4.1), Inches(8.2), Inches(2.75))


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
