"""
PowerPoint — Résultats étape 6 : Composé en développement interne (NHP)
FGFR2 inhibitor | PK/PD semi-mécaniste | Q3W × 3 cycles
"""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
import os

# ── Couleurs ─────────────────────────────────────────────
BLEU_FONCE = RGBColor(0x1A, 0x3A, 0x5C)
BLEU_CLAIR = RGBColor(0x2E, 0x75, 0xB6)
BLEU_PALE  = RGBColor(0xD6, 0xE4, 0xF0)
ORANGE     = RGBColor(0xC5, 0x50, 0x1A)
BLANC      = RGBColor(0xFF, 0xFF, 0xFF)
GRIS_CLAIR = RGBColor(0xF5, 0xF5, 0xF5)
GRIS_MED   = RGBColor(0xCC, 0xCC, 0xCC)
VERT       = RGBColor(0x37, 0x86, 0x44)
ROUGE      = RGBColor(0xC0, 0x00, 0x00)
NOIR       = RGBColor(0x1A, 0x1A, 0x1A)
JAUNE      = RGBColor(0xFF, 0xC0, 0x00)

prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]

# ── Chemins figures ───────────────────────────────────────
FIG_PK = "etape6_fgfr2/results/pk_fit_all_animals.png"
FIG_PD = "etape6_fgfr2/results/poster_PD_predicted_individual_noobs.png"
FIG_PD_OBS = "etape6_fgfr2/results/ind_pred_profiles.png"

def add_slide():
    return prs.slides.add_slide(BLANK)

def rect(slide, x, y, w, h, fill=None, line=None, line_w=Pt(1.5)):
    shp = slide.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.line.width = line_w
    if fill:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    else:
        shp.fill.background()
    if line:
        shp.line.color.rgb = line
    else:
        shp.line.fill.background()
    return shp

def txbox(slide, text, x, y, w, h,
          size=18, bold=False, color=NOIR, align=PP_ALIGN.LEFT,
          italic=False):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return tb

def add_image(slide, path, x, y, w, h):
    if os.path.exists(path):
        slide.shapes.add_picture(path, Inches(x), Inches(y),
                                 Inches(w), Inches(h))
        return True
    else:
        r = rect(slide, x, y, w, h, fill=GRIS_CLAIR, line=GRIS_MED, line_w=Pt(1))
        txbox(slide, f"[Image : {os.path.basename(path)}]",
              x + 0.1, y + h/2 - 0.2, w - 0.2, 0.4,
              size=10, italic=True, color=GRIS_MED, align=PP_ALIGN.CENTER)
        return False

def header_bar(slide, title, subtitle=None):
    rect(slide, 0, 0, 13.33, 1.15, fill=BLEU_FONCE)
    txbox(slide, title, 0.35, 0.08, 12.5, 0.65,
          size=24, bold=True, color=BLANC)
    if subtitle:
        txbox(slide, subtitle, 0.35, 0.72, 12.5, 0.38,
              size=14, color=BLEU_PALE, italic=True)

def footer(slide, txt="Confidentiel — Usage interne"):
    rect(slide, 0, 7.2, 13.33, 0.3, fill=BLEU_FONCE)
    txbox(slide, txt, 0.3, 7.21, 12.5, 0.28,
          size=9, color=GRIS_CLAIR, italic=True)

def bullet(slide, items, x, y, w, size=14, color=NOIR, spacing=0.38):
    for i, item in enumerate(items):
        prefix = "▸ " if not item.startswith("→") else ""
        txbox(slide, prefix + item, x, y + i * spacing, w, spacing + 0.05,
              size=size, color=color)

# ════════════════════════════════════════════════════════
# SLIDE 1 — Titre
# ════════════════════════════════════════════════════════
s1 = add_slide()
rect(s1, 0, 0, 13.33, 7.5, fill=BLEU_FONCE)
rect(s1, 0, 2.6, 13.33, 2.8, fill=BLEU_CLAIR)

txbox(s1, "Étape 6 — Composé en développement interne",
      1.0, 1.3, 11.33, 0.7, size=18, italic=True, color=BLEU_PALE, align=PP_ALIGN.CENTER)

txbox(s1, "Modélisation PK/PD de l'hématotoxicité",
      0.5, 2.8, 12.33, 1.0, size=36, bold=True, color=BLANC, align=PP_ALIGN.CENTER)

txbox(s1, "Données NHP précliniques — 8 singes cynomolgus — 4 niveaux de dose",
      1.0, 3.85, 11.33, 0.6, size=18, color=BLANC, align=PP_ALIGN.CENTER)

txbox(s1, "Modèle semi-mécaniste Fornari (2019)  |  Ajustement PK individuel  |  Calibration PD visuelle",
      1.0, 4.5, 11.33, 0.5, size=14, italic=True, color=BLEU_PALE, align=PP_ALIGN.CENTER)

rect(s1, 3.5, 5.3, 6.33, 0.06, fill=JAUNE)
txbox(s1, "Résultats confidentiels — ne pas diffuser",
      1.0, 5.5, 11.33, 0.45, size=13, italic=True, color=JAUNE, align=PP_ALIGN.CENTER)

txbox(s1, "M2 Sciences de la Donnée de Santé  |  Stage Pierre Fabre  |  2024–2025",
      1.0, 6.5, 11.33, 0.45, size=12, color=GRIS_CLAIR, align=PP_ALIGN.CENTER)

# ════════════════════════════════════════════════════════
# SLIDE 2 — Design expérimental
# ════════════════════════════════════════════════════════
s2 = add_slide()
header_bar(s2, "Design expérimental",
           "Étude PK/PD préclinique chez le singe cynomolgus (NHP)")
footer(s2)

# Tableau design
col_x = [0.4, 3.2, 6.0, 8.8, 11.0]
col_headers = ["Groupe", "Dose", "N animaux", "Schéma", "Observations PK/PD"]
col_w = [2.7, 2.7, 2.7, 2.2, 2.2]

# En-têtes tableau
for i, h in enumerate(col_headers):
    rect(s2, col_x[i], 1.3, col_w[i] - 0.05, 0.45, fill=BLEU_CLAIR)
    txbox(s2, h, col_x[i] + 0.05, 1.32, col_w[i] - 0.1, 0.4,
          size=13, bold=True, color=BLANC, align=PP_ALIGN.CENTER)

rows = [
    ["Groupe 1", "4 mg/kg",  "2 (1001, 1002)", "Q3W × 3", "PK multi-temps + hémogramme"],
    ["Groupe 2", "13 mg/kg", "2 (2001, 2002)", "Q3W × 3", "PK multi-temps + hémogramme"],
    ["Groupe 3", "26 mg/kg", "2 (4001, 4002)", "Q3W × 3", "PK multi-temps + hémogramme"],
    ["Groupe 4", "39 mg/kg", "2 (3101, 3002)", "Q3W × 3", "PK multi-temps + hémogramme"],
]
row_colors = [BLANC, GRIS_CLAIR, BLANC, GRIS_CLAIR]
for ri, row in enumerate(rows):
    y = 1.75 + ri * 0.48
    for ci, cell in enumerate(row):
        rect(s2, col_x[ci], y, col_w[ci] - 0.05, 0.45,
             fill=row_colors[ri], line=GRIS_MED, line_w=Pt(0.5))
        txbox(s2, cell, col_x[ci] + 0.08, y + 0.05, col_w[ci] - 0.15, 0.38,
              size=12, color=NOIR, align=PP_ALIGN.CENTER)

# Points clés
txbox(s2, "Points clés du pipeline", 0.4, 4.2, 6.0, 0.4,
      size=15, bold=True, color=BLEU_FONCE)
bullet(s2, [
    "Ajustement PK individuel : modèle 2 compartiments (CL, V1, Q, V2) par Nelder-Mead",
    "Concentrations prédites utilisées comme entrée du modèle PD",
    "Paramètres PD (Slope_MPP, Slope_CMP, Slope_MEP) ajustés visuellement / animal",
    "Baselines hématologiques individuelles dérivées des valeurs pré-dose (J-3)",
    "Cellules suivies : Neutrophiles, Plaquettes, Réticulocytes, RBC",
], 0.4, 4.65, 6.8, size=13)

# Schéma temporel
rect(s2, 7.3, 4.2, 5.6, 2.9, fill=BLEU_PALE, line=BLEU_CLAIR, line_w=Pt(1))
txbox(s2, "Schéma posologique", 7.5, 4.25, 5.0, 0.4,
      size=14, bold=True, color=BLEU_FONCE)
bullet(s2, [
    "Dose 1 : Jour 0",
    "Dose 2 : Jour 21",
    "Dose 3 : Jour 42",
    "Suivi PD jusqu'à J25 (post-dose 1)",
    "Récupération observée après J42",
], 7.5, 4.7, 5.0, size=13, spacing=0.35)

# ════════════════════════════════════════════════════════
# SLIDE 3 — Résultats PK
# ════════════════════════════════════════════════════════
s3 = add_slide()
header_bar(s3, "Résultats PK — Ajustement individuel 2 compartiments",
           "Modèle bi-exponentiel ajusté par Nelder-Mead (SSR log-scale) pour chaque animal")
footer(s3)

# Figure PK (à gauche)
add_image(s3, FIG_PK, 0.3, 1.25, 8.5, 5.9)

# Commentaires (à droite)
rect(s3, 9.0, 1.25, 4.0, 5.9, fill=GRIS_CLAIR, line=GRIS_MED, line_w=Pt(0.5))
txbox(s3, "Observations", 9.1, 1.35, 3.8, 0.45,
      size=15, bold=True, color=BLEU_FONCE)
bullet(s3, [
    "Excellent ajustement sur\nl'échelle log pour les 4 doses",
    "Cinétique bi-exponentielle\nbien capturée (phases α et β)",
    "Variabilité inter-animale\nmodérée sur CL et V1",
    "Dose-proportionnalité PK\nconfirmée (Cmax et AUC)",
    "Demi-vie terminale t½β\ncohérente entre animaux",
], 9.1, 1.85, 3.75, size=12, spacing=0.9)

txbox(s3, "Modèle :", 9.1, 6.35, 3.8, 0.3,
      size=12, bold=True, color=BLEU_FONCE)
txbox(s3, "dA1/dt = −(CL/V1 + Q/V1)·A1 + (Q/V2)·A2\nC(t) = A·e⁻ᵅᵗ + B·e⁻ᵝᵗ",
      9.1, 6.65, 3.8, 0.55, size=11, italic=True, color=NOIR)

# ════════════════════════════════════════════════════════
# SLIDE 4 — Résultats PD (avec observations)
# ════════════════════════════════════════════════════════
s4 = add_slide()
header_bar(s4, "Résultats PD — Profils individuels : modèle vs observations",
           "% de la baseline individuelle (J-3)  |  Lignes = modèle  |  Points = données observées")
footer(s4)

add_image(s4, FIG_PD_OBS, 0.3, 1.2, 9.5, 6.0)

rect(s4, 9.9, 1.2, 3.1, 6.0, fill=GRIS_CLAIR, line=GRIS_MED, line_w=Pt(0.5))
txbox(s4, "Lectures", 10.0, 1.3, 2.9, 0.4,
      size=15, bold=True, color=BLEU_FONCE)
bullet(s4, [
    "✓ Neutrophiles :\nbonne concordance\nmodèle–données",
    "✓ RBC & Plaquettes :\nsuppression modérée\nbien reproduite",
    "⚠ Réticulocytes :\npics D2–D8 non\ncaptés par le modèle",
    "→ Ces pics résultent\nde facteurs extérieurs\n(état préexistant,\nstress de manipulation)",
], 10.0, 1.8, 2.9, size=11, spacing=1.1)

# ════════════════════════════════════════════════════════
# SLIDE 5 — Résultats PD (prédictions seules, D0–D90)
# ════════════════════════════════════════════════════════
s5 = add_slide()
header_bar(s5, "Projections PD — Simulation Q3W × 3 cycles + récupération",
           "Profils prédits D0–D90 | % de la baseline individuelle | Sans données observées")
footer(s5)

add_image(s5, FIG_PD, 0.3, 1.2, 9.5, 6.0)

rect(s5, 9.9, 1.2, 3.1, 6.0, fill=GRIS_CLAIR, line=GRIS_MED, line_w=Pt(0.5))
txbox(s5, "Points clés", 10.0, 1.3, 2.9, 0.4,
      size=15, bold=True, color=BLEU_FONCE)
bullet(s5, [
    "Effet cumulatif sur\n3 cycles visible",
    "Nadir le plus profond\naprès la 3ᵉ dose (J42)",
    "Récupération partielle\nentre cycles",
    "Dose-réponse claire :\n39 mg/kg > 26 > 13 > 4",
    "Retour à la baseline\nattendu ~J80–D90",
], 10.0, 1.8, 2.9, size=12, spacing=0.9)

txbox(s5, "NOAEL", 10.0, 6.1, 2.9, 0.3,
      size=13, bold=True, color=VERT)
txbox(s5, "Aucun grade ≥3 prédit\nà 4 mg/kg", 10.0, 6.4, 2.9, 0.55,
      size=12, color=VERT)

# ════════════════════════════════════════════════════════
# SLIDE 6 — Prochaines étapes
# ════════════════════════════════════════════════════════
s6 = add_slide()
header_bar(s6, "Prochaines étapes",
           "De la caractérisation NHP à la prédiction clinique (FIH)")
footer(s6)

boxes = [
    ("1. Estimation formelle des paramètres PD",
     BLEU_CLAIR,
     ["Remplacer l'ajustement visuel par une optimisation numérique (Nelder-Mead ou NLME)",
      "Quantifier l'incertitude sur Slope_MPP, Slope_CMP, Slope_MEP",
      "Calculer des intervalles de confiance sur les nadirs prédits"]),
    ("2. Transposition PK NHP → Humain",
     VERT,
     ["Allométrie standard : CL_hum = CL_NHP × (70/BW_NHP)⁰·⁷⁵",
      "Validée sur T-DXd (facteur erreur < 2 sur CL)",
      "Générer les paramètres PK humains de départ"]),
    ("3. Simulation de la toxicité humaine attendue",
     ORANGE,
     ["Utiliser les Slopes NHP comme estimation initiale chez l'humain",
      "Simuler N=300 patients virtuels (IIV log-normal)",
      "Prédire la distribution des grades CTCAE v5 à différentes doses FIH"]),
    ("4. Support à la dose de départ FIH",
     RGBColor(0x60, 0x32, 0x8A),
     ["Identifier la dose associée à P(G≥3) < 10%",
      "Intégrer le NOAEL NHP comme contrainte de sécurité",
      "Documenter pour le dossier IND/CTA (FDA/EMA)"]),
]

for i, (title, color, items) in enumerate(boxes):
    col = i % 2
    row = i // 2
    x = 0.3 + col * 6.5
    y = 1.35 + row * 2.85

    rect(s6, x, y, 6.2, 2.65, fill=BLANC, line=color, line_w=Pt(2.5))
    rect(s6, x, y, 6.2, 0.45, fill=color)
    txbox(s6, title, x + 0.1, y + 0.04, 6.0, 0.38,
          size=13, bold=True, color=BLANC)
    for j, item in enumerate(items):
        txbox(s6, f"▸ {item}", x + 0.15, y + 0.52 + j * 0.62, 5.9, 0.58,
              size=11.5, color=NOIR)

# ════════════════════════════════════════════════════════
# SLIDE 7 — Synthèse
# ════════════════════════════════════════════════════════
s7 = add_slide()
header_bar(s7, "Synthèse — Ce que le modèle apporte",
           "Valeur ajoutée de la modélisation PK/PD pour le développement interne")
footer(s7)

cols = [
    ("Réalisé ✓", VERT, [
        "Ajustement PK 2-comp individuel\n(8 animaux, 4 doses)",
        "Calibration PD multi-lignées\n(Neut, Plt, Ret, RBC)",
        "Simulation 3 cycles + récupération\n(D0–D90)",
        "Prédiction nadir dose-dépendant\nconfirmée",
    ]),
    ("Limite connue ⚠", ORANGE, [
        "Ajustement PD visuel (non formel) :\nincertitude non quantifiée",
        "n=2 animaux/groupe : variabilité\ninter-individuelle mal estimée",
        "Pics réticulocytes non modélisés\n(facteurs extérieurs)",
        "Transposition clinique\nencore préliminaire",
    ]),
    ("Perspectives →", BLEU_CLAIR, [
        "Estimation formelle NLME\n(nlmixr2) sur données NHP",
        "Allométrie NHP → humain\n+ simulation FIH",
        "IC50 CFU cellules humaines\npour ancrer les Slopes",
        "Dossier réglementaire\nIND/CTA support M&S",
    ]),
]

for i, (title, color, items) in enumerate(cols):
    x = 0.3 + i * 4.3
    rect(s7, x, 1.25, 4.1, 5.9, fill=BLANC, line=color, line_w=Pt(2))
    rect(s7, x, 1.25, 4.1, 0.5, fill=color)
    txbox(s7, title, x + 0.1, 1.28, 3.9, 0.42,
          size=15, bold=True, color=BLANC, align=PP_ALIGN.CENTER)
    for j, item in enumerate(items):
        rect(s7, x + 0.15, 1.85 + j * 1.18, 3.8, 1.1,
             fill=GRIS_CLAIR, line=color, line_w=Pt(0.5))
        txbox(s7, item, x + 0.25, 1.9 + j * 1.18, 3.6, 1.0,
              size=12, color=NOIR)

# ── Export ────────────────────────────────────────────────
out = "FGFR2_NHP_resultats.pptx"
prs.save(out)
print(f"✓ PowerPoint généré : {out}")
print(f"  7 slides : Titre | Design | PK | PD obs | PD pred | Prochaines étapes | Synthèse")
print(f"\n  Note : insérer les figures manuellement si les PNG ne sont pas générés :")
print(f"    - Slide 3 : {FIG_PK}")
print(f"    - Slide 4 : {FIG_PD_OBS}")
print(f"    - Slide 5 : {FIG_PD}")
