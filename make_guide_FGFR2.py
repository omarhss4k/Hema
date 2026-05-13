"""
Génère guide_utilisation_FGFR2.pptx
Document destiné à la tutrice (non-codeur R) :
  - Quels modèles ont été utilisés et pourquoi
  - Quelles hypothèses ont été faites
  - Comment le modèle fonctionne biologiquement
  - Comment lancer les scripts après le stage
  - Comment lire les résultats
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
import os

OUT = "scripts/presentation/guide_utilisation_FGFR2.pptx"

# ── Palette ──────────────────────────────────────────────────────────────────
NAVY    = RGBColor(0x1A, 0x52, 0x76)
BLUE    = RGBColor(0x2E, 0x86, 0xC1)
WHITE   = RGBColor(0xFF, 0xFF, 0xFF)
GREY    = RGBColor(0x55, 0x55, 0x55)
GREEN   = RGBColor(0x1E, 0x84, 0x49)
RED     = RGBColor(0xC0, 0x39, 0x2B)
ORANGE  = RGBColor(0xE6, 0x7E, 0x22)
LGREY   = RGBColor(0xF4, 0xF6, 0xF7)
LBLUE   = RGBColor(0xD6, 0xEA, 0xF8)
LGREEN  = RGBColor(0xD5, 0xF5, 0xE3)
LYELLOW = RGBColor(0xFE, 0xF9, 0xE7)

prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)
blank = prs.slide_layouts[6]

# ── Helpers ───────────────────────────────────────────────────────────────────

def rect(slide, l, t, w, h, fill, border=None):
    s = slide.shapes.add_shape(1, Inches(l), Inches(t), Inches(w), Inches(h))
    s.fill.solid(); s.fill.fore_color.rgb = fill
    if border:
        s.line.color.rgb = border; s.line.width = Pt(1)
    else:
        s.line.fill.background()
    return s

def txt(slide, text, l, t, w, h,
        sz=13, bold=False, col=GREY, align=PP_ALIGN.LEFT,
        italic=False, wrap=True):
    tb = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = wrap
    p  = tf.paragraphs[0]; p.alignment = align
    r  = p.add_run(); r.text = text
    r.font.size = Pt(sz); r.font.bold = bold
    r.font.color.rgb = col; r.font.italic = italic
    return tb

def header(slide, title, subtitle="", color=NAVY):
    rect(slide, 0, 0, 13.33, 1.1, color)
    txt(slide, title, 0.25, 0.07, 12.5, 0.72,
        sz=24, bold=True, col=WHITE)
    if subtitle:
        txt(slide, subtitle, 0.25, 0.78, 12.5, 0.32,
            sz=11, col=RGBColor(0xAE,0xD6,0xF1))

def badge(slide, label, l, t, fill=LBLUE, border=BLUE, sz=11, bold=False, col=NAVY):
    rect(slide, l, t, len(label)*0.085+0.2, 0.35, fill, border)
    txt(slide, label, l+0.1, t+0.04, len(label)*0.085, 0.28,
        sz=sz, bold=bold, col=col)

def bullet_box(slide, title, items, l, t, w, h,
               fill=LGREY, title_col=NAVY, item_col=GREY):
    rect(slide, l, t, w, h, fill, BLUE)
    txt(slide, title, l+0.15, t+0.1, w-0.3, 0.38,
        sz=12, bold=True, col=title_col)
    step = (h - 0.5) / max(len(items), 1)
    for i, item in enumerate(items):
        txt(slide, item, l+0.2, t+0.48+i*step, w-0.35, step+0.05,
            sz=11, col=item_col)

def numbered_step(slide, num, text, l, t, w,
                  num_col=WHITE, num_bg=NAVY, txt_col=GREY):
    rect(slide, l, t, 0.42, 0.42, num_bg)
    txt(slide, str(num), l+0.01, t+0.02, 0.4, 0.38,
        sz=14, bold=True, col=num_col, align=PP_ALIGN.CENTER)
    txt(slide, text, l+0.5, t+0.03, w-0.55, 0.45, sz=11, col=txt_col)

def footer_note(slide, note):
    txt(slide, note, 0.2, 7.1, 12.9, 0.35,
        sz=9, col=RGBColor(0x99,0x99,0x99), italic=True)


# =============================================================================
# DIAPO 0 — Couverture
# =============================================================================
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, 13.33, 7.5, NAVY)
rect(sl, 0, 2.8, 13.33, 2.0, RGBColor(0x15,0x40,0x60))

txt(sl, "Guide d'utilisation", 1.0, 1.0, 11.5, 1.0,
    sz=36, bold=True, col=WHITE, align=PP_ALIGN.CENTER)
txt(sl, "Modélisation PK/PD — Fc-silent FGFR2-huBPA-LP1",
    1.0, 2.0, 11.5, 0.7,
    sz=20, col=RGBColor(0xAE,0xD6,0xF1), align=PP_ALIGN.CENTER)
txt(sl, "Ce document explique les modeles utilises, les hypotheses choisies,\n"
        "et comment relancer les analyses apres le stage.",
    1.5, 2.95, 10.5, 0.9,
    sz=14, col=RGBColor(0xD6,0xEA,0xF8), align=PP_ALIGN.CENTER)
txt(sl,
    "Modele PK : 2 compartiments   |   Modele PD : Simeoni-Emax   |   Langage : R",
    1.0, 4.2, 11.5, 0.5,
    sz=12, col=RGBColor(0x85,0xC1,0xE9), align=PP_ALIGN.CENTER)

rect(sl, 2.5, 5.1, 8.5, 0.06, RGBColor(0x2E,0x86,0xC1))
txt(sl, "Redige par : [votre nom]   |   Stage 2026   |   Encadrante : [nom tutrice]",
    1.0, 5.3, 11.5, 0.5,
    sz=11, col=RGBColor(0x85,0xC1,0xE9), align=PP_ALIGN.CENTER)


# =============================================================================
# DIAPO 1 — Table des matières
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "Sommaire du guide", "Ce document repond a 5 questions essentielles")

sections = [
    ("1", "Quelles donnees sont necessaires ?",
     "Format des fichiers Excel, unites, structure attendue", LBLUE),
    ("2", "Quel modele PK a ete utilise et pourquoi ?",
     "Modele 2 compartiments, hypotheses, parametres estimes", LBLUE),
    ("3", "Quel modele PD a ete utilise et pourquoi ?",
     "Simeoni 2004 + Emax, hypotheses biologiques, parametres", LGREEN),
    ("4", "Comment lancer les analyses ?",
     "Ordre des scripts, fichiers produits, que faire si erreur", LYELLOW),
    ("5", "Comment lire et interpreter les resultats ?",
     "Graphiques, tableaux, TGI, validation du modele", RGBColor(0xFD,0xED,0xEC)),
]

for i, (num, title, sub, fill) in enumerate(sections):
    y = 1.3 + i * 1.12
    rect(sl, 0.4, y, 12.5, 0.95, fill, BLUE)
    txt(sl, num, 0.5, y+0.12, 0.45, 0.7,
        sz=22, bold=True, col=NAVY, align=PP_ALIGN.CENTER)
    txt(sl, title, 1.1, y+0.08, 9.5, 0.42, sz=13, bold=True, col=NAVY)
    txt(sl, sub,   1.1, y+0.52, 9.5, 0.35, sz=10, col=GREY)


# =============================================================================
# DIAPO 2 — Données nécessaires
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "1 — Donnees necessaires", "Deux fichiers Excel a fournir avant de lancer les scripts")

# Fichier 1
rect(sl, 0.3, 1.2, 5.9, 5.5, LBLUE, BLUE)
txt(sl, "Fichier 1 : PK souris FGFR2.xlsx", 0.45, 1.28, 5.6, 0.45,
    sz=13, bold=True, col=NAVY)
txt(sl,
    "Contenu :\n"
    "  Colonne 1 : temps en heures\n"
    "  Colonne 2 : concentrations 10 mg/kg (µg/L)\n"
    "  Colonne 3 : concentrations  5 mg/kg (µg/L)\n"
    "  Colonne 4 : concentrations  1 mg/kg (µg/L)\n\n"
    "Valeurs manquantes : laisser la cellule vide\n"
    "(le script les ignore automatiquement)\n\n"
    "Unites importantes :\n"
    "  Temps    → heures\n"
    "  Doses    → mg/kg (converties en µg/kg par le script)\n"
    "  Conc.    → µg/L\n\n"
    "Nombre de points : au moins 6 par dose\n"
    "recommande pour un bon ajustement",
    0.45, 1.78, 5.6, 4.7, sz=11, col=GREY)

# Fichier 2
rect(sl, 6.8, 1.2, 6.2, 5.5, LGREEN, GREEN)
txt(sl, "Fichier 2 : TumorVolume_FGFR2.xlsx", 6.95, 1.28, 5.9, 0.45,
    sz=13, bold=True, col=RGBColor(0x1E,0x84,0x49))
txt(sl,
    "Contenu : 2 blocs separes par une ligne vide\n\n"
    "  Bloc 1 - Moyennes (mm3) :\n"
    "    Ligne 1 : 'Group' | j0 | j4 | j7 | ...\n"
    "    Ligne 2 : Group 01 (controle) | valeurs\n"
    "    Ligne 3 : Group 03 (10 mg/kg) | valeurs\n"
    "    Ligne 4 : Group 04 ( 3 mg/kg) | valeurs\n\n"
    "  Bloc 2 - SEM (mm3) :\n"
    "    Meme structure que le bloc 1\n\n"
    "Unites : volume tumoral en mm3\n"
    "(converti en grammes par le script : /1000)\n\n"
    "Protocole encode dans le script :\n"
    "  Q2W × 4 doses (j0, j14, j28, j42)\n"
    "  A adapter si le protocole change",
    6.95, 1.78, 5.9, 4.7, sz=11, col=GREY)

footer_note(sl, "Si les noms de groupes changent (ex: Group 02 au lieu de Group 01), "
                "modifier les patterns dans pkpd_FGFR2.R section 12")


# =============================================================================
# DIAPO 3 — Modèle PK : description
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "2 — Modele PK : 2 compartiments IV bolus",
       "Choix justifie par la cinetique bi-exponentielle observee dans les donnees")

# Schéma simplifié
rect(sl, 0.3, 1.2, 3.5, 1.4, LBLUE, BLUE)
txt(sl, "Dose injectee\n(IV bolus)", 0.4, 1.3, 3.2, 1.1,
    sz=12, bold=True, col=NAVY, align=PP_ALIGN.CENTER)

txt(sl, "→", 3.8, 1.7, 0.5, 0.5, sz=22, bold=True, col=NAVY)

rect(sl, 4.3, 1.2, 4.0, 1.4, LBLUE, BLUE)
txt(sl, "Compartiment central\n(sang + organes)\nC1 = dose / V1",
    4.4, 1.3, 3.7, 1.1, sz=11, col=NAVY, align=PP_ALIGN.CENTER)

txt(sl, "→ elimination", 8.3, 1.7, 1.8, 0.5, sz=11, col=RED)
rect(sl, 10.1, 1.45, 2.9, 0.9, RGBColor(0xFA,0xDB,0xD8), RED)
txt(sl, "CL × C1\n(urine, foie...)", 10.2, 1.5, 2.7, 0.8,
    sz=10, col=RED, align=PP_ALIGN.CENTER)

txt(sl, "↕ echange Q", 5.8, 2.6, 2.0, 0.4, sz=11, col=GREEN)

rect(sl, 4.3, 3.0, 4.0, 1.2, LGREEN, GREEN)
txt(sl, "Compartiment peripherique\n(tissus, muscles...)",
    4.4, 3.1, 3.7, 0.9, sz=11, col=GREEN, align=PP_ALIGN.CENTER)

# Hypothèses
rect(sl, 0.3, 4.4, 6.0, 2.8, LYELLOW, ORANGE)
txt(sl, "Hypotheses du modele PK", 0.45, 4.5, 5.7, 0.4,
    sz=13, bold=True, col=ORANGE)
hyp_pk = [
    "Injection instantanee dans le sang (IV bolus)",
    "Le medicament se repartit entre 2 zones : sang et tissus",
    "L'elimination ne se fait que depuis le compartiment central",
    "Les transferts sang<->tissus sont proportionnels aux concentrations",
    "La pharmacocinetique est LINEAIRE : doubler la dose double la concentration",
    "Les parametres sont les memes pour toutes les souris (pas de variabilite)",
]
for i, h in enumerate(hyp_pk):
    txt(sl, "✓  " + h, 0.45, 4.95+i*0.37, 5.7, 0.35, sz=10, col=GREY)

# Paramètres
rect(sl, 6.5, 4.4, 6.5, 2.8, LGREY, BLUE)
txt(sl, "Les 4 parametres estimes automatiquement",
    6.65, 4.5, 6.2, 0.4, sz=13, bold=True, col=NAVY)
params_pk = [
    ("CL",  "Clairance : vitesse d'elimination depuis le sang"),
    ("V1",  "Volume central : espace de distribution dans le sang"),
    ("V2",  "Volume peripherique : espace dans les tissus"),
    ("Q",   "Clairance intercompartimentale : echange sang<->tissus"),
]
for i, (p, d) in enumerate(params_pk):
    txt(sl, p, 6.65, 4.98+i*0.55, 0.7, 0.45, sz=13, bold=True, col=BLUE)
    txt(sl, d, 7.4,  4.98+i*0.55, 5.4, 0.45, sz=10, col=GREY)

footer_note(sl, "La linearite PK a ete verifiee : les courbes a 1, 5 et 10 mg/kg sont "
                "proportionnelles — confirme la validite du modele 2-compartiments")


# =============================================================================
# DIAPO 4 — Modèle PD : Simeoni
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "3 — Modele PD : Simeoni 2004 (croissance + transit)",
       "Modele publie, valide sur >50 molecules anticancereuses en preclinique")

# Explication biologique
rect(sl, 0.3, 1.2, 7.8, 2.5, LGREEN, GREEN)
txt(sl, "Ce que le modele represente biologiquement",
    0.45, 1.28, 7.5, 0.4, sz=13, bold=True, col=GREEN)
txt(sl,
    "Une tumeur contient des cellules dans differents etats :\n\n"
    "  x1 = cellules qui proliferent activement  (sensibles au medicament)\n"
    "  x2 = cellules endommagees par le medicament (ne proliferent plus)\n"
    "  x3 = cellules en cours de mort (transit)\n"
    "  x4 = cellules mourantes (bientot eliminées)\n\n"
    "Volume tumoral total = x1 + x2 + x3 + x4",
    0.45, 1.72, 7.5, 1.85, sz=11, col=GREY)

# Hypothèses PD
rect(sl, 0.3, 3.85, 7.8, 3.3, LYELLOW, ORANGE)
txt(sl, "Hypotheses du modele PD", 0.45, 3.93, 7.5, 0.4,
    sz=13, bold=True, col=ORANGE)
hyp_pd = [
    "Croissance exponentielle au debut (l0), lineaire quand la tumeur est grande (l1)",
    "La tumeur ne peut pas grossir indefiniment : elle est limitee par sa propre taille",
    "Seules les cellules x1 (proliferantes) sont tuees par le medicament",
    "Les cellules endommagees meurent progressivement (pas instantanement)",
    "Le temps moyen pour mourir apres dommage = 4/k1 (ici ~8 jours, Simeoni 2004)",
    "Les parametres de croissance (l0, l1) sont calibres sur le groupe controle seul",
    "Le coefficient de Hill p=1 (modele original Simeoni 2004)",
]
for i, h in enumerate(hyp_pd):
    txt(sl, "✓  " + h, 0.45, 4.38+i*0.40, 7.5, 0.38, sz=10, col=GREY)

# Paramètres PD
rect(sl, 8.3, 1.2, 4.7, 5.95, LGREY, BLUE)
txt(sl, "Parametres PD", 8.45, 1.28, 4.4, 0.4, sz=13, bold=True, col=NAVY)
params_pd = [
    ("l0",   "Taux de croissance\nexponentielle (/h)"),
    ("l1",   "Taux de croissance\nlineaire (g/h)"),
    ("p",    "Coeff. de Hill\n(fixe a 1)"),
    ("k1",   "Vitesse de transit\n(fixee, Simeoni 2004)"),
    ("Emax", "Destruction max\npar le medicament"),
    ("EC50", "Conc. a 50%\nde l'effet max"),
]
for i, (p, d) in enumerate(params_pd):
    bg = LGREEN if p in ("Emax","EC50") else LGREY
    rect(sl, 8.35, 1.72+i*0.75, 4.55, 0.68, bg)
    txt(sl, p, 8.45, 1.78+i*0.75, 0.8, 0.55,
        sz=14, bold=True, col=NAVY if p not in ("Emax","EC50") else GREEN)
    txt(sl, d, 9.35, 1.78+i*0.75, 3.4, 0.55, sz=10, col=GREY)

txt(sl, "Vert = estimes sur vos donnees\nBlanc = fixes (litterature)",
    8.35, 6.35, 4.55, 0.6, sz=9, col=GREY, italic=True)

footer_note(sl, "Reference : Simeoni M. et al., Cancer Research 2004, 64:1094-1101 "
                "| Modele le plus cite en oncologie preclinique")


# =============================================================================
# DIAPO 5 — Pourquoi Emax
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "3 (suite) — Pourquoi le modele Emax et pas un modele plus simple ?",
       "Le choix est justifie par les donnees experimentales elles-memes")

# Problème du modèle linéaire
rect(sl, 0.3, 1.2, 6.0, 2.7, RGBColor(0xFD,0xED,0xEC), RED)
txt(sl, "Modele lineaire (abandonne) : effet = k2 x C1",
    0.45, 1.28, 5.7, 0.4, sz=13, bold=True, col=RED)
txt(sl,
    "Ce modele suppose que doubler la dose\n"
    "double toujours l'effet.\n\n"
    "Probleme observe :\n"
    "  - 3 mg/kg  → TGI observe = 64%\n"
    "  - 10 mg/kg → TGI observe = 84%\n"
    "  Ratio reel : 84/64 = 1.3x\n\n"
    "  Si l'effet etait lineaire :\n"
    "  10/3 = 3.3x plus d'effet attendu\n\n"
    "→ Le modele convergeait vers k2=0\n"
    "   (aucun effet predit)",
    0.45, 1.72, 5.7, 2.1, sz=11, col=GREY)

# Emax : solution
rect(sl, 6.5, 1.2, 6.5, 2.7, LGREEN, GREEN)
txt(sl, "Modele Emax (choisi) : effet = Emax x C1 / (EC50 + C1)",
    6.65, 1.28, 6.2, 0.4, sz=12, bold=True, col=GREEN)
txt(sl,
    "Ce modele inclut un plafond d'effet :\n"
    "il y a un nombre limite de recepteurs\n"
    "FGFR2 a saturer sur les cellules.\n\n"
    "  - A 3 mg/kg  : effets deja forts\n"
    "                (recepteurs satures a ~X%)\n"
    "  - A 10 mg/kg : effet marginalement\n"
    "                plus fort (saturation totale)\n\n"
    "→ Permet de reproduire le ratio 1.3x\n"
    "   observe experimentalement\n\n"
    "→ Biologiquement coherent avec\n"
    "   le mecanisme d'action des Ac",
    6.65, 1.72, 6.2, 2.1, sz=11, col=GREY)

# Les 2 paramètres Emax expliqués simplement
rect(sl, 0.3, 4.05, 12.7, 3.1, LGREY, BLUE)
txt(sl, "Les 2 parametres Emax expliques simplement",
    0.45, 4.13, 12.4, 0.4, sz=13, bold=True, col=NAVY)

rect(sl, 0.4, 4.6, 5.9, 2.4, LGREEN, GREEN)
txt(sl, "Emax = taux de mort cellulaire maximal",
    0.55, 4.68, 5.6, 0.4, sz=12, bold=True, col=GREEN)
txt(sl,
    "Meme si on injecte 100 mg/kg, le medicament\n"
    "ne peut pas tuer les cellules plus vite que Emax.\n\n"
    "Pour que la tumeur regresse : Emax > l0\n"
    "(destruction > croissance)",
    0.55, 5.12, 5.6, 1.8, sz=11, col=GREY)

rect(sl, 6.7, 4.6, 6.2, 2.4, LBLUE, BLUE)
txt(sl, "EC50 = concentration a 50% de l'effet max",
    6.85, 4.68, 5.9, 0.4, sz=12, bold=True, col=NAVY)
txt(sl,
    "C'est la 'puissance' du medicament.\n\n"
    "  EC50 petit → medicament tres puissant\n"
    "               (effet fort meme a faible dose)\n\n"
    "  EC50 grand → medicament peu puissant\n"
    "               (il faut de fortes concentrations)",
    6.85, 5.12, 5.9, 1.8, sz=11, col=GREY)

footer_note(sl,
    "Emax et EC50 sont estimes automatiquement par un algorithme d'optimisation globale (DEoptim)")


# =============================================================================
# DIAPO 6 — Organisation des scripts
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "4 — Organisation des fichiers et scripts",
       "Tout est dans le depot GitHub — aucune installation manuelle requise (hors packages R)")

rect(sl, 0.3, 1.2, 12.7, 1.0, LGREY, BLUE)
txt(sl, "Structure des dossiers", 0.45, 1.28, 12.4, 0.35,
    sz=12, bold=True, col=NAVY)
txt(sl,
    "  FGFR/\n"
    "  ├── PK souris FGFR2.xlsx          ← donnees PK (a fournir)\n"
    "  ├── TumorVolume_FGFR2.xlsx        ← donnees tumorales (a fournir)\n"
    "  └── scripts/\n"
    "      ├── pkpd_FGFR2.R              ← script principal (PK + PD + Emax)\n"
    "      ├── presentation_FGFR2.R      ← script figures de presentation\n"
    "      └── presentation/             ← figures PNG + PPTX generes",
    0.45, 1.58, 12.4, 0.55, sz=10, col=GREY)

# Étapes
rect(sl, 0.3, 2.35, 12.7, 4.8, RGBColor(0xF8,0xF9,0xFA), BLUE)
txt(sl, "Comment lancer l'analyse — dans RStudio",
    0.45, 2.43, 12.4, 0.4, sz=13, bold=True, col=NAVY)

steps = [
    ("Installer les packages R (une seule fois)",
     'install.packages(c("rxode2","deSolve","DEoptim","ggplot2","readxl"))'),
    ("Ouvrir RStudio et definir le dossier de travail",
     'setwd("~/FGFR")   # chemin vers le dossier contenant les fichiers Excel'),
    ("Lancer le script principal (PK + calibration + Emax)",
     'source("scripts/pkpd_FGFR2.R")   # duree : ~5-15 min (optimisation)'),
    ("Lancer le script de presentation (9 figures + PowerPoint)",
     'source("scripts/presentation_FGFR2.R")   # duree : ~1-2 min'),
    ("Recuperer les resultats",
     "scripts/presentation/presentation_FGFR2.pptx   ← ouvrir ce fichier"),
]
for i, (desc, code) in enumerate(steps):
    y = 2.9 + i*0.82
    rect(sl, 0.4, y, 0.5, 0.45, NAVY)
    txt(sl, str(i+1), 0.41, y+0.03, 0.48, 0.38,
        sz=14, bold=True, col=WHITE, align=PP_ALIGN.CENTER)
    txt(sl, desc, 1.0, y+0.02, 4.8, 0.35, sz=11, bold=True, col=NAVY)
    rect(sl, 5.9, y, 7.2, 0.42, RGBColor(0x2C,0x3E,0x50))
    txt(sl, code, 6.0, y+0.06, 7.0, 0.32, sz=9,
        col=RGBColor(0xAE,0xD6,0xF1))

footer_note(sl,
    "Si un message d'erreur apparait : verifier que les fichiers Excel sont dans le bon dossier "
    "et que les noms de colonnes/groupes correspondent a ce qui est attendu")


# =============================================================================
# DIAPO 7 — Fichiers produits
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "4 (suite) — Fichiers produits par les scripts",
       "Tous les resultats sont sauvegardes automatiquement dans scripts/")

outputs = [
    ("PK", "resultats_PK2comp_rxode2_FGFR2.RData",
     "Parametres PK (CL, V1, V2, Q, demi-vies)\nReutilise automatiquement par le script PKPD",
     LBLUE),
    ("PD", "resultats_PKPD_Emax_FGFR2.RData",
     "Parametres PKPD (l0, l1, Emax, EC50)\nReutilise automatiquement par le script presentation",
     LGREEN),
    ("PNG", "plot_PK2comp_rxode2_FGFR2.png",
     "Courbes PK ajustees (echelle log)", LGREY),
    ("PNG", "plot_PKPD_Emax_FGFR2.png",
     "Courbes PKPD ajustees (modele Emax vs donnees)", LGREY),
    ("PNG", "plot_Emax_curve_FGFR2.png",
     "Courbe effet-concentration (% Emax atteint a chaque dose)", LGREY),
    ("CSV", "tableau_TGI_FGFR2.csv",
     "Tableau TGI observe a chaque point de temps (exportable Excel)", LYELLOW),
    ("PPTX","presentation/presentation_FGFR2.pptx",
     "Presentation complete 11 diapositives (prete a utiliser)", LGREEN),
]

for i, (ftype, fname, desc, fill) in enumerate(outputs):
    row, col = divmod(i, 2)
    x = 0.3 + col * 6.6
    y = 1.25 + row * 1.75
    rect(sl, x, y, 6.4, 1.6, fill, BLUE)
    rect(sl, x, y, 0.9, 0.45, NAVY)
    txt(sl, ftype, x+0.01, y+0.04, 0.88, 0.38,
        sz=11, bold=True, col=WHITE, align=PP_ALIGN.CENTER)
    txt(sl, fname, x+1.0, y+0.07, 5.2, 0.32,
        sz=10, bold=True, col=NAVY)
    txt(sl, desc,  x+0.15, y+0.52, 6.1, 0.95, sz=10, col=GREY)

footer_note(sl,
    "Les fichiers .RData sont des fichiers R internes — ne pas ouvrir dans Excel. "
    "Les .png et .csv s'ouvrent normalement.")


# =============================================================================
# DIAPO 8 — Lire les résultats
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "5 — Comment interpreter les resultats",
       "Guide de lecture des graphiques et du tableau TGI")

# Fit PK
rect(sl, 0.3, 1.2, 6.0, 2.8, LBLUE, BLUE)
txt(sl, "Graphique PK (plot_PK2comp)", 0.45, 1.28, 5.7, 0.4,
    sz=12, bold=True, col=NAVY)
txt(sl,
    "Points = mesures experimentales\n"
    "Lignes = courbes predites par le modele\n\n"
    "Un bon ajustement : les points sont\n"
    "proches des lignes sur toute la duree\n\n"
    "Les 3 couleurs = les 3 doses\n"
    "(1, 5, 10 mg/kg)\n\n"
    "Echelle log : une droite = decroissance\n"
    "exponentielle pure (1 compartiment)\n"
    "Une courbe = 2 compartiments (normal)",
    0.45, 1.72, 5.7, 2.15, sz=11, col=GREY)

# Fit PKPD
rect(sl, 6.5, 1.2, 6.5, 2.8, LGREEN, GREEN)
txt(sl, "Graphique PKPD (plot_PKPD_Emax)", 6.65, 1.28, 6.2, 0.4,
    sz=12, bold=True, col=GREEN)
txt(sl,
    "3 groupes : Controle / 3 mg/kg / 10 mg/kg\n"
    "Points = donnees observees\n"
    "Lignes = simulation du modele Emax\n\n"
    "Bon ajustement : lignes proches des points\n\n"
    "Les lignes verticales en pointilles\n"
    "indiquent les jours d'injection (j0, j14...)\n\n"
    "Si le controle est bien ajuste mais pas\n"
    "les groupes traites → changer les bornes\n"
    "de Emax/EC50 dans le script",
    6.65, 1.72, 6.2, 2.15, sz=11, col=GREY)

# TGI
rect(sl, 0.3, 4.15, 6.0, 3.1, LYELLOW, ORANGE)
txt(sl, "Tableau TGI (Tumour Growth Inhibition)", 0.45, 4.23, 5.7, 0.4,
    sz=12, bold=True, col=ORANGE)
txt(sl,
    "TGI = (1 - delta_traite / delta_controle) x 100\n\n"
    "Interpretation :\n"
    "  TGI = 100% → regression complete\n"
    "  TGI =  80% → tres bon effet\n"
    "  TGI =  50% → effet modere\n"
    "  TGI =   0% → aucun effet\n"
    "  TGI < 0%   → la tumeur traitee a cru\n"
    "               plus vite que le controle\n"
    "               (aux temps precoces, normal)\n\n"
    "Vert = TGI >= 60%\n"
    "Orange = TGI entre 30 et 60%",
    0.45, 4.67, 5.7, 2.45, sz=11, col=GREY)

# Résidus
rect(sl, 6.5, 4.15, 6.5, 3.1, LGREY, BLUE)
txt(sl, "Graphique des residus (plot_residus_Emax)", 6.65, 4.23, 6.2, 0.4,
    sz=12, bold=True, col=NAVY)
txt(sl,
    "Residus = (observe - predit) / predit x 100\n\n"
    "Idéalement :\n"
    "  - Points repartis aleatoirement\n"
    "    autour de 0 (ligne grise)\n"
    "  - Pas de tendance systematique\n"
    "  - Majorite des points entre\n"
    "    les tirets a ±30%\n\n"
    "Si tous les residus sont positifs :\n"
    "  le modele sous-estime la tumeur\n"
    "Si tendance croissante dans le temps :\n"
    "  le modele ne capture pas bien\n"
    "  la dynamique tardive",
    6.65, 4.67, 6.2, 2.45, sz=11, col=GREY)

footer_note(sl,
    "La console R affiche aussi les valeurs numeriques : TGI predit a j49, Emax/l0, EC50/Cmax")


# =============================================================================
# DIAPO 9 — Paramètres : que signifient-ils biologiquement
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "5 (suite) — Signification biologique des parametres",
       "Ce que chaque nombre vous dit sur le medicament et la tumeur")

params_bio = [
    ("t½ beta", "Demi-vie d'elimination",
     "Temps pour que la concentration diminue de moitie.\n"
     "Ici ~10.6 jours → coherent avec un anticorps IgG\n"
     "(les anticorps sont elimines lentement par le corps)",
     LBLUE, NAVY),
    ("V1 (volume central)", "Distribution dans le sang",
     "Proche du volume plasmatique de la souris (~0.07 L/kg)\n"
     "→ confirme que l'Ac reste majoritairement dans le sang\n"
     "(pas de distribution tissulaire massive)",
     LBLUE, NAVY),
    ("l0 (croissance expo.)", "Vitesse de doublement de la tumeur",
     "Plus l0 est grand, plus la tumeur croit vite.\n"
     "Calibre sur le groupe controle sans traitement.\n"
     "Temps de doublement ≈ ln(2)/l0",
     LGREEN, GREEN),
    ("Emax / l0 > 1", "Le medicament peut regresser la tumeur",
     "Si Emax > l0 : le medicament peut tuer\n"
     "  les cellules plus vite qu'elles ne se divisent\n"
     "  → regression tumorale possible\n"
     "Si Emax < l0 : inhibition seulement, pas regression",
     LGREEN, GREEN),
    ("EC50 vs C1_max", "Positionnement sur la courbe effet-dose",
     "EC50 << C1_max (a 10 mg/kg) :\n"
     "  → on est au plafond a 10 mg/kg, augmenter la dose\n"
     "     ne sert plus a rien\n"
     "EC50 > C1_max :\n"
     "  → on n'a pas encore atteint l'effet maximal",
     RGBColor(0xFD,0xED,0xEC), RED),
    ("k1 (transit cellulaire)", "Vitesse de mort des cellules endommagees",
     "Temps moyen entre dommage et mort = 4/k1.\n"
     "Fixe a 0.5/jour → mort en ~8 jours.\n"
     "Valeur de la publication originale Simeoni 2004.",
     LYELLOW, ORANGE),
]

for i, (name, title, desc, fill, col) in enumerate(params_bio):
    row, c = divmod(i, 2)
    x = 0.3 + c * 6.55
    y = 1.2 + row * 1.98
    rect(sl, x, y, 6.35, 1.85, fill, col)
    txt(sl, name,  x+0.15, y+0.08, 6.0, 0.38, sz=12, bold=True, col=col)
    txt(sl, title, x+0.15, y+0.47, 6.0, 0.3,  sz=10, bold=True, col=GREY)
    txt(sl, desc,  x+0.15, y+0.78, 6.0, 1.0,  sz=10, col=GREY)

footer_note(sl,
    "Toutes ces valeurs sont affichees dans la console R apres execution du script "
    "et resumees dans la figure 09_parametres_PKPD.png")


# =============================================================================
# DIAPO 10 — Limites et quand refaire l'analyse
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "Limites du modele et quand refaire l'analyse",
       "Conditions dans lesquelles les resultats restent valides")

rect(sl, 0.3, 1.2, 6.0, 5.9, RGBColor(0xFD,0xED,0xEC), RED)
txt(sl, "Limites actuelles", 0.45, 1.28, 5.7, 0.4,
    sz=13, bold=True, col=RED)
limits = [
    "Un seul Emax/EC50 pour les 2 doses\n  (modele identique pour 3 et 10 mg/kg)",
    "Pas de variabilite inter-individuelle\n  (modele deterministe, pas de population)",
    "Donnees moyennes (mean±SEM)\n  pas les donnees individuelles par souris",
    "k1 fixe a la valeur de Simeoni 2004\n  (non estime sur ces donnees)",
    "Censure des controles a j35\n  (souris euthanasees, chute artificielle)",
    "Extrapolation au-dela de j49 incertaine\n  (pas de donnees pour valider)",
]
for i, l in enumerate(limits):
    txt(sl, "⚠  " + l, 0.45, 1.75+i*0.85, 5.7, 0.78, sz=11, col=GREY)

rect(sl, 6.5, 1.2, 6.5, 5.9, LGREEN, GREEN)
txt(sl, "Quand refaire l'analyse ?", 6.65, 1.28, 6.2, 0.4,
    sz=13, bold=True, col=GREEN)
cases = [
    ("Nouveau produit",
     "Remplacer les 2 fichiers Excel.\nRelancer pkpd_FGFR2.R puis\npresentation_FGFR2.R"),
    ("Nouveau protocole de dosage",
     "Changer DOSE_DAYS dans pkpd_FGFR2.R\n(ex: Q3W → c(0,21,42,63))"),
    ("Groupes de dose differents",
     "Changer dose3_pd et dose10_pd,\net les patterns Group 03/04"),
    ("Durée d'experience differente",
     "La censure T_CENSOR_CTRL s'ajuste\nautomatiquement si controle plus court"),
    ("Ameliorer le fit 3 mg/kg",
     "Essayer des bornes plus larges pour\nEC50 dans obj_emax()"),
]
for i, (title, action) in enumerate(cases):
    y = 1.75 + i * 0.98
    rect(sl, 6.6, y, 6.3, 0.9, LGREY)
    txt(sl, "→  " + title, 6.75, y+0.05, 6.0, 0.32, sz=11, bold=True, col=NAVY)
    txt(sl, action, 6.75, y+0.38, 6.0, 0.48, sz=10, col=GREY)

footer_note(sl,
    "En cas de doute, contacter [votre nom] — les scripts sont commentes "
    "et le depot GitHub contient l'historique de toutes les modifications")


# =============================================================================
# DIAPO 11 — Résumé une page
# =============================================================================
sl = prs.slides.add_slide(blank)
header(sl, "Resume — Vue d'ensemble du workflow",
       "De la donnee brute au resultat en 2 scripts R", color=GREEN)

rect(sl, 0.3, 1.2, 12.7, 5.9, LGREY, BLUE)

boxes = [
    (0.5,  1.35, 2.8, "DONNEES\nD'ENTREE",
     "PK souris FGFR2.xlsx\nTumorVolume_FGFR2.xlsx", LBLUE, NAVY),
    (3.7,  1.35, 2.8, "MODELE PK\n2 compartiments",
     "CL, V1, V2, Q\nt½β = 10.6 j\nnlminb (log-espace)", LBLUE, NAVY),
    (6.9,  1.35, 2.8, "MODELE PD\nSimeoni-Emax",
     "l0, l1 (controle)\nEmax, EC50\nDEoptim global", LGREEN, GREEN),
    (10.1, 1.35, 2.8, "RESULTATS",
     "9 figures PNG\nTableau TGI\nPowerPoint pret", LYELLOW, ORANGE),
]
for x, y, w, title, desc, fill, col in boxes:
    rect(sl, x, y+0.1, w, 1.6, fill, col)
    txt(sl, title, x+0.1, y+0.18, w-0.2, 0.55,
        sz=12, bold=True, col=col, align=PP_ALIGN.CENTER)
    txt(sl, desc, x+0.1, y+0.78, w-0.2, 0.85,
        sz=10, col=GREY, align=PP_ALIGN.CENTER)

# Flèches entre boîtes
for ax in [3.3, 6.5, 9.7]:
    txt(sl, "→", ax, 1.85, 0.4, 0.5, sz=22, bold=True, col=BLUE)

# Tableau récapitulatif hypothèses clés
rect(sl, 0.4, 3.15, 12.5, 3.75, WHITE, BLUE)
txt(sl, "Hypotheses cles retenues — justification",
    0.55, 3.23, 12.2, 0.4, sz=13, bold=True, col=NAVY)

hyps_final = [
    ("Modele PK",    "2 compartiments",    "Decroissance bi-exponentielle visible dans les donnees"),
    ("Cinetique PK", "Lineaire en dose",   "Ratio concentrations proportionnel aux doses (verifie)"),
    ("Croissance",   "Simeoni (l0/l1)",    "Croissance expo au debut, lineaire a grand volume (standard)"),
    ("Effet drug",   "Emax sature",        "TGI 3 et 10 mg/kg non proportionnel a la dose → saturation"),
    ("Transit",      "k1 = 0.5/j fixe",   "Non identifiable sur ces donnees — valeur Simeoni 2004"),
    ("Hill",         "p = 1 fixe",         "Valeur originale Simeoni 2004, evite instabilite numerique"),
]
hdr_cols = ["Composante", "Choix", "Justification"]
for j, h in enumerate(hdr_cols):
    xp = [0.55, 2.65, 5.05][j]
    txt(sl, h, xp, 3.68, [1.9, 2.3, 7.5][j], 0.32,
        sz=10, bold=True, col=NAVY)
for i, (comp, choix, just) in enumerate(hyps_final):
    y2 = 4.05 + i * 0.47
    bg = LGREY if i % 2 == 0 else WHITE
    rect(sl, 0.45, y2, 12.4, 0.43, bg)
    txt(sl, comp,  0.55, y2+0.07, 1.9, 0.34, sz=10, bold=True, col=NAVY)
    txt(sl, choix, 2.65, y2+0.07, 2.3, 0.34, sz=10, bold=True, col=BLUE)
    txt(sl, just,  5.05, y2+0.07, 7.5, 0.34, sz=10, col=GREY)

footer_note(sl,
    "pkpd_FGFR2.R  →  resultats_PKPD_Emax_FGFR2.RData  →  presentation_FGFR2.R  →  presentation_FGFR2.pptx")


# =============================================================================
# Sauvegarde
# =============================================================================
prs.save(OUT)
print(f"Guide genere : {OUT}")
print(f"Diapositives : {len(prs.slides)}")
