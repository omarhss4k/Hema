"""
Génère presentation_FGFR2.pptx
9 diapositives narratives pour la tutrice
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import os

IMG_DIR = "scripts/presentation"
OUT     = "scripts/presentation/presentation_FGFR2.pptx"

# Palette couleurs
NAVY   = RGBColor(0x1A, 0x52, 0x76)
WHITE  = RGBColor(0xFF, 0xFF, 0xFF)
GREY   = RGBColor(0x64, 0x64, 0x64)
GREEN  = RGBColor(0x1E, 0x84, 0x49)

SLIDE_W = Inches(13.33)
SLIDE_H = Inches(7.5)

prs = Presentation()
prs.slide_width  = SLIDE_W
prs.slide_height = SLIDE_H

blank = prs.slide_layouts[6]   # layout vide


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def add_rect(slide, l, t, w, h, fill_rgb, alpha=None):
    shape = slide.shapes.add_shape(1, Inches(l), Inches(t), Inches(w), Inches(h))
    shape.line.fill.background()
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_rgb
    shape.line.color.rgb = fill_rgb
    return shape


def add_text(slide, text, l, t, w, h,
             font_size=18, bold=False, color=RGBColor(0x1A,0x52,0x76),
             align=PP_ALIGN.LEFT, italic=False):
    txb = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf  = txb.text_frame
    tf.word_wrap = True
    p   = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size  = Pt(font_size)
    run.font.bold  = bold
    run.font.color.rgb = color
    run.font.italic    = italic
    return txb


def add_image(slide, path, l, t, w, h=None):
    if h is None:
        pic = slide.shapes.add_picture(path, Inches(l), Inches(t), width=Inches(w))
    else:
        pic = slide.shapes.add_picture(path, Inches(l), Inches(t),
                                       width=Inches(w), height=Inches(h))
    return pic


def title_slide(slide, title, subtitle=""):
    """Bande de titre en haut."""
    add_rect(slide, 0, 0, 13.33, 1.15, NAVY)
    add_text(slide, title,
             0.25, 0.08, 12.5, 0.85,
             font_size=26, bold=True, color=WHITE, align=PP_ALIGN.LEFT)
    if subtitle:
        add_text(slide, subtitle,
                 0.25, 0.85, 12.5, 0.35,
                 font_size=13, bold=False, color=RGBColor(0xAE,0xD6,0xF1),
                 align=PP_ALIGN.LEFT)


def footer(slide, note):
    add_text(slide, note,
             0.2, 7.1, 12.9, 0.35,
             font_size=10, color=GREY, align=PP_ALIGN.LEFT, italic=True)


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 0 — Titre général
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
add_rect(sl, 0, 0, 13.33, 7.5, NAVY)

add_text(sl,
         "Modélisation PK/PD",
         1.5, 1.5, 10.5, 1.2,
         font_size=38, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

add_text(sl,
         "Fc-silent FGFR2-huBPA-LP1",
         1.5, 2.7, 10.5, 0.8,
         font_size=26, bold=False, color=RGBColor(0xAE,0xD6,0xF1),
         align=PP_ALIGN.CENTER)

add_text(sl,
         "Modele 2-compartiments PK  |  Simeoni 2004 PKPD  |  Souris xenogreffe",
         1.0, 3.6, 11.0, 0.6,
         font_size=16, bold=False, color=RGBColor(0x85,0xC1,0xE9),
         align=PP_ALIGN.CENTER)

add_text(sl,
         "Bloc 1 : Donnees d'entree     Bloc 2 : Structure des modeles\n"
         "Bloc 3 : Parametres estimes   Bloc 4 : Predit vs Observe     Bloc 5 : Efficacite (TGI)",
         0.8, 4.6, 12.0, 1.2,
         font_size=13, color=RGBColor(0xD6,0xEA,0xF8), align=PP_ALIGN.CENTER)

add_text(sl, "Mai 2026", 5.5, 6.7, 2.5, 0.4,
         font_size=11, color=RGBColor(0x85,0xC1,0xE9), align=PP_ALIGN.CENTER)


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 1 — Données d'entrée PK
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 1 — Donnees d'entree : Pharmacocinetique",
            "Concentrations plasmatiques observees — IV bolus | 3 doses | Souris")

add_text(sl,
         "Donnees brutes :",
         0.25, 1.25, 4.0, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         u"• 3 doses testees : 1, 5 et 10 mg/kg\n"
         u"• Decroissance bi-exponentielle visible\n"
         u"• Compatible avec un modele 2-compartiments\n"
         u"• Lineaire en dose (courbes proportionnelles)",
         0.25, 1.65, 4.0, 1.8,
         font_size=12, color=GREY)

add_image(sl, os.path.join(IMG_DIR, "01_donnees_PK_brutes.png"),
          0.2, 3.55, 12.9, 3.8)

footer(sl, "Echelle log (droite) : deux phases de decroissance distinctes — distribution rapide puis elimination lente")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 2 — Données tumorales brutes
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 1 — Donnees d'entree : Croissance tumorale vs Inhibition",
            "Volume tumoral (mean +/- SEM) | Xenogreffe souris | Protocole Q2W x4")

add_text(sl,
         "Observations cles :",
         0.25, 1.25, 4.5, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         u"• Controle : croissance rapide, plateau vers j21 (~1200 mm3)\n"
         u"• 3 mg/kg : inhibition partielle, tumeur stable apres j18\n"
         u"• 10 mg/kg : stase quasi-complete, legere regression\n"
         u"• Effet dose-dependant clairement visible",
         0.25, 1.65, 4.5, 1.8,
         font_size=12, color=GREY)

add_image(sl, os.path.join(IMG_DIR, "05_croissance_vs_inhibition.png"),
          0.15, 1.2, 8.7, 6.1)

add_text(sl, "Effet du traitement", 9.2, 1.5, 3.9, 0.4,
         font_size=13, bold=True, color=NAVY)

bullets_effect = [
    ("Controle",  "Croissance exponentielle", GREY),
    ("3 mg/kg",   "TGI ~ 64% a j25",          GREEN),
    ("10 mg/kg",  "TGI ~ 84% a j25",          RGBColor(0x1B,0x4F,0x9E)),
]
for i, (dose, txt, col) in enumerate(bullets_effect):
    add_rect(sl, 9.2, 2.1 + i*0.85, 3.9, 0.7, RGBColor(0xF4,0xF6,0xF7))
    add_text(sl, dose, 9.35, 2.12 + i*0.85, 1.5, 0.35,
             font_size=12, bold=True, color=col)
    add_text(sl, txt, 9.35, 2.45 + i*0.85, 3.7, 0.3,
             font_size=11, color=GREY)

footer(sl, "Les triangles rouges indiquent les jours d'administration (j0, j14, j28, j42)")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 3 — Schéma modèle PK
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 2 — Structure du modele PK : 2-compartiments",
            "IV bolus | Parametres : CL, V1, V2, Q | Optimisation nlminb en log-espace")

add_image(sl, os.path.join(IMG_DIR, "04_schema_modele_PK.png"),
          0.3, 1.2, 8.5, 5.8)

add_text(sl, "Equations du modele", 9.1, 1.5, 4.0, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         "dA1/dt = -(CL/V1 + Q/V1)*A1\n        + (Q/V2)*A2\n\n"
         "dA2/dt =  (Q/V1)*A1 - (Q/V2)*A2\n\n"
         "C1(t) = A1(t) / V1",
         9.1, 1.95, 4.0, 2.0,
         font_size=11, color=GREY)

add_text(sl, "Parametres estimes", 9.1, 4.2, 4.0, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         u"CL  Clairance systemique\n"
         u"V1  Volume central\n"
         u"V2  Volume peripherique\n"
         u"Q   Clairance intercompart.\n"
         u"C1(t) pilote l'effet PD",
         9.1, 4.65, 4.0, 1.8,
         font_size=11, color=GREY)

footer(sl, "Methode des residus pour les valeurs initiales → optimisation nlminb pour affiner les parametres")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 4 — Schéma Simeoni PKPD
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 2 — Structure du modele PKPD : Simeoni 2004",
            "Croissance tumorale + compartiments de transit pour les cellules endommagees")

add_image(sl, os.path.join(IMG_DIR, "06_schema_modele_Simeoni.png"),
          0.15, 1.15, 9.0, 6.2)

add_text(sl, "Logique du modele", 9.4, 1.5, 3.7, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         u"1. C1(t) issue du modele PK\n\n"
         u"2. C1 tue les cellules x1\n   via k2 * C1 * x1\n\n"
         u"3. Cellules endommagees\n   transitent x2->x3->x4\n\n"
         u"4. w(t) = x1+x2+x3+x4\n   = volume tumoral total\n\n"
         u"5. Retroaction : la masse w\n   limite la croissance de x1",
         9.4, 1.95, 3.7, 4.5,
         font_size=11, color=GREY)

footer(sl, "Reference : Simeoni et al., Cancer Research 2004 | p=1 (modele original)")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 5 — Paramètres PK estimés
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 3 — Sorties du modele PK : Parametres estimes",
            "Modele 2-compartiments | IV bolus | Fc-silent FGFR2-huBPA-LP1 | Souris")

add_image(sl, os.path.join(IMG_DIR, "03_parametres_PK.png"),
          0.3, 1.2, 8.5, 5.8)

add_text(sl, "Points cles", 9.1, 1.5, 4.0, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         u"t1/2 beta = 10.6 jours\n"
         u"-> Longue demi-vie\n"
         u"   typique d'un anticorps\n\n"
         u"V1 = 0.071 L/kg\n"
         u"-> Volume central proche\n"
         u"   du plasma (souris)\n\n"
         u"Vss = V1+V2 = 0.204 L/kg\n"
         u"-> Distribution moderee\n\n"
         u"CL faible = 0.017 L/j/kg\n"
         u"-> Elimination lente",
         9.1, 1.95, 4.0, 4.5,
         font_size=11, color=GREY)

footer(sl, "Optimisation nlminb en log-espace | Residus log egalement ponderes entre les 3 groupes de dose")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 6 — Paramètres PKPD
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 3 — Sorties du modele PKPD : Parametres estimes",
            "Simeoni 2004 | PK fixes | l0/l1 calibres sur controle | k2 optimise (DEoptim)")

add_image(sl, os.path.join(IMG_DIR, "09_parametres_PKPD.png"),
          0.3, 1.2, 12.7, 5.8)

footer(sl,
       "Vert = parametres PD estimes | Bleu = parametres PK fixes depuis le modele PK "
       "| k1 fixe a 0.5/j (Simeoni 2004) | k2 optimise par DEoptim")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 7 — Fit PK prédit vs observé
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 4 — Qualite du fit PK : Predit vs Observe",
            "3 doses | Echelle log | Lignes = modele | Points = donnees observees")

add_image(sl, os.path.join(IMG_DIR, "02_PK_predit_vs_observe.png"),
          0.3, 1.2, 8.8, 5.9)

add_text(sl, "Interpretation", 9.3, 1.5, 3.8, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         u"Bon ajustement global :\n\n"
         u"• Les 3 doses sont bien\n"
         u"  reproduites sur 650h\n\n"
         u"• Linearite PK confirmee\n"
         u"  (courbes proportionnelles)\n\n"
         u"• Deux phases distinctes :\n"
         u"  - Phase alpha (rapide)\n"
         u"    t1/2a ~ quelques heures\n"
         u"  - Phase beta (lente)\n"
         u"    t1/2b = 10.6 jours\n\n"
         u"• Legere deviation a 1 mg/kg\n"
         u"  en fin de courbe",
         9.3, 1.95, 3.8, 4.8,
         font_size=11, color=GREY)

footer(sl, "Objectif nlminb : somme des residus quadratiques log | Residus log egalement ponderes par groupe")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 8 — Fit PKPD prédit vs observé
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 4 — Qualite du fit PKPD : Predit vs Observe",
            "Simeoni 2004 | Q2W x4 doses | k2 = 1.00e-07 | Residus relatifs (droite)")

add_image(sl, os.path.join(IMG_DIR, "07_PKPD_predit_vs_observe.png"),
          0.15, 1.2, 12.9, 5.3)

add_text(sl, "Limite actuelle : k2 unique pour les deux doses traitees",
         0.25, 6.6, 10.0, 0.5,
         font_size=12, bold=True, color=RGBColor(0xC0,0x39,0x2B))
add_text(sl,
         u"-> 10 mg/kg bien capture | 3 mg/kg sous-estime -> necessite un modele Emax ou k2 par dose",
         0.25, 6.95, 12.5, 0.4,
         font_size=11, color=GREY)


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 9 — TGI tableau
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Bloc 5 — Efficacite : Tumour Growth Inhibition (TGI)",
            "TGI observe vs predit par le modele | Protocole Q2W | Fc-silent FGFR2-huBPA-LP1")

add_image(sl, os.path.join(IMG_DIR, "08_tableau_TGI.png"),
          0.3, 1.2, 8.5, 5.8)

add_text(sl, "Conclusions", 9.1, 1.5, 4.0, 0.4,
         font_size=13, bold=True, color=NAVY)
add_text(sl,
         u"10 mg/kg :\n"
         u"  TGI observe  : 83.5 %\n"
         u"  TGI predit   : 93.8 %\n"
         u"  -> Bon accord modele/donnees\n\n"
         u"3 mg/kg :\n"
         u"  TGI observe  : 64.2 %\n"
         u"  TGI predit   : 47.0 %\n"
         u"  -> Modele sous-estime l'effet\n\n"
         u"Piste d'amelioration :\n"
         u"  Modele Emax pour capturer\n"
         u"  la non-linearite de l'effet\n"
         u"  en fonction de la dose",
         9.1, 1.95, 4.0, 5.0,
         font_size=11, color=GREY)

footer(sl, "TGI = (1 - delta_W_traite / delta_W_controle) x 100 | Vert >= 60% | Orange 30-60%")


# ─────────────────────────────────────────────────────────────────────────────
# DIAPO 10 — Conclusions et perspectives
# ─────────────────────────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
title_slide(sl,
            "Conclusions et perspectives",
            "Fc-silent FGFR2-huBPA-LP1 | Modelisation PK/PD souris")

add_rect(sl, 0.3, 1.25, 5.9, 5.85, RGBColor(0xEA,0xF2,0xFB))
add_text(sl, "Ce qui fonctionne", 0.5, 1.35, 5.5, 0.45,
         font_size=14, bold=True, color=NAVY)

ok_items = [
    "Modele PK 2-compartiments : bon fit, 3 doses",
    "t1/2 beta = 10.6 j — coherent avec un Ac",
    "Linearite PK confirmee",
    "Croissance tumorale controle bien ajustee",
    "TGI predit 10 mg/kg ~ TGI observe (84 vs 94%)",
    "Tableau TGI : outil de communication clair",
]
for i, item in enumerate(ok_items):
    add_text(sl, u"✓  " + item,
             0.5, 1.9 + i*0.72, 5.5, 0.55,
             font_size=11, color=RGBColor(0x1E,0x84,0x49))

add_rect(sl, 6.8, 1.25, 6.2, 5.85, RGBColor(0xFD,0xED,0xEC))
add_text(sl, "Limites et prochaines etapes", 7.0, 1.35, 5.8, 0.45,
         font_size=14, bold=True, color=RGBColor(0xC0,0x39,0x2B))

next_items = [
    "k2 unique => sous-estime l'effet a 3 mg/kg",
    "Piste 1 : k2 libre par groupe de dose",
    "Piste 2 : modele Emax (Cmax, EC50)",
    "Extension : translabilite souris -> homme",
    "Simulation de nouveaux schemas posologiques",
]
for i, item in enumerate(next_items):
    add_text(sl, u"->  " + item,
             7.0, 1.9 + i*0.88, 5.8, 0.7,
             font_size=11, color=RGBColor(0x92,0x2B,0x21))

footer(sl, "Modele PK/PD developpe sous R (rxode2 + deSolve + DEoptim) | Simeoni 2004")


# ─────────────────────────────────────────────────────────────────────────────
# Sauvegarde
# ─────────────────────────────────────────────────────────────────────────────
prs.save(OUT)
print(f"PowerPoint genere : {OUT}")
print(f"Diapositives : {len(prs.slides)}")
