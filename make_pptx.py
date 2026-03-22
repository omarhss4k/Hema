"""
Génération du PowerPoint de réunion d'avancement
Projet : Modélisation PK/PD de la toxicité hématologique du carboplatine
Référence : Fornari 2019
"""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import os

# ── Couleurs ─────────────────────────────────────────────
BLEU_FONCE  = RGBColor(0x1A, 0x3A, 0x5C)
BLEU_CLAIR  = RGBColor(0x2E, 0x75, 0xB6)
ORANGE      = RGBColor(0xC5, 0x50, 0x1A)
BLANC       = RGBColor(0xFF, 0xFF, 0xFF)
GRIS_CLAIR  = RGBColor(0xF2, 0xF2, 0xF2)
VERT        = RGBColor(0x37, 0x86, 0x44)
ROUGE       = RGBColor(0xC0, 0x00, 0x00)
NOIR        = RGBColor(0x1A, 0x1A, 0x1A)

prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)

BLANK = prs.slide_layouts[6]  # layout vide

def add_slide():
    return prs.slides.add_slide(BLANK)

def rect(slide, x, y, w, h, fill=None, line=None, line_w=Pt(0)):
    from pptx.util import Pt
    shp = slide.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.line.width = line_w
    if fill:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    else:
        shp.fill.background()
    if line:
        shp.line.color.rgb = line
        shp.line.width = line_w if line_w else Pt(1)
    else:
        shp.line.fill.background()
    return shp

def txbox(slide, text, x, y, w, h,
          size=18, bold=False, color=NOIR, align=PP_ALIGN.LEFT,
          italic=False, wrap=True):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = wrap
    p  = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return tb

def header_bar(slide, titre, sous_titre=None):
    """Barre bleue foncée en haut + titre blanc"""
    rect(slide, 0, 0, 13.33, 1.1, fill=BLEU_FONCE)
    txbox(slide, titre, 0.3, 0.05, 12.5, 0.7,
          size=28, bold=True, color=BLANC, align=PP_ALIGN.LEFT)
    if sous_titre:
        txbox(slide, sous_titre, 0.3, 0.72, 12.5, 0.35,
              size=14, color=RGBColor(0xBD, 0xD7, 0xEE), align=PP_ALIGN.LEFT)
    # Ligne orange en bas de header
    rect(slide, 0, 1.05, 13.33, 0.06, fill=ORANGE)

def footer(slide, page, total):
    rect(slide, 0, 7.2, 13.33, 0.3, fill=BLEU_FONCE)
    txbox(slide, "Modélisation PK/PD — Carboplatine & Hématotoxicité",
          0.2, 7.2, 10, 0.28, size=9, color=BLANC)
    txbox(slide, f"{page}/{total}",
          12.8, 7.2, 0.5, 0.28, size=9, color=BLANC, align=PP_ALIGN.RIGHT)

def bullet_block(slide, lines, x, y, w, h,
                 size=14, indent=True, color=NOIR):
    """Bloc de bullets (liste de tuples (niveau, texte))"""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    first = True
    for (lvl, text) in lines:
        if first:
            p = tf.paragraphs[0]
            first = False
        else:
            p = tf.add_paragraph()
        p.level = lvl
        bullet = "• " if lvl == 0 else "  – " if lvl == 1 else "      · "
        run = p.add_run()
        run.text = bullet + text
        run.font.size = Pt(size if lvl == 0 else size - 1)
        run.font.bold = (lvl == 0)
        run.font.color.rgb = color if lvl == 0 else RGBColor(0x40, 0x40, 0x40)

N_SLIDES = 8

# ════════════════════════════════════════════════════════
# SLIDE 1 — Page de titre
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLEU_FONCE)
rect(sl, 0, 0, 13.33, 0.08, fill=ORANGE)
rect(sl, 0, 7.42, 13.33, 0.08, fill=ORANGE)

txbox(sl, "Modélisation de la toxicité hématologique",
      1, 1.6, 11.3, 1.0, size=34, bold=True, color=BLANC, align=PP_ALIGN.CENTER)
txbox(sl, "du Carboplatine — Réunion d'avancement",
      1, 2.55, 11.3, 0.8, size=28, bold=False, color=RGBColor(0xBD, 0xD7, 0xEE),
      align=PP_ALIGN.CENTER)

rect(sl, 3.5, 3.5, 6.3, 0.05, fill=ORANGE)

txbox(sl, "Référence : Fornari et al. 2019 — CPT:PSP",
      1, 3.7, 11.3, 0.5, size=16, color=RGBColor(0xBD, 0xD7, 0xEE),
      align=PP_ALIGN.CENTER)
txbox(sl, "Mars 2026",
      1, 6.4, 11.3, 0.5, size=14, color=RGBColor(0x9D, 0xB3, 0xCC),
      align=PP_ALIGN.CENTER)

footer(sl, 1, N_SLIDES)

# ════════════════════════════════════════════════════════
# SLIDE 2 — Contexte & Objectifs
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLANC)
header_bar(sl, "Contexte & Objectifs",
           "Reproduire et adapter le modèle Fornari 2019 chez l'humain")

# Colonne gauche — Contexte
rect(sl, 0.3, 1.25, 5.9, 5.6, fill=GRIS_CLAIR, line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Problème clinique", 0.5, 1.3, 5.5, 0.4,
      size=13, bold=True, color=BLEU_FONCE)
bullet_block(sl, [
    (0, "Carboplatine = chimiothérapie standard"),
    (1, "Dosage par formule de Calvert (AUC cible)"),
    (0, "Toxicité hématologique dose-limitante"),
    (1, "Neutropénie & thrombocytopénie"),
    (1, "Grades NCI-CTCAE : G1 → G4"),
    (0, "Besoin d'un modèle prédictif"),
    (1, "Anticiper les nadirs par patient"),
    (1, "Adapter la dose selon la fonction rénale"),
], 0.5, 1.75, 5.7, 4.8)

# Colonne droite — Objectifs
rect(sl, 6.5, 1.25, 6.5, 5.6, fill=RGBColor(0xE8, 0xF1, 0xFA),
     line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Objectifs du stage", 6.7, 1.3, 6.1, 0.4,
      size=13, bold=True, color=BLEU_FONCE)
bullet_block(sl, [
    (0, "1. Reproduire le modèle Fornari 2019"),
    (1, "ODEs 1–9 (hématopoïèse + PK/PD)"),
    (1, "Validation rat (Figure 3 papier)"),
    (0, "2. Transposer chez l'humain"),
    (1, "Scaling allométrique des paramètres"),
    (1, "PK carboplatine 2-compartiments"),
    (0, "3. Simuler la population humaine"),
    (1, "VPC — 1000 patients (GFR variable)"),
    (1, "Distribution des grades par AUC cible"),
    (0, "4. Calibrer le modèle"),
    (1, "Ajuster δ_Plt & γ pour reproduire"),
    (1, "les grades observés (Fornari Fig. 4c)"),
], 6.7, 1.75, 6.1, 4.8)

footer(sl, 2, N_SLIDES)

# ════════════════════════════════════════════════════════
# SLIDE 3 — Structure du modèle
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLANC)
header_bar(sl, "Structure du modèle PK/PD",
           "Système d'EDOs couplé — Hématopoïèse + Pharmacocinétique carboplatine")

# Schéma simplifié (boîtes et flèches textuelles)
boxes = [
    # (x, y, w, h, fill, texte, taille)
    (0.3,  1.3, 2.5, 0.7, BLEU_FONCE,  "PK Carboplatine\n2 compartiments", 11),
    (0.3,  2.3, 2.5, 0.7, BLEU_CLAIR,  "Damage (Eq. 3)\nadducts ADN", 11),
    (3.2,  1.3, 2.4, 0.7, RGBColor(0x37,0x86,0x44), "MPP\nProg. multipotente", 11),
    (3.2,  2.3, 2.4, 0.7, RGBColor(0x37,0x86,0x44), "CMP / MEP\nProg. engagées", 11),
    (6.0,  1.3, 2.4, 0.7, ORANGE,       "Neut / Mono\n(voie myéloïde)", 11),
    (6.0,  2.3, 2.4, 0.7, ORANGE,       "Ret / RBC\n(voie érythroïde)", 11),
    (8.8,  1.3, 2.4, 0.7, RGBColor(0xC0,0x00,0x00), "Plaquettes (Plt)\nThrombocytes", 11),
    (8.8,  2.3, 2.4, 0.7, BLEU_FONCE,  "Feedbacks\nEq. 11–13", 11),
]
for (bx, by, bw, bh, bf, bt, bs) in boxes:
    rect(sl, bx, by, bw, bh, fill=bf)
    txbox(sl, bt, bx+0.05, by+0.05, bw-0.1, bh-0.1,
          size=bs, color=BLANC, align=PP_ALIGN.CENTER, bold=True)

# Flèches (texte →)
for (ax, ay, atxt) in [
    (2.85, 1.58, "→"), (2.85, 2.58, "→"),
    (5.65, 1.58, "→"), (5.65, 2.58, "→"),
    (8.45, 1.58, "→"),
]:
    txbox(sl, atxt, ax, ay, 0.4, 0.3, size=18, color=BLEU_FONCE, bold=True)

# Effet drogue
rect(sl, 3.2, 3.3, 7.2, 0.6, fill=RGBColor(0xFF, 0xEB, 0xD6))
txbox(sl, "Effet drogue : Kill = Slope × Damage  (MPP, CMP, MEP)   +   δ × Damage  (Ret, Plt)",
      3.4, 3.35, 7.0, 0.5, size=12, color=ROUGE, bold=True)

# Paramètres clés
rect(sl, 0.3, 4.1, 12.5, 2.9, fill=GRIS_CLAIR, line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Paramètres clés implémentés", 0.5, 4.15, 6, 0.4,
      size=13, bold=True, color=BLEU_FONCE)

col1 = [
    (0, "PK : CL=6.18 L/h, V1=7.87 L, Q=1.98 L/h, V2=8.06 L"),
    (0, "Baselines : Plt0=345 · 10⁹/L, Neut0=4.5 · 10⁹/L"),
    (0, "MTT_Plt = 168 h, MTT_Neut = 210 h"),
]
col2 = [
    (0, "Slopes : MPP=0.79, CMP=0.57, MEP=0.66"),
    (0, "δ_Plt = 0.80  (calibré), δ_Ret = 2.8"),
    (0, "GFR ~ N(78, 20²) mL/min  pour la population"),
]
bullet_block(sl, col1, 0.5, 4.6, 6.1, 2.2, size=12)
bullet_block(sl, col2, 6.7, 4.6, 6.0, 2.2, size=12)

footer(sl, 3, N_SLIDES)

# ════════════════════════════════════════════════════════
# SLIDE 4 — Validation rat (Figure 3)
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLANC)
header_bar(sl, "Validation : modèle rat (Figure 3 Fornari 2019)",
           "Reproduction des courbes Neut, Plt, Ret — données rat carboplatine")

rect(sl, 0.3, 1.25, 12.4, 4.0, fill=RGBColor(0xE8, 0xF1, 0xFA),
     line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Résultats rat (Figure 3 du papier reproduite)", 0.5, 1.3, 10, 0.4,
      size=13, bold=True, color=BLEU_FONCE)

bullet_block(sl, [
    (0, "Simulation déterministe : patient rat typique, dose unique carboplatine"),
    (1, "Courbes Neut, Plt, Ret, RBC reproduites → nadirs & temps de récupération conformes au papier"),
    (0, "VPC (Supplementary Figure S1) : 1000 rats simulés avec variabilité inter-individuelle"),
    (1, "Bandes 5–95e percentile encadrent les données observées (rat)"),
    (0, "Grades NCI-CTCAE calculés sur la population simulée"),
    (1, "Distribution G0–G4 cohérente avec les données rat"),
    (0, "Validations techniques : convergence EDOs (LSODA, rtol=1e-7), diagnostics PK/PD"),
], 0.5, 1.75, 12.0, 3.4)

# Case résultat
rect(sl, 0.3, 5.4, 12.4, 1.7, fill=RGBColor(0xE2, 0xEF, 0xDA), line=VERT, line_w=Pt(1.5))
txbox(sl, "Conclusion — Étape rat",
      0.5, 5.45, 4, 0.35, size=12, bold=True, color=VERT)
txbox(sl,
      "Le modèle ODE (Éqs 1–9 + feedbacks 11–13) reproduit fidèlement les données rat de Fornari 2019. "
      "Les fichiers pkpd_model_FORNARI.R, parameters_rat.R, run_CORRECT.R, plots.R sont finalisés.",
      0.5, 5.85, 12.0, 1.1, size=12, color=NOIR)

footer(sl, 4, N_SLIDES)

# ════════════════════════════════════════════════════════
# SLIDE 5 — Simulation humaine déterministe (Figure 4)
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLANC)
header_bar(sl, "Simulation humaine — Patient type (Figure 4)",
           "Carboplatine AUC=5, schéma Q21D × 2 cycles, GFR=125 mL/min")

# 2 colonnes
rect(sl, 0.3, 1.25, 6.0, 5.7, fill=GRIS_CLAIR, line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Hypothèses", 0.5, 1.3, 5.6, 0.4, size=13, bold=True, color=BLEU_FONCE)
bullet_block(sl, [
    (0, "Dose Calvert : AUC × (GFR + 25)"),
    (1, "AUC=5, GFR=125 → Dose = 750 mg"),
    (0, "Infusion 1h, 2 cycles à J0 & J21"),
    (0, "Patient médian (paramètres Table 1 & 2)"),
    (0, "PK 2 compartiments (Zandvliet 2008)"),
    (0, "Fraction libre fu₀=fu∞=1 (carboplatine)"),
    (0, "Pas de δ_Neut (absent du modèle Fornari)"),
], 0.5, 1.75, 5.7, 3.5)

rect(sl, 6.7, 1.25, 6.3, 5.7, fill=RGBColor(0xE8, 0xF1, 0xFA),
     line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Résultats déterministes", 6.9, 1.3, 5.9, 0.4,
      size=13, bold=True, color=BLEU_FONCE)
bullet_block(sl, [
    (0, "Neut : nadir Cycle 1 ≈ 2.5 · 10⁹/L → Grade 1"),
    (1, "Nadir C2 légèrement inférieur (cumul)"),
    (0, "Plt : nadir C1 ≈ 160 · 10⁹/L → Grade 1"),
    (1, "Nadir C2 ≈ 146 · 10⁹/L → Grade 1"),
    (0, "Ret & RBC : légère suppression érythroïde"),
    (0, "Récupération complète entre les cycles"),
    (0, "Conforme à Figure 4 du papier Fornari"),
], 6.9, 1.75, 5.9, 3.5)

# Note grades
rect(sl, 0.3, 6.05, 12.4, 0.85, fill=RGBColor(0xFF, 0xF0, 0xCC), line=ORANGE, line_w=Pt(1))
txbox(sl, "Note : La courbe déterministe = patient médian (Grade 1 attendu). "
          "Les grades G3/G4 apparaissent dans la queue de distribution de la population simulée (voir VPC & Figure 4c).",
      0.5, 6.1, 12.0, 0.75, size=11, color=RGBColor(0x80, 0x40, 0x00))

footer(sl, 5, N_SLIDES)

# ════════════════════════════════════════════════════════
# SLIDE 6 — VPC & Grades population (Figure 4b/4c)
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLANC)
header_bar(sl, "Simulation population — VPC & Grades (Figure 4b/4c)",
           "1000 patients simulés, GFR ~ N(78, 20²) mL/min")

rect(sl, 0.3, 1.25, 6.0, 5.7, fill=GRIS_CLAIR, line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "VPC (Visual Predictive Check)", 0.5, 1.3, 5.6, 0.4,
      size=13, bold=True, color=BLEU_FONCE)
bullet_block(sl, [
    (0, "Variabilité inter-individuelle simulée"),
    (1, "GFR ~ N(78, 20²) → CL individuelle (Calvert)"),
    (1, "IIV sur baselines (CV 20–30%)"),
    (0, "Bandes 5e–95e percentile tracées"),
    (1, "Encadrent la médiane et les données Fornari"),
    (0, "AUC testées : 4, 5, 6 mg/mL·min"),
    (0, "Figures : Figure4_VPC_AUC4/5/6.pdf"),
], 0.5, 1.75, 5.7, 4.0)

rect(sl, 6.7, 1.25, 6.3, 5.7, fill=RGBColor(0xE8, 0xF1, 0xFA),
     line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Distribution des grades (Figure 4c)", 6.9, 1.3, 5.9, 0.4,
      size=13, bold=True, color=BLEU_FONCE)
bullet_block(sl, [
    (0, "% patients par grade NCI-CTCAE (Plt & Neut)"),
    (0, "AUC=4 : majorité G0–G1, très peu G3/G4"),
    (0, "AUC=5 : montée des G2–G3"),
    (1, "Plt : ~15% G2, ~5% G3"),
    (0, "AUC=6 : G3/G4 significatifs"),
    (1, "Nécessite ajustement de dose"),
    (0, "Calibration δ_Plt=0.80, γ=0.40"),
    (1, "Reproduit les proportions de Fornari Fig 4c"),
], 6.9, 1.75, 5.9, 4.0)

rect(sl, 0.3, 6.05, 12.4, 0.85, fill=RGBColor(0xE2, 0xEF, 0xDA), line=VERT, line_w=Pt(1))
txbox(sl, "Fichiers générés : Figure4_VPC.pdf, Figure4c_AUC4/5/6.pdf, Figure4c_grades_AUC5.pdf",
      0.5, 6.1, 12.0, 0.75, size=12, color=VERT, bold=True)

footer(sl, 6, N_SLIDES)

# ════════════════════════════════════════════════════════
# SLIDE 7 — Calibration & Difficultés rencontrées
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLANC)
header_bar(sl, "Calibration du modèle & Difficultés rencontrées",
           "Ajustement itératif des paramètres pour reproduire Fornari 2019")

rect(sl, 0.3, 1.25, 6.0, 5.7, fill=GRIS_CLAIR, line=BLEU_CLAIR, line_w=Pt(1))
txbox(sl, "Processus de calibration", 0.5, 1.3, 5.6, 0.4,
      size=13, bold=True, color=BLEU_FONCE)
bullet_block(sl, [
    (0, "Cible : nadirs Plt ≈ 160 (C1) & 146 (C2)"),
    (0, "Paramètre ajusté : δ_Plt (effet drogue Plt)"),
    (1, "Résultat : δ_Plt = 0.80"),
    (0, "γ_prolTrans = 0.40  (feedback Plt → MPP)"),
    (0, "GFR patient de référence : 125 mL/min (S11)"),
    (1, "Cohérence avec le papier (annexe S11)"),
    (0, "IIV ajustée pour reproduire les G3/G4"),
    (1, "CV 25% sur Plt0, 20% sur Neut0"),
    (0, "Procédure : script calibrate_plt.R"),
], 0.5, 1.75, 5.7, 4.8)

rect(sl, 6.7, 1.25, 6.3, 5.7, fill=RGBColor(0xFF, 0xEB, 0xEB),
     line=ROUGE, line_w=Pt(1))
txbox(sl, "Difficultés rencontrées", 6.9, 1.3, 5.9, 0.4,
      size=13, bold=True, color=ROUGE)
bullet_block(sl, [
    (0, "GFR de référence non explicite dans Table 1"),
    (1, "Trouvé dans le Supplementary S11 (125 mL/min)"),
    (0, "Ambiguïté courbe déterministe vs. grade G4"),
    (1, "Courbe = patient médian ≠ données observées"),
    (1, "G4 vient de la queue de la distribution pop."),
    (0, "Scaling allométrique des paramètres"),
    (1, "MTT, baselines, taux de circulation"),
    (1, "Références multiples (Table 1 + annexes)"),
    (0, "δ_Neut absent du modèle (confirmé papier)"),
    (1, "Tentative initiale erronée → corrigée"),
    (0, "Convergence numérique (stiff ODEs)"),
    (1, "LSODA + rtol/atol ajustés"),
], 6.9, 1.75, 5.9, 4.8)

footer(sl, 7, N_SLIDES)

# ════════════════════════════════════════════════════════
# SLIDE 8 — Bilan & Prochaines étapes
# ════════════════════════════════════════════════════════
sl = add_slide()
rect(sl, 0, 0, 13.33, 7.5, fill=BLANC)
header_bar(sl, "Bilan & Prochaines étapes", "")

# Bilan
rect(sl, 0.3, 1.25, 12.4, 2.8, fill=RGBColor(0xE2, 0xEF, 0xDA), line=VERT, line_w=Pt(1.5))
txbox(sl, "Ce qui est accompli", 0.5, 1.3, 4, 0.4, size=14, bold=True, color=VERT)

done = [
    (0, "Modèle PK/PD complet implémenté en R  (pkpd_model_FORNARI.R)"),
    (0, "Validation rat : Figure 3 + Figure S1 reproduites"),
    (0, "Simulation humaine déterministe : Figure 4 Neut & Plt"),
    (0, "VPC population (1000 patients) — AUC 4, 5, 6"),
    (0, "Grades NCI-CTCAE : barplots Figure 4c — AUC 4, 5, 6"),
    (0, "Calibration δ_Plt = 0.80 et γ_prolTrans = 0.40"),
]
bullet_block(sl, done, 0.5, 1.75, 12.0, 2.2, size=13, color=VERT)

# Prochaines étapes
rect(sl, 0.3, 4.2, 12.4, 2.8, fill=RGBColor(0xE8, 0xF1, 0xFA), line=BLEU_CLAIR, line_w=Pt(1.5))
txbox(sl, "Prochaines étapes envisagées", 0.5, 4.25, 6, 0.4, size=14, bold=True, color=BLEU_FONCE)

next_steps = [
    (0, "Optimisation bayésienne des paramètres IIV (NONMEM / nlmixr2)"),
    (0, "Intégration de données cliniques réelles (validation externe)"),
    (0, "Extension à d'autres schémas (Q28D, combinaisons)"),
    (0, "Outil de prédiction de grade par patient (Shiny app ou rapport automatisé)"),
]
bullet_block(sl, next_steps, 0.5, 4.7, 12.0, 2.2, size=13)

footer(sl, 8, N_SLIDES)

# ── Sauvegarde ───────────────────────────────────────────
out = "/home/user/Hema/Avancement_PK_PD_Carboplatine.pptx"
prs.save(out)
print(f"Fichier créé : {out}")
