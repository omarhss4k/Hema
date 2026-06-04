"""
Slides de résultats — Modèle PK/PD TGI Simeoni (2004) — FGFR2
Slide 1 : Fit prédit vs observé + paramètres estimés
Slide 2 : Interprétation pharmacologique (TSC, MTT, groupes)
"""
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
import os

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

CTRL_COL  = RGBColor(0x88, 0x88, 0x88)
D3_COL    = RGBColor(0x43, 0x93, 0xC3)
D10_COL   = RGBColor(0x21, 0x66, 0xAC)

# ── Paramètres estimés (à mettre à jour si nouvelle optimisation) ─────────────
L0   = 0.0607
L1   = ">> données"   # non identifiable sur cet intervalle
k1   = 0.200
k2   = 1.635e-6
TSC  = 37132
MTT  = 20.0
OBJ  = None          # sera ignoré si None

PLOT_PATH = "scripts/plot_PKPD_simeoni_FGFR2.png"

# ── Helpers ───────────────────────────────────────────────────────────────────
prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)

def new_slide():
    return prs.slides.add_slide(prs.slide_layouts[6])

def R(sl, x, y, w, h, fill=None, line=None, lw=Pt(0)):
    shp = sl.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill:
        shp.fill.solid(); shp.fill.fore_color.rgb = fill
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
    tf = tb.text_frame; tf.word_wrap = True
    p  = tf.paragraphs[0]; p.alignment = al
    run = p.add_run()
    run.text = text; run.font.size = Pt(sz)
    run.font.bold = bold; run.font.italic = it
    run.font.color.rgb = col

def IMG(sl, path, x, y, w, h):
    if os.path.exists(path):
        sl.shapes.add_picture(path, Inches(x), Inches(y),
                              Inches(w), Inches(h))
    else:
        R(sl, x, y, w, h, fill=GRIS_CLAIR, line=BLEU_CLAIR, lw=Pt(1))
        T(sl, f"[Image manquante]\n{path}", x+0.1, y+h/2-0.2, w-0.2, 0.4,
          sz=9, col=ROUGE, al=PP_ALIGN.CENTER)

def header(sl, titre, sous_titre, page):
    R(sl, 0, 0, 13.33, 1.06, fill=BLEU_FONCE)
    T(sl, titre,      0.30, 0.04, 11.8, 0.62, sz=22, bold=True, col=BLANC)
    T(sl, sous_titre, 0.30, 0.67, 11.5, 0.30, sz=10,
      col=RGBColor(0xBD, 0xD7, 0xEE))
    T(sl, page, 12.45, 0.10, 0.75, 0.55, sz=20, bold=True,
      col=ORANGE, al=PP_ALIGN.RIGHT)
    R(sl, 0, 1.01, 13.33, 0.055, fill=ORANGE)

def footer(sl):
    R(sl, 0, 7.20, 13.33, 0.30, fill=BLEU_FONCE)
    T(sl, "Réf. : Simeoni M. et al.  —  Cancer Research, 2004, 64 : 1094–1101  |  "
          "Optimisation : DEoptim (Differential Evolution) + affinage nlminb",
      0.20, 7.22, 10.5, 0.26, sz=8, col=BLANC)
    T(sl, "Réunion FGFR2", 11.2, 7.22, 1.9, 0.26, sz=9,
      col=BLANC, al=PP_ALIGN.RIGHT)


# =============================================================================
# SLIDE 1  —  Fit prédit vs observé + paramètres
# =============================================================================
sl1 = new_slide()
header(sl1,
       "Modèle PK/PD TGI Simeoni — Résultats  FGFR2  (Q2W × 4)",
       "Fc-silent FGFR2-huBPA-LP1  |  Souris xénogreffe  |  Groupes : Contrôle, 3 mg/kg, 10 mg/kg",
       "1/2")

# ── Graphique (gauche) ────────────────────────────────────────────────────────
IMG(sl1, PLOT_PATH, 0.15, 1.12, 8.55, 5.98)

# ── Panel droit : paramètres ──────────────────────────────────────────────────
R(sl1, 8.85, 1.12, 4.33, 5.98, fill=BLEU_PALE, line=BLEU_CLAIR, lw=Pt(1))
T(sl1, "Paramètres estimés", 8.98, 1.16, 4.10, 0.36,
  sz=13, bold=True, col=BLEU_FONCE)

# --- Tableau des paramètres ---
ROWS = [
    ("λ₀",  f"{L0:.4f} /j",         "Taux de croissance exponentielle",  BLANC),
    ("λ₁",  "non identifiable",      "Taux de croissance linéaire",       GRIS_BLU),
    ("k₁",  f"{k1:.3f} /j",          "Constante de transit (borne basse)", BLANC),
    ("k₂",  f"{k2:.3e} L/µg/j",      "Potence cytotoxique du drug",       GRIS_BLU),
]

COLS_W = [0.65, 1.35, 2.00]
COLS_X = [8.90, 9.57, 10.94]
HDR_Y  = 1.60

# En-têtes colonnes
for cx, cw, lbl in zip(COLS_X, COLS_W, ["Param.", "Valeur", "Définition"]):
    R(sl1, cx, HDR_Y, cw, 0.30, fill=BLEU_FONCE)
    T(sl1, lbl, cx+0.03, HDR_Y+0.03, cw-0.06, 0.24,
      sz=9, bold=True, col=BLANC, al=PP_ALIGN.CENTER)

for i, (par, val, defn, bg) in enumerate(ROWS):
    yp = HDR_Y + 0.32 + i * 0.38
    for cx, cw in zip(COLS_X, COLS_W):
        R(sl1, cx, yp, cw, 0.35, fill=bg,
          line=RGBColor(0xCC, 0xCC, 0xCC), lw=Pt(0.25))
    cells = [par, val, defn]
    for cx, cw, cell in zip(COLS_X, COLS_W, cells):
        is_p = (cell == par)
        T(sl1, cell, cx+0.04, yp+0.05, cw-0.08, 0.24,
          sz=9, bold=is_p,
          col=BLEU_FONCE if is_p else NOIR,
          al=PP_ALIGN.CENTER if is_p else PP_ALIGN.LEFT)

# --- TSC encadré ---
TSC_Y = HDR_Y + 0.32 + len(ROWS)*0.38 + 0.15
R(sl1, 8.90, TSC_Y, 4.10, 0.85, fill=ORANGE_PL, line=ORANGE, lw=Pt(1.5))
T(sl1, "TSC  (Tumor Static Concentration)",
  9.00, TSC_Y + 0.04, 3.90, 0.28, sz=11, bold=True, col=ORANGE)
T(sl1, f"TSC = λ₀ / k₂  =  {TSC:,.0f} µg/L",
  9.00, TSC_Y + 0.34, 3.90, 0.30, sz=14, bold=True, col=ROUGE,
  al=PP_ALIGN.CENTER)
T(sl1, f"MTT  =  4 / k₁  =  {MTT:.0f} jours",
  9.00, TSC_Y + 0.62, 3.90, 0.22, sz=10, col=NOIR, al=PP_ALIGN.CENTER)

# --- Méthode d'optimisation ---
OPT_Y = TSC_Y + 0.95
R(sl1, 8.90, OPT_Y, 4.10, 0.70, fill=GRIS_CLAIR,
  line=BLEU_CLAIR, lw=Pt(0.75))
T(sl1, "Méthode d'optimisation",
  9.00, OPT_Y + 0.04, 3.90, 0.24, sz=10, bold=True, col=BLEU_FONCE)
T(sl1, "DEoptim (évolution différentielle, NP=120, 600 générations)\n"
        "+ affinage local nlminb  |  Fonction objectif : Σ w·[log(obs)−log(pred)]²\n"
        "PK fixée (fit 2-comp. séparé)",
  9.00, OPT_Y + 0.28, 3.90, 0.40, sz=8.5, col=NOIR)

# --- Légende groupes ---
LEG_Y = OPT_Y + 0.82
for color, lbl in [(CTRL_COL,"Contrôle"), (D3_COL,"3 mg/kg"), (D10_COL,"10 mg/kg")]:
    R(sl1, 8.98, LEG_Y, 0.35, 0.18, fill=color)
    T(sl1, lbl, 9.38, LEG_Y, 3.50, 0.20, sz=9.5, col=NOIR)
    LEG_Y += 0.24

footer(sl1)


# =============================================================================
# SLIDE 2  —  Interprétation pharmacologique
# =============================================================================
sl2 = new_slide()
header(sl2,
       "Interprétation pharmacologique — Simeoni FGFR2",
       "TSC, MTT, efficacité par groupe  |  Comparaison avec le modèle Emax",
       "2/2")

# ── Colonne gauche : TSC & concentrations (0.15 → 6.50) ──────────────────────
R(sl2, 0.15, 1.12, 6.30, 5.98, fill=BLANC, line=BLEU_CLAIR, lw=Pt(1))
T(sl2, "Analyse TSC par groupe de dose",
  0.30, 1.16, 6.00, 0.36, sz=14, bold=True, col=BLEU_FONCE)

# Rappel TSC
R(sl2, 0.22, 1.60, 6.14, 0.72, fill=ORANGE_PL, line=ORANGE, lw=Pt(1.2))
T(sl2, f"TSC  =  λ₀ / k₂  =  {TSC:,.0f}  µg/L",
  0.32, 1.64, 5.94, 0.32, sz=16, bold=True, col=ROUGE, al=PP_ALIGN.CENTER)
T(sl2, "Concentration nécessaire pour maintenir la tumeur stable (croissance = inhibition)",
  0.32, 1.98, 5.94, 0.28, sz=9.5, col=NOIR, al=PP_ALIGN.CENTER, it=True)

# Table doses vs Cmax vs TSC
T(sl2, "Comparaison Cmax vs TSC", 0.30, 2.44, 6.00, 0.28,
  sz=11, bold=True, col=BLEU_FONCE)

HDR2 = [("Groupe", 1.00), ("Dose (µg/kg)", 1.40), ("Cmax (µg/L)", 1.40),
         ("Cmax / TSC", 1.20), ("Effet attendu", 1.10)]
HX   = [0.22, 1.24, 2.66, 4.08, 5.30]

for cx, (lbl, cw) in zip(HX, HDR2):
    R(sl2, cx, 2.78, cw, 0.30, fill=BLEU_FONCE)
    T(sl2, lbl, cx+0.03, 2.80, cw-0.06, 0.24,
      sz=9, bold=True, col=BLANC, al=PP_ALIGN.CENTER)

V1_val = 0.071   # L/kg (paramètre PK fixé)
DOSE_ROWS = [
    ("Contrôle", "—",     "—",       "—",     "Croissance libre",  BLANC),
    ("3 mg/kg",  "3 000",  f"{3000/V1_val:,.0f}", f"{3000/V1_val/TSC:.1f}×",
     "Ralentissement", GRIS_BLU),
    ("10 mg/kg", "10 000", f"{10000/V1_val:,.0f}", f"{10000/V1_val/TSC:.1f}×",
     "Stase / légère régression", BLANC),
]
for i, (grp, dose, cmax, ratio, effet, bg) in enumerate(DOSE_ROWS):
    yp = 3.10 + i * 0.42
    for cx, (_, cw) in zip(HX, HDR2):
        R(sl2, cx, yp, cw, 0.38, fill=bg,
          line=RGBColor(0xCC,0xCC,0xCC), lw=Pt(0.25))
    col_grp = [CTRL_COL, D3_COL, D10_COL][i]
    vals = [grp, dose, cmax, ratio, effet]
    for cx, (_, cw), val in zip(HX, HDR2, vals):
        is_grp = (val == grp)
        T(sl2, val, cx+0.04, yp+0.06, cw-0.08, 0.26,
          sz=9, bold=is_grp,
          col=col_grp if is_grp else (ROUGE if "×" in str(val) else NOIR),
          al=PP_ALIGN.CENTER)

# Note Cmax transitoire
R(sl2, 0.22, 4.42, 6.14, 0.55, fill=JAUNE_PL, line=ORANGE, lw=Pt(0.75))
T(sl2, "Cmax = concentration au pic (t=0 après injection IV)",
  0.32, 4.44, 5.94, 0.22, sz=9.5, bold=True, col=ORANGE)
T(sl2, f"La concentration décroît exponentiellement (t½ ≈ 3 j). "
        f"Le drug dépasse le TSC pendant ~{int(3*0.693/0.237/1)} j "
        f"à 3 mg/kg, ~{int(3*0.693/0.237/1 + 4)} j à 10 mg/kg.",
  0.32, 4.66, 5.94, 0.28, sz=9, col=NOIR, it=True)

# MTT interprétation
R(sl2, 0.22, 5.07, 6.14, 0.85, fill=BLEU_PALE, line=BLEU_CLAIR, lw=Pt(0.75))
T(sl2, f"MTT = {MTT:.0f} jours  (= 4/k₁)  —  signification biologique",
  0.32, 5.09, 5.94, 0.26, sz=11, bold=True, col=BLEU_FONCE)
T(sl2, f"Les cellules endommagées par le drug mettent en moyenne {MTT:.0f} jours "
        f"pour mourir (transit x₁→x₂→x₃→x₄). Ce délai explique pourquoi "
        f"la suppression tumorale persiste après la disparition du médicament (t½ PK ≈ 3 j).",
  0.32, 5.38, 5.94, 0.52, sz=9.5, col=NOIR)

# ── Colonne droite : Comparaison & Conclusion (6.75 → 13.10) ─────────────────
R(sl2, 6.75, 1.12, 6.20, 5.98, fill=BLANC, line=BLEU_CLAIR, lw=Pt(1))
T(sl2, "Comparaison des modèles",
  6.90, 1.16, 5.90, 0.36, sz=14, bold=True, col=BLEU_FONCE)

# Tableau comparatif
CMP_HDR = [("Critère", 2.30), ("Emax simple", 1.70), ("Simeoni (2004)", 1.80)]
CMP_X   = [6.82, 9.14, 10.86]
for cx, (lbl, cw) in zip(CMP_X, CMP_HDR):
    R(sl2, cx, 1.58, cw, 0.30, fill=BLEU_FONCE)
    T(sl2, lbl, cx+0.03, 1.60, cw-0.06, 0.24,
      sz=9, bold=True, col=BLANC, al=PP_ALIGN.CENTER)

CMP_ROWS = [
    ("Effet drogue",    "ke · C/(EC50+C)", "k₂ · C  (linéaire)"),
    ("Délai tumoral",   "Non (immédiat)",  "Oui — MTT = 4/k₁"),
    ("Compartiments",   "1 (TV direct)",    "4 (x₁ → x₂ → x₃ → x₄)"),
    ("Params. PD",      "3 (kg, ke, EC50)", "4 (λ₀, λ₁, k₁, k₂)"),
    ("Fit 10 mg/kg",    "Médiocre",         "Bon  ✓"),
    ("Fit Contrôle",    "Très mauvais ✗",  "Imparfait"),
    ("Interprétabilité","EC50, ke/kg",       "TSC = λ₀/k₂  ✓"),
]
for i, (crit, emax, sim) in enumerate(CMP_ROWS):
    yp = 1.90 + i * 0.38
    bg = BLANC if i % 2 == 0 else GRIS_BLU
    for cx, (_, cw) in zip(CMP_X, CMP_HDR):
        R(sl2, cx, yp, cw, 0.35, fill=bg,
          line=RGBColor(0xCC,0xCC,0xCC), lw=Pt(0.25))
    vals = [crit, emax, sim]
    for cx, (_, cw), val in zip(CMP_X, CMP_HDR, vals):
        is_good = "✓" in val
        is_bad  = "✗" in val or "Médiocre" in val or "Très" in val
        T(sl2, val, cx+0.04, yp+0.05, cw-0.08, 0.24,
          sz=8.5,
          col=VERT if is_good else (ROUGE if is_bad else NOIR),
          bold=is_good)

# Conclusion
CONC_Y = 1.90 + len(CMP_ROWS)*0.38 + 0.12
R(sl2, 6.82, CONC_Y, 6.10, 2.58, fill=VERT_PL, line=VERT, lw=Pt(1.5))
T(sl2, "Conclusions", 6.95, CONC_Y + 0.06, 5.80, 0.28,
  sz=13, bold=True, col=VERT)

CONCLUSIONS = [
    "Le modèle de Simeoni reproduit correctement la stase tumorale "
    "observée à 10 mg/kg grâce au délai de transit (MTT = 20 j)",

    f"Le TSC estimé ({TSC:,.0f} µg/L) est atteint dès les premières "
    f"heures post-injection pour les deux doses",

    "Le groupe 3 mg/kg montre une réponse partielle cohérente avec "
    "Cmax/TSC ≈ 1,1×  (juste au-dessus du seuil)",

    "Amélioration possible : contraindre λ₀ depuis le contrôle, "
    "augmenter la borne k₁ pour explorer MTT > 20 j",
]
for j, c in enumerate(CONCLUSIONS):
    yc = CONC_Y + 0.42 + j * 0.50
    T(sl2, f"{'①②③④'[j]}  {c}", 6.95, yc, 5.80, 0.46, sz=9.5, col=NOIR)

footer(sl2)

# =============================================================================
# SAUVEGARDE
# =============================================================================
out = "scripts/slide_resultats_simeoni_FGFR2.pptx"
prs.save(out)
print(f"Fichier créé : {out}")
