"""
Figure synthétique du pipeline PK/PD — 4 étapes
Export PNG 300 dpi, largeur 18 cm (A4 portrait)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ── Dimensions ────────────────────────────────────────────
FIG_W_CM = 18.0
FIG_H_CM = 22.0
DPI = 300
fig, ax = plt.subplots(figsize=(FIG_W_CM / 2.54, FIG_H_CM / 2.54))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# ── Palette ───────────────────────────────────────────────
BLEU      = "#2C3E50"
BLEU_CLAR = "#4A6FA5"
VERT      = "#1A6B3C"
ORANGE    = "#C0392B"
VIOLET    = "#6C3483"
GRIS_FOND = "#F4F6F7"
BLANC     = "#FFFFFF"
GRIS_TXT  = "#555555"

COULEURS = [BLEU_CLAR, VERT, ORANGE, VIOLET]

# ── Données des étapes ────────────────────────────────────
etapes = [
    {
        "num": "Étape 1",
        "titre": "Validation du cadre",
        "soustitre": "Rat / Carboplatine",
        "entree": [
            "Données : rat traité au carboplatine",
            "            (digitalisées — Fornari 2019)",
            "Modèle : Fornari 2019 — 8 compartiments ODE",
            "            MPP → CMP/MEP → Neut/Plt/Ret/RBC",
        ],
        "output": "RMSE < 15 % sur cellules circulantes\nValidation de l'implementation OK",
        "outils": "deSolve · rxode2 · ggplot2",
    },
    {
        "num": "Étape 2",
        "titre": "Preuve de concept humain",
        "soustitre": "T-DXd (trastuzumab déruxtécan)",
        "entree": [
            "PK : modèle 2-compartiments, paramètres FDA",
            "        (BLA 761139 — DESTINY-Breast01)",
            "PD : Slopes rat → transposés humain",
            "Simulation : N = 300 patients virtuels",
        ],
        "output": "Neutropénie tout grade : 29 % prédit vs 29 % FDA\nThrombocytopénie G3-4 : 3 % vs 3 % FDA (OK)",
        "outils": "optim() · CTCAE v5 · ggplot2",
    },
    {
        "num": "Étape 3",
        "titre": "Application NHP — composé interne",
        "soustitre": "Inhibiteur FGFR2 (confidentiel)",
        "entree": [
            "8 singes cynomolgus, 4 doses (4/13/26/39 mg/kg)",
            "PK : ajustement individuel Nelder-Mead",
            "        SSR log-scale, NCA comme initialisation",
            "PD : calibration visuelle itérative des Slopes",
        ],
        "output": "Profils PK/PD individuels concordants\nCL ~ 1,8 mL/h/kg · t½β ~ 55 h · dose-prop. (OK)",
        "outils": "rxode2 · optim() · PKNCA · ggplot2",
    },
    {
        "num": "Étape 4",
        "titre": "Perspectives cliniques (FIH)",
        "soustitre": "Transposition NHP → Humain",
        "entree": [
            "Allométrie PK : CL_hum = CL_NHP × (70/BW)⁰·⁷⁵",
            "Slopes NHP → point de départ humain",
            "Simulation hématotoxicité humaine préliminaire",
            "Comparaison NOAEL NHP comme contrainte",
        ],
        "output": "Dose FIH candidate : P(G≥3) < 10 %\nSupport dossier réglementaire IND/CTA",
        "outils": "rxode2 · ggplot2 · CTCAE v5",
    },
]

# ── Layout vertical ───────────────────────────────────────
BOX_X      = 0.04
BOX_W      = 0.92
BOX_H      = 0.175
GAP        = 0.038
ARROW_H    = 0.022
TOTAL_H    = len(etapes) * BOX_H + (len(etapes) - 1) * (GAP + ARROW_H)
Y_START    = 0.97

def draw_box(ax, etape, y_top, color):
    # Boîte principale
    box = FancyBboxPatch((BOX_X, y_top - BOX_H), BOX_W, BOX_H,
                         boxstyle="round,pad=0.008",
                         linewidth=1.5, edgecolor=color,
                         facecolor=GRIS_FOND)
    ax.add_patch(box)

    # Bande de titre colorée
    titre_h = 0.042
    titre_box = FancyBboxPatch((BOX_X, y_top - titre_h), BOX_W, titre_h,
                               boxstyle="round,pad=0.008",
                               linewidth=0, edgecolor=color,
                               facecolor=color)
    ax.add_patch(titre_box)

    # Numéro + titre
    ax.text(BOX_X + 0.012, y_top - titre_h / 2,
            etape["num"],
            ha="left", va="center", fontsize=8.5, fontweight="bold",
            color=BLANC, fontfamily="DejaVu Serif")
    ax.text(BOX_X + BOX_W / 2, y_top - titre_h / 2,
            f"{etape['titre']}  —  {etape['soustitre']}",
            ha="center", va="center", fontsize=8.5, fontweight="bold",
            color=BLANC, fontfamily="DejaVu Serif")

    # Corps : entrées (gauche) + output (droite)
    body_y_top = y_top - titre_h
    body_h     = BOX_H - titre_h
    mid_x      = BOX_X + BOX_W * 0.56
    sep_y_top  = body_y_top - 0.005
    sep_y_bot  = y_top - BOX_H + 0.012

    # Séparateur vertical
    ax.plot([mid_x, mid_x], [sep_y_bot, sep_y_top],
            color=color, lw=0.7, alpha=0.5)

    # En-tête colonnes
    col_y = body_y_top - 0.013
    ax.text(BOX_X + 0.022, col_y, "Données & méthodes",
            ha="left", va="top", fontsize=6.8, fontweight="bold",
            color=color, fontfamily="DejaVu Serif")
    ax.text(mid_x + 0.018, col_y, "Résultat clé",
            ha="left", va="top", fontsize=6.8, fontweight="bold",
            color=color, fontfamily="DejaVu Serif")

    # Entrées
    for i, line in enumerate(etape["entree"]):
        ax.text(BOX_X + 0.022, col_y - 0.018 - i * 0.022,
                line,
                ha="left", va="top", fontsize=6.2,
                color=GRIS_TXT, fontfamily="DejaVu Serif")

    # Output
    ax.text(mid_x + 0.018, col_y - 0.018,
            etape["output"],
            ha="left", va="top", fontsize=6.2,
            color="#1A1A1A", fontfamily="DejaVu Serif",
            linespacing=1.4)

    # Bande outils
    outil_y = y_top - BOX_H + 0.001
    ax.text(BOX_X + BOX_W / 2, outil_y + 0.01,
            f"Outils R :  {etape['outils']}",
            ha="center", va="bottom", fontsize=5.8, style="italic",
            color=color, fontfamily="DejaVu Serif")

def draw_arrow(ax, y_from, color_top, color_bot):
    mid_color = "#888888"
    ax.annotate("",
                xy=(0.5, y_from - ARROW_H),
                xytext=(0.5, y_from),
                arrowprops=dict(arrowstyle="-|>", color=mid_color,
                                lw=1.8, mutation_scale=14))

# ── Titre général ─────────────────────────────────────────
ax.text(0.5, 0.99,
        "Pipeline de modélisation PK/PD semi-mécaniste",
        ha="center", va="top", fontsize=11, fontweight="bold",
        color=BLEU, fontfamily="DejaVu Serif")
ax.text(0.5, 0.963,
        "Validation → Preuve de concept → Application NHP → Perspectives cliniques",
        ha="center", va="top", fontsize=7.5, style="italic",
        color=GRIS_TXT, fontfamily="DejaVu Serif")

# ── Dessin des boîtes et flèches ──────────────────────────
y = Y_START - 0.045
for i, etape in enumerate(etapes):
    draw_box(ax, etape, y, COULEURS[i])
    if i < len(etapes) - 1:
        y_arrow_from = y - BOX_H
        draw_arrow(ax, y_arrow_from, COULEURS[i], COULEURS[i + 1])
        y = y - BOX_H - GAP - ARROW_H

# ── Export ────────────────────────────────────────────────
out = "outputs/figure_pipeline_pkpd.png"
fig.savefig(out, dpi=DPI, bbox_inches="tight",
            facecolor="white", edgecolor="none")
print(f"✓ Figure exportée : {out}  ({DPI} dpi)")
print(f"  Dimensions : {FIG_W_CM} × {FIG_H_CM} cm")
