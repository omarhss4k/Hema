#!/usr/bin/env python3
"""Generate PK/PD model diagrams — simplified, very large fonts."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

def rbox(ax, x, y, w, h, text, fc, ec, tc='white', fs=22, fw='bold', lw=3, zorder=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder))
    ax.text(x+w/2, y+h/2, text, ha='center', va='center',
            fontsize=fs, fontweight=fw, color=tc, zorder=zorder+1,
            multialignment='center', linespacing=1.4)

def arr(ax, x1, y1, x2, y2, color='#2C2C2C', lw=4, rad=0):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

def darr(ax, x1, y1, x2, y2, color='white', lw=3, rad=0.35):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle='dashed',
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

BLUE  = '#2E5FAB'; BLUE_L = '#D6E4F7'; BLUE_M = '#4472C4'
ORG   = '#C0522A'; ORG_L  = '#FDEBD6'
GRN   = '#375623'; GRN_M  = '#70AD47'; GRN_L  = '#E2EFDA'; GRN_D = '#1F5C1F'
DARK  = '#1A1A1A'

# ═══════════════════════════════════════════════════════
# FIGURE 1 — Modèle PK/PD
# Espace coordonnées réduit = polices grandes visuellement
# ═══════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(22, 13))
ax.set_xlim(0, 11); ax.set_ylim(0, 6.5)
ax.axis('off'); fig.patch.set_facecolor('white')

# ── BLOC 1 : PK ──────────────────────
ax.add_patch(FancyBboxPatch((0.1, 1.0), 2.8, 5.3, boxstyle="round,pad=0.1",
    facecolor=BLUE_L, edgecolor=BLUE_M, linewidth=4, zorder=1))
ax.text(1.5, 6.05, "Bloc 1 — MODÈLE PK", ha='center',
        fontsize=24, fontweight='bold', color=BLUE, zorder=3)
ax.text(1.5, 5.62, "ADC 2-compartiments", ha='center',
        fontsize=18, color=DARK, style='italic', zorder=3)

rbox(ax, 0.25, 3.5, 1.1, 1.6, 'Central\nV₁',      BLUE_M, BLUE, fs=20)
rbox(ax, 1.6,  3.5, 1.1, 1.6, 'Périph.\nV₂',      BLUE_M, BLUE, fs=20)
ax.annotate('', xy=(1.6, 4.3), xytext=(1.35, 4.3),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=3), zorder=4)

ax.text(1.5, 2.8, "Élimination linéaire (CL)", ha='center', fontsize=17, color=DARK, zorder=3)
ax.text(1.5, 2.35, "↓ clairance plasmatique",   ha='center', fontsize=15,
        color='#555', style='italic', zorder=3)

arr(ax, -0.2, 4.3, 0.25, 4.3, color=DARK, lw=4)
ax.text(-0.1, 4.85, "Dose\nIV", ha='center', fontsize=17, fontweight='bold', color=DARK)

# ── BLOC 2 : DAMAGE ──────────────────
ax.add_patch(FancyBboxPatch((3.2, 2.5), 2.0, 3.7, boxstyle="round,pad=0.1",
    facecolor=ORG_L, edgecolor=ORG, linewidth=4, zorder=1))
ax.text(4.2, 5.95, "Bloc 2\nDOMMAGE", ha='center',
        fontsize=22, fontweight='bold', color=ORG, zorder=3)
rbox(ax, 3.4, 3.0, 1.6, 2.0, 'DAMAGE', ORG, '#8B3210', fs=22)

arr(ax, 2.9, 4.3, 3.2, 4.1, color=DARK, lw=4)
ax.text(3.05, 4.85, "C_ADC\n(µM)", ha='center', fontsize=16,
        fontweight='bold', color=DARK)

# ── BLOC 3 : PD / HÉMATOPOÏÈSE ───────
ax.add_patch(FancyBboxPatch((5.5, 0.2), 5.35, 6.1, boxstyle="round,pad=0.1",
    facecolor=GRN_L, edgecolor=GRN_M, linewidth=4, zorder=1))
ax.text(8.17, 6.1, "Bloc 3 — MODÈLE PD / HÉMATOPOÏÈSE", ha='center',
        fontsize=24, fontweight='bold', color=GRN, zorder=3)

# Niveau 1 : MPP
rbox(ax, 6.8, 4.9, 2.7, 1.0, 'MPP  (cellules souches)', GRN_M, GRN, fs=19)

# Niveau 2 : CMP + MEP
rbox(ax, 5.7, 3.3, 1.8, 0.9, 'CMP', GRN_M, GRN, fs=22)
rbox(ax, 8.8, 3.3, 1.8, 0.9, 'MEP', GRN_M, GRN, fs=22)
arr(ax, 7.5, 4.9, 6.6, 4.2, color='white', lw=3)
arr(ax, 8.2, 4.9, 9.7, 4.2, color='white', lw=3)

# Niveau 3 : cellules matures
rbox(ax, 5.6, 1.5, 1.8, 0.9, 'Neut.\n+ Mono.', GRN_D, GRN, fs=17)
rbox(ax, 7.7, 1.5, 1.5, 0.9, 'RBC',             GRN_D, GRN, fs=22)
rbox(ax, 9.5, 1.5, 1.2, 0.9, 'Plt',             GRN_D, GRN, fs=22)
arr(ax, 6.6, 3.3, 6.5, 2.4, color='white', lw=3)
arr(ax, 9.7, 3.3, 8.5, 2.4, color='white', lw=3)
arr(ax, 9.9, 3.3, 10.1, 2.4, color='white', lw=3)

# Feedbacks
darr(ax, 5.8, 2.4, 7.0, 5.2, color='white', lw=2.5, rad=0.35)
darr(ax, 10.5, 2.4, 9.4, 5.2, color='white', lw=2.5, rad=-0.35)

# Effets médicament
arr(ax, 5.2, 5.5, 6.8, 5.2, color=ORG, lw=3.5, rad=-0.1)
ax.text(5.8, 5.0, "Effet\nmédicament", ha='center',
        fontsize=14, color=ORG, fontweight='bold', zorder=6)
arr(ax, 5.2, 3.9, 5.7, 3.7, color=ORG, lw=3)
arr(ax, 5.2, 3.6, 8.8, 3.7, color=ORG, lw=3, rad=0.12)
arr(ax, 5.2, 3.8, 5.5, 3.8, color=DARK, lw=4)

plt.tight_layout(pad=0.3)
plt.savefig('/home/user/Hema/outputs/figure_pkpd_model.png', dpi=200,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pkpd_model.png")

# ═══════════════════════════════════════════════════════
# FIGURE 2 — Pipeline 4 étapes
# ═══════════════════════════════════════════════════════
PURPLE = '#5B2C8D'; TEAL = '#1A6B6B'

# 4 boîtes horizontales — texte minimal, polices très grandes
fig2, ax2 = plt.subplots(figsize=(20, 10))
ax2.set_xlim(0, 10); ax2.set_ylim(0, 5)
ax2.axis('off'); fig2.patch.set_facecolor('white')

steps = [
    (BLUE_M,   "ÉTAPE 1",  "Modèle\nhématopoïétique",  "Rat & Humain\nCarboplatine\nFornari 2019"),
    ('#B5510D', "ÉTAPE 2",  "Hématotoxicité\nT-DXd",    "FDA BLA 761139\nN = 300 patients\nRat → NHP → Humain"),
    (TEAL,     "ÉTAPE 3",  "Hématotoxicité\ncomposé",   "NHP interne\n8 singes — confidentiel\n4 niveaux de dose"),
    (PURPLE,   "ÉTAPE 4",  "Prédiction\nclinique",      "Transposition\nNHP → Humain\nAllométrie PK"),
]

BW = 2.2
BH = 4.2
GAP = 0.15

for i, (color, num, title, data) in enumerate(steps):
    x0 = i * (BW + GAP) + 0.1

    # Boîte principale
    ax2.add_patch(FancyBboxPatch((x0, 0.3), BW, BH,
        boxstyle="round,pad=0.06", facecolor='#F6F6F6',
        edgecolor=color, linewidth=4, zorder=1))

    # Bandeau titre
    ax2.add_patch(FancyBboxPatch((x0, 0.3 + BH - 1.4), BW, 1.4,
        boxstyle="round,pad=0.04", facecolor=color,
        edgecolor=color, linewidth=0, zorder=2))

    # Numéro étape
    ax2.text(x0 + BW/2, 0.3 + BH - 0.35, num,
             ha='center', va='center', fontsize=22,
             fontweight='bold', color='white', zorder=3)
    # Titre
    ax2.text(x0 + BW/2, 0.3 + BH - 0.95, title,
             ha='center', va='center', fontsize=19,
             fontweight='bold', color='white', zorder=3,
             multialignment='center', linespacing=1.3)

    # Données (centré verticalement dans la zone restante)
    ax2.text(x0 + BW/2, 0.3 + (BH - 1.4) / 2, data,
             ha='center', va='center', fontsize=18,
             color='#1A1A1A', zorder=3,
             multialignment='center', linespacing=1.6)

# Flèches entre étapes
for i in range(3):
    x_mid = 0.1 + (i+1)*(BW+GAP) - GAP/2
    ax2.annotate('', xy=(x_mid + 0.01, 2.7), xytext=(x_mid - 0.01, 2.7),
                 arrowprops=dict(arrowstyle='->', color='#999999', lw=6), zorder=6)

plt.tight_layout(pad=0.2)
plt.savefig('/home/user/Hema/outputs/figure_pipeline_steps.png', dpi=180,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pipeline_steps.png")
