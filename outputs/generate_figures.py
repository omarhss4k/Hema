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
# FIGURE 1 — Modèle PK/PD redesigné
# ═══════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(24, 14))
ax.set_xlim(0, 12); ax.set_ylim(0, 7)
ax.axis('off'); fig.patch.set_facecolor('white')

# ─── BLOC 1 : PK ─────────────────────────────────────────
ax.add_patch(FancyBboxPatch((0.15, 0.8), 3.1, 6.0, boxstyle="round,pad=0.1",
    facecolor=BLUE_L, edgecolor=BLUE_M, linewidth=4, zorder=1))
ax.text(1.7, 6.6, "Bloc 1", ha='center', fontsize=20, fontweight='bold', color=BLUE, zorder=3)
ax.text(1.7, 6.2, "MODÈLE PK", ha='center', fontsize=22, fontweight='bold', color=BLUE, zorder=3)
ax.text(1.7, 5.75, "ADC 2-compartiments", ha='center', fontsize=16, color=DARK, style='italic', zorder=3)

# Dose IV → V1
arr(ax, -0.15, 4.0, 0.35, 4.0, color=DARK, lw=4.5)
ax.text(-0.05, 4.55, "Dose IV", ha='center', fontsize=16, fontweight='bold', color=DARK)

# V1 et V2
rbox(ax, 0.35, 3.0, 1.2, 1.8, 'Central\nV₁', BLUE_M, BLUE, fs=19)
rbox(ax, 1.85, 3.0, 1.2, 1.8, 'Périph.\nV₂', BLUE_M, BLUE, fs=19)
ax.annotate('', xy=(1.85, 3.9), xytext=(1.55, 3.9),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=2.5), zorder=4)

# CL élimination
ax.text(1.7, 2.5, "↓  Élimination (CL)", ha='center', fontsize=16, color=DARK, zorder=3)
ax.annotate('', xy=(1.7, 1.5), xytext=(1.7, 2.3),
            arrowprops=dict(arrowstyle='->', color='#888', lw=3,
                            linestyle='dashed'), zorder=4)
ax.text(1.7, 1.2, "Clairance\nplasmatique", ha='center', fontsize=14,
        color='#666', style='italic', zorder=3)

# ─── FLÈCHE PK → DAMAGE ──────────────────────────────────
arr(ax, 3.25, 4.5, 3.75, 4.5, color=DARK, lw=4.5)
ax.text(3.5, 5.05, "C_ADC (µM)", ha='center', fontsize=16, fontweight='bold', color=DARK)

# ─── BLOC 2 : DAMAGE ─────────────────────────────────────
ax.add_patch(FancyBboxPatch((3.75, 2.8), 1.9, 3.5, boxstyle="round,pad=0.1",
    facecolor=ORG_L, edgecolor=ORG, linewidth=4, zorder=1))
ax.text(4.7, 6.15, "Bloc 2", ha='center', fontsize=18, fontweight='bold', color=ORG, zorder=3)
rbox(ax, 3.9, 3.3, 1.6, 2.3, 'DAMAGE', ORG, '#8B3210', fs=24)

# ─── FLÈCHE DAMAGE → BLOC 3 ──────────────────────────────
arr(ax, 5.65, 4.5, 6.1, 4.5, color=DARK, lw=4.5)

# ─── BLOC 3 : PD / HÉMATOPOÏÈSE ──────────────────────────
ax.add_patch(FancyBboxPatch((6.1, 0.3), 5.75, 6.5, boxstyle="round,pad=0.1",
    facecolor=GRN_L, edgecolor=GRN_M, linewidth=4, zorder=1))
ax.text(8.97, 6.65, "Bloc 3 — MODÈLE PD / HÉMATOPOÏÈSE",
        ha='center', fontsize=22, fontweight='bold', color=GRN, zorder=3)

# Niveau 1 : MPP (centré)
rbox(ax, 7.45, 5.1, 3.05, 1.0, 'MPP  (cellules souches)', GRN_M, GRN, fs=18)

# Niveau 2 : CMP (gauche) + MEP (droite)
rbox(ax, 6.4,  3.4, 1.9, 0.95, 'CMP', GRN_M, GRN, fs=22)
rbox(ax, 9.6,  3.4, 1.9, 0.95, 'MEP', GRN_M, GRN, fs=22)

# MPP → CMP et MPP → MEP avec label "différenciation"
arr(ax, 7.9, 5.1, 7.35, 4.35, color='white', lw=3)
arr(ax, 9.05, 5.1, 9.6, 4.35, color='white', lw=3)
ax.text(7.0, 4.85, "différenciation", ha='center', fontsize=12,
        color='#ccc', style='italic', rotation=55, zorder=3)
ax.text(10.0, 4.85, "différenciation", ha='center', fontsize=12,
        color='#ccc', style='italic', rotation=-55, zorder=3)

# Niveau 3 : cellules matures (avec label "maturation + transit")
rbox(ax, 6.2,  1.4, 2.0, 0.95, 'Neut. + Mono.', GRN_D, GRN, fs=16)
rbox(ax, 8.45, 1.4, 1.5, 0.95, 'RBC',            GRN_D, GRN, fs=22)
rbox(ax, 10.2, 1.4, 1.45,0.95, 'Plt',            GRN_D, GRN, fs=22)

arr(ax, 7.35, 3.4, 7.2, 2.35, color='white', lw=3)
arr(ax, 9.7,  3.4, 9.2, 2.35, color='white', lw=3)
arr(ax, 10.55,3.4,10.9, 2.35, color='white', lw=3)

ax.text(7.0, 2.95, "maturation\n(transit)", ha='center', fontsize=12,
        color='#ccc', style='italic', zorder=3)
ax.text(10.9, 2.95, "maturation\n(transit)", ha='center', fontsize=12,
        color='#ccc', style='italic', zorder=3)

# Feedbacks (rétroaction homéostatique)
darr(ax, 6.5, 2.35, 7.6, 5.45, color='white', lw=2.5, rad=0.4)
darr(ax, 11.4, 2.35, 10.3, 5.45, color='white', lw=2.5, rad=-0.4)
ax.text(6.25, 4.0, "rétroaction", ha='center', fontsize=12,
        color='#bbb', style='italic', rotation=80, zorder=3)

# Effets médicament (DAMAGE → MPP, CMP, MEP)
ax.annotate('', xy=(7.45, 5.6), xytext=(5.65, 5.6),
            arrowprops=dict(arrowstyle='->', color=ORG, lw=3.5,
                            connectionstyle='arc3,rad=-0.15'), zorder=6)
ax.annotate('', xy=(6.4, 3.9), xytext=(5.65, 4.2),
            arrowprops=dict(arrowstyle='->', color=ORG, lw=3,
                            connectionstyle='arc3,rad=0.1'), zorder=6)
ax.annotate('', xy=(9.6, 3.9), xytext=(5.65, 3.8),
            arrowprops=dict(arrowstyle='->', color=ORG, lw=3,
                            connectionstyle='arc3,rad=0.2'), zorder=6)
ax.text(5.88, 5.9, "Effet\nmédicament", ha='center', fontsize=14,
        fontweight='bold', color=ORG, zorder=6)

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
