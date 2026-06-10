#!/usr/bin/env python3
"""Generate PK/PD model diagrams for memoir — very large fonts."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib as mpl
mpl.rcParams['font.size'] = 14

def rbox(ax, x, y, w, h, text, fc, ec, tc='white', fs=16, fw='bold', lw=2.5, zorder=2):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fs, fontweight=fw, color=tc, zorder=zorder+1,
            multialignment='center', linespacing=1.35)

def arr(ax, x1, y1, x2, y2, color='#2C2C2C', lw=3, style='->', rad=0):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

def darr(ax, x1, y1, x2, y2, color='white', lw=2.5, rad=0.3):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle='dashed',
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

BLUE  = '#2E5FAB'; BLUE_L = '#D6E4F7'; BLUE_M = '#4472C4'
ORG   = '#C0522A'; ORG_L  = '#FDEBD6'
GRN   = '#375623'; GRN_M  = '#70AD47'; GRN_L  = '#E2EFDA'; GRN_D = '#1F5C1F'
DARK  = '#2C2C2C'; GREY   = '#666666'

# ═══════════════════════════════════════════════════════
# FIGURE 1 — Modèle PK/PD (3 blocs)
# Coordonnées réduites → polices plus grandes visuellement
# ═══════════════════════════════════════════════════════
W, H = 20, 12
fig, ax = plt.subplots(figsize=(26, 16))
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.axis('off'); fig.patch.set_facecolor('white')

# ── BLOC 1 : PK ──────────────────────────────────────
ax.add_patch(FancyBboxPatch((0.2, 3.5), 5.5, 8.2, boxstyle="round,pad=0.12",
    facecolor=BLUE_L, edgecolor=BLUE_M, linewidth=3.5, zorder=1))
ax.text(2.95, 11.35, "Bloc 1 — MODÈLE PK",
        ha='center', fontsize=20, fontweight='bold', color=BLUE, zorder=3)
ax.text(2.95, 10.8, "ADC 2-compartiments",
        ha='center', fontsize=15, color=DARK, style='italic', zorder=3)

rbox(ax, 0.5,  7.2, 2.2, 2.6, 'Central\nV₁',       BLUE_M, BLUE, fs=16)
rbox(ax, 3.2,  7.2, 2.2, 2.6, 'Périphérique\nV₂',  BLUE_M, BLUE, fs=16)
ax.annotate('', xy=(3.2, 8.5), xytext=(2.7, 8.5),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=2.5), zorder=4)

ax.text(2.95, 6.3, "Élimination linéaire (CL)",
        ha='center', fontsize=14, color=DARK, zorder=3)
ax.text(2.95, 5.7, "↓ clairance plasmatique",
        ha='center', fontsize=13, color=GREY, style='italic', zorder=3)

arr(ax, -0.3, 8.5, 0.5, 8.5, color=DARK, lw=3.5)
ax.text(-0.15, 9.3, "Dose\nIV", ha='center', fontsize=14,
        fontweight='bold', color=DARK)

# ── BLOC 2 : DAMAGE ──────────────────────────────────
ax.add_patch(FancyBboxPatch((6.3, 5.8), 3.0, 5.4, boxstyle="round,pad=0.12",
    facecolor=ORG_L, edgecolor=ORG, linewidth=3.5, zorder=1))
ax.text(7.8, 10.9, "Bloc 2 —\nDOMMAGE",
        ha='center', fontsize=18, fontweight='bold', color=ORG, zorder=3)
rbox(ax, 6.6, 6.8, 2.4, 2.8, 'DAMAGE', ORG, '#8B3210', fs=18)

arr(ax, 5.7, 8.5, 6.3, 8.0, color=DARK, lw=3.5)
ax.text(6.0, 9.3, "C_ADC\n(µM)", ha='center', fontsize=13,
        fontweight='bold', color=DARK)

# ── BLOC 3 : PD / HÉMATOPOÏÈSE ───────────────────────
ax.add_patch(FancyBboxPatch((9.8, 1.5), 9.95, 10.2, boxstyle="round,pad=0.12",
    facecolor=GRN_L, edgecolor=GRN_M, linewidth=3.5, zorder=1))
ax.text(14.75, 11.4, "Bloc 3 — MODÈLE PD / HÉMATOPOÏÈSE",
        ha='center', fontsize=20, fontweight='bold', color=GRN, zorder=3)

rbox(ax, 13.1, 9.6,  3.3, 1.1, 'MPP  (cellules souches)', GRN_M, GRN, fs=14)
rbox(ax, 10.2, 7.5,  2.8, 1.0, 'CMP', GRN_M, GRN, fs=16)
rbox(ax, 16.5, 7.5,  2.8, 1.0, 'MEP', GRN_M, GRN, fs=16)

arr(ax, 14.1, 9.6, 11.6,  8.5, color='white', lw=3)
arr(ax, 15.4, 9.6, 17.9,  8.5, color='white', lw=3)

rbox(ax, 10.2, 5.4, 2.8, 1.0, 'Neutrophiles\n+ Monocytes (transit)', GRN_M, GRN, fs=12, lw=2)
arr(ax, 11.6, 7.5, 11.6, 6.4, color='white', lw=3)

rbox(ax, 16.5, 5.4, 2.8, 1.0, 'Réticulocytes (transit)', GRN_M, GRN, fs=13, lw=2)
arr(ax, 17.9, 7.5, 17.9, 6.4, color='white', lw=3)

rbox(ax, 10.0, 2.8, 2.8, 1.0, 'Neutrophiles\n+ Monocytes', GRN_D, GRN, fs=13, lw=2)
rbox(ax, 13.4, 2.8, 2.4, 1.0, 'GR (RBC)',                  GRN_D, GRN, fs=14, lw=2)
rbox(ax, 16.6, 2.8, 2.8, 1.0, 'Plaquettes (Plt)',           GRN_D, GRN, fs=13, lw=2)

arr(ax, 11.6, 5.4, 11.4, 3.8, color='white', lw=3)
arr(ax, 17.9, 5.4, 16.0, 3.8, color='white', lw=3)
arr(ax, 18.4, 5.4, 18.8, 3.8, color='white', lw=3)

darr(ax, 10.4, 3.8, 13.2, 10.1, color='white', lw=2, rad= 0.35)
darr(ax, 19.2, 3.8, 16.3, 10.1, color='white', lw=2, rad=-0.35)

arr(ax, 9.3, 9.8, 13.1, 9.9,   color=ORG, lw=3, rad=-0.1)
ax.text(11.1, 9.4, "Effet\nmédicament", ha='center',
        fontsize=12, color=ORG, fontweight='bold', zorder=6)
arr(ax, 9.3, 8.2, 10.2, 8.0,   color=ORG, lw=3)
arr(ax, 9.3, 7.8, 16.5, 8.0,   color=ORG, lw=3, rad=0.15)
arr(ax, 9.3, 8.0, 9.8,  8.0,   color=DARK, lw=3.5)

# ── GLOSSAIRE ────────────────────────────────────────
ax.add_patch(FancyBboxPatch((0.2, 0.15), 9.3, 3.0, boxstyle="round,pad=0.1",
    facecolor='#F8F8F8', edgecolor='#AAAAAA', linewidth=2, zorder=1))
ax.text(4.85, 2.95, "GLOSSAIRE — Pharmacocinétique & Hématopoïèse",
        ha='center', fontsize=14, fontweight='bold', color=DARK, zorder=3)

pk_items = [
    ("Dose IV",    "Dose intraveineuse"),
    ("V₁ central", "Volume central (sang)"),
    ("V₂ périph.", "Volume périphérique (tissus)"),
    ("CL",         "Clairance — élimination"),
    ("C_ADC (µM)", "Concentration plasmatique"),
]
hema_items = [
    ("MPP",    "Progéniteur multipotent (cellule souche)"),
    ("CMP",    "Progéniteur myéloïde commun"),
    ("MEP",    "Progéniteur mégakaryocyte-érythroïde"),
    ("GR/RBC", "Globules rouges"),
    ("Plt",    "Plaquettes"),
]

ax.text(0.5,  2.58, "Pharmacocinétique (PK)", fontsize=13, fontweight='bold', color=BLUE, zorder=3)
ax.text(4.9,  2.58, "Hématopoïèse",           fontsize=13, fontweight='bold', color=GRN,  zorder=3)
ax.plot([4.65, 4.65], [0.3, 2.75], color='#CCCCCC', lw=1.5, zorder=3)

for i, (t, d) in enumerate(pk_items):
    y = 2.25 - i * 0.43
    ax.text(0.5,  y, f"• {t}",  fontsize=12, color=DARK,    fontweight='bold', zorder=3)
    ax.text(2.1,  y, d,         fontsize=12, color='#444444',                  zorder=3)

for i, (t, d) in enumerate(hema_items):
    y = 2.25 - i * 0.43
    ax.text(4.9,  y, f"• {t}",  fontsize=12, color=DARK,    fontweight='bold', zorder=3)
    ax.text(6.1,  y, d,         fontsize=12, color='#444444',                  zorder=3)

plt.tight_layout(pad=0.3)
plt.savefig('/home/user/Hema/outputs/figure_pkpd_model.png', dpi=180,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pkpd_model.png")

# ═══════════════════════════════════════════════════════
# FIGURE 2 — Pipeline 4 étapes
# ═══════════════════════════════════════════════════════
PURPLE = '#5B2C8D'; TEAL = '#1A6B6B'

fig2, ax2 = plt.subplots(figsize=(28, 14))
ax2.set_xlim(0, 20); ax2.set_ylim(0, 10)
ax2.axis('off'); fig2.patch.set_facecolor('white')

steps = [
    (0.2,
     "ÉTAPE 1\nModèle hématopoïétique\n(Fornari et al. 2019)",
     BLUE_M,
     ["Données & méthodes",
      "• Rat / carboplatine",
      "• Données digitalisées (Fornari 2019)",
      "• EDO 8 compartiments :",
      "  MPP→CMP/MEP→Neut/Mono/Ret/RBC/Plt",
      "• Transit + feedbacks",
      "• Outils R : deSolve · ggplot2"],
     ["Résultat clé",
      "• RMSE < 15 % (cellules circulantes)",
      "• Implémentation validée"],
     BLUE),

    (5.2,
     "ÉTAPE 2\nHématotoxicité T-DXd\n(Validation FDA BLA 761139)",
     '#B5510D',
     ["Données & méthodes",
      "• ADC 2-cpt PK + TMDD (NHP)",
      "• Transposition Rat → NHP → Humain",
      "• Slopes via IC50-scaling",
      "• N = 300 patients virtuels",
      "• Outils R : optim() · CTCAE v5"],
     ["Résultat clé",
      "• Neutropénie tout grade : 29 % vs 29 % FDA",
      "• Thrombocytopénie G3-4 : 3 % vs 3 % FDA"],
     '#B5510D'),

    (10.2,
     "ÉTAPE 3\nHématotoxicité composé\n(NHP interne — confidentiel)",
     TEAL,
     ["Données & méthodes",
      "• 8 singes cynomolgus",
      "• 4 doses : 4 / 13 / 26 / 39 mg/kg",
      "• Même cadre PK/PD que T-DXd",
      "• Ajustement PK individuel (Nelder-Mead)",
      "• Calibration visuelle des Slopes PD",
      "• Outils R : rxode2 · optim() · PKNCA"],
     ["Résultat clé",
      "• Profils PK/PD concordants",
      "• CL ≈ 1,8 mL/h/kg ;  t½ ≈ 55 h",
      "• Dose-proportionnalité validée"],
     TEAL),

    (15.2,
     "ÉTAPE 4\nPrédiction clinique\n(Transposition NHP → Humain)",
     PURPLE,
     ["Données & méthodes",
      "• Allométrie PK :",
      "  CL_hum = CL_NHP × (70/BW)^0,75",
      "• Slopes NHP → point de départ humain",
      "• Simulation hématotoxicité humaine",
      "• NOAEL NHP comme contrainte",
      "• Outils R : rxode2 · ggplot2 · CTCAE v5"],
     ["Résultat clé",
      "• Dose FIH : P(G≥3) < 10 %",
      "• Support dossier IND/CTA"],
     PURPLE),
]

W2 = 4.7
for x0, header, hbg, bullets, results, badge_c in steps:
    ax2.add_patch(FancyBboxPatch((x0, 0.5), W2, 9.2, boxstyle="round,pad=0.1",
        facecolor='#F7F7F7', edgecolor=badge_c, linewidth=3, zorder=1))

    ax2.add_patch(FancyBboxPatch((x0, 7.2), W2, 2.5, boxstyle="round,pad=0.05",
        facecolor=hbg, edgecolor=hbg, linewidth=0, zorder=2))
    ax2.text(x0 + W2/2, 8.45, header, ha='center', va='center',
             fontsize=13, fontweight='bold', color='white', zorder=3,
             multialignment='center', linespacing=1.5)

    for i, line in enumerate(bullets):
        y = 6.9 - i * 0.65
        fw = 'bold' if i == 0 else 'normal'
        fs = 13 if i == 0 else 11.5
        ax2.text(x0 + 0.18, y, line, fontsize=fs, fontweight=fw,
                 color='#333333', zorder=3)

    ax2.add_patch(FancyBboxPatch((x0 + 0.12, 0.62), W2 - 0.24, 2.5,
        boxstyle="round,pad=0.07", facecolor=badge_c, edgecolor=badge_c,
        linewidth=0, alpha=0.13, zorder=2))
    for i, line in enumerate(results):
        y = 2.95 - i * 0.62
        fw = 'bold' if i == 0 else 'normal'
        fc = badge_c if i == 0 else '#222222'
        ax2.text(x0 + 0.25, y, line, fontsize=11.5, fontweight=fw, color=fc, zorder=3)

# Arrows between steps
for x_arr in [4.9, 9.9, 14.9]:
    ax2.annotate('', xy=(x_arr + 0.3, 8.45), xytext=(x_arr - 0.0, 8.45),
                 arrowprops=dict(arrowstyle='->', color='#888888', lw=5), zorder=6)

# Step labels
for x0, lbl, lc in zip([0.2, 5.2, 10.2, 15.2],
                        ["Étape 1", "Étape 2", "Étape 3", "Étape 4"],
                        [BLUE_M, '#B5510D', TEAL, PURPLE]):
    rbox(ax2, x0 + 0.85, 0.05, 3.0, 0.42, lbl, lc, lc, fs=13, lw=2, zorder=3)

plt.tight_layout(pad=0.3)
plt.savefig('/home/user/Hema/outputs/figure_pipeline_steps.png', dpi=180,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pipeline_steps.png")
