#!/usr/bin/env python3
"""Generate PK/PD model diagrams for memoir."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe

# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────
def rbox(ax, x, y, w, h, text, fc, ec, tc='white', fs=9, fw='bold', lw=2, zorder=2):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.07",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fs, fontweight=fw, color=tc, zorder=zorder+1,
            multialignment='center', linespacing=1.35)

def arr(ax, x1, y1, x2, y2, color='#2C2C2C', lw=2, style='->', rad=0):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

def darr(ax, x1, y1, x2, y2, color='white', lw=1.5, rad=0.3):
    """Dashed feedback arrow."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle='dashed',
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

# ─────────────────────────────────────────────
# PALETTE
# ─────────────────────────────────────────────
BLUE      = '#2E5FAB'
BLUE_L    = '#D6E4F7'
BLUE_M    = '#4472C4'
ORG       = '#C0522A'
ORG_L     = '#FDEBD6'
GRN       = '#375623'
GRN_M     = '#70AD47'
GRN_L     = '#E2EFDA'
GRN_D     = '#1F5C1F'
DARK      = '#2C2C2C'
GREY      = '#666666'

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1 — Modèle PK/PD (3 blocs)
# ═══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(18, 11))
ax.set_xlim(0, 18); ax.set_ylim(0, 11)
ax.axis('off'); fig.patch.set_facecolor('white')

# ── BLOC 1 : PK ──────────────────────────────
ax.add_patch(FancyBboxPatch((0.2, 3.8), 5.0, 6.8, boxstyle="round,pad=0.12",
    facecolor=BLUE_L, edgecolor=BLUE_M, linewidth=2.5, zorder=1))
ax.text(2.7, 10.25, "Bloc 1 — MODÈLE PK", ha='center', fontsize=12, fontweight='bold', color=BLUE, zorder=3)
ax.text(2.7, 9.82, "ADC 2-compartiments", ha='center', fontsize=9, color=DARK, style='italic', zorder=3)

rbox(ax, 0.55, 7.0, 1.85, 2.0, 'Central\nV₁',        BLUE_M, BLUE, fs=10)
rbox(ax, 2.85, 7.0, 1.85, 2.0, 'Périphérique\nV₂',   BLUE_M, BLUE, fs=10)

# ↔ between V1-V2
ax.annotate('', xy=(2.85, 8.0), xytext=(2.4, 8.0),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=1.8), zorder=4)

ax.text(2.7, 6.2, "Élimination linéaire (CL)", ha='center', fontsize=8.5, color=DARK, zorder=3)
ax.text(2.7, 5.75, "↓ clairance plasmatique", ha='center', fontsize=8, color=GREY, style='italic', zorder=3)

# IV dose arrow
arr(ax, -0.1, 8.0, 0.55, 8.0, color=DARK, lw=2.5)
ax.text(-0.05, 8.5, "Dose\nIV", ha='center', fontsize=9, fontweight='bold', color=DARK)

# ── BLOC 2 : DAMAGE ──────────────────────────
ax.add_patch(FancyBboxPatch((5.8, 5.5), 2.6, 4.6, boxstyle="round,pad=0.12",
    facecolor=ORG_L, edgecolor=ORG, linewidth=2.5, zorder=1))
ax.text(7.1, 9.77, "Bloc 2 — DOMMAGE", ha='center', fontsize=11, fontweight='bold', color=ORG, zorder=3)
rbox(ax, 6.15, 6.4, 1.9, 2.1, 'DAMAGE', ORG, '#8B3210', fs=11)

# PK → DAMAGE
arr(ax, 5.2, 8.0, 5.8, 7.5, color=DARK, lw=2.5)
ax.text(5.5, 8.3, "C_ADC\n(µM)", ha='center', fontsize=8.5, fontweight='bold', color=DARK)

# ── BLOC 3 : PD / HÉMATOPOÏÈSE ───────────────
ax.add_patch(FancyBboxPatch((9.1, 1.5), 8.7, 9.2, boxstyle="round,pad=0.12",
    facecolor=GRN_L, edgecolor=GRN_M, linewidth=2.5, zorder=1))
ax.text(13.45, 10.4, "Bloc 3 — MODÈLE PD / HÉMATOPOÏÈSE", ha='center',
        fontsize=12, fontweight='bold', color=GRN, zorder=3)

# MPP
rbox(ax, 11.4, 8.7, 3.1, 0.95, 'MPP  (cellules souches)', GRN_M, GRN, fs=9)

# CMP + MEP
rbox(ax, 9.5,  6.65, 2.4, 0.9, 'CMP', GRN_M, GRN, fs=10)
rbox(ax, 13.8, 6.65, 2.4, 0.9, 'MEP', GRN_M, GRN, fs=10)

# MPP → CMP / MEP
arr(ax, 12.1, 8.7,  10.7, 7.55, color='white', lw=2)
arr(ax, 13.9, 8.7,  15.0, 7.55, color='white', lw=2)

# Transit / matures (myeloid)
rbox(ax, 9.5, 4.65, 2.4, 0.9, 'Neutrophiles\n+ Monocytes\n(transit)', GRN_M, GRN, fs=7.5, lw=1.5)
arr(ax, 10.7, 6.65, 10.7, 5.55, color='white', lw=2)

rbox(ax, 13.8, 4.65, 2.4, 0.9, 'Réticulocytes\n(transit)',            GRN_M, GRN, fs=8,   lw=1.5)
arr(ax, 15.0, 6.65, 15.0, 5.55, color='white', lw=2)

# Mature bottom row
rbox(ax, 9.3,  2.5, 2.4, 0.9, 'Neutrophiles\n+ Monocytes', GRN_D, GRN, fs=8, lw=1.5)
rbox(ax, 12.2, 2.5, 2.0, 0.9, 'GR (RBC)',                  GRN_D, GRN, fs=8.5, lw=1.5)
rbox(ax, 14.7, 2.5, 2.4, 0.9, 'Plaquettes (Plt)',           GRN_D, GRN, fs=8.5, lw=1.5)

arr(ax, 10.7, 4.65, 10.5, 3.4,  color='white', lw=2)
arr(ax, 15.0, 4.65, 13.2, 3.4,  color='white', lw=2)
arr(ax, 15.5, 4.65, 15.9, 3.4,  color='white', lw=2)

# Feedback dashed arrows (mature → MPP)
darr(ax, 9.8,  3.4, 11.5, 9.2,  color='white', lw=1.5, rad= 0.35)
darr(ax, 16.5, 3.4, 14.4, 9.2,  color='white', lw=1.5, rad=-0.35)

# Drug effect arrows (DAMAGE → progenitors)
arr(ax, 8.4, 8.6, 11.4, 8.9,   color=ORG, lw=2, rad=-0.15)
ax.text(10.1, 8.45, "Effet\nmédicament", ha='center', fontsize=7.5, color=ORG, fontweight='bold', zorder=6)

arr(ax, 8.4, 7.4, 9.5, 7.1,    color=ORG, lw=2)
arr(ax, 8.4, 7.0, 13.8, 7.1,   color=ORG, lw=2, rad=0.18)

# DAMAGE → bloc 3
arr(ax, 8.4, 7.2, 9.1, 7.2, color=DARK, lw=2.5)

# ── GLOSSAIRE (bas gauche) ────────────────────
ax.add_patch(FancyBboxPatch((0.2, 0.1), 8.7, 3.3, boxstyle="round,pad=0.1",
    facecolor='#F8F8F8', edgecolor='#AAAAAA', linewidth=1.5, zorder=1))

ax.text(4.55, 3.15, "GLOSSAIRE — Pharmacocinétique & Hématopoïèse",
        ha='center', fontsize=9, fontweight='bold', color=DARK, zorder=3)

pk_items = [
    ("Dose IV",     "Dose intraveineuse"),
    ("V₁ central",  "Volume central (compartiment sanguin)"),
    ("V₂ périph.",  "Volume périphérique (compartiment tissulaire)"),
    ("CL",          "Clairance (vitesse d'élimination)"),
    ("C_ADC (µM)",  "Concentration plasmatique en micromolaire"),
]
hema_items = [
    ("MPP", "Progéniteur multipotent (cellule souche auto-renouvelable)"),
    ("CMP", "Progéniteur myéloïde commun (précurseur granulocytes/monocytes)"),
    ("MEP", "Progéniteur mégakaryocyte-érythroïde (précurseur RBC/plaquettes)"),
    ("GR/RBC",  "Globules rouges (érythrocytes)"),
    ("Plt",     "Plaquettes (issues des mégakaryocytes)"),
]

ax.text(0.5,  2.78, "Pharmacocinétique (PK)", fontsize=8.5, fontweight='bold', color=BLUE, zorder=3)
ax.text(4.55, 2.78, "Hématopoïèse",           fontsize=8.5, fontweight='bold', color=GRN,  zorder=3)

for i, (t, d) in enumerate(pk_items):
    y = 2.42 - i * 0.44
    ax.text(0.5,  y, f"• {t}",  fontsize=7.5, color=DARK,    fontweight='bold', zorder=3)
    ax.text(1.85, y, d,         fontsize=7.5, color='#444444',                  zorder=3)

for i, (t, d) in enumerate(hema_items):
    y = 2.42 - i * 0.44
    ax.text(4.55, y, f"• {t}",  fontsize=7.5, color=DARK,    fontweight='bold', zorder=3)
    ax.text(5.55, y, d,         fontsize=7.5, color='#444444',                  zorder=3)

# Separator
ax.plot([4.3, 4.3], [0.25, 2.95], color='#CCCCCC', lw=1, zorder=3)

plt.tight_layout(pad=0.4)
plt.savefig('/home/user/Hema/outputs/figure_pkpd_model.png', dpi=150,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pkpd_model.png")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2 — Pipeline 4 étapes
# ═══════════════════════════════════════════════════════════════════════
PURPLE    = '#5B2C8D'
PURPLE_L  = '#E8DCEF'
TEAL      = '#1A6B6B'
TEAL_L    = '#D1ECEC'

fig2, ax2 = plt.subplots(figsize=(18, 9))
ax2.set_xlim(0, 18); ax2.set_ylim(0, 9)
ax2.axis('off'); fig2.patch.set_facecolor('white')

steps = [
    # (x, header_text, header_bg, bullet_lines, result_lines, badge_text, badge_color)
    (0.3,
     "ÉTAPE 1 — Modèle hématopoïétique\n(Fornari et al. 2019)",
     BLUE_M,
     ["Données & méthodes",
      "• Rat traité au carboplatine",
      "• Données digitalisées (Fornari 2019)",
      "• EDO 8-compartiments :",
      "  MPP → CMP/MEP → Neut/Mono/Ret/RBC/Plt",
      "• Compartiements de transit + feedbacks",
      "• Outils R : deSolve · ggplot2"],
     ["Résultat clé",
      "• RMSE < 15 % (cellules circulantes)",
      "• Implémentation validée"],
     "Modèle\nvalidé", BLUE),

    (4.7,
     "ÉTAPE 2 — Hématotoxicité T-DXd\n(Validation FDA BLA 761139)",
     '#B5510D',
     ["Données & méthodes",
      "• ADC 2-cpt PK + TMDD (modèle NHP)",
      "• Transposition croisée Rat → NHP → Humain",
      "• Slopes rat transposés via IC50-scaling",
      "• Sim. N = 300 patients virtuels",
      "• Outils R : optim() · CTCAE v5 · ggplot2"],
     ["Résultat clé",
      "• Neutropénie tout grade : 29 % prédit vs 29 % FDA",
      "• Thrombocytopénie G3-4 : 3 % vs 3 % FDA"],
     "Cadre\ntransféré", '#B5510D'),

    (9.1,
     "ÉTAPE 3 — Hématotoxicité composé\n(NHP interne, données confidentielles)",
     TEAL,
     ["Données & méthodes",
      "• 8 singes cynomolgus, 4 doses (4/13/26/39 mg/kg)",
      "• Même cadre PK/PD que T-DXd",
      "• Ajustement PK individuel (Nelder-Mead)",
      "• Calibration visuelle des Slopes PD",
      "• Outils R : rxode2 · optim() · PKNCA"],
     ["Résultat clé",
      "• Profils PK/PD individuels concordants",
      "• CL ≈ 1,8 mL/h/kg ; t½ ≈ 55 h",
      "• Dose-proportionnalité validée"],
     "Mise à l'échelle\ntranslationnelle", TEAL),

    (13.5,
     "ÉTAPE 4 — Prédiction clinique\n(Transposition NHP → Humain)",
     PURPLE,
     ["Données & méthodes",
      "• Allométrie PK : CL_hum = CL_NHP × (70/BW)^0,75",
      "• Slopes NHP → point de départ humain",
      "• Simulation hématotoxicité humaine préliminaire",
      "• NOAEL NHP comme contrainte",
      "• Outils R : rxode2 · ggplot2 · CTCAE v5"],
     ["Résultat clé",
      "• Dose FIH candidate : P(G≥3) < 10 %",
      "• Support dossier réglementaire IND/CTA"],
     "Prédiction\nclinique", PURPLE),
]

for x0, header, hbg, bullets, results, badge, badge_c in steps:
    W = 4.1

    # Outer box
    ax2.add_patch(FancyBboxPatch((x0, 0.6), W, 8.0, boxstyle="round,pad=0.1",
        facecolor='#F7F7F7', edgecolor=badge_c, linewidth=2, zorder=1))

    # Header band
    ax2.add_patch(FancyBboxPatch((x0, 6.8), W, 1.8, boxstyle="round,pad=0.05",
        facecolor=hbg, edgecolor=hbg, linewidth=0, zorder=2))
    ax2.text(x0 + W/2, 7.7, header, ha='center', va='center',
             fontsize=9, fontweight='bold', color='white', zorder=3,
             multialignment='center', linespacing=1.4)

    # Bullets
    for i, line in enumerate(bullets):
        y = 6.45 - i * 0.52
        fw = 'bold' if i == 0 else 'normal'
        fs = 8 if i == 0 else 7.5
        ax2.text(x0 + 0.15, y, line, fontsize=fs, fontweight=fw,
                 color='#333333', zorder=3)

    # Results block
    ax2.add_patch(FancyBboxPatch((x0 + 0.1, 0.75), W - 0.2, 2.2, boxstyle="round,pad=0.07",
        facecolor=badge_c, edgecolor=badge_c, linewidth=0, alpha=0.12, zorder=2))
    for i, line in enumerate(results):
        y = 2.75 - i * 0.48
        fw = 'bold' if i == 0 else 'normal'
        fc = badge_c if i == 0 else '#222222'
        ax2.text(x0 + 0.2, y, line, fontsize=7.5, fontweight=fw, color=fc, zorder=3)

# Arrows between steps
arrow_y = 7.7
for x_arr in [4.4, 8.8, 13.3]:
    ax2.annotate('', xy=(x_arr + 0.2, arrow_y), xytext=(x_arr - 0.15, arrow_y),
                 arrowprops=dict(arrowstyle='->', color='#888888', lw=3), zorder=6)

# Step labels at bottom
step_colors = [BLUE_M, '#B5510D', TEAL, PURPLE]
labels = ["Étape 1", "Étape 2", "Étape 3", "Étape 4"]
label_xs = [2.35, 6.75, 11.15, 15.55]
for lx, lbl, lc in zip(label_xs, labels, step_colors):
    rbox(ax2, lx - 0.7, 0.1, 1.4, 0.42, lbl, lc, lc, fs=8.5, lw=1.5, zorder=3)

plt.tight_layout(pad=0.4)
plt.savefig('/home/user/Hema/outputs/figure_pipeline_steps.png', dpi=150,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pipeline_steps.png")
