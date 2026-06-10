#!/usr/bin/env python3
"""Generate PK/PD model diagrams for memoir — large fonts."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

def rbox(ax, x, y, w, h, text, fc, ec, tc='white', fs=12, fw='bold', lw=2, zorder=2):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.07",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fs, fontweight=fw, color=tc, zorder=zorder+1,
            multialignment='center', linespacing=1.35)

def arr(ax, x1, y1, x2, y2, color='#2C2C2C', lw=2.5, style='->', rad=0):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

def darr(ax, x1, y1, x2, y2, color='white', lw=2, rad=0.3):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle='dashed',
                                connectionstyle=f'arc3,rad={rad}'), zorder=5)

BLUE   = '#2E5FAB'; BLUE_L  = '#D6E4F7'; BLUE_M  = '#4472C4'
ORG    = '#C0522A'; ORG_L   = '#FDEBD6'
GRN    = '#375623'; GRN_M   = '#70AD47'; GRN_L   = '#E2EFDA'; GRN_D = '#1F5C1F'
DARK   = '#2C2C2C'; GREY    = '#666666'

# ═══════════════════════════════════════════════════════
# FIGURE 1 — Modèle PK/PD (3 blocs)
# ═══════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(24, 14))
ax.set_xlim(0, 24); ax.set_ylim(0, 14)
ax.axis('off'); fig.patch.set_facecolor('white')

# ── BLOC 1 : PK ─────────────────────────────────────────
ax.add_patch(FancyBboxPatch((0.3, 4.5), 6.5, 9.0, boxstyle="round,pad=0.15",
    facecolor=BLUE_L, edgecolor=BLUE_M, linewidth=3, zorder=1))
ax.text(3.55, 13.1, "Bloc 1 — MODÈLE PK", ha='center',
        fontsize=16, fontweight='bold', color=BLUE, zorder=3)
ax.text(3.55, 12.5, "ADC 2-compartiments", ha='center',
        fontsize=12, color=DARK, style='italic', zorder=3)

rbox(ax, 0.7,  9.0, 2.4, 2.8, 'Central\nV₁',       BLUE_M, BLUE, fs=13)
rbox(ax, 3.7,  9.0, 2.4, 2.8, 'Périphérique\nV₂',  BLUE_M, BLUE, fs=13)

ax.annotate('', xy=(3.7, 10.4), xytext=(3.1, 10.4),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=2), zorder=4)

ax.text(3.55, 7.8, "Élimination linéaire (CL)", ha='center',
        fontsize=12, color=DARK, zorder=3)
ax.text(3.55, 7.2, "↓ clairance plasmatique", ha='center',
        fontsize=11, color=GREY, style='italic', zorder=3)

arr(ax, -0.2, 10.4, 0.7, 10.4, color=DARK, lw=3)
ax.text(-0.1, 11.1, "Dose\nIV", ha='center', fontsize=12,
        fontweight='bold', color=DARK)

# ── BLOC 2 : DAMAGE ─────────────────────────────────────
ax.add_patch(FancyBboxPatch((7.5, 7.0), 3.4, 6.0, boxstyle="round,pad=0.15",
    facecolor=ORG_L, edgecolor=ORG, linewidth=3, zorder=1))
ax.text(9.2, 12.7, "Bloc 2 — DOMMAGE", ha='center',
        fontsize=15, fontweight='bold', color=ORG, zorder=3)
rbox(ax, 7.9, 8.2, 2.6, 2.8, 'DAMAGE', ORG, '#8B3210', fs=15)

arr(ax, 6.8, 10.4, 7.5, 9.6, color=DARK, lw=3)
ax.text(7.15, 11.1, "C_ADC\n(µM)", ha='center', fontsize=11,
        fontweight='bold', color=DARK)

# ── BLOC 3 : PD / HÉMATOPOÏÈSE ──────────────────────────
ax.add_patch(FancyBboxPatch((11.5, 1.8), 12.2, 12.0, boxstyle="round,pad=0.15",
    facecolor=GRN_L, edgecolor=GRN_M, linewidth=3, zorder=1))
ax.text(17.6, 13.5, "Bloc 3 — MODÈLE PD / HÉMATOPOÏÈSE", ha='center',
        fontsize=16, fontweight='bold', color=GRN, zorder=3)

# MPP
rbox(ax, 15.2, 11.5, 4.2, 1.3, 'MPP  (cellules souches)', GRN_M, GRN, fs=12)

# CMP + MEP
rbox(ax, 12.0,  9.0, 3.2, 1.2, 'CMP', GRN_M, GRN, fs=13)
rbox(ax, 18.0,  9.0, 3.2, 1.2, 'MEP', GRN_M, GRN, fs=13)

arr(ax, 16.2, 11.5, 13.6,  10.2, color='white', lw=2.5)
arr(ax, 18.2, 11.5, 19.6,  10.2, color='white', lw=2.5)

# Transit
rbox(ax, 12.0, 6.5, 3.2, 1.2, 'Neutrophiles\n+ Monocytes (transit)', GRN_M, GRN, fs=10, lw=2)
arr(ax, 13.6,  9.0, 13.6,  7.7, color='white', lw=2.5)

rbox(ax, 18.0, 6.5, 3.2, 1.2, 'Réticulocytes (transit)', GRN_M, GRN, fs=11, lw=2)
arr(ax, 19.6,  9.0, 19.6,  7.7, color='white', lw=2.5)

# Mature
rbox(ax, 11.8, 3.8, 3.2, 1.2, 'Neutrophiles\n+ Monocytes', GRN_D, GRN, fs=11, lw=2)
rbox(ax, 15.8, 3.8, 2.8, 1.2, 'GR (RBC)',                  GRN_D, GRN, fs=12, lw=2)
rbox(ax, 19.4, 3.8, 3.2, 1.2, 'Plaquettes (Plt)',           GRN_D, GRN, fs=11, lw=2)

arr(ax, 13.6, 6.5, 13.4, 5.0, color='white', lw=2.5)
arr(ax, 19.6, 6.5, 18.0, 5.0, color='white', lw=2.5)
arr(ax, 20.2, 6.5, 21.0, 5.0, color='white', lw=2.5)

# Feedback dashed
darr(ax, 12.2, 5.0, 15.4, 12.0, color='white', lw=2, rad= 0.35)
darr(ax, 22.0, 5.0, 19.2, 12.0, color='white', lw=2, rad=-0.35)

# Drug effect arrows
arr(ax, 10.9, 11.8, 15.2, 11.8, color=ORG, lw=2.5, rad=-0.15)
ax.text(13.0, 11.4, "Effet\nmédicament", ha='center',
        fontsize=10, color=ORG, fontweight='bold', zorder=6)

arr(ax, 10.9,  9.8, 12.0,  9.5, color=ORG, lw=2.5)
arr(ax, 10.9,  9.2, 18.0,  9.5, color=ORG, lw=2.5, rad=0.18)

arr(ax, 10.9,  9.5, 11.5,  9.5, color=DARK, lw=3)

# ── GLOSSAIRE ────────────────────────────────────────────
ax.add_patch(FancyBboxPatch((0.3, 0.2), 11.0, 3.8, boxstyle="round,pad=0.12",
    facecolor='#F8F8F8', edgecolor='#AAAAAA', linewidth=2, zorder=1))
ax.text(5.8, 3.75, "GLOSSAIRE — Pharmacocinétique & Hématopoïèse",
        ha='center', fontsize=12, fontweight='bold', color=DARK, zorder=3)

pk_items = [
    ("Dose IV",    "Dose intraveineuse"),
    ("V₁ central", "Volume central (compartiment sanguin)"),
    ("V₂ périph.", "Volume périphérique (compartiment tissulaire)"),
    ("CL",         "Clairance (vitesse d'élimination)"),
    ("C_ADC (µM)", "Concentration plasmatique en micromolaire"),
]
hema_items = [
    ("MPP",    "Progéniteur multipotent (cellule souche auto-renouvelable)"),
    ("CMP",    "Progéniteur myéloïde commun (précurseur granulocytes/monocytes)"),
    ("MEP",    "Progéniteur mégakaryocyte-érythroïde (précurseur RBC/plaquettes)"),
    ("GR/RBC", "Globules rouges (érythrocytes)"),
    ("Plt",    "Plaquettes (issues des mégakaryocytes)"),
]

ax.text(0.6,  3.35, "Pharmacocinétique (PK)", fontsize=11, fontweight='bold', color=BLUE, zorder=3)
ax.text(5.9,  3.35, "Hématopoïèse",           fontsize=11, fontweight='bold', color=GRN,  zorder=3)

for i, (t, d) in enumerate(pk_items):
    y = 2.95 - i * 0.52
    ax.text(0.6,  y, f"• {t}",  fontsize=10, color=DARK,    fontweight='bold', zorder=3)
    ax.text(2.3,  y, d,         fontsize=10, color='#444444',                  zorder=3)

for i, (t, d) in enumerate(hema_items):
    y = 2.95 - i * 0.52
    ax.text(5.9,  y, f"• {t}",  fontsize=10, color=DARK,    fontweight='bold', zorder=3)
    ax.text(7.15, y, d,         fontsize=10, color='#444444',                  zorder=3)

ax.plot([5.6, 5.6], [0.4, 3.55], color='#CCCCCC', lw=1.5, zorder=3)

plt.tight_layout(pad=0.4)
plt.savefig('/home/user/Hema/outputs/figure_pkpd_model.png', dpi=150,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pkpd_model.png")

# ═══════════════════════════════════════════════════════
# FIGURE 2 — Pipeline 4 étapes
# ═══════════════════════════════════════════════════════
PURPLE = '#5B2C8D'; TEAL = '#1A6B6B'

fig2, ax2 = plt.subplots(figsize=(24, 12))
ax2.set_xlim(0, 24); ax2.set_ylim(0, 12)
ax2.axis('off'); fig2.patch.set_facecolor('white')

steps = [
    (0.3,
     "ÉTAPE 1 — Modèle\nhématopoïétique\n(Fornari et al. 2019)",
     BLUE_M,
     ["Données & méthodes",
      "• Rat / carboplatine",
      "• Données digitalisées (Fornari 2019)",
      "• EDO 8-compartiments :",
      "  MPP→CMP/MEP→Neut/Mono/Ret/RBC/Plt",
      "• Compartiments de transit + feedbacks",
      "• Outils R : deSolve · ggplot2"],
     ["Résultat clé",
      "• RMSE < 15 % (cellules circulantes)",
      "• Implémentation validée"],
     BLUE),

    (6.3,
     "ÉTAPE 2 — Hématotoxicité\nT-DXd\n(Validation FDA BLA 761139)",
     '#B5510D',
     ["Données & méthodes",
      "• ADC 2-cpt PK + TMDD (modèle NHP)",
      "• Transposition Rat → NHP → Humain",
      "• Slopes rat transposés via IC50-scaling",
      "• Sim. N = 300 patients virtuels",
      "• Outils R : optim() · CTCAE v5 · ggplot2"],
     ["Résultat clé",
      "• Neutropénie tout grade : 29 % vs 29 % FDA",
      "• Thrombocytopénie G3-4 : 3 % vs 3 % FDA"],
     '#B5510D'),

    (12.3,
     "ÉTAPE 3 — Hématotoxicité\ncomposé (NHP interne)\n(données confidentielles)",
     TEAL,
     ["Données & méthodes",
      "• 8 singes cynomolgus, 4 doses",
      "  (4 / 13 / 26 / 39 mg/kg)",
      "• Même cadre PK/PD que T-DXd",
      "• Ajustement PK individuel (Nelder-Mead)",
      "• Calibration visuelle des Slopes PD",
      "• Outils R : rxode2 · optim() · PKNCA"],
     ["Résultat clé",
      "• Profils PK/PD individuels concordants",
      "• CL ≈ 1,8 mL/h/kg ;  t½ ≈ 55 h",
      "• Dose-proportionnalité validée"],
     TEAL),

    (18.3,
     "ÉTAPE 4 — Prédiction\nclinique\n(Transposition NHP → Humain)",
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
      "• Support dossier réglementaire IND/CTA"],
     PURPLE),
]

W = 5.6
for x0, header, hbg, bullets, results, badge_c in steps:
    # Outer box
    ax2.add_patch(FancyBboxPatch((x0, 0.7), W, 11.0, boxstyle="round,pad=0.12",
        facecolor='#F7F7F7', edgecolor=badge_c, linewidth=2.5, zorder=1))

    # Header band
    ax2.add_patch(FancyBboxPatch((x0, 8.8), W, 2.9, boxstyle="round,pad=0.06",
        facecolor=hbg, edgecolor=hbg, linewidth=0, zorder=2))
    ax2.text(x0 + W/2, 10.25, header, ha='center', va='center',
             fontsize=11, fontweight='bold', color='white', zorder=3,
             multialignment='center', linespacing=1.5)

    # Bullets
    for i, line in enumerate(bullets):
        y = 8.45 - i * 0.67
        fw = 'bold' if i == 0 else 'normal'
        fs = 11 if i == 0 else 10
        ax2.text(x0 + 0.2, y, line, fontsize=fs, fontweight=fw,
                 color='#333333', zorder=3)

    # Results block
    ax2.add_patch(FancyBboxPatch((x0 + 0.15, 0.85), W - 0.3, 2.9, boxstyle="round,pad=0.08",
        facecolor=badge_c, edgecolor=badge_c, linewidth=0, alpha=0.13, zorder=2))
    for i, line in enumerate(results):
        y = 3.5 - i * 0.62
        fw = 'bold' if i == 0 else 'normal'
        fc = badge_c if i == 0 else '#222222'
        ax2.text(x0 + 0.3, y, line, fontsize=10, fontweight=fw, color=fc, zorder=3)

# Arrows between steps
for x_arr in [5.9, 11.9, 17.9]:
    ax2.annotate('', xy=(x_arr + 0.35, 10.25), xytext=(x_arr - 0.1, 10.25),
                 arrowprops=dict(arrowstyle='->', color='#888888', lw=4), zorder=6)

# Step labels
step_colors = [BLUE_M, '#B5510D', TEAL, PURPLE]
labels = ["Étape 1", "Étape 2", "Étape 3", "Étape 4"]
label_xs = [3.1, 9.1, 15.1, 21.1]
for lx, lbl, lc in zip(label_xs, labels, step_colors):
    rbox(ax2, lx - 1.0, 0.05, 2.0, 0.55, lbl, lc, lc, fs=11, lw=2, zorder=3)

plt.tight_layout(pad=0.4)
plt.savefig('/home/user/Hema/outputs/figure_pipeline_steps.png', dpi=150,
            bbox_inches='tight', facecolor='white')
plt.close()
print("✓ figure_pipeline_steps.png")
