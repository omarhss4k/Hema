import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# -- Paramètres --
t = np.linspace(0, 35, 1000)

# Concentration plasmatique : perfusion 1h puis décroissance biexponentielle
# alpha rapide (distribution) + beta lente (élimination)
A, alpha = 1.0, 0.6
B, beta  = 0.35, 0.055
C = A * np.exp(-alpha * t) + B * np.exp(-beta * t)
C[t < 0] = 0

# Neutrophiles : décroissance retardée via compartiments de transit (MTT ~7j)
# Simulation simplifiée : sigmoïde de décroissance + rebond
Neut0 = 1.0
nadir_t = 14.0
nadir_v = 0.28
recover_t = 28.0

def neut_profile(t):
    # Phase de déclin (j0 → nadir j14)
    decline = Neut0 - (Neut0 - nadir_v) * (1 - np.exp(-0.18 * (t - 2)))
    decline = np.clip(decline, nadir_v, Neut0)
    # Phase de rebond (j14 → j28)
    recover = nadir_v + (Neut0 - nadir_v) * (1 - np.exp(-0.22 * (t - nadir_t)))
    recover = np.clip(recover, nadir_v, Neut0 * 1.05)
    result = np.where(t < nadir_t, decline, recover)
    result[t < 1] = Neut0
    return result

Neut = neut_profile(t)

# -- Figure --
fig, ax1 = plt.subplots(figsize=(11, 6))
fig.patch.set_facecolor('white')

# Axe 1 : Concentration
color_pk = '#2166ac'
l1, = ax1.plot(t, C, color=color_pk, linewidth=2.8, label='Concentration plasmatique (ADC)', zorder=3)
ax1.set_xlabel('Temps (jours)', fontsize=13)
ax1.set_ylabel('Concentration plasmatique\n(unité normalisée)', color=color_pk, fontsize=12)
ax1.tick_params(axis='y', labelcolor=color_pk)
ax1.set_xlim(0, 35)
ax1.set_ylim(0, 1.45)
ax1.fill_between(t, C, alpha=0.08, color=color_pk)

# Axe 2 : Neutrophiles
ax2 = ax1.twinx()
color_pd = '#d6604d'
l2, = ax2.plot(t, Neut, color=color_pd, linewidth=2.8,
               linestyle='--', label='Neutrophiles circulants', zorder=3)
ax2.set_ylabel('Neutrophiles\n(fraction de la baseline)', color=color_pd, fontsize=12)
ax2.tick_params(axis='y', labelcolor=color_pd)
ax2.set_ylim(0, 1.45)
ax2.axhline(1.0, color=color_pd, linewidth=0.8, linestyle=':', alpha=0.5)

# -- Annotations --
# Pic concentration
cmax_t = 0.3
cmax_v = A * np.exp(-alpha * cmax_t) + B * np.exp(-beta * cmax_t)
ax1.annotate('Cmax\n(~30 min)', xy=(cmax_t, cmax_v), xytext=(3.5, 1.28),
             fontsize=10.5, color=color_pk, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color=color_pk, lw=1.5))

# Nadir neutrophiles
ax2.annotate(f'Nadir neutrophiles\n(J{int(nadir_t)})',
             xy=(nadir_t, nadir_v), xytext=(17.5, 0.13),
             fontsize=10.5, color=color_pd, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color=color_pd, lw=1.5))

# Ligne verticale Cmax
ax1.axvline(cmax_t, color=color_pk, linewidth=1.0, linestyle=':', alpha=0.6)

# Ligne verticale nadir
ax1.axvline(nadir_t, color=color_pd, linewidth=1.0, linestyle=':', alpha=0.6)

# Flèche de décalage entre Cmax et nadir
ax1.annotate('', xy=(nadir_t, 1.35), xytext=(cmax_t, 1.35),
             arrowprops=dict(arrowstyle='<->', color='#555555',
                             lw=1.8, mutation_scale=15))
ax1.text((nadir_t + cmax_t) / 2, 1.38, 'Décalage ~14 jours\n(hystérèse PK/PD)',
         ha='center', va='bottom', fontsize=10, color='#333333',
         fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#fffde7',
                   edgecolor='#aaa', alpha=0.9))

# Zone grisée "médicament éliminé mais toxicité maximale non atteinte"
elim_t = 12
ax1.axvspan(elim_t, nadir_t, alpha=0.07, color='red',
            label='Médicament quasi-éliminé,\ntoxicité pas encore au max')
ax1.text(13.1, 0.62, 'Médicament quasi-éliminé\nmais nadir pas encore atteint',
         fontsize=9, color='#b22222', style='italic', va='center')

# -- Légende --
lines = [l1, l2]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper right', fontsize=10.5,
           framealpha=0.95, edgecolor='#cccccc')

# -- Titre --
ax1.set_title('Hystérèse PK/PD — Décalage temporel entre exposition et toxicité hématologique',
              fontsize=13, fontweight='bold', pad=14)

ax1.set_xticks(range(0, 36, 7))
ax1.set_xticklabels([f'J{d}' for d in range(0, 36, 7)], fontsize=11)

plt.tight_layout()
plt.savefig('/home/user/Hema/outputs/hysteresis_pkpd.png', dpi=200,
            bbox_inches='tight', facecolor='white')
print("Sauvegardé : outputs/hysteresis_pkpd.png")
