import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ---- Physics-based InGaN solar cell efficiency vs temperature model ----
# Material: In0.43Ga0.57N p-n junction, AM1.5G illumination
# Based on published literature values (see Dahal et al., Wu et al.)

q = 1.602e-19
k = 1.381e-23
Pin = 100.0  # mW/cm2 (1-sun AM1.5G)

temps = np.arange(200, 501, 25)

def varshni_Eg(T):
    # In0.43Ga0.57N Varshni parameters interpolated between GaN and InN
    Eg0 = 1.85   # eV at 0K
    alpha = 5.8e-4
    beta = 600.0
    return Eg0 - alpha * T**2 / (T + beta)

def ni_sq(T):
    # ni^2 ~ Nc*Nv*exp(-Eg/kT), scaling from reference at 300K
    Eg = varshni_Eg(T)
    Eg300 = varshni_Eg(300)
    return (T/300.0)**3 * np.exp(-q*(Eg - Eg300)/(k*T))

def Jsc_T(T):
    Eg300 = varshni_Eg(300)
    EgT = varshni_Eg(T)
    # Slightly more photons absorbed as Eg decreases, Jsc weakly increases
    base_jsc = 16.5  # mA/cm2 at 300K (literature: Dahal 2010, InGaN 43% In)
    return base_jsc * (Eg300 / EgT) ** 0.5

def J0_T(T):
    # J0 scales strongly with temperature via ni^2
    J0_300 = 1e-12  # A/cm2 at 300K (typical for InGaN)
    return J0_300 * ni_sq(T)

def Voc_T(T):
    Vt = k * T / q
    jsc = Jsc_T(T)
    j0 = J0_T(T)
    if j0 <= 0:
        return 0.0
    return Vt * np.log(jsc * 1e-3 / j0 + 1)  # J0 in A/cm2, Jsc in mA/cm2

def FF_T(voc_v, T):
    # Green's empirical formula
    Vt = k * T / q
    if voc_v <= 0:
        return 0.0
    voc_norm = voc_v / Vt
    return (voc_norm - np.log(voc_norm + 0.72)) / (voc_norm + 1)

results = []
for T in temps:
    jsc = Jsc_T(T)
    voc = Voc_T(T)
    ff  = FF_T(voc, T) * 100.0
    eff = (jsc * voc * ff / 100.0) / Pin * 100.0
    results.append((T, varshni_Eg(T), jsc, voc, ff, eff))

results = np.array(results)
T_arr, Eg_arr, Jsc_arr, Voc_arr, FF_arr, Eff_arr = results.T

# Print table
print()
print("=== InGaN Solar Cell (In0.43Ga0.57N) - Efficiency vs Temperature ===")
print(f"{'Temp (K)':>9} {'Eg (eV)':>9} {'Jsc (mA/cm2)':>14} {'Voc (V)':>9} {'FF (%)':>8} {'Eff (%)':>9}")
print("-" * 62)
for row in results:
    print(f"{row[0]:>9.0f} {row[1]:>9.4f} {row[2]:>14.3f} {row[3]:>9.4f} {row[4]:>8.2f} {row[5]:>9.3f}")

# Plot
fig = plt.figure(figsize=(14, 9))
fig.patch.set_facecolor('#1a1a2e')
gs = GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.4)

def styled_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor('#16213e')
    ax.tick_params(colors='#e0e0e0', labelsize=9)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#555577')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title, color='#aab4d4', fontsize=10, fontweight='bold', pad=8)
    ax.set_xlabel(xlabel, color='#c0c8e0', fontsize=9)
    ax.set_ylabel(ylabel, color='#c0c8e0', fontsize=9)
    ax.grid(True, color='#2d3561', linewidth=0.6, linestyle='--')

# Main efficiency plot
ax1 = fig.add_subplot(gs[0, :2])
ax1.plot(T_arr, Eff_arr, color='#f7b731', lw=2.5, marker='o', markersize=5, label='Efficiency (%)')
ax1.fill_between(T_arr, Eff_arr, alpha=0.12, color='#f7b731')
ax1.axvline(300, color='#aaaaaa', ls=':', lw=1.3, label='Std. Temp (300K)')
idx_300 = np.argmin(np.abs(T_arr - 300))
ax1.annotate(f" {Eff_arr[idx_300]:.2f}%", xy=(300, Eff_arr[idx_300]),
             color='#f7b731', fontsize=10, fontweight='bold')
styled_ax(ax1, 'Power Conversion Efficiency vs Temperature', 'Temperature (K)', 'Efficiency (%)')
ax1.legend(facecolor='#16213e', edgecolor='#555577', labelcolor='#e0e0e0', fontsize=9)

ax2 = fig.add_subplot(gs[0, 2])
ax2.plot(T_arr, Voc_arr, color='#26de81', lw=2.0, marker='s', markersize=4)
styled_ax(ax2, 'Open-Circuit Voltage', 'Temperature (K)', 'Voc (V)')

ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(T_arr, Jsc_arr, color='#45aaf2', lw=2.0, marker='^', markersize=4)
styled_ax(ax3, 'Short-Circuit Current', 'Temperature (K)', 'Jsc (mA/cm²)')

ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(T_arr, FF_arr, color='#fd9644', lw=2.0, marker='D', markersize=4)
styled_ax(ax4, 'Fill Factor', 'Temperature (K)', 'FF (%)')

ax5 = fig.add_subplot(gs[1, 2])
ax5.plot(T_arr, Eg_arr, color='#fc5c65', lw=2.0, marker='v', markersize=4)
styled_ax(ax5, 'Bandgap (Varshni)', 'Temperature (K)', 'Eg (eV)')

fig.suptitle(
    'InGaN Solar Cell (In\u2080.\u2084\u2083Ga\u2080.\u2085\u2087N) \u2014 Analytical Model: Efficiency vs Temperature',
    color='#e8ecff', fontsize=13, fontweight='bold', y=0.98
)

plt.savefig('eff_vs_temp.png', dpi=150, bbox_inches='tight', facecolor='#1a1a2e')
print("\nPlot saved as eff_vs_temp.png")
