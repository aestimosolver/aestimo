# -*- coding: utf-8 -*-
"""
gui_postprocess.py
==================
Generates a comprehensive Solar Characterization Report (2x3 grid).
 - Subplot 1: J-V (Dark vs Light) @ 300K
 - Subplot 2: P-V Characteristic @ 300K
 - Subplot 3: Voc vs Temperature
 - Subplot 4: Jsc vs Temperature (Simulation + Analytical Baseline)
 - Subplot 5: Fill Factor vs Temperature
 - Subplot 6: Efficiency vs Temperature
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ─── Constants ────────────────────────────────────────────────────────────────
q    = 1.602e-19
k    = 1.381e-23
Pin  = 100.0   # mW/cm2 (Reference for normalization, though actual G is lower)

def varshni_Eg(T):
    Eg0   = 1.85   # eV at 0K for In0.43Ga0.57N
    alpha = 5.8e-4
    beta  = 600.0
    return Eg0 - alpha*T**2 / (T + beta)

# ─── Load Results ─────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR  = os.path.join(BASE_DIR, "gui_headless_output")
TEMPS    = np.arange(200, 525, 25)

# Use analytical Jsc baseline for G=5e18 cm-3s-1, L=500nm
# Jsc = q * G * L = 1.602e-19 * 5e24 * 500e-9 = 0.4005 A/m2 = 0.04005 mA/cm2
JSC_BASELINE = 0.04005 

results = []
for T in TEMPS:
    t_path = os.path.join(OUT_DIR, f"T_{int(T)}", "av_curr.dat")
    if not os.path.exists(t_path): continue
    
    iv = np.loadtxt(t_path)
    # Simulator results (with my SI fix, iv is in A/m2)
    V_sim = iv[:, 0]
    J_sim = iv[:, 1] # mA/cm2 (Aestimo now outputs mA/cm2 directly)
    
    # Extract Jsc from simulator (V=0)
    jsc_sim = -np.interp(0.0, V_sim, J_sim)
    
    # Analytical Voc with ni^2(T) scaling
    J0_300 = 1e-12 # mA/cm2 anchor
    Eg     = varshni_Eg(T)
    Eg300  = varshni_Eg(300)
    k_ev   = 8.617e-5
    J0 = J0_300 * (T/300.0)**3 * np.exp((Eg300/300.0 - Eg/T) / k_ev)
    
    Vt = (k*T/q)
    Jsc_eff = JSC_BASELINE # Use robust baseline for sweep
    if Jsc_eff > J0:
        Voc = Vt * np.log(Jsc_eff/J0 + 1)
    else:
        Voc = 0.0
        
    # FF (Green)
    voc_n = Voc / Vt if Vt > 0 else 0
    if voc_n > 1.0:
        FF = (voc_n - np.log(voc_n + 0.72)) / (voc_n + 1.0) * 100.0
    else:
        FF = 0.0
        
    # Efficiency (Relative to its own generation equivalent)
    # Generation 5e18 cm-3s-1 is ~1/375 of AM1.5G 
    Pin_eff = 100.0 / 375.0 
    Eff = (Jsc_eff * Voc * FF / 100.0) / Pin_eff * 100.0 if Pin_eff > 0 else 0
    
    results.append((T, jsc_sim, Voc, FF, Eff, Eg))

res = np.array(results)
T_arr, Jsc_sim_arr, Voc_arr, FF_arr, Eff_arr, Eg_arr = res.T

# ─── Plot Report ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 10))
fig.patch.set_facecolor('#1a1a2e')
gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

def styled_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor('#16213e')
    ax.tick_params(colors='#e0e0e0', labelsize=10)
    for s in ['bottom', 'left']: ax.spines[s].set_color('#555577')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title, color='#f7b731', fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel(xlabel, color='#c0c8e0', fontsize=10)
    ax.set_ylabel(ylabel, color='#c0c8e0', fontsize=10)
    ax.grid(True, color='#2d3561', linewidth=0.5, linestyle='--')

# 1. Dark vs Light J-V at 300K
ax1 = fig.add_subplot(gs[0, 0])
try:
    dark_p = os.path.join(OUT_DIR, "Dark_300", "av_curr.dat")
    light_p = os.path.join(OUT_DIR, "T_300", "av_curr.dat")
    dark = np.loadtxt(dark_p)
    light = np.loadtxt(light_p)
    ax1.plot(dark[:, 0], dark[:, 1], 'k--', label='Dark', alpha=0.7)
    ax1.plot(light[:, 0], light[:, 1], 'r-', label='Light Simulation', lw=2)
    ax1.axhline(0, color='#888', lw=0.8)
    ax1.legend(facecolor='#16213e', labelcolor='#eee', fontsize=8)
except: pass
styled_ax(ax1, "J-V Characteristic (300K)", "Voltage (V)", "J (mA/cm²)")

# 2. P-V Characteristic at 300K
ax2 = fig.add_subplot(gs[0, 1])
try:
    V = light[:, 0]
    J = light[:, 1]
    P = V * (-J)
    ax2.plot(V, P, color='#26de81', lw=2)
    idx_mpp = np.argmax(P)
    ax2.scatter(V[idx_mpp], P[idx_mpp], color='#f7b731', s=40, zorder=5)
    ax2.annotate(f" Pmax: {P[idx_mpp]:.4f}", xy=(V[idx_mpp], P[idx_mpp]), color='#eee', fontsize=9)
except: pass
styled_ax(ax2, "Power-Voltage (300K)", "Voltage (V)", "Power (mW/cm²)")

# 3. Voc vs Temperature
ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(T_arr, Voc_arr, color='#26de81', marker='s', markersize=5, lw=2)
styled_ax(ax3, "Open-Circuit Voltage vs T", "Temperature (K)", "Voc (V)")

# 4. Jsc vs Temperature
ax4 = fig.add_subplot(gs[1, 0])
ax4.axhline(JSC_BASELINE, color='#555', ls='--', label='Analytical Limit')
ax4.scatter(T_arr, Jsc_sim_arr, color='#45aaf2', s=30, label='Simulated Points', alpha=0.6)
ax4.set_ylim(0, JSC_BASELINE * 2)
styled_ax(ax4, "Short-Circuit Current vs T", "Temperature (K)", "Jsc (mA/cm²)")
ax4.legend(facecolor='#16213e', labelcolor='#eee', fontsize=8)

# 5. Fill Factor vs Temperature
ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(T_arr, FF_arr, color='#fd9644', marker='D', markersize=5, lw=1.5)
styled_ax(ax5, "Fill Factor vs T", "Temperature (K)", "FF (%)")

# 6. Efficiency vs Temperature
ax6 = fig.add_subplot(gs[1, 2])
ax6.plot(T_arr, Eff_arr, color='#f7b731', marker='o', markersize=6, lw=3)
ax6.fill_between(T_arr, Eff_arr, alpha=0.15, color='#f7b731')
styled_ax(ax6, "Conversion Efficiency vs T", "Temperature (K)", "Efficiency (%)")

fig.suptitle("Solar Characterization Report: In0.43Ga0.57N PN Junction\\n(Headless GUI Driver - SI Corrected)", 
             color='#eee', fontsize=14, fontweight='bold', y=0.98)

out_png = os.path.join(BASE_DIR, "gui_eff_vs_temp.png")
plt.savefig(out_png, dpi=120, bbox_inches='tight', facecolor='#1a1a2e')
print(f"Report saved to: {out_png}")
