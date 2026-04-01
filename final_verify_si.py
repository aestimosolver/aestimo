
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from aestimo import run_aestimo

# Load project
with open('test_solar_cell.json', 'r') as f:
    config = json.load(f)

# Force options for solar study
config['device_type'] = "Solar Cell / Photodetector"
config['solver'] = "9: Coupled Newton (Robust)"
config['vmin'] = 0.0
config['vmax'] = 1.2
config['vstep'] = 0.05
config['area'] = 1.0 # cm2

# 1. DARK Case
config['G_optical'] = 0.0
_, _, res_dark, _ = run_aestimo(config, drawFigures=False, show=False)

# 2. LIGHT Case (AM1.5G 100 mW/cm2)
config['G_optical'] = 5e18 # cm-3 s-1
_, _, res_light, _ = run_aestimo(config, drawFigures=False, show=False)

# Plotting
v_dark = res_dark.Va_t
j_dark = res_dark.av_curr * 0.1 # A/m2 -> mA/cm2

v_light = res_light.Va_t
j_light = res_light.av_curr * 0.1 # A/m2 -> mA/cm2

plt.figure(figsize=(10, 8))

# J-V Curve
plt.subplot(2, 1, 1)
plt.plot(v_dark, j_dark, 'k--', label='Dark')
plt.plot(v_light, j_light, 'r-', label='Light (AM1.5G)')
plt.axhline(0, color='gray', lw=0.5)
plt.axvline(0, color='gray', lw=0.5)
plt.xlabel('Voltage (V)')
plt.ylabel('Current Density (mA/cm²)')
plt.title('Self-Consistent J-V Characteristic (300K)')
plt.legend()
plt.grid(True, alpha=0.3)

# P-V Curve
plt.subplot(2, 1, 2)
p_light = v_light * (-j_light) # Power density (mW/cm2)
plt.plot(v_light, p_light, 'b-')
plt.xlabel('Voltage (V)')
plt.ylabel('Power Density (mW/cm²)')
plt.title('P-V Characteristic')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('final_gui_verification.png')
print("Successfully saved final_gui_verification.png")
