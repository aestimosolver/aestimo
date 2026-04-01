
import os
import numpy as np
import matplotlib.pyplot as plt

def load_dat(path):
    return np.loadtxt(path)

base = 'gui_headless_output'
light_path = os.path.join(base, 'T_300', 'av_curr.dat')
dark_path = os.path.join(base, 'Dark_300', 'av_curr.dat')

if not os.path.exists(light_path) or not os.path.exists(dark_path):
    print("Files not found.")
    exit(1)

light_data = load_dat(light_path)
dark_data = load_dat(dark_path)

# Aestimo saves: Volts, Current(A/m2)
v_l = light_data[:, 0]
j_l = light_data[:, 1] * 0.1 # A/m2 -> mA/cm2

v_d = dark_data[:, 0]
j_d = dark_data[:, 1] * 0.1 # A/m2 -> mA/cm2

plt.figure(figsize=(10, 8))

# J-V Curve
plt.subplot(2, 1, 1)
plt.plot(v_d, j_d, 'k--', label='Dark')
plt.plot(v_l, j_l, 'r-', label='Light (AM1.5G)')
plt.axhline(0, color='gray', lw=0.5)
plt.axvline(0, color='gray', lw=0.5)
plt.xlabel('Voltage (V)')
plt.ylabel('Current Density (mA/cm²)')
plt.title('Self-Consistent SI-Unified J-V Characteristic (300K)')
plt.legend()
plt.grid(True, alpha=0.3)

# P-V Curve
plt.subplot(2, 1, 2)
# Power quadrant should be V > 0 and J < 0 (using -J for positive power)
p_l = v_l * (-j_l)
plt.plot(v_l, p_l, 'b-')
plt.xlabel('Voltage (V)')
plt.ylabel('Power Density (mW/cm²)')
plt.title('P-V Characteristic (mW/cm²)')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('final_gui_si_result.png')
print("Successfully saved final_gui_si_result.png")
