import sys
import os
from os import path
import numpy as np
import matplotlib
matplotlib.use('Agg') # Non-interactive backend
import matplotlib.pyplot as plt

# Add root to sys.path
sys.path.append(os.getcwd())

import aestimo
from aeslibs.experimental_validation import (
    load_experimental_data, load_current_from_avcurr,
    compute_error_metrics, generate_validation_report,
    apply_parasitic_resistances, calculate_ideality_factor
)

# 1. Define inputs for Silicon p-n junction
material = [
    [2500.0, 'Si', 0.0, 0.0, 1e18, 'p', 'b'],
    [2500.0, 'Si', 0.0, 0.0, 1e18, 'n', 'b']
]
gridfactor = 10
inputfilename = "verify_fix"

class InputObj:
    def __init__(self):
        self.material = material
        self.vmin = 0.0
        self.vmax = 0.8
        self.Each_Step = 0.05
        self.T = 300.0
        self.computation_scheme = 9
        self.gridfactor = gridfactor
        self.maxgridpoints = 200000
        self.mat_type = 'Zincblende'
        self.Quantum_Regions = False
        self.Quantum_Regions_boundary = np.zeros((2, 2))
        self.surface = np.zeros(2)
        self.inputfilename = inputfilename
        self.device_area = 1e-4

input_obj = InputObj()
x_max = sum([layer[0] for layer in material])
n_max = int(x_max / gridfactor)
input_obj.dop_profile = np.zeros(n_max)

# Run simulation
print("Starting simulation to verify fix...")
aestimo.output_directory = os.path.join(os.getcwd(), "verify_fix_output")
if not os.path.isdir(aestimo.output_directory):
    os.makedirs(aestimo.output_directory, exist_ok=True)
aestimo.run_aestimo(input_obj, drawFigures=False, show=False)

# Extract and validate
output_dir = aestimo.output_directory
exp_file = "examples/experimental_data/si_pn_experimental_iv.csv"

if not os.path.exists(output_dir):
    print(f"Error: Output directory {output_dir} not created.")
    sys.exit(1)

exp_voltage, exp_current = load_experimental_data(exp_file)
calc_v_int, calc_i_int = load_current_from_avcurr(output_dir, device_area_cm2=input_obj.device_area)

# Apply default parasitics (Rs=5, Rsh=10k)
calc_v, calc_i = apply_parasitic_resistances(calc_v_int, calc_i_int, Rs=5.0, Rsh=10000.0)

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.plot(exp_voltage, exp_current, 'ko', label='Experimental', alpha=0.6, markersize=4)
ax1.plot(calc_v, calc_i, 'r-', label='Simulation (Fixed Eg=1.12eV)', linewidth=2)
ax1.set_xlabel("Voltage (V)")
ax1.set_ylabel("Current (A)")
ax1.set_title("I-V Comparison (Linear)")
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.semilogy(exp_voltage, np.abs(exp_current), 'ko', label='Experimental', alpha=0.6, markersize=4)
ax2.semilogy(calc_v, np.abs(calc_i), 'r-', label='Simulation (Fixed Eg=1.12eV)', linewidth=2)
ax2.set_xlabel("Voltage (V)")
ax2.set_ylabel("Current (log A)")
ax2.set_title("I-V Comparison (Semi-log)")
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("voltage_fix_verification.png")
print("Verification plot saved to voltage_fix_verification.png")

# Compute metrics
metrics = compute_error_metrics(exp_current[exp_voltage >= 0.1], 
                                np.interp(exp_voltage[exp_voltage >= 0.1], calc_v, calc_i))
print("\nNew Metrics:")
for k, v in metrics.items():
    print(f"{k}: {v}")
