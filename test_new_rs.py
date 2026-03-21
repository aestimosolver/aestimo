import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add aeslibs to path
sys.path.append(os.getcwd())
from aeslibs.experimental_validation import (
    load_experimental_data, 
    load_current_from_avcurr, 
    apply_parasitic_resistances, 
    compute_error_metrics
)

# Load data
exp_file = "examples/experimental_data/si_pn_experimental_iv.csv"
output_dir = "sample_pn_with_experimental_validation_output"
device_area = 1e-4

print(f"Loading experimental data from {exp_file}...")
exp_v, exp_i = load_experimental_data(exp_file)

print(f"Loading simulation data from {output_dir}...")
sim_v_int, sim_i_int = load_current_from_avcurr(output_dir, device_area_cm2=device_area)

print(f"Internal Current Range: {sim_i_int.min():.2e} to {sim_i_int.max():.2e} A")

# Params to test
rs_test = 2330.0
rsh_test = 5.0e7

print(f"\nTesting Parameters:")
print(f"Rs = {rs_test} Ohm")
print(f"Rsh = {rsh_test:.1e} Ohm")

# Apply
calc_v, calc_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=rs_test, Rsh=rsh_test)

# Interpolate for metrics
# Only use forward bias exp data
mask = exp_v >= 0.05
exp_v_fwd = exp_v[mask]
exp_i_fwd = exp_i[mask]

sim_i_interp = np.interp(exp_v_fwd, calc_v, calc_i)

metrics = compute_error_metrics(exp_i_fwd, sim_i_interp)

print("\nMetrics:")
for k, v in metrics.items():
    print(f"{k}: {v}")

# Plot
plt.figure(figsize=(10, 6))
plt.semilogy(exp_v, np.abs(exp_i), 'ko', label='Experimental', markersize=5, alpha=0.5)
plt.semilogy(calc_v, np.abs(calc_i), 'r-', label=f'Sim (Rs={rs_test}, Rsh={rsh_test:.1e})', linewidth=2)
# Original for comparison (Rs=14.5)
orig_v, orig_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=14.5, Rsh=1e12)
plt.semilogy(orig_v, np.abs(orig_i), 'b--', label='Original (Rs=14.5)', alpha=0.3)

plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.title(f"Diagnostic Fit (Log-RMSE: {metrics.get('log_rmse', 'N/A'):.4f})")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("diagnostic_fit.png")
print("\nPlot saved to diagnostic_fit.png")
