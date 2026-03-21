
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

# Add root to sys.path
sys.path.append(os.getcwd())

from aeslibs.experimental_validation import (
    load_experimental_data, load_current_from_avcurr,
    apply_parasitic_resistances, compute_error_metrics
)

# Configuration
output_dir = "verify_fix_output"
exp_file = "examples/experimental_data/si_pn_experimental_iv.csv"
area_cm2 = 1e-4

# Load Data
print(f"Loading experimental data from {exp_file}...")
exp_voltage, exp_current = load_experimental_data(exp_file)

print(f"Loading simulation data from {output_dir}...")
# Assuming simulation was run with Rs=0 internally
calc_v_int, calc_i_int = load_current_from_avcurr(output_dir, device_area_cm2=area_cm2)


# Optimization Loop
rs_values = np.arange(0, 15.0, 0.5)
errors = []

print("\nSweeping Rs values (0 to 15 Ohm)...")
best_rs = 0
min_error = float('inf')
metrics_at_best = {}

for rs in rs_values:
    # Apply external Rs
    calc_v, calc_i = apply_parasitic_resistances(calc_v_int, calc_i_int, Rs=rs, Rsh=1e12)
    
    # Interpolate to experimental voltage points
    # Filter for active region (V > 0.1)
    mask = exp_voltage >= 0.1
    sim_current_interp = np.interp(exp_voltage[mask], calc_v, calc_i)
    
    metrics = compute_error_metrics(exp_current[mask], sim_current_interp)
    
    # Use Log-RMSE as the primary optimization metric for I-V curves
    error = metrics['log_rmse']
    errors.append(error)
    
    if error < min_error:
        min_error = error
        best_rs = rs
        metrics_at_best = metrics

print("\n" + "="*40)
print(f"Optimization Complete")
print("="*40)
print(f"Best Rs: {best_rs:.2f} Ohm")
print(f"Minimum Log-RMSE: {min_error:.5f}")
print(f"R2 at Best Rs: {metrics_at_best['r2']:.5f}")
print("="*40)

# Plot Error Landscape
plt.figure(figsize=(8, 6))
plt.plot(rs_values, errors, 'b-o')
plt.axvline(x=best_rs, color='r', linestyle='--', label=rf'Best Rs={best_rs} $\Omega$')
plt.xlabel(r'Series Resistance ($\Omega$)')
plt.ylabel('Log-RMSE (decades)')
plt.title('Optimization of Series Resistance')
plt.legend()
plt.grid(True)
plt.savefig('rs_optimization_plot.png')
print("Saved optimization plot to rs_optimization_plot.png")

# Generate Comparison Plot at Best Rs
calc_v_opt, calc_i_opt = apply_parasitic_resistances(calc_v_int, calc_i_int, Rs=best_rs, Rsh=1e12)

plt.figure(figsize=(10, 6))
plt.semilogy(exp_voltage, np.abs(exp_current), 'ko', label='Experimental', markersize=4, alpha=0.6)
plt.semilogy(calc_v_opt, np.abs(calc_i_opt), 'r-', label=rf'Best Fit (Rs={best_rs}$\Omega$)', linewidth=2)
# Also show baseline for context
plt.semilogy(calc_v_int, np.abs(calc_i_int), 'b--', label=r'Baseline (Rs=0$\Omega$)', alpha=0.4)

plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.title(rf"I-V Comparison at Best Rs ({best_rs} $\Omega$)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('best_rs_fit.png')
print("Saved fit comparison to best_rs_fit.png")
