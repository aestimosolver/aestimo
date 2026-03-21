#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Systematic optimization of Rs and Rsh for InGaN p-n junction validation.
This script sweeps through different Rs/Rsh combinations to find the best fit.
"""

import numpy as np
import os
import sys
from os import path

# Setup paths
root_dir = os.getcwd()
if os.path.exists(os.path.join(root_dir, 'aeslibs')):
    sys.path.append(root_dir)

from aeslibs.experimental_validation import (
    load_experimental_data,
    load_current_from_avcurr,
    apply_parasitic_resistances,
    compute_error_metrics
)

# Load experimental data
exp_file = "examples/experimental_data/ingan_pn_experimental_iv.csv"
exp_v, exp_i = load_experimental_data(exp_file)

# Load simulation data (from last run with Rsh=1000)
output_dir = "pn_with_experimental_validation_ingan_output"
device_area = 5e-4  # cm^2

sim_v_int, sim_i_int = load_current_from_avcurr(output_dir, device_area_cm2=device_area)

# Define search grid
Rs_values = [30, 35, 40, 45, 50, 55]  # Ohm
Rsh_values = [800, 900, 1000, 1100, 1200, 1500]  # Ohm

print("="*70)
print("InGaN Rs/Rsh Optimization")
print("="*70)
print(f"Testing {len(Rs_values)} Rs values × {len(Rsh_values)} Rsh values = {len(Rs_values)*len(Rsh_values)} combinations\n")

best_r2 = -999
best_params = None
results = []

for Rs in Rs_values:
    for Rsh in Rsh_values:
        # Apply parasitics
        sim_v, sim_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=Rs, Rsh=Rsh)
        
        # Interpolate to experimental points
        sim_i_interp = np.interp(exp_v, sim_v, sim_i)
        
        # Compute metrics
        metrics = compute_error_metrics(exp_i, sim_i_interp)
        r2 = metrics['r2']
        log_rmse = metrics['log_rmse']
        
        results.append({
            'Rs': Rs,
            'Rsh': Rsh,
            'R2': r2,
            'Log_RMSE': log_rmse,
            'RMSE': metrics['rmse'],
            'MAE': metrics['mae']
        })
        
        # Track best
        if r2 > best_r2:
            best_r2 = r2
            best_params = (Rs, Rsh)
        
        print(f"Rs={Rs:4.0f} Ohm, Rsh={Rsh:5.0f} Ohm -> R2={r2:7.4f}, Log-RMSE={log_rmse:6.3f}")

print("\n" + "="*70)
print("OPTIMIZATION RESULTS")
print("="*70)
print(f"Best R2 = {best_r2:.4f}")
print(f"Optimal Rs  = {best_params[0]:.0f} Ohm")
print(f"Optimal Rsh = {best_params[1]:.0f} Ohm")
print("="*70)

# Save results to file
with open("ingan_optimization_results.txt", "w") as f:
    f.write("Rs(Ohm)\tRsh(Ohm)\tR2\tLog-RMSE\tRMSE\tMAE\n")
    for res in results:
        f.write(f"{res['Rs']}\t{res['Rsh']}\t{res['R2']:.4f}\t{res['Log_RMSE']:.4f}\t{res['RMSE']:.4e}\t{res['MAE']:.4e}\n")

print(f"\nDetailed results saved to: ingan_optimization_results.txt")
