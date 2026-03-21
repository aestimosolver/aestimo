#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
# Input File Description:  Si p/n junction with experimental validation.
# ----------------------------------------------------------------------
# This example demonstrates how to validate drift-diffusion simulation
# results against experimental I-V data.
# 
# The simulation uses:
# - Drift-Diffusion solver (scheme 9: Gummel & Newton map)
# - Si p-n junction with matched doping profile (loaded from exp file)
# - Voltage sweep to generate I-V characteristics
# - Automatic comparison with experimental data
# ----------------------------------------------------------------------

import numpy as np
import os
import sys
from os import path

# Detect if we are in root or examples dir to fix paths
root_dir = os.getcwd()
if os.path.exists(os.path.join(root_dir, 'aeslibs')):
    sys.path.append(root_dir)
else:
    # try moving up if in examples
    sys.path.append(path.join(path.dirname(__file__), '..'))

try:
    from aeslibs.experimental_validation import load_experimental_metadata
except ImportError:
    # Fallback if aeslibs not found in path
    def load_experimental_metadata(f): return {}

# ----------------
# GENERAL SETTINGS
# ----------------

# TEMPERATURE
T = 300.0  # Kelvin

# COMPUTATIONAL SCHEME
# 9: Schrodinger-Poisson-Drift_Diffusion using Gummel & Newton map
computation_scheme = 9

# QUANTUM
# Total subband number to be calculated
subnumber_h = 2
subnumber_e = 2

# VOLTAGE SWEEP PARAMETERS
# These should match the experimental data voltage range
Fapplied = 0.0  # Applied electric field (V/m)
vmax = 0.80     # Maximum voltage (V)
vmin = 0.0      # Minimum voltage (V)
Each_Step = 0.01  # Voltage step (V)

# --------------------------------
# EXPERIMENTAL CONFIGURATION
# --------------------------------
# Set to True to enable experimental data comparison
enable_experimental_validation = True

# Path to experimental I-V data file
# Default path assuming running from project root
experimental_iv_file = "examples/experimental_data/si_pn_experimental_iv.csv"

if not os.path.exists(experimental_iv_file):
     # Robust path finding
     script_dir = os.path.dirname(os.path.abspath(__file__))
     
     # Check 1: script_dir/experimental_data/... (if script is in examples folder)
     candidate1 = os.path.join(script_dir, "experimental_data", "si_pn_experimental_iv.csv")
     
     # Check 2: script_dir/examples/experimental_data/... (if script is in root)
     candidate2 = os.path.join(script_dir, "examples", "experimental_data", "si_pn_experimental_iv.csv")
     
     if os.path.exists(candidate1):
         experimental_iv_file = candidate1
     elif os.path.exists(candidate2):
         experimental_iv_file = candidate2

# DYNAMIC PARAMETER LOADING
# Validate that we match the experiment dimensions
print(f"Reading metadata from: {experimental_iv_file}")
metadata = load_experimental_metadata(experimental_iv_file)

# Default values if metadata missing
default_area = 1e-4
default_doping = 1e18

# Load values
device_area = metadata.get('area_cm2', default_area)
doping_p = metadata.get('doping_p', default_doping)
doping_n = metadata.get('doping_n', default_doping)

print(f"  - Device Area: {device_area:.2e} cm^2")
print(f"  - Doping (p): {doping_p:.2e} cm^-3")
print(f"  - Doping (n): {doping_n:.2e} cm^-3")

# Parasitic Resistances (Optimized)
Rs_ext = 2330.0  # Series Resistance (Ohm) - High value to match low Exp current
Rsh_ext = 5e7   # Shunt Resistance (Ohm) - To match subthreshold leakage

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
gridfactor = 10  # nm
maxgridpoints = 200000
mat_type = 'Zincblende'

# DEVICE STRUCTURE
# Si p-n junction with symmetric doping
# Doping levels matched to experimental device dynamically
material = [
    [2500.0, 'Si', 0.0, 0.0, doping_p, 'p', 'b'],  # p-side
    [2500.0, 'Si', 0.0, 0.0, doping_n, 'n', 'b']   # n-side
]

# ---------------------------------------- 
# STANDARD SETUP (DO NOT MODIFY)
# ----------------------------------------
x_max = sum([layer[0] for layer in material])

def round2int(x):
    return int(x + 0.5)

n_max = round2int(x_max / gridfactor)

dop_profile = np.zeros(n_max)
Quantum_Regions = False
Quantum_Regions_boundary = np.zeros((2, 2))
surface = np.zeros(2)
inputfilename = "sample_pn_with_experimental_validation"

# ----------------------------------------
# MAIN EXECUTION
# ----------------------------------------
if __name__ == "__main__":
    input_obj = vars()
    
    # Run the simulation
    import aestimo
    print("="*60)
    print("Running Si p-n Junction Simulation")
    print("="*60)
    print(f"Device structure: {len(material)} layers")
    print(f"Total device length: {x_max} nm")
    print(f"Voltage range: {vmin}V to {vmax}V (step: {Each_Step}V)")
    print("="*60)
    
    # Run simulation
    # Pass input_obj to run_aestimo
    results = aestimo.run_aestimo(input_obj)
    
    print("\nSimulation completed!")
    
    # ----------------------------------------
    # EXPERIMENTAL VALIDATION
    # ----------------------------------------
    if enable_experimental_validation:
        print("\n" + "="*60)
        print("EXPERIMENTAL VALIDATION")
        print("="*60)
        
        try:
            from aeslibs.experimental_validation import (
                load_experimental_data,
                load_current_from_avcurr,
                apply_parasitic_resistances,
                compute_error_metrics,
                plot_iv_comparison,
                generate_validation_report
            )
            
            # Load experimental data
            print(f"\nLoading experimental data...")
            exp_voltage, exp_current = load_experimental_data(experimental_iv_file)
            print(f"Loaded {len(exp_voltage)} experimental data points")
            print(f"Voltage range: {exp_voltage.min():.2f}V to {exp_voltage.max():.2f}V")
            
            # Extract simulation results
            print("\nExtracting simulation results...")
            output_dir = f"{inputfilename}_output"
            
            if os.path.exists(output_dir):
                # Extract actual currents from simulation output
                # Using load_current_from_avcurr to get reliable total current
                calc_v_int, calc_i_int = load_current_from_avcurr(output_dir, device_area_cm2=device_area)
                
                print(f"Extracted {len(calc_v_int)} simulation points.")
                
                # Apply Parasitic Resistances (External Rs mode)
                print(f"\nApplying Parasitic Resistances:")
                print(f"- Rs: {Rs_ext} Ohm")
                print(f"- Rsh: {Rsh_ext:.1e} Ohm")
                
                calc_v, calc_i = apply_parasitic_resistances(calc_v_int, calc_i_int, Rs=Rs_ext, Rsh=Rsh_ext)
                
                # Interpolate simulation to experimental voltage points for metrics
                sim_current_interp = np.interp(exp_voltage, calc_v, calc_i)
                
                # Compute error metrics (Focus on active region V > 0.1V)
                print("\nComputing error metrics (Active region V > 0.1V)...")
                mask = exp_voltage >= 0.1
                metrics = compute_error_metrics(exp_current[mask], sim_current_interp[mask])
                
                # Generate validation report
                report_path = path.join(output_dir, "validation_report.txt")
                generate_validation_report(metrics, output_path=report_path)
                
                # Plot comparison
                plot_path = path.join(output_dir, "iv_comparison.png")
                print(f"\nGenerating I-V comparison plot...")
                plot_iv_comparison(
                    exp_voltage, exp_current,
                    calc_v, calc_i,
                    output_path=plot_path,
                    show=False
                )
                
                print("\n" + "="*60)
                print("Validation complete!")
                print(f"- Metrics: Log-RMSE = {metrics.get('log_rmse', 'N/A'):.4f}")
                print(f"- Validation report: {report_path}")
                print(f"- Comparison plot: {plot_path}")
                print("="*60)
                
            else:
                print(f"Warning: Output directory not found: {output_dir}")
                print("Skipping validation.")
        
        except Exception as e:
            print(f"\nError during experimental validation: {e}")
            import traceback
            traceback.print_exc()
    
    print("\nDone!")
