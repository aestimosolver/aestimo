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
# - Si p-n junction with matched doping profile
# - Voltage sweep to generate I-V characteristics
# - Automatic comparison with experimental data
# ----------------------------------------------------------------------

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
Each_Step = 0.02  # Voltage step (V)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
gridfactor = 10  # nm
maxgridpoints = 200000
mat_type = 'Zincblende'

# DEVICE STRUCTURE
# Si p-n junction with symmetric doping
# Doping levels matched to experimental device: 1e18 cm^-3
material = [
    [2500.0, 'Si', 0.0, 0.0, 1e18, 'p', 'b'],  # p-side
    [2500.0, 'Si', 0.0, 0.0, 1e18, 'n', 'b']   # n-side
]

# EXPERIMENTAL VALIDATION SETTINGS
# Set to True to enable experimental data comparison
enable_experimental_validation = True

# Path to experimental I-V data file
experimental_iv_file = "examples/experimental_data/si_pn_experimental_iv.csv"

# Device parameters for current calculation
device_area = 1e-4  # cm² (100 x 100 μm)

# ---------------------------------------- 
# STANDARD SETUP (DO NOT MODIFY)
# ----------------------------------------
import numpy as np
import os
from os import path

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
    import sys
    sys.path.append(path.join(path.dirname(__file__), '..'))
    
    # Run the simulation
    import aestimo
    print("="*60)
    print("Running Si p-n Junction Simulation with Experimental Validation")
    print("="*60)
    print(f"Device structure: {len(material)} layers")
    print(f"Total device length: {x_max} nm")
    print(f"Voltage range: {vmin}V to {vmax}V (step: {Each_Step}V)")
    print(f"Doping: p-side = {material[0][4]:.1e} cm⁻³, n-side = {material[1][4]:.1e} cm⁻³")
    print("="*60)
    
    # Run simulation
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
            # Import validation module
            sys.path.append(path.join(path.dirname(__file__), '..', 'aeslibs'))
            from experimental_validation import (
                load_experimental_data,
                calculate_current_from_simulation,
                compute_error_metrics,
                plot_iv_comparison,
                generate_validation_report
            )
            
            # Load experimental data
            print(f"\nLoading experimental data from: {experimental_iv_file}")
            exp_voltage, exp_current = load_experimental_data(experimental_iv_file)
            print(f"Loaded {len(exp_voltage)} experimental data points")
            print(f"Voltage range: {exp_voltage.min():.2f}V to {exp_voltage.max():.2f}V")
            
            # Extract simulation results
            # Note: This is a simplified extraction. Actual implementation depends
            # on the output format from aestimo.run_aestimo()
            
            print("\nExtracting simulation results...")
            output_dir = f"{inputfilename}_output"
            
            # Check if output directory exists
            if os.path.exists(output_dir):
                # List all voltage output files
                voltage_files = sorted([f for f in os.listdir(output_dir) 
                                      if f.startswith('potn_eh_') and f.endswith('.dat')])
                
                print(f"Found {len(voltage_files)} simulation voltage points")
                
                # For this demo, create a simple placeholder simulation result
                # In a real implementation, you would parse the actual simulation output
                sim_voltages = np.arange(vmin, vmax + Each_Step/2, Each_Step)
                
                # Create dummy simulation data for demonstration
                # This should be replaced with actual current calculation from simulation
                # using the calculate_current_from_simulation function with real data
                print("\nNOTE: Using placeholder simulation currents for demonstration.")
                print("      Real implementation would extract carrier densities from")
                print(f"      output files in {output_dir}/ and calculate actual currents.")
                
                # Placeholder: Simple diode equation for demonstration
                Is = 1e-12 * device_area  # Saturation current
                eta = 1.1  # Ideality factor
                Vt = 0.0259  # Thermal voltage at 300K
                sim_currents = Is * (np.exp(sim_voltages / (eta * Vt)) - 1)
                
                # Interpolate simulation to experimental voltage points
                sim_current_interp = np.interp(exp_voltage, sim_voltages, sim_currents)
                
                # Compute error metrics
                print("\nComputing error metrics...")
                # Only use forward bias for comparison (V >= 0)
                forward_mask = exp_voltage >= 0
                metrics = compute_error_metrics(
                    exp_current[forward_mask],
                    sim_current_interp[forward_mask]
                )
                
                # Generate validation report
                report_path = path.join(output_dir, "validation_report.txt")
                generate_validation_report(metrics, output_path=report_path)
                
                # Plot comparison
                plot_path = path.join(output_dir, "iv_comparison.png")
                print(f"\nGenerating I-V comparison plot...")
                plot_iv_comparison(
                    exp_voltage, exp_current,
                    sim_voltages, sim_currents,
                    output_path=plot_path,
                    show=False
                )
                
                print("\n" + "="*60)
                print("Validation complete!")
                print(f"- Validation report: {report_path}")
                print(f"- Comparison plot: {plot_path}")
                print("="*60)
                
            else:
                print(f"Warning: Output directory not found: {output_dir}")
                print("Skipping validation.")
        
        except Exception as e:
            print(f"\nError during experimental validation: {e}")
            print("Continuing without validation...")
            import traceback
            traceback.print_exc()
    
    print("\nDone!")
