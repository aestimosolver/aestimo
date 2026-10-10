#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
# Input File Description: In(0.57)Ga(0.43)N p/n homojunction diode with experimental validation.
# ----------------------------------------------------------------------
# This benchmark validates drift-diffusion simulation against experimental I-V
# data for high-In InGaN homojunction (ingan_pn_experimental_iv.csv).
#
# Simulation uses:
# - Mode 10: Fully-Coupled Newton-Raphson drift-diffusion solver
# - Material: In(0.57)Ga(0.43)N (Wurtzite, Eg ~ 1.2-1.4 eV)
# - Doping: p-side 1e16 cm^-3 (Mg-compensated), n-side 5e17 cm^-3
# - Parasitics: Rs = 25 Ohm, Rsh = 2000 Ohm, Area = 5e-4 cm^2
# ----------------------------------------------------------------------

import numpy as np
import os
import sys
from os import path

# Ensure local workspace takes priority over site-packages
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
if os.getcwd() not in sys.path:
    sys.path.insert(0, os.getcwd())

try:
    from aeslibs.experimental_validation import load_experimental_metadata
except ImportError:
    def load_experimental_metadata(f): return {}

# ----------------
# GENERAL SETTINGS
# ----------------
T = 300.0  # Kelvin
computation_scheme = 10  # Fully-Coupled Newton-Raphson
comp_scheme = 10

subnumber_h = 2
subnumber_e = 2

# Voltage sweep
Fapplied = 0.0
vmin = 0.0
vmax = 0.80
Each_Step = 0.02

# --------------------------------
# EXPERIMENTAL CONFIGURATION
# --------------------------------
enable_experimental_validation = True
experimental_iv_file = os.path.join(repo_root, "examples", "experimental_data", "ingan_pn_experimental_iv.csv")
if not os.path.exists(experimental_iv_file):
    experimental_iv_file = os.path.join(script_dir, "experimental_data", "ingan_pn_experimental_iv.csv")

metadata = load_experimental_metadata(experimental_iv_file)
device_area = metadata.get('area_cm2', 5.0e-4)
doping_p = metadata.get('doping_p', 1.0e+16)
doping_n = metadata.get('doping_n', 5.0e+17)

Rs_ext = 50.0    # Series Resistance (Ohm) - Physical nitride contact/bulk resistance
Rsh_ext = 1.5e5  # Shunt Resistance (Ohm) - 150 kOhm defect-assisted leakage shunt

# --------------------------------
# REGIONAL SETTINGS
# --------------------------------
gridfactor = 2.0  # nm
maxgridpoints = 200000
mat_type = 'Wurtzite'

# Structure: In_0.57Ga_0.43N homojunction
material = [
    [200.0, 'InGaN', 0.57, 0.0, doping_p, 'p', 'b'],  # p-side
    [400.0, 'InGaN', 0.57, 0.0, doping_n, 'n', 'b']   # n-side
]

x_max = sum([layer[0] for layer in material])

def round2int(x):
    return int(x + 0.5)

n_max = round2int(x_max / gridfactor)
dop_profile = np.zeros(n_max)
Quantum_Regions = False
Quantum_Regions_boundary = np.zeros((2, 2))
surface = np.zeros(2)
inputfilename = "sample_pn_with_experimental_validation_ingan"

if __name__ == "__main__":
    input_obj = vars()
    import aestimo
    print("="*60)
    print("Running InGaN p-n Junction Simulation (Mode 10)")
    print("="*60)
    print(f"Device structure: {len(material)} layers, total length: {x_max} nm")
    print(f"Voltage range: {vmin}V to {vmax}V (step: {Each_Step}V)")
    print(f"Area: {device_area:.2e} cm^2, Rs: {Rs_ext} Ohm, Rsh: {Rsh_ext} Ohm")
    print("="*60)

    results = aestimo.run_aestimo(input_obj)
    print("\nSimulation completed!")

    if enable_experimental_validation:
        print("\n" + "="*60)
        print("EXPERIMENTAL VALIDATION (InGaN Homojunction)")
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
            exp_voltage, exp_current = load_experimental_data(experimental_iv_file)
            output_dir = f"{inputfilename}_output"

            if os.path.exists(output_dir):
                calc_v_int, calc_i_int = load_current_from_avcurr(output_dir, device_area_cm2=device_area)
                calc_v, calc_i = apply_parasitic_resistances(calc_v_int, calc_i_int, Rs=Rs_ext, Rsh=Rsh_ext)

                sim_current_interp = np.interp(exp_voltage, calc_v, calc_i)
                mask = exp_voltage >= 0.1
                metrics = compute_error_metrics(exp_current[mask], sim_current_interp[mask])

                report_path = path.join(output_dir, "validation_report.txt")
                generate_validation_report(metrics, output_path=report_path)

                plot_path = path.join(output_dir, "iv_comparison.png")
                plot_iv_comparison(exp_voltage, exp_current, calc_v, calc_i, output_path=plot_path, show=False)

                print(f"- Log-RMSE: {metrics.get('log_rmse', 'N/A'):.4f}")
                print(f"- Validation report: {report_path}")
                print(f"- Comparison plot: {plot_path}")
        except Exception as e:
            print(f"Error during InGaN validation: {e}")
            import traceback
            traceback.print_exc()

    print("\nDone!")
