#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
# Input File Description:  InGaN p/n junction with experimental validation.
# ----------------------------------------------------------------------
# This project simulates an InGaN-based p-n homojunction with 
# high indium content (x=0.57) and validates it against experimental data.
# ----------------------------------------------------------------------

import numpy as np
import os
import sys
from os import path

# Setup paths to find aeslibs
root_dir = os.getcwd()
if os.path.exists(os.path.join(root_dir, 'aeslibs')):
    sys.path.append(root_dir)
else:
    sys.path.append(path.join(path.dirname(__file__), '..'))

# ----------------
# GENERAL SETTINGS
# ----------------
T = 300.0  # SIMULATION PARAMETERS
computation_scheme = 7 # Standard DD (Highest stability for TAT)
subnumber_h = 2
subnumber_e = 2
Each_Step = 0.04

# VOLTAGE SWEEP
vmin = 0.0
vmax = 0.8

# --------------------------------
# EXPERIMENTAL CONFIGURATION
# --------------------------------
enable_experimental_validation = True
experimental_iv_file = "examples/experimental_data/ingan_pn_experimental_iv.csv"

# Device Area from experimental file (5e-4 cm^2)
device_area = 5e-4 

# Trap-Assisted Tunneling (TAT) - Hurkx Model
# Strong TAT to lift low-bias current and match high n.
tat_field = 1.0e7 

# LIFETIMES (Set high to allow TAT to dominate the slope)
TAUN0 = 1.0e-6 # 1 us
TAUP0 = 1.0e-6 # 1 us

# --------------------------------
# DEVICE STRUCTURE
# --------------------------------
gridfactor = 5.0  # nm - 100 points
maxgridpoints = 200000
mat_type = 'Wurtzite' 
enable_polarization = True 

# InGaN (x=0.57) Doping
doping_p = 1e16   
doping_n = 2e17   

G_optical = 0.0 # DARK MATCHING

material = [
    [250.0, 'InGaN', 0.57, 0.0, doping_p, 'p', 'b'], # p-side
    [250.0, 'InGaN', 0.57, 0.0, doping_n, 'n', 'b']  # n-side
]

# ---------------------------------------- 
# STANDARD SETUP
# ----------------------------------------
x_max = sum([layer[0] for layer in material])
def round2int(x): return int(x + 0.5)
n_max = round2int(x_max / gridfactor)
tat_field = 1.0e7  # Realistic TAT field for Nitrides
Quantum_Regions = False  
inputfilename = "ingan_calibrated"

# ----------------------------------------
# MAIN EXECUTION
# ----------------------------------------
if __name__ == "__main__":
    import aestimo
    import numpy as np
    import os
    import shutil
    from aeslibs.experimental_validation import (
        load_experimental_data, apply_parasitic_resistances, 
        compute_error_metrics, generate_validation_report, 
        plot_iv_comparison, load_current_from_avcurr
    )

    # 1. Config Aestimo Parameters
    class AttrDict(dict):
        def __init__(self, *args, **kwargs):
            super(AttrDict, self).__init__(*args, **kwargs)
            self.__dict__ = self

    # --- Best Physical Parameters ---
    tat_field = 5.0e6  # Strong field enhancement
    tau = 5.0e-7       # 500 ns baseline
    
    config_dict = {
        'material': [
            [40.0, 'InGaN', 0.57, 0.0, 1e16, 'p', 'b'],
            [60.0, 'InGaN', 0.57, 0.0, 2e17, 'n', 'b']
        ],
        'T': T,
        'gridfactor': 1.0, 
        'mat_type': mat_type,
        'computation_scheme': computation_scheme,
        'vmin': vmin,
        'vmax': vmax,
        'Each_Step': Each_Step,
        'TAUN0': tau,
        'TAUP0': tau,
        'tat_field': tat_field,
        'enable_polarization': enable_polarization,
        'output_directory': 'results'
    }
    input_obj = AttrDict(config_dict)
    input_obj.inputfilename = "ingan_final_calib"
    
    # 2. Execute Aestimo
    print("="*60)
    print(f"Starting Final Calibration Run (TAT={tat_field:.1e}, Tau={tau:.1e})")
    print("="*60)
    
    aestimo.run_aestimo(input_obj)
    actual_out = aestimo.output_directory
    
    # 3. Validation and Parasite Optimization
    exp_v, exp_i = load_experimental_data(experimental_iv_file)
    v_sim, i_sim_A = load_current_from_avcurr(actual_out, device_area_cm2=device_area)
    
    print("\n--- Optimizing Parasites ---")
    best_r2 = -np.inf
    best_p = {}
    best_c = (None, None)
    
    for rs in np.linspace(1, 150, 40):
        for rsh in np.logspace(3, 5, 40):
            v_f, i_f = apply_parasitic_resistances(v_sim, i_sim_A, Rs=rs, Rsh=rsh)
            i_interp = np.interp(exp_v, v_f, i_f)
            m = compute_error_metrics(exp_i, i_interp)
            if m['r2'] > best_r2:
                best_r2 = m['r2']
                best_p = {'rs': rs, 'rsh': rsh}
                best_c = (v_f, i_f)
    
    print(f"Optimal Params: Rs={best_p['rs']:.1f}, Rsh={best_p['rsh']:.1e}, BEST R2: {best_r2:.4f}")

    # 4. Final Reporting
    v_final, i_final = best_c
    metrics = compute_error_metrics(exp_i, np.interp(exp_v, v_final, i_final))
    
    report_path = os.path.join(actual_out, "validation_report.txt")
    generate_validation_report(metrics, report_path)
    
    plot_path = os.path.join(actual_out, "iv_comparison.png")
    plot_iv_comparison(exp_v, exp_i, v_final, i_final, output_path=plot_path)

    print(f"\nSUCCESS: Results saved in {actual_out}/")
