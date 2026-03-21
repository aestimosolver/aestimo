#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
# Input File Description: InGaAs p/n junction with experimental validation.
# ----------------------------------------------------------------------
# This project simulates an In(0.53)Ga(0.47)As p-n homojunction 
# and validates it against synthetic experimental data.
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
T = 300.0  # Kelvin
computation_scheme = 9 # SP-DD (Gummel-Newton)

# QUANTUM (Not strictly needed for bulk p-n but kept for consistency)
subnumber_h = 2
subnumber_e = 2

# VOLTAGE SWEEP
# Matching experimental range: 0.0V to 1.0V (Focus on forward bias)
Fapplied = 0.0
vmin = 0.0
vmax = 1.0
Each_Step = 0.05

# --------------------------------
# EXPERIMENTAL CONFIGURATION
# --------------------------------
enable_experimental_validation = True
experimental_iv_file = "examples/experimental_data/ingaas_pn_experimental_iv.csv"

# Device Area from experimental file (1e-4 cm^2)
device_area = 1e-4 

# Parasitic Resistances (Initial guesses, might need optimization)
Rs_ext = 10.0   # Series Resistance (Ohm)
Rsh_ext = 1.0e8 # Shunt Resistance (Ohm) - Optimized
rs_mode = "External (Fast)"

# --------------------------------
# DEVICE STRUCTURE
# --------------------------------
gridfactor = 5.0  # nm
maxgridpoints = 200000
mat_type = 'Zincblende'

# InGaAs (x=0.53) Doping
doping_p = 1e17   
doping_n = 1e17   

material = [
    [250.0, 'InGaAs', 0.53, 0.0, doping_p, 'p', 'b'], # p-side
    [250.0, 'InGaAs', 0.53, 0.0, doping_n, 'n', 'b']  # n-side
]

# ---------------------------------------- 
# STANDARD SETUP
# ----------------------------------------
x_max = sum([layer[0] for layer in material])
def round2int(x): return int(x + 0.5)
n_max = round2int(x_max / gridfactor)
dop_profile = np.zeros(n_max)
Quantum_Regions = False
Quantum_Regions_boundary = np.zeros((2, 2))
surface = np.zeros(2)
inputfilename = "sample_pn_with_experimental_validation_ingaas"

# ----------------------------------------
# MAIN EXECUTION
# ----------------------------------------
if __name__ == "__main__":
    input_obj = vars()
    import aestimo
    
    print("="*60)
    print("InGaAs p-n Junction: Standardized Validation Workflow")
    print("="*60)
    print(f"Material: In(0.53)Ga(0.47)As")
    print(f"Doping: {doping_p:e} (p) / {doping_n:e} (n)")
    print(f"Area: {device_area:e} cm^2")
    print("="*60)
    
    # Run simulation - Validation hook integrated in run_aestimo
    results = aestimo.run_aestimo(input_obj)
    
    print("\nProcess finished. Please check the output directory for results.")
