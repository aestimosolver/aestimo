#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
# Input File Description:  Si p/n junction with experimental validation.
# ----------------------------------------------------------------------
# This example demonstrates how to validate drift-diffusion simulation
# results against experimental I-V data.
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
computation_scheme = 9

# QUANTUM
subnumber_h = 2
subnumber_e = 2

# VOLTAGE SWEEP PARAMETERS
Fapplied = 0.0  # Applied electric field (V/m)
vmax = 0.80     # Maximum voltage (V)
vmin = 0.0      # Minimum voltage (V)
Each_Step = 0.01  # Voltage step (V)

# --------------------------------
# EXPERIMENTAL CONFIGURATION
# --------------------------------
enable_experimental_validation = True
experimental_iv_file = "examples/experimental_data/si_pn_experimental_iv.csv"

# Robust path finding for the data file
if not os.path.exists(experimental_iv_file):
     script_dir = os.path.dirname(os.path.abspath(__file__))
     candidates = [
         os.path.join(script_dir, "experimental_data", "si_pn_experimental_iv.csv"),
         os.path.join(script_dir, "examples", "experimental_data", "si_pn_experimental_iv.csv")
     ]
     for c in candidates:
         if os.path.exists(c):
             experimental_iv_file = c
             break

# DYNAMIC PARAMETER LOADING
print(f"Reading metadata from: {experimental_iv_file}")
metadata = load_experimental_metadata(experimental_iv_file)

# Default values if metadata missing
device_area = metadata.get('area_cm2', 1e-4) # cm2
doping_p = metadata.get('doping_p', 1e18)
doping_n = metadata.get('doping_n', 1e18)

print(f"  - Device Area: {device_area:.2e} cm^2")
print(f"  - Doping (p): {doping_p:.2e} cm^-3")
print(f"  - Doping (n): {doping_n:.2e} cm^-3")

# Parasitic Resistances (Optimized)
Rs_ext = 2330.0  # Series Resistance (Ohm)
Rsh_ext = 5e7    # Shunt Resistance (Ohm)
rs_mode = "External (Fast)"

# --------------------------------
# DEVICE STRUCTURE
# --------------------------------
gridfactor = 10  # nm
maxgridpoints = 200000
mat_type = 'Zincblende'

material = [
    [2500.0, 'Si', 0.0, 0.0, doping_p, 'p', 'b'],  # p-side
    [2500.0, 'Si', 0.0, 0.0, doping_n, 'n', 'b']   # n-side
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
inputfilename = "sample_pn_with_experimental_validation"

# ----------------------------------------
# MAIN EXECUTION
# ----------------------------------------
if __name__ == "__main__":
    input_obj = vars()
    import aestimo
    
    print("="*60)
    print("Running Standardized Experimental Validation Workflow")
    print("="*60)
    
    # Run simulation - Validation is now triggered automatically by the hook inside run_aestimo
    results = aestimo.run_aestimo(input_obj)
    
    print("\nProcess finished. Please check the output directory for results.")
