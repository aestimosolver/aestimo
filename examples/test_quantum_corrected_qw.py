#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import sys
from os import path

# Setup paths to find aeslibs
root_dir = os.getcwd()
sys.path.append(root_dir)

# ----------------
# GENERAL SETTINGS
# ----------------
T = 300.0 
computation_scheme = 9 # SP-DD (Gummel-Newton)

# QUANTUM
subnumber_h = 3
subnumber_e = 3

# VOLTAGE SWEEP
Fapplied = 0.0
vmin = 0.0
vmax = 0.2
Each_Step = 0.05

# DEVICE STRUCTURE
gridfactor = 1.0 # nm
maxgridpoints = 200000
mat_type = 'Zincblende'

# [Thickness (nm), Material, Alloy fraction x, Alloy fraction y, Doping(cm^-3), type, role]
material = [
    [250.0, 'AlGaAs', 0.3, 0.0, 1e18, 'p', 'b'],
    [50.0,  'AlGaAs', 0.3, 0.0, 0.0,  'n', 'b'],
    [15.0,  'GaAs',   0.0, 0.0, 0.0,  'n', 'w'], # GaAs Well
    [50.0,  'AlGaAs', 0.3, 0.0, 0.0,  'n', 'b'],
    [250.0, 'AlGaAs', 0.3, 0.0, 1e18, 'n', 'b']
]

# ---------------------------------------- 
x_max = sum([layer[0] for layer in material])
def round2int(x): return int(x + 0.5)
n_max = round2int(x_max / gridfactor)
dop_profile = np.zeros(n_max)

# ENABLE QUANTUM CORRECTED DD
Quantum_Regions = True
# Region 0: GaAs Well (starts at 250+50 = 300nm, lasts 15nm)
Quantum_Regions_boundary = np.array([[300.0, 315.0]])

surface = np.zeros(2)
inputfilename = "test_quantum_corrected_qw"

if __name__ == "__main__":
    # Ensure we use local aestimo
    sys.path.insert(0, root_dir)
    import aestimo
    input_obj = vars()
    print("="*60)
    print("Test Case: Simple Quantum Well with Quantum Corrected DD")
    print("="*60)
    aestimo.run_aestimo(input_obj)
    print("\nProcess finished.")
