
import sys
import os
import numpy as np

# Add current directory to path
sys.path.append(os.getcwd())

from aestimo import run_aestimo

class InputObject:
    T = 300.0
    computation_scheme = 2 
    subnumber_h = 3
    subnumber_e = 3
    gridfactor = 0.5 # nm
    maxgridpoints = 200000 
    mat_type= 'Zincblende'
    dx = 1e-10
    material = [[50.0, 'GaAs', 0.0, 0.0, 0.0, 'n', 'b'],
                [10.0, 'InGaAs', 0.15, 0.0, 0.0, 'n', 'w'],
                [50.0, 'GaAs', 0.0, 0.0, 0.0, 'n', 'b']]
    # Field parameters
    Fapplied = 0.0 
    vmax = 1.79
    vmin = 0.0
    Each_Step = 0.05

    surface = np.zeros(2)
    surface[0] = 0.0
    surface[1] = 0.6
    
    Quantum_Regions = False
    Quantum_Regions_boundary = np.zeros((2,2))
    
    total_thickness = sum(row[0] for row in material) * 1e-9
    n_max = int(total_thickness / dx)
    dop_profile = np.zeros(n_max)
    current_idx = 0
    for row in material:
        th = row[0] * 1e-9
        doping_val = row[4]
        if row[5] == 'p':
             doping_val = -doping_val
        steps = int(th / dx)
        end_idx = min(current_idx + steps, n_max)
        dop_profile[current_idx:end_idx] = doping_val
        current_idx = end_idx
    
InputObject.__file__ = os.path.abspath("repro_sim.py")

try:
    print("Running simulation...")
    run_aestimo(InputObject)
    print("Simulation success!")
except Exception as e:
    print(f"Simulation failed: {e}")
