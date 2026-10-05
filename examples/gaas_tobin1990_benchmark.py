#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tobin 1990 GaAs Benchmark Solar Cell Simulation (Mode 10 Newton-Raphson)
Reference: S. P. Tobin et al., IEEE Transactions on Electron Devices 37(2), 469-477 (1990).
Target Metrics: Jsc = 27.80 mA/cm2, Voc = 1.052 V, FF = 84.9%, Eff = 24.8% (1-sun AM1.5G, 298.15 K)
"""
import os, sys
import numpy as np

# General Simulation Settings
T = 298.15                         # Temperature [K]
computation_scheme = 10            # 10: Fully-Coupled Newton-Raphson
comp_scheme = 10
gridfactor = 5.0                   # Spatial grid spacing [nm]
dx = gridfactor * 1e-9             # [m]
maxgridpoints = 200000
mat_type = 'Zincblende'

# Voltage Sweep Parameters
vmin = 0.0                         # [V]
vmax = 1.14                        # [V]
Each_Step = 0.02                   # [V]
Fapplied = 0.0

# Solar & Physical Parameters
photovoltaic_mode = True
enable_polarization = False
G_optical = 6.759e+20              # Calibrated AM1.5G generation rate [cm^-3 s^-1]
device_area = 1.0                  # [cm^2]
device_area_m2 = device_area * 1e-4 # [m^2]
work_function_left = 5.2           # Anode / p-contact [eV]
work_function_right = 4.1          # Cathode / n-contact [eV]
surface = np.array([0.0, 0.0])     # Bulk/ohmic equilibrium potential reference
surface_recomb = (0, 0)
Rs = 0.0

# Subbands
subnumber_e = 5
subnumber_h = 5

# Quantum Regions
Quantum_Regions = False
Quantum_Regions_boundary = np.zeros((1, 2))

# Layer Structure: [thickness (nm), material, mole_x, mole_y, doping (cm^-3), type ('n'/'p'), role ('b'/'q')]
material = [
    [30.0,   'AlGaAs', 0.85, 0.0, 2.0e18, 'p', 'b'],  # Window
    [500.0,  'GaAs',   0.0,  0.0, 2.0e18, 'p', 'b'],  # Emitter
    [3000.0, 'GaAs',   0.0,  0.0, 2.0e17, 'n', 'q'],  # Base absorber
    [100.0,  'AlGaAs', 0.85, 0.0, 2.0e18, 'n', 'b'],  # BSF
]

inputfilename = 'gaas_tobin1990_benchmark'

if __name__ == '__main__':
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    import aestimo
    aestimo.output_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), 'gaas_tobin1990_benchmark_output'))
    aestimo.run_aestimo(vars())
