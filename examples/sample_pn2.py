#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
# Input File Description:  Si p/n junction.
# ----------------------------------------------------------------------
# ----------------
# GENERAL SETTINGS
# ----------------

# TEMPERATURE
T = 300.0 #Kelvin

# COMPUTATIONAL SCHEME
# 0: Schrodinger
# 1: Schrodinger + nonparabolicity
# 2: Schrodinger-Poisson
# 3: Schrodinger-Poisson with nonparabolicity
# 4: Schrodinger-Exchange interaction
# 5: Schrodinger-Poisson + Exchange interaction
# 6: Schrodinger-Poisson + Exchange interaction with nonparabolicity
# 7: Schrodinger-Poisson-Drift_Diffusion (Schrodinger solved with poisson then  poisson and DD)
# 8: Schrodinger-Poisson-Drift_Diffusion (Schrodinger solved with poisson and DD)
# 9: Schrodinger-Poisson-Drift_Diffusion (Schrodinger solved with poisson and DD) using Gummel & Newton map
# 10: Schrodinger-Poisson under testing using Newton iteration (will replace scheme 2)
computation_scheme = 2

# Non-parabolic effective mass function
# 0: no energy dependence
# 1: Nelson's effective 2-band model
# 2: k.p model from Vurgaftman's 2001 paper
#meff_method = 0

# Non-parabolic Dispersion Calculations for Fermi-Dirac
fermi_np_scheme = True

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_e = 1
subnumber_h = 1
# APPLIED ELECTRIC FIELD
Fapplied =  0.0# (V/m)2.5e7/50e-9
vmax= 2.0
vmin= 0.0
Each_Step=0.05# For 1D, z-axis is choosen
gridfactor = 1
maxgridpoints = 200000 #for controlling the size
mat_type='Zincblende'
# REGIONS
# Region input is a two-dimensional list input.
# An example:
# Si p-n diode. Firstly lets picturize the regional input.
#         | Thickness (nm) | Material | Alloy fraction | Doping(cm^-3) | n or p type |
# Layer 0 |      250.0     |   Si     |      0         |     1e16      |     n       |
# Layer 1 |      250.0     |   Si     |      0         |     1e16      |     p       |
#
# To input this list in Gallium, we use lists as:
material =[[500, 'Si', 0.1, 0.0, 1e+19, 'p','b'],
            [500, 'Si', 0.1, 0.0, 1e+19, 'n','b']]


material1 =[[350, 'InGaAs', 0.1, 0.0, 1e+16, 'p','b'],
            [15, 'InAs', 0.0, 0.0, 1e+16, 'n','w'],
            [600, 'InGaAs', 0.1, 0.0, 1e+17, 'n','b']] 
#---------------------------------------- 
import numpy as np
x_max = sum([layer[0] for layer in material])
def round2int(x):
    return int(x+0.5)
n_max=round2int(x_max/gridfactor)
#----------------------------------------
dop_profile=np.zeros(n_max)
#----------------------------------------
Quantum_Regions=False
Quantum_Regions_boundary=np.zeros((2,2))
#----------------------------------------  
surface=np.zeros(2)
#----------------------------------------
inputfilename = "sample_pn2"
from os import path
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import sys
    sys.path.append(path.join(path.dirname(__file__), '..'))
    import aestimo_eh
    aestimo_eh.run_aestimo(input_obj)