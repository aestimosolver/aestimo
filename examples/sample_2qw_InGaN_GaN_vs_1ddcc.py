#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Input File Description:  Barrier doped AlGaAs/GaAs heterostructure.
# -------------------------------------------------------------------
# ----------------
# GENERAL SETTINGS
# ----------------
import time
time0 = time.time() # timing audit
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

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_h = 1
subnumber_e = 1
# APPLIED ELECTRIC FIELD
Fapplied =  0.0# (V/m)-20e8
vmax= 3.3
vmin= 3.0
Each_Step=0.15
# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------
contact=0.0
# GRID
# For 1D, z-axis is choosen
gridfactor = 0.25#nm
maxgridpoints = 200000 #for controlling the size
mat_type='Wurtzite'
# REGIONS
# Region input is a two-dimensional list input.
# An example:
# Si p-n diode. Firstly lets picturize the regional input.
#         | Thickness (nm) | Material | Alloy fraction | Doping(cm^-3) | n or p type |
# Layer 0 |      250.0     |   Si     |      0         |     1e16      |     n       |
# Layer 1 |      250.0     |   Si     |      0         |     1e16      |     p       |
#
# To input this list in Gallium, we use lists as:

material  =[[ 200.0, 'GaN',  0.0, 0.0, 2e19, 'p','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.0e17, 'n','w'],
            [ 15.0, 'GaN',  0.0, 0.0, 0.0e17, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.0e17, 'n','w'],
            [ 300.0, 'GaN', 0.0, 0.0, 5e18, 'n','b']]


material1  =[[ 15.0, 'GaN',  0.0, 0.0, 2e19, 'p','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.01e17, 'n','w'],
            [ 14.0, 'GaN',  0.0, 0.0, 0.01e17, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.01e17, 'n','w'],
            [ 15.0, 'GaN', 0.0, 0.0, 5e18, 'n','b']]

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
Quantum_Regions_boundary[0,0]=10
Quantum_Regions_boundary[0,1]=25
Quantum_Regions_boundary[1,0]=26
Quantum_Regions_boundary[1,1]=38
#----------------------------------------
surface=np.zeros(2)
surface[0]=0.0
surface[1]=0.0
#----------------------------------------
inputfilename = "sample_2qw_InGaN_GaN_vs_1ddcc"
#this code allows you to run the input file directly
from os import path
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import sys
    sys.path.append(path.join(path.dirname(__file__), '..'))
    import aestimo_eh
    aestimo_eh.run_aestimo(input_obj)

time1 = time.time()
print("total running time=",time1-time0)
