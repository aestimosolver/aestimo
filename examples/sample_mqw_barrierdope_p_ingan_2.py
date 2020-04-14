#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Input File Description:  Barrier doped AlGaAs/GaAs heterostructure.
# -------------------------------------------------------------------
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
computation_scheme = 7

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_h = 1
subnumber_e = 1
# APPLIED ELECTRIC FIELD
Fapplied =  0.0# (V/m)2.5e7/50e-9
vmax= 3.2
vmin= 0.0
Each_Step=0.05# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.2#nm
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
material =[[ 300.0, 'AlGaN', 0.05, 0.0, 5e18, 'p','b'],
            [ 2.0, 'InGaN', 0.2, 0.0, 1e-10,'n','w'],
            [ 5.0, 'GaN', 0.0, 0.0, 1e-10, 'n','b'],
            [ 2.0, 'InGaN', 0.2, 0.0, 1e-10,'n','w'],
            [ 200.0, 'AlGaN', 0.05, 0.0, 5e18, 'n','b']]

material1 =[[ 20.0, 'AlGaN', 0.3, 0.0, 5e16, 'p','b'],
            [ 10.0, 'GaN', 0.2, 0.0, 1e18,'p','w'],
            [ 20.0, 'AlGaN', 0.3, 0.0, 5e16, 'p','b']]
#----------------------------------------
import numpy as np
x_max = sum([layer[0] for layer in material])
#----------------------------------------
def round2int(x):
    return int(x+0.5)
n_max=round2int(x_max/gridfactor)
dop_profile=np.zeros(n_max)  
#----------------------------------------
Quantum_Regions=False
Quantum_Regions_boundary=np.zeros((2,2))
#----------------------------------------  
surface=np.zeros(2)
#----------------------------------------
inputfilename = "sample_mqw_barrierdope_p_ingan_2"
from os import path
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import sys
    sys.path.append(path.join(path.dirname(__file__), '..'))
    import aestimo_eh
    aestimo_eh.run_aestimo(input_obj)
