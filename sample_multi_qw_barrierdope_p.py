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
computation_scheme = 7

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_h = 1
subnumber_e = 1
# APPLIED ELECTRIC FIELD
Fapplied = 0.00#/50e-9 # (V/m)
Vapplied=1.8# (V)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.2 #nm
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
material =[[ 200.0, 'AlGaAs', 0.3, 0.0, 5e17, 'p','b'],
            [ 2.0, 'AlGaAs', 0.3, 0.0, 5e17, 'p','b'],
            [ 2.0, 'AlGaAs', 0.3, 0.0, 0.0, 'p','b'],
            [ 3.0, 'GaAs', 0.0, 0.0, 0.0, 'p','w'],
            [ 8.0, 'AlGaAs', 0.3, 0.0, 0.0, 'p','b'],
            [ 3.0, 'GaAs', 0.0, 0.0, 0.0, 'p','w'],
            [ 8.0, 'AlGaAs', 0.3, 0.0, 0.0, 'p','b'],
            [ 3.0, 'GaAs', 0.0, 0.0, 0.0, 'p','w'],
            [ 8.0, 'AlGaAs', 0.3, 0.0, 0.0, 'p','b'],
            [ 3.0, 'GaAs', 0.0, 0.0, 0.0, 'p','w'],
            [ 2.0, 'AlGaAs', 0.3, 0.0, 0.0, 'p','b'],
            [ 2.0, 'AlGaAs', 0.3, 0.0, 5e17, 'n','b'],
            [ 200.0, 'AlGaAs', 0.3, 0.0, 5e17, 'n','b']]
 
import numpy as np
x_max = sum([layer[0] for layer in material])
def round2int(x):
    return int(x+0.5)
n_max=round2int(x_max/gridfactor)
dop_profile=np.zeros(n_max)  
surface=np.zeros(2)
#surface[0]=-3.0
#surface[1]=-0.3
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import aestimo_eh
    aestimo_eh.run_aestimo(input_obj)
