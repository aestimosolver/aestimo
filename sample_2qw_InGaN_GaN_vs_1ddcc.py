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
# 7: Schrodinger-Poisson-Drift_Diffusion

computation_scheme = 2

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_h = 1
subnumber_e = 1
# APPLIED ELECTRIC FIELD
Fapplied =  0.0# (V/m)-20e8
Vapplied=3.1
# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------
contact=0.0
# GRID
# For 1D, z-axis is choosen
gridfactor = 0.5#nm
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
            [ 3.0, 'InGaN', 0.2, 0.0, 0.01e17, 'n','w'],
            [ 15.0, 'GaN',  0.0, 0.0, 0.01e17, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.01e17, 'n','w'],
            [ 300.0, 'GaN', 0.0, 0.0, 5e18, 'n','b']]

material11  =[[ 15.0, 'GaN',  0.0, 0.0, 2e19, 'p','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.01e17, 'n','w'],
            [ 15.0, 'GaN',  0.0, 0.0, 0.01e17, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.01e17, 'n','w'],
            [ 15.0, 'GaN', 0.0, 0.0, 5e18, 'n','b']]
material1 =[[ 200.0, 'GaN',  0.0, 0.0, 2e19, 'p','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 1.0e17, 'n','w'],
            [ 300.0, 'GaN', 0.0, 0.0, 5e18, 'n','b']]

material2 =[[ 200.0, 'ZnO',  0.0, 0.0, 2e19, 'p','b'],
            [ 3.0, 'CdZnO', 0.2, 0.0, 1.0e17, 'n','w'],
            [ 15.0, 'ZnO',  0.0, 0.0, 1.0e17, 'n','b'],
            [ 3.0, 'CdZnO', 0.2, 0.0, 1.0e17, 'n','w'],
            [ 300.0, 'ZnO', 0.0, 0.0, 5e18, 'n','b']]

material3 =[[ 200.0, 'GaN',  0.0, 0.0, 1.2e18, 'p','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.7e17, 'n','w'],
            [ 15.0, 'GaN',  0.0, 0.0, 0.7e17, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.7e17, 'n','w'],
            [ 300.0, 'GaN', 0.0, 0.0, 4.2e17, 'n','b']]

material4  =[[ 200.0, 'GaN',  0.0, 0.0, 0.0, 'p','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.0, 'n','w'],
            [ 15.0, 'GaN',  0.0, 0.0, 0.0, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.0, 'n','w'],
            [ 300.0, 'GaN', 0.0, 0.0, 0.0, 'n','b']]

material5 =[[ 200.0, 'GaN',  0.0, 0.0, 1.8e19, 'p','b'],
            [ 3.0, 'AlGaInN', 0.0, 0.2, 0.8e17, 'n','w'],
            [ 15.0, 'GaN',  0.0, 0.0, 0.8e17, 'n','b'],
            [ 3.0, 'AlGaInN', 0.0, 0.2, 0.8e17, 'n','w'],
            [ 300.0, 'GaN', 0.0, 0.0, 4.8e18, 'n','b']]


material6 =[[ 300.0, 'GaN', 0.3, 0.0, 4.2e18, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.7e17, 'n','w'],
            [ 15.0,  'GaN', 0.3, 0.0, 0.7e17, 'n','b'],
            [ 3.0, 'InGaN', 0.2, 0.0, 0.7e17, 'n','w'],
            [ 200.0, 'GaN', 0.3, 0.0, 1.2e19, 'p','b']]


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

time1 = time.time()
print("total running time=",time1-time0)
