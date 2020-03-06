#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Input File Description:  Barrier doped AlGaAs/GaAs heterostructure.
# -------------------------------------------------------------------
# ----------------
# GENERAL SETTINGS
# ----------------
import numpy as np
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
computation_scheme = 2

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_h = 5
subnumber_e = 5
# APPLIED ELECTRIC FIELD
Fapplied =  0.0# (V/m)-20e8
vmax= 0.1
vmin= 0.0
Each_Step=0.05
# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------
contact=0.0
# GRID
# For 1D, z-axis is choosen
gridfactor = 0.1#nm
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


material =[[ 14.0, 'GaN', 0.0, 0.0,0.0, 'i','NA'],
           [ 30.0, 'AlGaN', 0.3, 0.0,0.0, 'i','NA'],
            [ 40.0, 'GaN', 0.0, 0.0,  0.0, 'i','NA']]

material1 =[[ 100.0, 'GaN', 0.0, 0.0,0.0, 'i','NA'],
           [ 17.0, 'AlGaN', 0.3, 0.0,0.0, 'i','NA'],
            [ 100.0, 'GaN', 0.0, 0.0,  0.0, 'i','NA']]

material2 =[[ 30.0, 'AlGaN', 0.3, 0.0,0.0, 'i','NA'],
            [ 40.0, 'GaN', 0.0, 0.0,  0.0, 'i','NA']]
#----------------------------------------
x_max = sum([layer[0] for layer in material])
def round2int(x):
    return int(x+0.5)
n_max=round2int(x_max/gridfactor)
#----------------------------------------
dop_profile=np.zeros(n_max)
#----------------------------------------
Quantum_Regions=True
#--------------------------------
Quantum_Regions_boundary=np.zeros((1,2))
Quantum_Regions_boundary[0,0]=2
Quantum_Regions_boundary[0,1]=55
#----------------------------------
#Quantum_Regions_boundary[0,0]=90
#Quantum_Regions_boundary[0,1]=127
#---------------------------------
#Quantum_Regions_boundary[0,0]=28
#Quantum_Regions_boundary[0,1]=50

#----------------------------------------
surface=np.zeros(2)
surface[0]=1.42
surface[1]=0.0
#----------------------------------------
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import aestimo_eh
    aestimo_eh.run_aestimo(input_obj)
