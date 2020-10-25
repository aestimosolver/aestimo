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
computation_scheme = 9
# QUANTUM
# Total subband number to be calculated for electrons
subnumber_h = 1
subnumber_e = 1
# APPLIED ELECTRIC FIELD
Fapplied = 0.#0.41348e8 (V/m)
vmax= 3.1
vmin= 0.0
Each_Step=0.5
# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.5#nm
maxgridpoints = 200000 #for controlling the size
#mat_type='Zincblende'
mat_type='Wurtzite'
# REGIONS
# Region input is a two-dimensional list input.
# An example:
# Si p-n diode. Firstly lets picturize the regional input.
#         | Thickness (nm) | Material | Alloy fraction | Doping(cm^-3) | n or p type |
# Layer 0 |      250.0     |   Si     |      0         |     1e16      |     n       |
# Layer 1 |      250.0     |   Si     |      0         |     1e16      |     p       |
# To input this list in Gallium, we use lists N:



material =[[ 30.0  , 'GaN'  , 0.3 , 0.0, 2e20, 'p','b'],           
           [ 100.0 , 'AlGaN', 0.14, 0.0, 1e20, 'p','b'],                      
           [ 50.0  , 'GaN'  , 0.3 , 0.0, 5e18, 'p','b'],
           [ 20.0  , 'AlGaN', 0.2 , 0.0, 1e19, 'p','b'],           
           [ 10.0  , 'AlGaN', 0.2 , 0.0, 7e16, 'n','b'],
           [ 4.0   , 'InGaN', 0.02, 0.0, 7e16, 'n','w'],
           [ 10.0  , 'AlGaN', 0.2 , 0.0, 7e16, 'n','b'],
           [ 4.0   , 'InGaN', 0.02, 0.0, 7e16, 'n','w'],
           [ 10.0  , 'AlGaN', 0.2 , 0.0, 7e16, 'n','b'],
           [ 4.0   , 'InGaN', 0.02, 0.0, 7e16, 'n','w'],
           [ 10.0  , 'AlGaN', 0.2 , 0.0, 7e16, 'n','b'],
           [ 50.0  , 'GaN'  , 0.3 , 0.0, 7e17, 'n','b'],
           [ 100.0 , 'AlGaN', 0.14, 0.0, 3e18, 'n','b'],
           [ 50.0  , 'InGaN', 0.1 , 0.0, 3e18, 'n','b'],
           [ 100.0 , 'GaN'  , 0.3 , 0.0, 3e18, 'n','b']]
#----------------------------------------
import numpy as np
x_max = sum([layer[0] for layer in material])
def round2int(x):
    return int(x+0.5)
n_max=round2int(x_max/gridfactor)
#----------------------------------------
#Doping profiles based on the LSS theory (ion implantation).
dop_n=np.zeros(n_max)
dop_p=np.zeros(n_max)
dop_profile=np.zeros(n_max)

"""   
import matplotlib.pyplot as pl
pl.plot(xaxis, dop_n*1e-6,'r',xaxis,dop_p*1e-6,'b')
#pl.plot(xaxis, dop_profile*1e-6,'k')
pl.xlabel('Position (m)')
pl.ylabel('electrons  and and holes concentrations (cm-3)' )
pl.title('electrons (red) and holes (blue)')
pl.grid(True)
khkhk
"""
#----------------------------------------
Quantum_Regions=False
Quantum_Regions_boundary=np.zeros((2,2))
#----------------------------------------
surface=np.zeros(2)
#surface[0]=-0.6
#----------------------------------------
inputfilename = "sample_qw_barrierdope_ingan_2"
from os import path
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import sys
    sys.path.append(path.join(path.dirname(__file__), '..'))
    import aestimo_eh
    aestimo_eh.run_aestimo(input_obj)
