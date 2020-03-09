#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------
# Input File Description:  Double Quantum well doped AlGaAs/GaAs heterostructure.
# ------------------------------------------------------------------------
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
computation_scheme = 2
# Non-parabolic effective mass function
# 0: no energy dependence
# 1: Nelson's effective 2-band model
# 2: k.p model from Vurgaftman's 2001 paper
meff_method = 2

# Non-parabolic Dispersion Calculations for Fermi-Dirac
fermi_np_scheme = True

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_e = 1
subnumber_h = 3
# APPLIED ELECTRIC FIELD
Fapplied = 0.0 # (V/m)
vmax= 1.7
vmin= 0.0
Each_Step=0.05
mat_type='Zincblende'
# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.5 #nm
maxgridpoints = 200000 #for controlling the size
# REGIONS
# Region input is a two-dimensional list input.
# An example:
# Si p-n diode. Firstly lets picturize the regional input.
#         | Thickness (nm)  | Material | Alloy fraction | Doping(cm^-3) | n or p type |
# Layer 0 |       250.0     |   Si     |      0         |     1e16      |     n       |
# Layer 1 |       250.0     |   Si     |      0         |     1e16      |     p       |
#
# To input this list in Gallium, we use lists as:

material =[[ 300.0, 'AlGaAs', 0.3, 0.3, 1e17, 'p','b'],
            [3.0, 'GaAs', 0, 0.3, 0.0, 'i','w'],
            [20.0, 'AlGaAs', 0.3, 0.3, 0.0, 'i','b'],
            [3.0, 'GaAs', 0, 0.3, 0.0, 'i','w'],
            [20.0, 'AlGaAs', 0.3, 0.3, 0.0, 'i','b'],
            [300.0, 'AlGaAs', 0.3, 0.3, 1e17, 'n','b']]

#----------------------------------------
#Doping profiles based on the LSS theory (ion implantation).
import numpy as np
x_max = sum([layer[0] for layer in material])
def round2int(x):
    return int(x+0.5)
n_max=round2int(x_max/gridfactor)
#----------------------------------------
dop_n=np.zeros(n_max)
dop_p=np.zeros(n_max)
dop_profile=np.zeros(n_max)
xaxis = np.arange(0,n_max)*gridfactor#[nm]
Q_n=2e12#implant dose [1/cm2]
Rp_n=86#projected range Rp [nm]
Delta_Rp_n=44#projected straggle Delta Rp [nm]
Q_p=1e11#implant dose [1/cm2]
Rp_p=75#projected range Rp [nm]
Delta_Rp_p=20#projected straggle Delta Rp [nm]
from math import sqrt, exp
def Lss_profile_dop(x,Q,Delta_Rp,Rp):   
    return Q/(sqrt(2*np.pi)*Delta_Rp*1e-7)*exp(-(x-Rp)**2/(2*Delta_Rp**2))
def Lss_profile_dop_diff(x,Q,Delta_Rp,Rp):   
    return Q/(2*sqrt(np.pi)*Delta_Rp*1e-7)*exp(-(x-Rp)**2/(4*Delta_Rp**2))
for i in range(n_max):   
    dop_n[i]=Lss_profile_dop(xaxis[n_max-1-i],Q_n,Delta_Rp_n,Rp_n)*1e6
    dop_p[i]=-Lss_profile_dop(xaxis[n_max-1-i],Q_p,Delta_Rp_p,Rp_p)*1e6
    #dop_profile[i]=dop_n[i]+dop_p[i] 
#----------------------------------------
Quantum_Regions=False
Quantum_Regions_boundary=np.zeros((2,2))
#----------------------------------------  
surface=np.zeros(2)
#surface[0]=-0.6
#----------------------------------------
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import aestimo_eh
    aestimo_eh.run_aestimo(input_obj)
