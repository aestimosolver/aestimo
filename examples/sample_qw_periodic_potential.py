#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Input File Description:  Barrier doped AlGaAs/GaAs heterostructure.
# -------------------------------------------------------------------
# ----------------
# GENERAL SETTINGS
# ----------------

# TEMPERATURE
T = 60.0 #Kelvin

# COMPUTATIONAL SCHEME
# 0: Schrodinger
# 1: Schrodinger + nonparabolicity
# 2: Schrodinger-Poisson
# 3: Schrodinger-Poisson with nonparabolicity
# 4: Schrodinger-Exchange interaction
# 5: Schrodinger-Poisson + Exchange interaction
# 6: Schrodinger-Poisson + Exchange interaction with nonparabolicity
computation_scheme = 6

# Non-parabolic effective mass function
# 0: no energy dependence
# 1: Nelson's effective 2-band model
# 2: k.p model from Vurgaftman's 2001 paper
meff_method = 2

# Non-parabolic Dispersion Calculations for Fermi-Dirac
fermi_np_scheme = True

# (optional) Dopant Poisson field boundary condition
# (only applies to electric field from dopants)
# 0: Efield = 0 at boundaries
# 1: Periodic Potential 
boundary_scheme = 1

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_e = 3
# Total subband number to be calculated for electrons (for aestimo_numpy_h)
subnumber_h = 1

# APPLIED ELECTRIC FIELD
Fapplied = 0.0/50e-9 # (V/m)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.1 #nm
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
material =[[ 10.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [ 5.0, 'AlGaAs', 0.3, 5e17, 'n'],
            [ 5.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [ 11.0, 'GaAs', 0, 0, 'n'],
            [ 5.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [ 5.0, 'AlGaAs', 0.3, 0e17, 'n'],
            [ 10.0, 'AlGaAs', 0.3, 0.0, 'n']]
 

inputfilename = "sample_qw_periodic_potential"
from os import path
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import sys
    sys.path.append(path.join(path.dirname(__file__), '..'))
    import aestimo
    aestimo.run_aestimo(input_obj)
    
