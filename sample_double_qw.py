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
subnumber_e = 6
subnumber_h = 4
# APPLIED ELECTRIC FIELD
Fapplied = 0.0 # (V/m)
mat_type='Zincblende'
# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.1 #nm
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
material =[[ 20.0, 'AlGaAs', 0.3, 1e14, 'n','b'],
            [10.0, 'GaAs', 0, 2e16, 'n','w'],
            [5.0, 'AlGaAs', 0.3, 1e14, 'n','b'],
            [10.0, 'GaAs', 0, 2e16, 'n','w'],
            [20.0, 'AlGaAs', 0.3, 1e14, 'n','b']]
 
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import aestimo
    aestimo.run_aestimo(input_obj)
"""    
if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import aestimo_h
    aestimo_h.run_aestimo(input_obj)
"""