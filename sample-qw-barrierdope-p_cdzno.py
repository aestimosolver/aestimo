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
computation_scheme = 2

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_h = 4
subnumber_e = 2
# APPLIED ELECTRIC FIELD
Fapplied = 0.00#/50e-9 # (V/m)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.2 #nm
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
material =[[ 1.0, 'ZnO', 0.0, 0.0, 'p','b'],
            [ 2.0, 'ZnO', 0.0, 0.0, 'p','b'],
            [ 2.0, 'ZnO', 0.0, 5e17, 'p','b'],
            [ 3.0, 'CdZnO', 0.15, 0, 'p','w'],
            [ 2.0, 'ZnO', 0.0, 5e17, 'p','b'],
            [ 2.0, 'ZnO', 0.0, 0.0, 'p','b'],
            [ 1.0, 'ZnO', 0.0, 0.0, 'p','b']]
 


if __name__ == "__main__": #this code allows you to run the input file directly
    input_obj = vars()
    import aestimo_h
    aestimo_h.run_aestimo(input_obj)