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
# For now, there is only one computational scheme
# 1: Schrodinger-Poisson
computation_scheme = 1

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_e = 3

# APPLIED ELECTRIC FIELD
Fapplied = 0.1/50e-9 # (V/m)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
z_coordinate_begin = 0.0 # nm
z_coordinate_end = 50.0 # nm  
gridfactor = 0.1 #nm
maxgridpoints = 200000 #for controlling the size

# REGIONS
# Region input is a two-dimensional list input.
# An example:
# Si p-n diode. Firstly lets picturize the regional input.
#         | layer # | z-coordinate-begin | z-coordinate-end | Material | Alloy fraction | Doping(cm^-3) | n or p type |
# Layer 0 |     0   |      0.0           |        250.0     |   Si     |      0         |     1e16      |     n       |
# Layer 1 |     1   |      250.0         |        500.0     |   Si     |      0         |     1e16      |     p       |
#
# To input this list in Gallium, we use lists as:
material =[[0, 0.0, 10.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [1, 10.0, 15.0, 'AlGaAs', 0.3, 5e17, 'n'],
            [2, 15.0, 20.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [3, 20.0, 30.0, 'GaAs', 0, 0, 'n'],
            [4, 30.0, 35.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [5, 35.0, 40.0, 'AlGaAs', 0.3, 5e17, 'n'],
            [6, 40.0, 50.0, 'AlGaAs', 0.3, 0.0, 'n']]
 

