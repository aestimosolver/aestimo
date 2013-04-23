#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------
# Input File Description:  Quantum well doped AlGaAs/GaAs heterostructure.
# ------------------------------------------------------------------------
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
subnumber_e = 2

# APPLIED ELECTRIC FIELD
Fapplied = 0.0 # (V/m)

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
material =[[0, 0.0, 20.0, 'AlGaAs', 0.2, 0, 'n'],
            [1, 20.0, 30.0, 'GaAs', 0, 2e18, 'n'],
            [1, 30.0, 50.0, 'AlGaAs', 0.2, 0, 'n']]
 

