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
# 2: Schrodinger-Poisson with non-parabolicity
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
material =[[ 20.0, 'AlGaAs', 0.2, 0, 'n'],
            [10.0, 'GaAs', 0, 2e18, 'n'],
            [20.0, 'AlGaAs', 0.2, 0, 'n']]
 

