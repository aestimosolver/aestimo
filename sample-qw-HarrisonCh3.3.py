#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Input File Description:  Trying to replicate some results given in sec. 3.3 
#               of Paul Harrison's book "Quantum Well's Wires and Dots".
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
computation_scheme = 0

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_e = 1

# APPLIED ELECTRIC FIELD
Fapplied = 0.00 # (V/m)

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
#         | Thickness (nm) | Material | Alloy fraction | Doping(cm^-3) | n or p type |
# Layer 0 |      250.0     |   Si     |      0         |     1e16      |     n       |
# Layer 1 |      250.0     |   Si     |      0         |     1e16      |     p       |
#

# The two structures defined in the text:
material =[[ 15.0, 'AlGaAs', 0.2, 0.0, 'n'],
            [ 10.0, 'GaAs',     0, 1e16, 'n'],
            [ 15.0, 'AlGaAs', 0.2, 0.0, 'n']]




if __name__=="__main__":
    import matplotlib.pyplot as pl
    import numpy as np
    #import config
    import aestimo_numpy as aestimo
    import database
    
    # Harrison's initial examples don't take account of different effective masses
    # in the different materials. Hence I need to adjust the database
    database.alloyproperty['AlGaAs']['eps_b'] = 0.0
    database.alloyproperty['AlGaAs']['cb_mass_b'] = 0.0
    
    """
    structure_param = {'Fapp': Fapplied,
                       'T': T,
                       'subnumber_e': subnumber_e,
                       'comp_scheme': computation_scheme,
                       'dx': gridfactor*1e-9, #grid in m
                       'maxgridpoints': maxgridpoints,
                       'material': material,
                       }
    # Initialise structure class
    model = aestimo.Structure(database,**structure_param)
    """
    # Initialise structure class
    model = aestimo.StructureFrom(globals(),database)
    model.create_structure_arrays()
    
    result= aestimo.Poisson_Schrodinger(model)
    
    # Perform the calculation
    results = []
    gridfactors = [0.2,0.1,0.05,0.02,0.01] #nm
    for gridfactor in gridfactors:
        model.dx = gridfactor*1e-9
        model.create_structure_arrays() #updates our structure object
        result= aestimo.Poisson_Schrodinger(model)
        results.append(result.E_state[0])
    
    print 'E[0]'
    for gridfactor,E in zip(gridfactors,results): print gridfactor,E
    # Write the simulation results in files
    #aestimo.save_and_plot(result,model)
    
    #Plot QW representation
    aestimo.QWplot(result)#,figno=None)
    