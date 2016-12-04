#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Input File Description:  Trying to replicate some results given in sec. 3.3 
#               of Paul Harrison's book "Quantum Well's Wires and Dots".
# -------------------------------------------------------------------

# ----------------
# HARRISON'S MATERIAL VALUES
# ----------------
import database

# Nb. Harrison's initial examples don't take account of different effective masses
# in the different materials.

database.materialproperty = {
    'GaAs':{
        'm_e':0.067,
        'm_e_alpha':0.0,
        'epsilonStatic':12.90,
        'Eg':0.0, #1.426, # set the energy scale origin to be at the GaAs condution band
        'Band_offset':0.67,
        },
    'AlAs':{
        'm_e':0.067, # normally Harrison would be using 0.15
        'm_e_alpha':0.0,
        'epsilonStatic':10.06,
        'Eg':2.673-1.426, #2.673, # set the energy scale origin to be at the GaAs condution band
        'Band_offset':0.67,
        },
    }

database.alloyproperty = {
    'AlGaAs':{
        'Bowing_param':0.0,
        'Band_offset':0.67,
        'Material1':'AlAs',
        'Material2':'GaAs',
        'm_e_alpha':0.0,
        },
    }

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

# Non-parabolic effective mass function
# 0: no energy dependence
# 1: Nelson's effective 2-band model
# 2: k.p model from Vurgaftman's 2001 paper
meff_method = 0

# Non-parabolic Dispersion Calculations for Fermi-Dirac
fermi_np_scheme = True

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
mat_type='Zincblende'
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
            [ 10.0, 'GaAs',     0, 1e10, 'n'],
            [ 15.0, 'AlGaAs', 0.2, 0.0, 'n']]




if __name__=="__main__":
    import matplotlib.pyplot as pl
    import numpy as np
    #import config
    import aestimo
    import database
    
    logger = aestimo.logger
    
    # Initialise structure class
    model = aestimo.StructureFrom(globals(),database)
    model.create_structure_arrays()
    
    result= aestimo.Poisson_Schrodinger(model)
    
    # Perform the calculation using different gridfactors
    results = []
    gridfactors = [0.2,0.1,0.05,0.02,0.01] #nm
    for gridfactor in gridfactors:
        model.dx = gridfactor*1e-9
        model.create_structure_arrays() #updates our structure object
        result= aestimo.Poisson_Schrodinger(model)
        results.append(result.E_state[0]) #-np.min(model.fi)*aestimo.J2meV)
    
    logger.info('gridfactor (nm),E[0]')
    for gridfactor,E in zip(gridfactors,results): logger.info('%g,%g',gridfactor,E)
    
    # Perform the calculation using different barrier widths
    results2 = []
    barrier_widths = [5,6,7,8,9,10,11,12,13,14,15,16,18,20,25,30,35] #nm
    for barrier in barrier_widths:
        material[0][0]=barrier;material[2][0]=barrier
        model = aestimo.StructureFrom(globals(),database)
        model.create_structure_arrays()
        result= aestimo.Poisson_Schrodinger(model)
        results2.append(result.E_state[0]) #-np.min(model.fi)*aestimo.J2meV)
    
    logger.info('barrier (nm),E[0]')
    for barrier,E in zip(barrier_widths,results2): logger.info('%g,%g',barrier,E)
        
    f1 = pl.figure()
    ax1 = f1.add_subplot(111)
    ax1.set_xlim(0,40)
    ax1.plot(barrier_widths,results2,'-o')
    
    # Write the simulation results in files
    #aestimo.save_and_plot(result,model)
    
    #Plot QW representation
    aestimo.QWplot(result)#,figno=None)
    
