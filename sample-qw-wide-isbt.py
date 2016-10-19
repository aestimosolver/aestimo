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
# Total subband number to be calculated for electrons (for aestimo_numpy_h)
subnumber_h = 1

# APPLIED ELECTRIC FIELD
Fapplied = 0.00/50e-9 # (V/m)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.1 #nm
maxgridpoints = 200000 #for controlling the size

# DOPING
Nd = 2e18

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
            [ 5.0, 'AlGaAs', 0.3, Nd/2.0, 'n'],
            [ 5.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [ 30.0, 'GaAs', 0, 0, 'n'],
            [ 5.0, 'AlGaAs', 0.3, 0.0, 'n'],
            [ 5.0, 'AlGaAs', 0.3, Nd/2.0, 'n'],
            [ 10.0, 'AlGaAs', 0.3, 0.0, 'n']]
 


if __name__=="__main__":
    import config
    import database
    import aestimo
    import intersubband_optical_transitions as isbt
    import numpy as np
    import os
    import time
    
    logger = aestimo.logger
    
    # Initialise structure class
    model = aestimo.StructureFrom(vars(),database)
    
    if True: #recalculate QW states
        # Perform the calculation
        result = aestimo.Poisson_Schrodinger(model)
        
        time4 = time.time() #timing audit
        logger.info("total running time (inc. loading libraries) %g s",(time4 - aestimo.time0))
        logger.info("total running time (exc. loading libraries) %g s",(time4 - aestimo.time1))
        
        # Write the simulation results in files
        fig1,fig2,fig3 = aestimo.save_and_plot(result,model)
        logger.info("Simulation is finished.")
    else: #load previously calculated results from output directory
        result = aestimo.load_results()
    
    #Set thickness of effective medium
    Lperiod = sum([layer[0] for layer in model.material])*1e-9 #m
    
    # set dielectric constants
    case = 2
    if case==1: #scalar dielectric constants
        eps_b = 12.90
        eps_z = 12.90
    
    elif case==2: #z-dependent dielectric constants
        eps_b = 10.364
        eps_gaas = 10.364 # @ 16um
        eps_algaas = 8.2067
        eps_z = isbt.eps_background_GaAs(model,eps_gaas,eps_algaas)
    
    # Linewidth
    def linewidth(freq): return 0.1*freq #define linewidth in THz
    
    # Optical Intersubband Transitions
    transitions_table,(hdr,units) = isbt.transitions(result,Lperiod,eps_z,linewidth)
    
    isbt.print_levels(result)
    isbt.print_transitions(transitions_table,hdr,units)
    isbt.print_multiplasmon_transitions(*isbt.calc_wR_multiplasmon(result,transitions_table,eps_z))
    
    fig4 = isbt.plotting_absorption(model,result,transitions_table,eps_b,eps_z,linewidth)

