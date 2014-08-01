#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Input File Description:  Trying to replicate some results given in sec.3.5 
#               of Paul Harrison's book "Quantum Well's Wires and Dots".
#               Modelling a parabolic quantum well.
# -------------------------------------------------------------------

# ----------------
# HARRISON'S MATERIAL VALUES
# ----------------
import database
import config
import numpy as np
import aestimo
import database
meV2J = aestimo.meV2J # conversion factor
q = aestimo.q # electron charge

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

# Non-parabolic Dispersion Calculations for Fermi-Dirac
fermi_np_scheme = True

# QUANTUM
# Total subband number to be calculated for electrons
subnumber_e = 10

# APPLIED ELECTRIC FIELD
Fapplied = 0.00 # (V/m)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.01 #nm
maxgridpoints = 200000 #for controlling the size

# STRUCTURE
# a finite parabolic well between two barriers
a = 10.0 #nm #maximum width of parabolic well (at the barrier's energy level)
b = 10.0 #nm #width of the first barrier layer
b2 = 5.0 #nm #width of the first barrier layer
xmin = 0.0 #minimum Al alloy in structure
xmax = 10.0 #maximum Al alloy in structure


def alloy_profile(z):
    """function of alloy profile for finite parabolic well"""
    well = xmin + (z - (b+a/2.0))**2 / (a/2.0)**2 * (xmax - xmin)    
    return np.where(well<xmax,well,xmax) #cuts off parabolic well at barrier height

def bandstructure_profile(x):
    mat1 = model.material_property['AlAs']
    mat2 = model.material_property['GaAs']
    alloy = model.alloy_property['AlGaAs']
    return alloy['Band_offset']*(x*mat1['Eg'] + (1-x)* mat2['Eg']-alloy['Bowing_param']*x*(1-x))*aestimo.q # for electron. Joule
    #return alloy['Band_offset']*(x*mat1['Eg'] + (1-x)* mat2['Eg'])*aestimo.q # for electron. Joule

# We can't declare our structure in the usual way, instead we will have to create
# the structure arrays ourselves.


structure_param = {'Fapp': Fapplied,
                       'T': T,
                       'subnumber_e': subnumber_e,
                       'comp_scheme': computation_scheme,
                       'dx': gridfactor*1e-9, #grid in m
                       'maxgridpoints': maxgridpoints,
                       }
# Initialise structure class
model = aestimo.Structure(database,**structure_param)

## Create structure arrays -----------------------------------------------------

# Calculate the required number of grid points
model.x_max = sum([b,a,b2])*1e-9 #total thickness (m)
model.n_max = n_max = int(model.x_max/model.dx)
# Check on n_max
if model.n_max > model.maxgridpoints:
    aestimo.logger(" Grid number is exceeding the max number of %d", maxgridpoints)
    exit()
#
meff = model.material_property['GaAs']['m_e']*aestimo.m_e
epsGaAs = model.material_property['GaAs']['epsilonStatic']

model.z = np.arange(0.0,model.x_max,model.dx)*1e9 #nm

model.cb_meff = np.ones(n_max)*meff	#conduction band effective mass
model.cb_meff_alpha = np.zeros(n_max)   #non-parabolicity constant.
model.eps = np.ones(n_max)*epsGaAs	#dielectric constant
model.dop = np.ones(n_max)*0.0          #doping
model.fi = bandstructure_profile(alloy_profile(model.z))   #Bandstructure potential

## -----------------------------------------------------------------------------

#config.d_E = 1e-5*meV2J
#config.Estate_convergence_test = 3e-10*meV2J
result= aestimo.Poisson_Schrodinger(model)

#Plot QW representation
config.wavefunction_scalefactor = 5000
aestimo.QWplot(result)#,figno=None)
    