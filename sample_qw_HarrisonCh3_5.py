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
#import database
meV2J = aestimo.meV2J # conversion factor
q = aestimo.q # electron charge

logger=aestimo.logger

# Nb. Harrison's initial examples don't take account of different effective masses
# in the different materials.

materialproperty = {
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

alloyproperty = {
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
s0 = {}


# TEMPERATURE
s0['T'] = 300.0 #Kelvin

# COMPUTATIONAL SCHEME
# 0: Schrodinger
# 1: Schrodinger + nonparabolicity
# 2: Schrodinger-Poisson
# 3: Schrodinger-Poisson with nonparabolicity
# 4: Schrodinger-Exchange interaction
# 5: Schrodinger-Poisson + Exchange interaction
# 6: Schrodinger-Poisson + Exchange interaction with nonparabolicity
s0['comp_scheme']= 0

# Non-parabolic effective mass function
# 0: no energy dependence
# 1: Nelson's effective 2-band model
# 2: k.p model from Vurgaftman's 2001 paper
s0['meff_method']= 0

# Non-parabolic Dispersion Calculations for Fermi-Dirac
s0['fermi_np_scheme'] = True

# QUANTUM
# Total subband number to be calculated for electrons
s0['subnumber_e'] = 10

# APPLIED ELECTRIC FIELD
s0['Fapp'] = 0.00 # (V/m)

# --------------------------------
# REGIONAL SETTINGS FOR SIMULATION
# --------------------------------

# GRID
# For 1D, z-axis is choosen
gridfactor = 0.01 #nm
s0['dx'] = gridfactor*1e-9
maxgridpoints = 200000 #for controlling the size
mat_type='Zincblind'
## Create structure arrays -----------------------------------------------------

# We can't declare our structure in the usual way using the StructureFrom class
# instead we will have to create and the structure arrays ourselves and work with
# the parent Structure class.

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
    mat1 = materialproperty['AlAs']
    mat2 = materialproperty['GaAs']
    alloy = alloyproperty['AlGaAs']
    return alloy['Band_offset']*(x*mat1['Eg'] + (1-x)* mat2['Eg']-alloy['Bowing_param']*x*(1-x))*aestimo.q # for electron. Joule
    #return alloy['Band_offset']*(x*mat1['Eg'] + (1-x)* mat2['Eg'])*aestimo.q # for electron. Joule

# Calculate the required number of grid points
s0['x_max'] = sum([b,a,b2])*1e-9 #total thickness (m)
s0['n_max'] = n_max = int(s0['x_max']/s0['dx'])
# Check on n_max
if s0['n_max']> maxgridpoints:
    aestimo.logger(" Grid number is exceeding the max number of %d", maxgridpoints)
    exit()
#
meff = materialproperty['GaAs']['m_e']*aestimo.m_e
epsGaAs = materialproperty['GaAs']['epsilonStatic']

z = np.arange(0.0,s0['x_max'],s0['dx'])*1e9 #nm

s0['cb_meff'] = np.ones(n_max)*meff	#conduction band effective mass
#s0['cb_meff_alpha'] = np.zeros(n_max)   #non-parabolicity constant.
s0['eps'] = np.ones(n_max)*epsGaAs	#dielectric constant
s0['dop'] = np.ones(n_max)*0.0          #doping
s0['fi'] = bandstructure_profile(alloy_profile(z))   #Bandstructure potential

#logger.info(s0.items())
#model = aestimo.Structure(T,Fapp,subnumber_e,dx,n_max, #parameters
                 #fi,eps,dop,cb_meff, #arrays
                 #comp_scheme,meff_scheme,fermi_np_scheme, #model choices
                 #cb_meff_alpha=None,Eg=None,Ep=None,F=None,delta_S0=None, #optional arrays
                 #**kwargs)
model = aestimo.Structure(**s0)

## -----------------------------------------------------------------------------

if __name__ == "__main__":
    #config.d_E = 1e-5*meV2J
    #config.Estate_convergence_test = 3e-10*meV2J
    result= aestimo.Poisson_Schrodinger(model)

    #Plot QW representation
    config.wavefunction_scalefactor = 5000
    fig = aestimo.QWplot(result)#,figno=None)

