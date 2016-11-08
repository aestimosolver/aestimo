#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This is the aestimo calculator of conduction band structure for quantum wells.

aestimo.py can be used as a script or a libary.

To use as a script, define the simulation in a python file (see the sample-*.py 
files for examples of the parameters needed) and then run aestimo on the command
line as 
  ./aestimo.py -i <input file>
Since we are abusing the python module system, the input 'file' needs to be 
importable by the aestimo script. Alternatively, define the input file in the
config module using the inputfilename parameter.

To use aestimo as a library, first create an instance of the Structure class; 
this can be more conveniently done using the StructureFrom class which builds 
the arrays describing a structure from a simple list format that describes the
structure's layers. See the class docstrings for details on the required 
parameters.

The simplest way to then calculate the structure is to use the Poisson_Schrodinger()
function. See the run_aestimo() function's source code for details on presenting or 
saving the results of this calculation using the returned object.

Calculations can be sped up by compiling the cythonised version of the psi_at_inf*
functions. This can be done using the setup_cython.py module.
"""
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Copyright (C) 2013-2016 Sefer Bora Lisesivdin and Aestimo group

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. See ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt .

    For the list of contributors, see ~/AUTHORS 
"""
__version__='1.1.0'

import time
time0 = time.time() # timing audit
#from scipy.optimize import fsolve
import matplotlib.pyplot as pl
import numpy as np
alen = np.alen
import os
import config,database
from math import log,exp
# --------------------------------------
import logging
logger = logging.getLogger('aestimo')
#File
hdlr = logging.FileHandler(config.logfile)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
#stderr
ch = logging.StreamHandler()
formatter2 = logging.Formatter('%(message)s')
ch.setFormatter(formatter2)
logger.addHandler(ch)
# LOG level can be INFO, WARNING, ERROR
logger.setLevel(logging.INFO)
# --------------------------------------

#Defining constants and material parameters
q = 1.602176e-19 #C
kb = 1.3806504e-23 #J/K
nii = 0.0
hbar = 1.054588757e-34
m_e= 9.1093826E-31 #kg
pi=np.pi
eps0= 8.8541878176e-12 #F/m

J2meV=1e3/q #Joules to meV
meV2J=1e-3*q #meV to Joules

time1 = time.time() # timing audit
#logger.info("Aestimo is starting...")

# Input Class
# -------------------------------------

def round2int(x):
    """int is sensitive to floating point numerical errors near whole numbers,
    this moves the discontinuity to the half interval. It is also equivalent
    to the normal rules for rounding positive numbers."""
    # int(x + (x>0) -0.5) # round2int for positive and negative numbers
    return int(x+0.5)

class Structure():
    def __init__(self,T,Fapp,subnumber_e,dx,n_max, #parameters
                 fi,eps,dop,cb_meff, #arrays
                 comp_scheme,meff_method,fermi_np_scheme, #model choices
                 cb_meff_alpha=None,Eg=None,Ep=None,F=None,delta_S0=None, #optional arrays
                 **kwargs):
        """This class holds details on the structure to be simulated.
        database is the module containing the material properties. Then
        this class should have the following attributes set
          Fapp - float - applied field (Vm**-1)
          T - float - Temperature (K)
          subnumber_e - int- number of subbands to look for.
          comp_scheme - int - computing scheme (see below)
          meff_method - int - choose effective mass function to model non-parabolicity (see below)
          fermi_np_scheme - bool - include nonparabolicity in calculation of Fermi level
        
          dx - float - grid step size (m)
          n_max - int - number of grid points
        
        #arrays - float/double (len = n_max)
          fi - Bandstructure potential (J) 
          eps - dielectric constant (including eps0)
          dop - doping distribution (m**-3)
          cb_meff - conduction band effective mass (kg)
        
        #optional arrays - float/double (len = n_max)
          cb_meff_alpha - non-parabolicity constant. Used for Nelson's empirical 2-band model)
          Eg - band gap energy (Used for k.p model found in Vurgaftman)
          Ep - (Used for k.p model found in Vurgaftman)
          F - (Used for k.p model found in Vurgaftman)
          delta_S0 - spin split-off energy (Used for k.p model found in Vurgaftman)
                
        COMPUTATIONAL SCHEMES (comp_scheme)
        0: Schrodinger
        1: Schrodinger + nonparabolicity
        2: Schrodinger-Poisson
        3: Schrodinger-Poisson with nonparabolicity
        4: Schrodinger-Exchange interaction
        5: Schrodinger-Poisson + Exchange interaction
        6: Schrodinger-Poisson + Exchange interaction with nonparabolicity
        
        EFFECTIVE MASS MODELS (meff_method) 
                        - for when nonparabolicity is modelled (see comp_scheme)
        0: no energy dependence
        1: Nelson's effective 2-band model
        2: k.p model from Vurgaftman's 2001 paper
        """
        #value attributes        
        self.T = T
        self.Fapp = Fapp
        self.subnumber_e = subnumber_e
        self.comp_scheme = comp_scheme
        self.fermi_np_scheme = fermi_np_scheme
        self.dx = dx
        self.n_max = n_max
        self.x_max = dx*n_max
        
        self.fi = fi
        self.eps = eps
        self.dop = dop
        self.cb_meff = cb_meff
        
        #Nelson's 2-band nonparabolicity model
        self.cb_meff_alpha= cb_meff_alpha
        #Vurgaftman's k.p nonparabolicity model
        self.Eg= Eg
        self.Ep= Ep
        self.F= F
        self.delta_S0= delta_S0
        
        # setting any extra parameters provided with initialisation
        for key,value in kwargs.items():
            setattr(self,key,value)
        
        #choosing effective mass function for non-parabolicity calculations
        self.meff_method = meff_method
        if self.meff_method == 1:
            self.cb_meff_E = self.cb_meff_E1
        elif self.meff_method == 2:
            self.cb_meff_E = self.cb_meff_E2
        
    def cb_meff_E(self,E,fi):
        """returns an array for the structure giving the effective mass for a particular
        energy. This version simply returns the conduction-band edge effective mass without
        any energy dependence. This default method should be overwritten in the class'
        __init__() method in order to model non-parabolicity.
        E - energy (J)
        fi - bandstructure potential (J) (numpy array)
        """
        return self.cb_meff*np.ones_like(E)
        
    def cb_meff_E1(self,E,fi):
        """returns an array for the structure giving the effective mass for a particular
        energy. This method implements Nelson's empirical 2-band model of non-parabolicity.
        This method needs to be aliased to the cb_meff_E() method in order for aestimo to
        use it.
        E - energy (J)
        fi - bandstructure potential (J) (numpy array)
        """
        return self.cb_meff*(1.0 + self._cb_meff_alpha*(E-fi))
        
    def cb_meff_E2(self,E,fi):
        """returns an array for the structure giving the effective mass for a particular
        energy. This uses the non-parabolicity calculation as given by Vurgaftman's 2001 
        paper on semiconductor properties. This method needs to be aliased to the 
        cb_meff_E() method in order for aestimo to use it.
        E - energy (J)
        fi - bandstructure potential (J) (numpy array)
        """
        EeV = (E - fi)/q
        return m_e/((1+2*self.F) + self.Ep/3.0*(2.0/(EeV+self.Eg) + 1.0/(EeV+self.Eg+self.delta_SO)))
 

class AttrDict(dict):
    """turns a dictionary into an object with attribute style lookups"""
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class StructureFrom(Structure):
    def __init__(self,inputfile,database):
        """database is a module or object containing the semiconductor material
        properties. Normally this is just aestimo's database.py module.
        
        inputfile is a dict or object with the following required parameters or
        attributes. These are exactly the parameters normally defined in the 
        input files (note that the parameter names are not exactly the same 
        as those for the Structure class) :
        
          Fapplied - applied field (Vm**-1) (float)
          T - Temperature (K) (float)
          subnumber_e - number of subbands to look for. (int)
          computation_scheme - computing scheme (see below) (int)
          meff_method - choose effective mass function to model non-parabolicity (see below) (int)
          fermi_np_scheme - include nonparabolicity in calculation of Fermi level (bool)

          gridfactor - grid step size (m) (float)
          maxgridpoints - number of grid points (int)
        
          material - a list describing a structure's layers. Each layer is a list
                     containing the following fields:
               thickness -  (nm)
               material type - any defined in the database ie. {GaAs,AlAs,AlGaAs,InGaAs...}  
               alloy fraction - value between 0.0 and 1.0
               doping (cm^-3) - density of dopants in the layer
               doping type - whether the dopants are n-type or p-type
              i.e.
              [[ 20.0, 'AlGaAs', 0.3, 1e17, 'n'],
               [ 50.0,   'GaAs',   0,    0, 'n']]
        

        COMPUTATIONAL SCHEMES (comp_scheme)
        0: Schrodinger
        1: Schrodinger + nonparabolicity
        2: Schrodinger-Poisson
        3: Schrodinger-Poisson with nonparabolicity
        4: Schrodinger-Exchange interaction
        5: Schrodinger-Poisson + Exchange interaction
        6: Schrodinger-Poisson + Exchange interaction with nonparabolicity

        EFFECTIVE MASS MODELS (meff_method) 
                        - for when nonparabolicity is modelled (see comp_scheme)
        0: no energy dependence
        1: Nelson's effective 2-band model
        2: k.p model from Vurgaftman's 2001 paper
        
        Nb. If changes are made to an instance's input parameters after 
        initialisation, the method create_structure_arrays() should be run to
        update the instance to the new parameters.
        """
        if type(inputfile)==dict:
            inputfile=AttrDict(inputfile)            
        # Parameters for simulation
        self.Fapp = inputfile.Fapplied
        self.T = inputfile.T
        self.subnumber_e = inputfile.subnumber_e
        self.comp_scheme = inputfile.computation_scheme
        self.fermi_np_scheme = inputfile.fermi_np_scheme
        self.dx = inputfile.gridfactor*1e-9 #grid in m
        self.maxgridpoints = inputfile.maxgridpoints
        
        # Loading material list
        self.material = inputfile.material
        totallayer = alen(self.material)
        logger.info("Total layer number: %s",totallayer)
        
        # Calculate the required number of grid points
        self.x_max = sum([layer[0] for layer in self.material])*1e-9 #total thickness (m)
        self.n_max = int(self.x_max/self.dx)
        
        # Check on n_max
        max_val = inputfile.maxgridpoints
        if self.n_max > max_val:
            logger.error(" Grid number is exceeding the max number of %d", max_val)
            exit()
                
        # Loading materials database
        self.material_property = database.materialproperty
        totalmaterial = alen(self.material_property)
        
        self.alloy_property = database.alloyproperty
        totalalloy = alen(self.alloy_property)
        
        logger.info("Total number of materials in database: %d",(totalmaterial+totalalloy))
        
        #choosing effective mass function for non-parabolicity calculations
        self.meff_method = inputfile.meff_method
        if self.meff_method == 1:
            self.cb_meff_E = self.cb_meff_E1
        elif self.meff_method == 2:
            self.cb_meff_E = self.cb_meff_E2
        
        # Initialise arrays
        
        #cb_meff #conduction band effective mass (array, len n_max)
        #fi #Bandstructure potential (array, len n_max)
        #eps #dielectric constant (array, len n_max)
        #dop #doping distribution (array, len n_max)
        self.create_structure_arrays()
        
    def create_structure_arrays(self):
        """ initialise arrays/lists for structure"""
        # Calculate the required number of grid points
        self.x_max = sum([layer[0] for layer in self.material])*1e-9 #total thickness (m)
        n_max = round2int(self.x_max/self.dx)
        # Check on n_max
        maxgridpoints = self.maxgridpoints
        if n_max > maxgridpoints:
            logger.error("Grid number is exceeding the max number of %d",maxgridpoints)
            exit()
        #
        self.n_max = n_max
        dx =self.dx
        material_property = self.material_property
        alloy_property = self.alloy_property
        
        cb_meff = np.zeros(n_max)   #conduction band effective mass
        cb_meff_alpha = np.zeros(n_max) #non-parabolicity constant.
        F = np.zeros(n_max)             #? (?)
        Eg = np.zeros(n_max)            #bandgap energy (?)
        delta_SO = np.zeros(n_max)      #spin-split off energy (?)
        Ep = np.zeros(n_max)            #? (?)
        fi = np.zeros(n_max)        #Bandstructure potential
        eps =np.zeros(n_max)        #dielectric constant
        dop = np.zeros(n_max)           #doping
        
        position = 0.0 # keeping in nanometres (to minimise errors)
        for layer in self.material:
            startindex = round2int(position*1e-9/dx)
            position += layer[0] # update position to end of the layer
            finishindex = round2int(position*1e-9/dx)
            #
            matType = layer[1]
            
            if matType in material_property:
                matprops = material_property[matType]
                fi[startindex:finishindex] = matprops['Band_offset']*matprops['Eg']*q #Joule
                eps[startindex:finishindex] = matprops['epsilonStatic']*eps0
                cb_meff[startindex:finishindex] = matprops['m_e']*m_e
                if self.meff_method == 1:
                    cb_meff_alpha[startindex:finishindex] = matprops['m_e_alpha']
                elif self.meff_method == 2:
                    F[startindex:finishindex] = matprops['F']
                    Eg[startindex:finishindex] = matprops['Eg'] #eV
                    delta_SO[startindex:finishindex] = matprops['delta']
                    Ep[startindex:finishindex] = matprops['Ep']

                
            elif matType in alloy_property:
                alloyprops = alloy_property[matType]
                mat1 = material_property[alloyprops['Material1']]
                mat2 = material_property[alloyprops['Material2']]
                x = layer[2] #alloy ratio
                fi[startindex:finishindex] = alloyprops['Band_offset']*(x*mat1['Eg'] + (1-x)* mat2['Eg']-alloyprops['Bowing_param']*x*(1-x))*q # for electron. Joule
                eps[startindex:finishindex] = (x*mat1['epsilonStatic'] + (1-x)* mat2['epsilonStatic'] )*eps0
                cb_meff_alloy = x*mat1['m_e'] + (1-x)* mat2['m_e']
                cb_meff[startindex:finishindex] = cb_meff_alloy*m_e
                if self.meff_method == 1:
                    cb_meff_alpha[startindex:finishindex] = alloyprops['m_e_alpha']*(mat2['m_e']/cb_meff_alloy) #non-parabolicity constant for alloy. THIS CALCULATION IS MOSTLY WRONG. MUST BE CONTROLLED. SBL
                elif self.meff_method == 2:
                    F[startindex:finishindex] = x*mat1['F'] + (1-x)* mat2['F']
                    Eg[startindex:finishindex] = x*mat1['Eg'] + (1-x)* mat2['Eg']-alloyprops['Bowing_param']*x*(1-x) #eV
                    delta_SO[startindex:finishindex] = x*mat1['delta'] + (1-x)* mat2['delta']-alloyprops['delta_bowing_param']*x*(1-x)
                    Ep[startindex:finishindex] = x*mat1['Ep'] + (1-x)* mat2['Ep']
            
            #doping
            if layer[4] == 'n':  
                chargedensity = layer[3]*1e6 #charge density in m**-3 (conversion from cm**-3)
            elif layer[4] == 'p': 
                chargedensity = -layer[3]*1e6 #charge density in m**-3 (conversion from cm**-3)
            else:
                chargedensity = 0.0
            
            dop[startindex:finishindex] = chargedensity
        
        self.fi = fi
        self.eps = eps
        self.dop = dop
        self.cb_meff = cb_meff
        if self.meff_method == 1: #only define optional arrays if we think that we need them
            self._cb_meff_alpha = cb_meff_alpha
        elif self.meff_method == 2:
            self.F = F
            self.Eg = Eg
            self.delta_SO = delta_SO
            self.Ep = Ep

        #return fi,cb_meff,eps,dop

# DO NOT EDIT UNDER HERE FOR PARAMETERS
# --------------------------------------

#Vegard's law for alloys
def vegard(first,second,mole):
    return first*mole+second*(1-mole)

# This function returns the value of the wavefunction (psi)
# at +infinity for a given value of the energy.  The solution
# to the energy occurs for psi(+infinity)=0.


# FUNCTIONS for SHOOTING ------------------

def psi_at_inf(E,fis,cb_meff,n_max,dx):
    """Shooting method for heterostructure as given in Harrison's book.
    This works much faster on lists than on numpy arrays for some reason 
    (probably type-casting related).
    E - energy (Joules)
    fis - potential energy (J)
    cb_meff - effective mass in conduction band (array)
    model - instance of Structure class
    n_max - number of points in arrays describing structure wrt z-axis
    dx - step size of distance quantisation (metres)
    """
    fis = fis.tolist() #lists are faster than numpy arrays for loops
    cb_meff = cb_meff.tolist() #lists are faster than numpy arrays for loops
    c0 = 2*(dx/hbar)**2
    # boundary conditions
    psi0 = 0.0                 
    psi1 = 1.0
    psi2 = None
    for j in range(1,n_max-1,1):
        # Last potential not used
        c1=2.0/(cb_meff[j]+cb_meff[j-1])
        c2=2.0/(cb_meff[j]+cb_meff[j+1])
        psi2=((c0*(fis[j]-E)+c2+c1)*psi1-c1*psi0)/c2
        psi0=psi1
        psi1=psi2
    return psi2
    
def psi_at_inf1(E,fis,model,n_max,dx):
    """shooting method with parabolic dispersions (energy independent effective mass).
    E - energy (Joules)
    fis - potential energy (J)
    cb_meff_f is a function that returns an array of effective mass at a given Energy E
    model - instance of Structure
    n_max - number of points in arrays describing structure wrt z-axis
    dx - step size of distance quantisation (metres)
    """
    return psi_at_inf(E,fis,model.cb_meff,n_max,dx)

def psi_at_inf2(E,fis,model,n_max,dx):
    """shooting method with non-parabolicity.
    E - energy (Joules)
    fis - potential energy (J)
    cb_meff_f is a function that returns an array of effective mass at a given Energy E
    model - instance of Structure
    n_max - number of points in arrays describing structure wrt z-axis
    dx - step size of distance quantisation (metres)
    """
    cb_meff = model.cb_meff_E(E,fis) # energy dependent mass
    return psi_at_inf(E,fis,cb_meff,n_max,dx)

try:
    from psi_at_inf_cython import psi_at_inf_numpy
    
    def psi_at_inf1_cython(E,fis,model,n_max,dx):
        return psi_at_inf_numpy(E,fis,model.cb_meff,n_max,dx)

    def psi_at_inf2_cython(E,fis,model,n_max,dx):
        """shooting method with non-parabolicity"""
        cb_meff_E = model.cb_meff_E(E,fis) # energy dependent mass
        return psi_at_inf_numpy(E,fis,cb_meff_E,n_max,dx)
        
    logger.info("using psi_at_inf_cython module")
except ImportError:
    logger.warning("psi_at_inf_cython module not found")


#nb. function was much slower when fi is a numpy array than a python list.
def calc_E_state(numlevels,fi,model,energyx0): # delta_E,d_E
    """Finds the Eigen-energies of any bound states of the chosen potential.
    numlevels - number of levels to find
    fi - Potential energy (Joules)
    model - any object with attributes: 
        cb_meff - array of effective mass (len n_max)
        n_max - length of arrays
        dx - step size (metres)
    energyx0 - minimum energy for starting subband search (Joules)"""
    # Shooting method parameters for Schrödinger Equation solution
    delta_E = config.delta_E #0.5*meV2J #Energy step (Joules) for initial search. Initial delta_E is 1 meV.
    d_E = config.d_E #1e-5*meV2J #Energy step (Joules) for Newton-Raphson method when improving the precision of the energy of a found level.
    Estate_convergence_test = config.Estate_convergence_test #1e-9*meV2J
    #
    E_state=[0.0]*numlevels #Energies of subbands (meV)
    cb_meff = model.cb_meff # effective mass of electrons in conduction band (kg)
    energyx = float(energyx0) #starting energy for subband search (Joules) + floats are faster than numpy.float64
    n_max = model.n_max
    dx = model.dx
    
    #choose shooting function
    if config.use_cython == True:
        if model.comp_scheme in (1,3,6): #then include non-parabolicity calculation
            psi_at_inf = psi_at_inf2_cython
        else:
            psi_at_inf = psi_at_inf1_cython
    else:
        if model.comp_scheme in (1,3,6): #then include non-parabolicity calculation
            psi_at_inf = psi_at_inf2
        else:
            psi_at_inf = psi_at_inf1
    
    #print 'energyx', energyx,type(energyx)
    #print 'cb_meff', cb_meff[0:10], type(cb_meff), type(cb_meff[0])
    #print 'n_max', n_max, type(n_max)
    #print 'fi', fi[0:10], type(fi), type(fi[0])
    #print 'dx', dx, type(dx)
    #exit()
    
    for i in range(0,numlevels,1):  
        #increment energy-search for f(x)=0
        y2=psi_at_inf(energyx,fi,model,n_max,dx)
        while True:
            y1=y2
            energyx += delta_E
            y2=psi_at_inf(energyx,fi,model,n_max,dx)
            if y1*y2 < 0:
                break
        # improve estimate using midpoint rule
        energyx -= abs(y2)/(abs(y1)+abs(y2))*delta_E
        #implement Newton-Raphson method
        while True:
            y = psi_at_inf(energyx,fi,model,n_max,dx)
            dy = (psi_at_inf(energyx+d_E,fi,model,n_max,dx)- psi_at_inf(energyx-d_E,fi,model,n_max,dx))/(2.0*d_E)
            energyx -= y/dy
            if abs(y/dy) < Estate_convergence_test:
                break
        E_state[i]=energyx*J2meV
        # clears x from solution
        energyx += delta_E # finish for i-th state.
    return E_state

# FUNCTIONS for ENVELOPE FUNCTION WAVEFUNCTION--------------------------------
def wf(E,fis,model):
    """This function returns the value of the wavefunction (psi)
    at +infinity for a given value of the energy.  The solution
    to the energy occurs for psi(+infinity)=0.
    psi[3] wavefunction at z-delta_z, z and z+delta_z 
    i index
    
    E - eigen-energy of state (Joules)
    fis - Potential energy of system (Joules)
    model - an object with atributes:
        cb_meff - array of effective mass (len n_max)
        n_max - length of arrays
        dx - step size (metres)"""
    #choosing effective mass function
    if model.comp_scheme in (1,3,6): #non-parabolicity calculation
        cb_meff_E = model.cb_meff_E(E,fis).tolist()
    else:
        cb_meff_E = model.cb_meff.tolist() #lists are faster than numpy arrays for loops
    fis = fis.tolist() #lists are faster than numpy arrays for loops
    n_max = model.n_max
    dx = model.dx
    #
    N = 0.0 # Normalization integral
    psi = []
    psi = [0.0]*3
    # boundary conditions
    psi[0] = 0.0                 
    psi[1] = 1.0
    b = [0.0]*n_max
    b[0] = psi[0]
    b[1] = psi[1]
    N += (psi[0])**2
    N += (psi[1])**2
    for j in range(1,n_max-1,1):
        # Last potential not used
        c1=2.0/(cb_meff_E[j]+cb_meff_E[j-1])
        c2=2.0/(cb_meff_E[j]+cb_meff_E[j+1])
        psi[2] = ((2*(dx/hbar)**2*(fis[j]-E)+c2+c1)*psi[1]-c1*psi[0])/c2
        b[j+1]=psi[2]
        N += (psi[2])**2
        psi[0]=psi[1]
        psi[1]=psi[2]
    b2 = np.array(b)
    b2/= N**0.5
    return b2 # units of dx**0.5

    
# FUNCTIONS for FERMI-DIRAC STATISTICS---SIMPLE---------------------------------   
def fd2(Ei,Ef,T):
    """integral of Fermi Dirac Equation for energy independent density of states.
    Ei [meV], Ef [meV], T [K]"""
    return kb*T*log(exp(meV2J*(Ef-Ei)/(kb*T))+1)

def calc_meff_state(wfe,cb_meff):
    """find subband effective masses"""
    tmp = 1.0/np.sum(wfe**2/cb_meff,axis=1)
    meff_state = tmp.tolist()
    return meff_state #kg

def calc_meff_state2(wfe,E_state,fi,model):
    """find subband effective masses including non-parabolicity
    (but stilling using a fixed effective mass for each subband dispersion)"""
    cb_meff_states = np.vstack([model.cb_meff_E(E*meV2J,fi) for E in E_state])
    tmp = 1.0/np.sum(wfe**2/cb_meff_states,axis=1)
    meff_state = tmp.tolist()
    return meff_state #kg
    
def fermilevel_0K(Ntotal2d,E_state,meff_state):
    Et,Ef=0.0,0.0
    for i,(Ei,csb_meff) in enumerate(zip(E_state,meff_state)):
        Et+=Ei
        Efnew=(Ntotal2d*hbar**2*pi/csb_meff*J2meV + Et)/(i+1)
        if Efnew>Ei:
            Ef=Efnew
        else:
            break #we have found Ef and so we should break out of the loop
    else: #exception clause for 'for' loop.
        logger.warning("Have processed all energy levels present and so can't be sure that Ef is below next higher energy level.")
    N_state=[0.0]*len(E_state)
    for i,(Ei,csb_meff) in enumerate(zip(E_state,meff_state)):
        Ni=(Ef - Ei)*csb_meff/(hbar**2*pi)*meV2J    # populations of levels
        Ni*=(Ni>0.0)
        N_state[i]=Ni
    return Ef,N_state #Fermi levels at 0K (meV), number of electrons in each subband at 0K
    
def fermilevel(Ntotal2d,T,E_state,meff_state):
    """Finds the Fermi level (meV)"""
    #parameters
    FD_d_E = config.FD_d_E #1e-9 Initial and minimum Energy step (meV) for derivative calculation for Newton-Raphson method to find E_F
    FD_convergence_test = config.FD_convergence_test #1e-6
    
    def func(Ef,E_state,meff_state,Ntotal2d,T):
        #return Ntotal2d - sum( [csb_meff*fd2(Ei,Ef,T) for Ei,csb_meff in zip(E_state,meff_state)] )/(hbar**2*pi)
        diff = Ntotal2d
        for Ei,csb_meff in zip(E_state,meff_state):
            diff -= csb_meff*fd2(Ei,Ef,T)/(hbar**2*pi)
        return diff
    Ef_0K,N_states_0K = fermilevel_0K(Ntotal2d,E_state,meff_state)
    #Ef=fsolve(func,Ef_0K,args=(E_state,meff_state,Ntotal2d,T))[0]
    #return float(Ef)
    #implement Newton-Raphson method
    Ef = Ef_0K
    d_E = FD_d_E #Energy step (meV)
    while True:
        y = func(Ef,E_state,meff_state,Ntotal2d,T)
        dy = (func(Ef+d_E,E_state,meff_state,Ntotal2d,T)- func(Ef-d_E,E_state,meff_state,Ntotal2d,T))/(2.0*d_E)
        if dy == 0.0: #increases interval size for derivative calculation in case of numerical error
            d_E*=2.0
            continue #goes back to start of loop, therefore d_E will increase until a non-zero derivative is found
        Ef -= y/dy
        if abs(y/dy) < FD_convergence_test:
            break
        #reduces the interval by a couple of notches ready for the next iteration
        for i in range(2):
            if d_E>FD_d_E: 
                d_E*=0.5
    return Ef #(meV)

def calc_N_state(Ef,T,Ns,E_state,meff_state):
    # Find the subband populations, taking advantage of step like d.o.s. and analytic integral of FD
    N_state=[fd2(Ei,Ef,T)*csb_meff/(hbar**2*pi) for Ei,csb_meff in zip(E_state,meff_state)]
    return N_state # number of carriers in each subband
    
# FUNCTIONS for FERMI-DIRAC STATISTICS---NON-PARABOLIC--------------------------

def calc_meff_state3(wfe,cb_meff):
    """find subband effective masses"""
    meff_states = 1.0/np.sum(wfe**2/cb_meff,axis=1)
    return meff_states #kg

def calc_dispersions(Emin,Emax,dE,wfe,E_state,fi,model):
    """Calculate dispersion curves and their effective masses for each subband."""
    output = [] #(Ea,cb_meff,k)
    for Ei,wfi in zip(E_state,wfe):
        Ea = np.arange(Ei,Emax,dE) #meV
        cb_meff_2darray = model.cb_meff_E(Ea[:,np.newaxis]*meV2J,fi)
        cb_meff_a = calc_meff_state3(wfi,cb_meff_2darray) # effective mass of dispersion wrt Ea
        ka = np.sqrt(2.0/hbar**2*meV2J*(Ea-Ei)*cb_meff_a) # k-space array wrt Ea
        output.append((Ea,cb_meff_a,ka))
    return output

def calc_N_state_np(Ef,T,level_dispersions):
    """numerical integral of Fermi Dirac Equation for 2d density of states with 
    energy dependent effective mass
    level_dispersions - list of tuples describing each level's dispersion. Each tuple
    is (energy array, effective mass array, k-vector array). This is non-coincidently
    the output of calc_dispersions.
    function.
    Ef - Fermi energy (meV)
    T - Temperature (K)."""
    N_state = []
    for Ea,cb_meff_a,ka in level_dispersions: #treat level separately
        Ea2 = Ea*meV2J
        tmp = cb_meff_a/(pi*hbar**2)/(np.exp((Ea2-meV2J*Ef)/(kb*T))+1.0)
        N_state.append(np.trapz(tmp,x=Ea2))
    return N_state

def fermilevel_np(Ntotal2d,T,wfe,E_state,fi,model):
    """Finds the Fermi level (meV) for non-parabolic subbands"""
    #parameters
    FD_d_E = config.FD_d_E #1e-9 Initial and minimum Energy step (meV) for derivative calculation for Newton-Raphson method to find E_F
    FD_convergence_test = config.FD_convergence_test #1e-6
    #level dispersions
    level_dispersions = calc_dispersions(Emin=0.0,Emax=E_state[-1]+kb*T*J2meV,dE=config.np_d_E,wfe=wfe,E_state=E_state,fi=fi,model=model)
    #error function
    def func(Ef,Ntotal2d,T,level_dispersions):
        return Ntotal2d - sum(calc_N_state_np(Ef,T,level_dispersions))
    #starting estimates
    Ef_0K,N_states_0K = fermilevel_0K(Ntotal2d,E_state,calc_meff_state2(wfe,E_state,fi,model))
    #Ef=fsolve(func,Ef_0K,args=(E_state,meff_state,Ntotal2d,T))[0]
    #return float(Ef)
    #implement Newton-Raphson method
    Ef = Ef_0K
    d_E = FD_d_E #Energy step (meV)
    while True:
        y = func(Ef,Ntotal2d,T,level_dispersions)
        dy = (func(Ef+d_E,Ntotal2d,T,level_dispersions)- func(Ef-d_E,Ntotal2d,T,level_dispersions))/(2.0*d_E)
        if dy == 0.0: #increases interval size for derivative calculation in case of numerical error
            d_E*=2.0
            continue #goes back to start of loop, therefore d_E will increase until a non-zero derivative is found
        Ef -= y/dy
        if abs(y/dy) < FD_convergence_test:
            break
        #reduces the interval by a couple of notches ready for the next iteration
        for i in range(2):
            if d_E>FD_d_E: 
                d_E*=0.5
    return Ef #(meV)
    
# FUNCTIONS for SELF-CONSISTENT POISSON----------------------------------------

def calc_sigma(wfe,N_state,model):
    """This function calculates `net' areal charge density
    n-type dopants lead to -ve charge representing electrons, and additionally 
    +ve ionised donors."""
    # note: model.dop is still a volume density, the delta_x converts it to an areal density
    sigma= model.dop*model.dx # The charges due to the dopant ions
    for j in range(0,model.subnumber_e,1): # The charges due to the electrons in the subbands
        sigma-= N_state[j]*(wfe[j])**2
    return sigma #charge per m**2 per dz (units of electronic charge)
    
##
def calc_field(sigma,eps):
    """calculate electric field as a function of z-
    sigma is a number density per unit area
    eps is dielectric constant"""
    # i index over z co-ordinates
    # j index over z' co-ordinates
    # Note: 
    F0 = -np.sum(q*sigma)/(2.0) #CMP'deki i ve j yer değişebilir - de + olabilir
    # is the above necessary since the total field due to the structure should be zero.
    # Do running integral
    tmp = np.hstack(([0.0],sigma[:-1])) + sigma
    tmp*= q/2.0 # Note: sigma is a number density per unit area, needs to be converted to Couloumb per unit area
    tmp[0] = F0 
    F = np.cumsum(tmp)/eps
    return F #electric field

def calc_field_convolve(sigma,eps):
    tmp = np.ones(len(sigma)-1)
    signstep = np.hstack((-tmp,[0.0],tmp)) # step function
    F = np.convolve(signstep,sigma,mode='valid')
    F*= q/(2.0*eps)
    return F

def calc_field_old(sigma,eps):
    """calculate F electric field as a function of z-
    sigma is a number density per unit area,
    eps is dielectric constant"""
    # i index over z co-ordinates
    # j index over z' co-ordinates
    n_max = len(sigma)
    # For wave function initialise F
    F = np.zeros(n_max)
    for i in range(0,n_max,1):
        for j in range(0,n_max,1):
           # Note sigma is a number density per unit area, needs to be converted to Couloumb per unit area
           F[i] = F[i] + q*sigma[j]*cmp(i,j)/(2*eps[i]) #CMP'deki i ve j yer değişebilir - de + olabilir
    return F

def calc_potn(F,dx):
    """This function calculates the potential (energy actually)"""
    # V electric field as a function of z-
    # i	index over z co-ordinates

    #Calculate the potential, defining the first point as zero
    tmp = q*F*dx
    V = np.cumsum(tmp) #+q -> electron -q->hole? 
    return V


# FUNCTIONS FOR EXCHANGE INTERACTION-------------------------------------------

def calc_Vxc(sigma,dop,eps,cb_meff,dx):
    """An effective field describing the exchange-interactions between the electrons
    derived from Kohn-Sham density functional theory. This formula is given in many
    papers, for example see Gunnarsson and Lundquist (1976), Ando, Taniyama, Ohtani 
    et al. (2003), or Ch.1 in the book 'Intersubband transitions in quantum wells' (edited
    by Liu and Capasso) by M. Helm.
    eps = dielectric constant array
    cb_meff = effective mass array
    sigma = charge carriers per m**2, however this includes the donor atoms and we are only
            interested in the electron density."""
    a_B = 4*pi*hbar**2/q**2 # Bohr radius.
    nz= -(sigma - dop*dx) # electron density per m**2
    nz_3 = nz**(1/3.) #cube root of charge density.
    #a_B_eff = eps/cb_meff*a_B #effective Bohr radius
    #r_s occasionally suffers from division by zero errors due to nz=0.
    #We will fix these by setting nz_3 = 1.0 for these points (a tiny charge in per m**2).
    nz_3 = nz_3.clip(1.0,max(nz_3))
    
    r_s = 1.0/((4*pi/3.0)**(1/3.)*nz_3*eps/cb_meff*a_B) #average distance between charges in units of effective Bohr radis.
    #A = q**4/(32*pi**2*hbar**2)*(9*pi/4.0)**(1/3.)*2/pi*(4*pi/3.0)**(1/3.)*4*pi*hbar**2/q**2 #constant factor for expression.
    A = q**2/(4*pi)*(3/pi)**(1/3.) #simplified constant factor for expression.
    #
    Vxc = -A*nz_3/eps * ( 1.0 + 0.0545 * r_s * np.log( 1.0 + 11.4/r_s) )
    return Vxc

# -----------------------------------------------------------------------------

def Poisson_Schrodinger(model):
    """Performs a self-consistent Poisson-Schrodinger calculation of a 1d quantum well structure.
    Model is an instance of the Structure class or an object with the following attributes:
      fi - float - Bandstructure potential (J) (array, len n_max)
      cb_meff - float - conduction band effective mass (kg)(array, len n_max)
      eps - float - dielectric constant (including eps0) (array, len n_max)
      dop - float - doping distribution (m**-3) ( array, len n_max)
      Fapp - float - Applied field (Vm**-1)
      T - float - Temperature (K)
      comp_scheme - int - simulation scheme
      subnumber_e - int - number of subbands for look for in the conduction band
      dx - float - grid spacing (m)
      n_max - int - number of points.
    """   
    fi = model.fi
    cb_meff = model.cb_meff # effective mass at band edge (effective mass with non-parabolicity is model.cb_meff_E(E,fi))
    eps = model.eps
    dop = model.dop
    Fapp = model.Fapp
    T = model.T
    comp_scheme = model.comp_scheme
    subnumber_e = model.subnumber_e
    dx = model.dx
    n_max = model.n_max
    
    #parameters
    E_start = config.E_start #0.0 #Energy to start shooting method from (if E_start = 0.0 uses minimum of energy of bandstructure)
    # Poisson Loop
    damping = config.damping #0.5 #averaging factor between iterations to smooth convergence.
    max_iterations= config.max_iterations #80 #maximum number of iterations.
    convergence_test= config.convergence_test #1e-6 #convergence is reached when the ground state energy (meV) is stable to within this number between iterations.
    
    # Check
    if comp_scheme ==6:
        scheme6warning = """The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere."""
        logger.warning(scheme6warning)
    
    # Preparing empty subband energy lists.
    E_state = [0.0]*subnumber_e     # Energies of subbands/levels (meV)
    N_state = [0.0]*subnumber_e     # Number of carriers in subbands  
    
    # Creating and Filling material arrays
    xaxis = np.arange(0,n_max)*dx       #metres
    fitot = np.zeros(n_max)             #Energy potential = Bandstructure + Coulombic potential
    #sigma = np.zeros(n_max)            #charge distribution (donors + free charges)
    #F = np.zeros(n_max)                #Electric Field
    #Vapp = np.zeros(n_max)             #Applied Electric Potential
    V = np.zeros(n_max)                 #Electric Potential
    
    # Subband wavefunction for electron list. 2-dimensional: [i][j] i:stateno, j:wavefunc
    wfe = np.zeros((subnumber_e,n_max))
    
    # Setup the doping
    Ntotal = sum(dop) # calculating total doping density m-3
    Ntotal2d = Ntotal*dx
    if not(config.messagesoff):
        #print "Ntotal ",Ntotal,"m**-3"
        logger.info("Ntotal2d %g m**-2", Ntotal2d)
    
    #Applied Field
    x0 = dx*n_max/2.0
    Vapp = q*Fapp*(xaxis-x0)

    # STARTING SELF CONSISTENT LOOP
    time2 = time.time() # timing audit
    iteration = 1   #iteration counter
    previousE0= 0   #(meV) energy of zeroth state for previous iteration(for testing convergence)
    fitot = fi + Vapp #For initial iteration sum bandstructure and applied field
    
    fi_min= min(fitot) #minimum potential energy of structure (for limiting the energy range when searching for states)
    if abs(E_start)>1e-3*meV2J: #energyx is the minimum energy (meV) when starting the search for bound states.
        energyx = E_start
    else:
        energyx = fi_min
        
    while True:
        if not(config.messagesoff) :
            logger.info("Iteration: %d", iteration)
        if iteration> 1:
            energyx=min(fi_min,min(fitot),)
        
        E_state=calc_E_state(subnumber_e,fitot,model,energyx)
        
        # Envelope Function Wave Functions
        #print 'wf'
        for j in range(0,subnumber_e,1):
            if not(config.messagesoff):
                logger.info("Working for subband no: %d",(j+1))
            wfe[j] = wf(E_state[j]*meV2J,fitot,model) #wavefunction units dx**0.5

        # Calculate the effective mass of each subband
        #print 'calc_meff_state'
        if model.comp_scheme in (1,3,6): #include non-parabolicity in calculation
            meff_state = calc_meff_state2(wfe,E_state,fitot,model)
        else:
            meff_state = calc_meff_state(wfe,cb_meff)
        
        ## Self-consistent Poisson
        
        # Calculate the Fermi energy and subband populations at 0K
        #E_F_0K,N_state_0K=fermilevel_0K(Ntotal2d,E_state,meff_state)
        # Calculate the Fermi energy at the temperature T (K)
        if model.comp_scheme in (1,3,6) and model.fermi_np_scheme == True: #include non-parabolicity in calculation
            #print 'fermilevel'
            E_F = fermilevel_np(Ntotal2d,T,wfe,E_state,fi,model)
            # Calculate the subband populations at the temperature T (K)
            #print 'calc_N_state'
            Emin_np = 0.0
            Emax_np = E_state[-1]+kb*T*J2meV
            dE_np = config.np_d_E
            level_dispersions = calc_dispersions(Emin_np,Emax_np,dE_np,wfe,E_state,fi,model)
            N_state = calc_N_state_np(E_F,T,level_dispersions)
        else:
            level_dispersions = None
            #print 'fermilevel'
            E_F = fermilevel(Ntotal2d,T,E_state,meff_state)
            # Calculate the subband populations at the temperature T (K)
            #print 'calc_N_state'
            N_state = calc_N_state(E_F,T,Ntotal2d,E_state,meff_state)
        # Calculate `net' areal charge density
        #print 'calc_sigma'
        sigma=calc_sigma(wfe,N_state,model) #one more instead of subnumber_e
        # Calculate electric field (Poisson/Hartree Effects)
        if comp_scheme != 4: #in (0,1,2,3,5,6):
            # Calculate electric field
            F=calc_field(sigma,eps)
            # Calculate potential due to charge distribution
            Vnew=calc_potn(F,dx)
        else:
            F=np.zeros(n_max)
            Vnew=0
        # Exchange interaction    
        if comp_scheme in (4,5,6):
            # Exchange Potential
            Vnew += calc_Vxc(sigma,dop,eps,cb_meff,dx)
            
        #
        #status
        if not(config.messagesoff):
            for i,level in enumerate(E_state):
                logger.info("E[%d]= %f meV",i,level)
            for i,meff in enumerate(meff_state):
                logger.info("meff[%d]= %f",i,meff/m_e)
            for i,Ni in enumerate(N_state):
                logger.info("N[%d]= %g m**-2",i,Ni)
            #print 'Efermi (at 0K) = ',E_F_0K,' meV'
            #for i,Ni in enumerate(N_state_0K):
            #    print 'N[',i,']= ',Ni
            logger.info('Efermi (at %gK) = %g meV',T, E_F)
            logger.info("total donor charge = %g m**-2",sum(dop)*dx)
            logger.info("total level charge = %g m**-2",sum(N_state))
            logger.info("total system charge = %g m**-2",sum(sigma))
        #
        if comp_scheme in (0,1): 
            #if we are not self-consistently including Poisson Effects then only do one loop
            break
        
        # Combine band edge potential with potential due to charge distribution
        # To increase convergence, we calculate a moving average of electric potential 
        #with previous iterations. By dampening the corrective term, we avoid oscillations.
        V+= damping*(Vnew - V)
        fitot = fi + V + Vapp
        
        if abs(E_state[0]-previousE0) < convergence_test: #Convergence test
            break
        elif iteration >= max_iterations: #Iteration limit
            logger.warning("Have reached maximum number of iterations")
            break
        else:
            iteration += 1
            previousE0 = E_state[0]
    
    # END OF SELF-CONSISTENT LOOP
    time3 = time.time() # timing audit
    logger.info("calculation time  %g s",(time3 - time2))
    
    class Results(): pass
    results = Results()
    
    results.xaxis = xaxis
    results.wfe = wfe
    results.fitot = fitot
    results.sigma = sigma
    results.F = F
    results.V = V
    results.E_state = E_state
    results.N_state = N_state
    results.meff_state = meff_state
    results.Fapp = Fapp
    results.T = T
    results.E_F = E_F
    results.dx = dx
    results.subnumber_e = subnumber_e
    results.level_dispersions = level_dispersions
    
    return results

def save_and_plot(result,model):
    xaxis = result.xaxis
    
    output_directory = config.output_directory
    
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
        
    def saveoutput(fname,datatuple,header=None):
        fname2 = os.path.join(output_directory,fname)
        fobj = file(fname2,'wb')
        if header: fobj.write(header+'\n')
        np.savetxt(fobj,np.column_stack(datatuple),fmt='%.6e', delimiter=' ')
        fobj.close()
        
    def saveoutput2(fname2,datatuple,header=None,fmt='%.6g',delimiter=', '):
        fname2 = os.path.join(output_directory,fname2)
        fobj = file(fname2,'wb')
        if header: fobj.write(header+'\n')
        np.savetxt(fobj,np.column_stack(datatuple),fmt=fmt, delimiter=delimiter)
        fobj.close()
    
    if config.parameters:
        saveoutput2("parameters.dat",header=('T (K), Fapp (V/m), E_F (meV)'),
                    datatuple=(result.T,result.Fapp,result.E_F))
    if config.sigma_out:
        saveoutput("sigma.dat",(xaxis,result.sigma))
    if config.electricfield_out:
        saveoutput("efield.dat",(xaxis,result.F))
    if config.potential_out:
        saveoutput("potn.dat",(xaxis,result.fitot))
    if config.states_out:
        rel_meff_state = [meff/m_e for meff in result.meff_state] #going to report relative effective mass.
        columns = range(model.subnumber_e), result.E_state, result.N_state, rel_meff_state
        #header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
        header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
        saveoutput("states.dat",columns, header = header )
    if config.probability_out:
        saveoutput("wavefunctions.dat",(xaxis,result.wfe.transpose()) )
            
    # Resultviewer
        
    if config.resultviewer:
        fig1 = pl.figure(figsize=(10,8))
        pl.suptitle('Aestimo Results')
        pl.subplots_adjust(hspace=0.4,wspace=0.4)
                            
        #Plotting Sigma
        #figure(0)
        pl.subplot(2,2,1)
        pl.plot(xaxis, result.sigma)
        pl.xlabel('Position (m)')
        pl.ylabel('Sigma (e/m^2)')
        pl.title('Sigma')
        pl.grid(True)
    
        #Plotting Efield
        #figure(1)
        pl.subplot(2,2,2)
        pl.plot(xaxis, result.F)
        pl.xlabel('Position (m)')
        pl.ylabel('Electric Field strength (V/m)')
        pl.title('Electric Field')
        pl.grid(True)
    
        #Plotting Potential
        #figure(2)
        pl.subplot(2,2,3)
        pl.plot(xaxis, result.fitot)
        pl.xlabel('Position (m)')
        pl.ylabel('E_c (J)')
        pl.title('Potential')
        pl.grid(True)
    
        #Plotting State(s)
        #figure(3)
        pl.subplot(2,2,4)
        for j,state in enumerate(result.wfe):
            pl.plot(xaxis, state, label='state %d' %j)
        pl.xlabel('Position (m)')
        pl.ylabel('Psi')
        pl.title('First state')
        pl.grid(True)
        
        #QW representation
        #figure(5)
        fig2 = pl.figure(figsize=(10,8))
        pl.suptitle('Aestimo Results')
        pl.subplot(1,1,1)
        pl.plot(xaxis, result.fitot*J2meV,'k')
        for level,state in zip(result.E_state,result.wfe): 
            pl.axhline(level,0.1,0.9,color='g',ls='--')
            pl.plot(xaxis, state*config.wavefunction_scalefactor+level,'b')
            #pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
        pl.axhline(result.E_F,0.1,0.9,color='r',ls='--')
        pl.xlabel('Position (m)')
        pl.ylabel('Energy (meV)')
        pl.grid(True)
        
        #dispersion plot
        fig3 = pl.figure(figsize=(10,8))
        pl.suptitle('Subband Dispersions')
        ax = pl.subplot(1,1,1)
        result.level_dispersions
        cb_meff0 = result.meff_state[0] #kg
        ka = np.linspace(0.0,np.sqrt(2.0*cb_meff0*result.fitot.ptp())/hbar,50) #m**-1
        kax = ka*1e-9
        for Ei,meffi in zip(result.E_state,result.meff_state):
            p1, = pl.plot(kax,Ei+J2meV*hbar**2*ka**2/(2*cb_meff0),'k')
            p2, = pl.plot(kax,Ei+J2meV*hbar**2*ka**2/(2*meffi),'b')
        if result.level_dispersions:
            for Ea,cb_meff_a,ka in result.level_dispersions:
                p3, = pl.plot(ka*1e-9,Ea,'g')
            ax.legend([p1,p2,p3],['parabolic dispersions','parabolic dispersion (subband meff)','non-parabolic dispersions'])
        else:
            ax.legend([p1,p2],['parabolic dispersions','parabolic dispersion (subband meff)'])
        pl.axhline(result.E_F,0.0,1.0,color='r',ls='--')
        pl.xlabel('k-space (nm**-1)')
        pl.ylabel('Energy (meV)')
        pl.grid(True)
        
        pl.show()
    return [fig1,fig2,fig3]

def QWplot(result,figno=None):
    """QW representation"""
    xaxis = result.xaxis
    fig = pl.figure(figno,figsize=(10,8))
    pl.suptitle('Aestimo Results')
    ax = pl.subplot(1,1,1)
    ax.plot(xaxis, result.fitot*J2meV,'k')
    for level,state in zip(result.E_state,result.wfe): 
        ax.axhline(level,0.1,0.9,color='g',ls='--')
        ax.plot(xaxis, state*config.wavefunction_scalefactor+level,'b')
        #pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
    ax.axhline(result.E_F,0.1,0.9,color='r',ls='--')
    pl.xlabel('Position (m)')
    pl.ylabel('Energy (meV)')
    ax.grid(True)
    pl.show()
    return fig
        
def dispersionplot(result,figno=None):
    """subband dispersion plot"""
    fig = pl.figure(figno,figsize=(10,8))
    pl.suptitle('Subband Dispersions')
    ax = pl.subplot(1,1,1)
    result.level_dispersions
    cb_meff0 = result.meff_state[0] #kg
    kmax = np.sqrt(2.0*cb_meff0*result.fitot.ptp()*meV2J)/hbar
    ka = np.linspace(0.0,kmax,50) #m**-1
    kax = ka*1e-9
    for Ei,meffi in zip(result.E_state,result.meff_state):
        p1, = ax.plot(kax,Ei+J2meV*hbar**2*ka**2/(2*cb_meff0),'k')
        p2, = ax.plot(kax,Ei+J2meV*hbar**2*ka**2/(2*meffi),'b')
    if result.level_dispersions:
        for Ea,cb_meff_a,ka in result.level_dispersions:
            p3, = ax.plot(ka*1e-9,Ea,'g')
        ax.legend([p1,p2,p3],['parabolic dispersions','parabolic dispersion (subband meff)','non-parabolic dispersions'])
    else:
        ax.legend([p1,p2],['parabolic dispersions','parabolic dispersion (subband meff)'])
    ax.axhline(result.E_F,0.0,1.0,color='r',ls='--')
    pl.xlabel('k-space (nm**-1)')
    pl.ylabel('Energy (meV)')
    ax.grid(True)
    pl.show()
    return fig


def load_results():
    """Loads the data stored in the output folder"""
    class Results(): pass
    results = Results()
    
    output_directory = config.output_directory
            
    def loadoutput(fname,header=False,unpack=True):
        fname2 = os.path.join(output_directory,fname)
        fobj = file(fname2,'rb')
        if header: header = fobj.readline()
        else: header = ''
        data = np.loadtxt(fobj,delimiter=' ',unpack=unpack)
        fobj.close()
        return data,header
    
    if config.parameters:
        results.T,results.Fapp,results.E_F = np.loadtxt(
                      open(os.path.join(output_directory,"parameters.dat"),'rb'),
                      unpack=True,delimiter=',',skiprows=1)
    if config.sigma_out:
        (results.xaxis,results.sigma),hdr = loadoutput("sigma.dat")
    if config.electricfield_out:
        (results.xaxis,results.F),hdr = loadoutput("efield.dat")
    if config.potential_out:
        (results.xaxis,results.fitot),hdr = loadoutput("potn.dat")
    if config.states_out:
        (states,results.E_state,results.N_state,rel_meff_state),hdr = loadoutput("states.dat", header=True)
        results.subnumber_e = max(states)
        results.meff_state = rel_meff_state*m_e
    if config.probability_out:
        _wfe,hdr = loadoutput("wavefunctions.dat",unpack=False)
        results.xaxis = _wfe[:,0]
        results.wfe = _wfe[:,1:].transpose()
    
    #missing variables
    #results.V
    results.dx = np.mean(results.xaxis[1:]-results.xaxis[:-1])
    #results.level_dispersions = level_dispersions
    
    return results

def run_aestimo(input_obj):
    """A utility function that performs the standard simulation run
    for 'normal' input files. Input_obj can be a dict, class, named tuple or 
    module with the attributes needed to create the StructureFrom class, see 
    the class implementation or some of the sample-*.py files for details."""
    logger.info("Aestimo is starting...")
        
    # Initialise structure class
    model = StructureFrom(input_obj,database)
         
    # Perform the calculation
    result = Poisson_Schrodinger(model)
    
    time4 = time.time() #timing audit
    logger.info("total running time (inc. loading libraries) %g s",(time4 - time0))
    logger.info("total running time (exc. loading libraries) %g s",(time4 - time1))

    
    # Write the simulation results in files
    save_and_plot(result,model)
    
    logger.info("""Simulation is finished. All files are closed. Please control the related files.
-----------------------------------------------------------------""")
    
    return input_obj, model, result


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-i","--inputfile",action="store", dest="inputfile", 
                  default=config.inputfilename,
                  help="chose input file to override default in config.py")
    (options, args) = parser.parse_args()
    
    # Import from config file
    inputfile = __import__(options.inputfile)
    logger.info("inputfile is %s",options.inputfile)
    
    run_aestimo(inputfile)
    
