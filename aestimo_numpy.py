#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Version v.0.8
 Copyright (C) 2013 Sefer Bora Lisesivdin and Aestimo group

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

  Description: This is the aestimo calculator conduction band calculations (Numpy version).
  
"""
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
logger = logging.getLogger('aestimo_numpy')
hdlr = logging.FileHandler(config.logfile)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
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
print "Aestimo_numpy is starting..."
logger.info("Aestimo_numpy is starting...")

# Input Class
# -------------------------------------

def round2int(x):
    """int is sensitive to floating point numerical errors near whole numbers,
    this moves the discontinuity to the half interval. It is also equivalent
    to the normal rules for rounding positive numbers."""
    # int(x + (x>0) -0.5) # round2int for positive and negative numbers
    return int(x+0.5)

class Structure():
    def __init__(self,database,**kwargs):
        """This class holds details on the structure to be simulated.
        database is the module containing the material properties. Then
        this class should have the following attributes set
        Fapp - applied field (Vm**-1)
        T - Temperature (K)
        subnumber_e - number of subbands to look for.
        comp_scheme - computing scheme
        dx - grid step size (m)
        n_max - number of grid points
        
        cb_meff #conduction band effective mass (kg) (array, len n_max)
        cb_meff_alpha #non-parabolicity constant.
        fi #Bandstructure potential (J) (array, len n_max)
        eps #dielectric constant (including eps0) (array, len n_max)
        dop #doping distribution (m**-3) (array, len n_max)
        
        These last 4 can be created by using the method 
        create_structure_arrays(material_list)
        """
        # setting any parameters provided with initialisation
        for key,value in kwargs.items():
            setattr(self,key,value)
                
        # Loading materials database
        self.material_property = database.materialproperty
        totalmaterial = alen(self.material_property)
        
        self.alloy_property = database.alloyproperty
        totalalloy = alen(self.alloy_property)
        
        print "Total material number in database: ",totalmaterial + totalalloy
        logger.info("Total material number in database: %d" %(totalmaterial + totalalloy))
        
    def create_structure_arrays(self):
        """ initialise arrays/lists for structure"""
        # Calculate the required number of grid points
        self.x_max = sum([layer[0] for layer in self.material])*1e-9 #total thickness (m)
        n_max = round2int(self.x_max/self.dx)
        # Check on n_max
        maxgridpoints = self.maxgridpoints
        if n_max > maxgridpoints:
            print " Grid number is exceeding the max number of ", maxgridpoints
            logger.error("Grid number is exceeding the max number of %d" %max_val)
            exit()
        #
        self.n_max = n_max
        dx =self.dx
        material_property = self.material_property
        alloy_property = self.alloy_property
        
        cb_meff = np.zeros(n_max)	#conduction band effective mass
        cb_meff_alpha = np.zeros(n_max) #non-parabolicity constant.
        fi = np.zeros(n_max)		#Bandstructure potential
        eps =np.zeros(n_max)		#dielectric constant
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
                cb_meff[startindex:finishindex] = matprops['m_e']*m_e
                cb_meff_alpha[startindex:finishindex] = matprops['m_e_alpha']
                fi[startindex:finishindex] = matprops['Band_offset']*matprops['Eg']*q #Joule
                eps[startindex:finishindex] = matprops['epsilonStatic']*eps0
                
            elif matType in alloy_property:
                alloyprops = alloy_property[matType] 
                x = layer[2] #alloy ratio
                cb_meff_alloy = x*material_property[alloyprops['Material1']]['m_e'] + (1-x)* material_property[alloyprops['Material2']]['m_e']
                cb_meff[startindex:finishindex] = cb_meff_alloy*m_e
                cb_meff_alpha[startindex:finishindex] = alloyprops['m_e_alpha']*(material_property[alloyprops['Material2']]['m_e']/cb_meff_alloy) #non-parabolicity constant for alloy. THIS CALCULATION IS MOSTLY WRONG. MUST BE CONTROLLED. SBL
                fi[startindex:finishindex] = alloyprops['Band_offset']*(x*material_property[alloyprops['Material1']]['Eg'] + (1-x)* material_property[alloyprops['Material2']]['Eg']-alloyprops['Bowing_param']*x*(1-x))*q # for electron. Joule
                eps[startindex:finishindex] = (x*material_property[alloyprops['Material1']]['epsilonStatic'] + (1-x)* material_property[alloyprops['Material2']]['epsilonStatic'] )*eps0
                
            #doping
            if layer[4] == 'n':  
                chargedensity = layer[3]*1e6 #charge density in m**-3 (conversion from cm**-3)
            elif layer[4] == 'p': 
                chargedensity = -layer[3]*1e6 #charge density in m**-3 (conversion from cm**-3)
            else:
                chargedensity = 0.0
            
            dop[startindex:finishindex] = chargedensity
        
        self.fi = fi
        self.cb_meff = cb_meff
        self.cb_meff_alpha = cb_meff_alpha
        self.eps = eps
        self.dop = dop
        #return fi,cb_meff,eps,dop

class AttrDict(dict):
    """turns a dictionary into an object with attribute style lookups"""
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class StructureFrom(Structure):
    def __init__(self,inputfile,database):
        if type(inputfile)==dict:
            inputfile=AttrDict(inputfile)            
        # Parameters for simulation
        self.Fapp = inputfile.Fapplied
        self.T = inputfile.T
        self.subnumber_e = inputfile.subnumber_e
        self.comp_scheme = inputfile.computation_scheme
        self.dx = inputfile.gridfactor*1e-9 #grid in m
        self.maxgridpoints = inputfile.maxgridpoints
        
        # Loading material list
        self.material = inputfile.material
        totallayer = alen(self.material)
        print "Total layer number: ",totallayer
        logger.info("Total layer number: %s" %totallayer)
        
        # Calculate the required number of grid points
        self.x_max = sum([layer[0] for layer in self.material])*1e-9 #total thickness (m)
        self.n_max = int(self.x_max/self.dx)
        
        # Check on n_max
        max_val = inputfile.maxgridpoints
        if self.n_max > max_val:
            print " Grid number is exceeding the max number of ", max_val
            exit()
                
        # Loading materials database #
        self.material_property = database.materialproperty
        totalmaterial = alen(self.material_property)
        
        self.alloy_property = database.alloyproperty
        totalalloy = alen(self.alloy_property)
        
        print "Total number of materials in database: ",totalmaterial+totalalloy
        logger.info("Total number of materials in database: %d" %(totalmaterial+totalalloy))
        
        # Initialise arrays
        
        #cb_meff #conduction band effective mass (array, len n_max)
        #fi #Bandstructure potential (array, len n_max)
        #eps #dielectric constant (array, len n_max)
        #dop #doping distribution (array, len n_max)
        self.create_structure_arrays()



# Shooting method parameters for Schrödinger Equation solution
delta_E = 1.0*meV2J #Energy step (Joules) for initial search. Initial delta_E is 1 meV. #This can be included in config as a setting?
d_E = 1e-5*meV2J #Energy step (Joules) for Newton-Raphson method when improving the precision of the energy of a found level.
E_start = 0.0    #Energy to start shooting method from #This can be included in config as a setting?
damping = 0.5    #averaging factor between iterations to smooth convergence.
max_iterations=80 #maximum number of iterations.
convergence_test=1e-6 #convergence is reached when the ground state energy (meV) is stable to within this number between iterations.


# DO NOT EDIT UNDER HERE FOR PARAMETERS
# --------------------------------------

#Vegard's law for alloys
def vegard(first,second,mole):
    return first*mole+second*(1-mole)

# This function returns the value of the wavefunction (psi)
# at +infinity for a given value of the energy.  The solution
# to the energy occurs for psi(+infinity)=0.

# FUNCTIONS for SHOOTING ------------------
if config.use_cython:
    from psi_at_inf_cython import psi_at_inf_numpy

    def psi_at_inf1(E,fis,cb_meff,cb_meff_alpha,n_max,dx):
        return psi_at_inf_numpy(E,fis,cb_meff,n_max,dx)

    def psi_at_inf2(E,fis,cb_meff,cb_meff_alpha,n_max,dx):
        """shooting method with non-parabolicity"""
        cb_meff_E = cb_meff*(1.0 + cb_meff_alpha*(E - fis))
        return psi_at_inf_numpy(E,fis,cb_meff_E,n_max,dx)
        
else:
    
    def psi_at_inf1(E,fis,cb_meff,cb_meff_alpha,n_max,dx): #inclusion of cb_meff_alpha only so both versions have the same calling arguments
        """Shooting method for heterostructure as given in Harrison's book"""
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
        
    def psi_at_inf2(E,fis,cb_meff,cb_meff_alpha,n_max,dx):
        """shooting method with non-parabolicity"""
        c0 = 2*(dx/hbar)**2
        # boundary conditions
        psi0 = 0.0                 
        psi1 = 1.0
        psi2 = None
        for j in range(1,n_max-1,1): # Last potential not used
            c1 = 2.0 / (cb_meff[j]*(1.0 + cb_meff_alpha[j]*(E - fis[j])) +\
                        cb_meff[j-1]*(1.0 + cb_meff_alpha[j-1]*(E - fis[j-1])) ) 
            c2 = 2.0 / (cb_meff[j]*(1.0 + cb_meff_alpha[j]*(E - fis[j])) +\
                        cb_meff[j+1]*(1.0 + cb_meff_alpha[j+1]*(E - fis[j+1])) )       
            psi2=((c0*(fis[j]-E)+c2+c1)*psi1-c1*psi0)/c2
            psi0=psi1
            psi1=psi2
        return psi2

    def psi_at_inf2b(E,fis,cb_meff,cb_meff_alpha,n_max,dx):
        """shooting method with non-parabolicity"""
        cb_meff_E = np.array(cb_meff)*(1.0 + np.array(cb_meff_alpha)*(E-np.array(fis))) # energy dependent mass
        cb_meff_E = cb_meff_E.tolist()
        return psi_at_inf1(E,fis,cb_meff_E,cb_meff_alpha,n_max,dx)

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
    E_state=[0.0]*numlevels #Energies of subbands (meV)
    cb_meff = model.cb_meff # effective mass of electrons in conduction band (kg)
    cb_meff_alpha = model.cb_meff_alpha # non-parabolicity constant for conduction band mass.
    energyx = float(energyx0) #starting energy for subband search (Joules) + floats are faster than numpy.float64
    n_max = model.n_max
    dx = model.dx
    
    if config.use_cython == False: #lists are faster than numpy arrays for loops
        fi = fi.tolist()    
        cb_meff = cb_meff.tolist()
        cb_meff_alpha = cb_meff_alpha.tolist()
    
    #print 'energyx', energyx,type(energyx)
    #print 'cb_meff', cb_meff[0:10], type(cb_meff), type(cb_meff[0])
    #print 'n_max', n_max, type(n_max)
    #print 'fi', fi[0:10], type(fi), type(fi[0])
    #print 'dx', dx, type(dx)
    #exit()
    #choose shooting function
    if model.comp_scheme in (1,3): #then include non-parabolicity calculation
        psi_at_inf = psi_at_inf2
    else:
        psi_at_inf = psi_at_inf1
    
    for i in range(0,numlevels,1):  
        #increment energy-search for f(x)=0
        y2=psi_at_inf(energyx,fi,cb_meff,cb_meff_alpha,n_max,dx)
        while True:
            y1=y2
            energyx += delta_E
            y2=psi_at_inf(energyx,fi,cb_meff,cb_meff_alpha,n_max,dx)
            if y1*y2 < 0:
                break
        # improve estimate using midpoint rule
        energyx -= abs(y2)/(abs(y1)+abs(y2))*delta_E
        #implement Newton-Raphson method
        while True:
            y = psi_at_inf(energyx,fi,cb_meff,cb_meff_alpha,n_max,dx)
            dy = (psi_at_inf(energyx+d_E,fi,cb_meff,cb_meff_alpha,n_max,dx)- psi_at_inf(energyx-d_E,fi,cb_meff,cb_meff_alpha,n_max,dx))/(2.0*d_E)
            energyx -= y/dy
            if abs(y/dy) < 1e-9*meV2J:
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
    fis = fis.tolist() #lists are faster than numpy arrays for loops
    cb_meff = model.cb_meff.tolist() #lists are faster than numpy arrays for loops
    cb_meff_alpha = model.cb_meff_alpha.tolist() #
    n_max = model.n_max
    dx = model.dx
    #choosing effective mass function
    if model.comp_scheme in (1,3): #non-parabolicity calculation
        cb_meff_E = lambda j : cb_meff[j]*(1.0 + cb_meff_alpha[j]*(E - fis[j])) # energy dependent mass
    else:
        cb_meff_E = lambda j : cb_meff[j]
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
        c1=2.0/(cb_meff_E(j)+cb_meff_E(j-1))
        c2=2.0/(cb_meff_E(j)+cb_meff_E(j+1))
        psi[2] = ((2*(dx/hbar)**2*(fis[j]-E)+c2+c1)*psi[1]-c1*psi[0])/c2
        b[j+1]=psi[2]
        N += (psi[2])**2
        psi[0]=psi[1]
        psi[1]=psi[2]
    b2 = np.array(b)
    b2/= N**0.5
    return b2 # units of dx**0.5

    
# FUNCTIONS for FERMI-DIRAC STATISTICS-----------------------------------------   
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
    cb_meff = model.cb_meff # effective mass of conduction band across structure
    cb_meff_alpha = model.cb_meff_alpha # non-parabolicity constant across structure
    cb_meff_states = np.array([cb_meff*(1.0 + cb_meff_alpha*(E*meV2J - fi)) for E in E_state])
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
        print "Have processed all energy levels present and so can't be sure that Ef is below next higher energy level."
        logger.warning("Have processed all energy levels present and so can't be sure that Ef is below next higher energy level.")
    N_state=[0.0]*len(E_state)
    for i,(Ei,csb_meff) in enumerate(zip(E_state,meff_state)):
        Ni=(Ef - Ei)*csb_meff/(hbar**2*pi)*meV2J    # populations of levels
        Ni*=(Ni>0.0)
        N_state[i]=Ni
    return Ef,N_state #Fermi levels at 0K (meV), number of electrons in each subband at 0K
    
def fermilevel(Ntotal2d,T,E_state,meff_state):
    #find the Fermi level (meV)
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
    d_E = 1e-9 #Energy step (meV)
    while True:
        y = func(Ef,E_state,meff_state,Ntotal2d,T)
        dy = (func(Ef+d_E,E_state,meff_state,Ntotal2d,T)- func(Ef-d_E,E_state,meff_state,Ntotal2d,T))/(2.0*d_E)
        Ef -= y/dy
        if abs(y/dy) < 1e-12:
            break
    return Ef #(meV)

def calc_N_state(Ef,T,Ns,E_state,meff_state):
    # Find the subband populations, taking advantage of step like d.o.s. and analytic integral of FD
    N_state=[fd2(Ei,Ef,T)*csb_meff/(hbar**2*pi) for Ei,csb_meff in zip(E_state,meff_state)]
    return N_state # number of carriers in each subband
    
# FUNCTIONS for SELF-CONSISTENT POISSON----------------------------------------

def calc_sigma(wfe,N_state,model):
    """This function calculates `net' areal charge density
    n-type dopants lead to -ve charge representing electrons, and additionally 
    +ve ionised donors."""
    # note: model.dop is still a volume density, the delta_x converts it to an areal density
    sigma= model.dop*model.dx # The charges due to the dopant ions
    for j in range(0,model.subnumber_e,1): # The charges due to the electrons in the subbands
        sigma-= N_state[j]*(wfe[j])**2
    return sigma #charge per m**2 (units of electronic charge)
    
##
def calc_field(sigma,eps):
    # F electric field as a function of z-
    # i index over z co-ordinates
    # j index over z' co-ordinates
    # Note: sigma is a number density per unit area, needs to be converted to Couloumb per unit area
    F0 = -np.sum(q*sigma)/(2.0) #CMP'deki i ve j yer değişebilir - de + olabilir
    # is the above necessary since the total field due to the structure should be zero.
    # Do running integral
    tmp = np.hstack(([0.0],sigma[:-1])) + sigma
    tmp*= q/2.0 # Note: sigma is a number density per unit area, needs to be converted to Couloumb per unit area
    tmp[0] = F0 
    F = np.cumsum(tmp)/eps
    return F

def calc_field_convolve(sigma,eps):
    tmp = np.ones(len(sigma)-1)
    signstep = np.hstack((-tmp,[0.0],tmp)) # step function
    F = np.convolve(signstep,sigma,mode='valid')
    F*= q/(2.0*eps)
    return F

def calc_field_old(sigma,eps):
    # F electric field as a function of z-
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
    # This function calculates the potential (energy actually)
    # V electric field as a function of z-
    # i	index over z co-ordinates

    #Calculate the potential, defining the first point as zero
    tmp = q*F*dx
    V = np.cumsum(tmp) #+q -> electron -q->hole? 
    return V


# FUNCTIONS FOR EXCHANGE INTERACTION-------------------------------------------

def calc_Vxc(sigma,eps,cb_meff):
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
    nz= -(sigma - model.dop*model.dx) # electron density per m**2
    nz_3 = nz**(1/3.) #cube root of charge density.
    #a_B_eff = eps/cb_meff*a_B #effective Bohr radius
    r_s = 1.0/((4*pi/3.0)**(1/3.)*nz_3*eps/cb_meff*a_B) #average distance between charges in units of effective Bohr radis.
    #A = q**4/(32*pi**2*hbar**2)*(9*pi/4.0)**(1/3.)*2/pi*(4*pi/3.0)**(1/3.)*4*pi*hbar**2/q**2 #constant factor for expression.
    A = q**2/(4*pi)*(3/pi)**(1/3.) #simplified constant factor for expression.
    #
    Vxc = -A*nz_3/eps * ( 1.0 + 0.0545 * r_s * np.log( 1.0 + 11.4/r_s) )
    return Vxc

# -----------------------------------------------------------------------------

def Poisson_Schrodinger(model):
    """Performs a self-consistent Poisson-Schrodinger calculation of a 1d quantum well structure.
    Model is an object with the following attributes:
    fi - Bandstructure potential (J) (array, len n_max)
    cb_meff - conduction band effective mass (kg)(array, len n_max)
    eps - dielectric constant (including eps0) (array, len n_max)
    dop - doping distribution (m**-3) ( array, len n_max)
    Fapp - Applied field (Vm**-1)
    T - Temperature (K)
    comp_scheme - simulation scheme (currently unused)
    subnumber_e - number of subbands for look for in the conduction band
    dx - grid spacing (m)
    n_max - number of points.
    """   
    fi = model.fi
    cb_meff = model.cb_meff
    eps = model.eps
    dop = model.dop
    Fapp = model.Fapp
    T = model.T
    comp_scheme = model.comp_scheme
    subnumber_e = model.subnumber_e
    dx = model.dx
    n_max = model.n_max
    
    # Check
    if comp_scheme ==6:
        print """The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere."""
        logger.warning("""The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere.""")
    
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
        print "Ntotal2d ",Ntotal2d," m**-2"
        logger.info("Ntotal2d %g m**-2" %Ntotal2d)
    
    fi_min= min(fi) #minimum potential energy of structure (for limiting the energy range when searching for states)
    
    if abs(E_start)<1e-3*meV2J: #energyx is the minimum energy (meV) when starting the search for bound states.
        energyx = fi_min
    else:
        energyx = E_start

    #Applied Field
    x0 = dx*n_max/2.0
    Vapp = q*Fapp*(xaxis-x0)

    # STARTING SELF CONSISTENT LOOP
    time2 = time.time() # timing audit
    iteration = 1   #iteration counter
    previousE0= 0   #(meV) energy of zeroth state for previous iteration(for testing convergence)
    fitot = fi + Vapp #For initial iteration sum bandstructure and applied field
    while True:
        if not(config.messagesoff) :
            print "Iteration:",iteration
            logger.info("Iteration: %d" %iteration)
        if iteration> 1:
            energyx=min(fi_min,min(fitot),)
        
        E_state=calc_E_state(subnumber_e,fitot,model,energyx)
        
        # Envelope Function Wave Functions
        for j in range(0,subnumber_e,1):
            if not(config.messagesoff) :
                print "Working for subband no:",j+1
                logger.info("Working for subband no: %d"%(j+1))
            wfe[j] = wf(E_state[j]*meV2J,fitot,model) #wavefunction units dx**0.5

        # Calculate the effective mass of each subband
        if model.comp_scheme in (1,3): #include non-parabolicity in calculation
            meff_state = calc_meff_state2(wfe,E_state,fitot,model)
        else:
            meff_state = calc_meff_state(wfe,cb_meff)
        
        ## Self-consistent Poisson
        
        # Calculate the Fermi energy and subband populations at 0K
        #E_F_0K,N_state_0K=fermilevel_0K(Ntotal2d,E_state,meff_state)
        # Calculate the Fermi energy at the temperature T (K)
        E_F = fermilevel(Ntotal2d,T,E_state,meff_state)
        # Calculate the subband populations at the temperature T (K)
        N_state=calc_N_state(E_F,T,Ntotal2d,E_state,meff_state)
        # Calculate `net' areal charge density
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
            Vnew += calc_Vxc(sigma,eps,cb_meff)
            
        #
        #status
        if not(config.messagesoff):
            for i,level in enumerate(E_state):
                print "E[",i,"]=",level,"meV" #can be written on file.
                logger.info("E[%d]= %f meV"%(i,level))
            for i,meff in enumerate(meff_state):
                print 'meff[',i,']= ',meff/m_e
                logger.info("meff[%d]= %f"%(i,meff/m_e))
            for i,Ni in enumerate(N_state):
                print 'N[',i,']= ',Ni,' m**-2'
                logger.info("N[%d]= %f m**-2"%(i,Ni))
            #print 'Efermi (at 0K) = ',E_F_0K,' meV'
            #for i,Ni in enumerate(N_state_0K):
            #    print 'N[',i,']= ',Ni
            print 'Efermi (at %gK) = ' %T, E_F,' meV'
            print "total donor charge = ",np.sum(dop)*dx,"m**-2"
            print "total level charge = ",sum(N_state),"m**-2"
            print "total system charge = ",np.sum(sigma),"m**-2"
            logger.info('Efermi (at %gK) = %g meV' %(T, E_F))
            logger.info("total donor charge = %g m**-2" %(sum(dop)*dx))
            logger.info("total level charge = %g m**-2" %(sum(N_state)))
            logger.info("total system charge = %g m**-2" %(sum(sigma)))
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
            print "Have reached maximum number of iterations"
            logger.warning("Have reached maximum number of iterations")
            break
        else:
            iteration += 1
            previousE0 = E_state[0]
            
    # END OF SELF-CONSISTENT LOOP
    time3 = time.time() # timing audit
    logger.info("calculation time  %g s" %(time3 - time2))
    
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
    
    return results

def save_and_plot(result,model):
    xaxis = result.xaxis
    
    output_directory = config.output_directory+"-numpy"
    
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
        
    def saveoutput(fname,datatuple,header=None):
        fname2 = os.path.join(output_directory,fname)
        fobj = file(fname2,'wb')
        if header: fobj.write(header+'\n')
        np.savetxt(fobj,np.column_stack(datatuple),fmt='%.6e', delimiter=' ')
        fobj.close()
        
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
        pl.figure(figsize=(10,8))
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
        pl.figure(figsize=(10,8))
        pl.suptitle('Aestimo Results')
        pl.subplot(1,1,1)
        pl.plot(xaxis, result.fitot*J2meV,'k')
        for level,state in zip(result.E_state,result.wfe): 
            pl.axhline(level,0.1,0.9,color='g',ls='--')
            pl.plot(xaxis, state*200.0+level,'b')
            #pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
        pl.axhline(result.E_F,0.1,0.9,color='r',ls='--')
        pl.xlabel('Position (m)')
        pl.ylabel('Energy (meV)')
        pl.grid(True)
        pl.show()

def QWplot(result,figno=None):
    #QW representation
    xaxis = result.xaxis
    pl.figure(figno,figsize=(10,8))
    pl.suptitle('Aestimo Results')
    pl.subplot(1,1,1)
    pl.plot(xaxis, result.fitot*J2meV,'k')
    for level,state in zip(result.E_state,result.wfe): 
        pl.axhline(level,0.1,0.9,color='g',ls='--')
        pl.plot(xaxis, state*200.0+level,'b')
        #pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
    pl.axhline(result.E_F,0.1,0.9,color='r',ls='--')
    pl.xlabel('Position (m)')
    pl.ylabel('Energy (meV)')
    pl.grid(True)
    pl.show()

if __name__=="__main__":
        
    # Import from config file
    inputfile = __import__(config.inputfilename)
    logger.info("inputfile is %s" %config.inputfilename)
    
    # Initialise structure class
    model = StructureFrom(inputfile,database)
         
    # Perform the calculation
    result = Poisson_Schrodinger(model)
    
    time4 = time.time() #timing audit
    logger.info("total running time (inc. loading libraries) %g s" %(time4 - time0))
    logger.info("total running time (exc. loading libraries) %g s" %(time4 - time1))

    
    # Write the simulation results in files
    save_and_plot(result,model)
    
    print "Simulation is finished. All files are closed."
    print "Please control the related files."
    logger.info("""Simulation is finished. All files are closed.Please control the related files.
        -----------------------------------------------------------------""")

