#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This is the 3x3 k.p aestimo calculator for valence band calculations 
   (Numpy version, there is no classic version for valence band calculations).

It can be used similarly to the aestimo.py module. aestimo_h.py can be used as 
a script or a libary.

To use as a script, define the simulation in a python file. See the following 
sample files for examples on usage and the required parameters:
    sample-qw-barrierdope-p.py
    sample-qw-barrierdope-p_cdzno.py
    sample-qw-barrierdope-p_ingran.py
    sample-multi-qw-barrierdope-p.py
    sample-multi-qw-barrierdope-p_ingran.py   
and then run aestimo on the command line as 
  ./aestimo.py -i <input file>
Since we are abusing the python module system, the input 'file' needs to be 
importable by the aestimo script. Alternatively, define the input file in the
config module using the inputfilename parameter.

To use aestimo_h.py as a library, first create an instance of the StructureFrom
class which builds the arrays describing a structure from the same input 
parameters that are found in the sample files. A simple list format is used to 
describes the structure's layers.
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
from math import log,exp,sqrt
import VBHM
from scipy import linalg
from VBHM import qsv,VBMAT1,VBMAT2,VBMAT_V,CBMAT,CBMAT_V
import config,database
# --------------------------------------
import logging
logger = logging.getLogger('aestimo')
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
        
        logger.info("Total material number in database: %d" ,(totalmaterial + totalalloy))
        
    def create_structure_arrays(self):
        """ initialise arrays/lists for structure"""
        self.N_wells_real0=sum(sum(np.char.count(self.material, 'w')))
        # Calculate the required number of grid points
        self.x_max = sum([layer[0] for layer in self.material])*1e-9 #total thickness (m)
        n_max = round2int(self.x_max/self.dx)
        # Check on n_max
        maxgridpoints = self.maxgridpoints
        mat_type= self.mat_type
        if n_max > maxgridpoints:
            logger.error("Grid number is exceeding the max number of %d",maxgridpoints)
            exit()
        #
        self.n_max = n_max
        dx =self.dx
        material_property = self.material_property
        alloy_property = self.alloy_property
        cb_meff = np.zeros(n_max)	#conduction band effective mass
        cb_meff_alpha = np.zeros(n_max) #non-parabolicity constant.
        m_hh = np.zeros(n_max)
        m_lh = np.zeros(n_max)
        m_so = np.zeros(n_max)
        # Elastic constants C11,C12
        C12 = np.zeros(n_max)		
        C11 = np.zeros(n_max)
        #Elastic constants Wurtzite C13,C33
        C13 = np.zeros(n_max)
        C33 = np.zeros(n_max)
        C44 = np.zeros(n_max)
        #Spontaneous and Piezoelectric Polarizations constants D15,D13,D33 and Psp
        D15 = np.zeros(n_max)
        D31 = np.zeros(n_max)
        D33 = np.zeros(n_max)
        Psp = np.zeros(n_max)
        # Luttinger Parameters γ1,γ2,γ3
        GA3 = np.zeros(n_max)		
        GA2 = np.zeros(n_max)            
        GA1 = np.zeros(n_max)
        #Hole eff. mass parameter  Wurtzite Semiconductors
        A1 = np.zeros(n_max)
        A2 = np.zeros(n_max)
        A3 = np.zeros(n_max)
        A4= np.zeros(n_max)
        A5 = np.zeros(n_max)
        A6 = np.zeros(n_max)
        # Lattice constant a0
        a0 = np.zeros(n_max)
        a0_wz = np.zeros(n_max)
        a0_sub =np.zeros(n_max)
        #  Deformation potentials ac,av,b		
        Ac = np.zeros(n_max)             
        Av = np.zeros(n_max)
        B = np.zeros(n_max)
        # Deformation potentials Wurtzite Semiconductors
        D1 = np.zeros(n_max)
        D2 = np.zeros(n_max)
        D3 = np.zeros(n_max)
        D4 = np.zeros(n_max)
        delta = np.zeros(n_max) #delta splitt off
        delta_so = np.zeros(n_max) #delta Spin–orbit split energy
        delta_cr = np.zeros(n_max) #delta Crystal-field split energy 
        # Strain related
        fi_h = np.zeros(n_max)    #Bandstructure potential
        fi = np.zeros(n_max)		#Bandstructure potential
        eps =np.zeros(n_max)		#dielectric constant
        dop = np.zeros(n_max)           #doping
        N_wells_real=0
        N_wells_real2=0
        N_wells_real0=self.N_wells_real0
        N_wells_virtual=N_wells_real0+2
        N_wells_virtual2=N_wells_real0+2
        Well_boundary=np.zeros((N_wells_virtual,2))
        Well_boundary2=np.zeros((N_wells_virtual,2))
        barrier_boundary=np.zeros((N_wells_virtual+1,2))
        n_max_general=np.zeros(N_wells_virtual)
        Well_boundary[N_wells_virtual-1,0]=n_max 
        Well_boundary[N_wells_virtual-1,1]=n_max
        Well_boundary2[N_wells_virtual-1,0]=n_max 
        Well_boundary2[N_wells_virtual-1,1]=n_max
        barrier_boundary[N_wells_virtual,0]=n_max
        barrier_len=np.zeros(N_wells_virtual+1)
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
                if mat_type=='Zincblende' :
                    a0_sub[startindex:finishindex]=matprops['a0']*1e-10
                    C11[startindex:finishindex] = matprops['C11'] 
                    C12[startindex:finishindex] = matprops['C12']
                    GA1[startindex:finishindex] = matprops['GA1']
                    GA2[startindex:finishindex] = matprops['GA2']
                    GA3[startindex:finishindex] = matprops['GA3']
                    Ac[startindex:finishindex] = matprops['Ac']*q
                    Av[startindex:finishindex] = matprops['Av']*q
                    B[startindex:finishindex] = matprops['B']*q
                    delta[startindex:finishindex] = matprops['delta']*q
                    fi_h[startindex:finishindex] =-(1-matprops['Band_offset'])*matprops['Eg']*q #Joule  #-0.8*q-(1-matprops['Band_offset'])*matprops['Eg']*q #Joule
                    eps[startindex:finishindex] = matprops['epsilonStatic']*eps0
                    a0[startindex:finishindex] = matprops['a0']*1e-10
                if mat_type=='Wurtzite' :
                    a0_sub[startindex:finishindex]=matprops['a0_wz']*1e-10
                    C11[startindex:finishindex] = matprops['C11']*1e10
                    C12[startindex:finishindex] = matprops['C12']*1e10
                    C13[startindex:finishindex] = matprops['C13']*1e10
                    C33[startindex:finishindex] = matprops['C33']*1e10
                    A1[startindex:finishindex] = matprops['A1']
                    A2[startindex:finishindex] = matprops['A2']
                    A3[startindex:finishindex] = matprops['A3']
                    A4[startindex:finishindex] = matprops['A4']
                    A5[startindex:finishindex] = matprops['A5']
                    A6[startindex:finishindex] = matprops['A6']
                    D1[startindex:finishindex] = matprops['D1']*q
                    D2[startindex:finishindex] = matprops['D2']*q
                    D3[startindex:finishindex] = matprops['D3']*q
                    D4[startindex:finishindex] = matprops['D4']*q                  
                    a0_wz[startindex:finishindex] = matprops['a0_wz']*1e-10 
                    delta_so[startindex:finishindex] = matprops['delta_so']*q
                    delta_cr[startindex:finishindex] = matprops['delta_cr']*q
                    eps[startindex:finishindex] = matprops['epsilonStatic']*eps0
                    fi_h[startindex:finishindex] =-(1-matprops['Band_offset'])*matprops['Eg']*q
                    Psp[startindex:finishindex]=matprops['Psp']
            elif matType in alloy_property:
                alloyprops = alloy_property[matType]
                mat1 = material_property[alloyprops['Material1']]
                mat2 = material_property[alloyprops['Material2']]
                x = layer[2] #alloy ratio
                cb_meff_alloy = x*mat1['m_e'] + (1-x)* mat2['m_e']
                cb_meff[startindex:finishindex] = cb_meff_alloy*m_e
                Eg = x*mat1['Eg'] + (1-x)* mat2['Eg']-alloyprops['Bowing_param']*x*(1-x) #eV
                fi[startindex:finishindex] = alloyprops['Band_offset']*Eg*q # for electron. Joule
                a0_sub[startindex:finishindex]=alloyprops['a0_sub']*1e-10
                if mat_type=='Zincblende':
                    C11[startindex:finishindex] = x*mat1['C11'] + (1-x)* mat2['C11']
                    C12[startindex:finishindex] = x*mat1['C12'] + (1-x)* mat2['C12']
                    GA1[startindex:finishindex] =x*mat1['GA1'] + (1-x)* mat2['GA1']
                    GA2[startindex:finishindex] = x*mat1['GA2'] + (1-x)* mat2['GA2']               
                    GA3[startindex:finishindex] = x*mat1['GA3'] + (1-x)* mat2['GA3']                
                    Ac_alloy = x*mat1['Ac'] + (1-x)* mat2['Ac']
                    Ac[startindex:finishindex] = Ac_alloy*q
                    Av_alloy = x*mat1['Av'] + (1-x)* mat2['Av']
                    Av[startindex:finishindex] = Av_alloy*q
                    B_alloy = x*mat1['B'] + (1-x)* mat2['B']
                    B[startindex:finishindex] = B_alloy*q
                    delta_alloy = x*mat1['delta'] + (1-x)* mat2['delta']
                    delta[startindex:finishindex] = delta_alloy*q
                    fi_h[startindex:finishindex] = -(1-alloyprops['Band_offset'])*Eg*q # -(-1.33*(1-x)-0.8*x)for electron. Joule-1.97793434e-20 #
                    eps[startindex:finishindex] = (x*mat1['epsilonStatic'] + (1-x)* mat2['epsilonStatic'] )*eps0
                    a0[startindex:finishindex] = ((1-x)*mat1['a0'] + x* mat2['a0'] )*1e-10
                    cb_meff_alpha[startindex:finishindex] = alloyprops['m_e_alpha']*(mat2['m_e']/cb_meff_alloy) #non-parabolicity constant for alloy. THIS CALCULATION IS MOSTLY WRONG. MUST BE CONTROLLED. SBL
                if mat_type=='Wurtzite' :
                    A1[startindex:finishindex] =x*material_property[alloyprops['Material1']]['A1'] + (1-x)* material_property[alloyprops['Material2']]['A1']
                    A2[startindex:finishindex] =x*material_property[alloyprops['Material1']]['A2'] + (1-x)* material_property[alloyprops['Material2']]['A2']
                    A3[startindex:finishindex] =x*material_property[alloyprops['Material1']]['A3'] + (1-x)* material_property[alloyprops['Material2']]['A3']
                    A4[startindex:finishindex] =x*material_property[alloyprops['Material1']]['A4'] + (1-x)* material_property[alloyprops['Material2']]['A4']
                    A5[startindex:finishindex] =x*material_property[alloyprops['Material1']]['A5'] + (1-x)* material_property[alloyprops['Material2']]['A5']
                    A6[startindex:finishindex] =x*material_property[alloyprops['Material1']]['A6'] + (1-x)* material_property[alloyprops['Material2']]['A6']
                    D1[startindex:finishindex] =(x*material_property[alloyprops['Material1']]['D1'] + (1-x)* material_property[alloyprops['Material2']]['D1'])*q
                    D2[startindex:finishindex] =(x*material_property[alloyprops['Material1']]['D2'] + (1-x)* material_property[alloyprops['Material2']]['D2'])*q
                    D3[startindex:finishindex] =(x*material_property[alloyprops['Material1']]['D3'] + (1-x)* material_property[alloyprops['Material2']]['D3'])*q
                    D4[startindex:finishindex] =(x*material_property[alloyprops['Material1']]['D4'] + (1-x)* material_property[alloyprops['Material2']]['D4'])*q
                    C13[startindex:finishindex] = (x*material_property[alloyprops['Material1']]['C13'] + (1-x)* material_property[alloyprops['Material2']]['C13'])*1e10# for newton/meter²
                    C33[startindex:finishindex] = (x*material_property[alloyprops['Material1']]['C33'] + (1-x)* material_property[alloyprops['Material2']]['C33'])*1e10
                    D31[startindex:finishindex] = x*material_property[alloyprops['Material1']]['D31'] + (1-x)* material_property[alloyprops['Material2']]['D31']
                    D33[startindex:finishindex] = x*material_property[alloyprops['Material1']]['D33'] + (1-x)* material_property[alloyprops['Material2']]['D33']
                    Psp[startindex:finishindex] = x*material_property[alloyprops['Material1']]['Psp'] + (1-x)* material_property[alloyprops['Material2']]['Psp']
                    C11[startindex:finishindex] = (x*material_property[alloyprops['Material1']]['C11'] + (1-x)* material_property[alloyprops['Material2']]['C11'])*1e10
                    C12[startindex:finishindex] = (x*material_property[alloyprops['Material1']]['C12'] + (1-x)* material_property[alloyprops['Material2']]['C12'])*1e10
                    a0_wz[startindex:finishindex] = ( x*material_property[alloyprops['Material1']]['a0_wz'] +(1-x)* material_property[alloyprops['Material2']]['a0_wz'])*1e-10
                    eps[startindex:finishindex] = (x*material_property[alloyprops['Material1']]['epsilonStatic'] + (1-x)* material_property[alloyprops['Material2']]['epsilonStatic'])*eps0
                    fi_h[startindex:finishindex] =-(1-alloyprops['Band_offset'])*Eg*q
                    delta_so[startindex:finishindex]  = (x*material_property[alloyprops['Material1']]['delta_so'] + (1-x)* material_property[alloyprops['Material2']]['delta_so'])*q
                    delta_cr[startindex:finishindex]  = (x*material_property[alloyprops['Material1']]['delta_cr'] + (1-x)* material_property[alloyprops['Material2']]['delta_cr'])*q
                    Ac_alloy = x*material_property[alloyprops['Material1']]['Ac'] + (1-x)* material_property[alloyprops['Material2']]['Ac']
                    Ac[startindex:finishindex] = Ac_alloy*q
            matRole= layer[5]
            if  matRole == 'w':
                N_wells_real2+=1
                Well_boundary2[N_wells_real2,0]=startindex 
                Well_boundary2[N_wells_real2,1]=finishindex
                
            for J in range(0,N_wells_virtual2):
                barrier_boundary[J,0]=Well_boundary2[J-1,1]
                barrier_boundary[J,1]=Well_boundary2[J,0]
                barrier_len[J]=barrier_boundary[J,1]-barrier_boundary[J,0]
            #doping
            if layer[4] == 'n':  
                chargedensity = layer[3]*1e6 #charge density in m**-3 (conversion from cm**-3)
            elif layer[4] == 'p': 
                chargedensity = -layer[3]*1e6 #charge density in m**-3 (conversion from cm**-3)
            else:
                chargedensity = 0.0
            dop[startindex:finishindex] = chargedensity
        #here we remove barriers who are less than anti_crossing_length
        #so we can constructe the new well boundary using the resulted barrier boundary
        brr=0
        anti_crossing_length=config.anti_crossing_length*1e-9
        for J in range(2,N_wells_virtual2-1):
            if (barrier_len[J]*dx <= anti_crossing_length ) :
                brr+=1
        brr_vec=np.zeros(brr)
        brr2=0
        for J in range(2,N_wells_virtual2-1):
            if (barrier_len[J]*dx <= anti_crossing_length ) :
                brr2+=1
                brr_vec[brr2-1]=J+1-brr2
        for I in range (0,brr):
            barrier_boundary=np.delete(barrier_boundary,brr_vec[I], 0)
        N_wells_virtual=N_wells_virtual-brr
        Well_boundary=np.resize(Well_boundary,(N_wells_virtual,2))
        for J in range(0,N_wells_virtual):
            Well_boundary[J-1,1]=barrier_boundary[J,0]
            Well_boundary[J,0]=barrier_boundary[J,1]
        self.fi = fi
        self.cb_meff = cb_meff
        self.cb_meff_alpha = cb_meff_alpha
        self.dop = dop
        #return fi,cb_meff,eps,dop
        self.C11 = C11
        self.C12 = C12
        self.GA1 = GA1
        self.GA2 = GA2
        self.GA3 = GA3
        self.Ac = Ac
        self.Av = Av
        self.B = B
        self.a0 = a0
        self.delta = delta
        self.fi_h = fi_h
        self.eps = eps
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.A4 = A4
        self.A5 = A5
        self.A6 = A6
        self.D1 = D1
        self.D2 = D2
        self.D3 = D3
        self.D4 = D4
        self.C13 = C13
        self.C33 = C33
        self.D31 = D31
        self.D33 = D33
        self.Psp = Psp
        self.a0_wz = a0_wz
        self.a0_sub =  a0_sub
        self.delta_so = delta_so
        self.delta_cr = delta_cr
        self.N_wells_virtual=N_wells_virtual
        self.N_wells_virtual2=N_wells_virtual2
        self.N_wells_real0=N_wells_real0
        self.Well_boundary=Well_boundary
        self.Well_boundary2=Well_boundary2
        self.barrier_boundary=barrier_boundary
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
        self.subnumber_h = inputfile.subnumber_h
        self.subnumber_e = inputfile.subnumber_e
        self.comp_scheme = inputfile.computation_scheme
        self.dx = inputfile.gridfactor*1e-9 #grid in m
        self.maxgridpoints = inputfile.maxgridpoints
        self.mat_type = inputfile.mat_type
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
                
        # Loading materials database #
        self.material_property = database.materialproperty
        totalmaterial = alen(self.material_property)
        
        self.alloy_property = database.alloyproperty
        totalalloy = alen(self.alloy_property)
        
        logger.info("Total number of materials in database: %d" %(totalmaterial+totalalloy))
        
        # Initialise arrays
        
        #cb_meff #conduction band effective mass (array, len n_max)
        #fi #Bandstructure potential (array, len n_max)
        #eps #dielectric constant (array, len n_max)
        #dop #doping distribution (array, len n_max)
        self.create_structure_arrays()


# No Shooting method parameters for Schrödinger Equation solution since we use a 3x3 KP solver
#delta_E = 1.0*meV2J #Energy step (Joules) for initial search. Initial delta_E is 1 meV. #This can be included in config as a setting?
#d_E = 1e-5*meV2J #Energy step (Joules) for Newton-Raphson method when improving the precision of the energy of a found level.
damping = 0.4    #averaging factor between iterations to smooth convergence.
max_iterations=80 #maximum number of iterations.
convergence_test=1e-6 #convergence is reached when the ground state energy (meV) is stable to within this number between iterations.

# DO NOT EDIT UNDER HERE FOR PARAMETERS
# --------------------------------------

#Vegard's law for alloys
def vegard(first,second,mole):
    return first*mole+second*(1-mole)

# FUNCTIONS for FERMI-DIRAC STATISTICS-----------------------------------------   
def fd1(Ei,Ef,model):#use
    """integral of Fermi Dirac Equation for energy independent density of states.
    Ei [meV], Ef [meV], T [K]"""
    T= model.T
    return kb*T*log(exp(meV2J*(Ei-Ef)/(kb*T))+1)

def fd2(Ei,Ef,model):
    """integral of Fermi Dirac Equation for energy independent density of states.
    Ei [meV], Ef [meV], T [K]"""
    T= model.T
    return kb*T*log(exp(meV2J*(Ef-Ei)/(kb*T))+1)
def calc_meff_state_general(wfh,wfe,model,fi,E_statec,list,m_hh,m_lh,m_so,n_max_general,j,Well_boundary):
    vb_meff= np.zeros((model.subnumber_h,n_max_general))
    #
    for i in range(0,model.subnumber_h,1):
        if list[i]=='hh1'  or  list[i]=='hh2' or  list[i]=='hh3':
           vb_meff[i]=m_hh[Well_boundary[j-1,1]:Well_boundary[j+1,0]]
        elif list[i] =='lh1'or list[i] =='lh2'or list[i] =='lh3':
           vb_meff[i]=m_lh[Well_boundary[j-1,1]:Well_boundary[j+1,0]]
        else:
           vb_meff[i]=m_so[Well_boundary[j-1,1]:Well_boundary[j+1,0]]
    tmp =1.0/np.sum(wfh[:,0:n_max_general]**2/vb_meff,axis=1)  #vb_meff[:,int(n_max/2)]
    meff_state = tmp.tolist()
    """find subband effective masses including non-parabolicity
    (but stilling using a fixed effective mass for each subband dispersion)"""
    cb_meff = model.cb_meff # effective mass of conduction band across structure
    cb_meff_alpha = model.cb_meff_alpha # non-parabolicity constant across structure
    cb_meff_states = np.array([cb_meff*(1.0 + cb_meff_alpha*(E*meV2J - fi)) for E in E_statec])
    tmp1 = 1.0/np.sum(wfe**2/cb_meff_states,axis=1)
    meff_statec = tmp1.tolist()
    return meff_statec,meff_state
def calc_meff_state(wfh,wfe,subnumber_h,subnumber_e,list,m_hh,m_lh,m_so):
    n_max=len(m_hh)
    vb_meff= np.zeros((subnumber_h,n_max))
    #
    for i in range(0,subnumber_h,1):
        if list[i]=='hh' :
           vb_meff[i]=m_hh
        elif list[i] =='lh' :
           vb_meff[i]=m_lh
        else:
           vb_meff[i]=m_so
    tmp = 1.0/np.sum(wfh**2/vb_meff,axis=1)
    meff_state = tmp.tolist()
    """find subband effective masses including non-parabolicity
    (but stilling using a fixed effective mass for each subband dispersion)"""
    cb_meff = model.cb_meff # effective mass of conduction band across structure
    #cb_meff_alpha = model.cb_meff_alpha # non-parabolicity constant across structure
    #cb_meff_states = np.array([cb_meff*(1.0 + cb_meff_alpha*(E*meV2J - fi)) for E in E_statec])
    #tmp1 = 1.0/np.sum(wfe**2/cb_meff_states,axis=1)
    tmp1 = 1.0/np.sum(wfe**2/cb_meff,axis=1)
    meff_statec = tmp1.tolist()
    return meff_statec,meff_state

def fermilevel_0Kc(Ntotal2d,E_statec,meff_statec,model):#use
    Et2,Ef=0.0,0.0
    meff_statec=np.array(meff_statec)
    E_statec=np.array(E_statec)
    for i in range (model.subnumber_e,0,-1): #,(Ei,vsb_meff) in enumerate(zip(E_state,meff_state)):
        Efnew2=(sum(E_statec[0:i]*meff_statec[0:i]))
        m2=(sum(meff_statec[0:i]))
        Et2+=E_statec[i-model.subnumber_e]
        Efnew=(Efnew2+Ntotal2d*hbar**2*pi*J2meV)/(m2)
        if  Efnew>Et2:
            Ef=Efnew
            #print 'Ef[',i-subnumber_h,']=',Ef
        else:
            break #we have found Ef and so we should break out of the loop
    else: #exception clause for 'for' loop.
        logger.warning("Have processed all energy levels present and so can't be sure that Ef is below next higher energy level.")

    #Ef1=(sum(E_state*meff_state)-Ntotal2d*hbar**2*pi)/(sum(meff_state))
    N_statec=[0.0]*len(E_statec)
    for i,(Ei,csb_meff) in enumerate(zip(E_statec,meff_statec)):
        Nic=(Ef - Ei)*csb_meff/(hbar**2*pi)*meV2J    # populations of levels
        Nic*=(Nic>0.0)
        N_statec[i]=Nic
    return Ef,N_statec #Fermi levels at 0K (meV), number of electrons in each subband at 0K

def fermilevel_0K(Ntotal2d,E_state,meff_state,model):#use
    Et1,Ef=0.0,0.0
    E_state=np.array(E_state)
    for i in range (model.subnumber_h,0,-1): #,(Ei,vsb_meff) in enumerate(zip(E_state,meff_state)):
        Efnew1=(sum(E_state[0:i]*meff_state[0:i]))
        m1=(sum(meff_state[0:i]))
        Et1+=E_state[i-model.subnumber_h]
        Efnew=(Efnew1+Ntotal2d*hbar**2*pi*J2meV)/(m1)
        if Efnew<Et1 :
            Ef=Efnew
            #print 'Ef[',i-subnumber_h,']=',Ef
        else:
            break #we have found Ef and so we should break out of the loop
    else: #exception clause for 'for' loop.
        logger.warning("Have processed all energy levels present and so can't be sure that Ef is below next higher energy level.")
    
    #Ef1=(sum(E_state*meff_state)-Ntotal2d*hbar**2*pi)/(sum(meff_state))    
    N_state=[0.0]*len(E_state)
    for i,(Ei,vsb_meff) in enumerate(zip(E_state,meff_state)):
        Ni=(Ei - Ef)*vsb_meff/(hbar**2*pi)*meV2J    # populations of levels
        Ni*=(Ni>0.0)
        N_state[i]=Ni
    return Ef,N_state  #Fermi levels at 0K (meV), number of electrons in each subband at 0K

def fermilevel(Ntotal2d,model,E_state,E_statec,meff_state,meff_statec):#use
    #find the Fermi level (meV)
    def func(Ef,E_state,meff_state,E_statec,meff_statec,Ntotal2d,model):
        #return Ntotal2d - sum( [vsb_meff*fd2(Ei,Ef,T) for Ei,vsb_meff in zip(E_state,meff_state)] )/(hbar**2*pi)
        diff,diff1,diff2=0.0,0.0,0.0
        diff = Ntotal2d
        for Ei,csb_meff in zip(E_statec,meff_statec):
            diff1 -= csb_meff*fd2(Ei,Ef,model)/(hbar**2*pi)
        for Ei,vsb_meff in zip(E_state,meff_state):
            diff2 += vsb_meff*fd1(Ei,Ef,model)/(hbar**2*pi)
        if Ntotal2d>0 :
            diff+=diff1
        else:
            diff+=diff2
        return diff
    if Ntotal2d>0 :
        Ef_0K,N_states_0K = fermilevel_0Kc(Ntotal2d,E_statec,meff_statec,model)
    else:
        Ef_0K,N_states_0K= fermilevel_0K(Ntotal2d,E_state,meff_state,model)
    #Ef=fsolve(func,Ef_0K,args=(E_state,meff_state,Ntotal2d,T))[0]
    #return float(Ef)
    #implement Newton-Raphson method
    Ef =Ef_0K
    #itr=0
    logger.info('Ef (at 0K)= %g',Ef)
    d_E = 1e-9 #Energy step (meV)
    while True:
        y = func(Ef,E_state,meff_state,E_statec,meff_statec,Ntotal2d,model)
        dy = (func(Ef+d_E,E_state,meff_state,E_statec,meff_statec,Ntotal2d,model)- func(Ef-d_E,E_state,meff_state,E_statec,meff_statec,Ntotal2d,model))/(2.0*d_E)
        if dy == 0.0: #increases interval size for derivative calculation in case of numerical error
            d_E*=2.0
            continue
        Ef -= y/dy
        if abs(y/dy) < 1e-12:
            break
        for i in range(2):
            if d_E>1e-9:
                d_E*=0.5
    return Ef #(meV)

def calc_N_state(Ef,model,E_state,meff_state,E_statec,meff_statec,Ntotal2d):#use
    # Find the subband populations, taking advantage of step like d.o.s. and analytic integral of FD
    N_statec,N_state=0.0,0.0
    if Ntotal2d>0 :
        N_statec=[fd2(Ei,Ef,model)*csb_meff/(hbar**2*pi) for Ei,csb_meff in zip(E_statec,meff_statec)]
    else:
        N_state=[fd1(Ei,Ef,model)*vsb_meff/(hbar**2*pi) for Ei,vsb_meff in zip(E_state,meff_state)]
    return N_state,N_statec  # number of carriers in each subband
    
# FUNCTIONS for SELF-CONSISTENT POISSON--------------------------------

def calc_sigma(wfh,wfe,N_state,N_statec,model,Ntotal2d): #use
    """This function calculates `net' areal charge density
    n-type dopants lead to -ve charge representing electrons, and additionally 
    +ve ionised donors."""
    # note: model.dop is still a volume density, the delta_x converts it to an areal density
    sigma= model.dop*model.dx # The charges due to the dopant ions
    if Ntotal2d>0 :
        for j in range(0,model.subnumber_e,1): # The charges due to the electrons in the subbands
            sigma-= N_statec[j]*(wfe[j])**2
    else:
        for i in range(0,model.subnumber_h,1): # The charges due to the electrons in the subbands
            sigma+= N_state[i]*(wfh[i])**2
    return sigma #charge per m**2 (units of electronic charge)
def calc_sigma_general(wfh,wfe,N_state,N_statec,model,Ntotal2d,j,Well_boundary): #use
    """This function calculates `net' areal charge density
    n-type dopants lead to -ve charge representing electrons, and additionally 
    +ve ionised donors."""
    # note: model.dop is still a volume density, the delta_x converts it to an areal density
    sigma= model.dop[Well_boundary[j-1,1]:Well_boundary[j+1,0]]*model.dx # The charges due to the dopant ions
    if Ntotal2d>0 :
        for j in range(0,model.subnumber_e,1): # The charges due to the electrons in the subbands
            sigma-= N_statec[j]*(wfe[j])**2
    else:
        for i in range(0,model.subnumber_h,1): # The charges due to the electrons in the subbands
            sigma+= N_state[i]*(wfh[i])**2
    return sigma #charge per m**2 (units of electronic charge)    
##
def calc_field(sigma,eps):
    # F electric field as a function of z-
    # i index over z co-ordinates
    # j index over z' co-ordinates
    # Note: sigma is a number density per unit area, needs to be converted to Couloumb per unit area
    sigma=sigma
    F0 = -np.sum(q*sigma)/(2.0) #CMP'deki i ve j yer değişebilir - de + olabilir
    # is the above necessary since the total field due to the structure should be zero.
    # Do running integral
    tmp = np.hstack(([0.0],sigma[:-1])) + sigma
    tmp*= q/2.0 # Note: sigma is a number density per unit area, needs to be converted to Couloumb per unit area
    tmp[0] = F0 
    F = np.cumsum(tmp)/eps
    return F

def calc_field_convolve(sigma,eps):#use
    tmp = np.ones(len(sigma)-1)
    signstep = np.hstack((-tmp,[0.0],tmp)) # step function
    F = np.convolve(signstep,sigma,mode='valid')
    F*= q/(2.0*eps)
    return F

def calc_field_old(sigma,eps):#use
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

def calc_potn(F,model):#use
    # This function calculates the potential (energy actually)
    # V electric field as a function of z-
    # i	index over z co-ordinates

    #Calculate the potential, defining the first point as zero
    tmp = q*F*model.dx
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

# valence band effective-masses hh,lh

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
    subnumber_h = model.subnumber_h
    N_wells_real0=model.N_wells_real0
    subnumber_e = model.subnumber_e
    dx = model.dx
    n_max = model.n_max
    mat_type= model.mat_type
    
    if comp_scheme in (4,5,6):
        logger.error("""aestimo_h doesn't currently include exchange interactions
        in its valence band calculations.""")
        exit()
    if comp_scheme in (1,3,6):
        logger.error("""aestimo_h doesn't currently include nonparabolicity effects in 
        its valence band calculations.""")
        exit()
    
    C11 = model.C11
    C12 = model.C12
    GA1 = model.GA1
    GA2 = model.GA2
    GA3 = model.GA3
    Ac = model.Ac
    Av = model.Av
    B = model.B
    a0 = model.a0
    delta = model.delta
    fi_h = model.fi_h
    A1 = model.A1
    A2 = model.A2
    A3 = model.A3
    A4 = model.A4
    A5 = model.A5
    A6 = model.A6
    D1 = model.D1
    D2 = model.D2
    D3 = model.D3
    D4 = model.D4
    C13 = model.C13
    C33 = model.C33
    D31 = model.D31
    D33 = model.D33
    Psp = model.Psp
    a0_wz = model.a0_wz
    a0_sub = model.a0_sub
    delta_so = model.delta_so
    delta_cr = model.delta_cr
    N_wells_virtual = model.N_wells_virtual
    N_wells_virtual2 = model.N_wells_virtual2
    Well_boundary=model.Well_boundary
    Well_boundary2=model.Well_boundary2
    barrier_boundary=model.barrier_boundary
    Ppz= np.zeros(n_max)
    HUPMAT1=np.zeros((n_max*3, n_max*3))
    HUPMATC1=np.zeros((n_max, n_max))
    UNIM = np.identity(n_max)
    EXX  = np.zeros(n_max)
    EZZ  = np.zeros(n_max)
    ZETA= np.zeros(n_max)
    CNIT= np.zeros(n_max)
    VNIT= np.zeros(n_max)
    S= np.zeros(n_max)
    k1= np.zeros(n_max)
    k2= np.zeros(n_max)
    k3= np.zeros(n_max)
    fp= np.ones(n_max)
    fm= np.ones(n_max)
    EPC= np.zeros(n_max)
    m_hh = np.zeros(n_max)
    m_lh = np.zeros(n_max)
    m_so = np.zeros(n_max)
    x_max=dx*n_max
    if config.strain :
        if mat_type=='Zincblende' :
            EXX= (a0_sub-a0)/a0
            EZZ= -2.0*C12/C11*EXX
            ZETA= -B/2.0*(EXX+EXX-2.0*EZZ)
            CNIT= Ac*(EXX+EXX+EZZ)
            VNIT= -Av*(EXX+EXX+EZZ)           
        if mat_type=='Wurtzite':
            EXX= (a0_sub-a0_wz)/a0_wz
            EZZ=-2.0*C13/C33*EXX
            CNIT= Ac*(EXX+EXX+EZZ)
            ZETA= (D2*(EXX+EXX)+D1*EZZ)
            VNIT= (D4*(EXX+EXX)+D3*EZZ)
            Ppz=((D31*(C11+C12)+D33*C13)*(EXX+EXX)+(2*D31*C13+D33*C33)*(EZZ))
            dx=x_max/n_max
            sum_1=0.0
            sum_2=0.0
            if config.piezo:
                """ spontaneous and piezoelectric polarization built-in field 
                [1] F. Bernardini and V. Fiorentini phys. stat. sol. (b) 216, 391 (1999)
                [2] Book 'Quantum Wells,Wires & Dots', Paul Harrison, pages 236-241"""
                for J in range(1,N_wells_virtual2-1):
                    BW=Well_boundary2[J,0]
                    WB=Well_boundary2[J,1]                    
                    Lw=(WB-BW)*dx
                    lb1=(BW-Well_boundary2[J-1,1])*dx
                    lb2=(Well_boundary2[J+1,0]-WB)*dx
                    sum_1+=(Psp[BW+1]+Ppz[BW+1])*Lw/eps[BW+1]+(Psp[BW-1]+Ppz[BW-1])*lb1/eps[BW-1]
                    sum_2+=Lw/eps[BW+1]+lb1/eps[BW-1]
                EPC=(sum_1-(Psp+Ppz)*sum_2)/(eps*sum_2)
                """
                sum_1=0.0
                sum_2=0.0
                sum_1=sum((Psp+Ppz)/eps)*dx
                sum_2=sum(1/eps)*dx
                EPC=(sum_1-(Psp+Ppz)*sum_2)/(eps*sum_2)
                """

    if mat_type=='Zincblende' :
        for i in range(0,n_max,1):
            if  EXX[i]!=0: 
                S[i]=ZETA[i]/delta[i]
                k1[i]=sqrt(1+2*S[i]+9*S[i]**2)
                k2[i]=S[i]-1+k1[i]
                k3[i]=S[i]-1-k1[i]
                fp[i]=(2*S[i]*(1+1.5*k2[i])+6*S[i]**2)/(0.75*k2[i]**2+k2[i]-3*S[i]**2)
                fm[i]=(2*S[i]*(1+1.5*k3[i])+6*S[i]**2)/(0.75*k3[i]**2+k3[i]-3*S[i]**2)
        m_hh = m_e/(GA1 -2*GA2 )
        m_lh = m_e/(GA1 +2*fp*GA2 )
        m_so = m_e/(GA1 +2*fm*GA2 )            
    if mat_type=='Wurtzite' :
        m_hh = -m_e/(A2 + A4 -A5)
        m_lh = -m_e/(A2 + A4 +A5 )
        m_so = -m_e/(A2)    
    RATIO=m_e/hbar**2*(x_max)**2
    AC1=(n_max+1)**2    
    AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,FSO,Pce,GDELM,DEL3,DEL1,DEL2=qsv(GA1,GA2,GA3,RATIO,VNIT,ZETA,CNIT,AC1,n_max,delta,A1,A2,A3,A4,A5,A6,delta_so,delta_cr,mat_type)
    KP=0.0
    KPINT=0.01
    if mat_type=='Zincblende' :        
        HUPMAT1=VBMAT1(KP,AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,FSO,GDELM,x_max,n_max,AC1,UNIM,KPINT)
    if mat_type=='Wurtzite' :
        HUPMAT1=-VBMAT2(KP,AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,x_max,n_max,AC1,UNIM,KPINT,DEL3,DEL1,DEL2)
    HUPMATC1=CBMAT(KP,Pce,cb_meff/m_e,x_max,n_max,AC1,UNIM,KPINT)
    def calc_E_state(HUPMAT1,HUPMATC1,subnumber_h,subnumber_e,fitot,fitotc):
        HUPMAT3=np.zeros((n_max*3, n_max*3))
        HUPMAT3=VBMAT_V(HUPMAT1,fitot,RATIO,n_max,UNIM)
        HUPMATC3=CBMAT_V(HUPMATC1,fitotc,RATIO,n_max,UNIM)
        #stop
        KPV1=[0.0]*subnumber_e
        la1,v1= linalg.eigh(HUPMATC3)
        tmp1=la1/RATIO*J2meV
        tmp1=tmp1.tolist()
        for i in range(0,subnumber_e,1):
            KPV1[i]=tmp1[i]
        KPV2=[0.0]*subnumber_h 
        la2,v2= linalg.eigh(HUPMAT3) 
        tmp=-la2/RATIO*J2meV
        tmp=tmp.tolist()
        for i in range(0,subnumber_h,1):
            KPV2[i]=tmp[i]
        return KPV1,v1,KPV2,v2
    KPV1=np.zeros((model.N_wells_virtual,subnumber_e))
    V1=np.zeros((model.N_wells_virtual,n_max,n_max))
    KPV2=np.zeros((model.N_wells_virtual,subnumber_h))
    V2=np.zeros((model.N_wells_virtual,n_max*3,n_max*3))
    n_max_general=np.zeros(model.N_wells_virtual)
    def calc_E_state_general(HUPMAT1,HUPMATC1,subnumber_h,subnumber_e,fitot,fitotc,model,Well_boundary):
        n_max_general=np.zeros(model.N_wells_virtual)
        HUPMAT3=np.zeros((n_max*3, n_max*3))
        HUPMAT3=VBMAT_V(HUPMAT1,fitot,RATIO,n_max,UNIM)
        HUPMATC3=CBMAT_V(HUPMATC1,fitotc,RATIO,n_max,UNIM)
        #stop
        tmp1=np.zeros((model.N_wells_virtual,n_max))
        KPV1=np.zeros((model.N_wells_virtual,subnumber_e))
        V1=np.zeros((model.N_wells_virtual,n_max,n_max))
        for J in range(1,model.N_wells_virtual-1):
            n_max_general[J]=Well_boundary[J+1,0]-Well_boundary[J-1,1]
            I1=int(Well_boundary[J-1,1])
            I2=int(Well_boundary[J+1,0])
            i1=I1-I1
            i2=I2-I1
            la1,v1= linalg.eigh(HUPMATC3[I1:I2,I1:I2])
            tmp1[J,i1:i2]=la1/RATIO*J2meV
            V1[J,i1:i2,i1:i2]=v1
        for j in range(1,model.N_wells_virtual-1):            
            for i in range(0,subnumber_e,1):
                KPV1[j,i]=tmp1[j,i]
        tmp=np.zeros((model.N_wells_virtual,n_max*3))
        KPV2=np.zeros((model.N_wells_virtual,subnumber_h))
        V2=np.zeros((model.N_wells_virtual,n_max*3,n_max*3))
        for k in range(1,model.N_wells_virtual-1):
            i_1=int(n_max_general[k])
            HUPMAT3_general=np.zeros((i_1*3,i_1*3))
            I1=int(Well_boundary[k-1,1])
            I2=int(Well_boundary[k+1,0])
            i1=I1-I1
            i2=I2-I1
            HUPMAT3_general[i1:i2,i1:i2]=HUPMAT3[I1:I2,I1:I2]
            HUPMAT3_general[i1+i_1:i2+i_1,i1:i2]=HUPMAT3[I1+n_max:I2+n_max,I1:I2]
            HUPMAT3_general[i1:i2,i1+i_1:i2+i_1]=HUPMAT3[I1:I2,I1+n_max:I2+n_max]
            HUPMAT3_general[i1+i_1:i2+i_1,i1+i_1:i2+i_1]=HUPMAT3[I1+n_max:I2+n_max,I1+n_max:I2+n_max]
            HUPMAT3_general[i1+i_1*2:i2+i_1*2,i1:i2]=HUPMAT3[I1+n_max*2:I2+n_max*2,I1:I2]
            HUPMAT3_general[i1:i2,i1+i_1*2:i2+i_1*2]=HUPMAT3[I1:I2,I1+n_max*2:I2+n_max*2]
            HUPMAT3_general[i1+i_1*2:i2+i_1*2,i1+i_1*2:i2+i_1*2]=HUPMAT3[I1+n_max*2:I2+n_max*2,I1+n_max*2:I2+n_max*2]
            HUPMAT3_general[i1+i_1:i2+i_1,i1+i_1*2:i2+i_1*2]=HUPMAT3[I1+n_max:I2+n_max,I1+n_max*2:I2+n_max*2]
            HUPMAT3_general[i1+i_1*2:i2+i_1*2,i1+i_1:i2+i_1]=HUPMAT3[I1+n_max*2:I2+n_max*2,I1+n_max:I2+n_max]
            la2,v2= linalg.eigh(HUPMAT3_general) 
            tmp[k,i1:i2*3]=-la2/RATIO*J2meV
            V2[k,i1:i2*3,i1:i2*3]=v2
        #tmp=tmp.tolist()
        for j in range(1,model.N_wells_virtual-1):            
            for i in range(0,subnumber_h,1):
                KPV2[j,i]=tmp[j,i]
        return KPV1,V1,KPV2,V2    
    
    # Check
    if comp_scheme ==6:
        logger.warning("""The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere.""")
    
    # Preparing empty subband energy lists.
    E_state = [0.0]*subnumber_h     # Energies of subbands/levels (meV)
    N_state = [0.0]*subnumber_h     # Number of carriers in subbands  
    E_statec = [0.0]*subnumber_e     # Energies of subbands/levels (meV)
    N_statec = [0.0]*subnumber_e     # Number of carriers in subbands
    # Preparing empty subband energy arrays for multiquantum wells.
    E_state_general = np.zeros((model.N_wells_virtual,subnumber_h))     # Energies of subbands/levels (meV)
    N_state_general = np.zeros((model.N_wells_virtual,subnumber_h))     # Number of carriers in subbands  
    E_statec_general = np.zeros((model.N_wells_virtual,subnumber_e))     # Energies of subbands/levels (meV)
    N_statec_general = np.zeros((model.N_wells_virtual,subnumber_e))     # Number of carriers in subbands
    meff_statec_general= np.zeros((model.N_wells_virtual,subnumber_e))
    meff_state_general= np.zeros((model.N_wells_virtual,subnumber_h))
    # Creating and Filling material arrays
    xaxis = np.arange(0,n_max)*dx   #metres
    fitot = np.zeros(n_max)         #Energy potential = Bandstructure + Coulombic potential
    fitotc = np.zeros(n_max)         #Energy potential = Bandstructure + Coulombic potentia
    #eps = np.zeros(n_max+2)	    #dielectric constant
    #dop = np.zeros(n_max+2)	    #doping distribution
    #sigma = np.zeros(n_max+2)      #charge distribution (donors + free charges)
    #F = np.zeros(n_max+2)          #Electric Field
    #Vapp = np.zeros(n_max+2)       #Applied Electric Potential
    V = np.zeros(n_max)             #Electric Potential

    # Subband wavefunction for holes list. 2-dimensional: [i][j] i:stateno, j:wavefunc
    wfh = np.zeros((subnumber_h,n_max))
    wfe = np.zeros((subnumber_e,n_max))
    wfh_general = np.zeros((model.N_wells_virtual,subnumber_h,n_max))
    wfe_general = np.zeros((model.N_wells_virtual,subnumber_e,n_max))
    E_F_general= np.zeros(model.N_wells_virtual)
    sigma_general = np.zeros(n_max)
    F_general = np.zeros(n_max)
    Vnew_general = np.zeros(n_max)
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
    fitot = fi_h + Vapp #For initial iteration sum bandstructure and applied field
    fitotc = fi + Vapp
    while True:
        if not(config.messagesoff) :
            logger.info("Iteration: %d", iteration)
        #HUPMAT2=np.zeros((n_max*3, n_max*3))
        E_statec_general,V1,E_state_general,V2=calc_E_state_general(HUPMAT1,HUPMATC1,subnumber_h,subnumber_e,fitot,fitotc,model,Well_boundary)
        for j in range(1,model.N_wells_virtual-1):
            #
            # Envelope Function Wave Functions
            n_max_general[j]=Well_boundary[j+1,0]-Well_boundary[j-1,1]
            wfh1s = np.zeros((subnumber_h,3,n_max_general[j]))
            maxwfh = np.zeros((subnumber_h,3))
            list = ['']*subnumber_h
            for i in range(0,subnumber_e,1):
                wfe_general[j,i,0:n_max_general[j]] = V1[j,0:n_max_general[j],i]                     
            wfh_pow=np.zeros(n_max)
            conter_hh,conter_lh,conter_so=0,0,0
            for jj in range(0,subnumber_h):
                for i in range(0,3):
                    wfh1s[jj,i,:] = V2[j,i*n_max_general[j]:(i+1)*n_max_general[j],jj]
                    wfh_pow=np.cumsum(wfh1s[jj,i,:]*wfh1s[jj,i,:])
                    maxwfh[jj,i]=wfh_pow[n_max_general[j]-1]
                if np.argmax(maxwfh[jj,:])==0 :
                    conter_hh+=1       
                    list[jj]='hh%d'%conter_hh                               
                    wfh_general[j,jj,0:n_max_general[j]]=wfh1s[jj,np.argmax(maxwfh[jj,:]),:]
                elif np.argmax(maxwfh[jj,:])==1:
                    conter_lh+=1       
                    list[jj]='lh%d'%conter_lh
                    wfh_general[j,jj,0:n_max_general[j]]=wfh1s[jj,np.argmax(maxwfh[jj,:]),:]
                else:
                    conter_so+=1       
                    list[jj]='so%d'%conter_so
                    wfh_general[j,jj,0:n_max_general[j]]=wfh1s[jj,np.argmax(maxwfh[jj,:]),:]
            meff_statec,meff_state = calc_meff_state_general(wfh_general[j,:,:],wfe_general[j,:,:],model,fitotc,E_statec_general[j,:],list,m_hh,m_lh,m_so,n_max_general[j],j,Well_boundary)
            meff_statec_general[j,:],meff_state_general[j,:] =meff_statec,meff_state
            E_F = fermilevel(Ntotal2d,model,E_state_general[j],E_statec_general[j],meff_state,meff_statec)
            E_F_general[j]=E_F
            # Calculate the subband populations at the temperature T (K)
            N_state,N_statec=calc_N_state(E_F,model,E_state_general[j,:],meff_state,E_statec_general[j,:],meff_statec,Ntotal2d)
            N_state_general[j,:],N_statec_general[j,:]=N_state,N_statec
            # Calculate `net' areal charge density
            sigma=calc_sigma_general(wfh_general[j,:,0:n_max_general[j]],wfe_general[j,:,0:n_max_general[j]],N_state,N_statec,model,Ntotal2d,j,Well_boundary) #one more instead of subnumber_h
            sigma_general[Well_boundary[j-1,1]:Well_boundary[j+1,0]]=sigma            
            #status
            if not(config.messagesoff):
                logger.info("-----------Starting calculation for Quantum Region number %d ------------------",j)
                for i,level in enumerate(E_state_general[j]):
                    logger.info("E[%d]= %f meV",i,level)
                for i,meff in enumerate(meff_state):
                    logger.info("meff[%d]= %f",i,meff/m_e)
                if Ntotal2d<0:
                    for i,Ni in enumerate(N_state):
                        logger.info("N[%d]= %g m**-2",i,Ni)
                for i,level in enumerate(E_statec_general[j]):
                    logger.info("Ec[%d]= %f meV"%(i,level))
                for i,meff in enumerate(meff_statec):
                    logger.info("meff[%d]= %f"%(i,meff/m_e))
                if Ntotal2d>0:
                    for i,Ni in enumerate(N_statec):
                        logger.info("N[%d]= %g m**-2"%(i,Ni))
                #print 'Efermi (at 0K) = ',E_F_0K,' meV'
                #for i,Ni in enumerate(N_state_0K):
                #    print 'N[',i,']= ',Ni
                logger.info('Efermi (at %gK) = %g meV',T, E_F)
                logger.info("total donor charge = %g m**-2",sum(dop)*dx)
                if Ntotal2d>0:
                    logger.info("total level chargec = %g m**-2" %(sum(N_statec)))
                else:
                    logger.info("total level charge = %g m**-2" %(sum(N_state)))
                logger.info("total system charge = %g m**-2" %(sum(sigma)))
        # Calculate electric field (Poisson/Hartree Effects)
        #F=calc_field(sigma,eps[Well_boundary[j-1,1]:Well_boundary[j+1,0]])
        #print EPC
        F_general=calc_field(sigma_general,eps)+EPC#[Well_boundary[j-1,1]:Well_boundary[j+1,0]]=F+EPC[Well_boundary[j-1,1]:Well_boundary[j+1,0]]
        # Calculate potential due to charge distribution
        #Vnew=calc_potn(F+EPC[Well_boundary[j-1,1]:Well_boundary[j+1,0]],model)
        Vnew_general=calc_potn(F_general,model)#[Well_boundary[j-1,1]:Well_boundary[j+1,0]]=Vnew        
        #
        #
        if comp_scheme in (0,1): 
            #if we are not self-consistently including Poisson Effects then only do one loop
            break
        
        # Combine band edge potential with potential due to charge distribution
        # To increase convergence, we calculate a moving average of electric potential 
        #with previous iterations. By dampening the corrective term, we avoid oscillations.
        #fi_h=np.resize(fi_h,n_max)
        V+= damping*(Vnew_general - V)
        fitot = fi_h + V + Vapp
        fitotc = fi + V + Vapp
        if abs(E_state_general[1,0]-previousE0) < convergence_test: #Convergence test
            break
        elif iteration >= max_iterations: #Iteration limit
            logger.warning("Have reached maximum number of iterations")
            break
        else:
            iteration += 1
            previousE0 = E_state_general[1,0]
            
    # END OF SELF-CONSISTENT LOOP
    time3 = time.time() # timing audit
    logger.info("calculation time  %g s",(time3 - time2))
    
    class Results(): pass
    results = Results()
    results.N_wells_virtual=N_wells_virtual
    results.Well_boundary=Well_boundary
    results.xaxis = xaxis
    results.wfh = wfh
    results.wfe = wfe
    results.wfh_general = wfh_general
    results.wfe_general = wfe_general
    results.fitot = fitot
    results.fitotc = fitotc
    results.sigma = sigma
    results.sigma_general = sigma_general
    #results.F = F
    results.V = V
    results.E_state = E_state
    results.N_state = N_state
    results.meff_state = meff_state
    results.E_statec = E_statec
    results.N_statec = N_statec
    results.meff_statec = meff_statec
    results.F_general = F_general
    results.E_state_general = E_state_general
    results.N_state_general = N_state_general
    results.meff_state_general = meff_state_general
    results.E_statec_general = E_statec_general
    results.N_statec_general = N_statec_general
    results.meff_statec_general = meff_statec_general
    results.Fapp = Fapp
    results.T = T
    results.E_F = E_F
    results.E_F_general = E_F_general
    results.dx = dx
    results.subnumber_h = subnumber_h
    results.subnumber_e = subnumber_e
    results.Ntotal2d = Ntotal2d
    return results

def save_and_plot(result,model):
    xaxis = result.xaxis
    
    output_directory = config.output_directory+"_h"
    
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    def saveoutput(fname,datatuple,header=None):
        fname2 = os.path.join(output_directory,fname)
        fobj = file(fname2,'wb')
        if header: fobj.write(header+'\n')
        np.savetxt(fobj,np.column_stack(datatuple),fmt='%.6e', delimiter=' ')
        fobj.close()
    if result.Ntotal2d<0:
        if config.sigma_out:
            saveoutput("sigma_h.dat",(xaxis,result.sigma_general))
        if config.electricfield_out:
            saveoutput("efield_h.dat",(xaxis,result.F_general))
        if config.potential_out:
            saveoutput("potn_h.dat",(xaxis,result.fitot))
        if config.states_out:
            for j in range(1,result.N_wells_virtual-1):                
                rel_meff_state = [meff/m_e for meff in result.meff_state_general[j]] #going to report relative effective mass.
                columns = range(model.subnumber_h), result.E_state_general[j], result.N_state_general[j], rel_meff_state
                #header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
                header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
                saveoutput("states_h_QWR%d.dat" %j,columns, header = header )
                if config.probability_out:
                    saveoutput("wavefunctions_h_QWR%d.dat" %j,(xaxis,result.wfh_general[j].transpose()))
    else:
        if config.sigma_out:
            saveoutput("sigma.dat",(xaxis,result.sigma_general))
        if config.electricfield_out:
            saveoutput("efield.dat",(xaxis,result.F_general))
        if config.potential_out:
            saveoutput("potn.dat",(xaxis,result.fitot))
        if config.states_out:                
            for j in range(1,result.N_wells_virtual-1):                
                rel_meff_state = [meff/m_e for meff in result.meff_state_general[j]] #going to report relative effective mass.
                columns = range(model.subnumber_h), result.E_state_general[j], result.N_state_general[j], rel_meff_state
                #header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
                header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
                saveoutput("states_h_QW%d.dat" %j,columns, header = header )
                if config.probability_out:
                    saveoutput("wavefunctions_h_QWR%d.dat" %j,(xaxis,result.wfh_general[j].transpose()))
    # Resultviewer
        
    if config.resultviewer:
        fig1 = pl.figure(figsize=(10,8))
        pl.suptitle('Aestimo Results')
        pl.subplots_adjust(hspace=0.4,wspace=0.4)
                            
        #Plotting Sigma
        #figure(0)
        pl.subplot(2,2,1)
        pl.plot(xaxis, result.sigma_general)
        pl.xlabel('Position (m)')
        pl.ylabel('Sigma (e/m^2)')
        pl.title('Sigma')
        pl.grid(True)
    
        #Plotting Efield
        #figure(1)
        pl.subplot(2,2,2)
        pl.plot(xaxis, result.F_general)
        pl.xlabel('Position (m)')
        pl.ylabel('Electric Field strength (V/m)')
        pl.title('Electric Field')
        pl.grid(True)
    
        #Plotting Potential
        #figure(2)
        if result.Ntotal2d<0:

            pl.subplot(2,2,3)
            pl.plot(xaxis, result.fitot)
            pl.xlabel('Position (m)')
            pl.ylabel('E_c (J)')
            pl.title('Potential')
            pl.grid(True)

            #Plotting State(s)
            #figure(3)
            pl.subplot(2,2,4)
            for j in range(1,result.N_wells_virtual-1):                
                for i,state in enumerate(result.wfh_general[j,:,:]):
                    pl.plot(xaxis[result.Well_boundary[j-1,1]:result.Well_boundary[j+1,0]], state[0:result.Well_boundary[j+1,0]-result.Well_boundary[j-1,1]], label='state %d' %i)
                    pl.xlabel('Position (m)')
                    pl.ylabel('Psi')
                    pl.title('First state')
                    pl.grid(True)
        else:
            pl.subplot(2,2,3)
            pl.plot(xaxis, result.fitotc)
            pl.xlabel('Position (m)')
            pl.ylabel('E_c (J)')
            pl.title('Potential')
            pl.grid(True)

            #Plotting State(s)
            #figure(3)
            pl.subplot(2,2,4)
            for j in range(1,result.N_wells_virtual-1):                
                for i,state in enumerate(result.wfe_general[j,:,:]):
                    pl.plot(xaxis[int(result.Well_boundary[j-1,1]):int(result.Well_boundary[j+1,0])], 
                                  state[0:int(result.Well_boundary[j+1,0])-int(result.Well_boundary[j-1,1])], label='state %d' %i)
                    pl.xlabel('Position (m)')
                    pl.ylabel('Psi')
                    pl.title('First state')
                    pl.grid(True)
        #QW representation
        #figure(5)
        span=np.ones(10000)
        fig2 = pl.figure(figsize=(10,8))
        pl.suptitle('Aestimo Results')
        pl.subplot(1,1,1)
        pl.plot(xaxis,result.fitot*J2meV,'k',xaxis,result.fitotc*J2meV,'k')
        for j in range(1,result.N_wells_virtual-1):
            I1=int(result.Well_boundary[j-1,1])
            I2=int(result.Well_boundary[j+1,0])
            i1=I1-I1
            i2=I2-I1            
            for levelc,statec in zip(result.E_statec_general[j,:],result.wfe_general[j,:,:]):
                #pl.axhline(levelc,0.1,0.9, hold=True,color='g',ls='--')
                pl.plot(xaxis[I1:I2], statec[i1:i2]*200.0+levelc,'b')
                pl.plot(xaxis[I1:I2],levelc*span[I1:I2],'g',ls='--')
            for level,state in zip(result.E_state_general[j,:],result.wfh_general[j,:,:]):
                #pl.axhline(level,0.1,0.9,color='g',ls='--')
                pl.plot(xaxis[I1:I2], state[i1:i2]*config.wavefunction_scalefactor+level,'b')
                pl.plot(xaxis[I1:I2],level*span[I1:I2],'g',ls='--')
            #pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
            pl.plot(xaxis[I1:I2],result.E_F_general[j]*span[I1:I2],'r',ls='--')
        #pl.axhline(result.E_F,0.1,0.9,color='r',ls='--')
        pl.xlabel('Position (m)')
        pl.ylabel('Energy (meV)')
        pl.grid(True)
        pl.show()
    return [fig1,fig2]

def QWplot(result,figno=None):
    #QW representation
    xaxis = result.xaxis
    fig = pl.figure(figno,figsize=(10,8))
    pl.suptitle('Aestimo Results')
    pl.subplot(1,1,1)
    pl.plot(xaxis,result.fitot*J2meV,'k',xaxis,result.fitotc*J2meV,'k')
    for levelc,statec in zip(result.E_statec,result.wfe):
        pl.axhline(levelc,0.1,0.9,color='g',ls='--')
        pl.plot(xaxis, statec*200.0+levelc,'b')
    for level,state in zip(result.E_state,result.wfh):
        pl.axhline(level,0.1,0.9,color='g',ls='--')
        pl.plot(xaxis, state*config.wavefunction_scalefactor+level,'b')
        #pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
    pl.axhline(result.E_F,0.1,0.9,color='r',ls='--')
    pl.xlabel('Position (m)')
    pl.ylabel('Energy (meV)')
    pl.grid(True)
    pl.show()
    return fig


def run_aestimo(input_obj):
    """A utility function that performs the standard simulation run
    for 'normal' input files. Input_obj can be a dict, class, named tuple or 
    module with the attributes needed to create the StructureFrom class, see 
    the class implementation or some of the sample-*.py files for details."""
    logger.info("Aestimo_h is starting...")
        
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
