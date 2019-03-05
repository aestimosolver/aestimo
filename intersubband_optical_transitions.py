#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Copyright (C) 2013-2018 Sefer Bora Lisesivdin and Aestimo group

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

File Information:
-----------------
This module calculates the optical intersubband transitions (ISBTs) for the 
conduction band levels.

The module can be run as a script, it will calculate the ISBTs for the inputfile
defined in the config.py module.

To learn how to apply the module to your needs, it is best to study the code for
the plotting_absorption() function and the code within the 'if __name__=="main"'
section at the bottom of the module. The module consists of a collection of 
functions, each of which is commented and so it should be relatively easy to work
out what each function does.

If you just want to calculate the intersubband absorption though, this is some 
example code:

    results = aestimo.Poisson_Schrodinger(model) 
    #see aestimo.py for defining the model object or the code at the bottom of this module
    
    transitions_table,(hdr,units)=transitions(results,Lperiod,eps_z,linewidths)
    
    wya,Ry2a = calc_wR_multiplasmon(results,transitions_table,eps_z)
    #print 'matrix method results'; print_multiplasmon_transitions(wya,Ry2a)
    
    inv_eps_zz = inv_eps_zz_multiplasmon(wya,Ry2a,transitions_table,linewidth,freqaxis,eps_z)
    #linewidth can be a value or a function dependent upon the transition frequency.
    eps_ratio = eps_b*inv_eps_zz
    absorption = uniaxial_layer_absorption(theta,freqaxis*f2w,eps_ratio,nk,d)



Theory Notes:

Intersubband transitions (ISBTs) in a quantum well occur between the well's 
different levels but stay within a single band. This is in contrast to
interband transitions between the valence and conduction bands, those transitions
mostly occur in the visible and near-infrared and are used for LEDs, lasers and 
detectors. Intersubband transitions occur in the mid- to far- infrared and are 
used for quantum cascade lasers and mid-infrared detectors etc.

The strongest transitions are electric dipole transitions and so we need to
calculate the dipole matrix elements of the transitions. However, ISBTs have some
additional complications to a standard dipole transition. Firstly, they are 
polarisation sensitive in that they can only couple to light polarised perpendicular
to the quantum well plane. Secondly, there are complications since there are 
many electrons in the quantum well layers; one such complication is that we can
consider the transition as a plasma and the effect of the depolarisation field
induced within the plasma when interacting with an electric field leads to a 
shift of the peak's position (called the depolarisation shift). The plasma 
behaviour of the transition is also referred to as the collective excitation of 
the electrons since it the effect of the electrons on each other. 

There are even some more advanced effects that can often be ignored (and have been 
in the models implemented here). So there is another frequency shift due to the exchange 
interaction between the electrons but this is much smaller that the depolarisation 
shift. We may also need to account for the non-parabolicity of the subbands and 
the different effective masses of the different layers (and so would need to integrate 
over k-space). Finally, we may need to consider the free electron absorption from 
electrons moving within the plane of the well (for instance using a Drude model). 
Very rarely we may need to include weaker transitions such as the magnetic dipole, 
electric quadrupole.

For the collective excitation, we need to calculate of the plasma density of the 
transition. This is not as simple as it first appears, since we shouldn't just use
the width of the well to define the density but calculate an effective length
using the two level's wavefunctions. The equation is a little complicated but can
be found in ?? and is given below.

Leff_ij = hbar**2/(2*S_ijij*meff*hbar*w_ij)

where meff is the effective mass and w_ij is the transition frequency (natural) and

S_ijij = integral_0_inf( Psi_i(z).Psi_j(z).integral_0_z( integral_0_z'( Psi_i(z'').Psi_j(z'') )dz'' )dz' )dz

We can account for the collective excitation effect as it is included automatically 
if we calculate the absorption of an anisotropic absorbing layer (where I mean
absorbing like an ISBT absorbs). This can be done using classical electromagnetism
and we find that the absorption is given by

abs ~ (epsilon_b/epsilon_zz).imag . n_b.(w/c).(sin(theta)**2/cos(theta))*d

when there is nothing absorbing in-plane and the dielectric constants apart from
the ISBTs are all epsilon_b. So we see that the absorption depends upon the inverse
of the dielectric constant rather than proportionally like usual. 

For a QW containing a single dominant ISBT transition, the dielectric constant
can be given by an effective medium approach. This means that we calculate the
effective dielectric constant for the QW + its barriers ie. one period of the
structure. Since the wavefunctions of the QW penetrate the barriers, this approach
seems more appropriate than considering the QW and the barriers separately. Assuming
that the QW amd the barriers have the same dielectric constants:

eps_b/eps_z01 = 1 - (Leff/LSQW).w_p01**2.f01 / ( w01**2 +w_p01**2 - w**2 + 2j.y01.w01 )

where
w01 - transition frequency (natural) from difference in energy levels
f01 - oscillator strength = 2.meff.w01.mu01 / (q**2.hbar**2 )
mu01 - dipole matrix element = integral_0_inf( Psi_0(z).z.Psi_1(z) )dz
w_p01 - plasma frequency = dN.q**2/(meff.eps_b.eps_0)
dN - population difference between subbands (m**-2)
y01 - transition broadening (same units as w01)

Since the effective mass of each level is different in aestimo, it is not immediately clear
what value to use when calculating the parameters above that characterise the transition. 
However, in the more advanced models it turns out not to be important for the final transition
since all of the effective mass terms cancel out in the final expression.

Note on Dielectric Constants
eps_z - array of dielectric constants wrt. z for layers withing the effective medium (excluding ISBTs)
eps_b - dielectric constant of media surrounding the effective medium/QW structure
      - or sometimes just the background dielectric constant.
eps_w - dielectric constant for well layer or rather an weighted average of the dielectric
        constants seen by the ground state.

Important:
To use the results from aestimo with these functions, we need to normalise the wavefunctions
    dx = results.dx #m
    wfe = results.wfe*dx**-0.5
    
Important:
Keep track of whether you are dealing with real or natural frequencies.

"""

import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as pl
from itertools import combinations,permutations
import types
from scipy.linalg import eigh,eig

if __package__: #explicit relative imports for using aestimo as a package (in python3)
    from .aestimo import round2int,logger
else:
    from aestimo import round2int,logger

sin,cos,log,exp = np.sin,np.cos,np.log,np.exp

#Defining constants and material parameters
q = 1.602176e-19 #C
kb = 1.3806504e-23 #J/K
h = 6.62606957e-34 #m2kg/s
hbar = 1.054588757e-34
m_e = 9.1093826E-31 #kg
pi = np.pi
eps0 = 8.8541878176e-12 #F/m
c = 299792458 #m/s

J2meV=1e3/q #Joules to meV
meV2J=1e-3*q #meV to Joules
f2w = 1e12*2*pi #THz to Hz (natural)

def eig_sorted(A):
    """returns the results from scipy.eig eigenvalue solver sorted in
    ascending eigenvalue order"""
    Adiag,U = eig(A,right=True)
    order = np.argsort(Adiag)
    return Adiag[order],U[:,order]    


# Electromagnetism
# -------------------------

# Standard absorption coefficient for a transition

def absorption_standard(w,epsilon,d):
    """Calculates the absorption for a material from its dielectric 
    'constant'. However, this is not really suitable for modelling
    the absorption of intersubband transitions in quantum wells.
    w - natural frequency (Hz)
    epsilon - dielectric constant
    d - distance (m)
    """
    return 2*w*np.sqrt(epsilon).imag*d/c

# Absorption of uniaxial absorbing layer

def uniaxial_layer_absorption(theta,w,eps_ratio,nk,d):
    """Approximately calculates the absorption of a uniaxial layer
    which absorbs along its extraordinary axis which is perpendicular
    to the plane of the layer. Basically, this is the situation for 
    an intersubband transition of a quantum well. More accurate modelling
    requires a transfer matrix model that can account for these uniaxial
    layers (i.e. see pyLuminous/pyFresnel module).
    theta - angle (rad)
    w - natural frequency(ies)
    eps_ratio  = eps_b/eps_zz = eps_background/eps_layer(out of plane component)
    nk - refractive index of media surrounding uniaxial layer
    d - thickness of layer (m).
    """
    return -eps_ratio.imag*nk*w/c*sin(theta)**2/cos(theta)*d

# Lorentzian oscillator model

def susceptibility_Losc(w,w0,f,w_p,y0,eps_b=1.0):
    """Dielectric Susceptibility for the Lorentzian oscillator model
    f - oscillator strength (unitless)
    The following should all use the same units i.e Hz or meV, natural frequency
    or real frequency ...
    w_p - plasma frequency 
    w - frequency axis
    w0 - transition frequency
    y0 - transition broadening
    
    eps_b is the background dielectric constant, it is included because w_p usually
    includes eps_b in its definition and so it is necessary to include this factor
    for this formula to calculate the susceptibility. Only really interesting when using
    a frequency dependent background dielectric constant, as it slightly affects 
    how we calculate the 'total' dielectric constant. It's the difference between 
    having to using the result from this function as 
    eps_b = eps_b0.(1 + Xi + Xi2...) 
    or (by including eps_b!=1.0 here) as
    eps_b = epsb0 + Xi + Xi2 +...  
    """
    chi = (eps_b*w_p**2)*f / ( w0**2 - w**2 - 2j*y0*w )
    return chi

# Effective medium

def eff_eps_z(layer_list,isbt_term=0.0):
    """calculates the effective dielectric constant for a stack
    of thin layers for the direction perpendicular to the plane
    of the layers, the layers' thicknesses should be less than
    the wavelengths of the light spectrum considered for this to
    work effectively.
    layer_list - list of (eps_z,f) tuples where f is the fraction
    of the total thickness filled by that dielectric constant.
    isbt_term - is the contribution of the intersubband transitions
    to the effective medium (see Zaluzny PRB1999, Ando1977 and Ando1982)
    """
    eps_z = 1.0/sum([f/eps_z for eps_z,f in layer_list],isbt_term)
    return eps_z
    
def eff_eps_z2(eps_z,isbt_term=0.0):
    """calculates the effective dielectric constant for a stack
    of thin layers for the direction perpendicular to the plane
    of the layers, the layers' thicknesses should be less than
    the wavelengths of the light spectrum considered for this to
    work effectively.
    eps_z - an array of values wrt z (uniformly spaced)
          - or a 2d array of values wrt z and w (frequency axis).
    isbt_term - is the contribution of the intersubband transitions
    to the effective medium (see Zaluzny PRB1999, Ando1977 and Ando1982)
    """
    eps = 1.0/sum(np.mean(np.atleast_1d(1.0/eps_z),axis=0),isbt_term)
    return eps
    
def eff_eps_x(layer_list):
    """calculates the effective dielectric constant for a stack
    of thin layers for the direction parallel to the plane of
    the layers, the layers' thicknesses should be less than
    the wavelengths of the light spectrum considered for this to
    work effectively.
    layer_list - list of (eps_x,f) tuples where f is the fraction
    of the total thickness filled by that dielectric constant.
    """
    eps_x = sum([f*eps_x for eps_x,f in layer_list])
    return eps_x
    
def eff_eps_x2(eps_x):
    """calculates the effective dielectric constant for a stack
    of thin layers for the direction parallel to the plane of
    the layers, the layers' thicknesses should be less than
    the wavelengths of the light spectrum considered for this to
    work effectively.
    eps_z - an array of values wrt z (uniformly spaced)
          - or a 2d array of values wrt z and w (frequency axis).
    layer_list - list of (eps_x,f) tuples where f is the fraction
    of the total thickness filled by that dielectric constant.
    """
    eps = np.mean(np.atleast_1d(eps_x),axis=0)
    return eps


# Singular Electric Dipole Transitions (each transition is treated separately)
# ------------------------------------

def dipole_matrix(z,wfe1,wfe2):
    """Calculates dipole matrix element numerically. Returns values in metres
    (electron charge is not included in calculation)"""
    return simps(z*wfe1*wfe2,z)

def oscStr(z_if,w_if,cb_meff):
    """Calculates oscillator strength. 
    w_if - frequency of transition (meV),
    cb_meff - relative effective electron mass in well layer, 
    z_if - dipole matrix element for position operator (m)."""
    return 2*cb_meff*(w_if*meV2J)*z_if**2/hbar**2

def calc_S(Psi0,Psi1,Psi2,Psi3,dz):
    """Calculates S, a quantity used to calculate the effective thickness of an
    intersubband transition of a quantum well.
    Psi0 - Psi3 are arrays describing the wavefunctions.
    dz is the step-size (metres) for the arrays (assumed to be constant).
    """
    i1=0.0; i2=0.0; i3=0.0
    for p0,p1,p2,p3 in zip(Psi0,Psi1,Psi2,Psi3):
        delta = p0*p1
        i1+=delta
        i2+=i1
        i3+=p2*p3*i2
    i3*=dz**3
    return -i3

def calc_S_b(Psi0,Psi1,Psi2,Psi3,zaxis):
    """Calculates S, a quantity used to calculate the effective thickness of an
    intersubband transition of a quantum well.
    Psi0 - Psi3 are arrays describing the wavefunctions.
    zaxis is an array of z-values for the wavefunctions (needn't be uniform) (metres)
    """
    dz_axis = zaxis[1:]-zaxis[:-1]
    dz_axis = np.hstack((dz_axis[0],dz_axis)) #preprend a value so that all values get used in calculation
    i1=0.0; i2=0.0; i3=0.0
    for p0,p1,p2,p3,z,dz in zip(Psi0,Psi1,Psi2,Psi3,zaxis,dz_axis):
        delta = p0*p1
        i1+=delta
        i2+=z*delta
        i3+=p2*p3*(z*i1 - i2)*dz**2
    return -i3

def L_eff(w_if,S_if,cb_meff):
    """Calculates the effective thickness of an intersubband transition of a quantum well.
    w_if is the transitions frequency in meV, 
    cb_meff is the relative effective electron mass in well layer, 
    S is a dimensionless quantity calculated via an integral. 
    Returns a value in metres"""
    return hbar**2/(2*abs(S_if)*cb_meff*w_if*meV2J)

def calc_w_p(dn_if,cb_meff,eps_w):
    """Calculates the plasma frequency of a transition. returns real Hz.
    dn_ij - (m**-3) population density difference between the initial and final levels
    cb_meff - effective electron mass for transition
    eps_w - background dielectric constant in well layer(approximately)
    """
    return np.sqrt(dn_if*q**2/(cb_meff*eps_w*eps0))/(2*pi) #real Hz
    
def calc_R2(w_if,z_if,dn_if,eps_w,Lperiod):
    """R**2 = f*Leff*w_p**2/Lperiod (oscillator strength * plasma frequency squared)
    This is more useful for quicky seeing the strength of a transition
    than either the plasma frequency or the oscillator strength alone.
    Also we see that there is no need to know the effective mass for this expression.
    Units of R are Hz (real).
    w_if - frequency of transition (meV)
    z_if - dipole matrix element for position operator (m)
    dn_if - (m**-3) population density difference between the initial and final levels
    eps_w - background dielectric constant in well layer (approximately)
    Lperiod - thickness of effective medium, this is the QW plus some barrier each side
            (normally this is the period length in a multiple QW stack). Whatever you use
            for Lperiod, you should use for 'd' when calculating the ISBT absorptions.
    """
    return 2*dn_if*(w_if*meV2J)*(q*z_if)**2/(hbar**2*(eps_w*eps0)*eps_w*Lperiod*(2*pi)**2) 

# Summary of Transitions

def transition_generator(seq):
    """returns the possible pairs in the input sequence. Each pair is
    only returned once and the ordering found in the input is maintained"""
    return combinations(seq,2)

def transitions(results,Lperiod,eps_z,linewidths):
    """Calculates the parameters needed to describe the intersubband transitions.
    Returns a list of dictionaries (one for each transition) with the following 
    keys:
    'ilevel','flevel','dE','freq','lambda','wavno','dN','z_if','f','Leff','S_if','S_if_b','wp'
    
    results - object created by aestimo containing results of the bandstructure simulation
    Lperiod - (m) length (m), thickness of effective medium containing heterostructure 
              i.e. should be equal or larger than extent of the structures wavefunctions.
    eps_z   - (unitless) dielectric constant array wrt position. Giving dielectric constant 
              of the structure's materials at the optical frequencies of interest.
    linewidths - (THz) a number or function that gives/returns the transition linewidths. If
              using a function, it should take the transition frequency (THz) as a parameter.
    """
    E_state = results.E_state #list of energy levels (meV)
    N_state = results.N_state #occupation of energy levels (m**-2)
    meff_state = results.meff_state #effective mass of each state
    T = results.T #K
    dx = results.dx #m
    xaxis = results.xaxis #m
    wfe = results.wfe*dx**-0.5 # normalising wavefunctions to be m**-0.5
    #reversethepolarities = np.ones(wfe.shape[0])
    #reversethepolarities[1::2]*=-1
    #for j,p in enumerate(reversethepolarities):wfe[j]*=p
    
    #calculate the mean dielectric constant as for the heterostructure
    #using the lowest subband of the system.
    eps_w = 1.0/(np.sum(wfe[0]**2/eps_z,axis=0)*dx)
    
    #create linewidth function
    if type(linewidths) is types.FunctionType:
        lw = linewidths
    else:
        lw = lambda freq: linewidths
    
    
    def transition(j,i,f): #Doing it this way would let me create a dielectric function for each transition using a function closure.
        """j - transition number (useful later)
           i - initial level
           f - final level
        """
        # Most of the standard transition parameters include effective mass in their
        # formulae, Should we use the effective mass of the initial or final state, 
        # or some combination? However, currently, in the more advanced version of 
        # the dielectric constant calculation it seems that all effective mass terms 
        # cancel out. Therefore, I will use the effective mass for lowest subband for
        # all my calculations.
        cb_meff= meff_state[0] # kg
        # Likewise, many quantities include a dielectric constant in their calculation,
        # but what value should we use when a wavefunction extends over many different
        # layers. To keep things simple, we will use the mean dielectric constant as
        # calculated using the lowest subband of the system.
        eps_wi = eps_w #= sum(eps_b*wfe[0]**2)*dx
        
        dE = E_state[f]-E_state[i] #meV
        dN = N_state[i]-N_state[f] #m**2
        z_if = dipole_matrix(xaxis,wfe[i],wfe[f]) #m
        S_if = calc_S(wfe[i],wfe[f],wfe[i],wfe[f],dx) #m
        S_if_b = calc_S_b(wfe[i],wfe[f],wfe[i],wfe[f],xaxis) #m
        L_if = L_eff(dE,S_if,cb_meff) #m
        
        col = {'j':j,
               'ilevel':i,
               'flevel':f,
               'dE':dE, #meV
               'freq':dE*1e-3*q/h/1e12, #THz (real)
               'lambda':1e6*h*c/(dE*1e-3*q), #(um)
               'wavno':dE*1e-3*q/h/c*1e-2, #(cm**-1)
               'dN':dN*10**(-4-11), #dN (m-2)
               'z_if':z_if*1e9, # z (dipole matrix element) (nm)
               'f':oscStr(z_if,dE,cb_meff), #f (oscillator strength)
               'Leff':L_if*1e9, #nm
               'S_if':S_if*1e9,  #nm
               'S_if_b':S_if_b*1e9,
               'wp':calc_w_p(dN/L_if,cb_meff,eps_wi)*1e-12, #real THz
               'R':np.sqrt(calc_R2(dE,z_if,dN,eps_wi,Lperiod))*1e-12, #real THz
               'Lperiod':Lperiod*1e9, #nm
               'eps_w':eps_wi,
               }
        col['y_if'] = lw(col['freq']) #(THz real) transition broadening
        return col
    
    transitions_table = []
    for j,(i,f) in enumerate(transition_generator(np.arange(len(E_state)))):
        col = transition(j,i,f)
        transitions_table.append(col)
    
    hdr=['j','ilevel','flevel','dE','freq','lambda','wavno','dN','z_if','f','Leff','S_if','S_if_b','wp','R','Lperiod','y_if','eps_w']
    units=['','','','meV','THz','um','cm-1','1e11cm-2 @%gK'%T,'nm','','nm','nm','nm','THz','THz','nm','THz','']
    
    return transitions_table,(hdr,units)


def print_levels(results):
    """prints out energy levels and their populations to the log. Also
    print out their gaps"""
    logger.info('the energy levels\population are (meV)\t(m**-2):')
    for Ei,Ni in zip(results.E_state,results.N_state): logger.info('%.5g\t%.5g',Ei,Ni)
    logger.info('T = %gK' %results.T)
    logger.info('the energy levels gaps are')
    logger.info('\t'.join(('(meV)','(THz)','(um)','(wavno)')))
    for leveli,levelj in transition_generator(results.E_state):
        gap=levelj-leveli
        freq=gap*1e-3*q/h/1e12
        wav=1e6*h*c/(gap*1e-3*q)
        wavno=gap*1e-3*q/h/c*1e-2
        logger.info('\t'.join('%.3g' %i for i in (gap,freq,wav,wavno)))
    
def print_transitions(transitions_table,hdr,units):
    """print out summary of transition values to the log""" 
    printwidth = np.get_printoptions()['linewidth']
    var_w = 8 #print width for variable names
    unit_w = 14 #print width for units
    data_w = 11 #print width for data
    # find number of repeats needed
    cols_per_repeat = (printwidth - var_w - unit_w)//data_w
    
    def repeat_generator(n,cols_per_repeat):
        startindex = 0
        while startindex < n:
            yield slice(startindex,startindex+cols_per_repeat)
            startindex += cols_per_repeat
    
    logger.info( "Summary of Intersubband Transitions")
    for selection in repeat_generator(len(transitions_table),cols_per_repeat):
        data = transitions_table[selection]
        for var,unit in zip(hdr,units):
            row = [var.rjust(var_w),unit.rjust(unit_w)]
            row += [('%.3g' %tr[var].real).rjust(data_w) for tr in data]
            logger.info( ''.join(row))

def get_Leff_est(transitions_table):
    """gets a value of Leff for the QW that will be applied to all transitions.
    We will use the Leff of the transition with the highest oscillator strength."""
    j = np.argmax([tra['R'] for tra in transitions_table])
    return transitions_table[j]['Leff']
    
## Calculating Dielectric Constants
# Below we have some different models for intersubband transitions

def inv_eps_zz_1(transitions_table,freqaxis,eps_z):
    """calculates eps_b/eps_zz using the analytically correct result for a single transition.
    If there are several active transitions that are close together then this will
    become increasingly incorrect. The dielectric constant is calculated for the 
    effective medium of QW + barrier. However, there is an assumption here that
    eps_background = eps_barrier = eps_well_layer"""
    isbt_term = 0.0
    for trn in transitions_table: #nb. first row of table describes the units of each variable
        w_if = np.sqrt(trn['freq']**2 + trn['wp']**2) #depolarisation shifted frequency
        Xi = susceptibility_Losc(freqaxis,w0=w_if,f=1.0,w_p=trn['R'],y0=trn['y_if'])
        #print trn['R'],np.sqrt(trn['f']*trn['wp']**2*trn['Leff']/trn['Lperiod'])
        isbt_term -= Xi
    inv_eps_zz = np.mean(np.atleast_1d(1.0/eps_z),axis=0) + isbt_term
    return inv_eps_zz

def eps_classical(transitions_table,freqaxis,eps_b=1.0):
    """Approximately calculates total dielectric constant epszz for QW by summing Lorentz 
    oscillator susceptibilities for each transition. This assumes that all transitions 
    share the same effective Length (which wasn't assumed when calculating the plasma 
    frequency values (wp) contained in the transitions_table). The average effective length 
    for all of the ISBTs of the QW will have to be used as a fitting parameter in order to
    get the best fit.
    
    If you leave eps_b=1.0 then the result should be multiplied by eps_b, if you use eps_b
    then you shouldn't need to do anything. In either case it should match the values used 
    for calculating the transition plasma frequencies.
        
    warning - This shouldn't be used on its own for modelling ISBTs using absorption_standard()
    to calculate the absorption since this doesn't take into account the anisotropic nature
    of the ISBTs and the resulting depolarisation shift.
    
    """
    eps = eps_b
    for trn in transitions_table: #nb. first row of table describes the units of each variable
        Xi = susceptibility_Losc(freqaxis,w0=trn['freq'],f=trn['f'],w_p=trn['wp'],y0=trn['y_if'],eps_b=eps_b)
        eps+=Xi
    return eps

def inv_eps_zz_classical(transitions_table,freqaxis,eps_z):
    """Calculates the dielectric constant for an effective medium containing 
    intersubband transitions using a classical approach (Lorentz oscillators &
    formula for effective medium).
    Slab model?
    """
    eps_w = transitions_table[0]['eps_w']
    if True:        
        eps_qw = eps_w*eps_classical(transitions_table,freqaxis) 
    else: 
        eps_qw = eps_classical(transitions_table,freqaxis,eps_b=eps_w) #?? what's the point??
    
    Lperiod = transitions_table[0]['Lperiod'] #nm
    Lqw = get_Leff_est(transitions_table) #(nm) 
    #using the effective length for the first transition as an estimate of the thickness of the 2d electron gas contained within the QW
    ff = Lqw/Lperiod
    #inv_eps_zz = ((1.0-ff)/eps_bb + ff/eps_qw) 
    inveps_bb_term = np.mean(np.atleast_1d(1.0/eps_z),axis=0) - ff/eps_w
    inv_eps_zz = inveps_bb_term + ff/eps_qw
    #eff_eps_z(layer_list,isbt_term=0.0)
    return inv_eps_zz

## Advanced Multilevel Model #################
## frequency independent dielectric constant

def dipole_matrix_b(z,wfe1,wfe2,eps_z):
    """Calculates something similar to the dipole matrix element numerically. 
    Returns values in metres(electron charge is not included in calculation)"""
    return simps(z*wfe1*wfe2/eps_z,z)

def calc_S_c(Psi0,Psi1,Psi2,Psi3,eps_z,zaxis):
    """Calculates S, a quantity used to calculate the effective thickness of an
    intersubband transition of a quantum well.
    Psi0 - Psi3 are arrays describing the wavefunctions.
    zaxis is an array of z-values for the wavefunctions (needn't be uniform) (metres)
    """
    eps_z = eps_z*np.ones_like(zaxis)
    dz_axis = zaxis[1:]-zaxis[:-1]
    dz_axis = np.hstack((dz_axis[0],dz_axis)) #preprend a value so that all values get used in calculation
    i1=0.0; i2=0.0; i3=0.0
    for p0,p1,p2,p3,eps,dz in zip(Psi0,Psi1,Psi2,Psi3,eps_z,dz_axis):
        delta = p0*p1
        i1+=delta
        i2+=i1/eps
        i3+=p2*p3*i2
    i3*=dz**3
    return -i3

def calc_interaction_matrix(results,transitions_table,eps_z):
    """calculates the matrix of describing collective interactions between the transitions
    and also the d vector"""
    #number of transitions
    ntr = len(transitions_table) #or n*(n-1)/2 where n is the number of energy levels
    # Note that following arrays are indexed by transition rather than energy level.
    #So [0,0] is the transition for levels 0->1
    #construct transition matrix S & B
    dx = results.dx #m
    wfe = results.wfe*dx**-0.5
    zaxis = results.xaxis #m
    # S matrix
    S = np.zeros((ntr,ntr),np.complex)
    for tra in transitions_table:
        a = tra['j']
        S[a,a] = calc_S_c(wfe[tra['ilevel']],wfe[tra['flevel']],wfe[tra['ilevel']],wfe[tra['flevel']],eps_z,zaxis)
    for tra,trb in combinations(transitions_table,2):
        a = tra['j']
        b = trb['j']
        S[a,b] = S[b,a] = calc_S_c(wfe[tra['ilevel']],wfe[tra['flevel']],wfe[trb['ilevel']],wfe[trb['flevel']],eps_z,zaxis)
    #print 'S';print S
    # T matrix
    const = 2*q**2/eps0*meV2J*1e15 # 1e15 converts dN values into carriers/m**2
    R = np.zeros((ntr,ntr),np.complex)
    for tra in transitions_table:
        a = tra['j']
        R[a,a] = const*S[a,a]*np.sqrt(tra['dN']*tra['dE']*tra['dN']*tra['dE'])
    for tra,trb in combinations(transitions_table,2):
        a = tra['j']
        b = trb['j']
        R[a,b] = R[b,a] = const*S[a,b]*np.sqrt(tra['dN']*tra['dE']*trb['dN']*trb['dE']) 
    # d vector
    d = np.zeros(ntr,np.complex)
    for tra in transitions_table:
        a = tra['j'] #find correct index
        i = tra['ilevel']
        f = tra['flevel']
        x_if = dipole_matrix_b(zaxis,wfe[i],wfe[f],eps_z)
        d[a] = np.sqrt(tra['dN']*1e15*tra['dE']*meV2J)*q*x_if
    return R,d

def calc_wR_multiplasmon(results,transitions_table,eps_z):
    """Uses a multilevel version of the mathematical formalism given in Ando 1977
    A matrix is constucted describing the transitions and the interactions between
    them which can be diagonalised to give a description of the system as a simple
    sequence of Lorentzian oscillators.
    eps_z is an array of the dielectric constant wrt z for the media in the barrier+QW+barrier
    structure.
    returns (w,R-squared) - (real frequency (THz), related to transition oscillator strength)
    """
    #Calculate transitions interactions matrix + rhs of system equation
    R,d = calc_interaction_matrix(results,transitions_table,eps_z)
    #Add transition energies to Transition interaction matrix
    B = R
    for tra in transitions_table:
        a = tra['j']
        B[a,a] += (tra['dE']*meV2J)**2
    
    #diagonalise
    if np.iscomplex(eps_z).any():
        logger.info('calc_wR_multiplasmon: using eig() solver for complex symmetric or general matrix')
        Bdiag,U = eig_sorted(B) #matrix will be complex symmetric but not Hermitian, this may be a problem with the theory...
    else:
        logger.info('calc_wR_multiplasmon: using eigh() solver for Hermitian matrix')
        Bdiag,U = eigh(B, lower=True, eigvals_only=False, turbo=True, type=1) #otherwise we can be sure that B is real symmetric
    #final values of R,w0
    Ry2a = np.dot(U.transpose(),d)**2 * 2.0/(eps0*tra['Lperiod']*1e-9)*(1e-12/h)**2#THz**2 (real)
    wya = np.sqrt(Bdiag)/h*1e-12 #THz (real)
    return wya,Ry2a

def print_multiplasmon_transitions(wya,Ry2a):
    """display the results from the multiplasmon matrix model of the intersubband
    transitions. These results are 'more accurate' than the transition table results.
        wya - transition frequencies (THz)
        Ry2a - R-squared - related to transition strength - numerator of a Lorentz 
               oscillator model of the transitions. 
    
    Ry2a = f*wp**2/(eps_b)*L_eff/L
        f - oscillator strength
        wp - plasma frequency
        eps_b - background frequency
        L_eff - effective width of the transition
        L - width of QW / effective medium period.
    """
    col_width = 10
    logger.info( "Optical transitions from multiplasmon matrix model")
    logger.info( "R^2 - related to oscillator strength of each transition (THz^2).")
    logger.info( "Other columns give the transition frequencies in various units.")
    logger.info( ''.join(s.rjust(col_width) for s in ('R','(meV)','(THz)','(um)','(wavno - cm^-1)')))
    for wy,Ry2 in zip(wya,Ry2a):
        gap=wy*1e12*h*J2meV
        freq=wy
        wav=c/(wy*1e12)*1e6
        wavno=wy*1e12/c*1e-2
        logger.info( ''.join(('%.4g' %i).rjust(col_width) for i in (np.sqrt(Ry2.real),gap,freq,wav,wavno)))

def inv_eps_zz_multiplasmon(wya,Ry2a,transitions_table,linewidth,freqaxis,eps_z):
    """calculate dielectric constant ratio - 1.0/eps_ISBT for results of matrix calculation.
    linewidth is either a function depending upon a transition frequency or a constant value (THz?).
    linewidths are calculated in an empirical fashion, currently using the undepolarisation shifted
    frequency if linewidth is a function."""
    inveps = np.mean(1.0/eps_z)
    ff0 = transitions_table[0]['Leff']/transitions_table[0]['Lperiod']
    w_if = np.sort([tra['dE'] for tra in transitions_table])*meV2J/h*1e-12 #(THz) initial transition frequencies
    #w_if = np.zeros(len(transitions_table))
    #for tra in transitions_table:
    #    w_if[tra['j']] = tra['dE']*meV2J/h*1e-12 #(THz) initial transition frequencies
    for wy,Ry2,wi in zip(wya,Ry2a,w_if):
        y_y = linewidth(wi) if callable(linewidth) else linewidth #(THz real?) guesstimate of transition broadening (written to get result as close as possible to other models)
        #y_y = linewidth(np.sqrt(wy**2-Ry2/ff0)) if callable(linewidth) else linewidth #(THz real?) guesstimate of transition broadening (written to get result as close as possible to other models)
        Xi = susceptibility_Losc(freqaxis,w0=wy,f=Ry2,w_p=1.0,y0=y_y)
        inveps-= Xi
    return inveps

def inv_eps_zz_multiplasmon_helper(results,transitions_table,linewidth,freqaxis,eps_z):
    """this calculates the dielectric constant ratio - 1.0/eps_ISBT for the ISBTs for a frequency independent
    dielectric constant.
    """
    wya,Ry2a = calc_wR_multiplasmon(results,transitions_table,eps_z)
    return inv_eps_zz_multiplasmon(wya,Ry2a,transitions_table,linewdith,freqaxis,eps_z)

## frequency dependent dielectric constant (but separable from position)

def inv_eps_zz_multiplasmon2(results,transitions_table,linewidth,freqaxis,eps_z,eps_w):
    """Uses a multilevel version of the mathematical formalism given in Ando 1977
    A matrix is constucted describing the transitions and the interactions between
    them which can be diagonalised to give a description of the system as a simple
    sequence of Lorentzian oscillators.
    
    This calculates the dielectric constant ratio - 1.0/eps_ISBT for the ISBTs for
    a background dielectric constant given by eps_z * eps_w where we can separate
    out the frequency dependent part, eps_w, which is an array wrt the frequency axis.
    eps_z is an array wrt z for the media in the barrier+QW+barrier structure.
    
    linewidth - either a function of the transition frequency or a value (THz)
    freqaxis is an array of frequencies (THz) to calculate the dielectric constant for.
    eps_zw is an array of dielectric constant values (exluding eps_0) for the frequencies
        corresponding to freqaxis. It has no z-dependence, the dielectric constants of
        the barrier and well layers are assumed to be the same
    """
    #Calculate transitions interactions matrix + rhs of system equation
    R,d = calc_interaction_matrix(results,transitions_table,eps_z)
    
    #Calculate the inverse dielectric constant ############
    
    #background inverse dielectric constant
    inveps_b = np.mean(1.0/eps_z)/eps_w + 0j
    
    #choose appropriate solver
    if np.iscomplex(eps_z).any() or np.iscomplex(eps_w).any():
        logger.info('calc_wR_multiplasmon2: using eig() solver for complex symmetric or general matrix')
        eigen = lambda B: eig_sorted(B) #matrix will be complex symmetric but not Hermitian, this may be a problem with the theory...
    else:
        logger.info('calc_wR_multiplasmon2: using eigh() solver for Hermitian matrix')
        eigen = lambda B: eigh(B, lower=True, eigvals_only=False, turbo=True, type=1) #otherwise we can be sure that B is real symmetric
    
    #transition energies
    E_if = np.zeros(len(transitions_table))
    for tra in transitions_table:
        E_if[tra['j']] = tra['dE']*meV2J
    E2_if = E_if**2 #transition energies squared
    
    diag_indices = np.diag_indices_from(R) # indices to access the diagonal of the transitions interaction matrix
    
    #linewidth
    #ff0 = transitions_table[0]['Leff']/transitions_table[0]['Lperiod']
    w_if = np.sort([tra['dE'] for tra in transitions_table])*meV2J/h*1e-12 #(THz) initial transition frequencies
    #w_if = E_if/h*1e-12 #(THz) initial transition frequencies
    y_y = linewidth(w_if) if callable(linewidth) else linewidth*np.ones_like(w_if)
    #y_y = linewidth(w_i) if callable(linewidth) else linewidth #(THz real?) guesstimate of transition broadening (written to get result as close as possible to other models)
    #y_y = linewidth(np.sqrt(wy**2-Ry2/ff0)) if callable(linewidth) else linewidth #(THz real?) guesstimate of transition broadening (written to get result as close as possible to other models)
    
    Lperiod = transitions_table[0]['Lperiod']
    const_factor = 2.0/(eps0*Lperiod*1e-9)*(1e-12/h)**2
    
    for i,(freq,eps_w_i) in enumerate(zip(freqaxis,eps_w)):
        inv_eps_w_i = 1.0/eps_w_i
        
        #Add transition energies to Transition interaction matrix
        B = R.copy()
        B[diag_indices] += eps_w_i*E2_if
            
        #diagonalise
        Bdiag,U = eigen(B)
        
        #final values of R,w0
        Ry2a = np.dot(U.transpose(),d)**2 * const_factor #THz**2 (real)
        wya = np.sqrt(inv_eps_w_i*Bdiag)/h*1e-12 #THz (real)
        
        #calculate the dielectric constant at this frequency
        Xi = susceptibility_Losc(freq,w0=wya,f=Ry2a,w_p=1.0,y0=y_y)
        inveps_b[i]-= np.sum(Xi)*inv_eps_w_i**2
        
    #import ipdb; ipdb.set_trace()
        
    return inveps_b

## frequency dependent dielectric constant by splitting structure into pieces


def inv_eps_zz_multiplasmon3(results,transitions_table,linewidth,freqaxis,dielectric_masks):
    """Uses a multilevel version of the mathematical formalism given in Ando 1977
    A matrix is constucted describing the transitions and the interactions between
    them which can be diagonalised to give a description of the system as a simple
    sequence of Lorentzian oscillators.
    
    This calculates the dielectric constant ratio - 1.0/eps_ISBT for the ISBTs for
    a background dielectric constant given by
    
    dielectric_masks - a sequence of (eps,mask_array) where
        eps - an array or function of dielectric constants wrt the frequency axis. If it is
              a function, it should accept an arguement for frequency in THz.
        mask_array - a bool or integer array wrt the z axis indicating were the eps applies.
    
    linewidth - either a function of the transition frequency or a value (THz)
    freqaxis is an array of frequencies (THz) to calculate the dielectric constant for.
    """
    #check dielectric_mask for completeness
    mask_check = np.zeros_like(results.xaxis,dtype=int)
    mask_check = np.sum((mask for eps,mask in dielectric_masks),mask_check)
    if not all(mask_check == 1): 
        logger.error('masks in dielectric_masks either overlap or do not cover entire structure')
    
    #split the model up in to pieced by material type
    Rs = []
    ds = []
    Epsilons = []
    for eps,mask in dielectric_masks:
        with np.errstate(divide='ignore'):
            maskB = 1.0/mask #infinite where mask==0 and unity where mask==1
        #Calculate transitions interactions matrix + rhs of system equation
        R,d = calc_interaction_matrix(results,transitions_table,eps_z=maskB)
        epsilon = eps(freqaxis) if callable(eps) else eps*np.ones_like(freqaxis)
        Rs.append(R)
        ds.append(d)
        Epsilons.append(epsilon)
    Epsilons = np.column_stack(Epsilons) #array of freqaxis vs structure-pieces-wrt-eps_b
    
    #Calculate the inverse dielectric constant ############
        
    #choose appropriate solver
    if np.iscomplex(Epsilons).any():
        logger.info('calc_wR_multiplasmon3: using eig() solver for complex symmetric or general matrix')
        eigen = lambda B: eig_sorted(B) #matrix will be complex symmetric but not Hermitian, this may be a problem with the theory...
    else:
        logger.info('calc_wR_multiplasmon3: using eigh() solver for Hermitian matrix')
        eigen = lambda B: eigh(B, lower=True, eigvals_only=False, turbo=True, type=1) #otherwise we can be sure that B is real symmetric
    
    #transition energies
    E_if = np.zeros(len(transitions_table))
    for tra in transitions_table:
        E_if[tra['j']] = tra['dE']*meV2J
    E2_if = E_if**2 #transition energies squared
    
    diag_indices = np.diag_indices_from(R) # indices to access the diagonal of the transitions interaction matrix

    #linewidth
    #ff0 = transitions_table[0]['Leff']/transitions_table[0]['Lperiod']
    w_if = np.sort([tra['dE'] for tra in transitions_table])*meV2J/h*1e-12 #(THz) initial transition frequencies
    #w_if = E_if/h*1e-12 #(THz) initial transition frequencies
    y_y = linewidth(w_if) if callable(linewidth) else linewidth*np.ones_like(w_if)
    #y_y = linewidth(w_i) if callable(linewidth) else linewidth #(THz real?) guesstimate of transition broadening (written to get result as close as possible to other models)
    #y_y = linewidth(np.sqrt(wy**2-Ry2/ff0)) if callable(linewidth) else linewidth #(THz real?) guesstimate of transition broadening (written to get result as close as possible to other models)
    
    const_factor = 2.0/(eps0*tra['Lperiod']*1e-9)*(1e-12/h)**2
    
    inveps_b = np.zeros_like(freqaxis) + 0j
    
    for i,freq in enumerate(freqaxis):
        #find dielectric constants for each subsection of the structure
        inv_eps_w_i = 1.0/Epsilons[i]
        #background inverse dielectric constant
        inveps_b[i] = np.mean(np.sum(mask*inve for inve,(_,mask) in zip(inv_eps_w_i,dielectric_masks)))
        #Add the pieces of B together
        B = np.sum(a*r for a,r in zip(inv_eps_w_i,Rs))
        #Add the pieces of d together
        d = np.sum(a*dl for a,dl in zip(inv_eps_w_i,ds))
        
        #Add transition energies to Transition interaction matrix
        B[diag_indices] += E2_if
            
        #diagonalise
        Bdiag,U = eigen(B)
        
        #final values of R,w0
        Ry2a = np.dot(U.transpose(),d)**2 * const_factor #THz**2 (real)
        wya = np.sqrt(Bdiag)/h*1e-12 #THz (real)
        
        #calculate the dielectric constant at this frequency
        Xi = susceptibility_Losc(freq,w0=wya,f=Ry2a,w_p=1.0,y0=y_y)
        inveps_b[i]-= np.sum(Xi)
        
    #import ipdb; ipdb.set_trace()
    
    return inveps_b


## Making plots of absorption

def plotting_absorption(model,results,transitions_table,eps_b,eps_z,linewidth):
    """plots an approximation to the ISBT absorptions of a QW.
    This is really a demo function, add a customised version to your
    script."""    
    f1 = pl.figure()
    ax1 = f1.add_subplot(111)
    ax1.set_xlabel('frequency (THz)')
    transitionfreq = [trn['freq'] for trn in transitions_table]
    frange = max(transitionfreq) - min(transitionfreq)
    freqaxis = np.linspace(0.5*min(transitionfreq),1.5*max(transitionfreq),800)
    
    for trn in transitions_table: #nb. first row of table describes the units of each variable
        #vertical lines for original transition energies
        ax1.axvline(trn['freq'])
    
    theta =pi/4
    nk = np.sqrt(np.mean(np.atleast_1d(eps_z),axis=0)) # should be eps_xx really
    d = transitions_table[0]['Lperiod']*1e-9
    f2w = 1e12*2*pi
    eps_z = np.real_if_close(eps_z)
    
    #model 0 # the slightly niave model usng the 'standard' absorption calculation and Lorentz oscillator model
    # this is only for comparison.
    eps_simple = eps_classical(transitions_table,freqaxis,np.mean(eps_z))#.conjugate()
    Leff0 = get_Leff_est(transitions_table)*1e-9
    absorption_simple = absorption_standard(freqaxis*f2w,eps_simple,Leff0).real
    #eps_b=1.0
    #ff = transitions_table[0]['Leff']/Lperiod
    #absorption_simple = uniaxial_layer_absorption(theta,freqaxis*f2w,eps_b/eps_simple,nk,ff*d)
    ax1.plot(freqaxis,absorption_simple,label='Naive Model')
    
    #model 1 # Uses the analytically correct result for a single transition but can be incorrect for multiple transitions
    eps_ratio1 = eps_b*inv_eps_zz_1(transitions_table,freqaxis,eps_z)
    absorption1 = uniaxial_layer_absorption(theta,freqaxis*f2w,eps_ratio1,nk,d).real
    ax1.plot(freqaxis,absorption1,label='Independent Transitions Model')
    
    #model 2 # A classical approach to modelling multiple transitions. Not exact but accounts for coupling between transitions in a physically intuitive way.
    #eps_ratio2 = eps_b*inv_eps_zz_classical(transitions_table,freqaxis,eps_z)
    #absorption2 = uniaxial_layer_absorption(theta,freqaxis*f2w,eps_ratio2,nk,d)
    #ax1.plot(freqaxis,absorption2,label='Classical Transitions Model')
    
    #model 3 # An accurate model for multiple transitions (neglecting non-parabolicity).  
    wya,Ry2a = calc_wR_multiplasmon(results,transitions_table,eps_z)
    #print 'matrix method results'; print_multiplasmon_transitions(wya,Ry2a)
    inv_eps_zz3 = inv_eps_zz_multiplasmon(wya,Ry2a,transitions_table,linewidth,freqaxis,eps_z)
    eps_ratio3 = eps_b*inv_eps_zz3
    absorption3 = uniaxial_layer_absorption(theta,freqaxis*f2w,eps_ratio3,nk,d).real
    ax1.plot(freqaxis,absorption3,label='Matrix Model')
    
    #model 4 # An accurate model for multiple transitions with frequency dependant dielectric constant
    #the frequency dependence is defined relative to the normal eps_z but needs to be the same for the
    #whole structure.
    eps_w = np.ones_like(freqaxis) #no actual frequency dependence here, this is just a demo.
    inv_eps_zz4 = inv_eps_zz_multiplasmon2(results,transitions_table,linewidth,freqaxis,eps_z,eps_w)
    eps_ratio4 = eps_b*inv_eps_zz4
    absorption4 = uniaxial_layer_absorption(theta,freqaxis*f2w,eps_ratio4,nk,d).real
    ax1.plot(freqaxis,absorption4,label='Matrix Model with eps(w)')

    #model 5 # An accurate model for multiple transitions with frequency dependant dielectric constant
    #The structure is divided up into a few pieces with respect to shared background dielectric constants
    #this enables a more rapid calculation of the relevant matrices at each frequency by appropriately
    #summing the pieces together.
    
    #eps_z -> pieces 
    eps_values = set(eps_z*np.ones(model.n_max))
    if len(eps_values) > 30: logger.warning('plotting_absorption:model5 eps_z has been than 10 pieces, calculation may be slow')
    dielectric_masks = [(eps,(eps_z*np.ones(model.n_max)==eps)) for eps in eps_values]
    
    #normally might calculate the dielectric_masks sequence manually via
    #dielectric_masks = [(eps_w_AlGaAs,np.sum(model.layer_mask(i) for i in [0,1,2,5,6])),
    #                    (eps_w_GaAs,np.sum(model.layer_mask(i) for i in [3,4])),
    #                    ...]
    #where eps_w_AlGaAs could be a number, an array (wrt freq_axis) or a function (that takes freq_axis 
    #as a parameter)
    inv_eps_zz5 = inv_eps_zz_multiplasmon3(results,transitions_table,linewidth,freqaxis,dielectric_masks)
    eps_ratio5 = eps_b*inv_eps_zz5
    absorption5 = uniaxial_layer_absorption(theta,freqaxis*f2w,eps_ratio5,nk,d).real
    ax1.plot(freqaxis,absorption5,label='Matrix Model with eps(z,w)')
        
    ax1.legend()
    if not pl.isinteractive(): pl.show()
    
    #import ipdb; ipdb.set_trace()
    
    return f1
    


def eps_background_GaAs(model,eps_gaas,eps_algaas):
    """Helper function for calculating background dielectric constant
    array for GaAs/AlGaAs structures"""
    eps_z = np.zeros(model.n_max,np.complex)
    
    position = 0.0 # keeping in nanometres (to minimise errors)
    for layer in model.material:
        startindex = round2int(position*1e-9/model.dx)
        position += layer[0] # update position to end of the layer
        finishindex = round2int(position*1e-9/model.dx)
        #
        matType = layer[1]
        if matType == 'GaAs':
            eps_z[startindex:finishindex] = eps_gaas
        elif matType == 'AlGaAs':
            eps_z[startindex:finishindex] = eps_algaas

    return eps_z    



if __name__ == "__main__":
    import config
    import database
    import aestimo
    import os
    import time
    
    np.set_printoptions(precision=3,linewidth=180)
    
    logger = aestimo.logger
    
    # Import from config file
    inputfile = __import__(config.inputfilename)
    logger.info("inputfile is %s",config.inputfilename)
    
    # Initialise structure class
    model = aestimo.StructureFrom(inputfile,database)
    
    if False: #recalculate QW states
        # Perform the calculation
        result = aestimo.Poisson_Schrodinger(model)
        
        time4 = time.time() #timing audit
        logger.info("total running time (inc. loading libraries) %g s",(time4 - aestimo.time0))
        logger.info("total running time (exc. loading libraries) %g s",(time4 - aestimo.time1))
        
        # Write the simulation results in files
        aestimo.save_and_plot(result,model)
        logger.info("Simulation is finished.")
    else: #load previously calculated results from output directory
        result = aestimo.load_results()
    
    #Set thickness of effective medium
    Lperiod = sum([layer[0] for layer in model.material])*1e-9 #m
    
    # set dielectric constants
    case = 2
    if case==1: #scalar dielectric constants
        eps_b = 12.90
        eps_z = 12.90 #+ 0.0j
    
    elif case==2: #z-dependent dielectric constants
        eps_b = 10.364
        eps_gaas = 10.364 # @ 16um
        eps_algaas = 8.2067
        eps_z = eps_background_GaAs(model,eps_gaas,eps_algaas)
        eps_z = np.real_if_close(eps_z)
    
    elif case==3: #w-dependent dielectric constants
        #because the zeroth axis is assumed to be the z-axis, our eps_z array must be 2d
        pass
        #currently the matrix model doesn't cope with frequency dependent dielectric constants
        #therefore the classical model would be the best approach (model2) although it seems
        #to over-estimate the coupling between the transitions.
        #Alternatively, we could resolve the matrix at each frequency (for each value of the
        #background dielectric constant which would be accurate but may be quite computationally
        #intensive.
    
    elif case==4: #z-dependent and w-dependent dielectric constants
        pass
        #currently the matrix model doesn't cope with frequency dependent dielectric constants
        #therefore the classical model is the best approach (model2) although it seems
        #to over-estimate the coupling between the transitions.
        #Alternatively, we could resolve the matrix at each frequency (for each value of the
        #background dielectric constant which would be accurate but may be quite computationally
        #intensive.
    
    # Linewidth
    def linewidth(freq): return 0.1*freq #define linewidth in THz

    #linewidth = 1.0 #THz
    
    # Optical Intersubband Transitions
    transitions_table,(hdr,units) = transitions(result,Lperiod,eps_z,linewidth)
    
    print_levels(result)
    print_transitions(transitions_table,hdr,units)
    
    plotting_absorption(model,result,transitions_table,eps_b,eps_z,linewidth)
    
"""
TO DO:
solve eps_b wrt freq issue
make compatible with pyLuminous/pyFresnel
chang
warburton
"""
