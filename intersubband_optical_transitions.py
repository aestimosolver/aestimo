#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Version v.0.9
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

 Description: This calculates the optical intersubband transitions for the conduction
            band levels.

Theory Notes:
Intersubband transitions (ISBTs) in a quantum well occur between the well's 
different levels but staying within a single band. This is in contrast to
interband transitions between the valence and conduction bands, those transitions
mostly occur in the visible and near-infrared and are used for LEDs, lasers and 
detectors. Intersuband transitions occur in the mid- to far- infrared and are 
used for quantum cascade lasers and mid-infrared detectors etc.

The strongest transitions are electric dipole transtions and so we need to
calculate the dipole matrix elements of the transitions. However, ISBTs have some
additional complications to a standard dipole transition. Firstly, they are 
polarisation sensitive in that they can only couple to light polarised perpendicular
to the quantum well plane. Secondly, there are complications since there are 
many electrons in the quantum well layers; ones such complication is that we can
consider the transtion as a plasma and the effect of the depolarisation field
induced within the plasma when interacting with an electric field leads to a 
shift of the peak's position (called the depolarisation shift). The plasma 
behaviour of the transition is also refered to as the collective excitation of 
the electrons since it the effect of the electrons on each other. In addition 
there is another shift due to the exchange interaction between the electrons. We
may also need to account for the non-parabolicty of the subbands and the different
effective masses of the different layers.

Finally, we may need to consider the free electron absorption from electrons
moving within the plane of the well (for instance using a Drude model). Very rarely
we may need to include weaker transtions such as the magnetic dipole, electric
quadrupole.

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

abs ~ (epsilon_b/epsilon_zz).n_b.(w/c).(sin(theta)**2/cos(theta))*d

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
f01 - oscillator strength = 
w_p01 - plasma frequency = 
y01 - transition broadening = 

Since the effective mass of each level is different in aestimo, it is not clear what
to use when calculating the parameters above that characterise the transition. However,
it is not important for the final transition since all of the effective mass terms
cancel out in the final expression.
"""
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as pl
from itertools import combinations
from scipy.linalg import eigh #,eig

sin,cos,log,exp = np.sin,np.cos,np.log,np.exp

#Defining constants and material parameters
q = 1.602176e-19 #C
kb = 1.3806504e-23 #J/K
nii = 0.0
h = 6.62606957e-34 #m2kg/s
hbar = 1.054588757e-34
m_e = 9.1093826E-31 #kg
pi = np.pi
eps0 = 8.8541878176e-12 #F/m
c = 299792458 #m/s

J2meV=1e3/q #Joules to meV
meV2J=1e-3*q #meV to Joules


# Electromagnetism
# -------------------------

# Standard absorption coefficient for a transition

def absorption_standard(w,epsilon,d):
    """Calculates the absorption for a material from its dielectric 
    'constant'. However, this is not really suitable for modelling
    the absorption of intersubband transitions in quantum wells.    
    """
    return w*np.sqrt(epsilon).imag/c

# Absorption of uniaxial absorbing layer

def uniaxial_layer_absorption(theta,w,eps_ratio,nk,d):
    """Approxiamately calculates the absorption of a uniaxial layer
    which absorbs along its extraordinary axis which is perpendicular
    to the plane of the layer. Basically, this is the situation for 
    an intersubband transition of a quantum well. More accurate modelling
    requires a transfer matrix model that can account for these uniaxial
    layers (i.e. see pyLuminous/pyFresnel module).
    theta - angle (rad)
    w - frequencies array
    eps_ratio = eps_in_plane/eps_out_of_place = eps_xx/eps_zz
    nk - refractive index of media surrounding uniaxial layer
    d - thickness of layer.
    """
    return eps_ratio.imag*nk*w/c*sin(theta)**2/cos(theta)*d

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
    includes eps_b in its definition and so it is necessary to remove it for this
    formula to calculate the susceptibility (and if you want to use a frequency 
    dependent dielectric constant later then it will be more slightly more accurate
    to do this). It's the difference between having to use the result here as 
    eps_b = eps_b0.(1 + Xi).(1 + Xi2)... 
    or (by including eps_b!=1.0 here) as
    eps_b = epsb0 + Xi + Xi2 +...  
    """
    chi = (eps_b*w_p**2)*f / ( w**2 - w0**2 + 2j*y0*w0 )
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


# Singular Electric Dipole Transitions (each transition is treated separately)
# ------------------------------------

def dipole_matrix(z,wfe1,wfe2):
    """Calculates dipole matrix element numerically. Returns values in metres
    (electron charge is not included in calculation)"""
    return simps(z*wfe1*wfe2,z)

def oscStr(mu_if,w_if,cb_meff):
    """Calculates oscillator strength. 
    w_if - frequency of transition (meV),
    cb_meff - relative effective electron mass in well layer, 
    mu_if - dipole matrix element for position operator (m)."""
    return 2*cb_meff*(w_if*meV2J)*mu_if**2/hbar**2

def calc_S(Psi0,Psi1,Psi2,Psi3,dx):
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
    i3*=dx**3
    return -i3

def calc_S_b(Psi0,Psi1,Psi2,Psi3,xaxis):
    """Calculates S, a quantity used to calculate the effective thickness of an
    intersubband transition of a quantum well.
    Psi0 - Psi3 are arrays describing the wavefunctions.
    dz is the step-size (metres) for the arrays (not assumed to be constant).
    """
    dx_axis = xaxis[1:]-xaxis[:-1]
    dx_axis = np.hstack((dx_axis[0],dx_axis)) #preprend a value so that all values get used in calculation
    i1=0.0; i2=0.0; i3=0.0
    for p0,p1,p2,p3,x,dx in zip(Psi0,Psi1,Psi2,Psi3,xaxis,dx_axis):
        delta = p0*p1
        i1+=delta
        i2+=x*delta
        i3+=p2*p3*(x*i1 - i2)*dx**2
    return -i3

def L_eff(w_if,S_if,cb_meff):
    """Calculates the effective thickness of an intersubband transition of a quantum
    well. w_if is the transitions frequency in meV, Mw is the relative effective 
    electron mass in well layer, S is a dimensionless quantity calculated via an
    integral. Returns a value in metres"""
    return hbar**2/(2*abs(S_if)*cb_meff*w_if*meV2J)

def calc_w_p(dn_if,cb_meff,eps_b):
    """Calculates the plasma frequency of a transition.
    dn_ij - (m**-3) population density difference between the initial and final levels
    cb_meff - effective electron mass for transition
    eps_b - background dielectric constant
    """
    return np.sqrt(dn_if*q**2/(cb_meff*eps_b*eps0))/(2*pi) #real Hz
    
def calc_R2(w_if,mu_if,dn_if,eps_b,Lperiod):
    """R**2 = f*Leff*w_p**2/Lperiod (oscillator strength * plasma frequency squared)
    This is more useful for quicky seeing the strength of a transition
    than either the plasma frequency or the oscillator strength alone.
    Also we see that there is need to know the effective mass in the final expression.
    Units of R are Hz (real).
    """
    return dn_if*q**2*2*(w_if*meV2J)*mu_if**2/(hbar**2*(eps_b*eps0)*Lperiod*(2*pi)**2) 

# Summary of Transitions

def transition_generator(seq):
    """returns the possible pairs in the input sequence. Each pair is
    only returned once and the ordering found in the input is maintained"""
    return combinations(seq,2)

def transitions(results,Lperiod):
    """Calculates the parameters needed to describe the intersubband transitions.
    Returns a list of dictionaries (one for each transition) with the following 
    keys:
    'ilevel','flevel','dE','freq','lambda','wavno','dN','z','f','Leff','S_if','S_if_b','wp'
    """
    E_state = results.E_state #list of energy levels (meV)
    N_state = results.N_state #occupation of energy levels (m**-2)
    T = results.T #K
    dx = results.dx #m
    xaxis = results.xaxis #m
    wfe = results.wfe*dx**-0.5 # normalising wavefunctions to be m**-0.5
    #reversethepolarities = np.ones(wfe.shape[0])
    #reversethepolarities[1::2]*=-1
    #for j,p in enumerate(reversethepolarities):wfe[j]*=p
    
    print 'the energy levels\population are (meV)\t(m**-2):'
    for Ei,Ni in zip(E_state,N_state): print Ei,'\t',Ni
    print
    print 'T = %gK' %T
    print
    #print 'the energy levels gaps are'
    #print '\t'.join(('(meV)','(THz)','(um)','(wavno)'))
    #for leveli,levelj in transition_generator(levels):
    #    gaps=levelj-leveli
    #    freq=gap*1e-3*eC/h/1e12
    #    wav=1e6*h*c/(gap*1e-3*eC)
    #    wavno=gap*1e-3*eC/h/c*1e-2
    #    print '\t'.join((gap,freq,wav,wavno))
    #print
    #Summary of intersubband transitions
    hdr=['j','ilevel','flevel','dE','freq','lambda','wavno','dN','z','f','Leff','S_if','S_if_b','wp','R','Lperiod','y_if']
    units=['','','','meV','THz','um','cm-1','1e11cm-2 @%gK'%T,'nm','','nm','nm','nm','THz','THz','nm','THz']
    table=[hdr,units]
    
    def transition(i,f): #Doing it this way would let me create a dielectric function for each transition using a function closure.
        meV=E_state[f]-E_state[i]
        THz=meV*1e-3*q/h/1e12 #THz (real)
        um=1e6*h*c/(meV*1e-3*q)
        wavno=meV*1e-3*q/h/c*1e-2
        dN=N_state[i]-N_state[f] #m-2
        mu=dipole_matrix(xaxis,wfe[i],wfe[f]) #m
        #
        cb_meff= 0.067*m_e
        fif=oscStr(mu,meV,cb_meff)
        #
        Sif=calc_S(wfe[i],wfe[f],wfe[i],wfe[f],dx)
        Sif_b=calc_S_b(wfe[i],wfe[f],wfe[i],wfe[f],xaxis)
        Lif=L_eff(meV,Sif,cb_meff)
        #
        eps_b = 12.90
        wp=calc_w_p(dN/Lif,cb_meff,eps_b) #real Hz
        #
        Rif=np.sqrt(calc_R2(w_if,mu_if,dN,eps_b,Lperiod)) #real Hz
        #
        yif = THz*0.1 #(THz real) guesstimate of transition broadening
        #
        col=i,f,meV,THz,um,wavno,dN*10**(-4-11),mu*1e9,fif,Lif*1e9,Sif*1e9,Sif_b*1e9,wp*1e-12,Rif*1e-12,Lperiod*1e9,yif
        return col
    
    for j,(i,f) in enumerate(transition_generator(np.arange(len(E_state)))):
        col=transition(i,f)
        table.append((j,)+col)
    
    print "Summary of Intersubband Transitions"
    wids=[8,14]+[11]*(len(table[0])-3)
    for row in zip(*table):
        print ''.join([row[0].rjust(wids[0]),row[1].rjust(wids[1])]+[('%.3g' %item).rjust(wid) for wid,item in zip(wids[2:],row[2:])])
    
    transitions_table = [dict((key,value) for key,value in zip(hdr,row)) for row in table[1:]]
    #nb. the first row of the transitions table consists of the units for each parameter
    return transitions_table

## Calculating Dielectric Constants
# Below we have some different models for intersubband transitions

def eps_zz_ratio(transitions_table,freqaxis):
    """calculates eps_b/eps_zz using the analytical result for a single transition.
    If there are several active transitions that are close together then this will
    become increasingly incorrect. The dielectric constant is calculated for the 
    effective medium of QW + barrier"""
    inveps = 1.0
    for trn in transitions_table[1:]: #nb. first row of table describes the units of each variable
        Xi = susceptibility_Losc(freqaxis,trn['freq'],trn['R']**2,1.0,trn['y_if'])
        inveps-= Xi
    return inveps

def eps_classical(transitions_table,freqaxis,eps_b=1.0):
    """calculates total dielectric constant epszz for QW by summing susceptibilities
    of each transition. If you leave eps_b=1.0 then the result should be multiplied
    by eps_b, if you use eps_b then it should match the values used for calculating
    the transition plasma frequencies.
    
    WARNING: This shouldn't be used on its own for modelling ISBTs."""
    epszz = eps_b
    for trn in transitions_table[1:]: #nb. first row of table describes the units of each variable
        Xi = susceptibility_Losc(freqaxis,trn['freq'],trn['f'],trn['wp'],trn['y_if'],eps_b=ep_b)
        epszz+=Xi
    return epszz

def eps_zz_ratio_classical(transitions_table,freqaxis):
    """Calculates the dielectric constant for an effective medium containing 
    intersubband transitions using a classical approach (Lorentz oscillators &
    formula for effective medium)
    """
    eps_classical()
    eff_eps_z(layer_list,isbt_term=0.0)

## Advanced Multilevel Model #################
    
def inv_eps_zz_Ando(results,transitions_table,freqaxis,eps_b):
    """Uses a multilevel version of the mathematical formalism given in Ando 197? ? ??
    A matrix is constucted describing the transitions and the interactions between
    them which can be diagonalised to give a description of the system as a simple
    sequence of Lorentzian oscillators.
    """
    #number of transitions
    ntr = len(transitions_table) #or n*(n-1)/2 where n is the number of energy levels   
    # Note that following arrays are indexed by transition rather than energy level.
    #So [0,0] is the transition for levels 0->1
    def trindex(i,j,ntr):
        """Takes energy level numbers for a transition and returns the chosen B index"""
        return j- i + ntr*i - i*(i+1)/2
    #construct transition matrix S & B
    S = np.zeros((ntr,ntr))
    for tra in transitions_table:
        a = trindex(tra['ilevel'],tra['flevel'])
        S[a,a] = tra['S_if']
    wfe = results.wfe
    xaxis = resuls.xaxis
    for tra,trb in combinations(transitions_table):
        a = trindex(tra['ilevel'],tra['flevel'])
        b = trindex(trb['ilevel'],trb['flevel'])
        S[a,b] = S[b,a] = calc_S_b(wfe[tra['ilevel']],wfe[tra['flevel']],wfe[tra['ilevel']],wfe[tra['flevel']],xaxis)
    #B = 2*S_ab*e**2/(eps*eps0) * sqrt((n_a - n_-a)*(n_b - n_b)*hbar*w_a*hbar*w_b) + delta_ab*hbar**2*w_a**2
    const = 2*q**2/(eps*eps0)*meV2J*1e15 # 1e15 converts dN values into carriers/m**2
    B = np.zeros((ntr,ntr))
    for tra,trb in permutations(transitions_table):
        a = trindex(tra['i'],tra['j'])
        b = trindex(trb['i'],trb['j'])
        B[a,b] = const*S[a,b]*np.sqrt(tra['dN']*trb['dN']*tra['dE']*tra['dE'])
        if a==b: B[a,a] += (tra['dE']*meV2J)**2
    #diagonalise
    Bdiag,U = eigh(B, lower=True, eigvals_only=False, turbo=True, type=1)
    #final values of R,w0
    rhs = np.zeros(ntr)
    for tra in transitions_table:
        a = trindex(tra['i'],tra['j']) #find correct index
        rhs[a] = np.sqrt(tra['dN']*tra['dE']/hbar)*tra['z']
    Ry = U*rhs
    wy = diag(Bdiag)    
    #calculate inveps
    
    return (wy,Ry),inveps

## pyFresnel
#class borrowed from pyFresnel so we can integrate better with that library.
class LayerUniaxial_eps():
    def __init__(self,epsxx,epszz,d):
        self.epsxx = epsxx # can be a function or array (if it has the same length as the spectral axis)
        self.epszz = epszz # can be a function or array (if it has the same length as the spectral axis)
        self.d = d
    
    def n(self,w):
        try: #Is self.eps a function of w?
            n = N.sqrt(self.epsxx(w))
        except:
            #Is self.eps an array with the same length as w? 
            if hasattr(self.epsxx,'__len__'):
                if len(self.epsxx)!=len(w):
                    raise NameError("w and epsilon arrays are not compatible")
                else:
                    n = N.sqrt(self.epsxx)
            else: #Is self.eps a number (integer/float/complex)
                n = N.repeat(N.sqrt(self.epsxx),len(w))

        return n

    def nzz(self,w):
        try: #Is self.eps a function of w?
            n = N.sqrt(self.epszz(w))
        except:
            #Is self.eps an array with the same length as w? 
            if hasattr(self.epszz,'__len__'):
                if len(self.epszz)!=len(w):
                    raise NameError("w and epsilon arrays are not compatible")
                else:
                    n = N.sqrt(self.epszz)
            else: #Is self.eps a number (integer/float/complex)
                n = N.repeat(N.sqrt(self.epszz),len(w))

        return n
                
    def __repr__(self):
        return "Layer"+"("+repr(self.epsxx)+", "+repr(self.epszz)+", "+repr(self.d)+" )"



## Making plots of absorption

def plotting_absorption(model,results,transitions_table):
    """plots an approximation to the ISBT absorptions of a QW"""    
    f1 = pl.figure()
    ax1 = f1.add_subplot(111)
    ax1.set_xlabel('frequency (THz)')
    transitionfreq = [trn['freq'] for trn in transitions_table[1:]]
    frange = max(transitionfreq) - min(transitionfreq)
    freqaxis = np.linspace(0.5*min(transitionfreq),1.5*max(transitionfreq),400)
    
    for trn in transitions_table[1:]: #nb. first row of table describes the units of each variable
        #vertical lines for original transition energies
        ax1.axvline(trn['freq'])
    
    theta =pi/4
    nk = np.sqrt(12.90)
    d = 20e-9
    #model 1
    eps_ratio1 = inv_eps_zz(transitions_table,freqaxis)
    absorption1 = uniaxial_layer_absorption(theta,freqaxis,eps_ratio1,nk,d)
    ax1.plot(freqaxis,absorption1,label='model 1')
    
    #model 2
    eps_ratio2 = 1.0/eps_zz(transitions_table,freqaxis,eps_b=1.0)
    absorption2 = uniaxial_layer_absorption(theta,freqaxis,eps_ratio2,nk,d)
    ax1.plot(freqaxis,absorption2,label='model2')
    
    ax1.legend()
    pl.show()

if __name__ == "__main__":
    import config
    import database
    import aestimo_numpy2 as aestimo
    import os
    import time
    
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
    
    # Optical Intersubband Transitions
    xaxis = result.xaxis
    dx = result.dx
    wfe = result.wfe*dx**-0.5 # normalising wavefunctions to be m**-0.5
    mu_if = dipole_matrix(xaxis,wfe[0],wfe[1])
    print mu_if
    print result.dx
    print simps(wfe[1]*wfe[1],xaxis)/result.dx
    
    w_if = result.E_state[1] - result.E_state[0]
    cb_meff = result.meff_state[0]
    print cb_meff
    print w_if
    print mu_if
    print oscStr(mu_if,w_if,cb_meff)
    
    print [i for i in transition_generator(range(4))]
        
    transitions_table = transitions(result,Lperiod)
    plotting_absorption(model,result,transitions_table)
    
"""
TO DO:
make compatible with pyLuminous/pyFresnel
chang
warburton
"""