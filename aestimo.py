#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
  Program:    Aestimo 1D Schrodinger-Poisson Solver
  Version:    v.0.6

  Description: This is the aestimo calculator.
  
"""
#from scipy.optimize import fsolve
from pylab import *
import numpy as np
import sys,config,database
from math import *

# Import from config file
inputfile = __import__(config.inputfilename)
# Flags for different outputs
E_out = config.electricfield_out
V_out = config.potential_out
S_out = config.sigma_out
P_out = config.probability_out
St_out = config.states_out
Resultview = config.resultviewer

print "Aestimo is starting..."
# Reading inputs and using local variables
max_val = inputfile.maxgridpoints
T = inputfile.T
comp_scheme = inputfile.computation_scheme
dx = inputfile.gridfactor*1e-9 #grid in m
x_max = (inputfile.z_coordinate_end - inputfile.z_coordinate_begin)*1e-9 # in m

#Making a seperate material list
material = []
for item in inputfile.material:
    material.append(item)
totallayer = alen(material)

material_property = []
for item in database.materialproperty:
    material_property.append(item)
totalmaterial = alen(material_property)

alloy_property = []
for item in database.alloyproperty:
    alloy_property.append(item)
totalalloy = alen(alloy_property)

# Changing material position info to meter
for j in range(0, totallayer,1):
    material[j][1] = float(inputfile.material[j][1])*1e-9
    material[j][2] = float(inputfile.material[j][2])*1e-9
print "Total layer number: ",totallayer
print "Total material number in database: ",totalmaterial

# Preparing empty subband energy lists. We will not use E_state[0]
E_state = []
N_state = []
E_stateresult = []
wfetot = []
E_state = [0.0]*(inputfile.subnumber_e+1)
N_state = [0.0]*(inputfile.subnumber_e+1)
N_stateresult = [0.0]*(inputfile.subnumber_e+1)
wfetot = [0.0]*(inputfile.subnumber_e+1)

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

manual_iterate = 3
previousE0=0
#subband_n_ratio = [0.8,0.2,0.0]
#subband_n_ratio = [0.8, 0.2] # This must be automized


# DO NOT EDIT UNDER HERE FOR PARAMETERS
# --------------------------------------

#Vegard's law for alloys
def vegard(first,second,mole):
    return first*mole+second*(1-mole)

# This function returns the value of the wavefunction (psi)
# at +infinity for a given value of the energy.  The solution
# to the energy occurs for psi(+infinity)=0.

# FUNCTIONS for SHOOTING ------------------
def psi_at_inf(E,fis,cb_meff):
    psi = []
    psi = [0.0]*3
    # boundary conditions
    psi[0] = 0.0                 
    psi[1] = 1.0
    if T_flag == True:
        for j in range(0,n_max,1):
            # Last potential not used
            c1=2.0/(cb_meff[j]+cb_meff[j-1])
            c2=2.0/(cb_meff[j]+cb_meff[j+1])
            psi[2]=((2*(dx/hbar)**2*(fis[j]-E)+c2+c1)*psi[1]-c1*psi[0])/c2
            psi[0]=psi[1]
            psi[1]=psi[2]
    else:
        for j in range(0,n_max,1): # Last potential not used
            psi[2]=(2*cb_meff[j]*(fis[j]-E)*(dx/hbar)**2+2)*psi[1]-psi[0]
            psi[0]=psi[1]
            psi[1]=psi[2]
    return psi[2]

#nb. function was much slower when fi is a numpy array than a python list.
def calc_E_state(numlevels,fi,cb_meff,energyx0): # delta_E,d_E
    energyx=energyx0
    E_state=[0.0]*(numlevels)
    for i in range(0,numlevels,1):  
        #increment energy-search for f(x)=0
        y2=psi_at_inf(energyx,fi,cb_meff)
        while True:
            y1=y2
            energyx += delta_E
            y2=psi_at_inf(energyx,fi,cb_meff)
            if y1*y2 < 0:
                break
        # improve estimate using midpoint rule
        energyx -= abs(y2)/(abs(y1)+abs(y2))*delta_E
        #implement Newton-Raphson method
        while True:
            y = psi_at_inf(energyx,fi,cb_meff)
            dy = (psi_at_inf(energyx+d_E,fi,cb_meff)- psi_at_inf(energyx-d_E,fi,cb_meff))/(2.0*d_E)
            energyx -= y/dy
            if abs(y/dy) < 1e-12*q:
                break
        E_state[i]=energyx/(1e-3*q)
        # clears x from solution
        if not(config.messagesoff) :
            print "E[",i,"]=",E_state[i],"meV" #can be written on file.
        energyx += delta_E # finish for i-th state.
    return E_state

# FUNCTIONS for ENVELOPE FUNCTION WAVEFUNCTION--------------------------------
def wf(E,fis,cb_meff):
    # This function returns the value of the wavefunction (psi)
    # at +infinity for a given value of the energy.  The solution
    #	to the energy occurs for psi(+infinity)=0.
    # psi[3] wavefunction at z-delta_z, z and z+delta_z 
    # i index
    N = 0.0 # Normalization integral
    psi = []
    psi = [0.0]*3
    # boundary conditions
    psi[0] = 0.0                 
    psi[1] = 1.0
    b[0] = psi[0]
    b[1] = psi[1]
    N += (psi[0])**2
    N += (psi[1])**2
    if T_flag == True:
        for j in range(0,n_max-1,1):
            # Last potential not used
            c1=2.0/(cb_meff[j]+cb_meff[j-1])
            c2=2.0/(cb_meff[j]+cb_meff[j+1])
            psi[2] = ((2*(dx/hbar)**2*(fis[j]-E)+c2+c1)*psi[1]-c1*psi[0])/c2
            b[j+1]=psi[2]
            N += (psi[2])**2
            psi[0]=psi[1]
            psi[1]=psi[2]
    else:
        for j in range(0,n_max,1): # Last potential not used
            psi[2] = (2*cb_meff[j]*(fis[j]-E)*(dx/hbar)**2+2)*psi[1]-psi[0]
            b[j+1]=psi[2]
            N += (psi[2])**2
            psi[0]=psi[1]
            psi[1]=psi[2]
    return b,float(N*dx)
    
# FUNCTIONS for FERMI-DIRAC STATISTICS-----------------------------------------   
def fd2(Ei,Ef,T):
    #integral of Fermi Dirac Equation for energy independent density of states.
    #Ei [meV], Ef [meV], T [K]"""
    return kb*T*log(exp(meV2J*(Ef-Ei)/(kb*T))+1)

def calc_meff_state(wfe,cb_meff):
    #find subband effective mass
    meff_state = [0.0]*alen(wfe)
    for j in range(0,inputfile.subnumber_e,1):
        total=0.0
        for b,meff in zip(wfe[j],cb_meff):
            total+=float(b)**2/meff
        meff_state[j] = 1.0/total
    return meff_state
    
def fermilevel_0K(Ntotal2d,E_state,meff_state):
    Et,Ef=0.0,0.0
    for i,(Ei,csb_meff) in enumerate(zip(E_state,meff_state)):
        Et+=Ei
        Efnew=(-Ntotal2d*hbar**2*pi/csb_meff*J2meV + Et)/(i+1)
        if Efnew>Ei:
            Ef=Efnew
        else:
            break #we have found Ef and so we should break out of the loop
    else: #exception clause for 'for' loop.
        print "Have processed all energy levels present and so can't be sure that Ef is below next higher energy level."
    
    N_state=[0.0]*len(E_state)
    for i,(Ei,csb_meff) in enumerate(zip(E_state,meff_state)):
        Ni=(Ef - Ei)*csb_meff/(hbar**2*pi)*meV2J    # populations of levels
        Ni*=(Ni>0.0)
        N_state[i]=Ni
    return Ef,N_state
    
def fermilevel(Ntotal2d,T,E_state,meff_state):
    #find the Fermi level
    def func(Ef,E_state,meff_state,Ntotal2d,T):
        #return Ntotal2d - sum( [csb_meff*fd2(Ei,Ef,T) for Ei,csb_meff in zip(E_state,meff_state)] )/(hbar**2*pi)
        diff = -Ntotal2d
        for Ei,csb_meff in zip(E_state,meff_state):
            diff -= csb_meff*fd2(Ei,Ef,T)/(hbar**2*pi)
        return diff
    Ef_0K,N_states_0K = fermilevel_0K(Ntotal2d,E_state,meff_state)
    #Ef=fsolve(func,Ef_0K,args=(E_state,meff_state,Ntotal2d,T))[0]
    #return float(Ef)
    #implement Newton-Raphson method
    Ef = Ef_0K
    d_E = 1e-9
    while True:
        y = func(Ef,E_state,meff_state,Ntotal2d,T)
        dy = (func(Ef+d_E,E_state,meff_state,Ntotal2d,T)- func(Ef-d_E,E_state,meff_state,Ntotal2d,T))/(2.0*d_E)
        Ef -= y/dy
        if abs(y/dy) < 1e-12:
            break
    return Ef

def calc_N_state(Ef,T,Ns,E_state,meff_state):
    # Find the subband populations, taking advantage of step like d.o.s. and analytic integral of FD
    N_state=[fd2(Ei,Ef,T)*csb_meff/(hbar**2*pi) for Ei,csb_meff in zip(E_state,meff_state)]
    return N_state
    
# FUNCTIONS for SELF-CONSISTENT POISSON--------------------------------
def calc_field(sigma,eps):
    # F electric field as a function of z-
    # i index over z co-ordinates
    # j index over z' co-ordinates

    # For wave function initialise F
    for i in range(0,n_max,1):
        F[i] = 0.0
    for i in range(0,n_max,1):
        for j in range(0,n_max,1):
            # Note sigma is a number density per unit area, needs to be converted to Couloumb per unit area
            F[i] = F[i] + q*sigma[j]*cmp(i,j)/(2*eps[j]) #CMP'deki i ve j yer değişebilir - de + olabilir
    return F

def calc_potn(F):
    # This function calculates the potential (energy actually)
    # V electric field as a function of z-
    # i	index over z co-ordinates

    #Calculate the potential, defining the first point as zero
    V = [0.0] * n_max
    for i in range(0,n_max,1):
        V[i]=V[i-1]+q*F[i]*dx #+q -> electron -q->hole? 
    return V

def calc_sigma(wfe,N_state,dop):
    # This function calculates `net' areal charge density
    # i index over z co-ordinates
    # is index over states
    for i in range(0,n_max,1):
        sigma[i] = 0.0
    for i in range(0,n_max,1):
        for j in range(0,inputfile.subnumber_e,1):
            sigma[i] = sigma[i] - N_state[j]*(float(wfe[j][i])**2)
            # n-type dopants give -ve *(N+j) representing electrons, hence 
            # addition of +ve ionised donors requires -*(Nda+i), note Nda is still a
            # volume density, the delta_z converts it to an areal density
        sigma[i] = sigma[i] - dop[i]*dx # This may be one tab indented.

    return sigma
    
def dop0(dop):
    posi = 0.0
    for i in range(0, n_max, 1):
        posi = i*dx
        for j in range(0, totallayer,1):
            if posi >= material[j][1] and posi <= material[j][2]:
                k=j
        if material[k][6] == 'n':
            if material[k][5] == 0:
                dop[i] = nii
            else:
                dop[i] = -material[k][5]*1e6
        else:
            if material[k][5] == 0:
                dop[i] = nii
            else:
                dop[i] = material[k][5]*1e6
    return dop
# ----------------------------------------------------

# Calculate the required number of grid points and renormalize dx
n_max = int(x_max/dx)
if n_max > max_val:
    print " Grid number is exceeding the max number of ", max_val
    exit()

# Creating and Filling material arrays
cb_meff = []
Nc = []
Nv = []
ni = []
fi = []
fis = []
fitot = []
eps =[]
cb_meff = [0.0]*(n_max+1)
Nc = [0.0]*(n_max+1)
Nv = [0.0]*(n_max+1)
ni = [0.0]*(n_max+1)
fi = [0.0]*(n_max+1)
fis = [0.0]* (n_max+1)
fitot = [0.0]*(n_max+1)
eps = [0.0]*(n_max+1)
posi = 0.0
fi_min= 0.0

dop = []
dop = [0.0]*(n_max+1)

Ntotal = 0.0
Ntotal2d = 0.0
sigma = []
sigma = [0.0] * (n_max+1)

F = []
F = [0.0] * (n_max+1)

b = []
b = [0.0]*(n_max+1)
# Subband wavefunction for electron list. 2-dimensional: [i][j] i:stateno, j:wavefunc
wfe = np.zeros((inputfile.subnumber_e,n_max),dtype = float)

for i in range(0, n_max+1, 1):
    posi = i*dx
    for j in range(0, totallayer,1):
        for m in range (0, totalmaterial,1):
            if posi >= material[j][1] and posi <= material[j][2]:
                k=j
            if material[k][3] == material_property[m][1]:
                cb_meff[i] = material_property[m][2]*m_e
                fi[i] = material_property[m][6]*material_property[m][5]*q #Joule
                eps[i] = material_property[m][4]*eps0
        for m in range (0, totalalloy,1):
            if posi >= material[j][1] and posi <= material[j][2]:
                k=j
            if material[k][3] == alloy_property[m][1]:
                cb_meff[i] = (alloy_property[m][2]+alloy_property[m][3]*material[k][4])*m_e
                fi[i] = alloy_property[m][6]*material[k][4]*q*alloy_property[m][7] # for electron. Joule
                eps[i] = (alloy_property[m][4]+alloy_property[m][5]*material[k][4])*eps0
    
    # Find fi-minimum
    if fi[i] < fi_min:
        fi_min= fi[i]

# Setup the doping

dop = dop0(dop)
Ntotal = sum(dop) # calculating total doping density m-3
Ntotal2d = Ntotal*dx
#print "Ntotal ",Ntotal,"m**-3"
print "Ntotal2d ",Ntotal2d," m**-2"

delta_acc = 1e-6

# Shooting method for Schrödinger Equation solution
delta_E = 1e-3*q #Initial delta_E is 1 meV. This can be included in config as a setting?
E_start = 0.0 #This can be included in config as a setting?
T_flag = True
d_E = 1e-8*q

if abs(E_start)<1e-6*q:
    energyx = fi_min
else:
    energyx = E_start

# Populize subbands. Must be rewritten, it is manual now!
#for j in range(0,inputfile.subnumber_e,1):
#    N_state[j] = subband_n_ratio[j]*Ntotal

# STARTING SELF CONSISTENT LOOP
iteration = 0
fitot = list(fi) #For initial iteration just copy fi. list(seq) returns a copy of the original rather than just an alias.

while True:
    if not(config.messagesoff) :
        print "Iteration:",iteration+1
    if iteration> 0:
        for i in range(0, n_max, 1):
            # Find fi-minimum --may got error.
            if fitot[i] < fi_min:
                energyx = fitot[i]
            else:
                energyx = fi_min
    
    E_state=calc_E_state(inputfile.subnumber_e,fitot,cb_meff,energyx)
    
    # Envelope Function Wave Functions
    for j in range(0,inputfile.subnumber_e,1):
        if not(config.messagesoff) :
            print "Working for subband no:",j+1
        b,Ntrial = wf(E_state[j]*1e-3*q,fitot,cb_meff)
        for i in range(0,n_max,1):
            wfe[j][i]=b[i]/(Ntrial/dx)**0.5 #Ntrial/dx?
    
    # Calculate the effective mass of each subband
    meff_state = calc_meff_state(wfe,cb_meff)
    
    for i,meff in enumerate(meff_state):
        print 'meff[',i,']= ',meff/m_e

    ## Self-consistent Poisson
    
    # Calculate the Fermi energy and subband populations at 0K
    E_F_0K,N_state_0K=fermilevel_0K(Ntotal2d,E_state,meff_state)
    print 'Efermi (at 0K) = ',E_F_0K,' J'
    #for i,Ni in enumerate(N_state_0K):
    #    print 'N[',i,']= ',Ni
    
    # Calculate the Fermi energy at the temperature T (K)
    E_F = fermilevel(Ntotal2d,T,E_state,meff_state)
    # Calculate the subband populations at the temperature T (K)
    N_state=calc_N_state(E_F,T,Ntotal2d,E_state,meff_state)
    
    for i,Ni in enumerate(N_state):
        print 'N[',i,']= ',Ni
    
    # Calculate `net' areal charge density and output to file
    sigma=calc_sigma(wfe,N_state,dop) #one more instead of inputfile.subnumber_e
    print "total donor charge = ",sum(dop)*dx,"m**-2"
    print "total system charge = ",sum(sigma),"m**-2"
    # Calculate electric field and output to file
    F=calc_field(sigma,eps)
    # Calculate potential due to charge distribution and output to file	*/
    V=calc_potn(F)   
    # Combine band edge potential with potential due to charge distribution */
    for i in range(0,n_max,1):
        fitot[i] = fi[i] + V[i]
    if abs(E_state[0]-previousE0) < 1e-6:
        break
    else:
        iteration += 1
        previousE0 = E_state[0]

# END OF SELF-CONSISTENT LOOP
# Write the simulation results in files
if S_out:
    out_sigma_file = open("outputs/sigma.dat", "w")
if E_out:
    out_efield_file = open("outputs/efield.dat", "w")
if V_out:
    out_potn_file = open("outputs/potn.dat", "w")
if St_out:
    out_states_file = open("outputs/states.dat", "w")
if P_out:
    out_first_file = open("outputs/firststate.dat", "w")

xx = 0.0
# Arrays for drawing, maybe not necessary. Will be deleted in future.
xaxis = []
drw_sigma = []
drw_efield = []
drw_potn = []
drw_state1 = []
drw_state2 = []
drw_first = []

# Filling them with zeros
xaxis = [0.0]*(n_max)
drw_sigma = [0.0]*(n_max)
drw_efield = [0.0]*(n_max)
drw_potn = [0.0]*(n_max)
drw_potn_meV = [0.0]*(n_max)
drw_state = [0.0]*(n_max)
drw_states = [list(drw_state) for i in range(inputfile.subnumber_e)] # need to create copies of original list rather than references
drw_first = [0.0]*(n_max)

for i in range(0,n_max,1):
    # Sigma
    if S_out:
        out_sigma_file.write(repr(xx) + " " + repr(sigma[i]) + "\n")
    # Electric Field
    if E_out:
        out_efield_file.write(repr(xx) + " " + repr(F[i]) + "\n")
    # Potential
    if V_out:
        out_potn_file.write(repr(xx) + " " + repr(fitot[i]) + "\n")
    # Probability of First state
    if P_out:
        out_first_file.write(repr(xx) + " " + repr(wfe[0][i]) + "\n")
    # for plotting
    xaxis[i] = xx
    drw_sigma[i] = sigma[i]
    drw_efield[i] = F[i]
    drw_potn[i] = fitot[i]
    drw_potn_meV[i] = fitot[i]*J2meV
    for state,wf in zip(drw_states,wfe):
        state[i] = wf[i]
    xx = xx+dx

if St_out:
    for j in range(0,inputfile.subnumber_e,1):
        out_states_file.write(repr(j) + " " + repr(N_state[j]) + " " + repr(E_state[j]) +"\n")

# Closing EQ files
if S_out:
    out_sigma_file.close() 
if E_out:
    out_efield_file.close()
if V_out:
    out_potn_file.close()
if St_out:
    out_states_file.close()
if P_out:
    out_first_file.close()

# Resultviewer
if Resultview:
    figure(figsize=(10,8))
    suptitle('Aestimo Results')
    subplots_adjust(hspace=0.4,wspace=0.4)
                          
    #Plotting Sigma
    #figure(0)
    subplot(2,2,1)
    plot(xaxis, drw_sigma)
    xlabel('Position (m)')
    ylabel('Sigma (e/m^2)')
    title('Sigma')
    grid(True)

    #Plotting Efield
    #figure(1)
    subplot(2,2,2)
    plot(xaxis, drw_efield)
    xlabel('Position (m)')
    ylabel('Electric Field strength (V/m)')
    title('Electric Field')
    grid(True)

    #Plotting Potential
    #figure(2)
    subplot(2,2,3)
    plot(xaxis, drw_potn)
    xlabel('Position (m)')
    ylabel('[V_cb + V_p] (J)')
    title('Potential')
    grid(True)

    #Plotting State(s)
    #figure(3)
    subplot(2,2,4)
    for j,drw_state in enumerate(drw_states):
        plot(xaxis, drw_state, label='state %d' %j)
    xlabel('Position (m)')
    ylabel('Psi')
    title('First state')
    grid(True)
    
    #QW representation
    #figure(5)
    figure(figsize=(10,8))
    suptitle('Aestimo Results')
    subplot(1,1,1)
    plot(xaxis,drw_potn_meV,'k')
    for state,drw_state in zip(E_state,drw_states): 
        axhline(state,0.1,0.9,color='g',ls='--')
        plot(xaxis, np.array(drw_state)*200.0+state,'b')
    xlabel('Position (m)')
    ylabel('Energy (meV)')
    grid(True)
    show()

print "Simulation is finished. All files are closed."
print "Please control the related files."
