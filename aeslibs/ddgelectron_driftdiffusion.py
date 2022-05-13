# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 13:59:47 2019



## Copyright (C) 2004-2008  Carlo de Falco
##
## SECS1D - A 1-D Drift--Diffusion Semiconductor Device Simulator
##
##  SECS1D is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  SECS1D is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with SECS1D; If not, see <http://www.gnu.org/licenses/>.
##
## author: Carlo de Falco <cdf _AT_ users.sourceforge.net>

## -*- texinfo -*-
##
## @deftypefn {Function File}@
## {@var{n}} = DDGelectron_driftdiffusion(@var{psi},@var{xaxis},@var{ng},@var{p},@var{ni},@var{TAUN0},@var{TAUP0},@var{mun})
##
## Solve the continuity equation for electrons
##
## Input:
## @itemize @minus
## @item psi: electric potential
## @item xaxis: integration domain
## @item ng: initial guess and BCs for electron density
## @item p: hole density (for SRH recombination)
## @end itemize
##
## Output:
## @itemize @minus
## @item n: updated electron density
## @end itemize
##
## @end deftypefn
"""

import numpy as np
from math import*
from scipy import sparse as sp

from .func_lib import DDGphin2n,DDGphip2p,Ucompmass,Ucomplap,Ucompconst,Ubernoulli
from .aestimo_poisson1d import equi_np_fi222
import config
    
def DDGelectron_driftdiffusion(psi,xaxis,ng,p,ni,TAUN0,TAUP0,mun,fi_e,fi_h,model,Vt,idata):
    
    nodes        = xaxis
    n_max     =len(nodes)
    """
    n=np.zeros(n_max)
    p=np.zeros(n_max)
    """
    fi_n=np.zeros(n_max)
    fi_p=np.zeros(n_max)    
    elements=np.zeros((n_max-1,2))
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    Nelements=np.size(elements[:,0])
    
    BCnodes= [0,n_max-1]
    
    nl = ng[0]
    nr = ng[n_max-1]
    h=nodes[1:n_max]-nodes[0:n_max-1]
    
    c=1/h
    """
    print("c=",c)
    print("h=",h)
    print("nr=",nr)
    print("nl=",nl)
    print("BCnodes=",BCnodes)
    print("Nelements=",Nelements)
    print("elements=",elements)
    print("n_max=",n_max)
    print("nodes=",nodes)
    check_point_15
    """
    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,psi,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,idata.n,p)    
    Bneg=Ubernoulli(-(psi[1:n_max]-psi[0:n_max-1])-(fi_n[1:n_max]-fi_n[0:n_max-1]),1)    
    Bpos=Ubernoulli( (psi[1:n_max]-psi[0:n_max-1])+(fi_p[1:n_max]-fi_p[0:n_max-1]),1)
    """
    print("Bneg=",Bneg)
    print("Bpos=",Bpos)
    check_point_16 
    """
    
    d0=np.zeros(n_max)
    d0[0]=c[0]*Bneg[0]
    d0[n_max-1]=c[len(c)-1]*Bpos[len(Bpos)-1]
    d0[1:n_max-1]=c[0:len(c)-1]*Bpos[0:len(Bpos)-1]+c[1:len(c)]*Bneg[1:len(Bneg)]    
    
    d1	= np.zeros(n_max)
    d1[0]=n_max
    d1[1:n_max]=-c* Bpos      
    dm1	= np.zeros(n_max)
    dm1[n_max-1]=n_max
    dm1[0:n_max-1]=-c* Bneg 
    """
    print("d0=",d0)
    print("d1=",d1)
    print("dm1=",dm1)
    check_point_17
    """
    A = sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense() 
    
    
    b = np.zeros(n_max)#%- A * ng
 
    ## SRH Recombination term
    SRHD = TAUP0 * (ng + ni) + TAUN0 * (p + ni)
    SRHL = p / SRHD
    SRHR = ni**2 / SRHD
    
    ASRH = Ucompmass (nodes,n_max,elements,Nelements,SRHL,np.ones(Nelements))
    bSRH = Ucompconst (nodes,n_max,elements,Nelements,SRHR,np.ones(Nelements))
    """
    print("ASRH=",ASRH)
    print("bSRH=",bSRH)
    check_point_18
    """    
    A = A + ASRH
    b = b + bSRH
    """ 
    print("A=",A)
    print("b=",b)
    check_point_19
    """  
    ## Boundary conditions
    b=np.delete(b, BCnodes, 0)
    b[0]         = - A[1,0] * nl
    b[len(b)-1]       =-A[n_max-2,n_max-1] * nr
    A=np.delete(A, BCnodes, 0)
    A=np.delete(A, BCnodes, 1)

    nn= np.linalg.solve(A, b)
    n=np.zeros(n_max)
    n[1:n_max-1]=nn
    n[0]=nl
    n[len(n)-1]=nr
    """
    print("BCnodes=",BCnodes)
    print("A=",A)
    print("b=",b)
    print("n=",n)
    check_point_20    
    """
    return n
