# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 20:39:48 2019


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
## {@var{p}} = DDGhole_driftdiffusio(@var{psi},@var{xaxis},@var{pg},@var{n},@var{ni},@var{TAUN0},@var{TAUP0},@var{mup})
##
## Solve the continuity equation for holes
##
## Input:
## @itemize @minus
## @item psi: electric potential
## @item xaxis: spatial grid
## @item ng: initial guess and BCs for electron density
## @item n: electron density (for SRH recombination)
## @end itemize
##
## Output:
## @itemize @minus
## @item p: updated hole density
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

def  DDGhole_driftdiffusion(psi,xaxis,pg,n,ni,TAUN0,TAUP0,mup,fi_e,fi_h,model,Vt,idata):
    
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
    
    pl = pg[0]
    pr = pg[n_max-1]
    h=nodes[1:len(nodes)]-nodes[0:len(nodes)-1]
    
    c=1/h
    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,psi,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,n,idata.p)
    Bneg=Ubernoulli(-(psi[1:n_max]-psi[0:n_max-1])-(fi_n[1:n_max]-fi_n[0:n_max-1]),1)
    Bpos=Ubernoulli( (psi[1:n_max]-psi[0:n_max-1])+(fi_p[1:n_max]-fi_p[0:n_max-1]),1)
    
    d0=np.zeros(n_max)
    d0[0]=c[0]*Bneg[0]
    d0[n_max-1]=c[len(c)-1]*Bpos[len(Bpos)-1]
    d0[1:n_max-1]=c[0:len(c)-1]*Bpos[0:len(Bpos)-1]+c[1:len(c)]*Bneg[1:len(Bneg)]    
    
    d1	= np.zeros(n_max)
    d1[0]=n_max
    d1[1:n_max]=-c* Bneg      
    dm1	= np.zeros(n_max)
    dm1[n_max-1]=n_max
    dm1[0:n_max-1]=-c* Bpos   
    A = sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense() 
    b = np.zeros(n_max)#%- A * ng
    
    ## SRH Recombination term
        
    SRHD = TAUP0 * (n + ni) + TAUN0 * (pg + ni)
    SRHL = n / SRHD
    SRHR = ni**2 / SRHD
    
    ASRH = Ucompmass (nodes,n_max,elements,Nelements,SRHL,np.ones(Nelements))
    bSRH = Ucompconst (nodes,n_max,elements,Nelements,SRHR,np.ones(Nelements))
    
    A = A + ASRH
    b = b + bSRH
    
    ## Boundary conditions
    b=np.delete(b, BCnodes, 0)
    b[0]         = - A[1,0] * pl
    b[len(b)-1]       = - A[n_max-2,n_max-1] * pr
    A=np.delete(A, BCnodes, 0)
    A=np.delete(A, BCnodes, 1)
    
    pp= np.linalg.solve(A,b)
    p=np.zeros(n_max)
    p[1:n_max-1]=pp    
    p[0]=pl
    p[len(p)-1]=pr
    return p
