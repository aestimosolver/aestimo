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
from scipy.sparse.linalg import spsolve

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
    
    if getattr(model, 'photovoltaic_mode', False):
        # Selective contact: fix electron density ONLY at n-type regions (cathode)
        dop = getattr(idata, 'dop', None)
        BCnodes = []
        if dop is not None:
            if dop[0] > 0: BCnodes.append(0)
            if dop[n_max-1] > 0: BCnodes.append(n_max-1)
        if not BCnodes: BCnodes = [0] # Fail-safe
    else:
        # standard ohmic: fix both
        BCnodes= [0,n_max-1]
    
    nl = ng[0]
    nr = ng[n_max-1]
    h=nodes[1:n_max]-nodes[0:n_max-1]
    
    mun_mid = (mun[0:n_max-1] + mun[1:n_max]) / 2.0
    c = mun_mid / h
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
    Bpos=Ubernoulli( (psi[1:n_max]-psi[0:n_max-1])+(fi_n[1:n_max]-fi_n[0:n_max-1]),1)
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
    A = sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).tocsc() 
    b = np.zeros(n_max)
    
    ## Trap-Assisted Tunneling (TAT) Modification
    dV = np.zeros(n_max)
    dV[1:-1] = (psi[2:] - psi[0:-2]) / (nodes[2:] - nodes[0:-2])
    dV[0] = (psi[1] - psi[0]) / (nodes[1] - nodes[0])
    dV[-1] = (psi[-1] - psi[-2]) / (nodes[-1] - nodes[-2])
    
    # Scaling for electric field (V/m)
    # len_ (xbar) is in meters. xaxis is normalized [0, 1]
    xs_val = getattr(model, 'dx', 1.0) * n_max * 1e-9 # meters
    E_field = np.abs(dV) * (Vt / (xs_val if xs_val > 0 else 1e-9))
    
    # Hurkx factor Gamma
    tat_field = getattr(model, 'tat_field', 1e10)
    trap_density_scale = max(getattr(model, 'trap_density_scale', 1.0), 1e-12)
    trap_energy_offset_ev = getattr(model, 'trap_energy_offset_ev', 0.0)
    Gamma = np.zeros(n_max)
    if tat_field < 1e9:
        mask_high_field = E_field > 1e4
        ratio = E_field[mask_high_field] / tat_field
        Gamma[mask_high_field] = 2.0 * np.sqrt(3.0 * np.pi) * ratio * np.exp(np.clip(ratio**2, 0, 20))
    Gamma *= trap_density_scale

    ## SRH Recombination term
    trap_arg = np.clip(trap_energy_offset_ev / max(Vt, 1e-12), -40.0, 40.0)
    n1 = ni * np.exp(trap_arg)
    p1 = ni * np.exp(-trap_arg)
    SRHD = (TAUP0 * (ng + n1) + TAUN0 * (p + p1)) / ((1.0 + Gamma) * trap_density_scale)
    SRHL = p / SRHD
    SRHR = ni**2 / SRHD
    
    ASRH = Ucompmass (nodes,n_max,elements,Nelements,SRHL,np.ones(Nelements))
    bSRH = Ucompconst (nodes,n_max,elements,Nelements,SRHR,np.ones(Nelements))
    
    ## Optical Generation
    G_opt = getattr(idata, 'G_optical', 0.0)
    bG = Ucompconst (nodes,n_max,elements,Nelements,np.ones(n_max) * G_opt,np.ones(Nelements))
    
    A = A + ASRH
    b = b + bSRH + bG
    
    ## Boundary conditions
    mask = np.ones(n_max, dtype=bool)
    mask[BCnodes] = False
    
    b_red = b[mask]
    # For Dirichlet BCs, move known boundary terms to the RHS
    # Only apply if the node is indeed a boundary node that we've masked out
    if 0 in BCnodes:
        # Node 1 is the first free node, it has a contribution from fixed node 0
        b_red[0] = b_red[0] - A[mask, 0].toarray().flatten()[0] * nl
    
    if (n_max-1) in BCnodes:
        # Node n_max-2 is the last free node if n_max-1 is fixed
        b_red[-1] = b_red[-1] - A[mask, n_max-1].toarray().flatten()[-1] * nr
    
    A_red = A[mask, :][:, mask]
    
    nn = spsolve(A_red, b_red)
    
    # Robustness: Check for NaN or Inf results
    if np.any(np.isnan(nn)) or np.any(np.isinf(nn)):
        # Return initial guess as a safe fallback
        return ng
        
    n = np.zeros(n_max)
    n[mask] = nn
    n[0]=nl
    n[len(n)-1]=nr
    return n
