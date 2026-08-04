# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:14:03 2019

"""

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
## {@var{odata},@var{it},@var{res}} = DDNnewtonmap(@var{xaxis},@var{idata},@var{toll},@var{maxit},@var{verbose})
##
## Solve the scaled stationary bipolar DD equation system using a
## coupled Newton algorithm
##
## Input:
## @itemize @minus
## @item xaxis: spatial grid
## @item idata.dop: doping profile
## @item idata.Ppz_Psp: Piezoelectric (Ppz) and Spontious (Psp) built-in polarization charge density profile
## @item idata.p: initial guess for hole concentration
## @item idata.n: initial guess for electron concentration
## @item idata.V: initial guess for electrostatic potential
## @item idata.Fn: initial guess for electron Fermi potential
## @item idata.Fp: initial guess for hole Fermi potential
## @item idata.l2: scaled electric permittivity (diffusion coefficient in Poisson equation)
## @item idata.mun: scaled electron mobility
## @item idata.mup: scaled electron mobility
## @item idata.nis: scaled intrinsic carrier density
## @item idata.TAUN0: scaled electron lifetime
## @item idata.TAUP0: scaled hole lifetime
## @item toll: tolerance for Newton iterarion convergence test
## @item maxit: maximum number of Newton iterarions
## @item verbose: verbosity level: 0,1,2
## @end itemize
##
## Output:
## @itemize @minus
## @item odata.n: electron concentration
## @item odata.p: hole concentration
## @item odata.V: electrostatic potential
## @item odata.Fn: electron Fermi potential
## @item odata.Fp: hole Fermi potential
## @item it: number of Newton iterations performed
## @item res: residual at each step
## @end itemize
##
## @end deftypefn

import numpy as np
from math import*
from scipy import sparse as sp
from scipy.sparse.linalg import spsolve

from .func_lib import Uscharfettergummel,Ucompmass,Ucomplap,Umediaarmonica
from .aestimo_poisson1d import equi_np_fi222
import config


def  DDNnewtonmap (ni,fi_e,fi_h,xaxis,idata,toll,maxit,verbose,model,Vt):
    odata     = idata
    n_max    = len(xaxis)
    fi_n=np.zeros(n_max)
    fi_p=np.zeros(n_max)
    Nelements=n_max-1
    elements=np.zeros((n_max-1,2))
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    BCnodesp = [0, n_max-1]
    if getattr(model, 'photovoltaic_mode', False):
        # Photovoltaic Mode: Selective Contacts
        # Fix majority carriers at their respective contacts based on doping
        BCnodesp1 = []
        BCnodesp2 = []
        dop = getattr(idata, 'dop', None)
        
        if dop is not None:
            if dop[0] > 0: # n-type on left
                BCnodesp1.append(n_max) # Fix Electron at Left
            else: # p-type on left
                BCnodesp2.append(2*n_max) # Fix Hole at Left
                
            if dop[n_max-1] > 0: # n-type on right
                BCnodesp1.append(2*n_max - 1) # Fix Electron at Right
            else: # p-type on right
                BCnodesp2.append(3*n_max - 1) # Fix Hole at Right
        else:
            BCnodesp1 = [n_max] # Fix Electron at Left
            BCnodesp2 = [3*n_max-1] # Fix Hole at Right
        
        BCnodes = np.array(BCnodesp + BCnodesp1 + BCnodesp2)
    else:
        # Standard Mode: Ohmic Contacts (Fix all at both ends)
        BCnodesp1 = [n_max, 2*n_max-1]
        BCnodesp2 = [2*n_max, 3*n_max-1]
        BCnodes_ = np.zeros((3, 2))
        BCnodes_[0, :] = BCnodesp
        BCnodes_[1, :] = BCnodesp1
        BCnodes_[2, :] = BCnodesp2
        BCnodes = BCnodes_.flatten()

    totaldofs= n_max-2
    dampcoef = 10
    maxdamp  = 10
    nrm_du_old=1.
    V = odata.V
    n = odata.n
    p = odata.p
    dop = idata.dop
    Ppz_Psp=idata.Ppz_Psp
    
    ## Create the complete unknown vector
    u = np.hstack(([V, n, p]))
    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,V,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,n,p)    
    ## Build fem matrices
    L = Ucomplap (xaxis,n_max,elements,Nelements,idata.l2*np.ones(Nelements))
    M = Ucompmass (xaxis,n_max,elements,Nelements,np.ones(n_max),np.ones(Nelements))
    DDn = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mun,1,V+fi_n)
    DDp = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mup,1,-V-fi_p)
    
    ## --- Trap-Assisted Tunneling (TAT) Modification ---
    # Calculate electric field for Hurkx model
    # xaxis is scaled (0 to 1), V is scaled (V/Vt)
    dV = np.zeros(n_max)
    dV[1:-1] = (V[2:] - V[0:-2]) / (xaxis[2:] - xaxis[0:-2])
    dV[0] = (V[1] - V[0]) / (xaxis[1] - xaxis[0])
    dV[-1] = (V[-1] - V[-2]) / (xaxis[-1] - xaxis[-2])
    
    # Scaling factors (from aestimo.py)
    # Physical field F [V/m] = (dV_scaled/dx_scaled) * (Vt / xs)
    # where xs is the total length in meters.
    xs = model.dx * n_max
    E_field = np.abs(dV) * (Vt / xs) # V/m
    
    # Hurkx field-enhancement factor Gamma
    # Use a high default for tat_field (e.g. 1e10) if not provided to disable TAT
    tat_field = float(getattr(model, 'tat_field', 1e10)) # V/m
    trap_density_scale = max(getattr(model, 'trap_density_scale', 1.0), 1e-12)
    trap_energy_offset_ev = getattr(model, 'trap_energy_offset_ev', 0.0)
    
    Gamma = np.zeros(n_max)
    mask_high_field = E_field > 1e4 # Threshold to avoid noise (~100 V/cm)
    if tat_field < 1e9:
        ratio = E_field[mask_high_field] / tat_field
        # Gamma = 2 * sqrt(3*pi) * ratio * exp(ratio**2)
        # Simplified phenomenological form often used for soft turn-on:
        Gamma[mask_high_field] = 2.0 * np.sqrt(3.0 * np.pi) * ratio * np.exp(np.clip(ratio**2, 0, 20))
    Gamma *= trap_density_scale
    
    ## Initialise RHS
    theta_n = getattr(idata, 'theta_n', idata.theta)
    theta_p = getattr(idata, 'theta_p', idata.theta)
    denomsrh   = idata.TAUN0 * (p + theta_p) + idata.TAUP0 * (n + theta_n)
    denomsrh   = denomsrh / trap_density_scale
    factauger  = idata.Cn * n + idata.Cp * p
    
    # Apply TAT enhancement to SRH term
    fact       = ((1.0 + Gamma) / denomsrh + factauger)
    
    # Optical Generation Rate (Pre-scaled in idata)
    G_opt = getattr(idata, 'G_optical', 0.0)
         
    r1  = L.dot(V) + M.dot(n - p - dop - Ppz_Psp)
    # Continuity: div(Jn) - R + G = 0  =>  div(Jn) - (R - G) = 0
    r2  = DDn.dot(n) + M.dot((p * n - idata.theta** 2) * fact) - M.dot(np.ones(n_max) * G_opt)
    # Continuity: div(Jp) + R - G = 0  =>  div(Jp) + (R - G) = 0
    r3  = DDp.dot(p) + M.dot((p * n - idata.theta** 2) * fact) - M.dot(np.ones(n_max) * G_opt)
    RHS=-np.hstack(( r1, r2, r3))

    ##  Apply BCs via masking
    mask = np.ones(3 * n_max, dtype=bool)
    mask[BCnodes.astype(int)] = False

    RHS_red = RHS[mask]
    nrm = np.linalg.norm(RHS_red, np.inf)
    res = np.zeros(maxit)
    res[0] = nrm
    ## Begin Newton Cycle
    for count in range (0, maxit):
        if verbose:
          print ("Newton Iteration Number:%d\n"%count)	
        # Ensure positive carrier densities for matrix stability
        n_safe = np.maximum(n, 1e-20)
        p_safe = np.maximum(p, 1e-20)
        
        Ln = Ucomplap (xaxis,n_max,elements,Nelements,Umediaarmonica(idata.mun*n_safe))
        Lp = Ucomplap (xaxis,n_max,elements,Nelements,Umediaarmonica(idata.mup*p_safe))
        Mn = Ucompmass (xaxis,n_max,elements,Nelements,np.ones(n_max),n_safe[0:n_max-1]*fact[0:n_max-1])
        Mp = Ucompmass (xaxis,n_max,elements,Nelements,np.ones(n_max),p_safe[0:n_max-1]*fact[0:n_max-1])
        Z  = np.zeros((n_max,n_max))   
        DDn = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mun,1,V+fi_n)
        DDp = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mup,1,-V-fi_p)
        A 	= L  #A11
        B	= M #A12
        C	=-M #A13
        DDD	=-Ln #A21
        E	= DDn+Mp#A22
        F	= Mn  #A23 (Removed Z+ since Z is sparse and zero)
        G	= Lp #A31
        H	= Mp #A32
        I	= DDp+Mn#A33
        ## Build LHS
        LHS = sp.bmat([(A, B, C), (DDD, E, F), (G, H, I)], format='csc')
        
        ## Apply BCs using masking
        A_red = LHS[mask, :][:, mask]
        # RHS_red is updated in the loop damping cycle if needed, 
        # but for count=0 we use the initial RHS_red
        if count > 0:
            RHS_red = RHS[mask]
        
        ## Solve the linearised system
        dutmp_reduced = spsolve(A_red, RHS_red) # Newton step is A*du = -R, and RHS_red is already -R
        
        # Reconstruct full dutmp
        dutmp = np.zeros(LHS.shape[0])
        dutmp[mask] = dutmp_reduced
        du = dutmp
        ## Check Convergence
        nrm_u = np.linalg.norm(u,np.inf)
        nrm_du = np.linalg.norm(du,np.inf)
    	
        ratio = nrm_du/nrm_u 
        if verbose:
          print ("ratio = %e\n"% ratio)		
        
        if (ratio <= toll):
            V 	 = u[0:n_max]
            n	    = u[n_max:2*n_max]
            p	    = u[2*n_max:len(u)]
            res[count]  = nrm
            break
        ## Begin damping cycle
        tj = 1
        if model.N_wells_virtual-2!=0 and config.quantum_effect:
            fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,V,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,n,p)
        for cc in range( 1,maxdamp):
          if verbose:
            print ("damping iteration number:%d\n"%cc)
            print ("reference residual norm:%f\n"%nrm)
          
          ## Update the unknown vector		
          utmp    = u + tj*du
          Vnew 	    = utmp[0:n_max]
          # Enforce carrier floor during damping
          nnew	    = np.maximum(utmp[n_max:2*n_max], 1e-20)
          pnew	    = np.maximum(utmp[2*n_max:3*n_max], 1e-20)
          ## Try a new RHS
          
          DDn = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mun,1,Vnew+fi_n)
          DDp = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mup,1,-Vnew-fi_p)
          
          r1  = L.dot(Vnew) + M.dot(nnew - pnew - dop - Ppz_Psp)
          r2  = DDn.dot(nnew) + M.dot((pnew * nnew - idata.theta** 2) * fact) - M.dot(np.ones(n_max) * G_opt)
          r3  = DDp.dot(pnew) + M.dot((pnew * nnew - idata.theta** 2) * fact) - M.dot(np.ones(n_max) * G_opt)
          RHS = -np.hstack(( r1, r2, r3))

          ## Apply BCs
          RHS_red_tmp = RHS[mask]
          nrmtmp = np.linalg.norm(RHS_red_tmp, np.inf)
          
          ## Update the damping coefficient
          if verbose:
              print("residual norm:%f\n\n"%nrmtmp)
            
          if (nrmtmp > nrm):
              tj = tj/(dampcoef*cc)
              if verbose:                  
                  print ("\ndamping coefficients = %f"%tj)
          else:
              RHS_red = RHS_red_tmp # Accept new residual
              break
        nrm_du = np.linalg.norm(tj*du,np.inf)
        u 	= utmp
        
        if (count>0):
            ratio = nrm_du/nrm_du_old
            if (ratio<.005):
                V 	    = u[0:n_max]
                n	    = np.maximum(u[n_max:2*n_max], 1e-20)
                p	    = np.maximum(u[2*n_max:3*n_max], 1e-20)            
                res[count]  = nrm
                break           
        nrm = nrmtmp
        res[count]  = nrm
        ## Convert result vector into distinct output vectors 
        V 	    = u[0:n_max]
        n	    = np.maximum(u[n_max:2*n_max], 1e-20)
        p	    = np.maximum(u[2*n_max:3*n_max], 1e-20)    
        nrm_du_old = nrm_du
    odata.V = V
    odata.n = n
    odata.p = p
    

    it   = count

    return [odata,it,res]
