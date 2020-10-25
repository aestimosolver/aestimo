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
from scipy.sparse import bsr_matrix

if __package__:  # explicit relative imports for using aestimo as a package (in python3)
    from . import func_lib
    from .func_lib import Uscharfettergummel,Ucompmass,Ucomplap,Umediaarmonica
    from .aestimo_poisson1d import equi_np_fi222
    from . import config
else:    
    import func_lib
    from func_lib import Uscharfettergummel,Ucompmass,Ucomplap,Umediaarmonica
    from aestimo_poisson1d import equi_np_fi222
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
    BCnodesp = [0,n_max-1]
    BCnodesp1 = [n_max,2*n_max-1]
    BCnodesp2 = [2*n_max,3*n_max-1]
    BCnodes_=np.zeros((3,2))
    BCnodes_[0,:]=BCnodesp
    BCnodes_[1,:]=BCnodesp1
    BCnodes_[2,:]=BCnodesp2
    BCnodes=BCnodes_.flatten()

    totaldofs= n_max-2
    dampcoef = 10
    maxdamp  = 2
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
    
    ## Initialise RHS
    denomsrh   = idata.TAUN0 * (p + idata.theta) + idata.TAUP0 * (n + idata.theta)
    factauger  = idata.Cn * n + idata.Cp * p
    fact       = (1 / denomsrh + factauger)
    r1  = np.dot(np.array(L) , V) + np.dot(np.array(M) , (n - p - dop-Ppz_Psp))
    r2  = np.dot(np.array(DDn) , n)+ np.dot(np.array(M) , (p * n - idata.theta** 2) * fact)
    r3  = np.dot(np.array(DDp) , p)+ np.dot(np.array(M) , (p * n - idata.theta** 2) * fact)
    RHS=-np.hstack(( r1, r2, r3))

    ##  Apply BCs
    RHS=np.delete(RHS, BCnodes, 0)
    nrm = np.linalg.norm(RHS,np.inf)
    res=np.zeros(maxit)
    res[0] = nrm
    ## Begin Newton Cycle
    for count in range (0, maxit):
        if verbose:
          print ("Newton Iteration Number:%d\n"%count)	
        Ln = Ucomplap (xaxis,n_max,elements,Nelements,Umediaarmonica(idata.mun*n))
        Lp = Ucomplap (xaxis,n_max,elements,Nelements,Umediaarmonica(idata.mup*p))
        Mn = Ucompmass (xaxis,n_max,elements,Nelements,np.ones(n_max),n[0:n_max-1]*fact[0:n_max-1])
        Mp = Ucompmass (xaxis,n_max,elements,Nelements,np.ones(n_max),p[0:n_max-1]*fact[0:n_max-1])
        Z  = np.zeros((n_max,n_max))   
        DDn = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mun,1,V+fi_n)
        DDp = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mup,1,-V-fi_p)
        A 	= L  #A11
        B	= M #A12
        C	=-M #A13
        DDD	=-Ln #A21
        E	= DDn+Mp#A22
        F	= Z+Mn  #A23
        G	= Lp #A31
        H	= Z+Mp #A32
        I	= DDp+Mn#A33
        ## Build LHS
        LHS= np.asarray(np.bmat([(A,	B, C),(DDD, E, F),(G, H, I)]))
        ## Apply BCs
        LHS=np.delete(LHS, BCnodes, 0)
        LHS=np.delete(LHS, BCnodes, 1)        
        ## Solve the linearised system
        dutmp= np.linalg.solve(LHS, RHS)#, rcond=None)[0]
        dv    = dutmp[0:totaldofs]
        dn    = dutmp[totaldofs:2*totaldofs]
        dp    = dutmp[2*totaldofs:3*totaldofs]
        du=np.hstack((0,dv,0,0,dn,0,0,dp,0))
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
          nnew	    = utmp[n_max:2*n_max]
          pnew	    = utmp[2*n_max:len(utmp)]
          ## Try a new RHS
          
          DDn = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mun,1,Vnew+fi_n)
          DDp = Uscharfettergummel(xaxis,n_max,elements,Nelements,idata.mup,1,-Vnew-fi_p)
          
          r1  = np.dot(np.array(L) , V) + np.dot(np.array(M) , (nnew - pnew - dop-Ppz_Psp))
          r2  = np.dot(np.array(DDn) , nnew)+ np.dot(np.array(M) , (pnew * nnew - idata.theta** 2) * fact)
          r3  = np.dot(np.array(DDp) , pnew) + np.dot(np.array(M) , (pnew * nnew - idata.theta** 2) * fact)
          RHS=-np.hstack(( r1, r2, r3))

          ## Apply BCs
          RHS=np.delete(RHS, BCnodes, 0)
          nrmtmp=np.linalg.norm(RHS,np.inf)
          
          ## Update the damping coefficient
          if verbose:
              print("residual norm:%f\n\n"%nrmtmp)
            
          if (nrmtmp>nrm):
              tj = tj/(dampcoef*cc)
              if verbose:                  
                  print ("\ndamping coefficients = %f"%tj)
          else:
              break
        nrm_du = np.linalg.norm(tj*du,np.inf)
        u 	= utmp
        
        if (count>0):
            ratio = nrm_du/nrm_du_old
            if (ratio<.005):
                V 	    = u[0:n_max]
                n	    = u[n_max:2*n_max]
                p	    = u[2*n_max:len(u)]            
                res[count]  = nrm
                break           
        nrm = nrmtmp
        res[count]  = nrm
        ## Convert result vector into distinct output vectors 
        V 	    = u[0:n_max]
        n	    = u[n_max:2*n_max]
        p	    = u[2*n_max:len(u)]    
        nrm_du_old = nrm_du
    odata.V = V
    odata.n = n
    odata.p = p
    

    it   = count

    return [odata,it,res]
