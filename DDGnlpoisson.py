# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 08:52:43 2019


## Copyright (C) 2004-2008  Carlo de Falco
##
## SECS1D - A 1-D Drift--Diffusion Semiconductor Device Simulator
##
##  SECS1D is free software you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation either version 2 of the License, or
##  (at your option) any later version.
##
##  SECS1D is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with SECS1D If not, see <http://www.gnu.org/licenses/>.
##
## author: Carlo de Falco <cdf _AT_ users.sourceforge.net>

## -*- texinfo -*-
##
## @deftypefn {Function File}@
## {@var{V},@var{n},@var{p},@var{res},@var{niter}} = @
## DDGnlpoisson(@var{xaxis},@var{sinodes},@var{Vin},@var{nin},@var{pin},@var{Fnin},@var{Fpin},@var{dop},@var{Ppz_Psp},@var{l2},@var{toll},@var{maxit},@var{verbose})
##
## Solve the non linear Poisson equation
## 
## - lamda^2 *V'' + (n(V,Fn) - p(V,Fp) -dop-Ppz_Psp) = 0 
##
## Input:
## @itemize @minus
## @item xaxis: spatial grid
## @item sinodes: index of the nodes of the grid which are in the
## semiconductor subdomain (remaining nodes are assumed to be in the oxide subdomain)
## @item Vin: initial guess for the electrostatic potential
## @item nin: initial guess for electron concentration
## @item pin: initial guess for hole concentration
## @item Fnin: initial guess for electron Fermi potential
## @item Fpin: initial guess for hole Fermi potential
## @item dop: doping profile
## @item Ppz_Psp: Piezoelectric (Ppz) and Spontious (Psp) built-in polarization charge density profile
## @item l2: scaled electric permittivity (diffusion coefficient)
## @item toll: tolerance for convergence test
## @item maxit: maximum number of Newton iterations
## @item verbose: verbosity level (0,1,2)
## @end itemize
##
## Output:
## @itemize @minus
## @item V: electrostatic potential
## @item n: electron concentration
## @item p: hole concentration
## @item res: residual norm at each step
## @item niter: number of Newton iterations
## @end itemize
##
## @end deftypefn
"""
import numpy as np
from math import*
from scipy import sparse as sp

if __package__:  # explicit relative imports for using aestimo as a package (in python3)
    from . import func_lib
    from .func_lib import DDGphin2n,DDGphip2p,Ucompmass,Ucomplap,Ucompconst
    from .aestimo_poisson1d import equi_np_fi222,equi_np_fi3,equi_np_fi
    from . import config
else:
    import func_lib
    from func_lib import DDGphin2n,DDGphip2p,Ucompmass,Ucomplap,Ucompconst
    from aestimo_poisson1d import equi_np_fi222,equi_np_fi3,equi_np_fi
    import config
    
def  DDGnlpoisson_new (idata,xaxis,sinodes,Vin,nin,pin,toll,maxit,verbose,fi_e,fi_h,model,Vt,surface,fi_stat,iteration,ns):
    ## Set some useful constants
    dampit = 10
    dampcoeff	= 5
    
    ## Convert grid info to FEM form
    #Ndiricheletnodes= 2
    nodes 	= xaxis
    sielements=np.zeros(len(sinodes)-1) 
    sielements[:] = sinodes[1:len(sinodes)]
    n_max	= len(nodes)
    #totdofs = n_max - Ndiricheletnodes
    elements=np.zeros((n_max-1,2))
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    Nelements=np.size(elements[:,0])
    BCnodes= n_max
    normr=np.zeros(maxit+1)
    """
    print("nodes=",nodes)
    print("sielements=",sielements)    
    print("n_max=",n_max)
    print("elements=",elements)
    print("Nelements=",Nelements)
    print("BCnodes=",BCnodes)
    check_point_8    
    """
    ## Initialization
    V = Vin
    EF = 0.0
    if iteration == 1:
        # Determination of the Fermi Level
        EF = 0.0
        V = np.zeros(n_max)
        n, p, V = equi_np_fi(
            iteration, idata.dop, idata.Ppz_Psp, n_max, idata.ni, model, Vt, surface
        )
        fi_stat = V
    else:
        if model.N_wells_virtual - 2 != 0:
            n, p, fi_non, EF = equi_np_fi3(
            V,
            idata.wfh_general,
            idata.wfe_general,
            model,
            idata.E_state_general,
            idata.E_statec_general,
            idata.meff_state_general,
            idata.meff_statec_general,
            n_max,
            idata.ni*ns,
        )
        else:
            n = np.exp(V)
            p = np.exp(-V)       
        n=n*idata.ni
        p=p*idata.ni
    """
    print("sinodes=",sinodes)
    print("n=",n)
    print("p=",p)
    check_point_9
    """
    if (sinodes[0]==0):
        n[1]=nin[0]
        p[1]=pin[0]
    if (sinodes[n_max-1]==n_max-1):
        n[n_max-1]=nin[n_max-1]
        p[n_max-1]=pin[n_max-1]
    
    ## Compute LHS matrices
    l22=idata.l2*np.ones(Nelements)
    L      = Ucomplap (nodes,n_max,elements,Nelements,l22)
    
    ## Compute Mv =  ( n + p)
    Mv            =  np.zeros(n_max)
    Mv[sinodes]   =  (n + p)
    Cv            =  np.ones(Nelements)
    M             =  Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
    
    ## Compute RHS vector
    Tv0            =  np.zeros(n_max)
    Tv0[sinodes]   = (n - p -idata.dop-idata.Ppz_Psp)
    Cv=  np.ones(Nelements)
    T0             =  Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
    
    """
    print('L=',L)
    print('M=',M)
    print('T0=',T0)
    check_point_10
    """
    ## Build LHS matrix and RHS of the linear system for 1st Newton step
    A=np.zeros((n_max,n_max))
    R=np.zeros(n_max)
    Anew=np.zeros((n_max,n_max))
    Rnew=np.zeros(n_max)
    
    A 		= L + M
    LV=np.dot(np.array(L) , V)
    R 		=  LV +T0
    
    ## Apply boundary conditions
    """
    print('A=',A)
    print('R=',R)
    check_point_11
    """
    
    A=np.delete(A, [0,BCnodes-1], 0)
    A=np.delete(A, [0,BCnodes-1], 1)
    R=np.delete(R, [0,BCnodes-1], 0)
    
    normr[0]		=  np.linalg.norm(R,np.inf)
    relresnorm 	= 1
    reldVnorm   = 1
    normrnew	= normr[0]
    
    ## Start of the newton cycle
    for newtit in range(1,maxit):
        if verbose:
            print("\n newton iteration: %d, reldVnorm = %f"%(newtit,reldVnorm))
        
        cc= np.linalg.solve(A, -R)#, rcond=None)[0]
        dV=np.zeros(n_max)
        dV[1:n_max-1] =cc
        ## Start of the damping procedure
        tk = 1
 
        for dit in range(1,dampit):
            if verbose:
                print("\n damping iteration: %d, residual norm = %f"%(dit,normrnew))
            Vnew   = V + tk * dV
            if iteration == 1:
                n = np.exp(Vnew)*idata.ni
                p = np.exp(-Vnew)*idata.ni
            else:
                if model.N_wells_virtual - 2 != 0:
                    n, p, fi_non, EF = equi_np_fi3(
                    Vnew,
                    idata.wfh_general,
                    idata.wfe_general,
                    model,
                    idata.E_state_general,
                    idata.E_statec_general,
                    idata.meff_state_general,
                    idata.meff_statec_general,
                    n_max,
                    idata.ni*ns,
                )
                else:
                    n = np.exp(Vnew)
                    p = np.exp(-Vnew)               
                n=n*idata.ni
                p=p*idata.ni            
            if (sinodes[0]==0):
                n[0]=nin[0]
                p[0]=pin[0]
            if (sinodes[n_max-1]==n_max-1):
                n[n_max-1]=nin[n_max-1]
                p[n_max-1]=pin[n_max-1]
            
            ## Compute LHS matrices
            Mv            =  np.zeros(n_max)
            Mv[sinodes]   =  (n + p)
            Cv            =  np.ones(Nelements)
            #Cv[sielements]=  1        
            M    = Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
            
            ## Compute RHS vector (-residual)
            Tv0            =  np.zeros(n_max)
            Tv0[sinodes]   =  (n - p -idata.dop-idata.Ppz_Psp)
            Cv=  np.ones(Nelements)
            Cv            =  np.ones(Nelements)
            #Cv[sielements]=  1
            T0     = Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
            """ 
            print('L=',L)
            print('M=',M)
            print('T0=',T0)
            check_point_12
            """           
            ## Build LHS matrix and RHS of the linear system for 1st Newton step
            Anew 		= L + M
            LVnew=np.dot(np.array(L) , Vnew)
            Rnew 		=  LVnew +T0
            """
            print('Anew=',Anew)
            print('Rnew=',Rnew)
            check_point_13
            """
            ## Apply boundary conditions
            Anew=np.delete(Anew, [0,BCnodes-1], 0)
            Anew=np.delete(Anew, [0,BCnodes-1], 1)
            Rnew=np.delete(Rnew, [0,BCnodes-1], 0)
            
            if ((dit>1) and (np.linalg.norm(Rnew,np.inf) >= np.linalg.norm(R,np.inf))):
                if verbose:
                    print("\nexiting damping cycle \n")
                break
            else:
                A = Anew
                R = Rnew
        
            ## Compute | R_{k+1} | for the convergence test
            normrnew= np.linalg.norm(R,np.inf)
            
            ## Check if more damping is needed
            if (normrnew > normr[newtit]):
                tk = tk/dampcoeff
            else:
                if verbose:
                    print("\nexiting damping cycle because residual norm = %f \n"%normrnew)
                break
    
        V= Vnew	
        normr[newtit+1] = normrnew
        dVnorm= np.linalg.norm(tk*dV,np.inf)
    
        ## Check if convergence has been reached
        reldVnorm           = dVnorm / np.linalg.norm(V,np.inf)
        if (reldVnorm <= toll):
            if(verbose):
                print("\nexiting newton cycle because reldVnorm= %f \n"%reldVnorm)
            break
    
    res = normr
    niter = newtit
    
    return [V,n,p,fi_stat]#,res,niter

def  DDGnlpoisson (idata,xaxis,sinodes,Vin,nin,pin,Fnin,Fpin,dop,Ppz_Psp,l2,toll,maxit,verbose,ni,fi_e,fi_h,model,Vt):
    ## Set some useful constants
    dampit = 10
    dampcoeff	= 5
    
    ## Convert grid info to FEM form
    #Ndiricheletnodes= 2
    nodes 	= xaxis
    sielements=np.zeros(len(sinodes)-1) 
    sielements[:] = sinodes[1:len(sinodes)]
    n_max	= len(nodes)
    fi_n=np.zeros(n_max)
    fi_p=np.zeros(n_max)
    #totdofs = n_max - Ndiricheletnodes
    elements=np.zeros((n_max-1,2))
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    Nelements=np.size(elements[:,0])
    BCnodes= n_max
    normr=np.zeros(maxit+1)
    """
    print("nodes=",nodes)
    print("sielements=",sielements)    
    print("n_max=",n_max)
    print("elements=",elements)
    print("Nelements=",Nelements)
    print("BCnodes=",BCnodes)
    check_point_8    
    """
    ## Initialization
    V = Vin
    Fn = Fnin
    Fp = Fpin
    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,V,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,idata.n,idata.p)

    n = DDGphin2n(V[sinodes]+fi_n[sinodes],Fn,idata.n)
    p = DDGphip2p(V[sinodes]+fi_p[sinodes],Fp,idata.p)
    """
    print("sinodes=",sinodes)
    print("n=",n)
    print("p=",p)
    check_point_9
    """
    if (sinodes[0]==0):
        n[1]=nin[0]
        p[1]=pin[0]
    if (sinodes[n_max-1]==n_max-1):
        n[n_max-1]=nin[n_max-1]
        p[n_max-1]=pin[n_max-1]
    
    ## Compute LHS matrices
    l22=l2*np.ones(Nelements)
    L      = Ucomplap (nodes,n_max,elements,Nelements,l22)
    
    ## Compute Mv =  ( n + p)
    Mv            =  np.zeros(n_max)
    Mv[sinodes]   =  (n + p)
    Cv            =  np.ones(Nelements)
    M             =  Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
    
    ## Compute RHS vector
    Tv0            =  np.zeros(n_max)
    Tv0[sinodes]   = (n - p -dop-Ppz_Psp)
    Cv=  np.ones(Nelements)
    T0             =  Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
    
    """
    print('L=',L)
    print('M=',M)
    print('T0=',T0)
    check_point_10
    """
    ## Build LHS matrix and RHS of the linear system for 1st Newton step
    A=np.zeros((n_max,n_max))
    R=np.zeros(n_max)
    Anew=np.zeros((n_max,n_max))
    Rnew=np.zeros(n_max)
    
    A 		= L + M
    LV=np.dot(np.array(L) , V)
    R 		=  LV +T0
    
    ## Apply boundary conditions
    """
    print('A=',A)
    print('R=',R)
    check_point_11
    """
    
    A=np.delete(A, [0,BCnodes-1], 0)
    A=np.delete(A, [0,BCnodes-1], 1)
    R=np.delete(R, [0,BCnodes-1], 0)
    
    normr[0]		=  np.linalg.norm(R,np.inf)
    relresnorm 	= 1
    reldVnorm   = 1
    normrnew	= normr[0]
    
    ## Start of the newton cycle
    for newtit in range(1,maxit):
        if verbose:
            print("\n newton iteration: %d, reldVnorm = %f"%(newtit,reldVnorm))
        
        cc= np.linalg.solve(A, -R)#, rcond=None)[0]
        dV=np.zeros(n_max)
        dV[1:n_max-1] =cc
        ## Start of the damping procedure
        tk = 1
 
        for dit in range(1,dampit):
            if verbose:
                print("\n damping iteration: %d, residual norm = %f"%(dit,normrnew))
            Vnew   = V + tk * dV
            if model.N_wells_virtual-2!=0 and config.quantum_effect:
                fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,Vnew,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,n,p)             
            n = DDGphin2n(Vnew[sinodes]+fi_n[sinodes],Fn,idata.n)
            p = DDGphip2p(Vnew[sinodes]+fi_p[sinodes],Fp,idata.p)
            if (sinodes[0]==0):
                n[0]=nin[0]
                p[0]=pin[0]
            if (sinodes[n_max-1]==n_max-1):
                n[n_max-1]=nin[n_max-1]
                p[n_max-1]=pin[n_max-1]
            
            ## Compute LHS matrices
            Mv            =  np.zeros(n_max)
            Mv[sinodes]   =  (n + p)
            Cv            =  np.ones(Nelements)
            #Cv[sielements]=  1        
            M    = Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
            
            ## Compute RHS vector (-residual)
            Tv0            =  np.zeros(n_max)
            Tv0[sinodes]   =  (n - p -dop-Ppz_Psp)
            Cv=  np.ones(Nelements)
            Cv            =  np.ones(Nelements)
            #Cv[sielements]=  1
            T0     = Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
            """ 
            print('L=',L)
            print('M=',M)
            print('T0=',T0)
            check_point_12
            """           
            ## Build LHS matrix and RHS of the linear system for 1st Newton step
            Anew 		= L + M
            LVnew=np.dot(np.array(L) , Vnew)
            Rnew 		=  LVnew +T0
            """
            print('Anew=',Anew)
            print('Rnew=',Rnew)
            check_point_13
            """
            ## Apply boundary conditions
            Anew=np.delete(Anew, [0,BCnodes-1], 0)
            Anew=np.delete(Anew, [0,BCnodes-1], 1)
            Rnew=np.delete(Rnew, [0,BCnodes-1], 0)
            
            if ((dit>1) and (np.linalg.norm(Rnew,np.inf) >= np.linalg.norm(R,np.inf))):
                if verbose:
                    print("\nexiting damping cycle \n")
                break
            else:
                A = Anew
                R = Rnew
        
            ## Compute | R_{k+1} | for the convergence test
            normrnew= np.linalg.norm(R,np.inf)
            
            ## Check if more damping is needed
            if (normrnew > normr[newtit]):
                tk = tk/dampcoeff
            else:
                if verbose:
                    print("\nexiting damping cycle because residual norm = %f \n"%normrnew)
                break
    
        V= Vnew	
        normr[newtit+1] = normrnew
        dVnorm= np.linalg.norm(tk*dV,np.inf)
    
        ## Check if convergence has been reached
        reldVnorm           = dVnorm / np.linalg.norm(V,np.inf)
        if (reldVnorm <= toll):
            if(verbose):
                print("\nexiting newton cycle because reldVnorm= %f \n"%reldVnorm)
            break
    
    res = normr
    niter = newtit
    
    return [V,n,p]#,res,niter
