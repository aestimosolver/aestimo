# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 08:52:43 2019

## Copyright (C) 2004-2008  Carlo de Falco
"""
import numpy as np
from math import*
from scipy import sparse as sp
from scipy.sparse.linalg import spsolve

from .func_lib import DDGphin2n,DDGphip2p,Ucompmass,Ucomplap,Ucompconst
from .aestimo_poisson1d import equi_np_fi222,equi_np_fi3,equi_np_fi
import config

def  DDGnlpoisson_new (idata,xaxis,sinodes,Vin,nin,pin,toll,maxit,verbose,fi_e,fi_h,model,Vt,surface,fi_stat,iteration,ns):
    ## Set some useful constants
    dampit = 10
    dampcoeff	= 5
    
    nodes 	= xaxis
    n_max	= len(nodes)
    elements=np.zeros((n_max-1,2))
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    Nelements=n_max-1
    BCnodes= n_max
    normr=np.zeros(maxit+1)

    ## Initialization
    V = Vin
    if iteration == 1:
        V = np.zeros(n_max)
        n, p, V = equi_np_fi(iteration, idata.dop, idata.Ppz_Psp, n_max, idata.ni, model, Vt, surface)
        fi_stat = V
    else:
        if model.N_wells_virtual - 2 != 0:
            n, p, fi_non, EF = equi_np_fi3(V, idata.wfh_general, idata.wfe_general, model, idata.E_state_general, idata.E_statec_general, idata.meff_state_general, idata.meff_statec_general, n_max, idata.ni*ns)
        else:
            n = np.exp(V)
            p = np.exp(-V)       
        n=n*idata.ni
        p=p*idata.ni

    if (sinodes[0]==0):
        n[1]=nin[0]
        p[1]=pin[0]
    if (sinodes[n_max-1]==n_max-1):
        n[n_max-1]=nin[n_max-1]
        p[n_max-1]=pin[n_max-1]
    
    ## Mask for Boundary Conditions
    mask = np.ones(n_max, dtype=bool)
    mask[[0, n_max-1]] = False

    ## Compute LHS matrices
    l22 = idata.l2*np.ones(Nelements)
    L = Ucomplap (nodes,n_max,elements,Nelements,l22)
    
    Mv = np.zeros(n_max)
    Mv[sinodes] = (n + p)
    Cv = np.ones(Nelements)
    M = Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
    
    Tv0 = np.zeros(n_max)
    Tv0[sinodes] = (n - p -idata.dop-idata.Ppz_Psp)
    T0 = Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
    
    A = L + M
    R = L.dot(V) + T0
    
    A_red = A[mask, :][:, mask]
    R_red = R[mask]
    
    normr[0] = np.linalg.norm(R_red, np.inf)
    reldVnorm = 1
    normrnew = normr[0]
    
    for newtit in range(1, maxit):
        if verbose: print("\n newton iteration: %d, reldVnorm = %f"%(newtit,reldVnorm))
        
        cc = spsolve(A_red, -R_red)
        dV = np.zeros(n_max)
        dV[mask] = cc
        tk = 1

        for dit in range(1, dampit):
            if verbose: print("\n damping iteration: %d, residual norm = %f"%(dit,normrnew))
            Vnew = V + tk * dV
            if iteration == 1:
                n = np.exp(Vnew)*idata.ni
                p = np.exp(-Vnew)*idata.ni
            else:
                if model.N_wells_virtual - 2 != 0:
                    n, p, fi_non, EF = equi_np_fi3(Vnew, idata.wfh_general, idata.wfe_general, model, idata.E_state_general, idata.E_statec_general, idata.meff_state_general, idata.meff_statec_general, n_max, idata.ni*ns)
                else:
                    n = np.exp(Vnew)
                    p = np.exp(-Vnew)               
                n=n*idata.ni
                p=p*idata.ni            
            
            if (sinodes[0]==0): n[0]=nin[0]; p[0]=pin[0]
            if (sinodes[n_max-1]==n_max-1): n[n_max-1]=nin[n_max-1]; p[n_max-1]=pin[n_max-1]
            
            Mv = np.zeros(n_max)
            Mv[sinodes] = (n + p)
            M = Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
            
            Tv0 = np.zeros(n_max)
            Tv0[sinodes] = (n - p -idata.dop-idata.Ppz_Psp)
            T0 = Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
            
            Anew = L + M
            Rnew = L.dot(Vnew) + T0
            
            Anew_red = Anew[mask, :][:, mask]
            Rnew_red = Rnew[mask]
            
            if (dit > 1 and np.linalg.norm(Rnew_red, np.inf) >= np.linalg.norm(R_red, np.inf)):
                if verbose: print("\nexiting damping cycle \n")
                break
            else:
                A_red = Anew_red
                R_red = Anew_red # Wait, should be Rnew_red
                R_red = Rnew_red # Corrected
                A = Anew
                R = Rnew
        
            normrnew = np.linalg.norm(R_red, np.inf)
            if (normrnew > normr[newtit]):
                tk = tk/dampcoeff
            else:
                if verbose: print("\nexiting damping cycle because residual norm = %f \n"%normrnew)
                break
    
        V = Vnew	
        normr[newtit+1] = normrnew
        dVnorm = np.linalg.norm(tk*dV, np.inf)
        reldVnorm = dVnorm / (np.linalg.norm(V, np.inf) + 1e-10)
        if (reldVnorm <= toll):
            if verbose: print("\nexiting newton cycle because reldVnorm= %f \n"%reldVnorm)
            break
    
    return [V,n,p,fi_stat]

def  DDGnlpoisson (idata,xaxis,sinodes,Vin,nin,pin,Fnin,Fpin,dop,Ppz_Psp,l2,toll,maxit,verbose,ni,fi_e,fi_h,model,Vt):
    dampit = 10
    dampcoeff = 5
    nodes = xaxis
    n_max = len(nodes)
    elements = np.zeros((n_max-1,2))
    elements[:,0] = np.arange(0,n_max-1)
    elements[:,1] = np.arange(1,n_max)
    Nelements = n_max-1
    BCnodes = n_max
    normr = np.zeros(maxit+1)

    V = Vin
    Fn = Fnin
    Fp = Fpin
    fi_n = np.zeros(n_max)
    fi_p = np.zeros(n_max)

    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,V,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,idata.n,idata.p)

    n = DDGphin2n(V[sinodes]+fi_n[sinodes],Fn,idata.n)
    p = DDGphip2p(V[sinodes]+fi_p[sinodes],Fp,idata.p)

    if (sinodes[0]==0): n[1]=nin[0]; p[1]=pin[0]
    if (sinodes[n_max-1]==n_max-1): n[n_max-1]=nin[n_max-1]; p[n_max-1]=pin[n_max-1]
    
    mask = np.ones(n_max, dtype=bool)
    mask[[0, n_max-1]] = False

    l22 = l2*np.ones(Nelements)
    L = Ucomplap (nodes,n_max,elements,Nelements,l22)
    
    Mv = np.zeros(n_max)
    Mv[sinodes] = (n + p)
    Cv = np.ones(Nelements)
    M = Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
    
    Tv0 = np.zeros(n_max)
    Tv0[sinodes] = (n - p -dop-Ppz_Psp)
    T0 = Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
    
    A = L + M
    R = L.dot(V) + T0
    
    A_red = A[mask, :][:, mask]
    R_red = R[mask]
    
    normr[0] = np.linalg.norm(R_red, np.inf)
    reldVnorm = 1
    normrnew = normr[0]
    
    for newtit in range(1, maxit):
        if verbose: print("\n newton iteration: %d, reldVnorm = %f"%(newtit,reldVnorm))
        
        cc = spsolve(A_red, -R_red)
        
        # Robustness: Check for NaN or Inf results
        if np.any(np.isnan(cc)) or np.any(np.isinf(cc)):
            if verbose: print("\n WARNING: Non-physical result (NaN/Inf) in NLPoisson. Jumping out of Newton cycle.")
            break
            
        dV = np.zeros(n_max)
        dV[mask] = cc
        tk = 1

        for dit in range(1, dampit):
            if verbose: print("\n damping iteration: %d, residual norm = %f"%(dit,normrnew))
            Vnew = V + tk * dV
            if model.N_wells_virtual-2!=0 and config.quantum_effect:
                fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,Vnew,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,n,p)
            
            n = DDGphin2n(Vnew[sinodes]+fi_n[sinodes],Fn,idata.n)
            p = DDGphip2p(Vnew[sinodes]+fi_p[sinodes],Fp,idata.p)
            
            if (sinodes[0]==0): n[0]=nin[0]; p[0]=pin[0]
            if (sinodes[n_max-1]==n_max-1): n[n_max-1]=nin[n_max-1]; p[n_max-1]=pin[n_max-1]
            
            Mv = np.zeros(n_max)
            Mv[sinodes] = (n + p)
            M = Ucompmass (nodes,n_max,elements,Nelements,Mv,Cv)
            
            Tv0 = np.zeros(n_max)
            Tv0[sinodes] = (n - p -dop-Ppz_Psp)
            T0 = Ucompconst (nodes,n_max,elements,Nelements,Tv0,Cv)
            
            Anew = L + M
            Rnew = L.dot(Vnew) + T0
            
            Anew_red = Anew[mask, :][:, mask]
            Rnew_red = Rnew[mask]
            
            if (dit > 1 and np.linalg.norm(Rnew_red, np.inf) >= np.linalg.norm(R_red, np.inf)):
                if verbose: print("\nexiting damping cycle \n")
                break
            else:
                A_red = Anew_red
                R_red = Rnew_red
                A = Anew
                R = Rnew
        
            normrnew = np.linalg.norm(R_red, np.inf)
            if (normrnew > normr[newtit]):
                tk = tk/dampcoeff
            else:
                if verbose: print("\nexiting damping cycle because residual norm = %f \n"%normrnew)
                break
    
        V = Vnew	
        normr[newtit+1] = normrnew
        dVnorm = np.linalg.norm(tk*dV, np.inf)
        reldVnorm = dVnorm / (np.linalg.norm(V, np.inf) + 1e-10)
        if (reldVnorm <= toll):
            if verbose: print("\nexiting newton cycle because reldVnorm= %f \n"%reldVnorm)
            break
            
    return [V,n,p]
