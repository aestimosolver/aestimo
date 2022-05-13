# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 07:27:07 2019

"""

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
## {@var{odata},@var{it},@var{res}} =  DDGgummelmap(@var{xaxis},@var{idata},@var{toll},@var{maxit},@var{ptoll},@var{pmaxit},@var{verbose})
##
## Solve the scaled stationary bipolar DD equation system using Gummel
## algorithm
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
## @item toll: tolerance for Gummel iterarion convergence test
## @item maxit: maximum number of Gummel iterarions
## @item ptoll: tolerance for Newton iterarion convergence test for non linear Poisson
## @item pmaxit: maximum number of Newton iterarions
## @item verbose: verbosity level (0,1,2)
## @end itemize
##
## Output:
## @itemize @minus
## @item odata.n: electron concentration
## @item odata.p: hole concentration
## @item odata.V: electrostatic potential
## @item odata.Fn: electron Fermi potential
## @item odata.Fp: hole Fermi potential
## @item it: number of Gummel iterations performed
## @item res: total potential increment at each step
## @end itemize
##
## @end deftypefn
import numpy as np
      
from .ddgnlpoisson import DDGnlpoisson
from .ddgelectron_driftdiffusion import DDGelectron_driftdiffusion
from .ddghole_driftdiffusion import DDGhole_driftdiffusion
from .func_lib import DDGp2phip,DDGn2phin

def DDGgummelmap (n_max,xaxis,idata,odata,toll,maxit,ptoll,pmaxit,verbose,ni,fi_e,fi_h,model,Vt):

    odata  = idata
    dop         = idata.dop
    Ppz_Psp= idata.Ppz_Psp
    vout=np.zeros((n_max,2))
    hole_density=np.zeros((n_max,3))
    electron_density=np.zeros((n_max,3))
    fermin=np.zeros((n_max,2))
    fermip=np.zeros((n_max,2))
    v_Nnodes=np.arange(n_max)
    
    vout[:,0] = idata.V
    hole_density [:,0] = idata.p
    electron_density [:,0]= idata.n
    fermin [:,0]=idata.Fn
    fermip [:,0]=idata.Fp
    nrm=np.zeros(maxit)
    for i in range (1,maxit):
        if (verbose>1):
          print(1,"*****************************************************************\n")  
          print(1,"****    start of gummel iteration number: %d\n",i)
          print(1,"*****************************************************************\n")        
        if (verbose>1):
            print(1,"solving non linear poisson equation\n\n")
        
        #print("here_1","i=",i,"/",maxit)                                
        [vout[:,1],electron_density[:,1],hole_density[:,1]] =DDGnlpoisson (idata,xaxis,v_Nnodes,vout[:,0],electron_density[:,0],hole_density[:,0],fermin[:,0],fermip[:,0],dop,Ppz_Psp,idata.l2,ptoll,pmaxit,verbose,ni,fi_e,fi_h,model,Vt)
        
        #print("here_2")
        
        """
        print("vout=",vout[:,1])
        print("electron_density=",electron_density[:,1])
        print("hole_density=",hole_density[:,1])
        
        
        check_point_7 
        """
                                                        	
        if (verbose>1):
          print (1,"\n\nupdating electron qfl\n\n")
        electron_density[:,2]=DDGelectron_driftdiffusion(vout[:,1], xaxis, electron_density[:,1],hole_density[:,1],idata.nis,idata.TAUN0,idata.TAUP0,idata.mun,fi_e,fi_h,model,Vt,idata)
        
        fermin[:,1] = DDGn2phin(vout[:,1],electron_density[:,2])
        fermin[0,1]   = idata.Fn[0]
        fermin[n_max-1,1] = idata.Fn[len(idata.Fn)-1]
        
       
        if (verbose>1):
          print("updating hole qfl\n\n")
        hole_density[:,2] = DDGhole_driftdiffusion(vout[:,1], xaxis, hole_density[:,1],electron_density[:,1],idata.nis,idata.TAUN0,idata.TAUP0,idata.mup,fi_e,fi_h,model,Vt,idata)
        
        fermip[:,1] = DDGp2phip(vout[:,1],hole_density[:,2])
        fermip[0,1]   = idata.Fp[0]
        """ 
        print("vout=",vout[:,1])
        print("electron_density=",electron_density[:,2])
        print("hole_density=",hole_density[:,2])
        
        
        check_point_14 
        """
        fermip[n_max-1,1] = idata.Fp[len(idata.Fp)-1]
        
        if (verbose>1):
          print("checking for convergence\n\n")
        nrfn= np.linalg.norm(fermin[:,1]-fermin[:,0],np.inf)
        nrfp= np.linalg.norm (fermip[:,1]-fermip[:,0],np.inf)
        nrv = np.linalg.norm (vout[:,1]-vout[:,0],np.inf)
        nrm[i] = max([nrfn,nrfp,nrv])
        if (verbose>1):
          print (" max(|phin_(k+1)-phinn_(k)| , |phip_(k+1)-phip_(k)| , |v_(k+1)- v_(k)| )= %d \n"%nrm[i])
        #print("norm=",nrm[i],toll)
        if (nrm[i]<toll):
        		break
        
        vout[:,0] = vout[:,1]
        hole_density [:,0] = hole_density [:,2] 
        electron_density [:,0]= electron_density [:,2]
        fermin [:,0]= fermin [:,1]
        fermip [:,0]= fermip [:,1]
        if(verbose and 1==2):
            DDGplotresults(xaxis,electron_density,hole_density,vout,fermin,fermip)		
    
    it = i
    res = nrm
    
    if (verbose>0):
        print("\n\nInitial guess computed by DD: # of Gummel iterations = %d \n\n"%it)
        
    odata.n     = electron_density[:,2]
    odata.p     = hole_density[:,2]
    odata.V     = vout[:,1]
    odata.Fn    = fermin[:,1]
    odata.Fp    = fermip[:,1]
 
    
    return [odata,it,res]
