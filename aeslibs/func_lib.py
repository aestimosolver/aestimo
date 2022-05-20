# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:09:15 2019


"""
import numpy as np
from math import*
from scipy import sparse as sp

def DDGphin2n (V,phin,n):
    """
    ## -*- texinfo -*-
    ##
    ## @deftypefn {Function File}@
    ## {@var{n}} = DDGphin2n(@var{V},@var{phin})
    ##
    ## Compute the electron density using Maxwell-Boltzmann statistics
    ##
    ## @end deftypefn
    """
    nmin = 0
    n=  np.exp ((V-phin))
    #n = n * (n>nmin) + nmin * (n<=nmin)
    return n

def DDGphip2p (V,phip,p):
    """
    ## -*- texinfo -*-
    ##
    ## @deftypefn {Function File}@
    ## {@var{p}} = DDGphip2p(@var{V},@var{phip})
    ##
    ## Compute the hole density using Maxwell-Boltzmann statistic
    ##
    ## @end deftypefn
    """
    ## Load constants
    pmin = 0
    p=  np.exp ((phip-V))
    #p = p * (p>pmin) + pmin * (p<=pmin) 
    return p

def DDGn2phin (V,n):
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
    ## {@var{phin}} = DDGn2phin(@var{V},@var{n})
    ##
    ## Compute the qfl for electrons using Maxwell-Boltzmann statistics.
    ##
    ## @end deftypefn
    """
    ## Load constants
    nmin = 0
    #n    = n * (n>nmin) + nmin * (n<=nmin)
    phin = V - np.log(abs(n))
    
    return phin

def DDGp2phip (V,p):
    """
    ## Copyright (C) 2004-2008  Carlo de Falco
    ##
    ## SECS1D - A 1-D Drift--Diffusion Semiconductor Device Simulator
    ##
    ## SECS1D is free software; you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation; either version 2 of the License, or
    ## (at your option) any later version.
    ##
    ## SECS1D is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.
    ##
    ## You should have received a copy of the GNU General Public License
    ## along with SECS1D; If not, see <http://www.gnu.org/licenses/>.
    ##
    ## author: Carlo de Falco <cdf _AT_ users.sourceforge.net>
    
    ## -*- texinfo -*-
    ##
    ## @deftypefn {Function File}@
    ## {@var{phip}} = DDGn2phin(@var{V},@var{p})
    ##
    ## Compute the qfl for holes using Maxwell-Boltzmann statistics
    ##
    ## @end deftypefn
    """
    ## Load constants
    pmin = 0
    #p    = p * (p>pmin) + pmin * (p<=pmin)
    phip = V + np.log(abs(p)) 
      
    return phip

def Ucompmass (nodes,n_max,elements,Nelements,Bvect,Cvect):
    """
    ## -*- texinfo -*-
    ##
    ## @deftypefn {Function File}@
    ## {@var{R}} = Ucompmass(@var{nodes},@var{n_max},@var{elements},@var{Nelements},@var{Bvect},@var{Cvect})
    ##
    ## Compute P1 finite element mass-matrix:
    ##
    ## @itemize @minus
    ## @item @var{nodes}: list of mesh nodes
    ## @item @var{n_max}: number of mesh nodes
    ## @item @var{elements}: list of mesh elements 
    ## @item @var{Nelements}: number of mesh elements
    ## @item @var{Bvect}: piecewise linear reaction coefficient
    ## @item @var{Cvect}: piecewise constant reaction coefficient
    ## @end itemize
    ##
    ## @end deftypefn
    """
    h 	= (nodes[1:len(nodes)]-nodes[0:len(nodes)-1])*Cvect  
    d0=np.zeros(n_max)
    d0[0]=Bvect[0]*h[0]/2
    d0[n_max-1]=Bvect[n_max-1]*h[len(h)-1]/2
    d0[1:n_max-1]=Bvect[1:n_max-1]*(h[0:len(h)-1]+h[1:len(h)])/2
    Bmat  = sp.spdiags(d0, [0], n_max,n_max).todense()      
    return Bmat

def Ucomplap (nodes,n_max,elements,Nelements,coeff):
    """
    ## -*- texinfo -*-
    ##
    ## @deftypefn {Function File}@
    ## {@var{R}} = Ucomplap(@var{nodes},@var{n_max},@var{elements},@var{Nelements},@var{coeff})
    ##
    ## Compute P1 finite element approximation of the differential operator:
    ## 
    ##  - d ( coeff d (.)\dx)\dx
    ##
    ## @itemize @minus
    ## @item @var{nodes}: list of mesh nodes
    ## @item @var{n_max}: number of mesh nodes
    ## @item @var{elements}: list of mesh elements 
    ## @item @var{Nelements}: number of mesh elements
    ## @item @var{coeff}: piecewise linear reaction coefficient
    ## @end itemize
    ##
    ## @end deftypefn
    """
    h 	= nodes[1:len(nodes)]-nodes[0:len(nodes)-1]
    d0=np.zeros(n_max)
    d0[0]=coeff[0]/h[0]
    d0[n_max-1]=coeff[len(coeff)-1]/h[len(h)-1]
    d0[1:n_max-1]=(coeff[0:len(coeff)-1]/h[0:len(h)-1])+(coeff[1:len(coeff)]/h[1:len(h)])
    d1	= np.zeros(n_max)
    d1[0]=n_max
    d1[1:n_max]=-coeff/h
    dm1	= np.zeros(n_max)
    dm1[n_max-1]=n_max
    dm1[0:n_max-1]=-coeff/h
    L	= sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense()      
    return L

def Ucompconst (nodes,n_max,elements,Nelements,D,C):
    """
    ## -*- texinfo -*-
    ##
    ## @deftypefn {Function File}@
    ## {@var{R}} = Ucompconst(@var{nodes},@var{n_max},@var{elements},@var{Nelements},@var{D},@var{C})
    ##
    ## Compute P1 finite element rhs:
    ##
    ## @itemize @minus
    ## @item @var{nodes}: list of mesh nodes
    ## @item @var{n_max}: number of mesh nodes
    ## @item @var{elements}: list of mesh elements 
    ## @item @var{Nelements}: number of mesh elements
    ## @item @var{D}: piecewise linear reaction coefficient
    ## @item @var{C}: piecewise constant reaction coefficient
    ## @end itemize
    ##
    ## @end deftypefn
    """ 
    h = (nodes[1:len(nodes)]-nodes[0:len(nodes)-1])*C
    R=np.zeros(n_max)
    R[0]=D[0]*h[0]/2
    R[n_max-1]=D[n_max-1]*h[len(h)-1]/2
    R [1:n_max-1]=D[1:n_max-1]*(h[0:len(h)-1]+h[1:len(h)])/2
    return R


def Ubernoulli(x,sg):
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
    ## @deftypefn {Function File} {@var{b}} = Ubernoulli(@var{x},@var{sg})
    ##
    ## Compute Bernoulli function for vector x:
    ##
    ## @itemize @minus
    ## @item @var{b} = @var{x}/(exp(@var{x})-1) if @var{sg} == 1
    ## @item @var{b} = @var{x} + B( @var{x} ) if @var{sg} == 0
    ## @end itemize
    ##
    ## @end deftypefn
    """
    bernp=np.zeros(len(x))
    bernn=np.zeros(len(x))
    
    for count in range (len(x)):
        [bp,bn] = Ubern(x[count])
        
        bernp[count]=bp
        bernn[count]=bn
      
    if (sg ==1):
        b=bernp
    elif (sg ==0):
        b=bernn
    return b


def Ubern(x):
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
    ## @deftypefn {Function File} {@var{bp},@var{bn}} = Ubern(@var{x})
    ##
    ## Compute Bernoulli function for scalar x:
    ##
    ## @itemize @minus
    ## @item @var{bp} = @var{x}/(exp(@var{x})-1)
    ## @item @var{bn} = @var{x} + B( @var{x} )
    ## @end itemize
    ##
    ## @end deftypefn
    """      
    xlim=1e-2
    ax=abs(x)
    ## Compute Bernoulli function for x = 0
    if (ax == 0):
      bp=1.
      bn=1.
      return [bp,bn]    
    ## Compute Bernoulli function for asymptotic values    
    if (ax > 80):
        if (x > 0):
          bp=0.
          bn=x
          return [bp,bn]
        else:
          bp=-x
          bn=0.
          return [bp,bn]    
    ## Compute Bernoulli function for intermediate values    
    if (ax > xlim):
       bp=x/(exp(x)-1)
       bn=x+bp
       return [bp,bn]
    else:
       ## Compute Bernoulli function for small x
       ## via Taylor expansion    
       ii=1
       fp=1.
       fn=1.
       df=1.
       segno=1.
       while (abs(df) >np.finfo(np.float).eps):#eps on octave 2.220446049250313e-16
         ii=ii+1
         segno=-segno
         df=df*x/ii
         fp=fp+df
         fn=fn+segno*df
         bp=1/fp
         bn=1/fn
       return [bp,bn]    
    return [bp,bn]


def Uscharfettergummel(nodes,n_max,elements,Nelements,acoeff,bcoeff,v):
    
    """
    ## Copyright (C) 2004-2008  Carlo de Falco
      ##
      ## SECS1D - A 1-D Drift--Diffusion Semiconductor Device Simulator
      ##
    ## SECS1D is free software; you can redistribute it and/or modify
      ##  it under the terms of the GNU General Public License as published by
      ##  the Free Software Foundation; either version 2 of the License, or
      ##  (at your option) any later version.
      ##
    ## SECS1D is distributed in the hope that it will be useful,
      ##  but WITHOUT ANY WARRANTY; without even the implied warranty of
      ##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      ##  GNU General Public License for more details.
      ##
      ##  You should have received a copy of the GNU General Public License
    ## along with SECS1D; If not, see <http://www.gnu.org/licenses/>.
    ##
    ## author: Carlo de Falco <cdf _AT_ users.sourceforge.net>
    
    ## -*- texinfo -*-
    ##
    ## @deftypefn {Function File}@
    ## {@var{R}} = Uscharfettergummel(@var{nodes},@var{n_max},@var{elements},@var{Nelements},@var{acoeff},@var{bcoeff},@var{v})
    ##
    ## Build the Scharfetter-Gummel matrix for the the discretization of
    ## the LHS of the Drift-Diffusion equation:
    ##
    ##  -(a(x) (u' - b v'(x) u))'= f
    ##
    ## @itemize @minus
    ## @item @var{nodes}: list of mesh nodes
    ## @item @var{n_max}: number of mesh nodes
    ## @item @var{elements}: list of mesh elements 
    ## @item @var{Nelements}: number of mesh elements
    ## @item @var{acoeff}: piecewise linear diffusion coefficient
    ## @item @var{bcoeff}: piecewise constant drift constant coefficient
    ## @item @var{v}: piecewise linear drift potential
    ## @end itemize
    ##
    ## @end deftypefn
    """     
    h=nodes[1:n_max]-nodes[0:n_max-1]     
    c=acoeff[0:n_max-1]/h
    
    Bneg=Ubernoulli(-(v[1:n_max]-v[0:n_max-1])*bcoeff,1)
    Bpos=Ubernoulli( (v[1:n_max]-v[0:n_max-1])*bcoeff,1)
    
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
    A = sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense()
    return A


def Umediaarmonica(w):
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
    ## {@var{m}} = Umediarmonica(@var{w})
    ##
    ## Return the harmonic mean value of @var{w}
    ##
    ## @end deftypefn
    """     
    dw = (1/w[0:len(w)-1])+(1/w[1:len(w)])
    m  = 2 / dw
    return m
