import cython
cimport cython
import numpy as np
cimport numpy as np
from math import exp,log
from scipy import sparse as sp
from aestimo_poisson1d import equi_np_fi222
import config
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
DTYPE = np.float64
#
cdef double q = 1.602176e-19 #C
cdef double kb = 1.3806504e-23 #J/K
cdef double hbar = 1.054588757e-34
cdef double m_e= 9.1093826e-31 #kg
pi=np.pi
cdef double eps0= 8.8541878176e-12 #F/m
cdef double Cc=2.0
cdef double accurcy=0.004
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
ctypedef np.float64_t DTYPE_t
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
@cython.boundscheck(False) # stops bounds checking on numpy array lookups
#--------func_lib---------#
#use:
#python aestimo_dd_lib_setup.py build_ext --inplace
def DDGphin2n_cython (np.ndarray[np.double_t, ndim=1] V,phin,n):
    n=  np.exp ((V-phin))
    #n = n * (n>nmin) + nmin * (n<=nmin)
    return n
def DDGphip2p_cython (np.ndarray[np.double_t, ndim=1] V,phip,p):
    p=  np.exp ((phip-V))
    #p = p * (p>pmin) + pmin * (p<=pmin) 
    return p
def  DDGn2phin_cython (np.ndarray[np.double_t, ndim=1] V,n):
    cdef np.ndarray[np.double_t, ndim=1] phin
    #n    = n * (n>nmin) + nmin * (n<=nmin)
    phin = V - np.log(abs(n))    
    return phin  
def DDGp2phip_cython (np.ndarray[np.double_t, ndim=1] V,p):
    cdef np.ndarray[np.double_t, ndim=1] phip
    #p    = p * (p>pmin) + pmin * (p<=pmin)
    phip = V + np.log(abs(p))       
    return phip
def  Ucompmass_cython (np.ndarray[np.double_t, ndim=1] nodes,int n_max,np.ndarray[np.double_t, ndim=1] Bvect,Cvect):
    h 	= (nodes[1:len(nodes)]-nodes[0:len(nodes)-1])*Cvect  
    d0=np.zeros(n_max, dtype=np.double)
    d0[0]=Bvect[0]*h[0]/2
    d0[n_max-1]=Bvect[n_max-1]*h[len(h)-1]/2
    d0[1:n_max-1]=Bvect[1:n_max-1]*(h[0:len(h)-1]+h[1:len(h)])/2
    Bmat  = sp.spdiags(d0, [0], n_max,n_max).todense()      
    return Bmat
def Ucomplap_cython (np.ndarray[np.double_t, ndim=1] nodes,int n_max,np.ndarray[np.double_t, ndim=1] coeff):
    h 	= nodes[1:len(nodes)]-nodes[0:len(nodes)-1]
    d0=np.zeros(n_max, dtype=np.double)
    d0[0]=coeff[0]/h[0]
    d0[n_max-1]=coeff[len(coeff)-1]/h[len(h)-1]
    d0[1:n_max-1]=(coeff[0:len(coeff)-1]/h[0:len(h)-1])+(coeff[1:len(coeff)]/h[1:len(h)])
    d1	= np.zeros(n_max, dtype=np.double)
    d1[0]=n_max
    d1[1:n_max]=-coeff/h
    dm1	= np.zeros(n_max, dtype=np.double)
    dm1[n_max-1]=n_max
    dm1[0:n_max-1]=-coeff/h
    L	= sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense()      
    return L
def Ucompconst_cython (np.ndarray[np.double_t, ndim=1] nodes,int n_max,np.ndarray[np.double_t, ndim=1] D,C):
    h = (nodes[1:len(nodes)]-nodes[0:len(nodes)-1])*C
    R=np.zeros(n_max, dtype=np.double)
    R[0]=D[0]*h[0]/2
    R[n_max-1]=D[n_max-1]*h[len(h)-1]/2
    R [1:n_max-1]=D[1:n_max-1]*(h[0:len(h)-1]+h[1:len(h)])/2
    return R
def   Ubern_cython(double x):  
    xlim=1e-2
    ax=abs(x)
    if (ax == 0):
      bp=1.
      bn=1.
      return [bp,bn]    
    if (ax > 80):
        if (x > 0):
          bp=0.
          bn=x
          return [bp,bn]
        else:
          bp=-x
          bn=0.
          return [bp,bn]    
    if (ax > xlim):
       bp=x/(exp(x)-1)
       bn=x+bp
       return [bp,bn]
    else:  
       ii=1
       fp=1.
       fn=1.
       df=1.
       segno=1.
       while (abs(df) >np.finfo(np.float).eps):
         ii=ii+1
         segno=-segno
         df=df*x/ii
         fp=fp+df
         fn=fn+segno*df
         bp=1/fp
         bn=1/fn
       return [bp,bn]    
    return [bp,bn]
def Ubernoulli_cython(np.ndarray[np.double_t, ndim=1] x,double sg):
    bernp=np.zeros(len(x), dtype=np.double)
    bernn=np.zeros(len(x), dtype=np.double)    
    for count in range (len(x)):
        [bp,bn] = Ubern_cython(x[count])        
        bernp[count]=bp
        bernn[count]=bn      
    if (sg ==1):
        b=bernp
    elif (sg ==0):
        b=bernn
    return b


#--------poisson---------#

def  DDGnlpoisson_cython (model,idata,np.ndarray[DTYPE_t,ndim=1] xaxis,sinodes,Vin,nin,pin,Fnin,Fpin,dop,Ppz_Psp,ni,fi_e,fi_h,l2,double toll,Vt,int maxit,verbose):
    cdef np.ndarray[np.double_t, ndim=2] elements,Anew,A
    cdef np.ndarray[np.double_t, ndim=1] sielements,normr,fi_n,fi_p,dV,Rnew,R,Tv0,Mv,res,nodes
    cdef int dampit,n_max,BCnodes,newtit,dit,niter
    cdef double dampcoeff,tk,normrnew,reldVnorm,relresnorm,dVnorm
    dampit = 10
    dampcoeff	= 5
    nodes 	= xaxis
    sielements=np.zeros(len(sinodes)-1, dtype=np.double) 
    sielements[:] = sinodes[1:len(sinodes)]
    n_max	= len(nodes)
    fi_n=np.zeros(n_max, dtype=np.double)
    fi_p=np.zeros(n_max, dtype=np.double)
    #totdofs = n_max - Ndiricheletnodes
    elements=np.zeros((n_max-1,2), dtype=np.double)
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    Nelements=np.size(elements[:,0])
    BCnodes= n_max
    normr=np.zeros(maxit+1, dtype=np.double)
    V = Vin
    Fn = Fnin
    Fp = Fpin
    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,V,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,idata.n,idata.p)
    n = DDGphin2n_cython(V[sinodes]+fi_n[sinodes],Fn,idata.n)
    p = DDGphip2p_cython(V[sinodes]+fi_p[sinodes],Fp,idata.p)
    if (sinodes[0]==0):
        n[1]=nin[0]
        p[1]=pin[0]
    if (sinodes[n_max-1]==n_max-1):
        n[n_max-1]=nin[n_max-1]
        p[n_max-1]=pin[n_max-1]
    l22=l2*np.ones(Nelements)
    L      = Ucomplap_cython (nodes,n_max,l22)
    Mv            =  np.zeros(n_max, dtype=np.double)
    Mv[sinodes]   =  (n + p)
    Cv            =  np.ones(Nelements)
    M             =  Ucompmass_cython (nodes,n_max,Mv,Cv)
    Tv0            =  np.zeros(n_max, dtype=np.double)
    Tv0[sinodes]   = (n - p -dop-Ppz_Psp)
    Cv=  np.ones(Nelements)
    T0             =  Ucompconst_cython (nodes,n_max,Tv0,Cv)
    A=np.zeros((n_max,n_max), dtype=np.double)
    R=np.zeros(n_max, dtype=np.double)
    Anew=np.zeros((n_max,n_max), dtype=np.double)
    Rnew=np.zeros(n_max, dtype=np.double)
    A 		= L + M
    LV=np.dot(np.array(L) , V)
    R 		=  LV +T0
    A=np.delete(A, [0,BCnodes-1], 0)
    A=np.delete(A, [0,BCnodes-1], 1)
    R=np.delete(R, [0,BCnodes-1], 0)    
    normr[0]		=  np.linalg.norm(R,np.inf)
    relresnorm 	= 1
    reldVnorm   = 1
    normrnew	= normr[0]
    for newtit in range(1,maxit):
        ccc= np.linalg.solve(A, -R)
        dV=np.zeros(n_max, dtype=np.double)
        dV[1:n_max-1] =ccc
        tk = 1
        for dit in range(1,dampit):
            Vnew   = V + tk * dV
            if model.N_wells_virtual-2!=0 and config.quantum_effect:
                fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,Vnew,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,n,p)             
            n = DDGphin2n_cython(Vnew[sinodes]+fi_n[sinodes],Fn,idata.n)
            p = DDGphip2p_cython(Vnew[sinodes]+fi_p[sinodes],Fp,idata.p)
            if (sinodes[0]==0):
                n[0]=nin[0]
                p[0]=pin[0]
            if (sinodes[n_max-1]==n_max-1):
                n[n_max-1]=nin[n_max-1]
                p[n_max-1]=pin[n_max-1]
            Mv            =  np.zeros(n_max, dtype=np.double)
            Mv[sinodes]   =  (n + p)
            Cv            =  np.ones(Nelements)
            M    = Ucompmass_cython (nodes,n_max,Mv,Cv)
            Tv0            =  np.zeros(n_max, dtype=np.double)
            Tv0[sinodes]   =  (n - p -dop-Ppz_Psp)
            Cv=  np.ones(Nelements)
            Cv            =  np.ones(Nelements)
            T0     = Ucompconst_cython (nodes,n_max,Tv0,Cv)
            Anew 		= L + M
            LVnew=np.dot(np.array(L) , Vnew)
            Rnew 		=  LVnew +T0
            Anew=np.delete(Anew, [0,BCnodes-1], 0)
            Anew=np.delete(Anew, [0,BCnodes-1], 1)
            Rnew=np.delete(Rnew, [0,BCnodes-1], 0)            
            if ((dit>1) and (np.linalg.norm(Rnew,np.inf) >= np.linalg.norm(R,np.inf))):
                break
            else:
                A = Anew
                R = Rnew
            normrnew= np.linalg.norm(R,np.inf)
            if (normrnew > normr[newtit]):
                tk = tk/dampcoeff
            else:
                break
        V= Vnew	
        normr[newtit+1] = normrnew
        dVnorm= np.linalg.norm(tk*dV,np.inf)
        reldVnorm           = dVnorm / np.linalg.norm(V,np.inf)
        if (reldVnorm <= toll):
            break
    res = normr
    niter = newtit
    return [V,n,p]
#--------endpoisson---------#

def Uscharfettergummel_cython(np.ndarray[np.double_t, ndim=1] nodes, acoeff,v,int n_max,bcoeff):
    h=nodes[1:n_max]-nodes[0:n_max-1]     
    c=acoeff[0:n_max-1]/h    
    Bneg=Ubernoulli_cython(-(v[1:n_max]-v[0:n_max-1])*bcoeff,1)
    Bpos=Ubernoulli_cython( (v[1:n_max]-v[0:n_max-1])*bcoeff,1)    
    d0=np.zeros(n_max, dtype=np.double)
    d0[0]=c[0]*Bneg[0]
    d0[n_max-1]=c[len(c)-1]*Bpos[len(Bpos)-1]
    d0[1:n_max-1]=c[0:len(c)-1]*Bpos[0:len(Bpos)-1]+c[1:len(c)]*Bneg[1:len(Bneg)]    
    d1	= np.zeros(n_max, dtype=np.double)
    d1[0]=n_max
    d1[1:n_max]=-c* Bpos     
    dm1	= np.zeros(n_max, dtype=np.double)
    dm1[n_max-1]=n_max
    dm1[0:n_max-1]=-c* Bneg 
    A = sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense()
    return A     
def Umediaarmonica_cython(np.ndarray[np.double_t, ndim=1] w):   
    dw = (1/w[0:len(w)-1])+(1/w[1:len(w)])
    m  = 2 / dw
    return m
#--------endfunc_lib---------#

#--------gummelmap---------#

def DDGgummelmap_cython (idata,odata,model,np.ndarray[DTYPE_t,ndim=1] xaxis,ni,fi_e,fi_h,double toll,Vt,ptoll,pmaxit,int n_max,verbose,maxit):
    cdef int i,it
    odata  = idata
    dop = idata.dop
    Ppz_Psp= idata.Ppz_Psp
    cdef np.ndarray[np.double_t, ndim=2] vout,fermin,fermip,hole_density,electron_density
    cdef np.ndarray[np.double_t, ndim=1] nrm,res
    vout=np.zeros((n_max,2), dtype=np.double)
    hole_density=np.zeros((n_max,3), dtype=np.double)
    electron_density=np.zeros((n_max,3), dtype=np.double)
    fermin=np.zeros((n_max,2), dtype=np.double)
    fermip=np.zeros((n_max,2), dtype=np.double)
    nrm=np.zeros(maxit, dtype=np.double)
    res=np.zeros(maxit, dtype=np.double)
    v_Nnodes=np.arange(n_max)
    vout[:,0] = idata.V
    hole_density [:,0] = idata.p
    electron_density [:,0]= idata.n
    fermin [:,0]=idata.Fn
    fermip [:,0]=idata.Fp
    for i in range (1,maxit):
        [vout[:,1],electron_density[:,1],hole_density[:,1]] =DDGnlpoisson_cython (model,idata,xaxis,v_Nnodes,vout[:,0],electron_density[:,0],hole_density[:,0],fermin[:,0],fermip[:,0],dop,Ppz_Psp,ni,fi_e,fi_h,idata.l2,ptoll,Vt,pmaxit,verbose)
        electron_density[:,2]=DDGelectron_driftdiffusion_cython(model,idata,vout[:,1], xaxis, electron_density[:,1],hole_density[:,1],idata.nis,idata.TAUN0,idata.TAUP0,idata.mun,fi_e,fi_h,Vt)
        fermin[:,1] = DDGn2phin_cython(vout[:,1],electron_density[:,2])
        fermin[0,1]   = idata.Fn[0]
        fermin[n_max-1,1] = idata.Fn[len(idata.Fn)-1]
        hole_density[:,2] = DDGhole_driftdiffusion_cython(model,idata,vout[:,1], xaxis, hole_density[:,1],electron_density[:,1],idata.nis,idata.TAUN0,idata.TAUP0,idata.mup,fi_e,fi_h,Vt)
        fermip[:,1] = DDGp2phip_cython(vout[:,1],hole_density[:,2])
        fermip[0,1]   = idata.Fp[0]
        fermip[n_max-1,1] = idata.Fp[len(idata.Fp)-1]
        nrfn= np.linalg.norm(fermin[:,1]-fermin[:,0],np.inf)
        nrfp= np.linalg.norm (fermip[:,1]-fermip[:,0],np.inf)
        nrv = np.linalg.norm (vout[:,1]-vout[:,0],np.inf)
        nrm[i] = max([nrfn,nrfp,nrv])
        if (nrm[i]<toll):
            break
        vout[:,0] = vout[:,1]
        hole_density [:,0] = hole_density [:,2] 
        electron_density [:,0]= electron_density [:,2]
        fermin [:,0]= fermin [:,1]
        fermip [:,0]= fermip [:,1]
    it = i
    res = nrm        
    odata.n     = electron_density[:,2]
    odata.p     = hole_density[:,2]
    odata.V     = vout[:,1]
    odata.Fn    = fermin[:,1]
    odata.Fp    = fermip[:,1]
    return [odata,it,res]
#--------endgummelmap---------#

#--------newtonmap ---------#
def  DDNnewtonmap_cython (idata,model,np.ndarray[DTYPE_t,ndim=1] ni,fi_e,fi_h,xaxis,double toll,Vt,int maxit,verbose):
    import time
    cdef np.ndarray[np.double_t, ndim=2] elements,BCnodes_,Z
    cdef np.ndarray[np.double_t, ndim=1] fi_n,fi_p,res,u
    cdef int n_max,totaldofs,it,cc,maxdamp
    cdef double dampcoef,nrm_u,ratio,tj,nrmtmp,nrm_du,nrm_du_old
    odata     = idata
    n_max    = len(xaxis)
    fi_n=np.zeros(n_max, dtype=np.double)
    fi_p=np.zeros(n_max, dtype=np.double)
    Nelements=n_max-1
    elements=np.zeros((n_max-1,2), dtype=np.double)
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    BCnodesp = [0,n_max-1]
    BCnodesp1 = [n_max,2*n_max-1]
    BCnodesp2 = [2*n_max,3*n_max-1]
    BCnodes_=np.zeros((3,2), dtype=np.double)
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
    L = Ucomplap_cython (xaxis,n_max,idata.l2*np.ones(Nelements))
    M = Ucompmass_cython (xaxis,n_max,np.ones(n_max),np.ones(Nelements))
    DDn = Uscharfettergummel_cython(xaxis,idata.mun,V+fi_n,n_max,1)
    DDp = Uscharfettergummel_cython(xaxis,idata.mup,-V-fi_p,n_max,1)
    
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
    res=np.zeros(maxit, dtype=np.double)
    res[0] = nrm
    ## Begin Newton Cycle
    for count in range (0, maxit):
        if verbose:
          print ("Newton Iteration Number:%d\n"%count)	
        Ln = Ucomplap_cython (xaxis,n_max,Umediaarmonica_cython(idata.mun*n))
        Lp = Ucomplap_cython (xaxis,n_max,Umediaarmonica_cython(idata.mup*p))
        Mn = Ucompmass_cython (xaxis,n_max,np.ones(n_max),n[0:n_max-1]*fact[0:n_max-1])
        Mp = Ucompmass_cython (xaxis,n_max,np.ones(n_max),p[0:n_max-1]*fact[0:n_max-1])
        Z  = np.zeros((n_max,n_max), dtype=np.double)   
        DDn = Uscharfettergummel_cython(xaxis,idata.mun,V+fi_n,n_max,1)
        DDp = Uscharfettergummel_cython(xaxis,idata.mup,-V-fi_p,n_max,1)
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
          
          DDn = Uscharfettergummel_cython(xaxis,idata.mun,Vnew+fi_n,n_max,1)
          DDp = Uscharfettergummel_cython(xaxis,idata.mup,-Vnew-fi_p,n_max,1)
          
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

#--------endnewtonmap ---------#

#--------hole_driftdiffusion---------#

def  DDGhole_driftdiffusion_cython(model,idata,np.ndarray[np.double_t, ndim=1] psi,xaxis,pg,n,ni,TAUN0,TAUP0,mup,fi_e,fi_h,double Vt):

    cdef np.ndarray[np.double_t, ndim=2] elements
    cdef np.ndarray[np.double_t, ndim=1] fi_n,fi_p,p,b,d0,d1,dm1
    cdef np.ndarray[np.int_t, ndim=1] BCnodes_h
    cdef int n_max
    cdef double pl,pr
    nodes        = xaxis
    n_max     =len(nodes)
    fi_n=np.zeros(n_max, dtype=np.double)
    fi_p=np.zeros(n_max, dtype=np.double)
    elements=np.zeros((n_max-1,2), dtype=np.double)
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    Nelements=np.size(elements[:,0])
    BCnodes_h=np.zeros(n_max, dtype=np.int)
    BCnodes_h[0]=0
    BCnodes_h[1]=n_max-1
    
    pl = pg[0]
    pr = pg[n_max-1]
    h=nodes[1:len(nodes)]-nodes[0:len(nodes)-1]
    
    c=1/h
    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,psi,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,n,idata.p)
    Bneg=Ubernoulli_cython(-(psi[1:n_max]-psi[0:n_max-1])-(fi_n[1:n_max]-fi_n[0:n_max-1]),1)
    Bpos=Ubernoulli_cython( (psi[1:n_max]-psi[0:n_max-1])+(fi_p[1:n_max]-fi_p[0:n_max-1]),1)
    
    d0=np.zeros(n_max, dtype=np.double)
    d0[0]=c[0]*Bneg[0]
    d0[n_max-1]=c[len(c)-1]*Bpos[len(Bpos)-1]
    d0[1:n_max-1]=c[0:len(c)-1]*Bpos[0:len(Bpos)-1]+c[1:len(c)]*Bneg[1:len(Bneg)]    
    
    d1	= np.zeros(n_max, dtype=np.double)
    d1[0]=n_max
    d1[1:n_max]=-c* Bneg      
    dm1	= np.zeros(n_max, dtype=np.double)
    dm1[n_max-1]=n_max
    dm1[0:n_max-1]=-c* Bpos   
    A = sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense() 
    b = np.zeros(n_max, dtype=np.double)#%- A * ng
    
    ## SRH Recombination term
        
    SRHD = TAUP0 * (n + ni) + TAUN0 * (pg + ni)
    SRHL = n / SRHD
    SRHR = ni**2 / SRHD
    
    ASRH = Ucompmass_cython (nodes,n_max,SRHL,np.ones(Nelements))
    bSRH = Ucompconst_cython (nodes,n_max,SRHR,np.ones(Nelements))
    
    A = A + ASRH
    b = b + bSRH
    
    ## Boundary conditions
    b=np.delete(b, BCnodes_h, 0)
    b[0]         = - A[1,0] * pl
    b[len(b)-1]       = - A[n_max-2,n_max-1] * pr
    A=np.delete(A, BCnodes_h, 0)
    A=np.delete(A, BCnodes_h, 1)
    
    pp= np.linalg.solve(A,b)
    p=np.zeros(n_max, dtype=np.double)
    p[1:n_max-1]=pp    
    p[0]=pl
    p[len(p)-1]=pr
    return p
#--------endhole_driftdiffusion---------#


#--------electron_driftdiffusion---------#
def DDGelectron_driftdiffusion_cython(model,idata,np.ndarray[np.double_t, ndim=1] psi,xaxis,ng,p,ni,TAUN0,TAUP0,mun,fi_e,fi_h,double Vt):
    cdef np.ndarray[np.double_t, ndim=2] elements
    cdef np.ndarray[np.double_t, ndim=1] fi_n,fi_p,n,b,d0,d1,dm1
    cdef np.ndarray[np.int_t, ndim=1] BCnodes_e
    cdef int n_max
    cdef double nl,nr    
    nodes        = xaxis
    n_max     =len(nodes)
    fi_n=np.zeros(n_max, dtype=np.double)
    fi_p=np.zeros(n_max, dtype=np.double)    
    elements=np.zeros((n_max-1,2), dtype=np.double)
    elements[:,0]= np.arange(0,n_max-1)
    elements[:,1]=np.arange(1,n_max)
    Nelements=np.size(elements[:,0])
    BCnodes_e=np.zeros(n_max, dtype=np.int)
    BCnodes_e[0]=0
    BCnodes_e[1]=n_max-1
    
    nl = ng[0]
    nr = ng[n_max-1]
    h=nodes[1:n_max]-nodes[0:n_max-1]    
    c=1/h
    if model.N_wells_virtual-2!=0 and config.quantum_effect:
        fi_n,fi_p =equi_np_fi222(ni,idata,fi_e,fi_h,psi,Vt,idata.wfh_general,idata.wfe_general,model,idata.E_state_general,idata.E_statec_general,idata.meff_state_general,idata.meff_statec_general,n_max,idata.n,p)    
    Bneg=Ubernoulli_cython(-(psi[1:n_max]-psi[0:n_max-1])-(fi_n[1:n_max]-fi_n[0:n_max-1]),1)    
    Bpos=Ubernoulli_cython( (psi[1:n_max]-psi[0:n_max-1])+(fi_p[1:n_max]-fi_p[0:n_max-1]),1)
    d0=np.zeros(n_max, dtype=np.double)
    d0[0]=c[0]*Bneg[0]
    d0[n_max-1]=c[len(c)-1]*Bpos[len(Bpos)-1]
    d0[1:n_max-1]=c[0:len(c)-1]*Bpos[0:len(Bpos)-1]+c[1:len(c)]*Bneg[1:len(Bneg)]    
    
    d1	= np.zeros(n_max, dtype=np.double)
    d1[0]=n_max
    d1[1:n_max]=-c* Bpos      
    dm1	= np.zeros(n_max, dtype=np.double)
    dm1[n_max-1]=n_max
    dm1[0:n_max-1]=-c* Bneg
    A = sp.spdiags([dm1, d0, d1],np.array([-1,0,1]),n_max,n_max).todense()   
    b = np.zeros(n_max, dtype=np.double)#%- A * ng
 
    ## SRH Recombination term
    SRHD = TAUP0 * (ng + ni) + TAUN0 * (p + ni)
    SRHL = p / SRHD
    SRHR = ni**2 / SRHD
    
    ASRH = Ucompmass_cython (nodes,n_max,SRHL,np.ones(Nelements))
    bSRH = Ucompconst_cython (nodes,n_max,SRHR,np.ones(Nelements))  
    A = A + ASRH
    b = b + bSRH
    ## Boundary conditions
    b=np.delete(b, BCnodes_e, 0)
    b[0]         = - A[1,0] * nl
    b[len(b)-1]       =-A[n_max-2,n_max-1] * nr
    A=np.delete(A, BCnodes_e, 0)
    A=np.delete(A, BCnodes_e, 1)

    nn= np.linalg.solve(A, b)
    n=np.zeros(n_max, dtype=np.double)
    n[1:n_max-1]=nn
    n[0]=nl
    n[len(n)-1]=nr
    return n

#--------endelectron_driftdiffusion---------#



