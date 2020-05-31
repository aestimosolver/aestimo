"""cythonised version of the psi_at_inf function. This function is the core of the 
shooting method and where the program spends the majority of its time, hence we
should be able to speed up aestimo considerably by compiling the function."""

hbar = 1.054588757e-34

#original function
def psi_at_inf_orig(E,fis,cb_meff,n_max,dx):
    """Shooting method for heterostructure as given in Harrison's book"""
    c0 = 2*(dx/hbar)**2
    # boundary conditions
    psi0 = 0.0                 
    psi1 = 1.0
    psi2 = None
    for j in xrange(1,n_max-1,1): # Last potential not used
        c1=2.0/(cb_meff[j]+cb_meff[j-1])
        c2=2.0/(cb_meff[j]+cb_meff[j+1])
        psi2=((c0*(fis[j]-E)+c2+c1)*psi1-c1*psi0)/c2
        psi0=psi1
        psi1=psi2
    return psi2

#cythonised version
def psi_at_inf(double E,list fis,list cb_meff,int n_max,dx):
    """Shooting method for heterostructure as given in Harrison's book.
    cythonised version of original psi_at_inf function."""
    #storing masses in temporary variables to limit number of list lookups
    cdef double cb_meff1 = cb_meff[0] 
    cdef double cb_meff2 = cb_meff[1]
    cdef double cb_meff3 = cb_meff[2]
    #
    cdef double fis1 = fis[1]
    #
    cdef double c0 = 2*(dx/hbar)**2
    cdef double c1,c2,ic2 # cython can guess this already
    # boundary conditions
    cdef double psi0 = 0.0                 
    cdef double psi1 = 1.0
    cdef double psi2 = 0.0
    cdef int j
    #for j from 1 <= j < n_max-1:
    #for j in xrange(1,n_max-1,1):
    for j in range(1,n_max-1,1): # Last potential not used
	# indexing the list in a separate step leads to a massive speed up of the code.
        cb_meff3 = cb_meff[j+1]
        fis1 = fis[j] 
        #
        c1=2.0/(cb_meff2+cb_meff1)
        c2=2.0/(cb_meff2+cb_meff3)
        ic2 = 0.5*(cb_meff2+cb_meff3) # avoids a division in the next line
        psi2=((c0*(fis1-E)+c2+c1)*psi1-c1*psi0)*ic2
        # preparing for the next iteration (shifting the values along by one for psi and cb_meff)
        psi0=psi1
        psi1=psi2
        cb_meff1 = cb_meff2
        cb_meff2 = cb_meff3
    return psi2
    
#cythonised version for numpy arrays
import cython
cimport cython
import numpy as np
cimport numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
DTYPE = np.float64
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
ctypedef np.float64_t DTYPE_t
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.

@cython.boundscheck(False) # stops bounds checking on numpy array lookups
def psi_at_inf_numpy(double E,np.ndarray[DTYPE_t,ndim=1] fis,np.ndarray[DTYPE_t,ndim=1] cb_meff,int n_max,dx):
    """Shooting method for heterostructure as given in Harrison's book.
    cythonised version of original psi_at_inf function."""
    assert fis.dtype == DTYPE and cb_meff.dtype == DTYPE
    #storing masses in temporary variables to limit number of list lookups
    cdef double cb_meff1 = cb_meff[0] 
    cdef double cb_meff2 = cb_meff[1]
    cdef double cb_meff3 = cb_meff[2]
    #
    cdef double fis1 = fis[1]
    #
    cdef double c0 = 2*(dx/hbar)**2
    cdef double c1,c2,ic2 # cython can guess this already
    # boundary conditions
    cdef double psi0 = 0.0                 
    cdef double psi1 = 1.0
    cdef double psi2 = 0.0
    cdef Py_ssize_t j # tells cython that indices will always be positive.
    #for j from 1 <= j < n_max-1:
    #for j in xrange(1,n_max-1,1):
    for j in range(1,n_max-1,1): # Last potential not used
	# indexing the list in a separate step leads to a massive speed up of the code.
        cb_meff3 = cb_meff[j+1]
        fis1 = fis[j] 
        #
        c1=2.0/(cb_meff2+cb_meff1)
        c2=2.0/(cb_meff2+cb_meff3)
        ic2 = 0.5*(cb_meff2+cb_meff3) # avoids a division in the next line
        psi2=((c0*(fis1-E)+c2+c1)*psi1-c1*psi0)*ic2
        # preparing for the next iteration (shifting the values along by one for psi and cb_meff)
        psi0=psi1
        psi1=psi2
        cb_meff1 = cb_meff2
        cb_meff2 = cb_meff3
    return psi2