#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
We consider the upper 3 × 3 Hamiltonian for a (001)-oriented
zinc blende (ZB) crystal see ref [1], after block diagonalization 
the KP-FDM (k.P theory + Finite Difference Method ) method is explained
in ref [2] to ensure the Hermitian property of Hamiltonian you have to
apply the we have to write all operators of the form presented in ref [3],
to understand Hermiticity property see ref [2]. page 110,same code in fortran
language presented in the index of ref [1]. Dirichlet boundary conditions
were applied [1].
[1]:D.Ahn & S-H.Park 'ENGINEERING QUANTUM MECHANICS' P 238
[2]:P.Harrison 'QUANTUM WELLS, WIRES AND DOTS' P 357-362
[3]: S-L.CHUANG 'physics of Optoelectronic Devices' P 183
"""
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Copyright (C) 2013-2016 Sefer Bora Lisesivdin and Aestimo group

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. See ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt .

    For the list of contributors, see ~/AUTHORS
"""
#from scipy.optimize import fsolve
import numpy as np
from math import *
#start
#material parametre for barrier and well
def qsv(GA1,GA2,GA3,RATIO,VNIT,ZETA,CNIT,AC1,n_max,delta,A1,A2,A3,A4,A5,A6,delta_so,delta_cr,mat_type):
    AP1= np.zeros(n_max+2)
    AP2= np.zeros(n_max+2)
    AP3= np.zeros(n_max+2)
    AP4= np.zeros(n_max+2)
    AP5= np.zeros(n_max+2)
    AP6= np.zeros(n_max+2)
    GDELM= np.zeros(n_max+2)
    FH= np.zeros(n_max+2)
    FL= np.zeros(n_max+2)
    FSO=  np.zeros(n_max+2)
    DEL3=  np.zeros(n_max+2)
    DEL1=  np.zeros(n_max+2)
    DEL2=  np.zeros(n_max+2)
    Pce=  np.zeros(n_max+2)
    if mat_type=='Zincblende':     
        AP1=GA1
        AP2=GA2
        AP3=GA1
        AP4=GA2
        AP5=(GA2+GA3)/2.0
        AP6=GA3
        GDELM=RATIO*(sqrt(2.0)*ZETA)/AC1
        FH=RATIO*(VNIT+ZETA)/AC1
        FL=RATIO*(VNIT-ZETA)/AC1
        FSO=RATIO*(VNIT+delta)/AC1
        Pce=RATIO*(CNIT)/AC1
    if mat_type=='Wurtzite':
        AP1=A1
        AP2=A2
        AP3=A3
        AP4=A4
        AP5=A5
        AP6=A6
        FH=RATIO*(ZETA)/AC1
        FL=RATIO*(VNIT)/AC1
        DEL3=RATIO*(delta_so/3)/AC1
        DEL1=RATIO*(delta_cr)/AC1
        DEL2=RATIO*(delta_so/3)/AC1
        Pce=RATIO*(CNIT)/AC1+DEL1+DEL2
    AP1=np.resize(AP1,n_max+2)
    AP2=np.resize(AP2,n_max+2)
    AP3=np.resize(AP3,n_max+2)
    AP4=np.resize(AP4,n_max+2)
    AP5=np.resize(AP5,n_max+2)
    AP6=np.resize(AP6,n_max+2)
    FH=np.resize(FH,n_max+2)
    FL=np.resize(FL,n_max+2)
    FSO=np.resize(FSO,n_max+2)
    GDELM=np.resize(GDELM,n_max+2)
    DEL3=np.resize(DEL3,n_max+2)
    DEL1=np.resize(DEL1,n_max+2)
    DEL2=np.resize(DEL2,n_max+2)
    Pce=np.resize(Pce,n_max+2)
    return AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,FSO,Pce,GDELM,DEL3,DEL1,DEL2
#
#define VB Hamiltonian

def VBMAT1(KP,AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,FSO,GDELM,x_max,n_max,AC1,UNIM,KPINT):
    AVH1= np.zeros(n_max+2)
    AVH2= np.zeros(n_max+2)
    BVH1= np.zeros(n_max+2)
    BVH2= np.zeros(n_max+2)
    GVH1= np.zeros(n_max+2)	
    GVH2= np.zeros(n_max+2)		
    CVH= np.zeros(n_max+2)
    DVH= np.zeros(n_max+2)		 
    FVH1= np.zeros(n_max+2)
    FVH2= np.zeros(n_max+2)
    AVHDF2=np.zeros((n_max+2, n_max+2))
    BVHDF2=np.zeros((n_max+2, n_max+2))
    FVHDF2=np.zeros((n_max+2, n_max+2))
    GVHDF2=np.zeros((n_max+2, n_max+2))
    DVHDF1=np.zeros((n_max+2, n_max+2))
    B=np.zeros((n_max*3, n_max*3)) 
#=====DEFINE MATRIX=================
    for I in range(0,n_max+2,1):
            AVH1[I]=0.5*(AP3[I]+AP4[I])*((x_max*KP*KPINT)**2)
            AVH2[I]=0.5*(AP1[I]-2.0*AP2[I])
            BVH1[I]=0.5*(AP3[I]-AP4[I])*((x_max*KP*KPINT)**2)
            BVH2[I]=0.5*(AP1[I]+2.0*AP2[I])
            GVH1[I]=0.5*AP3[I]*((x_max*KP*KPINT)**2)	
            GVH2[I]=0.5*AP1[I]			
            CVH[I]=0.5*sqrt(3.0)*AP5[I]*((x_max*KP*KPINT)**2) 
            DVH[I]=sqrt(3.0)*AP6[I]*x_max*KP*KPINT			 
            FVH1[I]=AP4[I]/sqrt(2.0)*((x_max*KP*KPINT)**2) 
            FVH2[I]=sqrt(2.0)*AP2[I]
            #boundary condition between barrier and well 
    for I in range (0,n_max,1):
                AVH2[I]=(AVH2[I-1]+AVH2[I+1])/2.0
                BVH2[I]=(BVH2[I-1]+BVH2[I+1])/2.0
                GVH2[I]=(GVH2[I-1]+GVH2[I+1])/2.0
                FVH2[I]=(FVH2[I-1]+FVH2[I+1])/2.0
                DVH[I]=(DVH[I-1]+DVH[I+1])/2.0
                CVH[I]=(CVH[I-1]+CVH[I+1])/2.0
                FVH1[I]=(FVH1[I-1]+FVH1[I+1])/2.0
                BVH1[I]=(BVH1[I-1]+BVH1[I+1])/2.0
                GVH1[I]=(GVH1[I-1]+GVH1[I+1])/2.0
                AVH1[I]=(AVH1[I-1]+AVH1[I+1])/2.0  
    for I in range (0,n_max,1):
            # filling i row
            AVHDF2[I,I]+=-AVH2[I+1]-(AVH2[I]+AVH2[I+2])/2.0
            BVHDF2[I,I]+=-BVH2[I+1]-(BVH2[I]+BVH2[I+2])/2.0
            FVHDF2[I,I]+=-FVH2[I+1]-(FVH2[I]+FVH2[I+2])/2.0
            GVHDF2[I,I]+=-GVH2[I+1]-(GVH2[I]+GVH2[I+2])/2.0
            # filling i-1 and i+1 rows
            AVHDF2[I,I+1]+=(AVH2[I+1]+AVH2[I+2])/2.0
            AVHDF2[I+1,I]+=(AVH2[I+2]+AVH2[I+1])/2.0
            BVHDF2[I,I+1]+=(BVH2[I+1]+BVH2[I+2])/2.0
            BVHDF2[I+1,I]+=(BVH2[I+2]+BVH2[I+1])/2.0
            FVHDF2[I,I+1]+=(FVH2[I+1]+FVH2[I+2])/2.0
            FVHDF2[I+1,I]+=(FVH2[I+2]+FVH2[I+1])/2.0
            GVHDF2[I,I+1]+=(GVH2[I+1]+GVH2[I+2])/2.0
            GVHDF2[I+1,I]+=(GVH2[I+2]+GVH2[I+1])/2.0
            DVHDF1[I,I+1]+=(DVH[I+2]+DVH[I+1])/2.0
            DVHDF1[I+1,I]+=-(DVH[I+2]+DVH[I+1])/2.0
    for I in range (0,n_max,1):
            for J in range (0,n_max,1):        
                AVHDF2[I,J]*=AC1
                BVHDF2[I,J]*=AC1
                FVHDF2[I,J]*=AC1
                GVHDF2[I,J]*=AC1
                DVHDF1[I,J]*=sqrt(AC1)/2.0
                #here where start filling the matrix
    for I in range (0,n_max,1):
        for J in range (0,n_max,1):
                B[I,J]+=AVH1[I+1]*UNIM[I,J]-AVHDF2[I,J]+FH[I+1]*UNIM[I,J]*AC1
                B[I+n_max,J+n_max]+=BVH1[I+1]*UNIM[I,J]-BVHDF2[I,J]+FL[I+1]*UNIM[I,J]*AC1
                B[I+n_max*2,J+n_max*2]+=GVH1[I+1]*UNIM[I,J]-GVHDF2[I,J]+FSO[I+1]*UNIM[I,J]*AC1
                B[I+n_max,J]+=CVH[I+1]*UNIM[I,J]+DVHDF1[I,J]
                B[I+n_max*2,J]+=sqrt(2.0)*CVH[I+1]*UNIM[I,J]-DVHDF1[I,J]/sqrt(2.0)
                B[I+n_max*2,J+n_max]+=FVH1[I+1]*UNIM[I,J]+FVHDF2[I,J]-sqrt(3.0/2.0)*DVHDF1[I,J]+GDELM[I+1]*UNIM[I,J]*AC1
                #filling the symitric part of the hamiltonian matrix
                B[I,J+n_max]=B[J+n_max,I]
                B[I,J+n_max*2]=B[J+n_max*2,I]
                B[I+n_max,J+n_max*2]=B[J+n_max*2,I+n_max]  
    return B

def VBMAT2(KP,AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,x_max,n_max,AC1,UNIM,KPINT,DEL3,DEL1,DEL2):
    B11V1= np.zeros(n_max+2)
    B11V2= np.zeros(n_max+2)
    B12V1= np.zeros(n_max+2)
    B13V1= np.zeros(n_max+2)	
    B21V1= np.zeros(n_max+2)		
    B22V1= np.zeros(n_max+2)
    B22V2= np.zeros(n_max+2)		 
    B23V1= np.zeros(n_max+2)
    B31V1= np.zeros(n_max+2)
    B32V1= np.zeros(n_max+2)		 
    B33V1= np.zeros(n_max+2)
    B33V2= np.zeros(n_max+2)
    B11=np.zeros((n_max+2, n_max+2))
    B22=np.zeros((n_max+2, n_max+2))
    B33=np.zeros((n_max+2, n_max+2))
    B13=np.zeros((n_max+2, n_max+2))
    B23=np.zeros((n_max+2, n_max+2))
    B31=np.zeros((n_max+2, n_max+2))
    B32=np.zeros((n_max+2, n_max+2))
    B=np.zeros((n_max*3, n_max*3)) 
#=====DEFINE MATRIX=================
    for I in range(0,n_max+2,1):
            B11V2[I]=0.5*(AP2[I]+AP4[I])*((KP*x_max*KPINT)**2)
            B12V1[I]=0.5*AP5[I]*((KP*x_max*KPINT)**2)
            B21V1[I]=0.5*AP5[I]*((KP*x_max*KPINT)**2)
            B22V2[I]=0.5*(AP2[I]+AP4[I])*((KP*x_max*KPINT)**2)
            B33V2[I]=0.5*AP2[I]*(KP*x_max*KPINT)**2
            ###################################################
            B33V1[I]=-0.5*AP1[I]
            B22V1[I]=-0.5*(AP1[I]+AP3[I])
            B11V1[I]=-0.5*(AP1[I]+AP3[I])
            ###################################################
            B13V1[I]=0.5*AP6[I]*(KP*x_max*KPINT)
            B23V1[I]=0.5*AP6[I]*(KP*x_max*KPINT)
            B31V1[I]=0.5*AP6[I]*(KP*x_max*KPINT)
            B32V1[I]=0.5*AP6[I]*(KP*x_max*KPINT)
            #boundary condition between barrier and well 
    for I in range (0,n_max,1):
                B33V2[I]=(B33V2[I-1]+B33V2[I+1])/2.0
                B33V1[I]=(B33V1[I-1]+B33V1[I+1])/2.0
                B32V1[I]=(B32V1[I-1]+B32V1[I+1])/2.0
                B31V1[I]=(B31V1[I-1]+B31V1[I+1])/2.0
                B23V1[I]=(B23V1[I-1]+B23V1[I+1])/2.0
                B22V2[I]=(B22V2[I-1]+B22V2[I+1])/2.0
                B22V1[I]=(B22V1[I-1]+B22V1[I+1])/2.0
                B21V1[I]=(B21V1[I-1]+B21V1[I+1])/2.0
                B13V1[I]=(B13V1[I-1]+B13V1[I+1])/2.0
                B12V1[I]=(B12V1[I-1]+B12V1[I+1])/2.0
                B11V2[I]=(B11V2[I-1]+B11V2[I+1])/2.0
                B11V1[I]=(B11V1[I-1]+B11V1[I+1])/2.0                
    for I in range (0,n_max,1):
            # filling i row
            B11[I,I]+=-B11V1[I+1]-(B11V1[I]+B11V1[I+2])/2.0
            B22[I,I]+=-B22V1[I+1]-(B22V1[I]+B22V1[I+2])/2.0
            B33[I,I]+=-B33V1[I+1]-(B33V1[I]+B33V1[I+2])/2.0
            # filling i-1 and i+1 rows
            B11[I,I+1]+=(B11V1[I+1]+B11V1[I+2])/2.0
            B11[I+1,I]+=(B11V1[I+2]+B11V1[I+1])/2.0
            B13[I,I+1]+=-(B13V1[I+1]+B13V1[I+2])/2.0
            B13[I+1,I]+=(B13V1[I+1]+B13V1[I+1])/2.0
            B22[I,I+1]+=(B22V1[I+1]+B22V1[I+2])/2.0
            B22[I+1,I]+=(B22V1[I+1]+B22V1[I+1])/2.0
            B23[I,I+1]+=-(B23V1[I+1]+B23V1[I+2])/2.0
            B23[I+1,I]+=(B23V1[I+1]+B23V1[I+1])/2.0
            B31[I,I+1]+=(B31V1[I+2]+B31V1[I+1])/2.0
            B31[I+1,I]+=-(B31V1[I+1]+B31V1[I+1])/2.0
            B32[I,I+1]+=(B32V1[I+1]+B32V1[I+2])/2.0
            B32[I+1,I]+=-(B32V1[I+1]+B32V1[I+1])/2.0
            B33[I,I+1]+=(B33V1[I+1]+B33V1[I+2])/2.0
            B33[I+1,I]+=(B33V1[I+2]+B33V1[I+1])/2.0           

    for I in range (0,n_max,1):
            for J in range (0,n_max,1):        
                B11[I,J]*=AC1
                B13[I,J]*=sqrt(AC1)/2
                B22[I,J]*=AC1
                B23[I,J]*=sqrt(AC1)/2
                B31[I,J]*=sqrt(AC1)/2
                B32[I,J]*=sqrt(AC1)/2
                B33[I,J]*=AC1
                #here where start filling the matrix
    for I in range (0,n_max,1):
        for J in range (0,n_max,1):
                B[I,J]+=B11[I,J]+B11V2[I+1]*UNIM[I,J]+(DEL1[I+1]+DEL2[I+1]+FH[I+1]+FL[I+1])*UNIM[I,J]*AC1
                B[I+n_max,J+n_max]+=B22[I,J]+B22V2[I+1]*UNIM[I,J]+(DEL1[I+1]-DEL2[I+1]+FH[I+1]+FL[I+1])*UNIM[I,J]*AC1
                B[I+n_max*2,J+n_max*2]+=B33[I,J]+B33V2[I+1]*UNIM[I,J]+FH[I+1]*UNIM[I,J]*AC1
                B[I+n_max,J]+=B21V1[I+1]*UNIM[I,J]
                B[I+n_max*2,J]+=B31[I,J]
                B[I+n_max*2,J+n_max]+=B32[I,J]+sqrt(2)*DEL3[I+1]*UNIM[I,J]*AC1       
                B[I,J+n_max]+=B12V1[I+1]*UNIM[I,J]
                B[I,J+n_max*2]+=B13[I,J]
                B[I+n_max,J+n_max*2]+=B23[I,J]+sqrt(2)*DEL3[I+1]*UNIM[I,J]*AC1
    return B
def VBMAT_V(B2,fi_h,RATIO,n_max,UNIM):
    """B2 - float?
    fi_h - 1d array, potential
    RATIO - float?
    n_max - int, number of points in grid.
    UNIM - 2d array
    """
    #tmp=np.resize(fi_h,n_max+1)[1:]
    tmp = -RATIO*np.roll(fi_h,-1)
    #tmp[-1] = 0.0 #?
    A = tmp[:,np.newaxis]*UNIM
    Z = np.zeros((n_max,n_max))
    #Creating HUPMAT4
    #1
    #HUPMAT4 = np.vstack((np.hstack((A,Z,Z)),
    #                     np.hstack((Z,A,Z)),
    #                     np.hstack((Z,Z,A))))
    #2 - using the strange c_ and r_ 'functions'
    #HUPMAT4 = np.r_[np.c_[A,Z,Z],
    #                np.c_[Z,A,Z],
    #                np.c_[Z,Z,A]]
    #3 - using block matrices (skips intermediates)
    HUPMAT4 = np.asarray(np.bmat([(A,Z,Z),(Z,A,Z),(Z,Z,A)]))
    
    HUPMAT4+=B2       
    return HUPMAT4

def VBMAT_V_2(B2,fi_h,RATIO,i_1,I1,UNIM):
    """B2 - float?
    fi_h - 1d array, potential
    RATIO - float?
    n_max - int, number of points in grid.
    UNIM - 2d array
    """
    I2=I1+i_1
    #tmp=np.resize(fi_h,n_max+1)[1:]
    tmp = -RATIO*np.roll(fi_h[I1:I2],-1)
    #tmp[-1] = 0.0 #?
    A = tmp[:,np.newaxis]*UNIM[I1:I2,I1:I2]
    Z = np.zeros((i_1,i_1))
    #Creating HUPMAT4
    #1
    #HUPMAT4 = np.vstack((np.hstack((A,Z,Z)),
    #                     np.hstack((Z,A,Z)),
    #                     np.hstack((Z,Z,A))))
    #2 - using the strange c_ and r_ 'functions'
    #HUPMAT4 = np.r_[np.c_[A,Z,Z],
    #                np.c_[Z,A,Z],
    #                np.c_[Z,Z,A]]
    #3 - using block matrices (skips intermediates)
    HUPMAT3_general_2 = np.asarray(np.bmat([(A,Z,Z),(Z,A,Z),(Z,Z,A)]))
    
    HUPMAT3_general_2+=B2       
    return HUPMAT3_general_2
 
def VBMAT_V_old(B2,fi_h,RATIO,n_max,UNIM):
    HUPMAT4=np.zeros((n_max*3, n_max*3))
    tmp1=np.zeros(n_max)
    tmp1=fi_h
    tmp1=np.resize(tmp1,n_max+2)    
    tmp=RATIO*(-tmp1)
    tmp=tmp.tolist()
    for I in range (0,n_max,1):
        for J in range (0,n_max,1):
            HUPMAT4[I,J]+=tmp[I+1]*UNIM[I,J]
            HUPMAT4[I+n_max,J+n_max]+=tmp[I+1]*UNIM[I,J]
            HUPMAT4[I+n_max*2,J+n_max*2]+=tmp[I+1]*UNIM[I,J]
    HUPMAT4+=B2       
    return HUPMAT4

def CBMAT(KP,Pce,EM,x_max,n_max,AC1,UNIM,KPINT):
    B11V1= np.zeros(n_max+2)
    B11V2= np.zeros(n_max+2)
    B11=np.zeros((n_max+2, n_max+2))
    B=np.zeros((n_max, n_max))
    EM=np.resize(EM,n_max+2)
#=====DEFINE MATRIX=================
    for I in range(0,n_max+2,1):
            B11V2[I]=0.5*((KP*x_max*KPINT)**2)/EM[I]
            B11V1[I]=-0.5/EM[I]
            #boundary condition between barrier and well
    for I in range (0,n_max,1):
                B11V2[I]=(B11V2[I-1]+B11V2[I+1])/2.0
                B11V1[I]=(B11V1[I-1]+B11V1[I+1])/2.0
    for I in range (0,n_max,1):
            # filling i row
            B11[I,I]+=-B11V1[I+1]-(B11V1[I]+B11V1[I+2])/2.0
            # filling i-1 and i+1 rows
            B11[I,I+1]+=(B11V1[I+1]+B11V1[I+2])/2.0
            B11[I+1,I]+=(B11V1[I+2]+B11V1[I+1])/2.0

    for I in range (0,n_max,1):
            for J in range (0,n_max,1):
                B11[I,J]*=AC1
                #here where start filling the matrix
    #UNIM=np.resize(UNIM,n_max+2)
    for I in range (0,n_max,1):
        for J in range (0,n_max,1):
                B[I,J]+=B11[I,J]+B11V2[I+1]*UNIM[I,J]+Pce[I+1]*UNIM[I,J]*AC1
    return B

def CBMAT_V(BC,fi,RATIO,n_max,UNIM):
    #tmp=np.resize(fi_h,n_max+1)[1:]
    tmp = RATIO*np.roll(fi,-1)
    #tmp[-1] = 0.0 #?
    HUPMAT7 = tmp[:,np.newaxis]*UNIM    
    HUPMAT7+=BC       
    return HUPMAT7

def CBMAT_V_old(BC,fi,RATIO,n_max,UNIM):
    HUPMAT7=np.zeros((n_max, n_max))
    tmp1=np.zeros(n_max)
    tmp=np.zeros(n_max+2)
    tmp1=fi
    tmp1=np.resize(tmp1,n_max+2)
    tmp=RATIO*(tmp1)
    tmp=tmp.tolist()
    for I in range (0,n_max,1):
        for J in range (0,n_max,1):
            HUPMAT7[I,J]+=tmp[I+1]*UNIM[I,J]
    HUPMAT7+=BC
    return HUPMAT7

