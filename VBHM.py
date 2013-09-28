#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Copyright (C) 2013 Sefer Bora Lisesivdin and Aestimo group

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

 Description:   We consider the upper 3 × 3 Hamiltonian for a (001)-oriented
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
#from scipy.optimize import fsolve
import numpy as np
from math import *
#start
#material parametre for barrier and well
def qsv(GA1,GA2,GA3,RATIO,VNIT,ZETA,AC1,n_max,delta):
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

    return AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,FSO,GDELM
#
#define VB Hamiltonian

def VBMAT1(KP,AP1,AP2,AP3,AP4,AP5,AP6,FH,FL,FSO,GDELM,x_max,n_max,AC1,UNIM,KPINT,WB,BW):
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
            if I== BW+1 or I == WB+1 :
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

def VBMAT_V(B2,fi_h,RATIO,n_max,UNIM):
    #tmp=np.resize(fi_h,n_max+1)[1:]
    tmp = -RATIO*np.roll(fi_h,-1)
    #tmp[-1] = 0.0 #?
    A = tmp[:np.newaxis]*UNIM
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
            


