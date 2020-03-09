# -*- coding: utf-8 -*-
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Copyright (C) 2013-2020 Aestimo Group

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

File Information:
-----------------
Created on Wed Aug  7 16:55:50 2019

"""
import time

time0 = time.time()  # timing audit
from math import *
import numpy as np
import constants
from constants import *
import DDGgummelmap
from DDGgummelmap import DDGgummelmap
import DDNnewtonmap
from DDNnewtonmap import DDNnewtonmap
import func_lib
from func_lib import Ubernoulli
import matplotlib.pyplot as pl

time0 = time.time()
len_ = 1e-6
n_max = 300

vmin = 0
vmax = 2

vstep = 0.15
Total_Steps = int((vmax - vmin) / vstep) + 1
istep = 0
va = vmin

xaxis = np.linspace(0, len_, n_max)
xm = np.mean(xaxis)

Nd = 1e25
Na = 1e25
dop = np.zeros(n_max)
TAUN0 = np.zeros(n_max)
TAUP0 = np.zeros(n_max)

nis = np.zeros(n_max)
mun = np.zeros(n_max)
mup = np.zeros(n_max)
l2 = np.zeros(n_max)
Fn = np.zeros(n_max)
Fp = np.zeros(n_max)
n = np.zeros(n_max)
p = np.zeros(n_max)
V = np.zeros(n_max)
vvect = np.zeros(Total_Steps)
n_ = np.zeros((Total_Steps, n_max))
p_ = np.zeros((Total_Steps, n_max))
Fn_ = np.zeros((Total_Steps, n_max))
Fp_ = np.zeros((Total_Steps, n_max))
V_ = np.zeros((Total_Steps, n_max))
Jn = np.zeros((Total_Steps, n_max - 1))
Jp = np.zeros((Total_Steps, n_max - 1))
J = np.zeros((Total_Steps, n_max - 1))
lambda2 = np.zeros(Total_Steps)
DV = np.zeros(Total_Steps)
Emax = np.zeros(Total_Steps)
for i in range(n_max):
    if xaxis[i] <= xm:
        dop[i] = -Na
    elif xaxis[i] > xm:
        dop[i] = Nd
nn = (Nd + sqrt(Nd ** 2 + 4 * ni_ ** 2)) / 2
pp = (Na + sqrt(Na ** 2 + 4 * ni_ ** 2)) / 2


xn = xm + 1e-7
xp = xm - 1e-7
## Scaling coefficients
xs = len_
ns = np.linalg.norm(dop, np.inf)
dop = dop / ns


class data:
    def __init__(self):
        self.dop = dop
        self.TAUN0 = TAUN0
        self.TAUP0 = TAUP0
        self.Fn = Fn
        self.Fp = Fp
        self.V = V
        self.n = n
        self.p = p
        self.nis = nis
        self.mun = mun
        self.mup = mup
        self.l2 = l2


idata = data()
odata = data()

Vs = Vt
us = mun0
Js = xs / (us * Vs * q * ns)
# print("Js=",Js,"Vt=",Vt,"nn=",nn,"pp=",pp,"dop[0]=",idata.dop[0],"dop[n_max-1]=",idata.dop[n_max-1])
# check_point_1
while va < vmax:
    #
    vvect[istep] = va
    #
    n_[istep, :] = nn * (xaxis >= xn) + (ni_ ** 2) / pp * (xaxis < xn)
    p_[istep, :] = (ni_ ** 2) / nn * (xaxis > xp) + pp * (xaxis <= xp)

    #
    Fn = va * (xaxis <= xm)
    Fp = Fn
    #
    V_[istep, :] = Fn - Vt * np.log(p_[istep, :] / ni_)
    #
    """
  print("n_=",n_[istep,:])
  print("p_=",p_[istep,:])
  print("V_=",V[istep,:])
  check_point_2
  """
    ## Scaling
    xin = xaxis / xs
    idata.n = n_[istep, :] / ns
    idata.p = p_[istep, :] / ns
    idata.V = V_[istep, :] / Vs
    idata.Fn = (Fn - Vs * log(ni_ / ns)) / Vs
    idata.Fp = (Fp + Vs * log(ni_ / ns)) / Vs
    #
    lambda2[istep] = idata.l2 = (Vs * eps) / (q * ns * xs ** 2)
    idata.nis = ni_ / ns
    idata.mun = mun0 / us
    idata.mup = mup0 / us
    #
    xbar = len_  # [m]
    Vbar = Vt  # [V]
    mubar = max(mun0, mup0)  # [m^2 V^{-1} s^{-1}]
    tbar = xbar ** 2 / (mubar * Vbar)  # [s]
    """
  print("n=",idata.n[:])
  print("p=",idata.p[:])
  print("V=",idata.V[:])
  check_point_3
  
  print("Fn=",idata.Fn[:])
  print("Fp=",idata.Fp[:])
  print("nis=",idata.nis)
  print("mun=",idata.mun)
  print("mup=",idata.mup)
  print("l2=",idata.l2)  
  check_point_4
  """
    ## Solution of DD system
    #
    ## Algorithm parameters
    toll = 1e-3
    maxit = 10
    ptoll = 1e-10
    pmaxit = 30
    verbose = 0
    sinodes = np.arange(len(xaxis))
    idata.TAUN0 = 1e-7 / tbar  # np.inf
    idata.TAUP0 = 1e-7 / tbar  # np.inf
    #
    """
  print("sinodes=",sinodes)
  print("TAUN0=",idata.TAUN0)  
  check_point_5
  """
    time00 = time.time()

    [odata, it, res] = DDGgummelmap(
        n_max, xin, idata, odata, toll, maxit, ptoll, pmaxit, verbose
    )

    """
  print("n=",odata.n[:])
  print("p=",odata.p[:])
  print("V=",odata.V[:])
  check_point_6  
  """
    time_DDGgummelmap = time.time()
    delta_time_DDGgummelmap = (time_DDGgummelmap - time00) / 60
    print("delta_time_DDGgummelmap=%fmn" % delta_time_DDGgummelmap)

    [odata, it, res] = DDNnewtonmap(xin, odata, toll, maxit, verbose)

    time_DDNnewtonmap = time.time()
    delta_time_DDNnewtonmap = (time_DDNnewtonmap - time_DDGgummelmap) / 60
    print("delta_time_DDNnewtonmap=%fmn" % delta_time_DDNnewtonmap)
    #
    n_[istep, :] = odata.n
    p_[istep, :] = odata.p
    V_[istep, :] = odata.V
    #
    Fn_[istep, :] = odata.Fn
    Fp_[istep, :] = odata.Fp
    DV[istep] = V_[istep, n_max - 1] - V_[0, istep]
    Emax[istep] = max(
        abs(
            (V_[istep, 1:n_max] - V_[istep, 0 : n_max - 1])
            / (xin[1:n_max] - xin[0 : n_max - 1])
        )
    )
    #
    Bp = Ubernoulli((V_[istep, 1:n_max] - V_[istep, 0 : n_max - 1]), 1)
    Bm = Ubernoulli((V_[istep, 1:n_max] - V_[istep, 0 : n_max - 1]), 0)
    Jn[istep, :] = (
        -odata.mun
        * (n_[istep, 1:n_max] * Bp - n_[istep, 0 : n_max - 1] * Bm)
        / (xin[1:n_max] - xin[0 : n_max - 1])
    )
    Jp[istep, :] = (
        odata.mup
        * (p_[istep, 1:n_max] * Bm - p_[istep, 0 : n_max - 1] * Bp)
        / (xin[1:n_max] - xin[0 : n_max - 1])
    )

    va = va + vstep
    istep = istep + 1
#
## Descaling
n_ = n_ * ns
p_ = p_ * ns
V_ = V_ * Vs
J = abs(Jp + Jn) * Js
Fn_ = Fn_ * Vs
Fp_ = Fp_ * Vs
#
time1 = time.time()
delta_t = (time1 - time0) / 60
print("time=%fmn" % delta_t)

pl.suptitle(
    "1D Drift Diffusion Model for pn Diodes Results - at Applied Bias (%.2f)V"
    % vvect[istep - 1],
    fontsize=12,
)
pl.subplots_adjust(hspace=0.4, wspace=0.4)
pl.subplot(2, 2, 1)
pl.plot(
    vvect[0 : istep - 1],
    abs(Jp[0 : istep - 1, 2] + Jn[0 : istep - 1, 2]) * us * q * ns,
    "b",
)
pl.xlabel("v [V]")
pl.ylabel("J [A]")
pl.title(" ", fontsize=10)
pl.grid(True)

pl.subplot(2, 2, 2)
pl.plot(xaxis * 1e6, n_[istep - 2, :])
pl.xlabel("xaxis(um)")
pl.ylabel("n_")
pl.title(" ", fontsize=10)
pl.grid(True)

pl.subplot(2, 2, 3)
pl.plot(xaxis * 1e6, p_[istep - 2, :])
pl.xlabel("xaxis(um)")
pl.ylabel("p_")
pl.title(" ", fontsize=10)
pl.grid(True)

pl.subplot(2, 2, 4)
pl.plot(xaxis * 1e6, Fn_[istep - 2, :], xaxis * 1e6, Fp_[istep - 2, :])
pl.xlabel("xaxis(um)")
pl.ylabel("Efn, Efp")
pl.title(" ", fontsize=10)
pl.grid(True)
pl.show()
