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
Most of the code in this file are converted from Matlab to python and adapted to fit aestimo's
needs from this book:
Computational electronics : semiclassical and quantum device modeling and simulation. by:
    [Dragica Vasileska; Stephen M Goodnick; Gerhard Klimeck]
"""
import numpy as np
import matplotlib.pyplot as pl
from math import exp, log, sqrt

if __package__:  # explicit relative imports for using aestimo as a package (in python3)
    from . import config
else:
    import config

# Defining constants and material parameters
q = 1.602176e-19  # C
kb = 1.3806504e-23  # J/K
nii = 0.0
hbar = 1.054588757e-34  # m2 kg / s
m_e = 9.1093826e-31  # kg
pi = np.pi
eps0 = 8.8541878176e-12  # F/m
# TEMPERATURE
T = 300.0  # Kelvin
Vt = kb * T / q  # [eV]
J2meV = 1e3 / q  # Joules to meV
meV2J = 1e-3 * q  # meV to Joules


def Poisson_equi2(
    ns,
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    n,
    p,
    dx,
    Ldi,
    dop,
    Ppz_Psp,
    pol_surf_char,
    ni,
    n_max,
    iteration,
    fi_old,
    Vt,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    surface,
    fi_stat,
):
    fi_out = np.zeros(n_max)
    dop_out = np.zeros(n_max)
    Ppz_Psp_out = np.zeros(n_max)
    d = np.zeros(n_max)
    v = np.zeros(n_max)
    f = np.zeros(n_max)
    c = np.zeros(n_max)
    b = np.zeros(n_max)
    a = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = fi_old
    dop_out = dop / (ns * ni)
    Ppz_Psp_out = Ppz_Psp / (ns * ni)
    delta_acc = 1.0e-5
    dx2 = dx * dx
    Ldi2 = Ldi * Ldi
    iteration0 = 1

    if iteration == iteration0:
        # Determination of the Fermi Level
        EF = 0.0
        fi_out = np.zeros(n_max)
        n, p, fi_out = equi_np_fi(
            iteration, dop, Ppz_Psp, n_max, ni, model, Vt, surface
        )
        fi_stat = fi_out
    else:
        n, p, fi_non, EF = equi_np_fi3(
            fi_out,
            wfh_general,
            wfe_general,
            model,
            E_state_general,
            E_statec_general,
            meff_state_general,
            meff_statec_general,
            n_max,
            ni,
        )
    for i in range(1, n_max - 1):
        a[i] = Ldi2[i] / (dx2)
        c[i] = Ldi2[i] / (dx2)
        b[i] = -(2 * Ldi2[i] / (dx2) + n[i] / ns + p[i] / ns)
        f[i] = (
            n[i] / ns
            - p[i] / ns
            - dop_out[i]
            - Ppz_Psp_out[i]
            - (fi_out[i] * (n[i] / ns + p[i] / ns))
        )
    # (B) Define the elements of the coefficient matrix and initialize the forcing
    # function at the ohmic contacts
    a[0] = 0.0
    c[0] = 0.0
    b[0] = 1.0
    f[0] = fi_out[0]
    a[n_max - 1] = 0.0
    c[n_max - 1] = 0.0
    b[n_max - 1] = 1.0
    f[n_max - 1] = fi_out[n_max - 1]
    # (C) Start the iterative procedure for the solution of the linearized Poisson
    # equation using LU decomposition method:
    flag_conv = True  # convergence of the Poisson loop
    k_iter = 0
    while flag_conv:
        k_iter = k_iter + 1
        d[0] = b[0]
        for i in range(1, n_max):
            d[i] = b[i] - a[i] * c[i - 1] / d[i - 1]
        # Solution of Lv=f:
        v[0] = f[0]
        for i in range(1, n_max):
            v[i] = f[i] - a[i] * v[i - 1] / d[i - 1]
        # Solution of U*fi=v:
        temp = v[n_max - 1] / d[n_max - 1]
        delta[n_max - 1] = temp - fi_out[n_max - 1]
        fi_out[n_max - 1] = temp
        for i in range(n_max - 2, -1, -1):
            temp = (v[i] - c[i] * fi_out[i + 1]) / d[i]
            delta[i] = temp - fi_out[i]
            fi_out[i] = temp
        # Test update in the outer iteration loop:
        delta_max = 0.0
        delta_max = max(abs(delta[:]))
        # print ('k_iter=',k_iter, 'delta_max=',delta_max)
        # Test convergence and recalculate forcing function and
        # central coefficient b if necessary:
        if delta_max < delta_acc:
            flag_conv = False
        else:
            if iteration == iteration0:
                n = np.exp(fi_out)
                p = np.exp(-fi_out)
            else:
                n, p, fi_non, EF = equi_np_fi3(
                    fi_out,
                    wfh_general,
                    wfe_general,
                    model,
                    E_state_general,
                    E_statec_general,
                    meff_state_general,
                    meff_statec_general,
                    n_max,
                    ni,
                )
            for i in range(1, n_max - 1):
                b[i] = -(2 * Ldi2[i] / (dx2) + n[i] / ns + p[i] / ns)
                f[i] = (
                    n[i] / ns
                    - p[i] / ns
                    - dop_out[i]
                    - Ppz_Psp_out[i]
                    - (fi_out[i] * (n[i] / ns + p[i] / ns))
                )
    return n, p, fi_out, EF, fi_stat


def Poisson_equi_non_2(
    vindex,
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    n,
    p,
    dx,
    Ldi,
    dop,
    Ppz_Psp,
    pol_surf_char,
    ni,
    n_max,
    iteration,
    fi_old,
    Vt,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    surface,
    fi_stat,
):
    fi_out = np.zeros(n_max)
    dop_out = np.zeros(n_max)
    Ppz_Psp_out = np.zeros(n_max)
    d = np.zeros(n_max)
    v = np.zeros(n_max)
    f = np.zeros(n_max)
    c = np.zeros(n_max)
    b = np.zeros(n_max)
    a = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = fi_old
    dop_out = dop / ni
    Ppz_Psp_out = Ppz_Psp / ni
    delta_acc = 1.0e-5
    dx2 = dx * dx
    Ldi2 = Ldi * Ldi
    if config.predic_correc:  # (check high concentrations)
        n = np.zeros(n_max)
        p = np.zeros(n_max)
        n, p, fi_non, EF = equi_np_fi3(
            fi_out,
            wfh_general,
            wfe_general,
            model,
            E_state_general,
            E_statec_general,
            meff_state_general,
            meff_statec_general,
            n_max,
            ni,
        )
    else:
        n, p, fi_non = equi_np_fi22(
            vindex,
            fitotc,
            fitot,
            Nc,
            Nv,
            fi_e,
            fi_h,
            iteration,
            fi_out,
            Vt,
            wfh_general,
            wfe_general,
            model,
            E_state_general,
            E_statec_general,
            meff_state_general,
            meff_statec_general,
            dop,
            Ppz_Psp,
            n_max,
            ni,
            n,
            p,
        )
    for i in range(0, n_max):
        a[i] = Ldi2[i] / (dx2)
        c[i] = Ldi2[i] / (dx2)
        b[i] = -(2 * Ldi2[i] / (dx2) + n[i] + p[i])
        f[i] = n[i] - p[i] - dop_out[i] - Ppz_Psp_out[i] - fi_out[i] * (n[i] + p[i])
    # (B) Define the elements of the coefficient matrix and initialize the forcing
    # function at the ohmic contacts
    a[0] = 0.0
    c[0] = 0.0
    b[0] = 1.0
    f[0] = fi_out[0]
    a[n_max - 1] = 0.0
    c[n_max - 1] = 0.0
    b[n_max - 1] = 1.0
    f[n_max - 1] = fi_out[n_max - 1]
    # (C) Start the iterative procedure for the solution of the linearized Poisson
    # equation using LU decomposition method:
    flag_conv = True  # convergence of the Poisson loop
    k_iter = 0
    while flag_conv:
        k_iter = k_iter + 1
        d[0] = b[0]
        for i in range(1, n_max):
            d[i] = b[i] - a[i] * c[i - 1] / d[i - 1]
        # Solution of Lv=f:
        v[0] = f[0]
        for i in range(1, n_max):
            v[i] = f[i] - a[i] * v[i - 1] / d[i - 1]
        # Solution of U*fi=v:
        temp = v[n_max - 1] / d[n_max - 1]
        delta[n_max - 1] = temp - fi_out[n_max - 1]
        fi_out[n_max - 1] = temp
        for i in range(n_max - 2, -1, -1):
            temp = (v[i] - c[i] * fi_out[i + 1]) / d[i]
            delta[i] = temp - fi_out[i]
            fi_out[i] = temp
        # Test update in the outer iteration loop:
        delta_max = 0.0
        delta_max = max(abs(delta[:]))
        # print ('k_iter=',k_iter, 'delta_max=',delta_max)
        # Test convergence and recalculate forcing function and
        # central coefficient b if necessary:
        if delta_max < delta_acc:
            flag_conv = False
        else:
            if config.predic_correc:
                n = np.zeros(n_max)
                p = np.zeros(n_max)
                n, p, fi_non, EF = equi_np_fi3(
                    fi_out,
                    wfh_general,
                    wfe_general,
                    model,
                    E_state_general,
                    E_statec_general,
                    meff_state_general,
                    meff_statec_general,
                    n_max,
                    ni,
                )
            else:
                n, p, fi_non = equi_np_fi22(
                    vindex,
                    fitotc,
                    fitot,
                    Nc,
                    Nv,
                    fi_e,
                    fi_h,
                    iteration,
                    fi_out,
                    Vt,
                    wfh_general,
                    wfe_general,
                    model,
                    E_state_general,
                    E_statec_general,
                    meff_state_general,
                    meff_statec_general,
                    dop,
                    Ppz_Psp,
                    n_max,
                    ni,
                    n,
                    p,
                )
            for i in range(1, n_max - 1):
                b[i] = -(2 * Ldi2[i] / (dx2) + n[i] + p[i])
                f[i] = (
                    n[i]
                    - p[i]
                    - dop_out[i]
                    - Ppz_Psp_out[i]
                    - fi_out[i] * (n[i] + p[i])
                )
    return n, p, fi_out, fi_stat


def fd1(Ei, Ef, model):  # use
    """integral of Fermi Dirac Equation for energy independent density of states.
    Ei [meV], Ef [meV], T [K]"""
    T = model.T
    return kb * T * log(exp(meV2J * (Ei - Ef) / (kb * T)) + 1)


def fd2(Ei, Ef, model):
    """integral of Fermi Dirac Equation for energy independent density of states.
    Ei [meV], Ef [meV], T [K]"""
    T = model.T
    return kb * T * log(exp(meV2J * (Ef - Ei) / (kb * T)) + 1)


def fd3(x):
    """
    Approximation used for the Fermi integral by:
    D. Bednarczyk and J. Bednarczyk, The approximation of the Fermi-Dirac
    integral f1/2(η), Phys. Lett., vol. 64A, pp. 409–410, 1978.
    """
    return 1 / (
        exp(-x)
        + 3
        / 4
        * sqrt(pi)
        * (x ** 4 + 50 + 33.6 * x * (1 - 0.68 * exp(-0.17 * (x + 1) ** 2))) ** (-3 / 8)
    )


def fd4(x):
    """
    Approximation used for the Fermi integral by:
    Ehrenberg W (1950) The electric conductivity of simple semiconductors. Proc Phys Soc A 63:75
    """
    return 2 * sqrt(pi) * (exp(x)) / (4 + exp(x))


def fd5(x):
    """
    Approximation used for the Fermi integral by:
    Shun_Lien_Chuang Physics of Photonic Devices, 2009 by John Wiley & Sons, p 36.
    """
    return 4 / 3 * sqrt(abs(x) ** 3 / pi)


def fd6(x):
    x = abs(x)
    """
    Approximation used for the Fermi integral by:

    """
    return (
        log(x) / (x ** 2 - 1)
        + (3 * sqrt(pi) * x / 4) ** (2 / 3)
        + (2 * 3 * sqrt(pi) * x / 4) / (3 + 3 * sqrt(pi) * x / 4) ** 2
    )


def fd7(x):
    """
    Approximation used for the Fermi integral by:
    D. Bednarczyk and J. Bednarczyk, The approximation of the Fermi-Dirac
    integral f1/2(η), Phys. Lett., vol. 64A, pp. 409–410, 1978.
    """
    return exp(-x) + 3 / 4 * sqrt(pi) * (
        x ** 4 + 50 + 33.6 * x * (1 - 0.68 * exp(-0.17 * (x + 1) ** 2))
    ) ** (-3 / 8)


def amort_wave(j, Well_boundary, n_max):
    import config

    I11 = Well_boundary[j, 0]
    I22 = Well_boundary[j, 1]
    amort_wave_0 = int(
        config.amort_wave_0 * (Well_boundary[j, 1] - Well_boundary[j, 0]) / 2
    )
    amort_wave_1 = int(
        config.amort_wave_1 * (Well_boundary[j, 1] - Well_boundary[j, 0]) / 2
    )

    I1 = I11 - amort_wave_0  # n_max-70#
    I2 = I22 + amort_wave_1  # n_max-5#
    return I1, I2, I11, I22


def equi_np_fi2(
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    iteration,
    fi_old,
    Vt,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    dop,
    Ppz_Psp,
    n_max,
    ni,
):  # use
    n = np.zeros(n_max)
    n_cl = np.zeros(n_max)
    n_qw = np.zeros(n_max)
    p = np.zeros(n_max)
    p_cl = np.zeros(n_max)
    p_qw = np.zeros(n_max)
    Ef_Ec = np.zeros(n_max)
    Ev_Ef = np.zeros(n_max)
    E3kbT = np.zeros(n_max)
    # Determination of the Fermi Level
    ###################################################
    """
    nn=20000
    n1=np.zeros(2*nn)
    n2=np.zeros(2*nn)
    p2=np.zeros(2*nn)    
    n3=np.zeros(2*nn)
    n7=np.zeros(2*nn)
    p1=np.zeros(2*nn)        
    for i5 in range(-nn,nn):
        p1[i5+nn]=-(i5)/(0.0001*nn)#*meV2J/(kb*T)#Nv[i1]*fd3((fi_h[i1]-Vt*q*fi_old[i1])/(kb*T))/ni[i1]
        p2[i5+nn]=-(i5)/(0.0001*nn)/1000
        n1[i5+nn]=fd3(p1[i5+nn]*meV2J/(kb*T))#Nc[i1]*fd3(-(fi_e[i1]-Vt*q*fi_old[i1])/(kb*T))/ni[i1]
        n2[i5+nn]=fd4(p1[i5+nn]*meV2J/(kb*T))
        n3[i5+nn]=exp(p1[i5+nn]*meV2J/(kb*T))
        
        #n7[i5+nn]=fd7(p1[i5+nn]*meV2J/(kb*T))
        
        #print(n1[i5+nn])

    xaxis = np.arange(0,n_max)*model.dx
    #pl.plot(p1, n7)
    pl.plot(p2, n1,'r',p2,n2,'b',p2,n3,'k')
    #pl.plot(xaxis, Ef_Ec,'r',xaxis,E3kbT,'b')
    pl.xlabel('Ef-Ec (meV)')
    pl.ylabel('Fermi integral F1/2' )
    pl.title('')
    pl.grid(True)
    pl.legend(('fd3','fd4','exp'),loc='best',fontsize=12)    
    pl.show()     
    sssss
    """
    #######################################################
    EF = 0.0
    for i1 in range(0, n_max):
        Ef_Ec[i1] = (EF - (fi_e[i1] - Vt * q * fi_old[i1])) / (kb * T)  # *J2meV
        Ev_Ef[i1] = ((fi_h[i1] - Vt * q * fi_old[i1]) - EF) / (kb * T)  # *J2meV
        E3kbT[i1] = -3 * kb * T / (kb * T)  # *J2meV
        # print('Ef_Ec[%d]='%i1,Ef_Ec[i1]*J2meV*(kb*T),'meV')
        # print('Eg[%d]='%i1,abs(Ef_Ec[i1]*J2meV*(kb*T)+Ev_Ef[i1]*J2meV*(kb*T)),'meV')
        if Ef_Ec[i1] * J2meV * (kb * T) > 1000000000000:
            print("Ef_Ec[%d]=" % i1, Ef_Ec[i1] * J2meV * (kb * T), "meV")
            print("fi_old=", fi_old[i1] * J2meV, "meV")
            exit()
        if EF * meV2J - (fi_e[i1] - Vt * q * fi_old[i1]) > -3 * kb * T:  # -3*kb*T
            n_cl[i1] = (
                Nc[i1]
                * fd3((EF * meV2J - (fi_e[i1] - Vt * q * fi_old[i1])) / (kb * T))
                / ni[i1]
            )
        else:
            n_cl[i1] = (
                Nc[i1]
                * exp((EF * meV2J - (fi_e[i1] - Vt * q * fi_old[i1])) / (kb * T))
                / ni[i1]
            )
            # print(n[i1])
        if (fi_h[i1] - Vt * q * fi_old[i1]) - EF * meV2J > -3 * kb * T:  # -3*kb*T
            p_cl[i1] = (
                Nv[i1]
                * fd3(((fi_h[i1] - Vt * q * fi_old[i1]) - EF * meV2J) / (kb * T))
                / ni[i1]
            )
        else:
            p_cl[i1] = (
                Nv[i1]
                * exp(((fi_h[i1] - Vt * q * fi_old[i1]) - EF * meV2J) / (kb * T))
                / ni[i1]
            )
    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        i1 = I1 - I1
        # See C. de Faco et al. / Journal of Computational Physics 204 (2005) page 538 for
        # the jump discontinuity in the electron density n across the interface
        n_cl[I1:I2] = 0.0
        p_cl[I1:I2] = 0.0
        # couter=0
        for j in range(0, model.subnumber_e, 1):
            for i in range(I1, I2):
                n_qw[i] += (
                    (
                        fd2(E_statec_general[k, j], EF, model)
                        * meff_statec_general[k, j]
                        / (hbar ** 2 * pi)
                    )
                    * (wfe_general[k, j, i - I1]) ** 2
                    / (ni[i] * model.dx)
                )
        for jj in range(0, model.subnumber_h, 1):
            for ii in range(I1, I2):
                p_qw[ii] += (
                    (
                        fd1(E_state_general[k, jj], EF, model)
                        * meff_state_general[k, jj]
                        / (hbar ** 2 * pi)
                    )
                    * (wfh_general[k, jj, ii - I1]) ** 2
                    / (ni[ii] * model.dx)
                )
    n = n_cl + n_qw
    p = p_cl + p_qw
    return n, p, fi_old, EF  # density of carriers


def equi_np_fi4(
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    iteration,
    fi_old,
    Vt,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    dop,
    Ppz_Psp,
    n_max,
    ni,
):  # use
    n = np.zeros(n_max)
    n_cl = np.zeros(n_max)
    n_qw = np.zeros(n_max)
    p = np.zeros(n_max)
    p_cl = np.zeros(n_max)
    p_qw = np.zeros(n_max)
    Ef_Ec = np.zeros(n_max)
    Ev_Ef = np.zeros(n_max)
    E3kbT = np.zeros(n_max)
    # Determination of the Fermi Level
    EF = 0.0
    for i1 in range(0, n_max):
        Ef_Ec[i1] = (EF - (fi_e[i1] - Vt * q * fi_old[i1])) / (kb * T)  # *J2meV
        Ev_Ef[i1] = ((fi_h[i1] - Vt * q * fi_old[i1]) - EF) / (kb * T)  # *J2meV
        E3kbT[i1] = -3 * kb * T / (kb * T)  # *J2meV
        if EF * meV2J - (fi_e[i1] - Vt * q * fi_old[i1]) > -3 * kb * T:  # -3*kb*T
            n_cl[i1] = (
                Nc[i1]
                * fd3((EF * meV2J - (fi_e[i1] - Vt * q * fi_old[i1])) / (kb * T))
                / ni[i1]
            )
        else:
            n_cl[i1] = (
                Nc[i1]
                * exp((EF * meV2J - (fi_e[i1] - Vt * q * fi_old[i1])) / (kb * T))
                / ni[i1]
            )
            # print(n[i1])
        if (fi_h[i1] - Vt * q * fi_old[i1]) - EF * meV2J > -3 * kb * T:  # -3*kb*T
            p_cl[i1] = (
                Nv[i1]
                * fd3(((fi_h[i1] - Vt * q * fi_old[i1]) - EF * meV2J) / (kb * T))
                / ni[i1]
            )
        else:
            p_cl[i1] = (
                Nv[i1]
                * exp(((fi_h[i1] - Vt * q * fi_old[i1]) - EF * meV2J) / (kb * T))
                / ni[i1]
            )
    if config.quantum_effect:

        for k in range(1, model.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
            i1 = I1 - I1
            """
            See C. de Faco et al. / Journal of Computational Physics 204 (2005) page 538 for 
            the jump discontinuity in the electron density n across the interface
            """
            n_cl[I1:I2] = 0.0
            p_cl[I1:I2] = 0.0

            # couter=0
            for j in range(0, model.subnumber_e, 1):
                for i in range(I1, I2):
                    n_qw[i] += (
                        (
                            fd2(E_statec_general[k, j], EF, model)
                            * meff_statec_general[k, j]
                            / (hbar ** 2 * pi)
                        )
                        * (wfe_general[k, j, i - I1]) ** 2
                        / (ni[i] * model.dx)
                    )
                    """
                    if (fi_e[i]-Vt*q*fi_old[i])*J2meV<E_statec_general[k,j] and couter==0:                    
                        n_cl[i]+=Nc[i]*fd4((EF*meV2J-(fi_e[i]-Vt*q*fi_old[i]))/(kb*T))*(max(wfe_general[k,j,0:I2])**2)/ni[i]
                        couter+=1
                    elif couter!=0:
                        if (fi_e[i]-Vt*q*fi_old[i])*J2meV<E_statec_general[k,j] and (fi_e[i]-Vt*q*fi_old[i])*J2meV>E_statec_general[k,j-1]:
                            n_cl[i]+=Nc[i]*fd4((EF*meV2J-(fi_e[i]-Vt*q*fi_old[i]))/(kb*T))*(max(wfe_general[k,j,0:I2])**2)/ni[i]
                    """
            for jj in range(0, model.subnumber_h, 1):
                for ii in range(I1, I2):
                    p_qw[ii] += (
                        (
                            fd1(E_state_general[k, jj], EF, model)
                            * meff_state_general[k, jj]
                            / (hbar ** 2 * pi)
                        )
                        * (wfh_general[k, jj, ii - I1]) ** 2
                        / (ni[ii] * model.dx)
                    )
    n = n_cl + n_qw
    p = p_cl + p_qw
    return n, p, fi_old, EF  # density of carriers


def equi_np_fi3(
    fi_old,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    n_max,
    ni,
):  # use
    n = np.exp(fi_old)
    p = np.exp(-fi_old)
    fi_stat = np.zeros(n_max)
    Delta_fi = np.zeros(n_max)
    if not (config.predic_correc):
        fi_stat = fi_old
    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        n[I1:I2] = 0.0
        p[I1:I2] = 0.0
        for i in range(I1, I2):
            Delta_fi[i] = -Vt * q * fi_stat[i] - (-Vt * q * fi_old[i])
            for j in range(0, model.subnumber_e, 1):
                #print("E_statec_general=",(E_statec_general[k, j] - Delta_fi[i] * J2meV))
                n[i] += (
                    (
                        fd2(
                            E_statec_general[k, j] - Delta_fi[i] * J2meV,
                            0.0,
                            model,
                        )
                        * meff_statec_general[k, j]
                        / (hbar ** 2 * pi)
                    )
                    * (wfe_general[k, j, i - I1]) ** 2
                    / (ni[i] * model.dx)
                )
            for jj in range(0, model.subnumber_h, 1):
                p[i] += (
                    (
                        fd1(
                            E_state_general[k, jj] - Delta_fi[i] * J2meV,
                            0.0,
                            model,
                        )
                        * meff_state_general[k, jj]
                        / (hbar ** 2 * pi)
                    )
                    * (wfh_general[k, jj, i - I1]) ** 2
                    / (ni[i] * model.dx)
                )
    return n, p, fi_old, 0.0


def Ber(x):
    flag_sum = True
    if x > 0.01:
        Ber_1 = x * exp(-x) / (1.0 - exp(-x))
    elif x < 0 and abs(x) > 0.01:
        Ber_1 = x / (exp(x) - 1.0)
    elif x == 0:
        Ber_1 = 1.0
    else:
        temp_term = 1.0
        sum_1 = temp_term
        j = 0.0
        while flag_sum:
            j = j + 1
            temp_term = temp_term * x / (j + 1)
            # print sum_1+temp_term,'=',sum_1
            if sum_1 + temp_term == sum_1:
                flag_sum = False
            sum_1 = sum_1 + temp_term
        Ber_1 = 1.0 / sum_1
    return Ber_1


def Ber2(x):
    if np.abs(x) < 1e-12:
        return 1.0
    else:
        return x / np.expm1(x)


def Ber3(x):
    if x < -0.7:
        return -x
    elif x < 0.009 and x > -0.009:
        return 1 - x / 2
    else:
        return x / (exp(x) + 1)


def Poisson_equi1(dx, dop, fi, n_max):
    fi_out = np.zeros(n_max)
    dop_out = np.zeros(n_max)
    d = np.zeros(n_max)
    v = np.zeros(n_max)
    f = np.zeros(n_max)
    c = np.zeros(n_max)
    b = np.zeros(n_max)
    a = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = fi
    dop_out = dop
    delta_acc = 1.0e-5
    dx2 = dx * dx
    for i in range(0, n_max):
        a[i] = 1.0 / (dx2)
        c[i] = 1.0 / (dx2)
        b[i] = -(2.0 / (dx2) + exp(fi_out[i]) + exp(-fi_out[i]))
        f[i] = (
            exp(fi_out[i])
            - exp(-fi_out[i])
            - dop_out[i]
            - fi_out[i] * (exp(fi_out[i]) + exp(-fi_out[i]))
        )
    # (B) Define the elements of the coefficient matrix and initialize the forcing
    # function at the ohmic contacts
    a[0] = 0.0
    c[0] = 0.0
    b[0] = 1.0
    f[0] = fi_out[0]
    a[n_max - 1] = 0.0
    c[n_max - 1] = 0.0
    b[n_max - 1] = 1.0
    f[n_max - 1] = fi_out[n_max - 1]
    # (C) Start the iterative procedure for the solution of the linearized Poisson
    # equation using LU decomposition method:
    flag_conv = True  # convergence of the Poisson loop
    k_iter = 0
    while flag_conv:
        k_iter = k_iter + 1
        d[0] = b[0]
        for i in range(1, n_max):
            d[i] = b[i] - a[i] * c[i - 1] / d[i - 1]
        # Solution of Lv=f:
        v[0] = f[0]
        for i in range(1, n_max):
            v[i] = f[i] - a[i] * v[i - 1] / d[i - 1]
        # Solution of U*fi=v:
        temp = v[n_max - 1] / d[n_max - 1]
        delta[n_max - 1] = temp - fi_out[n_max - 1]
        fi_out[n_max - 1] = temp
        for i in range(n_max - 2, -1, -1):
            temp = (v[i] - c[i] * fi_out[i + 1]) / d[i]
            delta[i] = temp - fi_out[i]
            fi_out[i] = temp
        # Test update in the outer iteration loop:
        delta_max = 0.0
        delta_max = max(abs(delta[:]))
        print("k_iter=", k_iter, "delta_max=", delta_max)
        # Test convergence and recalculate forcing function and
        # central coefficient b if necessary:
        if delta_max < delta_acc:
            flag_conv = False
        else:
            for i in range(1, n_max - 1):
                b[i] = -(2.0 / (dx2) + exp(fi_out[i]) + exp(-fi_out[i]))
                f[i] = (
                    exp(fi_out[i])
                    - exp(-fi_out[i])
                    - dop_out[i]
                    - fi_out[i] * (exp(fi_out[i]) + exp(-fi_out[i]))
                )
    return fi_out, a, b, c, d, f, v


def Mobility1(mun0, mup0, fi, Vt, Ldi, VSATN, VSATP, BETAN, BETAP, n_max, dx):
    #######################################################################
    #% 3.1 . Calculate Field Dependant Mobility for each value of 'fi'   ##
    #%       at each node point of the PN diode.                         ##
    #######################################################################
    Ef = np.zeros(n_max)
    mup = np.zeros(n_max)
    mun = np.zeros(n_max)

    """
    ### To test with Constant Mobility without field dependancy.        
    for i in range(0,n_max):           # Start Loop for Field Dep Mobility 
        mup[i] = mup0
        mun[i] = mun0
    # #           
    """
    # [0] Solution of electron current continuity equation:
    # .................
    # (1a) Define the elements of the coefficient matrix and
    # initialize the forcing function:
    ## Calculate the Electric Field at each Node
    for i in range(0, n_max - 1):
        Ef[i] = abs(fi[i] - fi[i + 1]) * Vt / (dx * Ldi)
    Ef[0] = Ef[1]
    Ef[n_max - 1] = Ef[n_max - 2]
    ## Calculate the Field Dependant Mobility at each Node
    for i in range(0, n_max):
        pdeno = (mup0 * Ef[i] / VSATP) ** BETAP
        mup[i] = mup0 * ((1 / (1 + pdeno)) ** (1 / BETAP))
        ndeno = (mun0 * Ef[i] / VSATN) ** BETAN
        mun[i] = mun0 * ((1 / (1 + ndeno)) ** (1 / BETAN))
    mup[0] = mup[1]
    mup[n_max - 1] = mup[n_max - 2]
    mun[0] = mun[1]
    mun[n_max - 1] = mun[n_max - 2]
    return mun, mup


def Continuity1(n, p, mun, mup, fi, Vt, Ldi, n_max, dx, TAUN0, TAUP0):
    #################################################################################
    ## 3.2 Solve Continuity Equation for Electron and Holes using LU Decomposition ##
    #################################################################################
    # (A) Define the elements of the coefficient matrix and initialize the forcing
    #    function at the ohmic contacts for ELECTRON and HOLE Continuity
    vp = np.zeros(n_max)
    dp = np.zeros(n_max)
    vn = np.zeros(n_max)
    fp = np.zeros(n_max)
    cp = np.zeros(n_max)
    bp = np.zeros(n_max)
    ap = np.zeros(n_max)
    fn = np.zeros(n_max)
    cn = np.zeros(n_max)
    bn = np.zeros(n_max)
    dn = np.zeros(n_max)
    an = np.zeros(n_max)
    betan = np.zeros(n_max)
    betap = np.zeros(n_max)
    dx2 = dx * dx
    an[0] = 0  # Co-ef for electron at Anode
    bn[0] = 1  # Co-ef for electron at Anode
    cn[0] = 0  # Co-ef for electron at Anode
    ap[0] = 0  # Co-ef for hole     at Anode
    bp[0] = 1  # Co-ef for hole     at Anode
    cp[0] = 0  # Co-ef for hole     at Anode
    # fnp[0] = (Ldi*Ldi*dx2/Vt) * ( p[0]*n[0] - 1 ) / ( TAUP0*(n[0] + 1 ) + TAUN0*(p[0] + 1 ) )
    fn[0] = n[0]
    fp[0] = p[0]
    an[n_max - 1] = 0  # Co-ef for electron at Cathode
    bn[n_max - 1] = 1  # Co-ef for electron at Cathode
    cn[n_max - 1] = 0  # Co-ef for electron at Cathode
    ap[n_max - 1] = 0  # Co-ef for hole     at Cathode
    bp[n_max - 1] = 1  # Co-ef for hole     at Cathode
    cp[n_max - 1] = 0  # Co-ef for hole     at Cathode
    # fnp[n_max-1] = (Ldi*Ldi*dx2/Vt) * ( p[n_max-1]*n[n_max-1] - 1 ) / ( TAUP0*(n[n_max-1] + 1) + TAUN0*(p[n_max-1] + 1) )
    fn[n_max - 1] = n[n_max - 1]
    fp[n_max - 1] = p[n_max - 1]
    # (B) Define the elements of the coefficient matrix for the internal nodes and
    #    initialize the forcing function
    for i in range(1, n_max - 1):
        munim1by2 = (mun[i - 1] + mun[i]) / 2
        munip1by2 = (mun[i] + mun[i + 1]) / 2
        mupim1by2 = (mup[i - 1] + mup[i]) / 2
        mupip1by2 = (mup[i] + mup[i + 1]) / 2
        ## Co-efficients for HOLE Continuity eqn
        ap[i] = mupim1by2 * Ber((fi[i] - fi[i - 1]))
        cp[i] = mupip1by2 * Ber((fi[i] - fi[i + 1]))
        bp[i] = -(
            mupim1by2 * Ber((fi[i - 1] - fi[i])) + mupip1by2 * Ber((fi[i + 1] - fi[i]))
        )
        ## Co-efficients for ELECTRON Continuity eqn
        an[i] = munim1by2 * Ber((fi[i - 1] - fi[i]))
        cn[i] = munip1by2 * Ber((fi[i + 1] - fi[i]))
        bn[i] = -(
            munim1by2 * Ber((fi[i] - fi[i - 1])) + munip1by2 * Ber(fi[i] - fi[i + 1])
        )
        ## Forcing Function for ELECTRON and HOLE Continuity eqns
        fn[i] = (
            (Ldi * Ldi * dx2 / Vt)
            * (p[i] * n[i] - 1)
            / (TAUP0 * (n[i] + 1) + TAUN0 * (p[i] + 1))
        )
        fp[i] = (
            (Ldi * Ldi * dx2 / Vt)
            * (p[i] * n[i] - 1)
            / (TAUP0 * (n[i] + 1) + TAUN0 * (p[i] + 1))
        )
    # (C)  Start the iterative procedure for the solution of the linearized Continuity
    #     equation for "ELECTRONS" using LU decomposition method:
    dn[0] = bn[0]
    for i in range(1, n_max):
        betan[i] = an[i] / dn[i - 1]
        dn[i] = bn[i] - betan[i] * cn[i - 1]
    # Solution of Lv = f #
    vn[0] = fn[0]
    for i in range(1, n_max):
        vn[i] = fn[i] - betan[i] * vn[i - 1]
    # Solution of U*fi = v #
    tempn = vn[n_max - 1] / dn[n_max - 1]
    # deltan[n_max-1] = tempn - n[n_max-1]
    n[n_max - 1] = tempn
    for i in range(n_max - 2, -1, -1):  # delta#
        tempn = (vn[i] - cn[i] * n[i + 1]) / dn[i]
        #  deltan[i] = tempn - n[i]
        n[i] = tempn
    ####################### END of ELECTRON Continuty Solver ###########
    # (D)  Start the iterative procedure for the solution of the linearized Continuity
    #     equation for "HOLES" using LU decomposition method:
    dp[0] = bp[0]
    for i in range(1, n_max):
        betap[i] = ap[i] / dp[i - 1]
        dp[i] = bp[i] - betap[i] * cp[i - 1]
    # Solution of Lv = f #
    vp[0] = fp[0]
    for i in range(1, n_max):
        vp[i] = fp[i] - betap[i] * vp[i - 1]
    # Solution of U*fi = v #
    tempp = vp[n_max - 1] / dp[n_max - 1]
    # deltap[n_max-1] = tempp - p[n_max-1]
    p[n_max - 1] = tempp
    for i in range(n_max - 2, -1, -1):  # delta#
        tempp = (vp[i] - cp[i] * p[i + 1]) / dp[i]
        #   deltap[i] = tempp - p[i]
        p[i] = tempp
    ####################### END of HOLE Continuty Solver ###########
    return n, p


def Poisson_non_equi1(n, p, dop, n_max, dx, fi, flag_conv_2):
    ####################################################################
    ## 3.3 Calculate potential fi again with new values of "n" and "p"##
    ##     and check for convergence                                  ##
    ####################################################################
    # Recalculate forcing function and central coefficient b for fi
    d = np.zeros(n_max)
    v = np.zeros(n_max)
    f = np.zeros(n_max)
    c = np.zeros(n_max)
    b = np.zeros(n_max)
    a = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = np.zeros(n_max)
    dop_out = np.zeros(n_max)
    d = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = fi
    dop_out = dop
    dx2 = dx * dx
    delta_acc = 1.0e-5
    for i in range(1, n_max - 1):
        a[i] = 1.0 / (dx2)
        c[i] = 1.0 / (dx2)
        b[i] = -(2 / (dx2) + n[i] + p[i])
        f[i] = n[i] - p[i] - dop_out[i] - (fi_out[i] * (n[i] + p[i]))
    a[0] = 0.0
    c[0] = 0.0
    b[0] = 1.0
    f[0] = fi_out[0]
    a[n_max - 1] = 0.0
    c[n_max - 1] = 0.0
    b[n_max - 1] = 1.0
    f[n_max - 1] = fi_out[n_max - 1]
    ## here values of n[i] and p[i] are used in place of exp(fi[i])
    # Solve for Updated potential given the new value of Forcing
    # Function using LU decomposition
    d[0] = b[0]
    for i in range(1, n_max):
        d[i] = b[i] - a[i] * c[i - 1] / d[i - 1]
    # Solution of Lv = f #
    v[0] = f[0]
    for i in range(1, n_max):
        v[i] = f[i] - a[i] * v[i - 1] / d[i - 1]
    # Solution of U*fi = v #
    temp = v[n_max - 1] / d[n_max - 1]
    delta[n_max - 1] = temp - fi_out[n_max - 1]
    fi_out[n_max - 1] = temp
    for i in range(n_max - 2, -1, -1):  # delta#
        temp = (v[i] - c[i] * fi_out[i + 1]) / d[i]
        delta[i] = temp - fi_out[i]
        fi_out[i] = temp
    delta_max = 0
    delta_max = max(abs(delta[:]))
    # Test convergence and start the loop if necessary else increase
    # the applied potential
    print("delta_max= ", delta_max)
    if delta_max < delta_acc:
        flag_conv_2 = False
    else:
        for i in range(1, n_max - 1):
            b[i] = -(2 / (dx2) + n[i] + p[i])
            f[i] = n[i] - p[i] - dop_out[i] - (fi_out[i] * (n[i] + p[i]))
    return fi_out, flag_conv_2


def Current1(
    vindex,
    n,
    p,
    mun,
    mup,
    fi,
    Vt,
    n_max,
    Total_Steps,
    q,
    dx,
    ni,
    Ldi,
    Jnip1by2,
    Jnim1by2,
    Jelec,
    Jpip1by2,
    Jpim1by2,
    Jhole,
):
    ##########################################################################
    ##                        CALCULATE CURRENT                             ##
    ##########################################################################
    for i in range(1, n_max - 1):

        # Electron Current
        Jnip1by2[vindex, i] = (
            (q * mun[i] * Vt / (dx * Ldi))
            * ni
            * (n[i + 1] * Ber((fi[i + 1] - fi[i])) - n[i] * Ber((fi[i] - fi[i + 1])))
        )
        Jnim1by2[vindex, i] = (
            (q * mun[i] * Vt / (dx * Ldi))
            * ni
            * (n[i] * Ber((fi[i] - fi[i - 1])) - n[i - 1] * Ber((fi[i - 1] - fi[i])))
        )
        Jelec[vindex, i] = (Jnip1by2[vindex, i] + Jnim1by2[vindex, i]) / 2
        # Hole Current
        Jpip1by2[vindex, i] = (
            (q * mup[i] * Vt / (dx * Ldi))
            * ni
            * (p[i + 1] * Ber((fi[i] - fi[i + 1])) - p[i] * Ber((fi[i + 1] - fi[i])))
        )
        Jpim1by2[vindex, i] = (
            (q * mup[i] * Vt / (dx * Ldi))
            * ni
            * (p[i] * Ber((fi[i - 1] - fi[i])) - p[i - 1] * Ber((fi[i] - fi[i - 1])))
        )
        Jhole[vindex, i] = (Jpip1by2[vindex, i] + Jpim1by2[vindex, i]) / 2
    ##         Jtotal(vindex) = Jelec
    return Jnip1by2, Jnim1by2, Jelec, Jpip1by2, Jpim1by2, Jhole


def Write_results_equi1(dEc, Vt, q, ni, n, p, dop, dx, Ldi, fi, n_max):
    ro = np.zeros(n_max)
    el_field1 = np.zeros(n_max)
    el_field2 = np.zeros(n_max)
    Ec = np.zeros(n_max)
    for i in range(1, n_max - 1):
        Ec[i] = dEc - Vt * fi[i]
        ro[i] = -q * ni * (exp(fi[i]) - exp(-fi[i]) - dop[i])
        el_field1[i] = -(fi[i + 1] - fi[i]) * Vt / (dx * Ldi)
        el_field2[i] = -(fi[i + 1] - fi[i - 1]) * Vt / (2.0 * dx * Ldi)
        n[i] = exp(fi[i])
        p[i] = exp(-fi[i])
    Ec[0] = Ec[1]
    Ec[n_max - 1] = Ec[n_max - 2]
    ro[0] = ro[1]
    ro[n_max - 1] = ro[n_max - 2]
    el_field1[0] = el_field1[1]
    el_field1[n_max - 1] = el_field1[n_max - 2]
    el_field2[0] = el_field2[1]
    el_field2[n_max - 1] = el_field2[n_max - 2]
    n[0] = n[1]
    n[n_max - 1] = n[n_max - 2]
    p[0] = p[1]
    p[n_max - 1] = p[n_max - 2]
    """
    fi[0]=fi[1]
    fi[n_max-1] = fi[n_max-2] 
    """
    Ec_result = np.zeros(n_max)
    ro_result = np.zeros(n_max)
    el_field1_result = np.zeros(n_max)
    el_field2_result = np.zeros(n_max)
    nf_result = np.zeros(n_max)
    pf_result = np.zeros(n_max)
    fi_result = np.zeros(n_max)
    Ec_result = Ec[0:n_max]
    ro_result = ro[0:n_max]
    el_field1_result = el_field1[0:n_max]
    el_field2_result = el_field2[0:n_max]
    nf_result = n[0:n_max] * ni
    pf_result = p[0:n_max] * ni
    fi_result = fi[0:n_max]
    return (
        Ec_result,
        ro_result,
        el_field1_result,
        el_field2_result,
        nf_result,
        pf_result,
        fi_result,
    )


def Write_results_equi2(ns, fitotc, fitot, Vt, q, ni, n, p, dop, dx, Ldi, fi, n_max):
    ro = np.zeros(n_max)
    el_field1 = np.zeros(n_max)
    el_field2 = np.zeros(n_max)
    Ec = np.zeros(n_max)
    Ev = np.zeros(n_max)
    for i in range(1, n_max - 1):
        Ec[i] = fitotc[i]
        Ev[i] = fitot[i]
        ro[i] = -q * (ni[i] * n[i] - ni[i] * p[i] - dop[i])
        el_field1[i] = -(fi[i + 1] - fi[i]) * Vt / (dx)
        el_field2[i] = -(fi[i + 1] - fi[i - 1]) * Vt / (2.0 * dx)
        n[i] = n[i]  # exp(fi[i])
        p[i] = p[i]  # exp(-fi[i])
    Ec[0] = Ec[1]
    Ec[n_max - 1] = Ec[n_max - 2]
    Ev[0] = Ev[1]
    Ev[n_max - 1] = Ev[n_max - 2]
    ro[0] = ro[1]
    ro[n_max - 1] = ro[n_max - 2]
    el_field1[0] = el_field1[1]
    el_field1[n_max - 1] = el_field1[n_max - 2]
    el_field2[0] = el_field2[1]
    el_field2[n_max - 1] = el_field2[n_max - 2]
    n[0] = n[1]
    n[n_max - 1] = n[n_max - 2]
    p[0] = p[1]
    p[n_max - 1] = p[n_max - 2]

    fi[0] = fi[1]
    fi[n_max - 1] = fi[n_max - 2]
    Ec_result = np.zeros(n_max)
    Ev_result = np.zeros(n_max)
    ro_result = np.zeros(n_max)
    el_field1_result = np.zeros(n_max)
    el_field2_result = np.zeros(n_max)
    nf_result = np.zeros(n_max)
    pf_result = np.zeros(n_max)
    fi_result = np.zeros(n_max)
    Ec_result = Ec[0:n_max]
    Ev_result = Ev[0:n_max]
    ro_result = ro[0:n_max]
    el_field1_result = el_field1[0:n_max]
    el_field2_result = el_field2[0:n_max]
    nf_result = n[0:n_max] * ni[0:n_max]
    pf_result = p[0:n_max] * ni[0:n_max]
    fi_result = fi[0:n_max]
    return (
        Ec_result,
        Ev_result,
        ro_result,
        el_field1_result,
        el_field2_result,
        nf_result,
        pf_result,
        fi_result,
    )


def Write_results_non_equi1(
    dEc,
    Vt,
    q,
    ni,
    n,
    p,
    dop,
    dx,
    Ldi,
    fi,
    n_max,
    Jnip1by2,
    Jnim1by2,
    Jelec,
    Jpip1by2,
    Jpim1by2,
    Jhole,
    Jtotal,
    Total_Steps,
):

    axis = np.zeros(n_max)
    ro = np.zeros(n_max)
    el_field1 = np.zeros(n_max)
    el_field2 = np.zeros(n_max)

    Ec = np.zeros(n_max)
    Ev = np.zeros(n_max)
    Ei = np.zeros(n_max)
    Efn = np.zeros(n_max)
    Efp = np.zeros(n_max)
    av_curr = np.zeros(Total_Steps)
    axis[0] = dx * 1e4
    for i in range(1, n_max - 1):
        Ec[i] = dEc - Vt * fi[i]  # Values from the second Node%
        ro[i] = -q * ni * (n[i] - p[i] - dop[i])
        el_field1[i] = -(fi[i + 1] - fi[i]) * Vt / (dx * Ldi)
        el_field2[i] = -(fi[i + 1] - fi[i - 1]) * Vt / (2 * dx * Ldi)
        axis[i] = axis[i - 1] + dx * Ldi * 1e4
    Jtotal[:, 0] = Jtotal[:, 1]
    Jelec[:, 0] = Jelec[:, 1]
    Jhole[:, 0] = Jhole[:, 1]
    Jtotal[:, n_max - 1] = Jtotal[:, n_max - 2]
    Jelec[:, n_max - 1] = Jelec[:, n_max - 2]
    Jhole[:, n_max - 1] = Jhole[:, n_max - 2]

    Ec[0] = Ec[1]
    Ec[n_max - 1] = Ec[n_max - 2]
    axis[n_max - 1] = axis[n_max - 2] + dx * Ldi * 1e4
    el_field1[0] = el_field1[1]
    el_field2[0] = el_field2[1]
    el_field1[n_max - 1] = el_field1[n_max - 2]
    el_field2[n_max - 1] = el_field2[n_max - 2]
    ro[0] = ro[1]
    ro[n_max - 1] = ro[n_max - 2]
    n[0] = n[1]
    n[n_max - 1] = n[n_max - 2]
    p[0] = p[1]
    p[n_max - 1] = p[n_max - 2]
    nf = n * ni
    pf = p * ni
    Ev = Ec - 1.12
    ## Calculate Quasi Fermi Level - Efn Efp
    for i in range(0, n_max):
        Ei[i] = Ec[i] - 0.56
        Efn[i] = Ei[i] + Vt * log(nf[i] / ni)
        Efp[i] = Ei[i] - Vt * log(pf[i] / ni)
    Efn[0] = Efn[1]
    Efn[n_max - 1] = Efn[n_max - 2]
    Efp[0] = Efp[1]
    Efp[n_max - 1] = Efp[n_max - 2]
    for j in range(0, Total_Steps):
        av_curr[j] = Jtotal[j, 0]
    Ec_result = np.zeros(n_max)
    Ev_result = np.zeros(n_max)
    Ei_result = np.zeros(n_max)
    Efp_result = np.zeros(n_max)
    Efn_result = np.zeros(n_max)
    fi_result = np.zeros(n_max)
    fi_result = fi[0:n_max]
    Efn_result = Efn[0:n_max]
    Efp_result = Efp[0:n_max]
    ro_result = np.zeros(n_max)
    el_field1_result = np.zeros(n_max)
    el_field2_result = np.zeros(n_max)
    nf_result = np.zeros(n_max)
    pf_result = np.zeros(n_max)
    Ec_result = Ec[0:n_max]
    Ev_result = Ev[0:n_max]
    Ei_result = Ei[0:n_max]
    ro_result = ro[0:n_max]
    el_field1_result = el_field1[0:n_max]
    el_field2_result = el_field2[0:n_max]
    nf_result = n[0:n_max] * ni
    pf_result = p[0:n_max] * ni
    return (
        fi_result,
        Efn_result,
        Efp_result,
        ro_result,
        el_field1_result,
        el_field2_result,
        nf_result,
        pf_result,
        Ec_result,
        Ev_result,
        Ei_result,
        axis,
        av_curr,
    )


def equi_np_fi(iteration, dop, Ppz_Psp, n_max, ni, model, Vt, surface):
    n = np.zeros(n_max)
    p = np.zeros(n_max)
    fi_old = np.zeros(n_max)
    for i2 in range(0, n_max):
        nn, pp = np_cal(dop[i2] + Ppz_Psp[i2], ni[i2])
        n[i2] = nn / ni[i2]
        p[i2] = pp / ni[i2]
        fi_old[i2] = log(n[i2])
    fi_old[0] = fi_old[1]
    fi_old[n_max - 1] = fi_old[n_max - 2]

    fi_old[0] -= surface[0] / Vt
    fi_old[n_max - 1] -= surface[1] / Vt
    n = np.exp(fi_old)
    p = np.exp(-fi_old)
    return n, p, fi_old


def np_cal(dop, ni):
    ni2 = ni ** 2
    dop2 = dop ** 2
    if dop > 0:
        tmp = dop + sqrt(dop2 + 4 * ni2)
        n = tmp / 2.0
        p = ni2 / n
    else:
        dop = -1 * dop
        tmp = dop + sqrt(dop2 + 4 * ni2)
        p = tmp / 2.0
        n = ni2 / p
    return n, p


def Mobility2(mun0, mup0, fi, Vt, Ldi, VSATN, VSATP, BETAN, BETAP, n_max, dx):
    #######################################################################
    #% 3.1 . Calculate Field Dependant Mobility for each value of 'fi'   ##
    #%       at each node point of the PN diode.                         ##
    #######################################################################
    Efield = np.zeros(n_max)
    mup = np.zeros(n_max)
    mun = np.zeros(n_max)

    """
    ### To test with Constant Mobility without field dependancy.        
    for i in range(0,n_max):           # Start Loop for Field Dep Mobility 
        mup[i] = mup0[i]
        mun[i] = mun0[i]
    # #           
    """
    # [0] Solution of electron current continuity equation:
    # .................
    # (1a) Define the elements of the coefficient matrix and
    # initialize the forcing function:
    ## Calculate the Electric Field at each Node
    for i in range(0, n_max - 1):
        Efield[i] = abs(fi[i] - fi[i + 1]) * Vt / (dx)
    Efield[0] = Efield[1]
    Efield[n_max - 1] = Efield[n_max - 2]
    ## Calculate the Field Dependant Mobility at each Node
    for i in range(0, n_max):
        pdeno = (mup0[i] * Efield[i] / VSATP[i]) ** BETAP[i]
        mup[i] = mup0[i] * ((1 / (1 + pdeno)) ** (1 / BETAP[i]))
        ndeno = (mun0[i] * Efield[i] / VSATN[i]) ** BETAN[i]
        mun[i] = mun0[i] * ((1 / (1 + ndeno)) ** (1 / BETAN[i]))
    mup[0] = mup[1]
    mup[n_max - 1] = mup[n_max - 2]
    mun[0] = mun[1]
    mun[n_max - 1] = mun[n_max - 2]

    return mun, mup


def Continuity2(n, p, mun, mup, fi, Vt, Ldi, n_max, dx, TAUN0, TAUP0):
    #################################################################################
    ## 3.2 Solve Continuity Equation for Electron and Holes using LU Decomposition ##
    #################################################################################
    # (A) Define the elements of the coefficient matrix and initialize the forcing
    #    function at the ohmic contacts for ELECTRON and HOLE Continuity
    vp = np.zeros(n_max)
    dp = np.zeros(n_max)
    vn = np.zeros(n_max)
    fp = np.zeros(n_max)
    cp = np.zeros(n_max)
    bp = np.zeros(n_max)
    ap = np.zeros(n_max)
    fn = np.zeros(n_max)
    cn = np.zeros(n_max)
    bn = np.zeros(n_max)
    dn = np.zeros(n_max)
    an = np.zeros(n_max)
    betan = np.zeros(n_max)
    betap = np.zeros(n_max)
    fnp = np.zeros(n_max)
    dx2 = dx * dx
    an[0] = 0  # Co-ef for electron at Anode
    bn[0] = 1  # Co-ef for electron at Anode
    cn[0] = 0  # Co-ef for electron at Anode
    ap[0] = 0  # Co-ef for hole     at Anode
    bp[0] = 1  # Co-ef for hole     at Anode
    cp[0] = 0  # Co-ef for hole     at Anode
    fnp[0] = (
        (dx2 / Vt) * (p[0] * n[0] - 1) / (TAUP0[0] * (n[0] + 1) + TAUN0[0] * (p[0] + 1))
    )
    fn[0] = n[0]
    fp[0] = p[0]
    an[n_max - 1] = 0  # Co-ef for electron at Cathode
    bn[n_max - 1] = 1  # Co-ef for electron at Cathode
    cn[n_max - 1] = 0  # Co-ef for electron at Cathode
    ap[n_max - 1] = 0  # Co-ef for hole     at Cathode
    bp[n_max - 1] = 1  # Co-ef for hole     at Cathode
    cp[n_max - 1] = 0  # Co-ef for hole     at Cathode
    fnp[n_max - 1] = (
        (dx2 / Vt)
        * (p[n_max - 1] * n[n_max - 1] - 1)
        / (
            TAUP0[n_max - 1] * (n[n_max - 1] + 1)
            + TAUN0[n_max - 1] * (p[n_max - 1] + 1)
        )
    )
    fn[n_max - 1] = n[n_max - 1]
    fp[n_max - 1] = p[n_max - 1]
    # (B) Define the elements of the coefficient matrix for the internal nodes and
    #    initialize the forcing function
    for i in range(1, n_max - 1):
        munim1by2 = (mun[i - 1] + mun[i]) / 2
        munip1by2 = (mun[i] + mun[i + 1]) / 2
        mupim1by2 = (mup[i - 1] + mup[i]) / 2
        mupip1by2 = (mup[i] + mup[i + 1]) / 2
        ## Co-efficients for HOLE Continuity eqn
        ap[i] = mupim1by2 * Ber((fi[i] - fi[i - 1]))
        cp[i] = mupip1by2 * Ber((fi[i] - fi[i + 1]))
        bp[i] = -(
            mupim1by2 * Ber((fi[i - 1] - fi[i])) + mupip1by2 * Ber((fi[i + 1] - fi[i]))
        )
        ## Co-efficients for ELECTRON Continuity eqn
        an[i] = munim1by2 * Ber((fi[i - 1] - fi[i]))
        cn[i] = munip1by2 * Ber((fi[i + 1] - fi[i]))
        bn[i] = -(
            munim1by2 * Ber((fi[i] - fi[i - 1])) + munip1by2 * Ber(fi[i] - fi[i + 1])
        )
        ## Forcing Function for ELECTRON and HOLE Continuity eqns
        fn[i] = (
            (dx2 / Vt)
            * (p[i] * n[i] - 1)
            / (TAUP0[i] * (n[i] + 1) + TAUN0[i] * (p[i] + 1))
        )
        fp[i] = (
            (dx2 / Vt)
            * (p[i] * n[i] - 1)
            / (TAUP0[i] * (n[i] + 1) + TAUN0[i] * (p[i] + 1))
        )
    # (C)  Start the iterative procedure for the solution of the linearized Continuity
    #     equation for "ELECTRONS" using LU decomposition method:
    dn[0] = bn[0]
    for i in range(1, n_max):
        betan[i] = an[i] / dn[i - 1]
        dn[i] = bn[i] - betan[i] * cn[i - 1]
    # Solution of Lv = f #
    vn[0] = fn[0]
    for i in range(1, n_max):
        vn[i] = fn[i] - betan[i] * vn[i - 1]
    # Solution of U*fi = v #
    tempn = vn[n_max - 1] / dn[n_max - 1]
    # deltan[n_max-1] = tempn - n[n_max-1]
    n[n_max - 1] = tempn
    for i in range(n_max - 2, -1, -1):  # delta#
        tempn = (vn[i] - cn[i] * n[i + 1]) / dn[i]
        #  deltan[i] = tempn - n[i]
        n[i] = tempn
    ####################### END of ELECTRON Continuty Solver ###########
    # (D)  Start the iterative procedure for the solution of the linearized Continuity
    #     equation for "HOLES" using LU decomposition method:
    # print(max(n[:]))
    dp[0] = bp[0]
    for i in range(1, n_max):
        betap[i] = ap[i] / dp[i - 1]
        dp[i] = bp[i] - betap[i] * cp[i - 1]
    # Solution of Lv = f #
    vp[0] = fp[0]
    for i in range(1, n_max):
        vp[i] = fp[i] - betap[i] * vp[i - 1]
    # Solution of U*fi = v #
    tempp = vp[n_max - 1] / dp[n_max - 1]
    # deltap[n_max-1] = tempp - p[n_max-1]
    p[n_max - 1] = tempp
    for i in range(n_max - 2, -1, -1):  # delta#
        tempp = (vp[i] - cp[i] * p[i + 1]) / dp[i]
        #   deltap[i] = tempp - p[i]
        p[i] = tempp
    ####################### END of HOLE Continuty Solver ###########
    return n, p


def Mobility3(
    mun0, mup0, fi, fi_n, fi_p, Vt, Ldi, VSATN, VSATP, BETAN, BETAP, n_max, dx
):
    #######################################################################
    #% 3.1 . Calculate Field Dependant Mobility for each value of 'fi'   ##
    #%       at each node point of the PN diode.                         ##
    #######################################################################
    Efield = np.zeros(n_max)
    mup = np.zeros(n_max)
    mun = np.zeros(n_max)

    """
    ### To test with Constant Mobility without field dependancy.        
    for i in range(0,n_max):           # Start Loop for Field Dep Mobility 
        mup[i] = mup0
        mun[i] = mun0
    # #           
    """
    # [0] Solution of electron current continuity equation:
    # .................
    # (1a) Define the elements of the coefficient matrix and
    # initialize the forcing function:
    ## Calculate the Electric Field at each Node
    for i in range(0, n_max - 1):
        Efield[i] = abs(fi[i] - fi[i + 1]) * Vt / (dx)
    Efield[0] = Efield[1]
    Efield[n_max - 1] = Efield[n_max - 2]
    ## Calculate the Field Dependant Mobility at each Node
    for i in range(0, n_max):
        pdeno = (mup0[i] * Efield[i] / VSATP[i]) ** BETAP[i]
        mup[i] = mup0[i] * ((1 / (1 + pdeno)) ** (1 / BETAP[i]))
        ndeno = (mun0[i] * Efield[i] / VSATN[i]) ** BETAN[i]
        mun[i] = mun0[i] * ((1 / (1 + ndeno)) ** (1 / BETAN[i]))
    mup[0] = mup[1]
    mup[n_max - 1] = mup[n_max - 2]
    mun[0] = mun[1]
    mun[n_max - 1] = mun[n_max - 2]
    return mun, mup


def Continuity3(n, p, mun, mup, fi, fi_n, fi_p, Vt, Ldi, n_max, dx, TAUN0, TAUP0):
    #################################################################################
    ## 3.2 Solve Continuity Equation for Electron and Holes using LU Decomposition ##
    #################################################################################
    # (A) Define the elements of the coefficient matrix and initialize the forcing
    #    function at the ohmic contacts for ELECTRON and HOLE Continuity
    vp = np.zeros(n_max)
    dp = np.zeros(n_max)
    vn = np.zeros(n_max)
    fp = np.zeros(n_max)
    cp = np.zeros(n_max)
    bp = np.zeros(n_max)
    ap = np.zeros(n_max)
    fn = np.zeros(n_max)
    cn = np.zeros(n_max)
    bn = np.zeros(n_max)
    dn = np.zeros(n_max)
    an = np.zeros(n_max)
    betan = np.zeros(n_max)
    betap = np.zeros(n_max)
    dx2 = dx * dx
    an[0] = 0  # Co-ef for electron at Anode
    bn[0] = 1  # Co-ef for electron at Anode
    cn[0] = 0  # Co-ef for electron at Anode
    ap[0] = 0  # Co-ef for hole     at Anode
    bp[0] = 1  # Co-ef for hole     at Anode
    cp[0] = 0  # Co-ef for hole     at Anode
    # fnp[0] = (Ldi*Ldi*dx2/Vt) * ( p[0]*n[0] - 1 ) / ( TAUP0*(n[0] + 1 ) + TAUN0*(p[0] + 1 ) )
    fn[0] = n[0]
    fp[0] = p[0]
    an[n_max - 1] = 0  # Co-ef for electron at Cathode
    bn[n_max - 1] = 1  # Co-ef for electron at Cathode
    cn[n_max - 1] = 0  # Co-ef for electron at Cathode
    ap[n_max - 1] = 0  # Co-ef for hole     at Cathode
    bp[n_max - 1] = 1  # Co-ef for hole     at Cathode
    cp[n_max - 1] = 0  # Co-ef for hole     at Cathode
    # fnp[n_max-1] = (Ldi*Ldi*dx2/Vt) * ( p[n_max-1]*n[n_max-1] - 1 ) / ( TAUP0*(n[n_max-1] + 1) + TAUN0*(p[n_max-1] + 1) )
    fn[n_max - 1] = n[n_max - 1]
    fp[n_max - 1] = p[n_max - 1]
    # (B) Define the elements of the coefficient matrix for the internal nodes and
    #    initialize the forcing function
    for i in range(1, n_max - 1):
        munim1by2 = (mun[i - 1] + mun[i]) / 2
        munip1by2 = (mun[i] + mun[i + 1]) / 2
        mupim1by2 = (mup[i - 1] + mup[i]) / 2
        mupip1by2 = (mup[i] + mup[i + 1]) / 2
        ## Co-efficients for HOLE Continuity eqn
        ap[i] = mupim1by2 * Ber((fi[i] - fi[i - 1]) + (fi_p[i] - fi_p[i - 1]))
        cp[i] = mupip1by2 * Ber((fi[i] - fi[i + 1]) + (fi_p[i] - fi_p[i + 1]))
        bp[i] = -(
            mupim1by2 * Ber((fi[i - 1] - fi[i]) + (fi_p[i - 1] - fi_p[i]))
            + mupip1by2 * Ber((fi[i + 1] - fi[i]) + (fi_p[i + 1] - fi_p[i]))
        )
        ## Co-efficients for ELECTRON Continuity eqn
        an[i] = munim1by2 * Ber((fi[i - 1] - fi[i]) + (fi_n[i - 1] - fi_n[i]))
        cn[i] = munip1by2 * Ber((fi[i + 1] - fi[i]) + (fi_n[i + 1] - fi_n[i]))
        bn[i] = -(
            munim1by2 * Ber((fi[i] - fi[i - 1]) + (fi_n[i] - fi_n[i - 1]))
            + munip1by2 * Ber((fi[i] - fi[i + 1]) + (fi_n[i] - fi_n[i + 1]))
        )
        ## Forcing Function for ELECTRON and HOLE Continuity eqns
        fn[i] = (
            (dx2 / Vt)
            * (p[i] * n[i] - 1)
            / (TAUP0[i] * (n[i] + 1) + TAUN0[i] * (p[i] + 1))
        )
        fp[i] = (
            (dx2 / Vt)
            * (p[i] * n[i] - 1)
            / (TAUP0[i] * (n[i] + 1) + TAUN0[i] * (p[i] + 1))
        )
    # (C)  Start the iterative procedure for the solution of the linearized Continuity
    #     equation for "ELECTRONS" using LU decomposition method:
    dn[0] = bn[0]
    for i in range(1, n_max):
        betan[i] = an[i] / dn[i - 1]
        dn[i] = bn[i] - betan[i] * cn[i - 1]
    # Solution of Lv = f #
    vn[0] = fn[0]
    for i in range(1, n_max):
        vn[i] = fn[i] - betan[i] * vn[i - 1]
    # Solution of U*fi = v #
    tempn = vn[n_max - 1] / dn[n_max - 1]
    # deltan[n_max-1] = tempn - n[n_max-1]
    n[n_max - 1] = tempn
    for i in range(n_max - 2, -1, -1):  # delta#
        tempn = (vn[i] - cn[i] * n[i + 1]) / dn[i]
        #  deltan[i] = tempn - n[i]
        n[i] = tempn
    ####################### END of ELECTRON Continuty Solver ###########
    # (D)  Start the iterative procedure for the solution of the linearized Continuity
    #     equation for "HOLES" using LU decomposition method:
    dp[0] = bp[0]
    for i in range(1, n_max):
        betap[i] = ap[i] / dp[i - 1]
        dp[i] = bp[i] - betap[i] * cp[i - 1]
    # Solution of Lv = f #
    vp[0] = fp[0]
    for i in range(1, n_max):
        vp[i] = fp[i] - betap[i] * vp[i - 1]
    # Solution of U*fi = v #
    tempp = vp[n_max - 1] / dp[n_max - 1]
    # deltap[n_max-1] = tempp - p[n_max-1]
    p[n_max - 1] = tempp
    for i in range(n_max - 2, -1, -1):  # delta#
        tempp = (vp[i] - cp[i] * p[i + 1]) / dp[i]
        #   deltap[i] = tempp - p[i]
        p[i] = tempp
    ####################### END of HOLE Continuty Solver ###########
    return n, p


def Poisson_non_equi2(
    fi_stat,
    n,
    p,
    dop,
    Ppz_Psp,
    pol_surf_char,
    n_max,
    dx,
    fi,
    flag_conv_2,
    Ldi,
    ni,
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    iteration,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
):
    ####################################################################
    ## 3.3 Calculate potential fi again with new values of "n" and "p"##
    ##     and check for convergence                                  ##
    ####################################################################
    # Recalculate forcing function and central coefficient b for fi
    d = np.zeros(n_max)
    v = np.zeros(n_max)
    f = np.zeros(n_max)
    c = np.zeros(n_max)
    b = np.zeros(n_max)
    a = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = np.zeros(n_max)
    dop_out = np.zeros(n_max)
    Ppz_Psp_out = np.zeros(n_max)
    d = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = fi
    dop_out = dop / ni
    Ppz_Psp_out = Ppz_Psp / ni
    pol_surf_char_out = pol_surf_char / ni
    dx2 = dx * dx
    Ldi2 = Ldi * Ldi
    delta_acc = 1.0e-4
    if model.N_wells_virtual - 2 != 0 and 1 == 2:
        n, p, fi_non, EF = equi_np_fi3(
            fi_out,
            wfh_general,
            wfe_general,
            model,
            E_state_general,
            E_statec_general,
            meff_state_general,
            meff_statec_general,
            n_max,
            ni,
        )
    for i in range(1, n_max - 1):
        a[i] = Ldi2[i] / (dx2)
        c[i] = Ldi2[i] / (dx2)
        b[i] = -(2 * Ldi2[i] / (dx2) + n[i] + p[i])
        f[i] = n[i] - p[i] - dop_out[i] - Ppz_Psp_out[i] - (fi_out[i] * (n[i] + p[i]))
        # f[i] = n[i] - p[i] - dop_out[i]-(pol_surf_char_out[i+1]-pol_surf_char_out[i-1])/(2*dx) - (fi_out[i]*(n[i] + p[i]))
    a[0] = 0.0
    c[0] = 0.0
    b[0] = 1.0
    f[0] = fi_out[0]
    a[n_max - 1] = 0.0
    c[n_max - 1] = 0.0
    b[n_max - 1] = 1.0
    f[n_max - 1] = fi_out[n_max - 1]
    ## here values of n[i] and p[i] are used in place of exp(fi[i])
    # Solve for Updated potential given the new value of Forcing
    # Function using LU decomposition
    d[0] = b[0]
    for i in range(1, n_max):
        d[i] = b[i] - a[i] * c[i - 1] / d[i - 1]
    # Solution of Lv = f #
    v[0] = f[0]
    for i in range(1, n_max):
        v[i] = f[i] - a[i] * v[i - 1] / d[i - 1]
    # Solution of U*fi = v #
    temp = v[n_max - 1] / d[n_max - 1]
    delta[n_max - 1] = temp - fi_out[n_max - 1]
    fi_out[n_max - 1] = temp
    for i in range(n_max - 2, -1, -1):  # delta#
        temp = (v[i] - c[i] * fi_out[i + 1]) / d[i]
        delta[i] = temp - fi_out[i]
        fi_out[i] = temp
    delta_max = 0
    delta_max = max(abs(delta[:]))
    # Test convergence and start the loop if necessary else increase
    # the applied potential
    # print ('delta_max= ',delta_max)
    if delta_max < delta_acc:
        flag_conv_2 = False
    else:
        fi_out0 = fi_out
        fi_out0 += 0.15 * delta
        fi_out = fi_out0
        """"""
        if model.N_wells_virtual - 2 != 0 and 1 == 2:
            n, p, fi_non, EF = equi_np_fi3(
                fi_out,
                wfh_general,
                wfe_general,
                model,
                E_state_general,
                E_statec_general,
                meff_state_general,
                meff_statec_general,
                n_max,
                ni,
            )
        for i in range(1, n_max - 1):
            b[i] = -(2 * Ldi2[i] / (dx2) + n[i] + p[i])
            f[i] = (
                n[i] - p[i] - dop_out[i] - Ppz_Psp_out[i] - (fi_out[i] * (n[i] + p[i]))
            )
            # f[i] = n[i] - p[i] - dop_out[i]-(pol_surf_char_out[i+1]-pol_surf_char_out[i-1])/(2*dx) - (fi_out[i]*(n[i] + p[i]))
    return fi_out, flag_conv_2


def equi_np_fi31(
    g,
    fi_stat31,
    n,
    p,
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    iteration,
    fi_old,
    Vt,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    dop,
    Ppz_Psp,
    n_max,
    ni,
):  # use
    EF = 0.0
    nf = n * ni
    pf = p * ni
    E_statec_general_pc = E_statec_general
    E_state_general_pc = E_state_general
    ## Calculate Quasi Fermi Level - Efn Efp
    Ec = np.zeros(n_max)
    Ev = np.zeros(n_max)
    Ei = np.zeros(n_max)
    Efn = np.zeros(n_max)
    Efp = np.zeros(n_max)
    for i in range(1, n_max - 1):
        Ec[i] = fi_e[i] / q - Vt * fi_old[i]  # Values from the second Node%
        Ev[i] = fi_h[i] / q - Vt * fi_old[i]  # Values from the second Node%
    Ec[0] = Ec[1]
    Ec[n_max - 1] = Ec[n_max - 2]
    Ev[0] = Ev[1]
    Ev[n_max - 1] = Ev[n_max - 2]
    """
    if not(config.quantum_effect) and 1==2:
        Ef_Ec=np.zeros(n_max)
        Ev_Ef=np.zeros(n_max)
        E3kbT=np.zeros(n_max)
        #Determination of the Fermi Level
        EF=0.0
        for i1 in range(0,n_max):
            Ef_Ec[i1]=(EF-(fi_e[i1]-Vt*q*fi_old[i1]))/(kb*T)#*J2meV
            Ev_Ef[i1]=((fi_h[i1]-Vt*q*fi_old[i1])-EF)/(kb*T)#*J2meV
            E3kbT[i1]=-3*kb*T/(kb*T)#*J2meV
            if (EF*meV2J-(fi_e[i1]-Vt*q*fi_old[i1])>-3*kb*T  ) :#-3*kb*T 
                n[i1]=Nc[i1]*fd3((EF*meV2J-(fi_e[i1]-Vt*q*fi_old[i1]))/(kb*T))/ni[i1]
            else:
                n[i1]=Nc[i1]*exp((EF*meV2J-(fi_e[i1]-Vt*q*fi_old[i1]))/(kb*T))/ni[i1]
                #print(n[i1])
            if ((fi_h[i1]-Vt*q*fi_old[i1])-EF*meV2J>-3*kb*T ):#-3*kb*T 
                
                p[i1]=Nv[i1]*fd3(((fi_h[i1]-Vt*q*fi_old[i1])-EF*meV2J)/(kb*T))/ni[i1]
            else:
                p[i1]=Nv[i1]*exp(((fi_h[i1]-Vt*q*fi_old[i1])-EF*meV2J)/(kb*T))/ni[i1]        
        #fi_stat31=Ei=Efn=Efp
        
    else:    
        for i in range(0,n_max):
    """
    Ei[i] = Ec[i] - ((fi_e[i] - fi_h[i]) / (2 * q))
    Delta_fi = np.zeros(n_max)
    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        n[I1:I2] = 0.0
        p[I1:I2] = 0.0
        for j in range(0, model.subnumber_e, 1):
            for i in range(I1, I2):
                Efn[i] = Ei[i] + Vt * log(abs((nf[i] + 1e-1) / ni[i]))
                Delta_fi[i] = -Vt * q * fi_stat31[i] - (-Vt * q * fi_old[i])
                n[i] += (
                    (
                        fd2(
                            E_statec_general[k, j] - Delta_fi[i] * J2meV,
                            Efn[i] * q * J2meV,
                            model,
                        )
                        * meff_statec_general[k, j]
                        / (hbar ** 2 * pi)
                    )
                    * (wfe_general[k, j, i - I1]) ** 2
                    / (ni[i] * model.dx)
                )
        for jj in range(0, model.subnumber_h, 1):
            for ii in range(I1, I2):
                # print(Ei[ii],' -', Vt,'*log(abs(',pf[ii],'/',ni[ii],'))')
                Efp[ii] = Ei[ii] - Vt * log(abs((pf[ii] + 1e-1) / ni[ii]))
                Delta_fi[ii] = -Vt * q * fi_stat31[ii] - (
                    -Vt * q * fi_old[ii]
                )  # predictor-corrector-type approach.
                # print(Delta_fi[ii],'=', -Vt*q*fi_stat31[ii],'-',(-Vt*q*fi_old[ii]))
                p[ii] += (
                    (
                        fd1(
                            E_state_general[k, jj] - Delta_fi[ii] * J2meV,
                            Efp[ii] * q * J2meV,
                            model,
                        )
                        * meff_state_general[k, jj]
                        / (hbar ** 2 * pi)
                    )
                    * (wfh_general[k, jj, ii - I1]) ** 2
                    / (ni[ii] * model.dx)
                )
    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        for j in range(0, model.subnumber_e, 1):
            E_statec_general_pc[k, j] = (
                E_statec_general[k, j] - Delta_fi[I1:I2].mean() * J2meV
            )
        for jj in range(0, model.subnumber_h, 1):
            E_state_general_pc[k, jj] = (
                E_state_general[k, jj] - Delta_fi[I1:I2].mean() * J2meV
            )
        # print(Delta_fi[I1:I2])
        # print(fi_old[I1:I2])
    xaxis = np.arange(0, n_max) * model.dx
    span = np.ones(100000000)
    pl.plot(
        xaxis, Ev * 1e3, "k", xaxis, Ec * 1e3, "k"
    )  # ,xaxis,Efn*1e3,'r',xaxis,Efp*1e3,'b')
    for j in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(j, model.Well_boundary, model.n_max)
        i1 = I1 - I1
        i2 = I2 - I1
        for levelc, statec in zip(E_statec_general_pc[j, :], wfe_general[j, :, :]):
            # pl.axhline(levelc,0.1,0.9, hold=True,color='g',ls='--')
            pl.plot(
                xaxis[I1:I2],
                statec[i1:i2] * config.wavefunction_scalefactor + levelc,
                "b",
            )
            pl.plot(xaxis[I1:I2], levelc * span[I1:I2], "g", ls="--")
        for level, state in zip(E_state_general_pc[j, :], wfh_general[j, :, :]):
            # pl.axhline(level,0.1,0.9,color='g',ls='--')
            pl.plot(
                xaxis[I1:I2],
                state[i1:i2] * config.wavefunction_scalefactor + level,
                "b",
            )
            pl.plot(xaxis[I1:I2], level * span[I1:I2], "g", ls="--")
        # pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
    # pl.plot(xaxis,EF*span[0:model.n_max],'r',ls='--')
    # pl.axhline(E_F,0.1,0.9,color='r',ls='--')
    pl.xlabel("Position (m)")
    pl.ylabel("Energy (meV)")
    pl.grid(True)
    """"""
    return n, p, fi_old, EF  # density of carriers


def equi_np_fi22(
    vindex,
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    iteration,
    fi_old,
    Vt,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    dop,
    Ppz_Psp,
    n_max,
    ni,
    n,
    p,
):  # use
    n_q = np.zeros(n_max)
    p_q = np.zeros(n_max)
    Ec = np.zeros(n_max)
    Ev = np.zeros(n_max)
    Ei = np.zeros(n_max)
    Efn = np.zeros(n_max)
    Efp = np.zeros(n_max)
    fi_n = np.zeros(n_max)
    fi_p = np.zeros(n_max)
    for i in range(1, n_max - 1):
        Ec[i] = fi_e[i] / q - Vt * fi_old[i]  # Values from the second Node%
        Ev[i] = fi_h[i] / q - Vt * fi_old[i]  # Values from the second Node%
    Ec[0] = Ec[1]
    Ec[n_max - 1] = Ec[n_max - 2]
    Ev[0] = Ev[1]
    Ev[n_max - 1] = Ev[n_max - 2]
    for i in range(0, n_max):
        Ei[i] = Ec[i] - ((fi_e[i] - fi_h[i]) / (2 * q))
    #######################################################
    # EF=0.0

    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        # n[I1:I2]=0.0
        # p[I1:I2]=0.0
        # couter=0
        for j in range(0, model.subnumber_e, 1):
            for i in range(I11, I22):
                Efn[i] = Ei[i] + Vt * log(abs(n[i] + 1))
                n_q[i] += (
                    (
                        fd2(E_statec_general[k, j], Efn[i] * q * J2meV, model)
                        * meff_statec_general[k, j]
                        / (hbar ** 2 * pi)
                    )
                    * (wfe_general[k, j, i - I1]) ** 2
                    / (ni[i] * model.dx)
                )
        for jj in range(0, model.subnumber_h, 1):
            for ii in range(I11, I22):
                Efp[ii] = Ei[ii] - Vt * log(abs(p[ii] + 1))
                p_q[ii] += (
                    (
                        fd1(E_state_general[k, jj], Efp[ii] * q * J2meV, model)
                        * meff_state_general[k, jj]
                        / (hbar ** 2 * pi)
                    )
                    * (wfh_general[k, jj, ii - I1]) ** 2
                    / (ni[ii] * model.dx)
                )
    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        for i in range(I11, I22):
            fi_n[i] = Vt * log(abs((n_q[i] + 1) / (n[i] + 1)))
            fi_p[i] = -Vt * log(abs((p_q[i] + 1) / (p[i] + 1)))
    if vindex == 1000000:
        # print(max(n[:])/1e18,max(p[:])/1e18)
        xaxis = np.arange(0, n_max) * model.dx
        span = np.ones(100000000)
        pl.plot(
            xaxis,
            Ev * 1e3,
            "k",
            xaxis,
            Ec * 1e3,
            "k",
            xaxis,
            Efn * 1e3,
            "r",
            xaxis,
            Efp * 1e3,
            "b",
        )
        for j in range(1, model.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(j, model.Well_boundary, model.n_max)
            i1 = I1 - I1
            i2 = I2 - I1
            for levelc, statec in zip(E_statec_general[j, :], wfe_general[j, :, :]):
                # pl.axhline(levelc,0.1,0.9, hold=True,color='g',ls='--')
                pl.plot(
                    xaxis[I1:I2],
                    statec[i1:i2] * config.wavefunction_scalefactor + levelc,
                    "b",
                )
                pl.plot(xaxis[I1:I2], levelc * span[I1:I2], "g", ls="--")
            for level, state in zip(E_state_general[j, :], wfh_general[j, :, :]):
                # pl.axhline(level,0.1,0.9,color='g',ls='--')
                pl.plot(
                    xaxis[I1:I2],
                    state[i1:i2] * config.wavefunction_scalefactor + level,
                    "b",
                )
                pl.plot(xaxis[I1:I2], level * span[I1:I2], "g", ls="--")
            # pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
        # pl.plot(xaxis,EF*span[0:model.n_max],'r',ls='--')
        # pl.axhline(E_F,0.1,0.9,color='r',ls='--')
        pl.xlabel("Position (m)")
        pl.ylabel("Energy (meV)")
        pl.grid(True)
    return n_q, p_q, fi_n, fi_p  # density of carriers


def equi_np_fi222(
    ni,
    idata,
    fi_e,
    fi_h,
    fi_old,
    Vt,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
    n_max,
    n,
    p,
):  # use
    n_q = np.zeros(n_max)
    p_q = np.zeros(n_max)
    Ec = np.zeros(n_max)
    Ev = np.zeros(n_max)
    Ei = np.zeros(n_max)
    Efn = np.zeros(n_max)
    Efp = np.zeros(n_max)
    fi_n = np.zeros(n_max)
    fi_p = np.zeros(n_max)
    for i in range(1, n_max - 1):
        Ec[i] = fi_e[i] / q - Vt * fi_old[i]  # Values from the second Node%
        Ev[i] = fi_h[i] / q - Vt * fi_old[i]  # Values from the second Node%
    Ec[0] = Ec[1]
    Ec[n_max - 1] = Ec[n_max - 2]
    Ev[0] = Ev[1]
    Ev[n_max - 1] = Ev[n_max - 2]
    for i in range(0, n_max):
        Ei[i] = Ec[i] - ((fi_e[i] - fi_h[i]) / (2 * q))
    #######################################################
    # EF=0.0

    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        # n[I1:I2]=0.0
        # p[I1:I2]=0.0
        # couter=0
        for j in range(0, model.subnumber_e, 1):
            for i in range(I11, I22):
                Efn[i] = Ei[i] + Vt * log(abs(n[i] + 1))
                n_q[i] += (
                    (
                        fd2(E_statec_general[k, j], Efn[i] * q * J2meV, model)
                        * meff_statec_general[k, j]
                        / (hbar ** 2 * pi)
                    )
                    * (wfe_general[k, j, i - I1]) ** 2
                    / (ni[i] * model.dx)
                )
        for jj in range(0, model.subnumber_h, 1):
            for ii in range(I11, I22):
                Efp[ii] = Ei[ii] - Vt * log(abs(p[ii] + 1))
                p_q[ii] += (
                    (
                        fd1(E_state_general[k, jj], Efp[ii] * q * J2meV, model)
                        * meff_state_general[k, jj]
                        / (hbar ** 2 * pi)
                    )
                    * (wfh_general[k, jj, ii - I1]) ** 2
                    / (ni[ii] * model.dx)
                )
    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, model.Well_boundary, n_max)
        for i in range(I11, I22):
            fi_n[i] = Vt * log(abs((n_q[i] + 1) / (n[i] + 1)))
            fi_p[i] = -Vt * log(abs((p_q[i] + 1) / (p[i] + 1)))
    return fi_n, fi_p  # density of carriers


def Poisson_non_equi3(
    vindex,
    fi_stat3,
    n,
    p,
    dop,
    Ppz_Psp,
    pol_surf_char,
    n_max,
    dx,
    fi,
    flag_conv_2,
    Ldi,
    ni,
    fitotc,
    fitot,
    Nc,
    Nv,
    fi_e,
    fi_h,
    iteration,
    wfh_general,
    wfe_general,
    model,
    E_state_general,
    E_statec_general,
    meff_state_general,
    meff_statec_general,
):
    ####################################################################
    ## 3.3 Calculate potential fi again with new values of "n" and "p"##
    ##     and check for convergence                                  ##
    ####################################################################
    # Recalculate forcing function and central coefficient b for fi
    d = np.zeros(n_max)
    v = np.zeros(n_max)
    f = np.zeros(n_max)
    c = np.zeros(n_max)
    b = np.zeros(n_max)
    a = np.zeros(n_max)
    n_q = np.zeros(n_max)
    p_q = np.zeros(n_max)
    fi_n = np.zeros(n_max)
    fi_p = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = np.zeros(n_max)
    dop_out = np.zeros(n_max)
    Ppz_Psp_out = np.zeros(n_max)
    d = np.zeros(n_max)
    delta = np.zeros(n_max)
    fi_out = fi
    fi_stat30 = np.zeros(n_max)
    fi_stat30 = fi_stat3
    dop_out = dop / ni
    Ppz_Psp_out = Ppz_Psp / ni
    pol_surf_char_out = pol_surf_char / ni
    dx2 = dx * dx
    Ldi2 = Ldi * Ldi
    delta_acc = 1.0e-4
    if model.N_wells_virtual - 2 != 0 and config.quantum_effect:
        n_q, p_q, fi_n, fi_p = equi_np_fi22(
            vindex,
            fitotc,
            fitot,
            Nc,
            Nv,
            fi_e,
            fi_h,
            iteration,
            fi_out,
            Vt,
            wfh_general,
            wfe_general,
            model,
            E_state_general,
            E_statec_general,
            meff_state_general,
            meff_statec_general,
            dop,
            Ppz_Psp,
            n_max,
            ni,
            n,
            p,
        )
    # n,p,fi_non,EF =equi_np_fi31(g,fi_stat30,n,p,fitotc,fitot,Nc,Nv,fi_e,fi_h,iteration,fi_out,Vt,wfh_general,wfe_general,model,E_state_general,E_statec_general,meff_state_general,meff_statec_general,dop,Ppz_Psp,n_max,ni)
    for i in range(1, n_max - 1):
        a[i] = Ldi2[i] / (dx2)
        c[i] = Ldi2[i] / (dx2)
        b[i] = -(2 * Ldi2[i] / (dx2) + n[i] + p[i])
        f[i] = n[i] - p[i] - dop_out[i] - Ppz_Psp_out[i] - (fi_out[i] * (n[i] + p[i]))
        # f[i] = n[i] - p[i] - dop_out[i]-(pol_surf_char_out[i+1]-pol_surf_char_out[i-1])/(2*dx) - (fi_out[i]*(n[i] + p[i]))
    a[0] = 0.0
    c[0] = 0.0
    b[0] = 1.0
    f[0] = fi_out[0]
    a[n_max - 1] = 0.0
    c[n_max - 1] = 0.0
    b[n_max - 1] = 1.0
    f[n_max - 1] = fi_out[n_max - 1]
    ## here values of n[i] and p[i] are used in place of exp(fi[i])
    # Solve for Updated potential given the new value of Forcing
    # Function using LU decomposition
    # flag_conv =True  # convergence of the Poisson loop
    k_iter = 0
    # while flag_conv:
    k_iter = k_iter + 1
    d[0] = b[0]
    for i in range(1, n_max):
        d[i] = b[i] - a[i] * c[i - 1] / d[i - 1]
    # Solution of Lv = f #
    v[0] = f[0]
    for i in range(1, n_max):
        v[i] = f[i] - a[i] * v[i - 1] / d[i - 1]
    # Solution of U*fi = v #
    temp = v[n_max - 1] / d[n_max - 1]
    delta[n_max - 1] = temp - fi_out[n_max - 1]
    fi_out[n_max - 1] = temp
    for i in range(n_max - 2, -1, -1):  # delta#
        temp = (v[i] - c[i] * fi_out[i + 1]) / d[i]
        delta[i] = temp - fi_out[i]
        fi_out[i] = temp
    delta_max = 0
    delta_max = max(abs(delta[:]))
    # Test convergence and start the loop if necessary else increase
    # the applied potential
    print("delta_max= ", delta_max)
    # print(fi_stat30)
    if delta_max < delta_acc:
        flag_conv_2 = False
    else:
        if model.N_wells_virtual - 2 != 0 and config.quantum_effect:
            n_q, p_q, fi_n, fi_p = equi_np_fi22(
                vindex,
                fitotc,
                fitot,
                Nc,
                Nv,
                fi_e,
                fi_h,
                iteration,
                fi_out,
                Vt,
                wfh_general,
                wfe_general,
                model,
                E_state_general,
                E_statec_general,
                meff_state_general,
                meff_statec_general,
                dop,
                Ppz_Psp,
                n_max,
                ni,
                n,
                p,
            )
        # n,p,fi_non,EF =equi_np_fi31(g,fi_stat30,n,p,fitotc,fitot,Nc,Nv,fi_e,fi_h,iteration,fi_out,Vt,wfh_general,wfe_general,model,E_state_general,E_statec_general,meff_state_general,meff_statec_general,dop,Ppz_Psp,n_max,ni)
        for i in range(1, n_max - 1):
            b[i] = -(2 * Ldi2[i] / (dx2) + n[i] + p[i])
            f[i] = (
                n[i] - p[i] - dop_out[i] - Ppz_Psp_out[i] - (fi_out[i] * (n[i] + p[i]))
            )
            # f[i] = n[i] - p[i] - dop_out[i+1]-(pol_surf_char_out[i+1]-pol_surf_char_out[i-1])/(2*dx) - (fi_out[i]*(n[i] + p[i]))
    return fi_out, flag_conv_2, n_q, p_q, fi_n, fi_p


def Current2(
    vindex,
    n,
    p,
    mun,
    mup,
    fi,
    Vt,
    n_max,
    Total_Steps,
    q,
    dx,
    ni,
    Ldi,
    Jnip1by2,
    Jnim1by2,
    Jelec,
    Jpip1by2,
    Jpim1by2,
    Jhole,
):
    ##########################################################################
    ##                        CALCULATE CURRENT                             ##
    ##########################################################################
    for i in range(1, n_max - 1):

        # Electron Current
        Jnip1by2[vindex, i] = (
            (q * mun[i] * Vt / (dx))
            * ni[i]
            * (n[i + 1] * Ber((fi[i + 1] - fi[i])) - n[i] * Ber((fi[i] - fi[i + 1])))
        )
        Jnim1by2[vindex, i] = (
            (q * mun[i] * Vt / (dx))
            * ni[i]
            * (n[i] * Ber((fi[i] - fi[i - 1])) - n[i - 1] * Ber((fi[i - 1] - fi[i])))
        )
        Jelec[vindex, i] = (Jnip1by2[vindex, i] + Jnim1by2[vindex, i]) / 2
        # Hole Current
        Jpip1by2[vindex, i] = (
            (q * mup[i] * Vt / (dx))
            * ni[i]
            * (p[i + 1] * Ber((fi[i] - fi[i + 1])) - p[i] * Ber((fi[i + 1] - fi[i])))
        )
        Jpim1by2[vindex, i] = (
            (q * mup[i] * Vt / (dx))
            * ni[i]
            * (p[i] * Ber((fi[i - 1] - fi[i])) - p[i - 1] * Ber((fi[i] - fi[i - 1])))
        )
        Jhole[vindex, i] = (Jpip1by2[vindex, i] + Jpim1by2[vindex, i]) / 2
    ##         Jtotal(vindex) = Jelec
    return Jnip1by2, Jnim1by2, Jelec, Jpip1by2, Jpim1by2, Jhole


def Write_results_non_equi2(
    Nc,
    Nv,
    fi_e,
    fi_h,
    Vt,
    q,
    ni,
    n,
    p,
    dop,
    dx,
    Ldi,
    fi,
    n_max,
    Jnip1by2,
    Jnim1by2,
    Jelec,
    Jpip1by2,
    Jpim1by2,
    Jhole,
    Jtotal,
    Total_Steps,
):

    ro = np.zeros(n_max)
    el_field1 = np.zeros(n_max)
    el_field2 = np.zeros(n_max)

    Ec = np.zeros(n_max)
    Ev = np.zeros(n_max)
    Ei = np.zeros(n_max)
    Efn = np.zeros(n_max)
    Efp = np.zeros(n_max)
    av_curr = np.zeros(Total_Steps)
    for i in range(1, n_max - 1):
        Ec[i] = fi_e[i] / q - Vt * fi[i]  # Values from the second Node%
        Ev[i] = fi_h[i] / q - Vt * fi[i]  # Values from the second Node%
        ro[i] = -q * (ni[i] * n[i] - ni[i] * p[i] - dop[i])
        el_field1[i] = -(fi[i + 1] - fi[i]) * Vt / (dx)
        el_field2[i] = -(fi[i + 1] - fi[i - 1]) * Vt / (2 * dx)
    Jtotal[:, 0] = Jtotal[:, 1]
    Jelec[:, 0] = Jelec[:, 1]
    Jhole[:, 0] = Jhole[:, 1]
    Jtotal[:, n_max - 1] = Jtotal[:, n_max - 2]
    Jelec[:, n_max - 1] = Jelec[:, n_max - 2]
    Jhole[:, n_max - 1] = Jhole[:, n_max - 2]

    Ec[0] = Ec[1]
    Ec[n_max - 1] = Ec[n_max - 2]
    Ev[0] = Ev[1]
    Ev[n_max - 1] = Ev[n_max - 2]
    el_field1[0] = el_field1[1]
    el_field2[0] = el_field2[1]
    el_field1[n_max - 1] = el_field1[n_max - 2]
    el_field2[n_max - 1] = el_field2[n_max - 2]
    ro[0] = ro[1]
    ro[n_max - 1] = ro[n_max - 2]
    n[0] = n[1]
    n[n_max - 1] = n[n_max - 2]
    p[0] = p[1]
    p[n_max - 1] = p[n_max - 2]
    nf = n * ni
    pf = p * ni
    ## Calculate Quasi Fermi Level - Efn Efp
    for i in range(0, n_max):
        Ei[i] = Ec[i] - ((fi_e[i] - fi_h[i]) / (2 * q))
        Efn[i] = Ei[i] + Vt * log(nf[i] / ni[i] + 1)
        Efp[i] = Ei[i] - Vt * log(pf[i] / ni[i] + 1)
    Efn[0] = Efn[1]
    Efn[n_max - 1] = Efn[n_max - 2]
    Efp[0] = Efp[1]
    Efp[n_max - 1] = Efp[n_max - 2]
    for j in range(0, Total_Steps):
        av_curr[j] = Jtotal[j, 0]
    Ec_result = np.zeros(n_max)
    Ev_result = np.zeros(n_max)
    Ei_result = np.zeros(n_max)
    Efp_result = np.zeros(n_max)
    Efn_result = np.zeros(n_max)
    fi_result = np.zeros(n_max)
    fi_result = fi[0:n_max]
    Efn_result = Efn[0:n_max]
    Efp_result = Efp[0:n_max]
    ro_result = np.zeros(n_max)
    el_field1_result = np.zeros(n_max)
    el_field2_result = np.zeros(n_max)
    nf_result = np.zeros(n_max)
    pf_result = np.zeros(n_max)
    Ec_result = Ec[0:n_max]
    Ev_result = Ev[0:n_max]
    Ei_result = Ei[0:n_max]
    ro_result = ro[0:n_max]
    el_field1_result = el_field1[0:n_max]
    el_field2_result = el_field2[0:n_max]
    nf_result = n[0:n_max] * ni[0:n_max]
    pf_result = p[0:n_max] * ni[0:n_max]
    return (
        fi_result,
        Efn_result,
        Efp_result,
        ro_result,
        el_field1_result,
        el_field2_result,
        nf_result,
        pf_result,
        Ec_result,
        Ev_result,
        Ei_result,
        av_curr,
    )
