#!/usr/bin/env python
# -*- coding: utf-8 -*-
Description = f'''
Aestimo 1D Schrodinger-Poisson Solver
Copyright (C) 2013-2022

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. This program is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public
License along with this program. See ~/COPYING file or
http://www.gnu.org/copyleft/gpl.txt

INFORMATION:
This is the effective mass calculator for conduction band and
3x3 k.p Numpy calculator for valence band calculations.

Usage:
$ ./aestimo.py  <args>
'''
import os, sys, getopt, time, logging
import importlib.util # Fixed import
import matplotlib.pyplot as pl
import numpy as np
from math import log, exp, sqrt
from scipy import linalg
from argparse import ArgumentParser, HelpFormatter
from pathlib import Path
import textwrap
from aeslibs.VBHM import qsv, VBMAT1, VBMAT2, VBMAT_V, CBMAT, CBMAT_V, VBMAT_V_2
import config, database
from aeslibs.plotting import save_and_plot, save_and_plot2 # Fixed import
# Import aestimo_numpy_h lazily to avoid import conflicts
from aeslibs.aestimo_poisson1d import (
    Poisson_equi2,
    equi_np_fi,
    Write_results_equi2,
    equi_np_fi2,
    equi_np_fi3,
    Poisson_non_equi3,
    Poisson_equi_non_2,
    equi_np_fi22,
    equi_np_fi222,
)
from aeslibs.aestimo_poisson1d import (
    Poisson_equi1,
    Mobility2,
    Continuity2,
    Mobility3,
    Continuity3,
    Poisson_non_equi2,
    Current2,
    Write_results_non_equi2,
    Write_results_equi1,
    amort_wave,
)
from aeslibs.ddggummelmap import DDGgummelmap, apply_photovoltaic_BCs
from aeslibs.ddnnewtonmap import DDNnewtonmap
from aeslibs.func_lib import Ubernoulli
from aeslibs.ddgnlpoisson import DDGnlpoisson_new

time0 = time.time()  # timing audit

# Because alen is not used anymore
def alen(x):
    return 1 if np.isscalar(x) else len(x)

#alen = np.alen 

# Version
__version__ = "4.0.0"

def main():
    import runpy
    runpy.run_module('aestimo', run_name='__main__')

drawFigures = False

# Logger
logger = logging.getLogger('aestimo')

def initialize_logger():
    hdlr = logging.FileHandler(os.path.abspath(os.path.join(output_directory, 'aestimo.log')))
    hdlr_formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s')
    hdlr.setFormatter(hdlr_formatter)
    logger.addHandler(hdlr)
    # stderr
    ch = logging.StreamHandler()
    ch_formatter = logging.Formatter('%(levelname)s %(message)s')
    ch.setFormatter(ch_formatter)
    logger.addHandler(ch)
    # LOG level can be INFO, WARNING, ERROR
    logger.setLevel(logging.INFO)



# Defining constants and material parameters
q = 1.602176e-19  # C
kb = 1.3806504e-23  # J/K
hbar = 1.054588757e-34  # Js
m_e = 9.1093826e-31  # kg
pi = np.pi
eps0 = 8.8541878176e-12  # F/m
# TEMPERATURE
T = 300.0  # Kelvin
Vt = kb * T / q  # [eV]
J2meV = 1e3 / q  # Joules to meV
meV2J = 1e-3 * q  # meV to Joules

time1 = time.time()  # timing audit

# To print Description variable with argparse
class RawFormatter(HelpFormatter):
    def _fill_text(self, text, width, indent):
        return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])


def round2int(x):
    """int is sensitive to floating point numerical errors near whole numbers,
    this moves the discontinuity to the half interval. It is also equivalent
    to the normal rules for rounding positive numbers."""
    # int(x + (x>0) -0.5) # round2int for positive and negative numbers
    return int(x + 0.5)


def vegard1(first, second, mole):
    return first * mole + second * (1 - mole)


def get_varshni_Eg(matprops, T):
    """Calculates Eg at temperature T using Varshni's law from Eg at 300K."""
    Eg_300 = matprops.get("Eg", 0.0)
    alpha = matprops.get("alpha_var", 0.0)
    beta = matprops.get("beta_var", 0.0)
    if alpha == 0 or Eg_300 == 0:
        return Eg_300
    
    def dEg(temp):
        return -(alpha * temp ** 2) / (temp + beta)
    
    # Eg(T) = Eg(0) + dEg(T)
    # Eg(300) = Eg(0) + dEg(300) => Eg(0) = Eg(300) - dEg(300)
    # Eg(T) = Eg(300) - dEg(300) + dEg(T)
    return Eg_300 - dEg(300.0) + dEg(float(T))


class Structure:
    def __init__(self, database, **kwargs):
        """This class holds details on fthe structure to be simulated.
        database is the module containing the material properties. Then
        this class should have the following attributes set
        Fapp - applied field (Vm**-1)
        T - Temperature (K)
        subnumber_e - number of subbands to look for.
        comp_scheme - computing scheme
        dx - grid step size (m)
        n_max - number of grid points
        
        cb_meff #conduction band effective mass (kg) (array, len n_max)
        cb_meff_alpha #non-parabolicity constant.
        fi #Bandstructure potential (J) (array, len n_max)
        eps #dielectric constant (including eps0) (array, len n_max)
        dop #doping distribution (m**-3) (array, len n_max)
        
        These last 4 can be created by using the method
        create_structure_arrays(material_list)
        """
        # setting any parameters provided with initialisation
        for key, value in kwargs.items():
            setattr(self, key, value)
        # Loading materials database
        self.material_property = database.materialproperty
        totalmaterial = alen(self.material_property)

        self.alloy_property = database.alloyproperty
        totalalloy = alen(self.alloy_property)

        self.alloy_property_4 = database.alloyproperty4
        totalalloy += alen(self.alloy_property_4)

        if not hasattr(self, 'mat_crys_strc'):
             self.mat_crys_strc = getattr(self, 'mat_type', 'Zincblende')
        
        if not hasattr(self, 'material'):
             # If material list is missing, we assume arrays are provided manually
             # We just need to ensure standard attributes exist
             self.n_max = getattr(self, 'n_max', 0)
             if hasattr(self, 'fi') and not hasattr(self, 'fi_e'):
                 self.fi_e = self.fi
             return

        self.create_structure_arrays()

    def create_structure_arrays(self):
        """ initialise arrays/lists for structure"""
        # self.N_wells_real0=sum(sum(np.char.count(self.material,'w')))
        self.N_wells_real0 = sum(
            [1 for layer in self.material if len(layer) > 6 and layer[6] == "w"]
        )
        self.N_layers_real0 = len(
            self.material
        )  # sum(np.char.count([layer[6] for layer in self.material],'w'))+sum(np.char.count([layer[6] for layer in self.material],'b'))

        # Calculate the required number of grid points
        self.x_max = (
            sum([layer[0] for layer in self.material]) * 1e-9
        )  # total thickness (m)
        n_max = round2int(self.x_max / self.dx)
        # Check on n_max
        maxgridpoints = self.maxgridpoints
        mat_crys_strc = self.mat_crys_strc
        if n_max > maxgridpoints:
            logger.error("Grid number is exceeding the max number of %d", maxgridpoints)
            sys.exit()
        #
        self.n_max = n_max
        dx = self.dx
        material_property = self.material_property
        alloy_property = self.alloy_property
        alloy_property_4 = self.alloy_property_4
        cb_meff = np.zeros(n_max)  # conduction band effective mass
        cb_meff_alpha = np.zeros(n_max)  # non-parabolicity constant.
        m_hh = np.zeros(n_max)
        m_lh = np.zeros(n_max)
        m_so = np.zeros(n_max)
        # Elastic constants C11,C12
        C12 = np.zeros(n_max)
        C11 = np.zeros(n_max)
        # Elastic constants Wurtzite C13,C33
        C13 = np.zeros(n_max)
        C33 = np.zeros(n_max)
        C44 = np.zeros(n_max)
        # Spontaneous and Piezoelectric Polarizations constants D15,D13,D33 and Psp
        D15 = np.zeros(n_max)
        D31 = np.zeros(n_max)
        D33 = np.zeros(n_max)
        Psp = np.zeros(n_max)
        # Luttinger Parameters γ1,γ2,γ3
        GA3 = np.zeros(n_max)
        GA2 = np.zeros(n_max)
        GA1 = np.zeros(n_max)
        # Hole eff. mass parameter  Wurtzite Semiconductors
        A1 = np.zeros(n_max)
        A2 = np.zeros(n_max)
        A3 = np.zeros(n_max)
        A4 = np.zeros(n_max)
        A5 = np.zeros(n_max)
        A6 = np.zeros(n_max)
        # Lattice constant a0
        a0 = np.zeros(n_max)
        a0_wz = np.zeros(n_max)
        a0_sub = np.zeros(n_max)
        #  Deformation potentials ac,av,b
        Ac = np.zeros(n_max)
        Av = np.zeros(n_max)
        B = np.zeros(n_max)
        # Deformation potentials Wurtzite Semiconductors
        D1 = np.zeros(n_max)
        D2 = np.zeros(n_max)
        D3 = np.zeros(n_max)
        D4 = np.zeros(n_max)
        delta = np.zeros(n_max)  # delta splitt off
        delta_so = np.zeros(n_max)  # delta Spin–orbit split energy
        delta_cr = np.zeros(n_max)  # delta Crystal-field split energy
        # Strain related
        fi_h = np.zeros(n_max)  # Bandstructure potential
        fi_e = np.zeros(n_max)  # Bandstructure potential
        eps = np.zeros(n_max)  # dielectric constant
        dop = np.zeros(n_max)  # doping
        pol_surf_char = np.zeros(n_max)
        N_wells_real = 0
        N_wells_real2 = 0
        N_layers_real2 = 0
        N_wells_real0 = self.N_wells_real0
        N_layers_real0 = self.N_layers_real0
        N_wells_virtual = N_wells_real0 + 2
        N_wells_virtual2 = N_wells_real0 + 2
        N_layers_virtual = N_layers_real0 + 2
        Well_boundary = np.zeros((N_wells_virtual, 2), dtype=int)
        Well_boundary2 = np.zeros((N_wells_virtual, 2), dtype=int)
        barrier_boundary = np.zeros((N_wells_virtual + 1, 2), dtype=int)
        layer_boundary = np.zeros((N_layers_virtual, 2), dtype=int)
        n_max_general = np.zeros(N_wells_virtual, dtype=int)
        Well_boundary[N_wells_virtual - 1, 0] = n_max - 1
        Well_boundary[N_wells_virtual - 1, 1] = n_max - 1
        Well_boundary2[N_wells_virtual - 1, 0] = n_max - 1
        Well_boundary2[N_wells_virtual - 1, 1] = n_max - 1
        barrier_boundary[N_wells_virtual, 0] = n_max - 1
        barrier_len = np.zeros(N_wells_virtual + 1)
        n = np.zeros(n_max)
        p = np.zeros(n_max)
        TAUN0 = np.zeros(n_max)
        TAUP0 = np.zeros(n_max)
        mun0 = np.zeros(n_max)
        mup0 = np.zeros(n_max)

        Cn0 = np.zeros(n_max)
        Cp0 = np.zeros(n_max)
        BETAN = np.zeros(n_max)
        BETAP = np.zeros(n_max)
        VSATN = np.zeros(n_max)
        VSATP = np.zeros(n_max)
        position = 0.0  # keeping in nanometres (to minimise errors)
        for layer in self.material:
            startindex = round2int(position * 1e-9 / dx)
            z0 = round2int(position * 1e-9 / dx)
            position += layer[0]  # update position to end of the layer
            finishindex = round2int(position * 1e-9 / dx)
            z1 = round2int(position * 1e-9 / dx)
            #
            matType = layer[1]
            if matType in material_property:
                matprops = material_property[matType]
                cb_meff[startindex:finishindex] = matprops["m_e"] * m_e
                cb_meff_alpha[startindex:finishindex] = matprops["m_e_alpha"]
                Eg_T = get_varshni_Eg(matprops, self.T)
                fi_e[startindex:finishindex] = (
                    matprops["Band_offset"] * Eg_T * q
                )  # Joule
                is_wz_mat = "A1" in matprops
                if (mat_crys_strc == "Zincblende" or not is_wz_mat) and "a0_sub" in matprops:
                    a0_sub[startindex:finishindex] = matprops.get("a0_sub", 5.6533) * 1e-10
                    C11[startindex:finishindex] = matprops.get("C11", 11.879) * 1e10
                    C12[startindex:finishindex] = matprops.get("C12", 5.376) * 1e10
                    GA1[startindex:finishindex] = matprops.get("GA1", 6.8)
                    GA2[startindex:finishindex] = matprops.get("GA2", 1.9)
                    GA3[startindex:finishindex] = matprops.get("GA3", 2.73)
                    Ac[startindex:finishindex] = matprops.get("Ac", -7.17) * q
                    Av[startindex:finishindex] = matprops.get("Av", 1.16) * q
                    B[startindex:finishindex] = matprops.get("B", -1.7) * q
                    delta[startindex:finishindex] = matprops.get("delta", 0.28) * q
                    fi_h[startindex:finishindex] = (
                        -(1 - matprops["Band_offset"]) * Eg_T * q
                    )  # Joule
                    eps[startindex:finishindex] = matprops["epsilonStatic"] * eps0
                    a0[startindex:finishindex] = matprops.get("a0", 5.6533) * 1e-10
                    TAUN0[startindex:finishindex] = matprops.get("TAUN0", 1e-8)
                    TAUP0[startindex:finishindex] = matprops.get("TAUP0", 1e-8)
                    # Mobility is in m^2/Vs in database
                    mun0[startindex:finishindex] = matprops.get("mun0", 0.1)
                    mup0[startindex:finishindex] = matprops.get("mup0", 0.02)

                    Cn0[startindex:finishindex] = matprops.get("Cn0", 2.8e-31) * 1e-12
                    Cp0[startindex:finishindex] = matprops.get("Cp0", 2.8e-32) * 1e-12
                    BETAN[startindex:finishindex] = matprops.get("BETAN", 2.0)
                    BETAP[startindex:finishindex] = matprops.get("BETAP", 1.0)
                    VSATN[startindex:finishindex] = matprops.get("VSATN", 3e5)
                    VSATP[startindex:finishindex] = matprops.get("VSATP", 6e5)
                elif mat_crys_strc == "Wurtzite" and is_wz_mat:
                    a0_sub[startindex:finishindex] = matprops["a0_sub"] * 1e-10
                    C11[startindex:finishindex] = matprops["C11"] * 1e10
                    C12[startindex:finishindex] = matprops["C12"] * 1e10
                    C13[startindex:finishindex] = matprops["C13"] * 1e10
                    C33[startindex:finishindex] = matprops["C33"] * 1e10
                    A1[startindex:finishindex] = matprops["A1"]
                    A2[startindex:finishindex] = matprops["A2"]
                    A3[startindex:finishindex] = matprops["A3"]
                    A4[startindex:finishindex] = matprops["A4"]
                    A5[startindex:finishindex] = matprops["A5"]
                    A6[startindex:finishindex] = matprops["A6"]
                    Ac[startindex:finishindex] = matprops["Ac"] * q
                    D1[startindex:finishindex] = matprops["D1"] * q
                    D2[startindex:finishindex] = matprops["D2"] * q
                    D3[startindex:finishindex] = matprops["D3"] * q
                    D4[startindex:finishindex] = matprops["D4"] * q
                    D31[startindex:finishindex] = matprops["D31"]
                    D33[startindex:finishindex] = matprops["D33"]
                    a0_wz[startindex:finishindex] = matprops["a0_wz"] * 1e-10
                    delta_so[startindex:finishindex] = matprops["delta_so"] * q
                    delta_cr[startindex:finishindex] = matprops["delta_cr"] * q
                    eps[startindex:finishindex] = matprops["epsilonStatic"] * eps0
                    fi_h[startindex:finishindex] = (
                        -(1 - matprops["Band_offset"]) * Eg_T * q
                    )
                    # Mobility is in m^2/Vs in database
                    mun0[startindex:finishindex] = matprops["mun0"]
                    mup0[startindex:finishindex] = matprops["mup0"]
                    
                    # Apply global tau override if provided
                    if getattr(self, 'tau', None) is not None:
                        logger.info(f"Applying global tau override: {self.tau} s")
                        TAUN0[startindex:finishindex] = float(self.tau)
                        TAUP0[startindex:finishindex] = float(self.tau)

                    Cn0[startindex:finishindex] = matprops["Cn0"] * 1e-12
                    Cp0[startindex:finishindex] = matprops["Cp0"] * 1e-12
                    BETAN[startindex:finishindex] = matprops["BETAN"]
                    BETAP[startindex:finishindex] = matprops["BETAP"]
                    VSATN[startindex:finishindex] = matprops["VSATN"]
                    VSATP[startindex:finishindex] = matprops["VSATP"]
            elif matType in alloy_property:
                alloyprops = alloy_property[matType]
                mat1 = material_property[alloyprops["Material1"]]
                mat2 = material_property[alloyprops["Material2"]]
                x = layer[2]  # alloy ratio
                cb_meff_alloy = x * mat1["m_e"] + (1 - x) * mat2["m_e"]
                cb_meff[startindex:finishindex] = cb_meff_alloy * m_e
                Eg1 = get_varshni_Eg(mat1, self.T)
                Eg2 = get_varshni_Eg(mat2, self.T)
                Eg = (
                    x * Eg1
                    + (1 - x) * Eg2
                    - alloyprops["Bowing_param"] * x * (1 - x)
                )  # eV
                fi_e[startindex:finishindex] = (
                    alloyprops["Band_offset"] * Eg * q
                )  # for electron. Joule
                a0_sub[startindex:finishindex] = alloyprops["a0_sub"] * 1e-10
                # Apply global tau override if provided
                if getattr(self, 'tau', None) is not None:
                    TAUN0[startindex:finishindex] = float(self.tau)
                    TAUP0[startindex:finishindex] = float(self.tau)
                else:
                    TAUN0[startindex:finishindex] = alloyprops["TAUN0"]
                    TAUP0[startindex:finishindex] = alloyprops["TAUP0"]
                
                # Mobility is in m^2/Vs in database
                mun0[startindex:finishindex] = alloyprops["mun0"]
                mup0[startindex:finishindex] = alloyprops["mup0"]
                Cn0[startindex:finishindex] = alloyprops["Cn0"] * 1e-12
                Cp0[startindex:finishindex] = alloyprops["Cp0"] * 1e-12

                BETAN[startindex:finishindex] = alloyprops["BETAN"]
                BETAP[startindex:finishindex] = alloyprops["BETAP"]
                VSATN[startindex:finishindex] = alloyprops["VSATN"]
                VSATP[startindex:finishindex] = alloyprops["VSATP"]
                is_wz_alloy = ("A1" in mat1 and "A1" in mat2)
                if mat_crys_strc == "Zincblende" or not is_wz_alloy:
                    C11[startindex:finishindex] = (
                        x * mat1.get("C11", 11.879) + (1 - x) * mat2.get("C11", 11.879)
                    ) * 1e10
                    C12[startindex:finishindex] = (
                        x * mat1.get("C12", 5.376) + (1 - x) * mat2.get("C12", 5.376)
                    ) * 1e10
                    GA1[startindex:finishindex] = (
                        x * mat1.get("GA1", 6.8) + (1 - x) * mat2.get("GA1", 6.8)
                    )
                    GA2[startindex:finishindex] = (
                        x * mat1.get("GA2", 1.9) + (1 - x) * mat2.get("GA2", 1.9)
                    )
                    GA3[startindex:finishindex] = (
                        x * mat1.get("GA3", 2.73) + (1 - x) * mat2.get("GA3", 2.73)
                    )
                    Ac_alloy = x * mat1.get("Ac", -7.17) + (1 - x) * mat2.get("Ac", -7.17)
                    Ac[startindex:finishindex] = Ac_alloy * q
                    Av_alloy = x * mat1.get("Av", 1.16) + (1 - x) * mat2.get("Av", 1.16)
                    Av[startindex:finishindex] = Av_alloy * q
                    B_alloy = x * mat1.get("B", -1.7) + (1 - x) * mat2.get("B", -1.7)
                    B[startindex:finishindex] = B_alloy * q
                    delta_alloy = x * mat1.get("delta", 0.28) + (1 - x) * mat2.get("delta", 0.28)
                    delta[startindex:finishindex] = delta_alloy * q
                    fi_h[startindex:finishindex] = (
                        -(1 - alloyprops["Band_offset"]) * Eg * q
                    )  # -(-1.33*(1-x)-0.8*x)for electron. Joule-1.97793434e-20 #
                    eps[startindex:finishindex] = (
                        x * mat1["epsilonStatic"] + (1 - x) * mat2["epsilonStatic"]
                    ) * eps0
                    a0[startindex:finishindex] = (
                       x  * mat1.get("a0", 5.6533) + (1 - x) * mat2.get("a0", 5.6533)
                    ) * 1e-10
                    cb_meff_alpha[startindex:finishindex] = alloyprops.get("m_e_alpha", 0.0) * (
                        mat2["m_e"] / cb_meff_alloy
                    )  # non-parabolicity constant for alloy. THIS CALCULATION IS MOSTLY WRONG. MUST BE CONTROLLED. SBL

                    mun0[startindex:finishindex] = (
                        x * mat1.get("mun0", 0.1) + (1 - x) * mat2.get("mun0", 0.1)
                    )
                    mup0[startindex:finishindex] = (
                        x * mat1.get("mup0", 0.02) + (1 - x) * mat2.get("mup0", 0.02)
                    )

                    Cn0[startindex:finishindex] = (
                        x * mat1.get("Cn0", 2.8e-31) + (1 - x) * mat2.get("Cn0", 2.8e-31)
                    ) * 1e-12
                    Cp0[startindex:finishindex] = (
                        x * mat1.get("Cp0", 2.8e-32) + (1 - x) * mat2.get("Cp0", 2.8e-32)
                    ) * 1e-12
                elif mat_crys_strc == "Wurtzite" and is_wz_alloy:
                    # A1[startindex:finishindex] =vegard1(mat1['A1'],mat1['A1'],x)
                    A1[startindex:finishindex] = x * mat1["A1"] + (1 - x) * mat2["A1"]
                    A2[startindex:finishindex] = x * mat1["A2"] + (1 - x) * mat2["A2"]
                    A3[startindex:finishindex] = x * mat1["A3"] + (1 - x) * mat2["A3"]
                    A4[startindex:finishindex] = x * mat1["A4"] + (1 - x) * mat2["A4"]
                    A5[startindex:finishindex] = x * mat1["A5"] + (1 - x) * mat2["A5"]
                    A6[startindex:finishindex] = x * mat1["A6"] + (1 - x) * mat2["A6"]
                    D1[startindex:finishindex] = (
                        x * mat1["D1"] + (1 - x) * mat2["D1"]
                    ) * q
                    D2[startindex:finishindex] = (
                        x * mat1["D2"] + (1 - x) * mat2["D2"]
                    ) * q
                    D3[startindex:finishindex] = (
                        x * mat1["D3"] + (1 - x) * mat2["D3"]
                    ) * q
                    D4[startindex:finishindex] = (
                        x * mat1["D4"] + (1 - x) * mat2["D4"]
                    ) * q
                    C13[startindex:finishindex] = (
                        x * mat1["C13"] + (1 - x) * mat2["C13"]
                    ) * 1e10  # for newton/meter²
                    C33[startindex:finishindex] = (
                        x * mat1["C33"] + (1 - x) * mat2["C33"]
                    ) * 1e10
                    D31[startindex:finishindex] = (
                        x * mat1["D31"] + (1 - x) * mat2["D31"]
                    )
                    D33[startindex:finishindex] = (
                        x * mat1["D33"] + (1 - x) * mat2["D33"]
                    )
                    Psp[startindex:finishindex] = (
                        x * mat1["Psp"] + (1 - x) * mat2["Psp"]
                    )
                    C11[startindex:finishindex] = (
                        x * mat1["C11"] + (1 - x) * mat2["C11"]
                    ) * 1e10
                    C12[startindex:finishindex] = (
                        x * mat1["C12"] + (1 - x) * mat2["C12"]
                    ) * 1e10
                    a0_wz[startindex:finishindex] = (
                        x * mat1["a0_wz"] + (1 - x) * mat2["a0_wz"]
                    ) * 1e-10
                    eps[startindex:finishindex] = (
                        x * mat1["epsilonStatic"] + (1 - x) * mat2["epsilonStatic"]
                    ) * eps0
                    fi_h[startindex:finishindex] = (
                        -(1 - alloyprops["Band_offset"]) * Eg * q
                    )
                    delta_so[startindex:finishindex] = (
                        x * mat1["delta_so"] + (1 - x) * mat2["delta_so"]
                    ) * q
                    delta_cr[startindex:finishindex] = (
                        x * mat1["delta_cr"] + (1 - x) * mat2["delta_cr"]
                    ) * q
                    Ac_alloy = x * mat1["Ac"] + (1 - x) * mat2["Ac"]
                    Ac[startindex:finishindex] = Ac_alloy * q
                    mun0[startindex:finishindex] = (
                        x * mat1["mun0"] + (1 - x) * mat2["mun0"]
                    )
                    mup0[startindex:finishindex] = (
                        x * mat1["mup0"] + (1 - x) * mat2["mup0"]
                    )

                    Cn0[startindex:finishindex] = (
                        x * mat1["Cn0"] + (1 - x) * mat2["Cn0"]
                    ) * 1e-12
                    Cp0[startindex:finishindex] = (
                        x * mat1["Cp0"] + (1 - x) * mat2["Cp0"]
                    ) * 1e-12
                    #############################################
            elif matType in alloy_property_4:
                alloyprops = alloy_property_4[matType]
                TAUN0[startindex:finishindex] = alloyprops["TAUN0"]
                TAUP0[startindex:finishindex] = alloyprops["TAUP0"]
                BETAN[startindex:finishindex] = alloyprops["BETAN"]
                BETAP[startindex:finishindex] = alloyprops["BETAP"]
                VSATN[startindex:finishindex] = alloyprops["VSATN"]
                VSATP[startindex:finishindex] = alloyprops["VSATP"]
                if mat_crys_strc == "Zincblende":

                    alloyprops = alloy_property_4[matType]
                    mat1 = material_property[alloyprops["Material1"]]
                    mat2 = material_property[alloyprops["Material2"]]
                    mat3 = material_property[alloyprops["Material3"]]
                    mat4 = material_property[alloyprops["Material4"]]
                    # mat1:InAs
                    # mat2:GaAs
                    # mat3:InP
                    # mat4:GaP
                    # This is accourding to interpolated Vegard’s law for quaternary AxB(1-x)CyD(1-y)=InxGa(1-x)AsyP(1-y)
                    x = layer[2]  # alloy ratio x
                    y = layer[3]  # alloy ratio y
                    cb_meff_alloy_ABC_x = x * mat1["m_e"] + (1 - x) * mat2["m_e"]
                    cb_meff_alloy_ABD_x = x * mat3["m_e"] + (1 - x) * mat4["m_e"]
                    cb_meff_alloy_ACD_y = y * mat1["m_e"] + (1 - y) * mat3["m_e"]
                    cb_meff_alloy_BCD_y = y * mat2["m_e"] + (1 - y) * mat4["m_e"]
                    cb_meff_alloy = (
                        x
                        * (1 - x)
                        * (y * cb_meff_alloy_ABC_x + (1 - y) * cb_meff_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * cb_meff_alloy_ACD_y + (1 - x) * cb_meff_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))
                    cb_meff[startindex:finishindex] = cb_meff_alloy * m_e

                    Eg_alloy_ABC_x = (
                        x * mat1["Eg"]
                        + (1 - x) * mat2["Eg"]
                        - alloyprops["Bowing_param_ABC"] * x * (1 - x)
                    )  # eV InGaAs
                    Eg_alloy_ABD_x = (
                        x * mat3["Eg"]
                        + (1 - x) * mat4["Eg"]
                        - alloyprops["Bowing_param_ABD"] * x * (1 - x)
                    )  # eV InGaP
                    Eg_alloy_ACD_y = (
                        y * mat1["Eg"]
                        + (1 - y) * mat3["Eg"]
                        - alloyprops["Bowing_param_ACD"] * y * (1 - y)
                    )  # eV InAsP
                    Eg_alloy_BCD_y = (
                        y * mat2["Eg"]
                        + (1 - y) * mat4["Eg"]
                        - alloyprops["Bowing_param_BCD"] * y * (1 - y)
                    )  # eV GaAsP
                    Eg = (
                        x * (1 - x) * (y * Eg_alloy_ABC_x + (1 - y) * Eg_alloy_ABD_x)
                        + y * (1 - y) * (x * Eg_alloy_ACD_y + (1 - x) * Eg_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))

                    fi_e[startindex:finishindex] = (
                        alloyprops["Band_offset"] * Eg * q
                    )  # for electron. Joule
                    a0_sub[startindex:finishindex] = alloyprops["a0_sub"] * 1e-10
                    C11_alloy_ABC_x = x * mat1["C11"] + (1 - x) * mat2["C11"]
                    C11_alloy_ABD_x = x * mat3["C11"] + (1 - x) * mat4["C11"]
                    C11_alloy_ACD_y = y * mat1["C11"] + (1 - y) * mat3["C11"]
                    C11_alloy_BCD_y = y * mat2["C11"] + (1 - y) * mat4["C11"]
                    C11[startindex:finishindex] = (
                        (
                            x
                            * (1 - x)
                            * (y * C11_alloy_ABC_x + (1 - y) * C11_alloy_ABD_x)
                            + y
                            * (1 - y)
                            * (x * C11_alloy_ACD_y + (1 - x) * C11_alloy_BCD_y)
                        )
                        / (x * (1 - x) + y * (1 - y))
                    ) * 1e10

                    C12_alloy_ABC_x = x * mat1["C12"] + (1 - x) * mat2["C12"]
                    C12_alloy_ABD_x = x * mat3["C12"] + (1 - x) * mat4["C12"]
                    C12_alloy_ACD_y = y * mat1["C12"] + (1 - y) * mat3["C12"]
                    C12_alloy_BCD_y = y * mat2["C12"] + (1 - y) * mat4["C12"]
                    C12[startindex:finishindex] = (
                        (
                            x
                            * (1 - x)
                            * (y * C12_alloy_ABC_x + (1 - y) * C12_alloy_ABD_x)
                            + y
                            * (1 - y)
                            * (x * C12_alloy_ACD_y + (1 - x) * C12_alloy_BCD_y)
                        )
                        / (x * (1 - x) + y * (1 - y))
                    ) * 1e10

                    GA1_alloy_ABC_x = x * mat1["GA1"] + (1 - x) * mat2["GA1"]
                    GA1_alloy_ABD_x = x * mat3["GA1"] + (1 - x) * mat4["GA1"]
                    GA1_alloy_ACD_y = y * mat1["GA1"] + (1 - y) * mat3["GA1"]
                    GA1_alloy_BCD_y = y * mat2["GA1"] + (1 - y) * mat4["GA1"]
                    GA1[startindex:finishindex] = (
                        x * (1 - x) * (y * GA1_alloy_ABC_x + (1 - y) * GA1_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * GA1_alloy_ACD_y + (1 - x) * GA1_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))

                    GA2_alloy_ABC_x = x * mat1["GA2"] + (1 - x) * mat2["GA2"]
                    GA2_alloy_ABD_x = x * mat3["GA2"] + (1 - x) * mat4["GA2"]
                    GA2_alloy_ACD_y = y * mat1["GA2"] + (1 - y) * mat3["GA2"]
                    GA2_alloy_BCD_y = y * mat2["GA2"] + (1 - y) * mat4["GA2"]
                    GA2[startindex:finishindex] = (
                        x * (1 - x) * (y * GA2_alloy_ABC_x + (1 - y) * GA2_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * GA2_alloy_ACD_y + (1 - x) * GA2_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))

                    GA3_alloy_ABC_x = x * mat1["GA3"] + (1 - x) * mat2["GA3"]
                    GA3_alloy_ABD_x = x * mat3["GA3"] + (1 - x) * mat4["GA3"]
                    GA3_alloy_ACD_y = y * mat1["GA3"] + (1 - y) * mat3["GA3"]
                    GA3_alloy_BCD_y = y * mat2["GA3"] + (1 - y) * mat4["GA3"]
                    GA3[startindex:finishindex] = (
                        x * (1 - x) * (y * GA3_alloy_ABC_x + (1 - y) * GA3_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * GA3_alloy_ACD_y + (1 - x) * GA3_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))

                    Ac_alloy_ABC_x = x * mat1["Ac"] + (1 - x) * mat2["Ac"]
                    Ac_alloy_ABD_x = x * mat3["Ac"] + (1 - x) * mat4["Ac"]
                    Ac_alloy_ACD_y = y * mat1["Ac"] + (1 - y) * mat3["Ac"]
                    Ac_alloy_BCD_y = y * mat2["Ac"] + (1 - y) * mat4["Ac"]
                    Ac_alloy = (
                        x * (1 - x) * (y * Ac_alloy_ABC_x + (1 - y) * Ac_alloy_ABD_x)
                        + y * (1 - y) * (x * Ac_alloy_ACD_y + (1 - x) * Ac_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))
                    Ac[startindex:finishindex] = Ac_alloy * q

                    Av_alloy_ABC_x = x * mat1["Av"] + (1 - x) * mat2["Av"]
                    Av_alloy_ABD_x = x * mat3["Av"] + (1 - x) * mat4["Av"]
                    Av_alloy_ACD_y = y * mat1["Av"] + (1 - y) * mat3["Av"]
                    Av_alloy_BCD_y = y * mat2["Av"] + (1 - y) * mat4["Av"]
                    Av_alloy = (
                        x * (1 - x) * (y * Av_alloy_ABC_x + (1 - y) * Av_alloy_ABD_x)
                        + y * (1 - y) * (x * Av_alloy_ACD_y + (1 - x) * Av_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))
                    Av[startindex:finishindex] = Av_alloy * q

                    B_alloy_ABC_x = x * mat1["B"] + (1 - x) * mat2["B"]
                    B_alloy_ABD_x = x * mat3["B"] + (1 - x) * mat4["B"]
                    B_alloy_ACD_y = y * mat1["B"] + (1 - y) * mat3["B"]
                    B_alloy_BCD_y = y * mat2["B"] + (1 - y) * mat4["B"]
                    B_alloy = (
                        x * (1 - x) * (y * B_alloy_ABC_x + (1 - y) * B_alloy_ABD_x)
                        + y * (1 - y) * (x * B_alloy_ACD_y + (1 - x) * B_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))
                    B[startindex:finishindex] = B_alloy * q

                    delta_alloy_ABC_x = x * mat1["delta"] + (1 - x) * mat2["delta"]
                    delta_alloy_ABD_x = x * mat3["delta"] + (1 - x) * mat4["delta"]
                    delta_alloy_ACD_y = y * mat1["delta"] + (1 - y) * mat3["delta"]
                    delta_alloy_BCD_y = y * mat2["delta"] + (1 - y) * mat4["delta"]
                    delta_alloy = (
                        x
                        * (1 - x)
                        * (y * delta_alloy_ABC_x + (1 - y) * delta_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * delta_alloy_ACD_y + (1 - x) * delta_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))
                    delta[startindex:finishindex] = delta_alloy * q

                    fi_h[startindex:finishindex] = (
                        -(1 - alloyprops["Band_offset"]) * Eg * q
                    )  # -(-1.33*(1-x)-0.8*x)for electron. Joule-1.97793434e-20 #

                    eps_alloy_ABC_x = (
                        x * mat1["epsilonStatic"] + (1 - x) * mat2["epsilonStatic"]
                    )
                    eps_alloy_ABD_x = (
                        x * mat3["epsilonStatic"] + (1 - x) * mat4["epsilonStatic"]
                    )
                    eps_alloy_ACD_y = (
                        y * mat1["epsilonStatic"] + (1 - y) * mat3["epsilonStatic"]
                    )
                    eps_alloy_BCD_y = (
                        y * mat2["epsilonStatic"] + (1 - y) * mat4["epsilonStatic"]
                    )
                    eps_alloy = (
                        x * (1 - x) * (y * eps_alloy_ABC_x + (1 - y) * eps_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * eps_alloy_ACD_y + (1 - x) * eps_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))
                    eps[startindex:finishindex] = eps_alloy * eps0

                    a0_alloy_ABC_x = x * mat1["a0"] + (1 - x) * mat2["a0"]
                    a0_alloy_ABD_x = x * mat3["a0"] + (1 - x) * mat4["a0"]
                    a0_alloy_ACD_y = y * mat1["a0"] + (1 - y) * mat3["a0"]
                    a0_alloy_BCD_y = y * mat2["a0"] + (1 - y) * mat4["a0"]
                    a0_alloy = (
                        x * (1 - x) * (y * a0_alloy_ABC_x + (1 - y) * a0_alloy_ABD_x)
                        + y * (1 - y) * (x * a0_alloy_ACD_y + (1 - x) * a0_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))
                    a0[startindex:finishindex] = a0_alloy * 1e-10

                    mun0_alloy_ABC_x = x * mat1["mun0"] + (1 - x) * mat2["mun0"]
                    mun0_alloy_ABD_x = x * mat3["mun0"] + (1 - x) * mat4["mun0"]
                    mun0_alloy_ACD_y = y * mat1["mun0"] + (1 - y) * mat3["mun0"]
                    mun0_alloy_BCD_y = y * mat2["mun0"] + (1 - y) * mat4["mun0"]
                    mun0[startindex:finishindex] = (
                        x
                        * (1 - x)
                        * (y * mun0_alloy_ABC_x + (1 - y) * mun0_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * mun0_alloy_ACD_y + (1 - x) * mun0_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))

                    mup0_alloy_ABC_x = x * mat1["mup0"] + (1 - x) * mat2["mup0"]
                    mup0_alloy_ABD_x = x * mat3["mup0"] + (1 - x) * mat4["mup0"]
                    mup0_alloy_ACD_y = y * mat1["mup0"] + (1 - y) * mat3["mup0"]
                    mup0_alloy_BCD_y = y * mat2["mup0"] + (1 - y) * mat4["mup0"]
                    mup0[startindex:finishindex] = (
                        x
                        * (1 - x)
                        * (y * mup0_alloy_ABC_x + (1 - y) * mup0_alloy_ABD_x)
                        + y
                        * (1 - y)
                        * (x * mup0_alloy_ACD_y + (1 - x) * mup0_alloy_BCD_y)
                    ) / (x * (1 - x) + y * (1 - y))

                    Cn0_alloy_ABC_x = x * mat1["Cn0"] + (1 - x) * mat2["Cn0"]
                    Cn0_alloy_ABD_x = x * mat3["Cn0"] + (1 - x) * mat4["Cn0"]
                    Cn0_alloy_ACD_y = y * mat1["Cn0"] + (1 - y) * mat3["Cn0"]
                    Cn0_alloy_BCD_y = y * mat2["Cn0"] + (1 - y) * mat4["Cn0"]
                    Cn0[startindex:finishindex] = (
                        (
                            x
                            * (1 - x)
                            * (y * Cn0_alloy_ABC_x + (1 - y) * Cn0_alloy_ABD_x)
                            + y
                            * (1 - y)
                            * (x * Cn0_alloy_ACD_y + (1 - x) * Cn0_alloy_BCD_y)
                        )
                        / (x * (1 - x) + y * (1 - y))
                        * 1e-12
                    )

                    Cp0_alloy_ABC_x = x * mat1["Cp0"] + (1 - x) * mat2["Cp0"]
                    Cp0_alloy_ABD_x = x * mat3["Cp0"] + (1 - x) * mat4["Cp0"]
                    Cp0_alloy_ACD_y = y * mat1["Cp0"] + (1 - y) * mat3["Cp0"]
                    Cp0_alloy_BCD_y = y * mat2["Cp0"] + (1 - y) * mat4["Cp0"]
                    Cp0[startindex:finishindex] = (
                        (
                            x
                            * (1 - x)
                            * (y * Cp0_alloy_ABC_x + (1 - y) * Cp0_alloy_ABD_x)
                            + y
                            * (1 - y)
                            * (x * Cp0_alloy_ACD_y + (1 - x) * Cp0_alloy_BCD_y)
                        )
                        / (x * (1 - x) + y * (1 - y))
                        * 1e-12
                    )

                    cb_meff_alpha[startindex:finishindex] = alloyprops["m_e_alpha"] * (
                        mat2["m_e"] / cb_meff_alloy
                    )  # non-parabolicity constant for alloy. THIS CALCULATION IS MOSTLY WRONG. MUST BE CONTROLLED. SBL
                if mat_crys_strc == "Wurtzite":
                    alloyprops = alloy_property_4[matType]
                    mat1 = material_property[alloyprops["Material1"]]  # GaN
                    mat2 = material_property[alloyprops["Material2"]]  # InN
                    mat3 = material_property[alloyprops["Material3"]]  # AlN
                    # This is accourding to interpolated Vegard’s law for quaternary BxCyD1-x-yA=AlxInyGa1-x-yN
                    # I. Vurgaftman, J.R. Meyer, L.R. RamMohan, J. Appl. Phys. 89 (2001) 5815.
                    # C. K. Williams, T. H. Glisson, J. R. Hauser, and M. A. Littlejohn, J. Electron. Mater. 7, 639 (1978).
                    x = layer[2]  # alloy ratio x
                    y = layer[3]  # alloy ratio y
                    u_4 = (1 - x + y) / 2
                    v_4 = (2 - x - 2 * y) / 2
                    w_4 = (2 - 2 * x - y) / 2
                    cb_meff_alloy_ABC = (
                        u_4 * mat2["m_e"] + (1 - u_4) * mat3["m_e"]
                    )  # AlInN
                    cb_meff_alloy_ACD = (
                        v_4 * mat1["m_e"] + (1 - v_4) * mat2["m_e"]
                    )  # InGaN
                    cb_meff_alloy_ABD = (
                        w_4 * mat1["m_e"] + (1 - w_4) * mat3["m_e"]
                    )  # AlGaN
                    cb_meff_alloy = (
                        x * y * cb_meff_alloy_ABC
                        + y * (1 - x - y) * cb_meff_alloy_ACD
                        + x * (1 - x - y) * cb_meff_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    cb_meff[startindex:finishindex] = cb_meff_alloy * m_e

                    Eg_alloy_ABC = (
                        u_4 * mat2["Eg"]
                        + (1 - u_4) * mat3["Eg"]
                        - alloyprops["Bowing_param_ABC"] * u_4 * (1 - u_4)
                    )  # eV AlInN
                    Eg_alloy_ACD = (
                        v_4 * mat1["Eg"]
                        + (1 - v_4) * mat2["Eg"]
                        - alloyprops["Bowing_param_ACD"] * v_4 * (1 - v_4)
                    )  # eV InGaN
                    Eg_alloy_ABD = (
                        w_4 * mat1["Eg"]
                        + (1 - w_4) * mat3["Eg"]
                        - alloyprops["Bowing_param_ABD"] * w_4 * (1 - w_4)
                    )  # eV AlGaN
                    Eg = (
                        x * y * Eg_alloy_ABC
                        + y * (1 - x - y) * Eg_alloy_ACD
                        + x * (1 - x - y) * Eg_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    fi_e[startindex:finishindex] = (
                        alloyprops["Band_offset"] * Eg * q
                    )  # for electron. Joule
                    a0_sub[startindex:finishindex] = alloyprops["a0_sub"] * 1e-10
                    A1_alloy_ABC = u_4 * mat2["A1"] + (1 - u_4) * mat3["A1"]
                    A1_alloy_ACD = v_4 * mat1["A1"] + (1 - v_4) * mat2["A1"]
                    A1_alloy_ABD = w_4 * mat1["A1"] + (1 - w_4) * mat3["A1"]
                    A1[startindex:finishindex] = (
                        x * y * A1_alloy_ABC
                        + y * (1 - x - y) * A1_alloy_ACD
                        + x * (1 - x - y) * A1_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    A2_alloy_ABC = u_4 * mat2["A2"] + (1 - u_4) * mat3["A2"]
                    A2_alloy_ACD = v_4 * mat1["A2"] + (1 - v_4) * mat2["A2"]
                    A2_alloy_ABD = w_4 * mat1["A2"] + (1 - w_4) * mat3["A2"]
                    A2[startindex:finishindex] = (
                        x * y * A2_alloy_ABC
                        + y * (1 - x - y) * A2_alloy_ACD
                        + x * (1 - x - y) * A2_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    A3_alloy_ABC = u_4 * mat2["A3"] + (1 - u_4) * mat3["A3"]
                    A3_alloy_ACD = v_4 * mat1["A3"] + (1 - v_4) * mat2["A3"]
                    A3_alloy_ABD = w_4 * mat1["A3"] + (1 - w_4) * mat3["A3"]
                    A3[startindex:finishindex] = (
                        x * y * A3_alloy_ABC
                        + y * (1 - x - y) * A3_alloy_ACD
                        + x * (1 - x - y) * A3_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    A4_alloy_ABC = u_4 * mat2["A4"] + (1 - u_4) * mat3["A4"]
                    A4_alloy_ACD = v_4 * mat1["A4"] + (1 - v_4) * mat2["A4"]
                    A4_alloy_ABD = w_4 * mat1["A4"] + (1 - w_4) * mat3["A4"]
                    A4[startindex:finishindex] = (
                        x * y * A4_alloy_ABC
                        + y * (1 - x - y) * A4_alloy_ACD
                        + x * (1 - x - y) * A4_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    A5_alloy_ABC = u_4 * mat2["A5"] + (1 - u_4) * mat3["A5"]
                    A5_alloy_ACD = v_4 * mat1["A5"] + (1 - v_4) * mat2["A5"]
                    A5_alloy_ABD = w_4 * mat1["A5"] + (1 - w_4) * mat3["A5"]
                    A5[startindex:finishindex] = (
                        x * y * A5_alloy_ABC
                        + y * (1 - x - y) * A5_alloy_ACD
                        + x * (1 - x - y) * A5_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    A6_alloy_ABC = u_4 * mat2["A6"] + (1 - u_4) * mat3["A6"]
                    A6_alloy_ACD = v_4 * mat1["A6"] + (1 - v_4) * mat2["A6"]
                    A6_alloy_ABD = w_4 * mat1["A6"] + (1 - w_4) * mat3["A6"]
                    A6[startindex:finishindex] = (
                        x * y * A6_alloy_ABC
                        + y * (1 - x - y) * A6_alloy_ACD
                        + x * (1 - x - y) * A6_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    D1_alloy_ABC = u_4 * mat2["D1"] + (1 - u_4) * mat3["D1"]
                    D1_alloy_ACD = v_4 * mat1["D1"] + (1 - v_4) * mat2["D1"]
                    D1_alloy_ABD = w_4 * mat1["D1"] + (1 - w_4) * mat3["D1"]
                    D1_alloy = (
                        x * y * D1_alloy_ABC
                        + y * (1 - x - y) * D1_alloy_ACD
                        + x * (1 - x - y) * D1_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    D1[startindex:finishindex] = D1_alloy * q

                    D2_alloy_ABC = u_4 * mat2["D2"] + (1 - u_4) * mat3["D2"]
                    D2_alloy_ACD = v_4 * mat1["D2"] + (1 - v_4) * mat2["D2"]
                    D2_alloy_ABD = w_4 * mat1["D2"] + (1 - w_4) * mat3["D2"]
                    D2_alloy = (
                        x * y * D2_alloy_ABC
                        + y * (1 - x - y) * D2_alloy_ACD
                        + x * (1 - x - y) * D2_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    D2[startindex:finishindex] = D2_alloy * q

                    D3_alloy_ABC = u_4 * mat2["D3"] + (1 - u_4) * mat3["D3"]
                    D3_alloy_ACD = v_4 * mat1["D3"] + (1 - v_4) * mat2["D3"]
                    D3_alloy_ABD = w_4 * mat1["D3"] + (1 - w_4) * mat3["D3"]
                    D3_alloy = (
                        x * y * D3_alloy_ABC
                        + y * (1 - x - y) * D3_alloy_ACD
                        + x * (1 - x - y) * D3_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    D3[startindex:finishindex] = D3_alloy * q

                    D4_alloy_ABC = u_4 * mat2["D4"] + (1 - u_4) * mat3["D4"]
                    D4_alloy_ACD = v_4 * mat1["D4"] + (1 - v_4) * mat2["D4"]
                    D4_alloy_ABD = w_4 * mat1["D4"] + (1 - w_4) * mat3["D4"]
                    D4_alloy = (
                        x * y * D4_alloy_ABC
                        + y * (1 - x - y) * D4_alloy_ACD
                        + x * (1 - x - y) * D4_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    D4[startindex:finishindex] = D4_alloy * q

                    D31_alloy_ABC = u_4 * mat2["D31"] + (1 - u_4) * mat3["D31"]
                    D31_alloy_ACD = v_4 * mat1["D31"] + (1 - v_4) * mat2["D31"]
                    D31_alloy_ABD = w_4 * mat1["D31"] + (1 - w_4) * mat3["D31"]
                    D31[startindex:finishindex] = (
                        x * y * D31_alloy_ABC
                        + y * (1 - x - y) * D31_alloy_ACD
                        + x * (1 - x - y) * D31_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    D33_alloy_ABC = u_4 * mat2["D33"] + (1 - u_4) * mat3["D33"]
                    D33_alloy_ACD = v_4 * mat1["D33"] + (1 - v_4) * mat2["D33"]
                    D33_alloy_ABD = w_4 * mat1["D33"] + (1 - w_4) * mat3["D33"]
                    D33[startindex:finishindex] = (
                        x * y * D33_alloy_ABC
                        + y * (1 - x - y) * D33_alloy_ACD
                        + x * (1 - x - y) * D33_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    Psp_alloy_ABC = u_4 * mat2["Psp"] + (1 - u_4) * mat3["Psp"]
                    Psp_alloy_ACD = v_4 * mat1["Psp"] + (1 - v_4) * mat2["Psp"]
                    Psp_alloy_ABD = w_4 * mat1["Psp"] + (1 - w_4) * mat3["Psp"]
                    Psp[startindex:finishindex] = (
                        x * y * Psp_alloy_ABC
                        + y * (1 - x - y) * Psp_alloy_ACD
                        + x * (1 - x - y) * Psp_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    C11_alloy_ABC = u_4 * mat2["C11"] + (1 - u_4) * mat3["C11"]
                    C11_alloy_ACD = v_4 * mat1["C11"] + (1 - v_4) * mat2["C11"]
                    C11_alloy_ABD = w_4 * mat1["C11"] + (1 - w_4) * mat3["C11"]
                    C11[startindex:finishindex] = (
                        (
                            x * y * C11_alloy_ABC
                            + y * (1 - x - y) * C11_alloy_ACD
                            + x * (1 - x - y) * C11_alloy_ABD
                        )
                        / (x * y + y * (1 - x - y) + x * (1 - x - y))
                        * 1e10
                    )  # for newton/meter²

                    C12_alloy_ABC = u_4 * mat2["C12"] + (1 - u_4) * mat3["C12"]
                    C12_alloy_ACD = v_4 * mat1["C12"] + (1 - v_4) * mat2["C12"]
                    C12_alloy_ABD = w_4 * mat1["C12"] + (1 - w_4) * mat3["C12"]
                    C12[startindex:finishindex] = (
                        (
                            x * y * C12_alloy_ABC
                            + y * (1 - x - y) * C12_alloy_ACD
                            + x * (1 - x - y) * C12_alloy_ABD
                        )
                        / (x * y + y * (1 - x - y) + x * (1 - x - y))
                        * 1e10
                    )

                    C13_alloy_ABC = u_4 * mat2["C13"] + (1 - u_4) * mat3["C13"]
                    C13_alloy_ACD = v_4 * mat1["C13"] + (1 - v_4) * mat2["C13"]
                    C13_alloy_ABD = w_4 * mat1["C13"] + (1 - w_4) * mat3["C13"]
                    C13[startindex:finishindex] = (
                        (
                            x * y * C13_alloy_ABC
                            + y * (1 - x - y) * C13_alloy_ACD
                            + x * (1 - x - y) * C13_alloy_ABD
                        )
                        / (x * y + y * (1 - x - y) + x * (1 - x - y))
                        * 1e10
                    )

                    C33_alloy_ABC = u_4 * mat2["C33"] + (1 - u_4) * mat3["C33"]
                    C33_alloy_ACD = v_4 * mat1["C33"] + (1 - v_4) * mat2["C33"]
                    C33_alloy_ABD = w_4 * mat1["C33"] + (1 - w_4) * mat3["C33"]
                    C33[startindex:finishindex] = (
                        (
                            x * y * C33_alloy_ABC
                            + y * (1 - x - y) * C33_alloy_ACD
                            + x * (1 - x - y) * C33_alloy_ABD
                        )
                        / (x * y + y * (1 - x - y) + x * (1 - x - y))
                        * 1e10
                    )

                    fi_h[startindex:finishindex] = (
                        -(1 - alloyprops["Band_offset"]) * Eg * q
                    )  # -(-1.33*(1-x)-0.8*x)for electron. Joule-1.97793434e-20 #

                    eps_alloy_ABC = (
                        u_4 * mat2["epsilonStatic"] + (1 - u_4) * mat3["epsilonStatic"]
                    )
                    eps_alloy_ACD = (
                        v_4 * mat1["epsilonStatic"] + (1 - v_4) * mat2["epsilonStatic"]
                    )
                    eps_alloy_ABD = (
                        w_4 * mat1["epsilonStatic"] + (1 - w_4) * mat3["epsilonStatic"]
                    )
                    eps_alloy = (
                        x * y * eps_alloy_ABC
                        + y * (1 - x - y) * eps_alloy_ACD
                        + x * (1 - x - y) * eps_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    eps[startindex:finishindex] = eps_alloy * eps0

                    a0_wz_alloy_ABC = u_4 * mat2["a0_wz"] + (1 - u_4) * mat3["a0_wz"]
                    a0_wz_alloy_ACD = v_4 * mat1["a0_wz"] + (1 - v_4) * mat2["a0_wz"]
                    a0_wz_alloy_ABD = w_4 * mat1["a0_wz"] + (1 - w_4) * mat3["a0_wz"]
                    a0_wz_alloy = (
                        x * y * a0_wz_alloy_ABC
                        + y * (1 - x - y) * a0_wz_alloy_ACD
                        + x * (1 - x - y) * a0_wz_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    a0_wz[startindex:finishindex] = a0_wz_alloy * 1e-10

                    delta_so_alloy_ABC = (
                        u_4 * mat2["delta_so"] + (1 - u_4) * mat3["delta_so"]
                    )
                    delta_so_alloy_ACD = (
                        v_4 * mat1["delta_so"] + (1 - v_4) * mat2["delta_so"]
                    )
                    delta_so_alloy_ABD = (
                        w_4 * mat1["delta_so"] + (1 - w_4) * mat3["delta_so"]
                    )
                    delta_so_alloy = (
                        x * y * delta_so_alloy_ABC
                        + y * (1 - x - y) * delta_so_alloy_ACD
                        + x * (1 - x - y) * delta_so_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    delta_so[startindex:finishindex] = delta_so_alloy * q

                    delta_cr_alloy_ABC = (
                        u_4 * mat2["delta_cr"] + (1 - u_4) * mat3["delta_cr"]
                    )
                    delta_cr_alloy_ACD = (
                        v_4 * mat1["delta_cr"] + (1 - v_4) * mat2["delta_cr"]
                    )
                    delta_cr_alloy_ABD = (
                        w_4 * mat1["delta_cr"] + (1 - w_4) * mat3["delta_cr"]
                    )
                    delta_cr_alloy = (
                        x * y * delta_cr_alloy_ABC
                        + y * (1 - x - y) * delta_cr_alloy_ACD
                        + x * (1 - x - y) * delta_cr_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    delta_cr[startindex:finishindex] = delta_cr_alloy * q

                    Ac_alloy_ABC = u_4 * mat2["Ac"] + (1 - u_4) * mat3["Ac"]
                    Ac_alloy_ACD = v_4 * mat1["Ac"] + (1 - v_4) * mat2["Ac"]
                    Ac_alloy_ABD = w_4 * mat1["Ac"] + (1 - w_4) * mat3["Ac"]
                    Ac_alloy = (
                        x * y * Ac_alloy_ABC
                        + y * (1 - x - y) * Ac_alloy_ACD
                        + x * (1 - x - y) * Ac_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))
                    Ac[startindex:finishindex] = Ac_alloy * q

                    mun0_alloy_ABC = u_4 * mat2["mun0"] + (1 - u_4) * mat3["mun0"]
                    mun0_alloy_ACD = v_4 * mat1["mun0"] + (1 - v_4) * mat2["mun0"]
                    mun0_alloy_ABD = w_4 * mat1["mun0"] + (1 - w_4) * mat3["mun0"]
                    mun0[startindex:finishindex] = (
                        x * y * mun0_alloy_ABC
                        + y * (1 - x - y) * mun0_alloy_ACD
                        + x * (1 - x - y) * mun0_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    mup0_alloy_ABC = u_4 * mat2["mup0"] + (1 - u_4) * mat3["mup0"]
                    mup0_alloy_ACD = v_4 * mat1["mup0"] + (1 - v_4) * mat2["mup0"]
                    mup0_alloy_ABD = w_4 * mat1["mup0"] + (1 - w_4) * mat3["mup0"]
                    mup0[startindex:finishindex] = (
                        x * y * mup0_alloy_ABC
                        + y * (1 - x - y) * mup0_alloy_ACD
                        + x * (1 - x - y) * mup0_alloy_ABD
                    ) / (x * y + y * (1 - x - y) + x * (1 - x - y))

                    Cn0_alloy_ABC = u_4 * mat2["Cn0"] + (1 - u_4) * mat3["Cn0"]
                    Cn0_alloy_ACD = v_4 * mat1["Cn0"] + (1 - v_4) * mat2["Cn0"]
                    Cn0_alloy_ABD = w_4 * mat1["Cn0"] + (1 - w_4) * mat3["Cn0"]
                    Cn0[startindex:finishindex] = (
                        (
                            x * y * Cn0_alloy_ABC
                            + y * (1 - x - y) * Cn0_alloy_ACD
                            + x * (1 - x - y) * Cn0_alloy_ABD
                        )
                        / (x * y + y * (1 - x - y) + x * (1 - x - y))
                        * 1e-12
                    )

                    Cp0_alloy_ABC = u_4 * mat2["Cp0"] + (1 - u_4) * mat3["Cp0"]
                    Cp0_alloy_ACD = v_4 * mat1["Cp0"] + (1 - v_4) * mat2["Cp0"]
                    Cp0_alloy_ABD = w_4 * mat1["Cp0"] + (1 - w_4) * mat3["Cp0"]
                    Cp0[startindex:finishindex] = (
                        (
                            x * y * Cp0_alloy_ABC
                            + y * (1 - x - y) * Cp0_alloy_ACD
                            + x * (1 - x - y) * Cp0_alloy_ABD
                        )
                        / (x * y + y * (1 - x - y) + x * (1 - x - y))
                        * 1e-12
                    )
            # wells and barriers boundaries
            if len(layer) == 5:
                # [thickness, material, alloy, doping, type]
                matDope = layer[3]
                matType = layer[4]
                matRole = "b"
            elif len(layer) == 6:
                # [thickness, material, alloy1, alloy2, doping, type]
                matDope = layer[4]
                matType = layer[5]
                matRole = "b"
            else:
                matDope = layer[4] if len(layer) > 4 else 0.0
                matType = layer[5] if len(layer) > 5 else "n"
                matRole = layer[6] if len(layer) > 6 else "b"

            if matRole == "w":
                N_wells_real2 += 1
                Well_boundary2[N_wells_real2, 0] = startindex
                Well_boundary2[N_wells_real2, 1] = finishindex
            N_layers_real2 += 1
            layer_boundary[N_layers_real2, 0] = startindex
            layer_boundary[N_layers_real2, 1] = finishindex
            for J in range(0, N_wells_virtual2):
                barrier_boundary[J, 0] = Well_boundary2[J - 1, 1]
                barrier_boundary[J, 1] = Well_boundary2[J, 0]
                barrier_len[J] = barrier_boundary[J, 1] - barrier_boundary[J, 0]
            # doping
            if len(self.dop_profile) != self.n_max:
                self.dop_profile = np.zeros(self.n_max)
            dop_profile = self.dop_profile
            if matType == "n":
                dop[startindex:finishindex] = (
                    matDope * 1e6 + dop_profile[startindex:finishindex] + 1
                )  # charge density in m**-3 (conversion from cm**-3)
            elif matType == "p":
                dop[startindex:finishindex] = (
                    -matDope * 1e6 + dop_profile[startindex:finishindex] - 1
                )  # charge density in m**-3 (conversion from cm**-3)
            else:
                dop[startindex:finishindex] = dop_profile[startindex:finishindex] + 1
        
        # Here we remove barriers that are less than the anti_crossing_length
        # so we can constructe the new well boundary using the resulted barrier boundary

        brr = 0
        anti_crossing_length = config.anti_crossing_length * 1e-9
        if not (self.Quantum_Regions):
            for J in range(2, N_wells_virtual2 - 1):
                if barrier_len[J] * dx <= anti_crossing_length:
                    brr += 1
            brr_vec = np.zeros(brr, dtype=int)
            brr2 = 0
            for J in range(2, N_wells_virtual2 - 1):
                if barrier_len[J] * dx <= anti_crossing_length:
                    brr2 += 1
                    brr_vec[brr2 - 1] = int(J + 1 - brr2)
            for I in range(0, brr):
                barrier_boundary = np.delete(barrier_boundary, int(brr_vec[I]), 0)
            N_wells_virtual = N_wells_virtual - brr
            Well_boundary = np.resize(Well_boundary, (N_wells_virtual, 2))
            for J in range(0, N_wells_virtual):
                Well_boundary[J - 1, 1] = barrier_boundary[J, 0]
                Well_boundary[J, 0] = barrier_boundary[J, 1]
        else:
            # setup of independent quantum regions
            # ratio of half well's width for wavefunction  to penetration into the the left adjacent barrier
            config.amort_wave_0 = 0.0
            config.amort_wave_1 = 0.0
            N_wells_real0 = len(self.Quantum_Regions_boundary[:, 0])
            N_wells_virtual = N_wells_real0 + 2
            N_wells_virtual2 = N_wells_real0 + 2
            N_layers_virtual = N_layers_real0 + 2
            Well_boundary = np.zeros((N_wells_virtual, 2), dtype=int)
            Well_boundary2 = np.zeros((N_wells_virtual, 2), dtype=int)
            barrier_boundary = np.zeros((N_wells_virtual + 1, 2), dtype=int)
            layer_boundary = np.zeros((N_layers_virtual, 2), dtype=int)
            n_max_general = np.zeros(N_wells_virtual, dtype=int)
            Well_boundary[N_wells_virtual - 1, 0] = n_max - 1
            Well_boundary[N_wells_virtual - 1, 1] = n_max - 1
            Well_boundary2[N_wells_virtual - 1, 0] = n_max - 1
            Well_boundary2[N_wells_virtual - 1, 1] = n_max - 1
            barrier_boundary[N_wells_virtual, 0] = n_max - 1
            barrier_len = np.zeros(N_wells_virtual + 1)
            for i in range(len(self.Quantum_Regions_boundary[:, 0])):
                for j in range(2):
                    Well_boundary[i + 1, j] = round2int(
                        self.Quantum_Regions_boundary[i, j] * 1e-9 / dx
                    )
        self.fi_e = fi_e
        self.fi_h = fi_h
        self.cb_meff = cb_meff
        self.cb_meff_alpha = cb_meff_alpha
        self.dop = dop
        self.pol_surf_char = pol_surf_char
        # return fi_e,cb_meff,eps,dop
        self.C11 = C11
        self.C12 = C12
        self.GA1 = GA1
        self.GA2 = GA2
        self.GA3 = GA3
        self.Ac = Ac
        self.Av = Av
        self.B = B
        self.n = n
        self.p = p
        self.a0 = a0
        self.delta = delta
        self.eps = eps
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.A4 = A4
        self.A5 = A5
        self.A6 = A6
        self.D1 = D1
        self.D2 = D2
        self.D3 = D3
        self.D4 = D4
        self.C13 = C13
        self.C33 = C33
        self.D31 = D31
        self.D33 = D33
        self.Psp = Psp
        self.a0_wz = a0_wz
        self.a0_sub = a0_sub
        self.delta_so = delta_so
        self.delta_cr = delta_cr
        self.N_wells_virtual = N_wells_virtual
        self.N_wells_virtual2 = N_wells_virtual2
        self.N_wells_real0 = N_wells_real0
        self.Well_boundary = Well_boundary
        self.Well_boundary2 = Well_boundary2
        self.barrier_boundary = barrier_boundary
        self.N_layers_real2 = N_layers_real2
        self.layer_boundary = layer_boundary
        self.TAUN0 = TAUN0
        self.TAUP0 = TAUP0
        self.mun0 = mun0
        self.mup0 = mup0
        self.Cn0 = Cn0
        self.Cp0 = Cp0
        self.BETAN = BETAN
        self.BETAP = BETAP
        self.VSATN = VSATN
        self.VSATP = VSATP


class AttrDict(dict):
    """turns a dictionary into an object with attribute style lookups"""

    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class StructureFrom(Structure):
    def __init__(self, inputfile, database):
        if type(inputfile) == dict:
            inputfile = AttrDict(inputfile)
        # Parameters for simulation
        defaults = {
            'Fapplied': 0.0,
            'vmax': 0.0,
            'vmin': 0.0,
            'Each_Step': 0.1,
            'surface': [0.0, 0.0],
            'T': 300.0,
            'subnumber_h': 1,
            'subnumber_e': 1,
            'computation_scheme': 0,
            'gridfactor': 0.1,
            'maxgridpoints': 200000,
            'max_iterations': 120,
            'dd_max_iterations': 25,
            'dd_residual_tolerance': 0.02,
            'dd_current_atol': 1e-8,  # mA/cm^2, total-current spatial span
            'dd_current_rtol': 1e-3,
            'mat_type': 'Zincblende',
            'dop_profile': np.zeros(1),
            'Quantum_Regions_boundary': np.zeros((1, 2)),
            'Quantum_Regions': False,
            'device_area': 1.0e-4, # cm^2 (default)
            'tat_field': 1.0e10, # V/m (default - disabled)
            'enable_polarization': True, # Default enabled
            'G_optical': 0.0, # cm^-3 s^-1 (default - dark)
            'photovoltaic_mode': False,
            'work_function_left': 4.5,
            'work_function_right': 5.2,
            'surface_recomb_val': [1e7, 1e7],
            'tau': None,  # Global lifetime override (s)
            'use_newton_solver': False, # Toggle fully-coupled Newton solver
            'enable_qw_solver': False, # Toggle QW Confined-State Solver
            'num_electron_states': 3, # Conduction subbands count
            'num_hole_states': 3, # Valence subbands count
            'qw_self_consistent': False, # Toggle self-consistent QW-Poisson
            'qw_max_iterations': 20,
            'qw_tolerance': 1e-4,
            'qw_damping': 0.2,
            'qw_coupling_mode': 'Coupled MQW'
        }
        for key, default in defaults.items():
            val = getattr(inputfile, key, default)
            setattr(self, key if key != 'Fapplied' else 'Fapp', val)
            
        # Additional PV settings support
        if not hasattr(self, 'surface_recomb'):
             self.surface_recomb = getattr(inputfile, 'surface_recomb', self.surface_recomb_val)
        
        self.Vt = Vt # Thermal voltage
        # Mapping compatibility
        # If comp_scheme is provided in input (and not None), use it; otherwise use computation_scheme
        cs = getattr(inputfile, 'comp_scheme', None)
        if cs is None:
            cs = getattr(inputfile, 'computation_scheme', getattr(self, 'computation_scheme', 0))
        self.comp_scheme = int(cs)
        self.computation_scheme = self.comp_scheme
        print(f"DEBUG: Structure initialized. photovoltaic_mode={self.photovoltaic_mode}, comp_scheme={self.comp_scheme}")
        self.dx = getattr(inputfile, 'gridfactor', 0.1) * 1e-9  # grid in m
        self.mat_crys_strc = self.mat_type
        # Area in m^2 (input is in cm^2)
        self.device_area_m2 = getattr(inputfile, 'device_area', 1.0e-4) * 1e-4
        
        # Loading material list
        self.material = inputfile.material
        self.inputfilename = inputfile
        totallayer = alen(self.material)

        # Add to log
        logger.info("Total layer number: %s", totallayer)

        # Calculate the required number of grid points
        self.x_max = (
            sum([layer[0] for layer in self.material]) * 1e-9
        )  # total thickness (m)
        self.n_max = int(self.x_max / self.dx)
        # Check on n_max
        max_val = self.maxgridpoints

        if self.n_max > max_val:
            logger.error("Grid number is exceeding the max number of %d", max_val)
            sys.exit()
        # Loading materials database #
        self.material_property = database.materialproperty
        totalmaterial = alen(self.material_property)

        self.alloy_property = database.alloyproperty
        totalalloy = alen(self.alloy_property)

        self.alloy_property_4 = database.alloyproperty4
        totalalloy += alen(self.alloy_property_4)

        # Extract Series Resistance (Rs) from first layer material
        try:
            first_mat = self.material[0][1]
            if first_mat in self.material_property:
                self.Rs = self.material_property[first_mat].get('Rs', 0.0)
                logger.info(f"Model internal Rs set to {self.Rs} Ohm (from {first_mat})")
            else:
                self.Rs = 0.0
        except Exception as e:
            logger.warning(f"Could not extract Rs: {e}")
            self.Rs = 0.0
        # Add to log
        logger.info("Total number of materials in database: %d" % (totalmaterial + totalalloy))
        # Initialise arrays

        # cb_meff #conduction band effective mass (array, len n_max)
        # fi_e #Bandstructure potential (array, len n_max)
        # eps #dielectric constant (array, len n_max)
        # dop #doping distribution (array, len n_max)
        self.create_structure_arrays()
        
        # Override doping if provided in input (Validation fix)
        # Override doping if provided in input (Validation fix)
        # Check self.dop_profile (already loaded from input)
        if hasattr(self, 'dop_profile'):
             logger.info(f"DEBUG: Found self.dop_profile (unconditional) with len {len(self.dop_profile)}. n_max: {self.n_max}")
        else:
             logger.info("DEBUG: self.dop_profile NOT FOUND")
             
        if hasattr(inputfile, 'dop_profile') and len(inputfile.dop_profile) > 1 and np.any(inputfile.dop_profile != 0):
            if len(inputfile.dop_profile) == self.n_max:
                self.dop = inputfile.dop_profile
                logger.info("Overriding doping profile from input configuration.")
            else:
                 logger.warning(f"Input dop_profile length {len(inputfile.dop_profile)} does not match n_max {self.n_max}. Ignoring.")
        else:
             logger.info("Using doping profile calculated from material layers.")
        
        self.ionization_efficiency = getattr(inputfile, 'ionization_efficiency', 1.0)
        self.poisson_damping = getattr(inputfile, 'poisson_damping', 0.1)
        self.continuity_damping = getattr(inputfile, 'continuity_damping', 0.7)
        
        # Apply ionization efficiency to p-type dopants (deep acceptors)
        if self.ionization_efficiency != 1.0:
            for i in range(len(self.dop)):
                if self.dop[i] < 0:
                    self.dop[i] *= self.ionization_efficiency


# No Shooting method parameters for Schrödinger Equation solution since we use a 3x3 KP solver
# delta_E = 1.0*meV2J #Energy step (Joules) for initial search. Initial delta_E is 1 meV. #This can be included in config as a setting?
# d_E = 1e-5*meV2J #Energy step (Joules) for Newton-Raphson method when improving the precision of the energy of a found level.
"""damping:An adjustable parameter  (0 < damping < 1) is typically set to 0.5 at low carrier densities. With increasing
carrier densities, a smaller value of it is needed for rapid convergence."""
damping = 0.1  # averaging factor between iterations to smooth convergence.
max_iterations = 120  # maximum number of iterations.
convergence_test = 1e-5  # convergence is reached when the ground state energy (eV) is stable to within this number between iterations.
convergence_test0 = 1e-5
# DO NOT EDIT UNDER HERE FOR PARAMETERS
# --------------------------------------

# Vegard's law for alloys
def vegard(first, second, mole):
    return first * mole + second * (1 - mole)


# FUNCTIONS for FERMI-DIRAC STATISTICS-----------------------------------------
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


def calc_meff_state_general(
    wfh,
    wfe,
    model,
    fi_e,
    E_statec,
    list,
    m_hh,
    m_lh,
    m_so,
    n_max_general,
    j,
    Well_boundary,
    n_max,
):
    vb_meff = np.zeros((model.subnumber_h, n_max_general))
    #
    I1, I2, I11, I22 = amort_wave(j, Well_boundary, n_max)
    i2 = I2 - I1
    for i in range(0, model.subnumber_h, 1):
        if list[i] == "hh1" or list[i] == "hh2" or list[i] == "hh3":
            vb_meff[i] = m_hh[I1:I2]
        elif list[i] == "lh1" or list[i] == "lh2" or list[i] == "lh3":
            vb_meff[i] = m_lh[I1:I2]
        else:
            vb_meff[i] = m_so[I1:I2]
    tmp = 1.0 / np.sum(wfh[:, 0:i2] ** 2 / vb_meff, axis=1)  # vb_meff[:,int(n_max/2)]
    meff_state = tmp.tolist()
    """find subband effective masses including non-parabolicity
    (but stilling using a fixed effective mass for each subband dispersion)"""
    cb_meff = model.cb_meff  # effective mass of conduction band across structure
    cb_meff_alpha = model.cb_meff_alpha  # non-parabolicity constant across structure
    cb_meff_states = np.array(
        [cb_meff * (1.0 + cb_meff_alpha * (E * meV2J - fi_e)) for E in E_statec]
    )
    tmp1 = 1.0 / np.sum(wfe[:, 0:i2] ** 2 / cb_meff_states[:, I1:I2], axis=1)
    meff_statec = tmp1.tolist()
    return meff_statec, meff_state


def calc_meff_state(wfh, wfe, subnumber_h, subnumber_e, list, m_hh, m_lh, m_so, model):
    n_max = len(m_hh)
    vb_meff = np.zeros((subnumber_h, n_max))
    for i in range(0, subnumber_h, 1):
        if list[i] == "hh":
            vb_meff[i] = m_hh
        elif list[i] == "lh":
            vb_meff[i] = m_lh
        else:
            vb_meff[i] = m_so
    tmp = 1.0 / np.sum(wfh ** 2 / vb_meff, axis=1)
    meff_state = tmp.tolist()
    """find subband effective masses including non-parabolicity
    (but stilling using a fixed effective mass for each subband dispersion)"""
    cb_meff = model.cb_meff  # effective mass of conduction band across structure
    # cb_meff_alpha = model.cb_meff_alpha # non-parabolicity constant across structure
    # cb_meff_states = np.array([cb_meff*(1.0 + cb_meff_alpha*(E*meV2J - fi_e)) for E in E_statec])
    # tmp1 = 1.0/np.sum(wfe**2/cb_meff_states,axis=1)
    tmp1 = 1.0 / np.sum(wfe ** 2 / cb_meff, axis=1)
    meff_statec = tmp1.tolist()
    return meff_statec, meff_state


def fermilevel_0Kc(Ntotal2d, E_statec, meff_statec, model):  # use
    Et2, Ef = 0.0, 0.0
    meff_statec = np.array(meff_statec)
    E_statec = np.array(E_statec)
    for i in range(
        model.subnumber_e, 0, -1
    ):  # ,(Ei,vsb_meff) in enumerate(zip(E_state,meff_state)):
        Efnew2 = sum(E_statec[0:i] * meff_statec[0:i])
        m2 = sum(meff_statec[0:i])
        Et2 += E_statec[i - model.subnumber_e]
        Efnew = (Efnew2 + Ntotal2d * hbar ** 2 * pi * J2meV) / (m2)
        if Efnew > Et2:
            Ef = Efnew
            # print 'Ef[',i-subnumber_h,']=',Ef
        else:
            break  # we have found Ef and so we should break out of the loop
    else:  # exception clause for 'for' loop.
        # Add to log
        logger.warning("Have processed all energy levels present and so can't be sure that Ef is below next higher energy level.")
    # Ef1=(sum(E_state*meff_state)-Ntotal2d*hbar**2*pi)/(sum(meff_state))
    N_statec = [0.0] * len(E_statec)
    for i, (Ei, csb_meff) in enumerate(zip(E_statec, meff_statec)):
        Nic = (Ef - Ei) * csb_meff / (hbar ** 2 * pi) * meV2J  # populations of levels
        Nic *= Nic > 0.0
        N_statec[i] = Nic
    return (
        Ef,
        N_statec,
    )  # Fermi levels at 0K (meV), number of electrons in each subband at 0K


def fermilevel_0K(Ntotal2d, E_state, meff_state, model):  # use
    Et1, Ef = 0.0, 0.0
    E_state = np.array(E_state)
    for i in range(
        model.subnumber_h, 0, -1
    ):  # ,(Ei,vsb_meff) in enumerate(zip(E_state,meff_state)):
        Efnew1 = sum(E_state[0:i] * meff_state[0:i])
        m1 = sum(meff_state[0:i])
        Et1 += E_state[i - model.subnumber_h]
        Efnew = (Efnew1 + Ntotal2d * hbar ** 2 * pi * J2meV) / (m1)
        if Efnew < Et1:
            Ef = Efnew
            # print 'Ef[',i-subnumber_h,']=',Ef
        else:
            break  # we have found Ef and so we should break out of the loop
    else:  # exception clause for 'for' loop.
        # Add to log
        logger.warning("Have processed all energy levels present and so can't be sure that Ef is below next higher energy level.")
    # Ef1=(sum(E_state*meff_state)-Ntotal2d*hbar**2*pi)/(sum(meff_state))
    N_state = [0.0] * len(E_state)
    for i, (Ei, vsb_meff) in enumerate(zip(E_state, meff_state)):
        Ni = (Ei - Ef) * vsb_meff / (hbar ** 2 * pi) * meV2J  # populations of levels
        Ni *= Ni > 0.0
        N_state[i] = Ni
    return (
        Ef,
        N_state,
    )  # Fermi levels at 0K (meV), number of electrons in each subband at 0K


def fermilevel(Ntotal2d, model, E_state, E_statec, meff_state, meff_statec):  # use
    # find the Fermi level (meV)
    def func(Ef, E_state, meff_state, E_statec, meff_statec, Ntotal2d, model):
        # return Ntotal2d - sum( [vsb_meff*fd2(Ei,Ef,T) for Ei,vsb_meff in zip(E_state,meff_state)] )/(hbar**2*pi)
        diff, diff1, diff2 = 0.0, 0.0, 0.0
        diff = Ntotal2d
        for Ei, csb_meff in zip(E_statec, meff_statec):
            diff1 -= csb_meff * fd2(Ei, Ef, model) / (hbar ** 2 * pi)
        for Ei, vsb_meff in zip(E_state, meff_state):
            diff2 += vsb_meff * fd1(Ei, Ef, model) / (hbar ** 2 * pi)
        if Ntotal2d > 0:
            diff += diff1
        else:
            diff += diff2
        return diff

    if Ntotal2d > 0:
        Ef_0K, N_states_0K = fermilevel_0Kc(Ntotal2d, E_statec, meff_statec, model)
    else:
        Ef_0K, N_states_0K = fermilevel_0K(Ntotal2d, E_state, meff_state, model)
    # Ef=fsolve(func,Ef_0K,args=(E_state,meff_state,Ntotal2d,T))[0]
    # return float(Ef)
    # implement Newton-Raphson method
    Ef = Ef_0K
    # itr=0
    # logger.info('Ef (at 0K)= %g',Ef)
    d_E = 1e-9  # Energy step (meV)
    while True:
        y = func(Ef, E_state, meff_state, E_statec, meff_statec, Ntotal2d, model)
        dy = (
            func(Ef + d_E, E_state, meff_state, E_statec, meff_statec, Ntotal2d, model)
            - func(
                Ef - d_E, E_state, meff_state, E_statec, meff_statec, Ntotal2d, model
            )
        ) / (2.0 * d_E)
        if (
            dy == 0.0
        ):  # increases interval size for derivative calculation in case of numerical error
            d_E *= 2.0
            continue
        Ef -= y / dy
        if abs(y / dy) < 1e-12:
            break
        for i in range(2):
            if d_E > 1e-9:
                d_E *= 0.5
    return Ef  # (meV)


def calc_N_state(
    Ef, model, E_state, meff_state, E_statec, meff_statec, Ntotal2d
):  # use
    # Find the subband populations, taking advantage of step like d.o.s. and analytic integral of FD
    N_statec, N_state = 0.0, 0.0
    if Ntotal2d > 0:
        N_statec = [
            fd2(Ei, Ef, model) * csb_meff / (hbar ** 2 * pi)
            for Ei, csb_meff in zip(E_statec, meff_statec)
        ]
    else:
        N_state = [
            fd1(Ei, Ef, model) * vsb_meff / (hbar ** 2 * pi)
            for Ei, vsb_meff in zip(E_state, meff_state)
        ]
    return N_state, N_statec  # number of carriers in each subband


# FUNCTIONS for SELF-CONSISTENT POISSON--------------------------------


def calc_sigma(wfh, wfe, N_state, N_statec, model, Ntotal2d):  # use
    """This function calculates `net' areal charge density
    n-type dopants lead to -ve charge representing electrons, and additionally 
    +ve ionised donors."""
    # note: model.dop is still a volume density, the delta_x converts it to an areal density
    sigma = model.dop * model.dx  # The charges due to the dopant ions
    if Ntotal2d > 0:
        for j in range(
            0, model.subnumber_e, 1
        ):  # The charges due to the electrons in the subbands
            sigma -= N_statec[j] * (wfe[j]) ** 2
    else:
        for i in range(
            0, model.subnumber_h, 1
        ):  # The charges due to the electrons in the subbands
            sigma += N_state[i] * (wfh[i]) ** 2
    return sigma  # charge per m**2 (units of electronic charge)


def calc_sigma_general2(n_max, dopi, n, p):  # use
    """This function calculates `net' areal charge density
    n-type dopants lead to -ve charge representing electrons, and additionally 
    +ve ionised donors."""
    sigma = np.zeros(len(dopi))
    sigma = sigma + dopi  # The charges due to the dopant ions
    for i in range(0, n_max):  # The charges due to the electrons in the subbands
        sigma[i] += p[i] - n[i]
    return sigma  # charge per m**3 (units of electronic charge)


def calc_sigma_general(
    pol_surf_char, wfh, wfe, N_state, N_statec, model, Ntotal2d, j, Well_boundary
):  # use
    """This function calculates `net' areal charge density
    n-type dopants lead to -ve charge representing electrons, and additionally 
    +ve ionised donors."""
    # note: model.dop is still a volume density, the delta_x converts it to an areal density
    sigma = (
        model.dop[Well_boundary[j - 1, 1] : Well_boundary[j + 1, 0]] * model.dx
    )  # +pol_surf_char[Well_boundary[j-1,1]:Well_boundary[j+1,0]]  The charges due to the dopant ions
    if Ntotal2d > 0:
        for j in range(
            0, model.subnumber_e, 1
        ):  # The charges due to the electrons in the subbands
            sigma -= N_statec[j] * (wfe[j]) ** 2
    else:
        for i in range(
            0, model.subnumber_h, 1
        ):  # The charges due to the electrons in the subbands
            sigma += N_state[i] * (wfh[i]) ** 2
    return sigma  # charge per m**2 (units of electronic charge)


def calc_field(sigma, eps):
    # F electric field as a function of z-
    # i index over z co-ordinates
    # j index over z' co-ordinates
    # Note: sigma is a number density per unit area, needs to be converted to Couloumb per unit area
    sigma = sigma
    F0 = -np.sum(q * sigma) / (2.0)  # CMP'deki i ve j yer değişebilir - de + olabilir
    # is the above necessary since the total field due to the structure should be zero.
    # Do running integral
    tmp = np.hstack(([0.0], sigma[:-1])) + sigma
    tmp *= (
        q / 2.0
    )  # Note: sigma is a number density per unit area, needs to be converted to Couloumb per unit area
    tmp[0] = F0
    F = np.cumsum(tmp) / eps
    return F


def calc_field_convolve(sigma, eps):  # use
    tmp = np.ones(len(sigma) - 1)
    signstep = np.hstack((-tmp, [0.0], tmp))  # step function
    F = np.convolve(signstep, sigma, mode="valid")
    F *= q / (2.0 * eps)
    return F


def calc_field_old(sigma, eps):  # use
    # F electric field as a function of z-
    # i index over z co-ordinates
    # j index over z' co-ordinates
    n_max = len(sigma)
    # For wave function initialise F
    F = np.zeros(n_max)
    for i in range(0, n_max, 1):
        for j in range(0, n_max, 1):
            # Note sigma is a number density per unit area, needs to be converted to Couloumb per unit area
            F[i] = F[i] + q * sigma[j] * cmp(i, j) / (
                2 * eps[i]
            )  # CMP'deki i ve j yer değişebilir - de + olabilir
    return F


def calc_potn(F, model):  # use
    # This function calculates the potential (energy actually)
    # V electric field as a function of z-
    # i	index over z co-ordinates

    # Calculate the potential, defining the first point as zero
    tmp = q * F * model.dx
    V = np.cumsum(tmp)  # +q -> electron -q->hole?
    return V


# FUNCTIONS FOR EXCHANGE INTERACTION-------------------------------------------


def calc_Vxc(sigma, eps, cb_meff, model):
    """An effective field describing the exchange-interactions between the electrons
    derived from Kohn-Sham density functional theory. This formula is given in many
    papers, for example see Gunnarsson and Lundquist (1976), Ando, Taniyama, Ohtani 
    et al. (2003), or Ch.1 in the book 'Intersubband transitions in quantum wells' (edited
    by Liu and Capasso) by M. Helm.
    eps = dielectric constant array
    cb_meff = effective mass array
    sigma = charge carriers per m**2, however this includes the donor atoms and we are only
            interested in the electron density."""
    a_B = 4 * pi * hbar ** 2 / q ** 2  # Bohr radius.
    nz = -(sigma - model.dop * model.dx)  # electron density per m**2
    nz_3 = nz ** (1 / 3.0)  # cube root of charge density.
    # a_B_eff = eps/cb_meff*a_B #effective Bohr radius
    # r_s occasionally suffers from division by zero errors due to nz=0.
    # We will fix these by setting nz_3 = 1.0 for these points (a tiny charge in per m**2).
    nz_3 = nz_3.clip(1.0, max(nz_3))

    r_s = 1.0 / (
        (4 * pi / 3.0) ** (1 / 3.0) * nz_3 * eps / cb_meff * a_B
    )  # average distance between charges in units of effective Bohr radis.
    # A = q**4/(32*pi**2*hbar**2)*(9*pi/4.0)**(1/3.)*2/pi*(4*pi/3.0)**(1/3.)*4*pi*hbar**2/q**2 #constant factor for expression.
    A = (
        q ** 2 / (4 * pi) * (3 / pi) ** (1 / 3.0)
    )  # simplified constant factor for expression.
    #
    Vxc = -A * nz_3 / eps * (1.0 + 0.0545 * r_s * np.log(1.0 + 11.4 / r_s))
    ionization_efficiency = 1.0  # Fraction of p-type dopants that are active
    use_newton_solver = False # Toggle fully-coupled Newton solver
    return Vxc


# -----------------------------------------------------------------------------


def wave_func_tri(j, Well_boundary, n_max, V1, V2, subnumber_h, subnumber_e, model):
    # Envelope Function Wave Functions
    wfh_general = np.zeros((model.N_wells_virtual, subnumber_h, n_max))
    wfe_general = np.zeros((model.N_wells_virtual, subnumber_e, n_max))
    # n_max_general = np.zeros(model.N_wells_virtual,int)
    n_max_general2 = np.zeros(model.N_wells_virtual, int)
    I1, I2, I11, I22 = amort_wave(j, Well_boundary, n_max)
    i_1 = I2 - I1
    n_max_general2[j] = int(I2 - I1)
    # n_max_general[j]=int(Well_boundary[j+1,0]-Well_boundary[j-1,1])
    wfh1s2 = np.zeros((subnumber_h, 3, i_1))
    maxwfh = np.zeros((subnumber_h, 3))
    list = [""] * subnumber_h
    for i in range(0, subnumber_e, 1):
        wfe_general[j, i, 0:i_1] = V1[j, 0:i_1, i] + 1e-20
    wfh_pow = np.zeros(n_max)
    conter_hh, conter_lh, conter_so = 0, 0, 0
    for jj in range(0, subnumber_h):
        for i in range(0, 3):
            wfh1s2[jj, i, :] = V2[j, i * i_1 : (i + 1) * i_1, jj]
            wfh_pow = np.cumsum(wfh1s2[jj, i, :] * wfh1s2[jj, i, :])
            maxwfh[jj, i] = wfh_pow[i_1 - 1]
        if np.argmax(maxwfh[jj, :]) == 0:
            conter_hh += 1
            list[jj] = "hh%d" % conter_hh
            wfh_general[j, jj, 0:i_1] = wfh1s2[jj, np.argmax(maxwfh[jj, :]), :] + 1e-20
        elif np.argmax(maxwfh[jj, :]) == 1:
            conter_lh += 1
            list[jj] = "lh%d" % conter_lh
            wfh_general[j, jj, 0:i_1] = wfh1s2[jj, np.argmax(maxwfh[jj, :]), :] + 1e-20
        else:
            conter_so += 1
            list[jj] = "so%d" % conter_so
            wfh_general[j, jj, 0:i_1] = wfh1s2[jj, np.argmax(maxwfh[jj, :]), :] + 1e-20
    return wfh_general, wfe_general, list, n_max_general2


def Strain_and_Masses(model):
    n_max = model.n_max
    EXX = np.zeros(n_max)
    EZZ = np.zeros(n_max)
    ZETA = np.zeros(n_max)
    CNIT = np.zeros(n_max)
    VNIT = np.zeros(n_max)
    S = np.zeros(n_max)
    k1 = np.zeros(n_max)
    k2 = np.zeros(n_max)
    k3 = np.zeros(n_max)
    fp = np.ones(n_max)
    fm = np.ones(n_max)
    EPC = np.zeros(n_max)
    m_hh = np.zeros(n_max)
    m_lh = np.zeros(n_max)
    m_so = np.zeros(n_max)
    Ppz = np.zeros(n_max)
    Ppz_Psp = np.zeros(n_max)
    Ppz_Psp0 = np.zeros(n_max)
    pol_surf_char = np.zeros(n_max)
    pol_surf_char1 = np.zeros(n_max)
    x_max = model.dx * n_max
    if config.strain:
        if model.mat_crys_strc == "Zincblende":
            EXX = (model.a0_sub - model.a0) / model.a0
            EZZ = -2.0 * model.C12 / model.C11 * EXX
            ZETA = -model.B / 2.0 * (EXX + EXX - 2.0 * EZZ)
            CNIT = model.Ac * (EXX + EXX + EZZ)
            VNIT = -model.Av * (EXX + EXX + EZZ)
        if model.mat_crys_strc == "Wurtzite":
            EXX = (model.a0_sub - model.a0_wz) / model.a0_wz
            # EXX= (4.189*1e-10-model.a0_wz)/model.a0_wz
            EZZ = -2.0 * model.C13 / model.C33 * EXX
            CNIT = model.Ac * (EXX + EXX + EZZ)
            ZETA = model.D2 * (EXX + EXX) + model.D1 * EZZ
            VNIT = model.D4 * (EXX + EXX) + model.D3 * EZZ
            Ppz = (model.D31 * (model.C11 + model.C12) + model.D33 * model.C13) * (
                EXX + EXX
            ) + (2 * model.D31 * model.C13 + model.D33 * model.C33) * (EZZ)
            dx = x_max / n_max
            sum_1 = 0.0
            sum_2 = 0.0
            if config.piezo:
                """ Spontaneous and piezoelectric polarization built-in field 
                [1] F. Bernardini and V. Fiorentini phys. stat. sol. (b) 216, 391 (1999)
                [2] Book 'Quantum Wells,Wires & Dots', Paul Harrison, pages 236-241"""
                for J in range(1, model.N_wells_virtual2 - 1):
                    BW = model.Well_boundary2[J, 0]
                    WB = model.Well_boundary2[J, 1]
                    Lw = (WB - BW) * dx
                    lb1 = (BW - model.Well_boundary2[J - 1, 1]) * dx
                    # lb2=(Well_boundary2[J+1,0]-WB)*dx
                    sum_1 += (model.Psp[BW + 1] + Ppz[BW + 1]) * Lw / model.eps[
                        BW + 1
                    ] + (model.Psp[BW - 1] + Ppz[BW - 1]) * lb1 / model.eps[BW - 1]
                    sum_2 += Lw / model.eps[BW + 1] + lb1 / model.eps[BW - 1]
                EPC = (sum_1 - (model.Psp + Ppz) * sum_2) / (model.eps * sum_2)
            if config.piezo1:
                pol_surf_char = np.zeros(n_max)
                pol_surf_char1 = np.zeros(n_max)
                for i in range(0, n_max):
                    pol_surf_char[i] = (model.Psp[i] + Ppz[i]) / (q)
                for i in range(1, n_max - 1):
                    pol_surf_char1[i] = (
                        (model.Psp[i - 1] + Ppz[i - 1])
                        - (model.Psp[i + 1] + Ppz[i + 1])
                    ) / (q)
                for i in range(1, n_max - 1):
                    Ppz_Psp0[i] = (pol_surf_char[i] - pol_surf_char[i - 1]) / (dx)
                for I in range(1, model.N_wells_virtual2 - 1):
                    BW = model.Well_boundary2[I, 0]
                    WB = model.Well_boundary2[I, 1]
                    Ppz_Psp[WB] = (pol_surf_char[WB + 1] - pol_surf_char[WB - 1]) / (dx)
                    Ppz_Psp[BW] = (pol_surf_char[BW + 1] - pol_surf_char[BW - 1]) / (dx)
                Ppz_Psp0[0] = (pol_surf_char[0] - 0.0) / (dx)

                Ppz_Psp0[n_max - 1] = (0.0 - pol_surf_char[n_max - 1]) / (dx)

    if config.piezo1 and not (config.strain):
        if model.mat_crys_strc == "Zincblende":
            EXX = (model.a0_sub - model.a0) / model.a0
            EZZ = -2.0 * model.C12 / model.C11 * EXX
        if model.mat_crys_strc == "Wurtzite":
            EXX = (model.a0_sub - model.a0_wz) / model.a0_wz
            EZZ = -2.0 * model.C13 / model.C33 * EXX
        dx = x_max / n_max
        Ppz = (model.D31 * (model.C11 + model.C12) + model.D33 * model.C13) * (
            EXX + EXX
        ) + (2 * model.D31 * model.C13 + model.D33 * model.C33) * (EZZ)
        if config.piezo1:
            pol_surf_char = np.zeros(n_max)
            pol_surf_char1 = np.zeros(n_max)
            for i in range(0, n_max):
                pol_surf_char[i] = (model.Psp[i] + Ppz[i]) / (q)
            for i in range(1, n_max - 1):
                Ppz_Psp0[i] = (pol_surf_char[i] - pol_surf_char[i - 1]) / (dx)
            for i in range(1, n_max - 1):
                pol_surf_char1[i] = (
                    (model.Psp[i - 1] + Ppz[i - 1]) - (model.Psp[i + 1] + Ppz[i + 1])
                ) / (q)
        for I in range(1, model.N_wells_virtual2 - 1):
            BW = model.Well_boundary2[I, 0]
            WB = model.Well_boundary2[I, 1]
            Ppz_Psp[WB] = (pol_surf_char[WB + 1] - pol_surf_char[WB - 1]) / (dx)
            Ppz_Psp[BW] = (pol_surf_char[BW + 1] - pol_surf_char[BW - 1]) / (dx)
    if model.mat_crys_strc == "Zincblende":
        for i in range(0, n_max, 1):
            if EXX[i] != 0:
                S[i] = ZETA[i] / model.delta[i]
                k1[i] = sqrt(1 + 2 * S[i] + 9 * S[i] ** 2)
                k2[i] = S[i] - 1 + k1[i]
                k3[i] = S[i] - 1 - k1[i]
                fp[i] = (2 * S[i] * (1 + 1.5 * k2[i]) + 6 * S[i] ** 2) / (
                    0.75 * k2[i] ** 2 + k2[i] - 3 * S[i] ** 2
                )
                fm[i] = (2 * S[i] * (1 + 1.5 * k3[i]) + 6 * S[i] ** 2) / (
                    0.75 * k3[i] ** 2 + k3[i] - 3 * S[i] ** 2
                )
        m_hh = m_e / (model.GA1 - 2 * model.GA2)
        m_lh = m_e / (model.GA1 + 2 * fp * model.GA2)
        m_so = m_e / (model.GA1 + 2 * fm * model.GA2)
    if model.mat_crys_strc == "Wurtzite":
        m_hh = -m_e / (model.A2 + model.A4 - model.A5)
        m_lh = -m_e / (model.A2 + model.A4 + model.A5)
        m_so = -m_e / (model.A2)
    return m_hh, m_lh, m_so, VNIT, ZETA, CNIT, Ppz_Psp0, EPC, pol_surf_char


def calc_E_state_general(
    HUPMAT3_reduced_list,
    HUPMATC1,
    subnumber_h,
    subnumber_e,
    fitot,
    fitotc,
    model,
    Well_boundary,
    UNIM,
    RATIO,
):
    n_max = model.n_max
    n_max_general = np.zeros(model.N_wells_virtual, dtype=int)
    # HUPMAT3=np.zeros((n_max*3, n_max*3))
    # HUPMAT3=VBMAT_V(HUPMAT1,fitot,RATIO,n_max,UNIM)
    HUPMATC3 = CBMAT_V(HUPMATC1, fitotc, RATIO, n_max, UNIM)
    # stop
    tmp1 = np.zeros((model.N_wells_virtual, n_max))
    KPV1 = np.zeros((model.N_wells_virtual, subnumber_e))
    V1 = np.zeros((model.N_wells_virtual, n_max, n_max))
    V11 = np.zeros((model.N_wells_virtual, n_max, n_max))
    for J in range(1, model.N_wells_virtual - 1):
        n_max_general[J] = Well_boundary[J + 1, 0] - Well_boundary[J - 1, 1]
        I1, I2, I11, I22 = amort_wave(J, Well_boundary, n_max)
        i_1 = I2 - I1
        i1 = I1 - I1
        i2 = I2 - I1
        la1, v1 = linalg.eigh(HUPMATC3[I1:I2, I1:I2])
        tmp1[J, i1:i2] = la1 / RATIO * J2meV
        V1[J, i1:i2, i1:i2] = v1
        if max(tmp1[J, 0:subnumber_e]) > max(fitotc[I11:I22]) * J2meV and 1 == 2:
            logger.warning(
                ":You may experience convergence problem due to unconfined states."
            )
    """ 
    for j in range(1,model.N_wells_virtual-1):            
        for i in range(0,subnumber_e,1):
            KPV1[j,i]=tmp1[j,i]
    """
    for j in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(j, Well_boundary, n_max)
        i_1 = I2 - I1
        i1 = I1 - I1
        i2 = I2 - I1
        i11 = I11 - I1
        i22 = I22 - I1
        couter = 0
        for i in range(i1, i2):
            wfe_pow1 = np.cumsum(
                V1[j, i11:i22, i] * V1[j, i11:i22, i]
            )  # and tmp1[j,i]<max(fitotc[I11-1:I22+1])*J2meV
            if (
                (tmp1[j, i] > min(fitotc[I11 - 1 : I22 + 1]) * J2meV)
                and couter + 1 <= subnumber_e
                and (wfe_pow1[i22 - i11 - 1] > 1e-1)
            ):
                KPV1[j, couter] = tmp1[j, i]
                V11[j, i1:i2, couter] = V1[j, i1:i2, i]
                couter += 1
        if couter > subnumber_e:
            print("For this QW, the number confined states of e-levels is: ", couter)
    KPV2 = np.zeros((model.N_wells_virtual, subnumber_h))

    V22 = np.zeros((model.N_wells_virtual, 3 * n_max, 3 * n_max))
    n_max_general3 = np.zeros(model.N_wells_virtual, int)
    wfh_general3 = np.zeros((model.N_wells_virtual, n_max, n_max))
    for k in range(1, model.N_wells_virtual - 1):
        I1, I2, I11, I22 = amort_wave(k, Well_boundary, n_max)
        i_1 = I2 - I1

        V2 = np.zeros((model.N_wells_virtual, i_1 * 3, i_1 * 3))
        tmp = np.zeros((model.N_wells_virtual, i_1 * 3))
        # HUPMAT3_general=np.zeros((i_1*3,i_1*3))
        HUPMAT3_general_2 = np.zeros((i_1 * 3, i_1 * 3))
        i1 = I1 - I1
        i2 = I2 - I1
        HUPMAT3_general_2 = HUPMAT3_reduced_list[k - 1]
        HUPMAT3_general_2 = VBMAT_V_2(HUPMAT3_general_2, fitot, RATIO, i_1, I1, UNIM)
        la2, v2 = linalg.eigh(HUPMAT3_general_2)
        tmp[k, i1 : i2 * 3] = -la2 / RATIO * J2meV
        V2[k, i1 : i2 * 3, i1 : i2 * 3] = v2

        if max(tmp[k, 0:subnumber_h]) > max(fitot[I11:I22]) * J2meV and 1 == 2:
            logger.warning(
                ":You may experience convergence problem due to unconfined states."
            )
        i11 = I11 - I1
        i22 = I22 - I1
        n_max_general3[k] = int(I2 - I1)
        wfh1s3 = np.zeros((i2, 3, n_max_general3[k]))
        maxwfh = np.zeros((i2, 3))
        couter1 = 0
        for i in range(i1, i2):
            for kk in range(0, 3):
                wfh1s3[i, kk, :] = V2[
                    k, kk * n_max_general3[k] : (kk + 1) * n_max_general3[k], i
                ]
                wfh_pow = np.cumsum(wfh1s3[i, kk, :] * wfh1s3[i, kk, :])
                maxwfh[i, kk] = wfh_pow[n_max_general3[k] - 1]
            if np.argmax(maxwfh[i, :]) == 0:
                wfh_general3[k, i, 0 : n_max_general3[k]] = wfh1s3[
                    i, np.argmax(maxwfh[i, :]), :
                ]
            elif np.argmax(maxwfh[i, :]) == 1:
                wfh_general3[k, i, 0 : n_max_general3[k]] = wfh1s3[
                    i, np.argmax(maxwfh[i, :]), :
                ]
            else:
                wfh_general3[k, i, 0 : n_max_general3[k]] = wfh1s3[
                    i, np.argmax(maxwfh[i, :]), :
                ]
            wfh_pow1 = np.cumsum(
                wfh_general3[k, i, i11:i22] * wfh_general3[k, i, i11:i22]
            )  # tmp[k,i]<max(fitot[I11-1:I22+1])*J2meV and

            if (
                (tmp[k, i] > min(fitot[I11 - 1 : I22 + 1]) * J2meV)
                and couter1 + 1 <= subnumber_h
                and (wfh_pow1[i22 - i11 - 1] > 1e-1)
            ):
                # print(wfh_pow1[i22-i11-1],'!=0')
                # print(max(fitot[I11-1:I22+1])*J2meV ,'>',tmp[j,i],'>',min(fitot[I11-1:I22+1])*J2meV)
                KPV2[k, couter1] = tmp[k, i]
                V22[k, i1 : i2 * 3, couter1] = V2[k, i1 : i2 * 3, i]
                couter1 += 1
        if couter1 > subnumber_h:
            print("For this QW, the number confined states of h-levels is: ", couter1)
    """
    for j in range(1,model.N_wells_virtual-1):            
        for i in range(0,subnumber_h,1):
            KPV2[j,i]=tmp[j,i]    
    """
    return KPV1, V11, KPV2, V22


def Main_Str_Array(model):
    n_max = model.n_max
    HUPMAT1 = np.zeros((n_max * 3, n_max * 3))
    HUPMATC1 = np.zeros((n_max, n_max))
    x_max = model.dx * n_max
    m_hh, m_lh, m_so, VNIT, ZETA, CNIT, Ppz_Psp, EPC, pol_surf_char = Strain_and_Masses(
        model
    )
    UNIM = np.identity(n_max)
    RATIO = m_e / hbar ** 2 * (x_max) ** 2
    AC1 = (n_max + 1) ** 2
    AP1, AP2, AP3, AP4, AP5, AP6, FH, FL, FSO, Pce, GDELM, DEL3, DEL1, DEL2 = qsv(
        model.GA1,
        model.GA2,
        model.GA3,
        RATIO,
        VNIT,
        ZETA,
        CNIT,
        AC1,
        n_max,
        model.delta,
        model.A1,
        model.A2,
        model.A3,
        model.A4,
        model.A5,
        model.A6,
        model.delta_so,
        model.delta_cr,
        model.mat_crys_strc,
    )
    KP = 0.0
    KPINT = 0.01
    mat_crys = str(model.mat_crys_strc).lower()
    if "zincblende" in mat_crys and (model.N_wells_virtual - 2 != 0):
        HUPMAT1 = VBMAT1(
            KP,
            AP1,
            AP2,
            AP3,
            AP4,
            AP5,
            AP6,
            FH,
            FL,
            FSO,
            GDELM,
            x_max,
            n_max,
            AC1,
            UNIM,
            KPINT,
        )
        HUPMATC1 = CBMAT(KP, Pce, model.cb_meff / m_e, x_max, n_max, AC1, UNIM, KPINT)
    elif "wurtzite" in mat_crys and (model.N_wells_virtual - 2 != 0):
        HUPMAT1 = -VBMAT2(
            KP,
            AP1,
            AP2,
            AP3,
            AP4,
            AP5,
            AP6,
            FH,
            FL,
            x_max,
            n_max,
            AC1,
            UNIM,
            KPINT,
            DEL3,
            DEL1,
            DEL2,
        )
        HUPMATC1 = CBMAT(KP, Pce, model.cb_meff / m_e, x_max, n_max, AC1, UNIM, KPINT)
    return HUPMAT1, HUPMATC1, m_hh, m_lh, m_so, Ppz_Psp, pol_surf_char


def Schro(
    HUPMAT3_reduced_list,
    HUPMATC1,
    subnumber_h,
    subnumber_e,
    fitot,
    fitotc,
    model,
    Well_boundary,
    UNIM,
    RATIO,
    m_hh,
    m_lh,
    m_so,
    n_max,
):
    # V1=np.zeros((model.N_wells_virtual,n_max,n_max))
    # V2=np.zeros((model.N_wells_virtual,n_max*3,n_max*3))
    n_max_general = np.zeros(model.N_wells_virtual, dtype=int)
    wfh_general = np.zeros((model.N_wells_virtual, subnumber_h, n_max))
    wfe_general = np.zeros((model.N_wells_virtual, subnumber_e, n_max))
    meff_statec_general = np.zeros((model.N_wells_virtual, subnumber_e))
    meff_state_general = np.zeros((model.N_wells_virtual, subnumber_h))
    E_statec_general, V1, E_state_general, V2 = calc_E_state_general(
        HUPMAT3_reduced_list,
        HUPMATC1,
        subnumber_h,
        subnumber_e,
        fitot,
        fitotc,
        model,
        Well_boundary,
        UNIM,
        RATIO,
    )
    for j in range(1, model.N_wells_virtual - 1):
        wfh_general_tmp = np.zeros((model.N_wells_virtual, subnumber_h, n_max))
        wfe_general_tmp = np.zeros((model.N_wells_virtual, subnumber_e, n_max))
        wfh_general_tmp, wfe_general_tmp, list, n_max_general = wave_func_tri(
            j, Well_boundary, n_max, V1, V2, subnumber_h, subnumber_e, model
        )
        wfh_general[j, :, :] += wfh_general_tmp[j, :, :]
        wfe_general[j, :, :] += wfe_general_tmp[j, :, :]
        meff_statec, meff_state = calc_meff_state_general(
            wfh_general[j, :, :],
            wfe_general[j, :, :],
            model,
            fitotc,
            E_statec_general[j, :],
            list,
            m_hh,
            m_lh,
            m_so,
            int(n_max_general[j]),
            j,
            Well_boundary,
            n_max,
        )
        meff_statec_general[j, :], meff_state_general[j, :] = meff_statec, meff_state
    return (
        E_statec_general,
        E_state_general,
        wfe_general,
        wfh_general,
        meff_statec_general,
        meff_state_general,
    )


def Poisson_Schrodinger(model):
    """Performs a self-consistent Poisson-Schrodinger calculation of a 1d quantum well structure.
    Model is an object with the following attributes:
    fi_e - Bandstructure potential (J) (array, len n_max)
    cb_meff - conduction band effective mass (kg)(array, len n_max)
    eps - dielectric constant (including eps0) (array, len n_max)
    dop - doping distribution (m**-3) ( array, len n_max)
    Fapp - Applied field (Vm**-1)
    T - Temperature (K)
    comp_scheme - simulation scheme (currently unused)
    subnumber_e - number of subbands for look for in the conduction band
    dx - grid spacing (m)
    n_max - number of points.
    """
    fi_e = model.fi_e
    cb_meff = model.cb_meff
    eps = model.eps
    dop = model.dop
    Fapp = model.Fapp
    vmax = model.vmax
    vmin = model.vmin
    Each_Step = model.Each_Step
    surface = model.surface
    T = model.T
    comp_scheme = model.comp_scheme
    subnumber_h = model.subnumber_h
    subnumber_e = model.subnumber_e
    dx = model.dx
    n_max = model.n_max
    if comp_scheme in (4, 5, 6):
        logger.error(
            """Aestimo doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        sys.exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """Aestimo doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        sys.exit()
    fi_h = model.fi_h
    N_wells_virtual = model.N_wells_virtual
    Well_boundary = model.Well_boundary
    Ppz_Psp = np.zeros(n_max)
    """
    HUPMAT1=np.zeros((n_max*3, n_max*3))
    """
    HUPMATC1 = np.zeros((n_max, n_max))

    UNIM = np.identity(n_max)
    x_max = dx * n_max
    RATIO = m_e / hbar ** 2 * (x_max) ** 2
    HUPMAT3_reduced_list = []
    has_quantum = getattr(config, 'quantum_effect', True) and (not getattr(model, 'photovoltaic_mode', False) or getattr(model, 'Quantum_Regions', False))
    if (model.N_wells_virtual - 2 != 0) and has_quantum:
        HUPMAT1, HUPMATC1, m_hh, m_lh, m_so, Ppz_Psp, pol_surf_char = Main_Str_Array(
            model
        )
        for k in range(1, model.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(k, Well_boundary, n_max)
            i_1 = I2 - I1
            HUPMAT3_reduced = np.zeros((i_1 * 3, i_1 * 3))
            i1 = I1 - I1
            i2 = I2 - I1
            HUPMAT3_reduced[i1:i2, i1:i2] = HUPMAT1[I1:I2, I1:I2]
            HUPMAT3_reduced[i1 + i_1 : i2 + i_1, i1:i2] = HUPMAT1[
                I1 + n_max : I2 + n_max, I1:I2
            ]
            HUPMAT3_reduced[i1:i2, i1 + i_1 : i2 + i_1] = HUPMAT1[
                I1:I2, I1 + n_max : I2 + n_max
            ]
            HUPMAT3_reduced[i1 + i_1 : i2 + i_1, i1 + i_1 : i2 + i_1] = HUPMAT1[
                I1 + n_max : I2 + n_max, I1 + n_max : I2 + n_max
            ]
            HUPMAT3_reduced[i1 + i_1 * 2 : i2 + i_1 * 2, i1:i2] = HUPMAT1[
                I1 + n_max * 2 : I2 + n_max * 2, I1:I2
            ]
            HUPMAT3_reduced[i1:i2, i1 + i_1 * 2 : i2 + i_1 * 2] = HUPMAT1[
                I1:I2, I1 + n_max * 2 : I2 + n_max * 2
            ]
            HUPMAT3_reduced[
                i1 + i_1 * 2 : i2 + i_1 * 2, i1 + i_1 * 2 : i2 + i_1 * 2
            ] = HUPMAT1[
                I1 + n_max * 2 : I2 + n_max * 2, I1 + n_max * 2 : I2 + n_max * 2
            ]
            HUPMAT3_reduced[i1 + i_1 : i2 + i_1, i1 + i_1 * 2 : i2 + i_1 * 2] = HUPMAT1[
                I1 + n_max : I2 + n_max, I1 + n_max * 2 : I2 + n_max * 2
            ]
            HUPMAT3_reduced[i1 + i_1 * 2 : i2 + i_1 * 2, i1 + i_1 : i2 + i_1] = HUPMAT1[
                I1 + n_max * 2 : I2 + n_max * 2, I1 + n_max : I2 + n_max
            ]
            HUPMAT3_reduced_list.append(HUPMAT3_reduced)
    else:
        (
            m_hh,
            m_lh,
            m_so,
            VNIT,
            ZETA,
            CNIT,
            Ppz_Psp,
            EPC,
            pol_surf_char,
        ) = Strain_and_Masses(model)
    # Check
    if comp_scheme == 6:
        logger.warning(
            """The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere."""
        )
    # Preparing empty subband energy lists.
    E_state = [0.0] * subnumber_h  # Energies of subbands/levels (meV)
    N_state = [0.0] * subnumber_h  # Number of carriers in subbands
    E_statec = [0.0] * subnumber_e  # Energies of subbands/levels (meV)
    N_statec = [0.0] * subnumber_e  # Number of carriers in subbands
    # Preparing empty subband energy arrays for multiquantum wells.
    E_state_general = np.zeros(
        (model.N_wells_virtual, subnumber_h)
    )  # Energies of subbands/levels (meV)
    N_state_general = np.zeros(
        (model.N_wells_virtual, subnumber_h)
    )  # Number of carriers in subbands
    E_statec_general = np.zeros(
        (model.N_wells_virtual, subnumber_e)
    )  # Energies of subbands/levels (meV)
    N_statec_general = np.zeros(
        (model.N_wells_virtual, subnumber_e)
    )  # Number of carriers in subbands
    meff_statec_general = np.zeros((model.N_wells_virtual, subnumber_e))
    meff_state_general = np.zeros((model.N_wells_virtual, subnumber_h))
    # Creating and Filling material arrays
    xaxis = np.arange(0, n_max) * dx  # metres
    fitot = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potential
    fitotc = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potentia
    # eps = np.zeros(n_max+2)	    #dielectric constant
    # dop = np.zeros(n_max+2)	    #doping distribution
    # sigma = np.zeros(n_max+2)      #charge distribution (donors + free charges)
    # F = np.zeros(n_max+2)          #Electric Field
    # Vapp = np.zeros(n_max+2)       #Applied Electric Potential
    V = np.zeros(n_max)  # Electric Potential

    # Subband wavefunction for holes list. 2-dimensional: [i][j] i:stateno, j:wavefunc
    wfh = np.zeros((subnumber_h, n_max))
    wfe = np.zeros((subnumber_e, n_max))
    wfh_general = np.zeros((model.N_wells_virtual, subnumber_h, n_max))
    wfe_general = np.zeros((model.N_wells_virtual, subnumber_e, n_max))
    (
        E_statec_general0,
        E_state_general0,
        wfe_general0,
        wfh_general0,
        meff_statec_general0,
        meff_state_general0,
    ) = (
        E_statec_general,
        E_state_general,
        wfe_general,
        wfh_general,
        meff_statec_general,
        meff_state_general,
    )
    E_F_general = np.zeros(model.N_wells_virtual)
    sigma_general = np.zeros(n_max)
    F_general = np.zeros(n_max)
    Vnew_general = np.zeros(n_max)
    fi = np.zeros(n_max)
    fi_stat = np.zeros(n_max)
    # Setup the doping
    Ntotal = sum(dop)  # calculating total doping density m-3
    Ntotal2d = Ntotal * dx
    # Add to log
    logger.info("Ntotal2d %g m**-2", Ntotal2d)
    # Applied Field
    Vapp = calc_potn(Fapp * eps0 / eps, model)
    Vapp[n_max - 1] -= Vapp[
        n_max // 2
    ]  # Offsetting the applied field's potential so that it is zero in the centre of the structure.
    # s
    # setting up Ldi and Ld p and n
    Ld_n_p = np.zeros(n_max)
    Ldi = np.zeros(n_max)
    Nc = np.zeros(n_max)
    Nv = np.zeros(n_max)
    vb_meff = np.zeros(n_max)
    ni = np.zeros(n_max)
    n = np.zeros(n_max)
    p = np.zeros(n_max)
    hbark = hbar * 2 * pi
    for i in range(n_max):
        vb_meff[i] = (m_hh[i] ** (3 / 2) + m_lh[i] ** (3 / 2)) ** (2 / 3)
    Nc = 2 * (2 * pi * cb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Nv = 2 * (2 * pi * vb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Half_Eg = np.zeros(n_max)
    Eg_ = np.zeros(n_max)
    ns1 = np.linalg.norm(dop, np.inf)
    ns2 = np.linalg.norm(Ppz_Psp, np.inf)
    ns = max(ns1, ns2)

    offset0 = 0.0
    offset1 = 0.0
    for i in range(n_max):
        ni[i] = sqrt(
            Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T))
        )  # Intrinsic carrier concentration [1/m^3]
        if dop[i] == 1:
            dop[i] *= ni[i]
        dop_val = max(abs(dop[i]), 1e6)
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * dop_val))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ns * ni[i]))
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
        Eg_[i] = fi_e[i] - fi_h[i]

        # fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_e scaled
        # fi_h scaled

    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )

    # STARTING SELF CONSISTENT LOOP
    time2 = time.time()  # timing audit
    iteration = 1  # iteration counter
    # previousE0= 0   (meV) energy of zeroth state for previous iteration(for testing convergence)
    previousfi0 = 0  # (meV) energy of  for previous iteration(for testing convergence)
    fitot = fi_h  # + Vapp #For initial iteration sum bandstructure and applied field
    fitotc = fi_e  # + Vapp
    # initializing Stern damping method variables
    r = 0.0
    w_n_minus_max = 1.0
    w_n_max = 0.0
    w_n = np.zeros(n_max)
    damping_n_plus = 0.1
    damping_n = 0.1
    Ppz_Psp0 = Ppz_Psp
    EF = 0.0

    if config.predic_correc:
        print("Predictor–corrector method is activated")
    while True:
        if model.comp_scheme == 9:
            break
        print("Iteration:", iteration)
        # Add to log
        logger.info("Iteration: %d", iteration)
        if (model.N_wells_virtual - 2 != 0) and has_quantum:
            if config.predic_correc and iteration == 1:
                (
                    E_statec_general,
                    E_state_general,
                    wfe_general,
                    wfh_general,
                    meff_statec_general,
                    meff_state_general,
                ) = Schro(
                    HUPMAT3_reduced_list,
                    HUPMATC1,
                    subnumber_h,
                    subnumber_e,
                    fitot,
                    fitotc,
                    model,
                    Well_boundary,
                    UNIM,
                    RATIO,
                    m_hh,
                    m_lh,
                    m_so,
                    n_max,
                )
            elif not(config.predic_correc):
                (
                    E_statec_general,
                    E_state_general,
                    wfe_general,
                    wfh_general,
                    meff_statec_general,
                    meff_state_general,
                ) = Schro(
                    HUPMAT3_reduced_list,
                    HUPMATC1,
                    subnumber_h,
                    subnumber_e,
                    fitot,
                    fitotc,
                    model,
                    Well_boundary,
                    UNIM,
                    RATIO,
                    m_hh,
                    m_lh,
                    m_so,
                    n_max,
                )
            damping = 0.15  # 0.1 works between high and low doping
        else:
            damping = 1
        n, p, fi, EF, fi_stat = Poisson_equi2(
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
            Ppz_Psp0,
            pol_surf_char,
            ni,
            n_max,
            iteration,
            fi,
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
        )
        #
        if comp_scheme in (0, 1):
            # if we are not self-consistently including Poisson Effects then only do one loop
            break
        """
        # Combine band edge potential with potential due to charge distribution
        # To increase convergence, we calculate a moving average of electric potential 
        #with previous iterations. By dampening the corrective term, we avoid oscillations.
        #tryng new dmping method 
        F. Stern, J. Computational Physics 6, 56 (1970).
        #the extrapolated-convergence-factor method instead of the fixed-convergence-factor method
        """
        Vnew_general = -Vt * q * fi
        w_n = Vnew_general - V
        w_n_max = max(abs(w_n[:])) * J2meV
        r = w_n_max / w_n_minus_max
        w_n_minus_max = w_n_max
        damping_n_plus = damping_n / (1 - abs(r))
        damping_n = damping_n_plus
        if config.Stern_damping:
            V += damping_n_plus * (w_n)
        else:
            V += damping * (w_n)
        fitot = fi_h + V + Vapp
        fitotc = fi_e + V + Vapp
        xaxis = np.arange(0, n_max) * dx
        delta0 = V - Vnew_general
        delta_max0 = max(abs(delta0[:]))
        # print('w_n_max=',w_n_max)
        # print('r=',r)
        # print('damping_n=',damping_n)
        # print('damping_n_plus=',damping_n_plus)
        # print('w_n_minus_max=',w_n_minus_max)
        print("error_potential=", delta_max0 * J2meV, "meV")
        if config.predic_correc:
            delta1 = Vnew_general - previousfi0
            delta_max1 = max(abs(delta1[:]))
            if delta_max1 / q < convergence_test0:  # Convergence test
                # print('error=',abs(E_state_general[1,0]-previousE0)/1e3)
                # if abs(E_state_general[1,0]-previousE0)/1e3 < convergence_test: #Convergence test
                if (model.N_wells_virtual - 2 != 0) and has_quantum:
                    (
                        E_statec_general,
                        E_state_general,
                        wfe_general,
                        wfh_general,
                        meff_statec_general,
                        meff_state_general,
                    ) = Schro(
                        HUPMAT3_reduced_list,
                        HUPMATC1,
                        subnumber_h,
                        subnumber_e,
                        fitot,
                        fitotc,
                        model,
                        Well_boundary,
                        UNIM,
                        RATIO,
                        m_hh,
                        m_lh,
                        m_so,
                        n_max,
                    )

                break
            elif iteration >= max_iterations:  # Iteration limit
                logger.warning("Have reached maximum number of iterations")
                break
            else:
                iteration += 1
                previousfi0 = V
        else:
            delta1 = Vnew_general - previousfi0
            delta_max1 = max(abs(delta1[:]))
            if delta_max1 / q < convergence_test0:  # Convergence test
                if (model.N_wells_virtual - 2 != 0) and has_quantum:
                    (
                        E_statec_general,
                        E_state_general,
                        wfe_general,
                        wfh_general,
                        meff_statec_general,
                        meff_state_general,
                    ) = Schro(
                        HUPMAT3_reduced_list,
                        HUPMATC1,
                        subnumber_h,
                        subnumber_e,
                        fitot,
                        fitotc,
                        model,
                        Well_boundary,
                        UNIM,
                        RATIO,
                        m_hh,
                        m_lh,
                        m_so,
                        n_max,
                    )

                break
            elif iteration >= max_iterations:  # Iteration limit
                logger.warning("Have reached maximum number of iterations")
                break
            else:
                iteration += 1
                previousfi0 = V
                # END OF SELF-CONSISTENT LOOP
    (
        Ec_result,
        Ev_result,
        ro_result,
        el_field1_result,
        el_field2_result,
        nf_result,
        pf_result,
        fi_result,
    ) = Write_results_equi2(ns, fitotc, fitot, Vt, q, ni, n, p, dop, dx, Ldi, fi, n_max)
    time3 = time.time()  # timing audit
    # Add to log
    logger.info("calculation time  %g s", (time3 - time2))

    class Results:
        pass

    results = Results()
    results.N_wells_virtual = N_wells_virtual
    results.Well_boundary = Well_boundary
    results.xaxis = xaxis
    results.wfh = wfh
    results.wfe = wfe
    results.fitot = fitot
    results.fitotc = fitotc
    results.fi_e = fi_e
    results.fi_h = fi_h
    # results.sigma = sigma
    results.sigma_general = sigma_general
    # results.F = F
    results.V = V
    results.E_state = E_state
    results.N_state = N_state
    # results.meff_state = meff_state
    results.E_statec = E_statec
    results.N_statec = N_statec
    # results.meff_statec = meff_statec
    results.F_general = F_general
    results.E_state_general = E_state_general
    results.N_state_general = N_state_general
    results.meff_state_general = meff_state_general
    results.E_statec_general = E_statec_general
    results.N_statec_general = N_statec_general
    results.meff_statec_general = meff_statec_general
    results.wfh_general = wfh_general
    results.wfe_general = wfe_general

    results.E_state_general0 = E_state_general0
    results.E_statec_general0 = E_statec_general0
    results.meff_state_general0 = meff_state_general0
    results.meff_statec_general0 = meff_statec_general0
    results.wfh_general0 = wfh_general0
    results.wfe_general0 = wfe_general0
    results.Fapp = Fapp
    results.T = T
    # results.E_F = E_F
    results.E_F_general = E_F_general
    results.dx = dx
    results.subnumber_h = subnumber_h
    results.subnumber_e = subnumber_e
    results.Ntotal2d = Ntotal2d
    ########################
    results.Ec_result = Ec_result
    results.Ev_result = Ev_result
    results.ro_result = ro_result
    results.el_field1_result = el_field1_result
    results.el_field2_result = el_field2_result
    results.nf_result = nf_result
    results.pf_result = pf_result
    results.fi_result = fi_result
    results.EF = EF
    results.HUPMAT3_reduced_list = HUPMAT3_reduced_list
    results.m_hh = m_hh
    results.m_lh = m_lh
    results.m_so = m_so
    results.Ppz_Psp = Ppz_Psp
    results.pol_surf_char = pol_surf_char
    results.HUPMATC1 = HUPMATC1
    ##########################
    return results

def Poisson_Schrodinger_new(model):
    """Performs a self-consistent Poisson-Schrodinger calculation of a 1d quantum well structure.
    Model is an object with the following attributes:
    fi_e - Bandstructure potential (J) (array, len n_max)
    cb_meff - conduction band effective mass (kg)(array, len n_max)
    eps - dielectric constant (including eps0) (array, len n_max)
    dop - doping distribution (m**-3) ( array, len n_max)
    Fapp - Applied field (Vm**-1)
    T - Temperature (K)
    comp_scheme - simulation scheme (currently unused)
    subnumber_e - number of subbands for look for in the conduction band
    dx - grid spacing (m)
    n_max - number of points.
    """
    fi_e = model.fi_e
    cb_meff = model.cb_meff
    eps = model.eps
    dop = model.dop
    Fapp = model.Fapp
    surface = model.surface
    T = model.T
    comp_scheme = model.comp_scheme
    subnumber_h = model.subnumber_h
    subnumber_e = model.subnumber_e
    dx = model.dx
    n_max = model.n_max
    if comp_scheme in (4, 5, 6):
        logger.error(
            """Aestimo doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        sys.exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """Aestimo doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        sys.exit()
    fi_h = model.fi_h
    N_wells_virtual = model.N_wells_virtual
    Well_boundary = model.Well_boundary
    Ppz_Psp = np.zeros(n_max)
    HUPMATC1 = np.zeros((n_max, n_max))
    UNIM = np.identity(n_max)
    x_max = dx * n_max
    RATIO = m_e / hbar ** 2 * (x_max) ** 2
    HUPMAT3_reduced_list = []
    if model.N_wells_virtual - 2 != 0:
        HUPMAT1, HUPMATC1, m_hh, m_lh, m_so, Ppz_Psp, pol_surf_char = Main_Str_Array(
            model
        )
        for k in range(1, model.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(k, Well_boundary, n_max)
            i_1 = I2 - I1
            HUPMAT3_reduced = np.zeros((i_1 * 3, i_1 * 3))
            i1 = I1 - I1
            i2 = I2 - I1
            HUPMAT3_reduced[i1:i2, i1:i2] = HUPMAT1[I1:I2, I1:I2]
            HUPMAT3_reduced[i1 + i_1 : i2 + i_1, i1:i2] = HUPMAT1[
                I1 + n_max : I2 + n_max, I1:I2
            ]
            HUPMAT3_reduced[i1:i2, i1 + i_1 : i2 + i_1] = HUPMAT1[
                I1:I2, I1 + n_max : I2 + n_max
            ]
            HUPMAT3_reduced[i1 + i_1 : i2 + i_1, i1 + i_1 : i2 + i_1] = HUPMAT1[
                I1 + n_max : I2 + n_max, I1 + n_max : I2 + n_max
            ]
            HUPMAT3_reduced[i1 + i_1 * 2 : i2 + i_1 * 2, i1:i2] = HUPMAT1[
                I1 + n_max * 2 : I2 + n_max * 2, I1:I2
            ]
            HUPMAT3_reduced[i1:i2, i1 + i_1 * 2 : i2 + i_1 * 2] = HUPMAT1[
                I1:I2, I1 + n_max * 2 : I2 + n_max * 2
            ]
            HUPMAT3_reduced[
                i1 + i_1 * 2 : i2 + i_1 * 2, i1 + i_1 * 2 : i2 + i_1 * 2
            ] = HUPMAT1[
                I1 + n_max * 2 : I2 + n_max * 2, I1 + n_max * 2 : I2 + n_max * 2
            ]
            HUPMAT3_reduced[i1 + i_1 : i2 + i_1, i1 + i_1 * 2 : i2 + i_1 * 2] = HUPMAT1[
                I1 + n_max : I2 + n_max, I1 + n_max * 2 : I2 + n_max * 2
            ]
            HUPMAT3_reduced[i1 + i_1 * 2 : i2 + i_1 * 2, i1 + i_1 : i2 + i_1] = HUPMAT1[
                I1 + n_max * 2 : I2 + n_max * 2, I1 + n_max : I2 + n_max
            ]
            HUPMAT3_reduced_list.append(HUPMAT3_reduced)
    else:
        (
            m_hh,
            m_lh,
            m_so,
            VNIT,
            ZETA,
            CNIT,
            Ppz_Psp,
            EPC,
            pol_surf_char,
        ) = Strain_and_Masses(model)
    # Check
    if comp_scheme == 6:
        logger.warning(
            """The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere."""
        )
    # Preparing empty subband energy lists.
    E_state = [0.0] * subnumber_h  # Energies of subbands/levels (meV)
    N_state = [0.0] * subnumber_h  # Number of carriers in subbands
    E_statec = [0.0] * subnumber_e  # Energies of subbands/levels (meV)
    N_statec = [0.0] * subnumber_e  # Number of carriers in subbands
    # Preparing empty subband energy arrays for multiquantum wells.
    E_state_general = np.zeros(
        (model.N_wells_virtual, subnumber_h)
    )  # Energies of subbands/levels (meV)
    N_state_general = np.zeros(
        (model.N_wells_virtual, subnumber_h)
    )  # Number of carriers in subbands
    E_statec_general = np.zeros(
        (model.N_wells_virtual, subnumber_e)
    )  # Energies of subbands/levels (meV)
    N_statec_general = np.zeros(
        (model.N_wells_virtual, subnumber_e)
    )  # Number of carriers in subbands
    meff_statec_general = np.zeros((model.N_wells_virtual, subnumber_e))
    meff_state_general = np.zeros((model.N_wells_virtual, subnumber_h))
    # Creating and Filling material arrays
    xaxis = np.arange(0, n_max) * dx  # metres
    fitot = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potential
    fitotc = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potentia
    # eps = np.zeros(n_max+2)	    #dielectric constant
    # dop = np.zeros(n_max+2)	    #doping distribution
    # sigma = np.zeros(n_max+2)      #charge distribution (donors + free charges)
    # F = np.zeros(n_max+2)          #Electric Field
    # Vapp = np.zeros(n_max+2)       #Applied Electric Potential
    V = np.zeros(n_max)  # Electric Potential

    # Subband wavefunction for holes list. 2-dimensional: [i][j] i:stateno, j:wavefunc
    wfh = np.zeros((subnumber_h, n_max))
    wfe = np.zeros((subnumber_e, n_max))
    wfh_general = np.zeros((model.N_wells_virtual, subnumber_h, n_max))
    wfe_general = np.zeros((model.N_wells_virtual, subnumber_e, n_max))
    (
        E_statec_general0,
        E_state_general0,
        wfe_general0,
        wfh_general0,
        meff_statec_general0,
        meff_state_general0,
    ) = (
        E_statec_general,
        E_state_general,
        wfe_general,
        wfh_general,
        meff_statec_general,
        meff_state_general,
    )
    E_F_general = np.zeros(model.N_wells_virtual)
    sigma_general = np.zeros(n_max)
    F_general = np.zeros(n_max)
    Vnew_general = np.zeros(n_max)
    fi = np.zeros(n_max)
    fi_stat = np.zeros(n_max)
    # Setup the doping
    Ntotal = sum(dop)  # calculating total doping density m-3
    Ntotal2d = Ntotal * dx
    # Add to log
    logger.info("Ntotal2d %g m**-2", Ntotal2d)
    # Applied Field
    Vapp = calc_potn(Fapp * eps0 / eps, model)
    Vapp[n_max - 1] -= Vapp[
        n_max // 2
    ]  # Offsetting the applied field's potential so that it is zero in the centre of the structure.
    # s
    # setting up Ldi and Ld p and n
    Ld_n_p = np.zeros(n_max)
    Ldi = np.zeros(n_max)
    Nc = np.zeros(n_max)
    Nv = np.zeros(n_max)
    vb_meff = np.zeros(n_max)
    ni = np.zeros(n_max)
    n = np.zeros(n_max)
    p = np.zeros(n_max)
    hbark = hbar * 2 * pi
    for i in range(n_max):
        vb_meff[i] = (m_hh[i] ** (3 / 2) + m_lh[i] ** (3 / 2)) ** (2 / 3)
    Nc = 2 * (2 * pi * cb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Nv = 2 * (2 * pi * vb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Half_Eg = np.zeros(n_max)
    Eg_ = np.zeros(n_max)
    ns1 = np.linalg.norm(dop, np.inf)
    ns2 = np.linalg.norm(Ppz_Psp, np.inf)
    ns = max(ns1, ns2)

    for i in range(n_max):
        ni[i] = sqrt(
            Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T))
        )  # Intrinsic carrier concentration [1/m^3]
        if dop[i] == 1:
            dop[i] *= ni[i]
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * abs(dop[i])))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ns * ni[i]))
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
        Eg_[i] = fi_e[i] - fi_h[i]

        # fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_e scaled
        # fi_h scaled

    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )

    # STARTING SELF CONSISTENT LOOP
    time2 = time.time()  # timing audit
    iteration = 1  # iteration counter
    # previousE0= 0   (meV) energy of zeroth state for previous iteration(for testing convergence)
    previousfi0 = 0  # (meV) energy of  for previous iteration(for testing convergence)
    fitot = fi_h  # + Vapp #For initial iteration sum bandstructure and applied field
    fitotc = fi_e  # + Vapp
    # initializing Stern damping method variables
    r = 0.0
    w_n_minus_max = 1.0
    w_n_max = 0.0
    w_n = np.zeros(n_max)
    damping_n_plus = 0.1
    damping_n = 0.1
    
    # Global Switch to disable polarization
    if getattr(model, 'enable_polarization', True) == False:
        if 'Ppz_Psp' in locals():
            Ppz_Psp = np.zeros_like(Ppz_Psp)
            logger.info("Polarization effects disabled by configuration (Global Switch).")
            
    Ppz_Psp0 = Ppz_Psp
    EF = 0.0

    l2 = (Vt * eps[0 : n_max - 1]) / (q * ns * xaxis[n_max - 1] ** 2)
    class data:
        def __init__(self):
            self.l2 = l2
            self.dop = dop
            self.V = V
            self.n = n
            self.p = p
            self.Ppz_Psp = Ppz_Psp
            self.E_state_general = E_state_general
            self.meff_state_general = meff_state_general
            self.E_statec_general = E_statec_general
            self.meff_statec_general = meff_statec_general
            self.wfh_general = wfh_general
            self.wfe_general = wfe_general
            self.ni=ni

    idata = data()
    odata = data()
    toll = 1e-3
    maxit = 10
    ptoll = 1e-10
    pmaxit = 30
    verbose = 0
    v_Nnodes=np.arange(n_max)
    idata.l2 = (Vt * eps[0 : n_max - 1]) / (q * ns )#* xaxis[n_max - 1] ** 2
    idata.ni = ni / ns
    idata.dop = dop / ns
    idata.Ppz_Psp = Ppz_Psp / ns
    if config.predic_correc:
        print("Predictor–corrector method is activated")
    while True:
        if model.comp_scheme == 9:
            break
        print("Iteration:", iteration)
        # Add to log
        logger.info("Iteration: %d", iteration)
        if model.N_wells_virtual - 2 != 0:
            if config.predic_correc and iteration == 1:
                (
                    idata.E_statec_general,
                    idata.E_state_general,
                    idata.wfe_general,
                    idata.wfh_general,
                    idata.meff_statec_general,
                    idata.meff_state_general,
                ) = Schro(
                    HUPMAT3_reduced_list,
                    HUPMATC1,
                    subnumber_h,
                    subnumber_e,
                    fitot,
                    fitotc,
                    model,
                    Well_boundary,
                    UNIM,
                    RATIO,
                    m_hh,
                    m_lh,
                    m_so,
                    n_max,
                )
            elif not(config.predic_correc):
                (
                    idata.E_statec_general,
                    idata.E_state_general,
                    idata.wfe_general,
                    idata.wfh_general,
                    idata.meff_statec_general,
                    idata.meff_state_general,
                ) = Schro(
                    HUPMAT3_reduced_list,
                    HUPMATC1,
                    subnumber_h,
                    subnumber_e,
                    fitot,
                    fitotc,
                    model,
                    Well_boundary,
                    UNIM,
                    RATIO,
                    m_hh,
                    m_lh,
                    m_so,
                    n_max,
                )
            damping = 0.15  # 0.1 works between high and low doping
        else:
            damping = 1
        [fi,n,p,fi_stat] =DDGnlpoisson_new (idata,xaxis,v_Nnodes,fi,n,p,ptoll,pmaxit,verbose,fi_e,fi_h,model,Vt,surface,fi_stat,iteration,ns)
        #
        if comp_scheme in (0, 1):
            # if we are not self-consistently including Poisson Effects then only do one loop
            break
        """
        # Combine band edge potential with potential due to charge distribution
        # To increase convergence, we calculate a moving average of electric potential 
        #with previous iterations. By dampening the corrective term, we avoid oscillations.
        #tryng new dmping method 
        F. Stern, J. Computational Physics 6, 56 (1970).
        #the extrapolated-convergence-factor method instead of the fixed-convergence-factor method
        """
        Vnew_general = -Vt * q * fi
        w_n = Vnew_general - V
        w_n_max = max(abs(w_n[:])) * J2meV
        r = w_n_max / w_n_minus_max
        w_n_minus_max = w_n_max
        damping_n_plus = damping_n / (1 - abs(r))
        damping_n = damping_n_plus
        if config.Stern_damping:
            V += damping_n_plus * (w_n)
        else:
            V += damping * (w_n)
        fitot = fi_h + V + Vapp
        fitotc = fi_e + V + Vapp
        xaxis = np.arange(0, n_max) * dx
        delta0 = V - Vnew_general
        delta_max0 = max(abs(delta0[:]))
        # print('w_n_max=',w_n_max)
        # print('r=',r)
        # print('damping_n=',damping_n)
        # print('damping_n_plus=',damping_n_plus)
        # print('w_n_minus_max=',w_n_minus_max)
        print("error_potential=", delta_max0 * J2meV, "meV")
        if config.predic_correc:
            delta1 = Vnew_general - previousfi0
            delta_max1 = max(abs(delta1[:]))
            if delta_max1 / q < convergence_test0 :  # Convergence test
                # print('error=',abs(E_state_general[1,0]-previousE0)/1e3)
                # if abs(E_state_general[1,0]-previousE0)/1e3 < convergence_test: #Convergence test
                if model.N_wells_virtual - 2 != 0:
                    (
                        E_statec_general,
                        E_state_general,
                        wfe_general,
                        wfh_general,
                        meff_statec_general,
                        meff_state_general,
                    ) = Schro(
                        HUPMAT3_reduced_list,
                        HUPMATC1,
                        subnumber_h,
                        subnumber_e,
                        fitot,
                        fitotc,
                        model,
                        Well_boundary,
                        UNIM,
                        RATIO,
                        m_hh,
                        m_lh,
                        m_so,
                        n_max,
                    )

                break
            elif iteration >= max_iterations:  # Iteration limit
                logger.warning("Have reached maximum number of iterations")
                break
            else:
                iteration += 1
                previousfi0 = V
        else:
            delta1 = Vnew_general - previousfi0
            delta_max1 = max(abs(delta1[:]))
            if delta_max1 / q < convergence_test0:  # Convergence test
                if model.N_wells_virtual - 2 != 0:
                    (
                        E_statec_general,
                        E_state_general,
                        wfe_general,
                        wfh_general,
                        meff_statec_general,
                        meff_state_general,
                    ) = Schro(
                        HUPMAT3_reduced_list,
                        HUPMATC1,
                        subnumber_h,
                        subnumber_e,
                        fitot,
                        fitotc,
                        model,
                        Well_boundary,
                        UNIM,
                        RATIO,
                        m_hh,
                        m_lh,
                        m_so,
                        n_max,
                    )

                break
            elif iteration >= max_iterations:  # Iteration limit
                logger.warning("Have reached maximum number of iterations")
                break
            else:
                iteration += 1
                previousfi0 = V
                # END OF SELF-CONSISTENT LOOP
    (
        Ec_result,
        Ev_result,
        ro_result,
        el_field1_result,
        el_field2_result,
        nf_result,
        pf_result,
        fi_result,
    ) = Write_results_equi2(ns, fitotc, fitot, Vt, q, ni, n*ns/ni, p*ns/ni, dop, dx, Ldi, fi, n_max)
    time3 = time.time()  # timing audit
    # Add to log
    logger.info("calculation time  %g s", (time3 - time2))

    class Results:
        pass

    results = Results()
    results.N_wells_virtual = N_wells_virtual
    results.Well_boundary = Well_boundary
    results.xaxis = xaxis
    results.wfh = wfh
    results.wfe = wfe
    results.fitot = fitot
    results.fitotc = fitotc
    results.fi_e = fi_e
    results.fi_h = fi_h
    # results.sigma = sigma
    results.sigma_general = sigma_general
    # results.F = F
    results.V = V
    results.E_state = E_state
    results.N_state = N_state
    # results.meff_state = meff_state
    results.E_statec = E_statec
    results.N_statec = N_statec
    # results.meff_statec = meff_statec
    results.F_general = F_general
    results.E_state_general = E_state_general
    results.N_state_general = N_state_general
    results.meff_state_general = meff_state_general
    results.E_statec_general = E_statec_general
    results.N_statec_general = N_statec_general
    results.meff_statec_general = meff_statec_general
    results.wfh_general = wfh_general
    results.wfe_general = wfe_general

    results.E_state_general0 = E_state_general0
    results.E_statec_general0 = E_statec_general0
    results.meff_state_general0 = meff_state_general0
    results.meff_statec_general0 = meff_statec_general0
    results.wfh_general0 = wfh_general0
    results.wfe_general0 = wfe_general0
    results.Fapp = Fapp
    results.T = T
    # results.E_F = E_F
    results.E_F_general = E_F_general
    results.dx = dx
    results.subnumber_h = subnumber_h
    results.subnumber_e = subnumber_e
    results.Ntotal2d = Ntotal2d
    ########################
    results.Ec_result = Ec_result
    results.Ev_result = Ev_result
    results.ro_result = ro_result
    results.el_field1_result = el_field1_result
    results.el_field2_result = el_field2_result
    results.nf_result = nf_result
    results.pf_result = pf_result
    results.fi_result = fi_result
    results.EF = EF
    results.HUPMAT3_reduced_list = HUPMAT3_reduced_list
    results.m_hh = m_hh
    results.m_lh = m_lh
    results.m_so = m_so
    results.Ppz_Psp = Ppz_Psp
    results.pol_surf_char = pol_surf_char
    results.HUPMATC1 = HUPMATC1
    ##########################
    return results

def Poisson_Schrodinger_DD(result, model):
    # Initialize all potential result variables to avoid UnboundLocalError
    n_max = int(model.n_max)
    Va_t = np.zeros(1)
    Efn_result = Efp_result = Ei_result = Ec_result = Ev_result = np.zeros(n_max)
    ro_result = el_field1_result = el_field2_result = nf_result = pf_result = np.zeros(n_max)
    fi_result = np.zeros(n_max)
    av_curr = np.zeros(1)
    EF = fi_va = Ec_result_ = Ev_result_ = None
    Total_Steps = 1

    fi = result.fi_result
    E_state_general = result.E_state_general
    meff_state_general = result.meff_state_general
    E_statec_general = result.E_statec_general
    meff_statec_general = result.meff_statec_general
    wfh_general = result.wfh_general
    wfe_general = result.wfe_general
    n_max = model.n_max
    dx = model.dx
    HUPMAT3_reduced_list = result.HUPMAT3_reduced_list
    HUPMATC1 = result.HUPMATC1
    m_hh = result.m_hh
    m_lh = result.m_lh
    m_so = result.m_so
    Ppz_Psp = result.Ppz_Psp
    pol_surf_char = result.pol_surf_char
    fi_e = model.fi_e
    cb_meff = model.cb_meff
    eps = model.eps
    dop = model.dop
    Fapp = model.Fapp
    vmax = model.vmax
    vmin = model.vmin
    Each_Step = model.Each_Step
    surface = model.surface
    T = model.T
    comp_scheme = model.comp_scheme
    subnumber_h = model.subnumber_h
    subnumber_e = model.subnumber_e
    TAUN0 = model.TAUN0
    TAUP0 = model.TAUP0
    mun0 = model.mun0
    mup0 = model.mup0
    BETAN = model.BETAN
    BETAP = model.BETAP
    VSATN = model.VSATN
    VSATP = model.VSATP
    Cn0 = model.Cn0
    Cp0 = model.Cp0
    # Setup Optical Generation Profile
    gen_type = getattr(model, 'generation_type', 'uniform')
    G_optical_val = float(getattr(model, 'G_optical', 0.0))
    G_optical = np.zeros(n_max)
    
    if gen_type == 'uniform':
        G_optical[:] = G_optical_val * 1e6
        logger.info("Using uniform G_optical = %g cm^-3 s^-1", G_optical_val)
    elif gen_type == 'exponential':
        alpha = getattr(model, 'alpha', 1e5) # cm^-1
        alpha_m = alpha * 1e2 # m^-1
        xaxis_local = np.arange(0, n_max) * dx # now in meters
        G_optical = (G_optical_val * 1e6) * np.exp(-alpha_m * xaxis_local)
        logger.info("Using exponential G_optical (alpha = %g cm^-1)", alpha)
    else:
        G_optical[:] = G_optical_val * 1e6
    if comp_scheme in (4, 5, 6):
        logger.error(
            """Aestimo doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        sys.exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """Aestimo doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        sys.exit()
    fi_h = model.fi_h
    N_wells_virtual = model.N_wells_virtual
    Well_boundary = model.Well_boundary
    x_max = dx * n_max
    UNIM = np.identity(n_max)
    RATIO = m_e / hbar ** 2 * (x_max) ** 2

    # Check
    if comp_scheme == 6:
        logger.warning(
            """The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere."""
        )
    # Preparing empty subband energy lists.
    E_state = [0.0] * subnumber_h  # Energies of subbands/levels (meV)
    N_state = [0.0] * subnumber_h  # Number of carriers in subbands
    E_statec = [0.0] * subnumber_e  # Energies of subbands/levels (meV)
    N_statec = [0.0] * subnumber_e  # Number of carriers in subbands
    # Preparing empty subband energy arrays for multiquantum wells.
    N_state_general = np.zeros(
        (model.N_wells_virtual, subnumber_h)
    )  # Number of carriers in subbands
    N_statec_general = np.zeros(
        (model.N_wells_virtual, subnumber_e)
    )  # Number of carriers in subbands

    # Creating and Filling material arrays
    xaxis = np.arange(0, n_max) * dx  # metres
    fitot = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potential
    fitotc = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potentia
    # eps = np.zeros(n_max+2)	    #dielectric constant
    # dop = np.zeros(n_max+2)	    #doping distribution
    # sigma = np.zeros(n_max+2)      #charge distribution (donors + free charges)
    # F = np.zeros(n_max+2)          #Electric Field
    # Vapp = np.zeros(n_max+2)       #Applied Electric Potential
    V = np.zeros(n_max)  # Electric Potential

    # Subband wavefunction for holes list. 2-dimensional: [i][j] i:stateno, j:wavefunc

    wfh = np.zeros((subnumber_h, n_max))
    wfe = np.zeros((subnumber_e, n_max))
    """
    wfh_general = np.zeros((model.N_wells_virtual,subnumber_h,n_max))
    wfe_general = np.zeros((model.N_wells_virtual,subnumber_e,n_max))
    """
    E_F_general = np.zeros(model.N_wells_virtual)
    sigma_general = np.zeros(n_max)
    F_general = np.zeros(n_max)
    Vnew_general = np.zeros(n_max)
    # fi = np.zeros(n_max)
    # Setup the doping
    Ntotal = sum(dop)  # calculating total doping density m-3
    Ntotal2d = Ntotal * dx
    # Add to log
    logger.info("Ntotal2d %g m**-2", Ntotal2d)
    # Applied Field
    # Vapp = calc_potn(Fapp*eps0/eps,model)
    # Vapp[n_max-1] -= Vapp[n_max//2] #Offsetting the applied field's potential so that it is zero in the centre of the structure.
    # s
    # setting up Ldi and Ld p and n
    Ld_n_p = np.zeros(n_max)
    Ldi = np.zeros(n_max)
    Nc = np.zeros(n_max)
    Nv = np.zeros(n_max)
    vb_meff = np.zeros(n_max)
    ni = np.zeros(n_max)
    ni_phys = np.zeros(n_max)
    hbark = hbar * 2 * pi
    Ppz_Psp_tmp = Ppz_Psp
    Ppz_Psp = np.zeros(n_max)
    for i in range(n_max):
        vb_meff[i] = (m_hh[i] ** (3 / 2) + m_lh[i] ** (3 / 2) + m_so[i] ** (3 / 2)) ** (
            2 / 3
        )
    Nc = 2 * (2 * pi * cb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Nv = 2 * (2 * pi * vb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Half_Eg = np.zeros(n_max)
    for i in range(n_max):
        val_ni = sqrt(Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T)))
        ni_phys[i] = val_ni
        # We use a stable reference density (ni_ref) for all dimensionless normalization.
        # 1e18 m^-3 (1e12 cm^-3) is a robust scaling unit for wide-bandgap DD.
        ni[i] = np.maximum(val_ni, 1e18)
        
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * abs(dop[i] + 1e-20)))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ni[i]))
        if dop[i] == 1:
            dop[i] *= ni[i]
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
    
    # Store real ni for physics solvers
    model.ni_phys = ni_phys
        # fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_e scaled
        # fi_h scaled
    n = result.nf_result / ni
    p = result.pf_result / ni
    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )

    # STARTING SELF CONSISTENT LOOP
    time2 = time.time()  # timing audit
    iteration = 1  # iteration counter
    # previousE0= 0   #(meV) energy of zeroth state for previous iteration(for testing convergence)
    # fitot = fi_h + Vapp #For initial iteration sum bandstructure and applied field
    # fitotc = fi_e + Vapp

    #vmin = 0.0
    Total_Steps = int(((vmax - vmin) ) / (Each_Step))+1
    xaxis = np.arange(0, n_max) * dx  # metres
    mup = np.zeros(n_max)
    mun = np.zeros(n_max)
    EF = 0.0
    av_curr = np.zeros(Total_Steps)
    Va_t = np.zeros(Total_Steps)
    Jnim1by2 = np.zeros((Total_Steps, n_max))
    Jnip1by2 = np.zeros((Total_Steps, n_max))
    Jelec = np.zeros((Total_Steps, n_max))
    Jpim1by2 = np.zeros((Total_Steps, n_max))
    Jpip1by2 = np.zeros((Total_Steps, n_max))
    Jhole = np.zeros((Total_Steps, n_max))
    Jtotal = np.zeros((Total_Steps, n_max))
    fi_va = np.zeros((Total_Steps, n_max))
    Ec_result_ = np.zeros((Total_Steps, n_max))
    Ev_result_ = np.zeros((Total_Steps, n_max))
    fi_stat = fi.copy()
    fi[0] -= vmin / Vt
    if Total_Steps < 2 and not getattr(model, 'use_newton_solver', False):
        print("Equilibrium only (Total_Steps < 2)")
    else:
        print("Convergence of the Gummel cycles")
        vindex = 0
        for vindex in range(0, Total_Steps):
            if vindex>int(Total_Steps*4/5):
                Ppz_Psp = Ppz_Psp_tmp
            # Start Va increment loop
            Va = Each_Step * vindex
            if vindex > 0:
                fi[0] -= Each_Step / Vt
            flag_conv_2 = True  # Convergence of the Poisson loop
            #% Initialize the First and Last Node for Poisson's eqn

            Va_t[vindex] = Va+vmin
            logger.info("Voltage Step %d/%d: Va = %g V", vindex + 1, Total_Steps, Va_t[vindex])
            iteration = 1 # Reset iteration counter for each voltage step
            # previousE0 = 2
            while flag_conv_2:
                if iteration % 20 == 0:
                    logger.info("  Iteration %d...", iteration)
                    sys.stdout.flush()
                
                # Hard iteration cap for stability
                max_iter_val = getattr(model, 'dd_max_iterations', 25)
                curr_p_damp = getattr(model, 'poisson_damping', 0.4)
                curr_c_damp = getattr(model, 'continuity_damping', 0.7)                
                if not getattr(model, 'use_newton_solver', False) and iteration > max_iter_val:
                    flag_conv_2 = False
                    break
                    
                if getattr(model, 'use_newton_solver', False):
                    if not hasattr(model, 'newton_solver'):
                        from aeslibs.newton_raphson import CoupledNewtonSolver
                        model.newton_solver = CoupledNewtonSolver(
                            model, n_max, dx, ni, dop, Ldi, Ppz_Psp, pol_surf_char, Nc, Nv, fi_stat, n_stat=n, p_stat=p
                        )
                    
                    # Update mobility for Newton solver
                    mun, mup = Mobility2(
                        mun0, mup0, fi, Vt, Ldi, VSATN, VSATP, BETAN, BETAP, n_max, dx
                    )
                    
                    fi, n, p, newton_ok = model.newton_solver.solve(
                        fi, n, p, mun, mup, TAUN0, TAUP0, Cn0, Cp0, G_optical, iteration, Va=Va_t[vindex]
                    )
                    flag_conv_2 = False
                    
                    model.newton_solver.require_convergence(newton_ok, fi, n, p, Va_t[vindex])

                if not getattr(model, 'use_newton_solver', False):
                    fi, flag_conv_2 = Poisson_non_equi2(
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
                    damping=curr_p_damp,
                )
                
                    #
                    mun, mup = Mobility2(
                        mun0, mup0, fi, Vt, Ldi, VSATN, VSATP, BETAN, BETAP, n_max, dx
                    )
                    
                    ########### END of FIELD Dependant Mobility Calculation ###########
                    n, p = Continuity2(n, p, mun, mup, fi, Vt, Ldi, n_max, dx, TAUN0, TAUP0, ni, G_optical, iteration, model=model, dop=dop, Cn0=Cn0, Cp0=Cp0, damping=curr_c_damp)
                
                # Check for numerical instability
                if not np.all(np.isfinite(n)) or not np.all(np.isfinite(p)):
                    logger.error("  Numerical Instability: Carrier densities reached non-finite values at Va = %g V. Terminating Gummel loop.", Va_t[vindex])
                    flag_conv_2 = False
                    break
                    
                iteration += 1
                ####################### END of HOLE Continuty Solver ###########
                # End of WHILE Loop for Poisson's eqn solver
            if getattr(model, 'use_newton_solver', False) and hasattr(model, 'newton_solver') and model.newton_solver.Jtot is not None:
                Jelec[vindex, :n_max-1] = model.newton_solver.Jn
                Jhole[vindex, :n_max-1] = model.newton_solver.Jp
                Jelec[vindex, -1] = Jelec[vindex, -2]
                Jhole[vindex, -1] = Jhole[vindex, -2]
                Jtotal[vindex, :] = Jelec[vindex, :] + Jhole[vindex, :]
                av_curr[vindex] = model.newton_solver.last_Jtot
            else:
                Jnip1by2, Jnim1by2, Jelec, Jpip1by2, Jpim1by2, Jhole = Current2(
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
                )
                
                # Note: Current2 now returns physical current density in A/m^2 (SI)
                # because dx and ni are SI, and mun is converted internally.
                # Convert to mA/cm^2 for Aestimo GUI/Reports (1 A/m2 = 0.1 mA/cm2)
                Jelec[vindex, :] *= 0.1
                Jhole[vindex, :] *= 0.1
                Jtotal[vindex, :] = Jelec[vindex, :] + Jhole[vindex, :]

            # End of main FOR loop for Va increment.
            fi_va[vindex, :] = fi
            
            # No early stopping — always run the full sweep from vmin to vmax for accurate Voc/Pmax extraction

        for vindex in range(Total_Steps):
            Ec_result_[vindex, :] = fi_e / q - Vt * fi_va[vindex, :]
            Ev_result_[vindex, :] = fi_h / q - Vt * fi_va[vindex, :]

        # Compute av_curr for sequential solver if not already computed by Newton solver
        for vindex in range(Total_Steps):
            if not getattr(model, 'use_newton_solver', False):
                idx_lo = int(0.9 * n_max)
                idx_hi = n_max - 1
                av_curr[vindex] = np.median(Jtotal[vindex, idx_lo:idx_hi])
            
        av_curr = av_curr[:Total_Steps]
        Ec_result_ = Ec_result_[:Total_Steps, :]
        Ev_result_ = Ev_result_[:Total_Steps, :]
        ##########################################################################
        ##                 END OF NON-EQUILIBRIUM  SOLUTION PART                ##
        ##########################################################################
        # Write the results of the simulation in files #
        (
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
        ) = Write_results_non_equi2(
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
        )
        fitot = fi_h - Vt * q * fi
        fitotc = fi_e - Vt * q * fi
        has_quantum = getattr(config, 'quantum_effect', True) and (not getattr(model, 'photovoltaic_mode', False) or getattr(model, 'Quantum_Regions', False))
        if (model.N_wells_virtual - 2 != 0) and has_quantum:

            (
                E_statec_general,
                E_state_general,
                wfe_general,
                wfh_general,
                meff_statec_general,
                meff_state_general,
            ) = Schro(
                HUPMAT3_reduced_list,
                HUPMATC1,
                subnumber_h,
                subnumber_e,
                fitot,
                fitotc,
                model,
                Well_boundary,
                UNIM,
                RATIO,
                m_hh,
                m_lh,
                m_so,
                n_max,
            )
    time3 = time.time()  # timing audit
    # Add to log
    logger.info("calculation time  %g s", (time3 - time2))

    class Results:
        pass

    results = Results()
    results.N_wells_virtual = N_wells_virtual
    results.Well_boundary = Well_boundary
    results.xaxis = xaxis
    results.wfh = wfh
    results.wfe = wfe
    results.wfh_general = wfh_general
    results.wfe_general = wfe_general
    results.fitot = fitot
    results.fitotc = fitotc
    results.fi_e = fi_e
    results.fi_h = fi_h
    # results.sigma = sigma
    results.sigma_general = sigma_general
    # results.F = F
    results.V = V
    results.E_state = E_state
    results.N_state = N_state
    # results.meff_state = meff_state
    results.E_statec = E_statec
    results.N_statec = N_statec
    # results.meff_statec = meff_statec
    results.F_general = F_general
    results.E_state_general = E_state_general
    results.N_state_general = N_state_general
    results.meff_state_general = meff_state_general
    results.E_statec_general = E_statec_general
    results.N_statec_general = N_statec_general
    results.meff_statec_general = meff_statec_general
    results.Fapp = Fapp
    results.T = T
    # results.E_F = E_F
    results.E_F_general = E_F_general
    results.dx = dx
    results.subnumber_h = subnumber_h
    results.subnumber_e = subnumber_e
    results.Ntotal2d = Ntotal2d
    ########################
    results.Va_t = Va_t
    results.Efn_result = Efn_result
    results.Efp_result = Efp_result
    results.Ei_result = Ei_result
    results.av_curr = av_curr
    results.Ec_result = Ec_result
    results.Ev_result = Ev_result
    results.ro_result = ro_result
    results.el_field1_result = el_field1_result
    results.el_field2_result = el_field2_result
    results.nf_result = nf_result
    results.pf_result = pf_result
    results.fi_result = fi_result
    results.EF = EF
    results.Total_Steps = Total_Steps
    results.fi_va = fi_va
    results.Ec_result_ = Ec_result_
    results.Ev_result_ = Ev_result_
    

    
    return results


################################################
def Poisson_Schrodinger_DD_test(result, model):
    """Performs a self-consistent Poisson-Schrodinger calculation of a 1d quantum well structure.
    Model is an object with the following attributes:
    fi_e - Bandstructure potential (J) (array, len n_max)
    cb_meff - conduction band effective mass (kg)(array, len n_max)
    eps - dielectric constant (including eps0) (array, len n_max)
    dop - doping distribution (m**-3) ( array, len n_max)
    Fapp - Applied field (Vm**-1)
    T - Temperature (K)
    comp_scheme - simulation scheme (currently unused)
    subnumber_e - number of subbands for look for in the conduction band
    dx - grid spacing (m)
    n_max - number of points.
    """
    fi = result.fi_result
    E_state_general = result.E_state_general
    meff_state_general = result.meff_state_general
    E_statec_general = result.E_statec_general
    meff_statec_general = result.meff_statec_general
    wfh_general = result.wfh_general
    wfe_general = result.wfe_general
    HUPMAT3_reduced_list = result.HUPMAT3_reduced_list
    E_state_general0 = result.E_state_general0
    meff_state_general0 = result.meff_state_general0
    E_statec_general0 = result.E_statec_general0
    meff_statec_general0 = result.meff_statec_general0
    wfh_general0 = result.wfh_general0
    wfe_general0 = result.wfe_general0
    n_max = model.n_max
    dx = model.dx
    HUPMAT3_reduced_list = result.HUPMAT3_reduced_list
    m_hh = result.m_hh
    m_lh = result.m_lh
    m_so = result.m_so
    Ppz_Psp = result.Ppz_Psp
    pol_surf_char = result.pol_surf_char
    HUPMATC1 = result.HUPMATC1
    fi_e = model.fi_e
    cb_meff = model.cb_meff
    eps = model.eps
    dop = model.dop
    Fapp = model.Fapp
    vmax = model.vmax
    vmin = model.vmin
    Each_Step = model.Each_Step
    surface = model.surface
    T = model.T
    comp_scheme = model.comp_scheme
    subnumber_h = model.subnumber_h
    subnumber_e = model.subnumber_e
    dx = model.dx
    n_max = model.n_max
    TAUN0 = model.TAUN0
    TAUP0 = model.TAUP0
    mun0 = model.mun0
    mup0 = model.mup0
    BETAN = model.BETAN
    BETAP = model.BETAP
    VSATN = model.VSATN
    VSATP = model.VSATP
    if comp_scheme in (4, 5, 6):
        logger.error(
            """Aestimo doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        sys.exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """Aestimo doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        sys.exit()
    fi_h = model.fi_h
    N_wells_virtual = model.N_wells_virtual
    Well_boundary = model.Well_boundary
    

    x_max = dx * n_max
    # Check
    if comp_scheme == 6:
        logger.warning(
            """The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere."""
        )
    # Preparing empty subband energy lists.
    E_state = [0.0] * subnumber_h  # Energies of subbands/levels (meV)
    N_state = [0.0] * subnumber_h  # Number of carriers in subbands
    E_statec = [0.0] * subnumber_e  # Energies of subbands/levels (meV)
    N_statec = [0.0] * subnumber_e  # Number of carriers in subbands
    # Preparing empty subband energy arrays for multiquantum wells.
    N_state_general = np.zeros(
        (model.N_wells_virtual, subnumber_h)
    )  # Number of carriers in subbands
    N_statec_general = np.zeros(
        (model.N_wells_virtual, subnumber_e)
    )  # Number of carriers in subbands

    # Optical Generation Rate
    G_optical_cm3 = getattr(model, 'G_optical', 0.0)
    G_optical = G_optical_cm3 * 1e6 # Convert from cm^-3 s^-1 to m^-3 s^-1 for internal solver
    logger.info("Using G_optical = %g cm^-3 s^-1 (Internal: %g m^-3 s^-1)", G_optical_cm3, G_optical)
    

    # Creating and Filling material arrays
    xaxis = np.arange(0, n_max) * dx  # metres
    fitot = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potential
    fitotc = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potentia
    # eps = np.zeros(n_max+2)	    #dielectric constant
    # dop = np.zeros(n_max+2)	    #doping distribution
    # sigma = np.zeros(n_max+2)      #charge distribution (donors + free charges)
    # F = np.zeros(n_max+2)          #Electric Field
    # Vapp = np.zeros(n_max+2)       #Applied Electric Potential
    V = np.zeros(n_max)  # Electric Potential

    # Subband wavefunction for holes list. 2-dimensional: [i][j] i:stateno, j:wavefunc

    wfh = np.zeros((subnumber_h, n_max))
    wfe = np.zeros((subnumber_e, n_max))
    E_F_general = np.zeros(model.N_wells_virtual)
    sigma_general = np.zeros(n_max)
    F_general = np.zeros(n_max)
    Vnew_general = np.zeros(n_max)
    # fi = np.zeros(n_max)
    # Setup the doping
    Ntotal = sum(dop)  # calculating total doping density m-3
    Ntotal2d = Ntotal * dx
    # Add to log
    logger.info("Ntotal2d %g m**-2", Ntotal2d)
    # Applied Field
    # Vapp = calc_potn(Fapp*eps0/eps,model)
    # Vapp[n_max-1] -= Vapp[n_max//2] #Offsetting the applied field's potential so that it is zero in the centre of the structure.
    # s
    # setting up Ldi and Ld p and n
    Ld_n_p = np.zeros(n_max)
    Ldi = np.zeros(n_max)
    Nc = np.zeros(n_max)
    Nv = np.zeros(n_max)
    vb_meff = np.zeros(n_max)
    ni = np.zeros(n_max)
    hbark = hbar * 2 * pi
    # m_hh,m_lh,m_so,VNIT,ZETA,CNIT,Ppz_Psp,EPC,pol_surf_char=Strain_and_Masses(model)

    Ppz_Psp_tmp = Ppz_Psp
    Ppz_Psp = np.zeros(n_max)
    
    # Define scaling factors for normalization (consistent with Aestimo's internal convention)
    xs = dx
    Vs = Vt
    us = np.max(np.abs(mun0)) if np.max(np.abs(mun0)) > 1e-12 else 0.1

    UNIM = np.identity(n_max)
    x_max = dx * n_max
    RATIO = m_e / hbar ** 2 * (x_max) ** 2
    for i in range(n_max):
        vb_meff[i] = (m_hh[i] ** (3 / 2) + m_lh[i] ** (3 / 2) + m_so[i] ** (3 / 2)) ** (
            2 / 3
        )
    Nc = 2 * (2 * pi * cb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Nv = 2 * (2 * pi * vb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Half_Eg = np.zeros(n_max)
    for i in range(n_max):
        val_ni = sqrt( Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T)) )
        ni[i] = max(val_ni, 1e18)  # Consistent ni_ref scaling
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * abs(dop[i])))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ni[i]))
        # fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
    n = result.nf_result / ni
    p = result.pf_result / ni
    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )
        sys.exit()
    # STARTING SELF CONSISTENT LOOP
    time2 = time.time()  # timing audit
    iteration = 1  # iteration counter
    previousE0 = 0  # (meV) energy of zeroth state for previous iteration(for testing convergence)
    previousfi0 = 0  # (meV) energy of  for previous iteration(for testing convergence)
    fitot = fi_h  # + Vapp #For initial iteration sum bandstructure and applied field
    fitotc = fi_e  # + Vapp
    Total_Steps = int(((vmax - vmin) ) / (Each_Step))+1
    xaxis = np.arange(0, n_max) * dx  # metres
    mup = np.zeros(n_max)
    mun = np.zeros(n_max)
    n_q = np.zeros(n_max)
    p_q = np.zeros(n_max)
    fi_n = np.zeros(n_max)
    fi_p = np.zeros(n_max)
    EF = 0.0
    av_curr = np.zeros(Total_Steps)
    Va_t = np.zeros(Total_Steps)
    Jnim1by2 = np.zeros((Total_Steps, n_max))
    Jnip1by2 = np.zeros((Total_Steps, n_max))
    Jelec = np.zeros((Total_Steps, n_max))
    Jpim1by2 = np.zeros((Total_Steps, n_max))
    Jpip1by2 = np.zeros((Total_Steps, n_max))
    Jhole = np.zeros((Total_Steps, n_max))
    Jtotal = np.zeros((Total_Steps, n_max))
    fi_va = np.zeros((Total_Steps, n_max))
    Ec_result_ = np.zeros((Total_Steps, n_max))
    Ev_result_ = np.zeros((Total_Steps, n_max))
    fi_stat = fi
    fi+=vmin/Vt
    if Total_Steps < 2:
        print("Equilibrium only (Total_Steps < 2)")
    else:
        print("vindex=0")
        print("Convergence of the Gummel cycles")
        vindex = 0
        for vindex in range(0, Total_Steps):
            iteration = 1  # iteration counter
            Ppz_Psp = Ppz_Psp_tmp
            # Start Va increment loop
            Va = Each_Step * vindex
            if vindex == 0:
                fi[0] += 0.0  # Apply potential to Anode (1st node)
            else:
                fi[0] += Each_Step/Vt
            flag_conv_2 = True  # Convergence of the Poisson loop
            Va_t[vindex] = Va+vmin
            logger.info("Voltage Step %d/%d: Va = %g V", vindex + 1, Total_Steps, Va_t[vindex])
            max_poisson_iter = 40
            while flag_conv_2:
                if iteration > max_poisson_iter:
                    logger.warning(f"  [WARN] Voltage Step {vindex + 1}: Poisson-Gummel loop exceeded {max_poisson_iter} iterations. Proceeding anyway.")
                    flag_conv_2 = False
                    break
                    
                fitot = fi_h - Vt * q * fi
                fitotc = fi_e - Vt * q * fi
                if model.N_wells_virtual - 2 != 0:
                    (
                        E_statec_general,
                        E_state_general,
                        wfe_general,
                        wfh_general,
                        meff_statec_general,
                        meff_state_general,
                    ) = Schro(
                        HUPMAT3_reduced_list,
                        HUPMATC1,
                        subnumber_h,
                        subnumber_e,
                        fitot,
                        fitotc,
                        model,
                        Well_boundary,
                        UNIM,
                        RATIO,
                        m_hh,
                        m_lh,
                        m_so,
                        n_max,
                    )
                fi, flag_conv_2, n_q, p_q, fi_n, fi_p = Poisson_non_equi3(
                    vindex,
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
                )
                iteration += 1

                mun, mup = Mobility3(
                    mun0,
                    mup0,
                    fi,
                    fi_n,
                    fi_p,
                    Vt,
                    Ldi,
                    VSATN,
                    VSATP,
                    BETAN,
                    BETAP,
                    n_max,
                    dx,
                    ni,
                    n,
                    p,
                )
                ########### END of FIELD Dependant Mobility Calculation ###########
                # Time scaling for normalized continuity equation
                ts_scaling = (xs**2) / (us * Vs)
                n_old, p_old = n.copy(), p.copy()
                n, p = Continuity3(
                    n, p, mun, mup, fi, fi_n, fi_p, Vt, Ldi, n_max, dx, TAUN0, TAUP0, G_opt=G_optical, ni=ni, ts=ts_scaling, Cn0=Cn0, Cp0=Cp0, model=model
                )
                
            Jnip1by2, Jnim1by2, Jelec, Jpip1by2, Jpim1by2, Jhole = Current2(
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
            )
            # End of main FOR loop for Va increment.
            # Current2 now returns Jelec and Jhole in physical units (A/m^2)
            # No further scaling by Js is required.

            Jtotal = Jelec + Jhole
            fi_va[vindex,:] =fi

        for vindex in range(Total_Steps):
            Ec_result_[vindex, :] = fi_e / q - Vt * fi_va[vindex, :]  # Values from the all Node%
            Ev_result_[vindex, :] = fi_h / q - Vt * fi_va[vindex, :]  # Values from the all Node%
        ##########################################################################
        ##                 END OF NON-EQUILIBRIUM  SOLUTION PART                ##
        ##########################################################################
        # Write the results of the simulation in files #
        (
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
        ) = Write_results_non_equi2(
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
        )
        fitot = fi_h - Vt * q * fi
        fitotc = fi_e - Vt * q * fi
    time3 = time.time()  # timing audit
    # Add to log
    logger.info("calculation time  %g s", (time3 - time2))

    class Results:
        pass

    results = Results()
    results.N_wells_virtual = N_wells_virtual
    results.Well_boundary = Well_boundary
    results.xaxis = xaxis
    results.wfh = wfh
    results.wfe = wfe
    results.wfh_general = wfh_general
    results.wfe_general = wfe_general
    results.fitot = fitot
    results.fitotc = fitotc
    results.fi_e = fi_e
    results.fi_h = fi_h
    # results.sigma = sigma
    results.sigma_general = sigma_general
    # results.F = F
    results.V = V
    results.E_state = E_state
    results.N_state = N_state
    # results.meff_state = meff_state
    results.E_statec = E_statec
    results.N_statec = N_statec
    # results.meff_statec = meff_statec
    results.F_general = F_general
    results.E_state_general = E_state_general
    results.N_state_general = N_state_general
    results.meff_state_general = meff_state_general
    results.E_statec_general = E_statec_general
    results.N_statec_general = N_statec_general
    results.meff_statec_general = meff_statec_general
    results.Fapp = Fapp
    results.T = T
    # results.E_F = E_F
    results.E_F_general = E_F_general
    results.dx = dx
    results.subnumber_h = subnumber_h
    results.subnumber_e = subnumber_e
    results.Ntotal2d = Ntotal2d
    ########################
    results.Va_t = Va_t
    results.Efn_result = Efn_result
    results.Efp_result = Efp_result
    results.Ei_result = Ei_result
    results.av_curr = av_curr
    results.Ec_result = Ec_result
    results.Ev_result = Ev_result
    results.ro_result = ro_result
    results.el_field1_result = el_field1_result
    results.el_field2_result = el_field2_result
    results.nf_result = nf_result
    results.pf_result = pf_result
    results.fi_result = fi_result
    results.EF = EF
    results.Total_Steps = Total_Steps
    results.fi_va = fi_va
    results.Ec_result_ = Ec_result_
    results.Ev_result_ = Ev_result_
    return results


def Poisson_Schrodinger_DD_test_2(result, model):
    # Initialize all potential result variables to avoid UnboundLocalError
    # We must be careful to define them before any potential access
    n_max = model.n_max
    Va_t = np.zeros(1)
    Efn_result = Efp_result = Ei_result = Ec_result = Ev_result = np.zeros(n_max)
    ro_result = el_field1_result = el_field2_result = nf_result = pf_result = np.zeros(n_max)
    fi_result = np.zeros(n_max)
    av_curr = np.zeros(1)
    EF = fi_va = Ec_result_ = Ev_result_ = None
    Total_Steps = 1
    
    fi = result.fi_result
    E_state_general = result.E_state_general
    meff_state_general = result.meff_state_general
    E_statec_general = result.E_statec_general
    meff_statec_general = result.meff_statec_general
    wfh_general = result.wfh_general
    wfe_general = result.wfe_general
    n_max = model.n_max
    dx = model.dx
    HUPMAT3_reduced_list = result.HUPMAT3_reduced_list
    HUPMATC1 = result.HUPMATC1
    m_hh = result.m_hh
    m_lh = result.m_lh
    m_so = result.m_so
    Ppz_Psp = result.Ppz_Psp
    pol_surf_char = result.pol_surf_char
    """Performs a self-consistent Poisson-Schrodinger calculation of a 1d quantum well structure.
    Model is an object with the following attributes:
    fi_e - Bandstructure potential (J) (array, len n_max)
    cb_meff - conduction band effective mass (kg)(array, len n_max)
    eps - dielectric constant (including eps0) (array, len n_max)
    dop - doping distribution (m**-3) ( array, len n_max)
    Fapp - Applied field (Vm**-1)
    T - Temperature (K)
    comp_scheme - simulation scheme (currently unused)
    subnumber_e - number of subbands for look for in the conduction band
    dx - grid spacing (m)
    n_max - number of points.
    """
    fi_e = model.fi_e
    cb_meff = model.cb_meff
    eps = model.eps
    dop = model.dop
    Fapp = model.Fapp
    vmax = model.vmax
    vmin = model.vmin
    Each_Step = model.Each_Step
    surface = model.surface
    T = model.T
    comp_scheme = model.comp_scheme
    subnumber_h = model.subnumber_h
    subnumber_e = model.subnumber_e
    dx = model.dx
    n_max = model.n_max
    TAUN0 = model.TAUN0
    TAUP0 = model.TAUP0
    mun0 = model.mun0
    mup0 = model.mup0
    
    # Scaling factors for DD system (Corrected to SI units)
    # Scaling factors for DD system (Corrected to SI units)
    # dx is already in meters, n_max is number of points
    xbar = dx * n_max  # Total device length in meters
    Vbar = Vt
    # Convert mubar from cm2/Vs to m2/Vs
    mubar_raw = max(max(mun0), max(mup0)) if max(max(mun0), max(mup0)) > 1e-12 else 0.1
    mubar = mubar_raw * 1e-4  # m2/Vs
    tbar = xbar ** 2 / (mubar * Vbar)
    ns_scale = np.linalg.norm(dop, np.inf) if np.linalg.norm(dop, np.inf) > 1e15 else 1e18
    # Rbar is the normalization for generation/recombination rate [m^-3 s^-1]
    Rbar = ns_scale / tbar
    
    # Pass physical G_optical to solver (Convert cm^-3 s^-1 to m^-3 s^-1)
    G_opt_phys = float(getattr(model, 'G_optical', 0.0)) * 1e6 
    model.G_optical_scaled = G_opt_phys / Rbar # For solver normalization (dimensionless)
    print(f"DEBUG: xbar={xbar:.2e} m, tbar={tbar:.2e} s, Rbar={Rbar:.2e} m^-3/s")
    print(f"DEBUG: G_opt_phys={G_opt_phys:.2e}, G_scaled={model.G_optical_scaled:.2e}")
    
    Cn0 = model.Cn0
    Cp0 = model.Cp0
    BETAN = model.BETAN
    BETAP = model.BETAP
    VSATN = model.VSATN
    VSATP = model.VSATP

    if comp_scheme in (4, 5, 6):
        logger.error(
            """Aestimo doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        sys.exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """Aestimo doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        sys.exit()
    fi_h = model.fi_h
    N_wells_virtual = model.N_wells_virtual
    Well_boundary = model.Well_boundary
    x_max = dx * n_max
    UNIM = np.identity(n_max)
    RATIO = m_e / hbar ** 2 * (x_max) ** 2

    # Check
    if comp_scheme == 6:
        logger.warning(
            """The calculation of Vxc depends upon m*, however when non-parabolicity is also 
                 considered m* becomes energy dependent which would make Vxc energy dependent.
                 Currently this effect is ignored and Vxc uses the effective masses from the 
                 bottom of the conduction bands even when non-parabolicity is considered 
                 elsewhere."""
        )
    # Preparing empty subband energy lists.
    E_state = [0.0] * subnumber_h  # Energies of subbands/levels (meV)
    N_state = [0.0] * subnumber_h  # Number of carriers in subbands
    meff_state = [0.0] * subnumber_h # Effective mass of subbands
    E_statec = [0.0] * subnumber_e  # Energies of subbands/levels (meV)
    N_statec = [0.0] * subnumber_e  # Number of carriers in subbands
    meff_statec = [0.0] * subnumber_e # Effective mass of subbands
    # Preparing empty subband energy arrays for multiquantum wells.
    """
    E_state_general = np.zeros((model.N_wells_virtual,subnumber_h))     # Energies of subbands/levels (meV)
    E_statec_general = np.zeros((model.N_wells_virtual,subnumber_e))     # Energies of subbands/levels (meV)
    meff_statec_general= np.zeros((model.N_wells_virtual,subnumber_e))
    meff_state_general= np.zeros((model.N_wells_virtual,subnumber_h))
    """
    N_state_general = np.zeros(
        (model.N_wells_virtual, subnumber_h)
    )  # Number of carriers in subbands
    N_statec_general = np.zeros(
        (model.N_wells_virtual, subnumber_e)
    )  # Number of carriers in subbands

    # Creating and Filling material arrays
    xaxis = np.arange(0, n_max) * dx  # metres
    fitot = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potential
    fitotc = np.zeros(n_max)  # Energy potential = Bandstructure + Coulombic potentia
    # eps = np.zeros(n_max+2)	    #dielectric constant
    # dop = np.zeros(n_max+2)	    #doping distribution
    # sigma = np.zeros(n_max+2)      #charge distribution (donors + free charges)
    # F = np.zeros(n_max+2)          #Electric Field
    # Vapp = np.zeros(n_max+2)       #Applied Electric Potential
    V = np.zeros(n_max)  # Electric Potential

    # Subband wavefunction for holes list. 2-dimensional: [i][j] i:stateno, j:wavefunc

    wfh = np.zeros((subnumber_h, n_max))
    wfe = np.zeros((subnumber_e, n_max))
    """
    wfh_general = np.zeros((model.N_wells_virtual,subnumber_h,n_max))
    wfe_general = np.zeros((model.N_wells_virtual,subnumber_e,n_max))
    """
    E_F_general = np.zeros(model.N_wells_virtual)
    sigma_general = np.zeros(n_max)
    F_general = np.zeros(n_max)
    Vnew_general = np.zeros(n_max)
    # fi = np.zeros(n_max)
    # Setup the doping
    Ntotal = sum(dop)  # calculating total doping density m-3
    Ntotal2d = Ntotal * dx
    # Add to log
    logger.info("Ntotal2d %g m**-2", Ntotal2d)
    # Applied Field
    # Vapp = calc_potn(Fapp*eps0/eps,model)
    # Vapp[n_max-1] -= Vapp[n_max//2] #Offsetting the applied field's potential so that it is zero in the centre of the structure.
    # s
    # setting up Ldi and Ld p and n
    Ld_n_p = np.zeros(n_max)
    Ldi = np.zeros(n_max)
    Nc = np.zeros(n_max)
    Nv = np.zeros(n_max)
    vb_meff = np.zeros(n_max)
    ni = np.zeros(n_max)
    hbark = hbar * 2 * pi
    Ppz_Psp_tmp = Ppz_Psp
    Ppz_Psp = np.zeros(n_max)
    for i in range(n_max):
        vb_meff[i] = (m_hh[i] ** (3 / 2) + m_lh[i] ** (3 / 2) + m_so[i] ** (3 / 2)) ** (
            2 / 3
        )
    Nc = 2 * (2 * pi * cb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Nv = 2 * (2 * pi * vb_meff * kb * T / hbark ** 2) ** (3 / 2)
    Half_Eg = np.zeros(n_max)
    for i in range(n_max):
        val_ni = sqrt( Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T)) )
        ni[i] = max(val_ni, 1e18) # Unified ni_ref scaling
        # No longer using Ldi here as it is recomputed scaled
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
        # fi_e scaled
        # fi_h scaled
        # fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        # fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
    n = result.nf_result / ni
    p = result.pf_result / ni
    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )

    # STARTING SELF CONSISTENT LOOP
    time2 = time.time()  # timing audit
    iteration = 1  # iteration counter
    # previousE0= 0   #(meV) energy of zeroth state for previous iteration(for testing convergence)
    # fitot = fi_h + Vapp #For initial iteration sum bandstructure and applied field
    # fitotc = fi_e + Vapp
    xaxis = np.arange(0, n_max) * dx  # metres
    mup = np.zeros(n_max)
    mun = np.zeros(n_max)
    n_q = np.zeros(n_max)
    p_q = np.zeros(n_max)
    fi_n = np.zeros(n_max)
    fi_p = np.zeros(n_max)
    EF = 0.0
    Total_Steps = int((vmax - vmin) / Each_Step) + 1
    vindex = 0
    Va_t = np.zeros(Total_Steps)

    Jtotal = np.zeros((Total_Steps, n_max))
    J_Tunnling= np.zeros((Total_Steps, n_max))
    ###############################################################
    len_ = xaxis[n_max - 1]  # xaxis is already in meters

    #
    xm = np.mean(xaxis)
    nis = np.zeros(n_max)
    mun = np.zeros(n_max)
    mup = np.zeros(n_max)
    Cn = np.zeros(n_max)
    Cp = np.zeros(n_max)
    l2 = np.zeros(n_max - 1)
    Fn = np.zeros(n_max)
    Fp = np.zeros(n_max)
    n = np.zeros(n_max)
    p = np.zeros(n_max)
    nn = np.zeros(n_max)
    pp = np.zeros(n_max)
    V = np.zeros(n_max)
    vvect = np.zeros(Total_Steps)
    n_ = np.zeros((Total_Steps, n_max))
    p_ = np.zeros((Total_Steps, n_max))
    Fn_ = np.zeros((Total_Steps, n_max))
    Fp_ = np.zeros((Total_Steps, n_max))
    V_ = np.zeros((Total_Steps, n_max))
    Jn = np.zeros((Total_Steps, n_max))
    Jp = np.zeros((Total_Steps, n_max))
    fi_va = np.zeros((Total_Steps, n_max))
    Ec_result_ = np.zeros((Total_Steps, n_max))
    Ev_result_ = np.zeros((Total_Steps, n_max))
    nf_result_ = np.zeros((Total_Steps, n_max))
    pf_result_ = np.zeros((Total_Steps, n_max))
    el_field1_result_ = np.zeros((Total_Steps, n_max))
    el_field2_result_ = np.zeros((Total_Steps, n_max))
    ro_result_ = np.zeros((Total_Steps, n_max))
    # J=np.zeros((Total_Steps,n_max-1))
    DV = np.zeros(Total_Steps)
    Emax = np.zeros(Total_Steps)


    nn, pp, fi_out = equi_np_fi(iteration, dop, Ppz_Psp, n_max, ni, model, Vt, surface)

    # xn = xm+1e-7
    # xp = xm-1e-7
    ## Scaling coefficients
    ## SI Scaling coefficients (Unified)
    xs = len_ # [m]
    ns = np.linalg.norm(dop, np.inf) if np.linalg.norm(dop, np.inf) > 1e15 else 1e24 # [m-3]
    Vs = Vt # [V]
        
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
            self.Ppz_Psp = Ppz_Psp
            self.Cn = Cn
            self.Cp = Cp
            self.E_state_general = result.E_state_general
            self.meff_state_general = result.meff_state_general
            self.E_statec_general = result.E_statec_general
            self.meff_statec_general = result.meff_statec_general
            self.wfh_general = result.wfh_general
            self.wfe_general = result.wfe_general

    idata = data()
    odata = data()
    idata.n = nn * ni / ns # Both in m-3
    idata.p = pp * ni / ns
    idata.V = fi_out
    
    Vs = Vt
    us_raw = max(max(mun0), max(mup0)) if max(max(mun0), max(mup0)) > 1e-12 else 0.1
    us = us_raw * 1e-4  # m2/Vs
    xbar = len_  # [m]
    Vbar = Vt  # [V]
    mubar = us  # [m^2 V^{-1} s^{-1}]
    tbar = xbar ** 2 / (mubar * Vbar)  # [s]
    # Rbar is the normalization for generation/recombination rate [m^-3 s^-1]
    # ns is m-3, tbar is s
    Rbar = ns / tbar
    # CAubar is Auger normalization [m^6 s^-1]
    # R_phys = Cn * n^3 => R_norm = (Cn * ns^2 / (1/tbar)) * n_norm^3
    CAubar = Rbar / ns ** 2  
    idata.Cn = Cn0 / CAubar
    idata.Cp = Cp0 / CAubar
    # Using unified SI G_optical (m^-3 s^-1)
    idata.G_optical = model.G_optical_scaled
    ###############################################################
    if Total_Steps < 2:
        print("Equilibrium only (Total_Steps < 2)")
    else:
        print("Convergence of the Gummel cycles")
        vindex = 0
        for vindex in range(0, Total_Steps):
            # introducing piezo spont effect with increasing ratio till 33.33%
            Ppz_Psp = Ppz_Psp_tmp / (Total_Steps + 2 - vindex)
            # piezo_ratio=100*np.linalg.norm(Ppz_Psp,np.inf)/np.linalg.norm(Ppz_Psp_tmp,np.inf)
            # print("ratio of piezo=%.2f"%piezo_ratio," %")
            # Start Va increment loop
            Va = vmin
            Va += Each_Step * vindex
            Va_t[vindex] = Va

            print("Va_t[", vindex, "]=%.2f" % Va_t[vindex])
            #####################################################################################
            vvect[vindex] = Va
            # z
            xin = xaxis / xs
            
            # Seed each step with equilibrium solution (normalized)
            n_[vindex, :] = nn # Normalized to ns
            p_[vindex, :] = pp # Normalized to ns
            # Non-equilibrium seeds (normalized to Vs = Vt)
            V_app_norm = (Va / Vs) * (xaxis <= xm)
            V_[vindex, :] = fi_out + V_app_norm
            Fn_[vindex, :] = (V_app_norm) - np.log(ni / ns)
            Fp_[vindex, :] = (V_app_norm) + np.log(ni / ns)

            idata.l2 = (Vs * eps[0 : n_max - 1]) / (q * ns * xs ** 2)
            idata.nis = ni / ns
            idata.dop = dop / ns
            idata.Ppz_Psp = Ppz_Psp / ns
            # mun,mup=Mobility2(mun0,mup0,fi,Vt,Ldi,VSATN,VSATP,BETAN,BETAP,n_max,dx)
            from aeslibs.func_lib import CaugheyThomasMobility
            mun_ct, mup_ct = CaugheyThomasMobility(n_[vindex, :], p_[vindex, :], material_type='Si')
            idata.mun = mun_ct / us
            idata.mup = mup_ct / us

            # sinodes = np.arange(len(xaxis))
            idata.TAUN0 = TAUN0 / tbar  # np.inf
            idata.TAUP0 = TAUP0 / tbar  # np.inf
            idata.theta = ni / ns

            idata.n = n_[vindex, :] # Already normalized to ns
            idata.p = p_[vindex, :] # Already normalized to ns
            idata.V = V_[vindex, :] / Vs
            idata.Fn = Fn_[vindex, :] / Vs
            idata.Fp = Fp_[vindex, :] / Vs
            idata.V_applied = Va # For selective contact BCs
            fitot = fi_h - Vt * q * idata.V
            fitotc = fi_e - Vt * q * idata.V
            if model.N_wells_virtual - 2 != 0 and config.quantum_effect:
                (
                    idata.E_statec_general,
                    idata.E_state_general,
                    idata.wfe_general,
                    idata.wfh_general,
                    idata.meff_statec_general,
                    idata.meff_state_general,
                ) = Schro(
                    HUPMAT3_reduced_list,
                    HUPMATC1,
                    subnumber_h,
                    subnumber_e,
                    fitot,
                    fitotc,
                    model,
                    Well_boundary,
                    UNIM,
                    RATIO,
                    m_hh,
                    m_lh,
                    m_so,
                    n_max,
                )
            ## Solution of DD system
            #
            ## Algorithm parameters
            toll = 1e-6  # Gummel convergence tolerance
            maxit = 50   # Max Gummel iterations
            ptoll = 1e-10  # Poisson solver tolerance
            pmaxit = 20   # Poisson solver max iterations
            verbose = 0   # Quiet mode for performance
               
            [odata, it, res] = DDGgummelmap(
                n_max,
                xin,
                idata,
                odata,
                toll,
                maxit,
                ptoll,
                pmaxit,
                verbose,
                ni,
                fi_e,
                fi_h,
                model,
                Vt,
            )
            if getattr(model, 'photovoltaic_mode', False):
                 # Sync Fermi levels with PV boundary conditions after solver
                 # This ensures Fn/Fp at contacts are consistent with selective contacts
                 fermin = np.vstack((odata.V, odata.Fn)).T
                 fermip = np.vstack((odata.V, odata.Fp)).T
                 fermin, fermip = apply_photovoltaic_BCs(fermin, fermip, odata.V, n_max, model, Va, idata)
                 odata.Fn = fermin[:, 1]
                 odata.Fp = fermip[:, 1]
                 
                 # Ensure positive carrier densities
                 odata.n = np.maximum(np.exp(odata.V - odata.Fn), 1e-20)
                 odata.p = np.maximum(np.exp(odata.Fp - odata.V), 1e-20)
            
            
            # --- Series Resistance (Rs) Iterative Solver ---
            from aeslibs.func_lib import Ubernoulli
            
            Rs_val = getattr(model, 'Rs', 0.0)
            initial_V_bc = idata.V[n_max-1] 
            rs_converged = False
            
            Device_Area = getattr(model, 'device_area_m2', 1e-8) 
            rs_iters = 200 if Rs_val > 1e-6 else 1
            
            newton_toll = 1e-5
            newton_maxit = 100
            
            # --- Best State Tracker ---
            best_diff = 1e20
            best_odata = None
            damp = 0.1 # Reduced damping for high-bias stability
            
            # Scale before Newton solve
            idata.n = n_[vindex, :]
            idata.p = p_[vindex, :]
            idata.V = V_[vindex, :]
            
            for rs_it in range(rs_iters):
                [odata, it, res] = DDNnewtonmap(
                    ni, fi_e, fi_h, xin, odata, newton_toll, newton_maxit, verbose, model, Vs
                )
                
                if Rs_val <= 1e-6:
                    best_odata = odata
                    rs_converged = True
                    break
                
                # Manual current calculation for Rs adjustment
                v_curr = odata.V
                n_curr = odata.n
                p_curr = odata.p
                arg = -(v_curr[1:] - v_curr[:-1])
                Bp_vec = Ubernoulli(arg, 1)
                Bm_vec = Ubernoulli(arg, 0)
                dx_vec = xin[1:] - xin[:-1]
                
                Jn_scaled = -odata.mun[0:n_max-1] * (n_curr[1:]*Bp_vec - n_curr[:-1]*Bm_vec) / dx_vec
                
                arg_p = v_curr[1:] - v_curr[:-1]
                Jp_scaled = odata.mup[0:n_max-1] * (p_curr[1:]*Ubernoulli(arg_p, 0) - p_curr[:-1]*Ubernoulli(arg_p, 1)) / dx_vec
                
                # Aestimo library Bernoulli sign convention check - ensure consistency
                J_tot_scaled = np.abs(Jn_scaled + Jp_scaled)
                J_tot_val_scaled = np.median(J_tot_scaled)
                
                J_physical = J_tot_val_scaled * (us * q * ns * Vs / xs)
                I_physical = J_physical * Device_Area
                
                V_drop_scaled = (I_physical * Rs_val) / Vs
                new_V_bc = initial_V_bc - V_drop_scaled
                
                diff = abs(new_V_bc - odata.V[n_max-1])
                
                # Track best state
                if diff < best_diff:
                    best_diff = diff
                    import copy
                    best_odata = copy.deepcopy(odata)
                
                if diff < 1e-4:
                    rs_converged = True
                    break
                    
                # Damping with clipping
                max_step = 1.0 # 25mV limit
                delta_V = np.clip(new_V_bc - odata.V[n_max-1], -max_step, max_step)
                odata.V[n_max-1] += damp * delta_V
                
                # Sync for next solve
                idata.V = odata.V.copy()
            
            # Use the best state found during iterations
            odata = best_odata
            # ---------------------------------------------

            n_[vindex, :] = odata.n
            p_[vindex, :] = odata.p
            V_[vindex, :] = odata.V
            fi_va[vindex, :] = odata.V
            # print("n_newt=",odata.n[:])
            Fn_[vindex, :] = odata.Fn
            Fp_[vindex, :] = odata.Fp
            DV[vindex] = V_[vindex, n_max - 1] - V_[0, vindex]
            Emax[vindex] = max(
                abs(
                    (V_[vindex, 1:n_max] - V_[vindex, 0 : n_max - 1])
                    / (xin[1:n_max] - xin[0 : n_max - 1])
                )
            )
            #

            # Band offsets for Bernoulli current (Normalized by Vs)
            # fi_e, fi_h are in J. Convert to normalized potential.
            fi_n_norm = -fi_e / (kb * T)
            fi_p_norm = -fi_h / (kb * T)

            Bp = Ubernoulli(
                (V_[vindex, 1:n_max] - V_[vindex, 0 : n_max - 1])
                + (fi_n_norm[1:n_max] - fi_n_norm[0 : n_max - 1]),
                1,
            )
            Bm = Ubernoulli(
                (V_[vindex, 1:n_max] - V_[vindex, 0 : n_max - 1])
                + (fi_p_norm[1:n_max] - fi_p_norm[0 : n_max - 1]),
                0,
            )

            # Use SI units for current density calculation
            # xin is in m, mun is cm2/Vs, n_ is normalized by ns
            # Jn = (mun * 1e-4) * q * (ns * n_) * Vt * (Bp - Bm) / dx
            # However, aestimo uses a slightly different normalized form.
            # We standardize to SI: J = q * mu * n * E + q * D * grad(n)
            
            # Physical Current Density Calculation (SI Pure [A/m^2])
            # Formula: J = (q * mu * Vt * ni / dx) * [n_norm_{i+1} * B(dv) - n_norm_i * B(-dv)]
            # We use physical mobility [m2/Vs] and physical density n_phys = n_norm * ni
            
            # Local dx [m] and Vt [V]
            dx_eff = (xin[1:n_max] - xin[0 : n_max - 1]) * xs
            
            # Local mobilities converted to m2/Vs
            mun_m2 = odata.mun[0 : n_max - 1] * us
            mup_m2 = odata.mup[0 : n_max - 1] * us
            
            # Current components with Scharfetter-Gummel discretization
            # n_norm here is n_phys / ni (from equi_np_fi or Solver)
            Jn[vindex, 0 : n_max - 1] = (
                (q * mun_m2 * Vt * ns / dx_eff) 
                * (n_[vindex, 1:n_max] * Bp - n_[vindex, 0 : n_max - 1] * Bm)
            )
            Jp[vindex, 0 : n_max - 1] = (
                (q * mup_m2 * Vt * ns / dx_eff) 
                * (p_[vindex, 0 : n_max - 1] * Bp - p_[vindex, 1:n_max] * Bm)
            )
            
            # No early stopping - full sweep required for accurate Voc extraction
            
        ## Descaling to physical SI units
        # Restore carrier and potential scaling for GUI displays
        # Potential is normalized to Vt, densities to ns (m^-3)
        n_ = n_ * ns
        p_ = p_ * ns
        V_ = V_ * Vs
        Fn_ = Fn_ * Vs
        Fp_ = Fp_ * Vs

        # Jtotal is in A/m^2 (SI). Convert to mA/cm^2 for Aestimo GUI (1 A/m2 = 0.1 mA/cm2)
        Jtotal = (Jp + Jn) * 0.1
        #Fn = V_ / Vs - np.log(n_)
        #Fp = V_ / Vs + np.log(p_)
        # Fn_=Fn_*Vs
        # Fp_=Fp_*Vs
        #
        time1 = time.time()
        delta_t = (time1 - time0) / 60
        print("time=%.2fmn" % delta_t)

        ro_result = np.zeros(n_max)
        el_field1_result = np.zeros(n_max)
        el_field2_result = np.zeros(n_max)
        Ec_result = np.zeros(n_max)
        Ev_result = np.zeros(n_max)
        Ei_result = np.zeros(n_max)
        Efn_result = np.zeros(n_max)
        Efp_result = np.zeros(n_max)
        av_curr = np.zeros(Total_Steps)
        fi_result = V_[vindex, :]
        # Efn_result,Efp_result=Fn_[vindex,:],Fp_[vindex,:]
        nf_result, pf_result = n_[vindex, :], p_[vindex, :]
        # Use median of p-side (0-20% of device) with sign negated for photovoltaic convention.
        # In the Newton-Krylov path, Jtotal in the n-side (80-100%) is positive and INCREASES with
        # forward bias because the dark current adds in the same direction as photocurrent.
        # The p-side (0-20%) has the correct sign: photocurrent is negative (flows right-to-left
        # in the n→p conventional direction). Negating gives the standard convention:
        #   av_curr < 0 at V=0 (= -Jsc), rises toward 0 at Voc, positive for V > Voc.
        idx_lo = 1
        idx_hi = max(2, int(0.2 * n_max))
        for k in range(Total_Steps):
            av_curr[k] = -np.median(Jtotal[k, idx_lo:idx_hi])
        for i in range(1, n_max - 1):
            Ec_result[i] = fi_e[i] / q - V_[vindex, i]  # Values from the second Node%
            Ev_result[i] = fi_h[i] / q - V_[vindex, i]  # Values from the second Node%
            Ei_result[i] = Ec_result[i] - ((fi_e[i] - fi_h[i]) / (2 * q))
            ro_result[i] = -q * (n_[vindex, i] - p_[vindex, i] - dop[i])
            el_field1_result[i] = -(V_[vindex, i + 1] - V_[vindex, i]) / (dx)
            el_field2_result[i] = -(V_[vindex, i + 1] - V_[vindex, i - 1]) / (2 * dx)
            Efn_result[i] = Ei_result[i] + Vt * np.log(np.maximum(n_[vindex, i]/ni[i], 1e-20))
            Efp_result[i] = Ei_result[i] - Vt * np.log(np.maximum(p_[vindex, i]/ni[i], 1e-20))
        Ec_result[0] = Ec_result[1]
        Ec_result[n_max - 1] = Ec_result[n_max - 2]
        Ev_result[0] = Ev_result[1]
        Ev_result[n_max - 1] = Ev_result[n_max - 2]

        Ei_result[0] = Ei_result[1]
        Ei_result[n_max - 1] = Ei_result[n_max - 2]
        Efn_result[0] = Efn_result[1]
        Efn_result[n_max - 1] = Efn_result[n_max - 2]

        Efp_result[0] = Efp_result[1]
        Efp_result[n_max - 1] = Efp_result[n_max - 2]
        el_field1_result[0] = el_field1_result[1]
        el_field2_result[0] = el_field2_result[1]
        el_field1_result[n_max - 1] = el_field1_result[n_max - 2]
        el_field2_result[n_max - 1] = el_field2_result[n_max - 2]
        ro_result[0] = ro_result[1]
        ro_result[n_max - 1] = ro_result[n_max - 2]
        nf_result[0] = nf_result[1]
        nf_result[n_max - 1] = nf_result[n_max - 2]
        pf_result[0] = pf_result[1]
        pf_result[n_max - 1] = pf_result[n_max - 2]
        Va_t = vvect
        fitot = fi_h - Vt * q * odata.V
        fitotc = fi_e - Vt * q * odata.V
        for vindex in range(Total_Steps):
            Ec_result_[vindex, :] = fi_e / q - V_[vindex, :]  # Values from the all Node%
            Ev_result_[vindex, :] = fi_h / q - V_[vindex, :]  # Values from the all Node%
            nf_result_[vindex, :] = n_[vindex, :]
            pf_result_[vindex, :] = p_[vindex, :]
            ro_result_[vindex, 1:n_max-1] = -q * (n_[vindex, 1:n_max-1] - p_[vindex, 1:n_max-1] - dop[1:n_max-1])
            ro_result_[vindex, 0] = ro_result_[vindex, 1]
            ro_result_[vindex, n_max-1] = ro_result_[vindex, n_max-2]
            el_field1_result_[vindex, 1:n_max-1] = -(V_[vindex, 2:n_max] - V_[vindex, 1:n_max-1]) / (dx)
            el_field1_result_[vindex, 0] = el_field1_result_[vindex, 1]
            el_field1_result_[vindex, n_max-1] = el_field1_result_[vindex, n_max-2]
            el_field2_result_[vindex, 1:n_max-1] = -(V_[vindex, 2:n_max] - V_[vindex, 0:n_max-2]) / (2 * dx)
            el_field2_result_[vindex, 0] = el_field2_result_[vindex, 1]
            el_field2_result_[vindex, n_max-1] = el_field2_result_[vindex, n_max-2]
        if model.N_wells_virtual - 2 != 0 and config.quantum_effect:
            (
                idata.E_statec_general,
                idata.E_state_general,
                idata.wfe_general,
                idata.wfh_general,
                idata.meff_statec_general,
                idata.meff_state_general,
            ) = Schro(
                HUPMAT3_reduced_list,
                HUPMATC1,
                subnumber_h,
                subnumber_e,
                fitot,
                fitotc,
                model,
                Well_boundary,
                UNIM,
                RATIO,
                m_hh,
                m_lh,
                m_so,
                n_max,
            )
    time3 = time.time()  # timing audit
    # Add to log
    logger.info("calculation time  %g s", (time3 - time2))

    class Results:
        pass

    results = Results()
    results.N_wells_virtual = N_wells_virtual
    results.Well_boundary = Well_boundary
    results.xaxis = xaxis
    results.wfh = wfh
    results.wfe = wfe
    results.wfh_general = idata.wfh_general
    results.wfe_general = idata.wfe_general
    results.fitot = fitot
    results.fitotc = fitotc
    results.fi_e = fi_e
    results.fi_h = fi_h
    # results.sigma = sigma
    results.sigma_general = sigma_general
    # results.F = F
    results.V = V
    results.E_state = E_state
    results.N_state = N_state
    results.meff_state = meff_state
    results.E_statec = E_statec
    results.N_statec = N_statec
    results.meff_statec = meff_statec
    results.F_general = F_general
    results.E_state_general = idata.E_state_general
    results.N_state_general = N_state_general
    results.meff_state_general = idata.meff_state_general
    results.E_statec_general = idata.E_statec_general
    results.N_statec_general = N_statec_general
    results.meff_statec_general = idata.meff_statec_general
    results.Fapp = Fapp
    results.T = T
    # results.E_F = E_F
    results.E_F_general = E_F_general
    results.dx = dx
    results.subnumber_h = subnumber_h
    results.subnumber_e = subnumber_e
    results.Ntotal2d = Ntotal2d
    ########################
    results.Va_t = Va_t
    results.Efn_result = Efn_result
    results.Efp_result = Efp_result
    results.Ei_result = Ei_result
    results.av_curr = av_curr
    results.Ec_result = Ec_result
    results.Ev_result = Ev_result
    results.ro_result = ro_result
    results.el_field1_result = el_field1_result
    results.el_field2_result = el_field2_result
    results.nf_result = nf_result
    results.pf_result = pf_result
    results.fi_result = fi_result
    results.EF = EF
    results.Total_Steps = Total_Steps
    results.fi_va = fi_va
    results.Ec_result_ = Ec_result_
    results.Ev_result_ = Ev_result_
    results.nf_result_ = nf_result_
    results.pf_result_ = pf_result_
    results.el_field1_result_ = el_field1_result_
    results.el_field2_result_ = el_field2_result_
    results.ro_result_ = ro_result_
    return results





def run_aestimo(input_obj, drawFigures=drawFigures, show=True):
    """A utility function that performs the standard simulation run
    for 'normal' input files. Input_obj can be a dict, class, named tuple or 
    module with the attributes needed to create the StructureFrom class, see 
    the class implementation or some of the sample-*.py files for details."""
    
    global output_directory
    # Hack: If output_directory is currently generic 'output', and we know the input filename, switch it.
    # This supports running examples directly like `python examples/sample.py` which import aestimo.
    if os.path.basename(output_directory) == 'output':
        # Try to find meaningful name
        name = None
        if hasattr(input_obj, '__file__'):
            name = Path(input_obj.__file__).stem
        elif isinstance(input_obj, dict) and '__file__' in input_obj:
            name = Path(input_obj['__file__']).stem
        elif hasattr(input_obj, 'inputfilename'):
            name = input_obj.inputfilename
        
        if name:
             # Repoint output directory
             new_out = os.path.join(os.getcwd(), name + "_output")
             if not os.path.isdir(new_out):
                 os.makedirs(new_out, exist_ok=True)
             output_directory = new_out

    # Add to log
    # Note: If we changed output_directory, the logger is still pointing to the old file 
    # if it was already initialized. However, usually initialize_logger is called before this.
    # If we want logs in the new directory, we'd need to re-init logger. 
    # But initialize_logger uses the global output_directory.
    # For now, we accept logs might be in 'output' or we should technically re-init logger here.
    # Let's leave logger as is to avoid complex side effects, as users mainly care about data results.
    
    # Initialise structure class
    model = StructureFrom(input_obj, database)

    # Perform the calculation
    
    print(f"DEBUG: run_aestimo called. model.comp_scheme={model.comp_scheme}")
    if model.comp_scheme == 11:
        # Scheme 11: 8-band k·p with arbitrary crystal orientation
        result = Poisson_Schrodinger_new(model)
    else:
        result = Poisson_Schrodinger(model)
    
    print(f"DEBUG: Poisson_Schrodinger done. checking comp_scheme for DD routing: {model.comp_scheme}")
    if model.comp_scheme in (7, 10):
        if model.comp_scheme == 10:
            model.use_newton_solver = True
        print("DEBUG: Routing to Poisson_Schrodinger_DD (Fully-Coupled Newton-Raphson enabled for scheme 10)")
        result_dd = Poisson_Schrodinger_DD(result, model)
    if model.comp_scheme == 8:
        result_dd = Poisson_Schrodinger_DD_test(result, model)
    if model.comp_scheme == 9:
        result_dd = Poisson_Schrodinger_DD_test_2(result, model)
    time4 = time.time()  # timing audit
    # Add to log
    logger.info("total running time (inc. loading libraries) %g s", (time4 - time0))
    logger.info("total running time (exc. loading libraries) %g s", (time4 - time1))
    # Write the simulation results in files

    figs_out = []
    if model.comp_scheme in (0, 1, 2, 7, 8, 10):
        res_figs = save_and_plot(result, model, output_directory, drawFigures=drawFigures, show=show)
        if isinstance(res_figs, list): figs_out.extend(res_figs)
    if model.comp_scheme in (7,8,9,10):
        res_figs2 = save_and_plot2(result_dd, model, output_directory, drawFigures=drawFigures, show=show)
        if isinstance(res_figs2, list): figs_out.extend(res_figs2)
    figures = figs_out
    
    # Experimental Validation Hook
    # Check if input_obj or model has experimental validation enabled
    enable_val = False
    if isinstance(input_obj, dict):
        enable_val = input_obj.get('enable_experimental_validation', False)
    else:
        enable_val = getattr(input_obj, 'enable_experimental_validation', False)

    # If DD was performed, return that result as it contains more info
    final_res = result
    if model.comp_scheme in (7, 8, 9, 10) and 'result_dd' in locals():
        final_res = result_dd

    # Quantum-Well Confined States Solver Hook (Advanced Physics Mode)
    enable_qw = getattr(model, 'enable_qw_solver', False)
    if isinstance(input_obj, dict):
        enable_qw = enable_qw or input_obj.get('enable_qw_solver', False)
    else:
        enable_qw = enable_qw or getattr(input_obj, 'enable_qw_solver', False)

    if enable_qw:
        try:
            from aeslibs.quantum_well import solve_quantum_well, solve_self_consistent_qw_poisson
            logger.info("Running Advanced Quantum-Well Confined-State Solver module...")
            
            n_pts = model.n_max
            dx_nm = model.dx * 1e9
            z_grid = np.arange(n_pts) * dx_nm

            # Get Ec profile in eV
            if hasattr(final_res, 'fitotc') and final_res.fitotc is not None and len(final_res.fitotc) == n_pts:
                ec_raw = np.asarray(final_res.fitotc, dtype=float)
            elif hasattr(final_res, 'Ec_result') and final_res.Ec_result is not None and len(final_res.Ec_result) == n_pts and np.any(final_res.Ec_result != 0):
                ec_raw = np.asarray(final_res.Ec_result, dtype=float)
            elif hasattr(model, 'fi_e') and len(model.fi_e) == n_pts:
                ec_raw = np.asarray(model.fi_e, dtype=float)
            else:
                ec_raw = np.zeros(n_pts)

            if np.max(np.abs(ec_raw)) < 1e-10:
                ec_ev = ec_raw / 1.602176634e-19
            else:
                ec_ev = ec_raw

            # Get Ev profile in eV
            if hasattr(final_res, 'fitot') and final_res.fitot is not None and len(final_res.fitot) == n_pts:
                ev_raw = np.asarray(final_res.fitot, dtype=float)
            elif hasattr(final_res, 'Ev_result') and final_res.Ev_result is not None and len(final_res.Ev_result) == n_pts and np.any(final_res.Ev_result != 0):
                ev_raw = np.asarray(final_res.Ev_result, dtype=float)
            elif hasattr(model, 'fi_h') and len(model.fi_h) == n_pts:
                ev_raw = np.asarray(model.fi_h, dtype=float)
            else:
                ev_raw = ec_raw - 1.424 * 1.602176634e-19

            if np.max(np.abs(ev_raw)) < 1e-10:
                ev_ev = ev_raw / 1.602176634e-19
            else:
                ev_ev = ev_raw

            layer_dicts = []
            for l in model.material:
                th = float(l[0])
                mat = str(l[1])
                x = float(l[2]) if len(l) > 2 else 0.0
                y = float(l[3]) if len(l) > 3 else 0.0
                dop = float(l[4]) if len(l) > 4 else 0.0
                dtype = str(l[5]) if len(l) > 5 else "n"
                ltype = "well" if (len(l) > 6 and str(l[6]).lower() == 'w') else "barrier"
                layer_dicts.append({
                    "thickness": th,
                    "material": mat,
                    "mole": x,
                    "mole_y": y,
                    "doping": dop,
                    "doping_type": dtype,
                    "type": ltype
                })

            num_e = getattr(model, 'num_electron_states', 3)
            num_h = getattr(model, 'num_hole_states', 3)
            self_consistent = getattr(model, 'qw_self_consistent', False)

            if self_consistent:
                dop_arr = model.dop * 1e-6 if hasattr(model, 'dop') else np.zeros(n_pts)
                eps_arr = model.eps / 8.8541878128e-12 if hasattr(model, 'eps') else np.full(n_pts, 12.9)
                qw_result = solve_self_consistent_qw_poisson(
                    z_nm=z_grid,
                    initial_ec=ec_ev,
                    initial_ev=ev_ev,
                    dielectric_rel=eps_arr,
                    doping_profile_cm3=dop_arr,
                    layers=layer_dicts,
                    temperature_k=getattr(model, 'T', 300.0),
                    max_iterations=getattr(model, 'qw_max_iterations', 20),
                    tolerance_ev=getattr(model, 'qw_tolerance', 1e-4),
                    damping_factor=getattr(model, 'qw_damping', 0.2),
                    num_e_states=num_e,
                    num_h_states=num_h
                )
            else:
                qw_result = solve_quantum_well(
                    band_profile={"z": z_grid, "ec": ec_ev, "ev": ev_ev},
                    layers=layer_dicts,
                    temperature_k=getattr(model, 'T', 300.0),
                    num_electron_states=num_e,
                    num_hole_states=num_h,
                    coupling_mode=getattr(model, 'qw_coupling_mode', "Coupled MQW"),
                    mat_system=getattr(model, 'mat_type', 'Zincblende')
                )

            final_res.qw_result = qw_result
            logger.info("QW Solver: Found %d electron states and %d hole states.", len(qw_result.electron_energies), len(qw_result.hole_energies))
            if qw_result.dominant_transitions:
                top_trans = qw_result.dominant_transitions[0]
                logger.info("QW Ground Optical Transition: %s | E = %.4f eV | lambda = %.1f nm | Overlap Gamma = %.3f",
                            top_trans['name'], top_trans['energy_ev'], top_trans['wavelength_nm'], top_trans['overlap'])
        except Exception as qw_err:
            logger.warning("Quantum-Well Confined-State solver encountered an error: %s", qw_err)

    # Add to log
    logger.info("Simulation is finished. All files are closed. Please control the related files.")
    return input_obj, model, final_res, figures

def calculate_eqe(model, wavelengths):
    """
    Calculates the External Quantum Efficiency (EQE) for a given set of wavelengths.
    This performs a spectral sweep, calculating current response for each wavelength.
    """
    eqe_results = []
    logger.info("Starting EQE Calculation for %d wavelengths", len(wavelengths))
    
    # Save original generation state
    orig_type = getattr(model, 'generation_type', 'Uniform')
    orig_alpha = getattr(model, 'alpha', 0.0)
    orig_g = getattr(model, 'G_optical', 0.0)
    
    # Photon flux for EQE (typically 1e17 cm^-2 s^-1 or similar low injection)
    photon_flux = 1e17 # cm^-2 s^-1
    
    for wl in wavelengths:
        logger.info(f"  Calculating EQE at {wl} nm...")
        # Update model for this wavelength
        # In a real scenario, alpha(lambda) would come from a database.
        # Here we use a simple placeholder or the user-provided alpha.
        model.generation_type = 'Exponential (Beer-Lambert)'
        # For demo, use current alpha or a simple 1/wl scaling if alpha not provided
        model.G_optical = photon_flux 
        
        # Run one-point simulation (usually at V=0 for EQE/IQE)
        orig_vmin, orig_vmax, orig_step = model.vmin, model.vmax, model.Each_Step
        model.vmin, model.vmax, model.Each_Step = 0.0, 0.0, 0.1
        
        _, _, res, _ = run_aestimo(model, drawFigures=False)
        
        # Restore voltages
        model.vmin, model.vmax, model.Each_Step = orig_vmin, orig_vmax, orig_step
        
        # Calculate EQE
        # EQE = (Jsc / q) / PhotonFlux
        # Jsc is in A/m^2. q = 1.6e-19 C. 
        # Jsc / q -> electrons/m^2/s. 
        # PhotonFlux in cm^-2s^-1 -> multiply by 1e4 for m^-2s^-1
        jsc = abs(res.av_curr[0])
        flux_m2 = photon_flux * 1e4
        eqe = (jsc / q) / flux_m2 if flux_m2 > 0 else 0.0
        eqe_results.append(eqe)
        
    # Restore model
    model.generation_type = orig_type
    model.alpha = orig_alpha
    model.G_optical = orig_g
    
    return np.array(eqe_results)

if __name__ == "__main__":
    # Arguments parsing
    parser = ArgumentParser(prog ='aestimo.py', description=Description, formatter_class=RawFormatter)

    parser.add_argument("-i", "--input", dest = "inputfile", help="Input filename (will open output directory with same name)")
    parser.add_argument("-v", "--version", dest="version", action='store_true')
    parser.add_argument("-d", "--drawfigures", dest="drawfigures", action='store_true', help="Draws all data at the end of calculation.")

    args = None

    args = parser.parse_args()
    
    #Exit if no argument is provided
    if args is None:
        print("Please provide at least -i argument")  
        sys.exit()

    try:
        if args.version == True:
            import requests
            try:
                response = requests.get("https://api.github.com/repos/aestimosolver/aestimo/releases/latest", timeout=5)
                print('-------------------------------------------------------------------------------------------------------')
                print('\033[91mAestimo\033[0m 1D Version '+str(__version__))
                print('-------------------------------------------------------------------------------------------------------')
                print('The latest STABLE release was '+response.json()["tag_name"]+', which is published at '+response.json()["published_at"])
                print('Download the latest STABLE tarball release at: '+response.json()["tarball_url"])
                print('Download the latest STABLE zipball release at: '+response.json()["zipball_url"])
                print('Download the latest DEV zipball release at: https://github.com/aestimosolver/aestimo/archive/refs/heads/master.zip')
            except (requests.ConnectionError, requests.Timeout) as exception:
                print('-------------------------------------------------------------------------------------------------------')
                print('\033[91Aestimo\033[0m 1D Version '+str(__version__))
                print('-------------------------------------------------------------------------------------------------------')
                print('No internet connection available.')
            sys.exit()

        if args.inputfile is not None:
            # Get the folder information and use
            inputFile = os.path.abspath(args.inputfile)
            input_dir = os.path.dirname(inputFile)
            sys.path.append(input_dir)
            
            # Load input file
            module_name = Path(inputFile).stem
            if inputFile.endswith('.json'):
                import json
                class InputObject:
                    def __init__(self, data):
                        for key, value in data.items():
                            setattr(self, key, value)
                        
                        # Normalize JSON data to match what the core expects
                        # 1. Map 'layers' to 'material' list format
                        if hasattr(self, 'layers'):
                            self.material = []
                            for layer in self.layers:
                                # Solver expects: [thickness, mat, mole, mole_y, doping, doping_type, layer_type]
                                self.material.append([
                                    float(layer.get('thickness', 0)),
                                    layer.get('material', ''),
                                    float(layer.get('mole', 0)),
                                    float(layer.get('mole_y', 0)),
                                    float(layer.get('doping', 0)),
                                    layer.get('doping_type', 'n'),
                                    layer.get('type', 'barrier')[:1] # 'b' or 'w'
                                ])
                        
                        # Handle generation type case/string mapping
                        if hasattr(self, 'generation_type'):
                            gt = str(self.generation_type).lower()
                            if 'exponential' in gt:
                                self.generation_type = 'exponential'
                            else:
                                self.generation_type = 'uniform'
                        
                        # 2. Map other common keys if they differ
                        if hasattr(self, 'temp'): self.T = float(self.temp)
                        if hasattr(self, 'grid_step'): self.gridfactor = float(self.grid_step)
                        if hasattr(self, 'max_pts'): self.maxgridpoints = int(self.max_pts)
                        if hasattr(self, 'mat_sys'): self.mat_type = self.mat_sys
                        if hasattr(self, 'vstep'): self.Each_Step = float(self.vstep)
                        if hasattr(self, 'area'): self.device_area = float(self.area)
                        if hasattr(self, 'tat_field'): self.tat_field = float(self.tat_field)
                        
                        # 3. Solver and Mode Mapping
                        if hasattr(self, 'solver'):
                            # Extract numeric index from strings like "7: SP-Drift Diffusion"
                            import re
                            match = re.search(r'(\d+)', str(self.solver))
                            if match:
                                self.comp_scheme = int(match.group(1))
                                self.computation_scheme = self.comp_scheme
                        
                        if hasattr(self, 'device_type'):
                            if "Solar Cell" in str(self.device_type):
                                self.photovoltaic_mode = True
                            else:
                                self.photovoltaic_mode = False
                        
                        # Ensure numeric types for critical params
                        for attr in ['vmin', 'vmax', 'G_optical', 'alpha']:
                            if hasattr(self, attr):
                                try: setattr(self, attr, float(getattr(self, attr)))
                                except: pass
                
                with open(inputFile, 'r') as f:
                    data = json.load(f)
                inputfile_import = InputObject(data)
                inputfile_import.__file__ = inputFile
            else:
                # Load input file with importlib
                spec = importlib.util.spec_from_file_location(module_name, inputFile)
                if spec and spec.loader:
                    inputfile_import = importlib.util.module_from_spec(spec)
                    sys.modules[module_name] = inputfile_import
                    spec.loader.exec_module(inputfile_import)
                else:
                    print(f"Could not load input file: {inputFile}")
                    sys.exit(1)
                
            # Add to log
            logger.info("Inputfile is %s", inputFile)
        else:
            print("Please provide input file with -i argument")
            sys.exit()

        if args.drawfigures == True:
            # Control the drawing data at the end of the calculation
            drawFigures = True
        

    except getopt.error as err:
        # output error, and return with an error code
        print (str(err))

    output_directory = os.path.join(os.getcwd(), Path(inputFile).stem + "_output")

    #If output directory is not available, make one.
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    initialize_logger()

    os.sys.stderr.write("WARNING: Aestimo 1D logs in the output directory.\n")

    run_aestimo(inputfile_import)

else:
    # When imported as a module or default run without arguments?
    # Actually this else block runs if __name__ != "__main__", which means it's imported.
    # But this code block is inside `if __name__ == "__main__":` ?
    # Wait, looking at file... lines 4643 is `if __name__ == "__main__":`
    # The `else` at 4723 is matched to `if args.inputfile is not None:` ?
    # Let's check indentation.
    # Line 4679: if args.inputfile is not None:
    # Line 4723: else:
    #     output_directory = os.path.join(os.getcwd(), 'output')
    
    # Yes. This else handles the case where no input file is provided but it survived the earlier check?
    # Line 4656: if args is None: ... sys.exit()
    # But argparse handles this. 
    # Actually if args.inputfile is None, we print "Please provide..." and exit at 4700.
    # So the `else` block at 4723 is effectively dead code or unreachable given current logic?
    # Or maybe it was intended for something else.
    # Regardless, let's leave it as 'output' or maybe 'aestimo_output'.
    # User asked for "each example file should have output folder in its name".
    
    output_directory = os.path.join(os.getcwd(), 'output')

    #If output directory is not available, make one.
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    initialize_logger()

    os.sys.stderr.write("WARNING: Aestimo 1D logs automatically to aestimo.log in the output directory.\n")

