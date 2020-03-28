#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This is the 3x3 k.p aestimo calculator for valence band calculations 
   (Numpy version, there is no classic version for valence band calculations).

It can be used similarly to the aestimo.py module. aestimo_eh.py can be used as 
a script or a library.

To use as a script, define the simulation in a python file. See the following 
sample files for examples on usage and the required parameters:
    sample-qw-barrierdope-p.py
    sample-qw-barrierdope-p_cdzno.py
    sample-qw-barrierdope-p_ingran.py
    sample-multi-qw-barrierdope-p.py
    sample-multi-qw-barrierdope-p_ingran.py   
and then run aestimo on the command line as
  ./aestimo.py -i <input file>
Since we are abusing the python module system, the input 'file' needs to be 
importable by the aestimo script. Alternatively, define the input file in the
config module using the inputfilename parameter.

To use aestimo_eh.py as a library, first create an instance of the StructureFrom
class which builds the arrays describing a structure from the same input 
parameters that are found in the sample files. A simple list format is used to 
describes the structure's layers.
"""
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
"""
__version__ = "2.0"
import time

time0 = time.time()  # timing audit
# from scipy.optimize import fsolve
import matplotlib.pyplot as pl
import numpy as np

alen = np.alen
import os,sys

#importing examples directory
examplesdir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'examples'))
sys.path.append(examplesdir)

from math import log, exp, sqrt
import VBHM
from scipy import linalg
from VBHM import qsv, VBMAT1, VBMAT2, VBMAT_V, CBMAT, CBMAT_V, VBMAT_V_2
import config, database
from aestimo_poisson1d import (
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
from aestimo_poisson1d import (
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
import DDGgummelmap
from DDGgummelmap import DDGgummelmap
import DDNnewtonmap
from DDNnewtonmap import DDNnewtonmap
import func_lib
from func_lib import Ubernoulli

# --------------------------------------
import logging

logger = logging.getLogger("aestimo")

if not os.path.isdir(os.path.abspath(os.path.join(examplesdir, config.output_directory))):
    os.makedirs(os.path.abspath(os.path.join(examplesdir, config.output_directory)))

hdlr = logging.FileHandler(os.path.abspath(os.path.join(examplesdir, os.path.join(config.output_directory,config.logfile))))
formatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
# stderr
ch = logging.StreamHandler()
formatter2 = logging.Formatter("%(levelname)s %(message)s")
ch.setFormatter(formatter2)
logger.addHandler(ch)
# LOG level can be INFO, WARNING, ERROR
logger.setLevel(logging.INFO)
if not (config.messagesoff):
    os.sys.stderr.write(
        "WARNING aestimo_eh logs automatically to aestimo.log in the example's directory.\n"
    )
# --------------------------------------

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
# logger.info("Aestimo is starting...")
# Input Class
# -------------------------------------


def round2int(x):
    """int is sensitive to floating point numerical errors near whole numbers,
    this moves the discontinuity to the half interval. It is also equivalent
    to the normal rules for rounding positive numbers."""
    # int(x + (x>0) -0.5) # round2int for positive and negative numbers
    return int(x + 0.5)


def vegard1(first, second, mole):
    return first * mole + second * (1 - mole)


class Structure:
    def __init__(self, database, **kwargs):
        """This class holds details on the structure to be simulated.
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

        logger.info(
            "Total material number in database: %d", (totalmaterial + totalalloy)
        )

    def create_structure_arrays(self):
        """ initialise arrays/lists for structure"""
        # self.N_wells_real0=sum(sum(np.char.count(self.material,'w')))
        self.N_wells_real0 = sum(
            np.char.count([layer[6] for layer in self.material], "w")
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
            exit()
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
                fi_e[startindex:finishindex] = (
                    matprops["Band_offset"] * matprops["Eg"] * q
                )  # Joule
                if mat_crys_strc == "Zincblende":
                    a0_sub[startindex:finishindex] = matprops["a0"] * 1e-10
                    C11[startindex:finishindex] = matprops["C11"] * 1e10
                    C12[startindex:finishindex] = matprops["C12"] * 1e10
                    GA1[startindex:finishindex] = matprops["GA1"]
                    GA2[startindex:finishindex] = matprops["GA2"]
                    GA3[startindex:finishindex] = matprops["GA3"]
                    Ac[startindex:finishindex] = matprops["Ac"] * q
                    Av[startindex:finishindex] = matprops["Av"] * q
                    B[startindex:finishindex] = matprops["B"] * q
                    delta[startindex:finishindex] = matprops["delta"] * q
                    fi_h[startindex:finishindex] = (
                        -(1 - matprops["Band_offset"]) * matprops["Eg"] * q
                    )  # Joule  #-0.8*q-(1-matprops['Band_offset'])*matprops['Eg']*q #Joule
                    eps[startindex:finishindex] = matprops["epsilonStatic"] * eps0
                    a0[startindex:finishindex] = matprops["a0"] * 1e-10
                    TAUN0[startindex:finishindex] = matprops["TAUN0"]
                    TAUP0[startindex:finishindex] = matprops["TAUP0"]
                    mun0[startindex:finishindex] = matprops["mun0"]
                    mup0[startindex:finishindex] = matprops["mup0"]

                    Cn0[startindex:finishindex] = matprops["Cn0"] * 1e-12
                    Cp0[startindex:finishindex] = matprops["Cp0"] * 1e-12
                    BETAN[startindex:finishindex] = matprops["BETAN"]
                    BETAP[startindex:finishindex] = matprops["BETAP"]
                    VSATN[startindex:finishindex] = matprops["VSATN"]
                    VSATP[startindex:finishindex] = matprops["VSATP"]
                if mat_crys_strc == "Wurtzite":
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
                        -(1 - matprops["Band_offset"]) * matprops["Eg"] * q
                    )
                    Psp[startindex:finishindex] = matprops["Psp"]
                    TAUN0[startindex:finishindex] = matprops["TAUN0"]
                    TAUP0[startindex:finishindex] = matprops["TAUP0"]
                    mun0[startindex:finishindex] = matprops["mun0"]
                    mup0[startindex:finishindex] = matprops["mup0"]

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
                Eg = (
                    x * mat1["Eg"]
                    + (1 - x) * mat2["Eg"]
                    - alloyprops["Bowing_param"] * x * (1 - x)
                )  # eV
                fi_e[startindex:finishindex] = (
                    alloyprops["Band_offset"] * Eg * q
                )  # for electron. Joule
                a0_sub[startindex:finishindex] = alloyprops["a0_sub"] * 1e-10
                TAUN0[startindex:finishindex] = alloyprops["TAUN0"]
                TAUP0[startindex:finishindex] = alloyprops["TAUP0"]

                BETAN[startindex:finishindex] = alloyprops["BETAN"]
                BETAP[startindex:finishindex] = alloyprops["BETAP"]
                VSATN[startindex:finishindex] = alloyprops["VSATN"]
                VSATP[startindex:finishindex] = alloyprops["VSATP"]
                if mat_crys_strc == "Zincblende":
                    C11[startindex:finishindex] = (
                        x * mat1["C11"] + (1 - x) * mat2["C11"]
                    ) * 1e10
                    C12[startindex:finishindex] = (
                        x * mat1["C12"] + (1 - x) * mat2["C12"]
                    ) * 1e10
                    GA1[startindex:finishindex] = (
                        x * mat1["GA1"] + (1 - x) * mat2["GA1"]
                    )
                    GA2[startindex:finishindex] = (
                        x * mat1["GA2"] + (1 - x) * mat2["GA2"]
                    )
                    GA3[startindex:finishindex] = (
                        x * mat1["GA3"] + (1 - x) * mat2["GA3"]
                    )
                    Ac_alloy = x * mat1["Ac"] + (1 - x) * mat2["Ac"]
                    Ac[startindex:finishindex] = Ac_alloy * q
                    Av_alloy = x * mat1["Av"] + (1 - x) * mat2["Av"]
                    Av[startindex:finishindex] = Av_alloy * q
                    B_alloy = x * mat1["B"] + (1 - x) * mat2["B"]
                    B[startindex:finishindex] = B_alloy * q
                    delta_alloy = x * mat1["delta"] + (1 - x) * mat2["delta"]
                    delta[startindex:finishindex] = delta_alloy * q
                    fi_h[startindex:finishindex] = (
                        -(1 - alloyprops["Band_offset"]) * Eg * q
                    )  # -(-1.33*(1-x)-0.8*x)for electron. Joule-1.97793434e-20 #
                    eps[startindex:finishindex] = (
                        x * mat1["epsilonStatic"] + (1 - x) * mat2["epsilonStatic"]
                    ) * eps0
                    a0[startindex:finishindex] = (
                        (1 - x) * mat1["a0"] + x * mat2["a0"]
                    ) * 1e-10
                    cb_meff_alpha[startindex:finishindex] = alloyprops["m_e_alpha"] * (
                        mat2["m_e"] / cb_meff_alloy
                    )  # non-parabolicity constant for alloy. THIS CALCULATION IS MOSTLY WRONG. MUST BE CONTROLLED. SBL

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
                if mat_crys_strc == "Wurtzite":
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
                    """
                        I. Vurgaftman, J.R. Meyer, L.R. RamMohan, J. Appl. Phys. 89 (2001) 5815.
                        C. K. Williams, T. H. Glisson, J. R. Hauser, and M. A. Littlejohn, J. Electron. Mater. 7, 639 (1978).                       
                        """
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
            matRole = layer[6]
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

            dop_profile = self.dop_profile
            if layer[5] == "n":
                dop[startindex:finishindex] = (
                    layer[4] * 1e6 + dop_profile[startindex:finishindex] + 1
                )  # charge density in m**-3 (conversion from cm**-3)
            elif layer[5] == "p":
                dop[startindex:finishindex] = (
                    -layer[4] * 1e6 + dop_profile[startindex:finishindex] - 1
                )  # charge density in m**-3 (conversion from cm**-3)
            else:
                dop[startindex:finishindex] = dop_profile[startindex:finishindex] + 1
        """
        Here we remove barriers that are less than the anti_crossing_length
        so we can constructe the new well boundary using the resulted barrier boundary
        """
        brr = 0
        anti_crossing_length = config.anti_crossing_length * 1e-9
        if not (self.Quantum_Regions):
            for J in range(2, N_wells_virtual2 - 1):
                if barrier_len[J] * dx <= anti_crossing_length:
                    brr += 1
            brr_vec = np.zeros(brr)
            brr2 = 0
            for J in range(2, N_wells_virtual2 - 1):
                if barrier_len[J] * dx <= anti_crossing_length:
                    brr2 += 1
                    brr_vec[brr2 - 1] = J + 1 - brr2
            for I in range(0, brr):
                barrier_boundary = np.delete(barrier_boundary, brr_vec[I], 0)
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
        self.Fapp = inputfile.Fapplied
        self.vmax = inputfile.vmax
        self.vmin = inputfile.vmin
        self.Each_Step = inputfile.Each_Step
        self.surface = inputfile.surface
        self.T = inputfile.T
        self.subnumber_h = inputfile.subnumber_h
        self.subnumber_e = inputfile.subnumber_e
        self.comp_scheme = inputfile.computation_scheme
        self.dx = inputfile.gridfactor * 1e-9  # grid in m
        self.maxgridpoints = inputfile.maxgridpoints
        self.mat_crys_strc = inputfile.mat_type
        # Loading material list
        self.material = inputfile.material

        totallayer = alen(self.material)
        if not (config.messagesoff):
            logger.info("Total layer number: %s", totallayer)
        # Calculate the required number of grid points
        self.x_max = (
            sum([layer[0] for layer in self.material]) * 1e-9
        )  # total thickness (m)
        self.n_max = int(self.x_max / self.dx)
        # Check on n_max
        max_val = inputfile.maxgridpoints

        self.dop_profile = inputfile.dop_profile
        self.Quantum_Regions_boundary = inputfile.Quantum_Regions_boundary
        self.Quantum_Regions = inputfile.Quantum_Regions
        if self.n_max > max_val:
            logger.error(" Grid number is exceeding the max number of %d", max_val)
            exit()
        # Loading materials database #
        self.material_property = database.materialproperty
        totalmaterial = alen(self.material_property)

        self.alloy_property = database.alloyproperty
        totalalloy = alen(self.alloy_property)

        self.alloy_property_4 = database.alloyproperty4
        totalalloy += alen(self.alloy_property_4)
        if not (config.messagesoff):
            logger.info(
                "Total number of materials in database: %d"
                % (totalmaterial + totalalloy)
            )
        # Initialise arrays

        # cb_meff #conduction band effective mass (array, len n_max)
        # fi_e #Bandstructure potential (array, len n_max)
        # eps #dielectric constant (array, len n_max)
        # dop #doping distribution (array, len n_max)
        self.create_structure_arrays()


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
        if not (config.messagesoff):
            logger.warning(
                "Have processed all energy levels present and so can't be sure that Ef is below next higher energy level."
            )
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
        if not (config.messagesoff):
            logger.warning(
                "Have processed all energy levels present and so can't be sure that Ef is below next higher energy level."
            )
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

            """
            E31=(C11+C12)*D31+C13*D33
            E33=2*C13*D31+C33*D33
            print('E31=',E31,'E33=',E33)
            Ppz2=2*EXX*(E31-E33*C13/C33)
            print('Pps2=',Ppz2)
            """
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
                """
                xaxis = np.arange(0,n_max)*dx
                pl.plot(xaxis*1e6,Ppz_Psp0,'r')#,xaxis*1e6,Ppz_Psp0,'b'
                pl.xlabel('x [um]')
                pl.ylabel('Energy [eV]')
                pl.title('pze (0 (red) & 1 (bleu)) vs Position', fontsize=12)
                pl.legend(('Efn','Efp'),loc='best',fontsize=12)
                pl.grid(True)
                pl.show()
                ggggggggggggg
                """
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
    # HUPMAT1=np.zeros((n_max*3, n_max*3))
    # HUPMATC1=np.zeros((n_max, n_max))
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
    if model.mat_crys_strc == "Zincblende" and (model.N_wells_virtual - 2 != 0):
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
    if model.mat_crys_strc == "Wurtzite" and (model.N_wells_virtual - 2 != 0):
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
            """aestimo_eh doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """aestimo_eh doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        exit()
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
    if not (config.messagesoff):
        # print "Ntotal ",Ntotal,"m**-3"
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
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * abs(dop[i])))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ns * ni[i]))
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
        Eg_[i] = fi_e[i] - fi_h[i]

        fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
    """
    fi_e-=fi_e[0]
    fi_h-=fi_e[0]
    
    fi_e+=offset0+offset1#+kb*T*np.log(Nv/Nc)/2
    fi_h+=offset0+offset1#+kb*T*np.log(Nv/Nc)/2
    
    
    pl.plot(xaxis, Half_Eg/q,'k')
    pl.xlabel('Position (m)')
    pl.ylabel('electrons  and and holes concentrations (cm-3)' )
    pl.title('electrons (red) and holes (blue)')
    pl.grid(True)
    """
    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )
        # exit()
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
        if not (config.messagesoff):
            logger.info("Iteration: %d", iteration)
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
            if iteration == 1:
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
            damping = 0.15  # 0.1 works between high and low doping
            if config.predic_correc:
                (
                    E_statec_general,
                    E_state_general,
                    wfe_general,
                    wfh_general,
                    meff_statec_general,
                    meff_state_general,
                ) = (
                    E_statec_general0,
                    E_state_general0,
                    wfe_general0,
                    wfh_general0,
                    meff_statec_general0,
                    meff_state_general0,
                )
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
        else:
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
            damping = 1
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
    if not (config.messagesoff):
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
    BETAN = model.BETAN
    BETAP = model.BETAP
    VSATN = model.VSATN
    VSATP = model.VSATP
    if comp_scheme in (4, 5, 6):
        logger.error(
            """aestimo_eh doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """aestimo_eh doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        exit()
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
    if not (config.messagesoff):
        # print "Ntotal ",Ntotal,"m**-3"
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
        ni[i] = sqrt(
            Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T))
        )  # Intrinsic carrier concentration [1/m^3] kb*T/q
        # print("%.3E" % (ni[i]*1e-6))
        # print(fi_e[i]-fi_h[i])
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * abs(dop[i])))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ni[i]))
        if dop[i] == 1:
            dop[i] *= ni[i]
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
        fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
    n = result.nf_result / ni
    p = result.pf_result / ni
    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )
        # exit()
    # STARTING SELF CONSISTENT LOOP
    time2 = time.time()  # timing audit
    iteration = 1  # iteration counter
    # previousE0= 0   #(meV) energy of zeroth state for previous iteration(for testing convergence)
    # fitot = fi_h + Vapp #For initial iteration sum bandstructure and applied field
    # fitotc = fi_e + Vapp
    Va_max = vmax  # 1.8#input()0.625
    # Va_max=0.625#input()0.625
    dVa = 0.5 * Vt  # input()0.01
    dVa = dVa / Vt
    Each_Step = dVa
    vmin = 0.0
    Total_Steps = int(((Va_max - vmin) / Vt) / (Each_Step))
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

    fi_stat = fi
    if Va_max == 0:
        print("Va_max=0")
    else:
        print("Convergence of the Gummel cycles")
        vindex = 0
        for vindex in range(0, Total_Steps):
            # if vindex>int(Total_Steps*4/5):
            Ppz_Psp = Ppz_Psp_tmp
            # Start Va increment loop
            Va = Each_Step * vindex
            if vindex == 0:
                fi[0] += 0.0  # Apply potential to Anode (1st node)
            else:
                fi[0] += Each_Step
            flag_conv_2 = True  # Convergence of the Poisson loop
            #% Initialize the First and Last Node for Poisson's eqn

            Va_t[vindex] = Va
            print("Va_t[", vindex, "]=", Va_t[vindex] * Vt)
            print("vindex=", vindex)
            # previousE0= 2   #(meV) energy of zeroth state for previous iteration(for testing convergence)
            while flag_conv_2:
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
                )
                #
                mun, mup = Mobility2(
                    mun0, mup0, fi, Vt, Ldi, VSATN, VSATP, BETAN, BETAP, n_max, dx
                )
                ########### END of FIELD Dependant Mobility Calculation ###########
                n, p = Continuity2(n, p, mun, mup, fi, Vt, Ldi, n_max, dx, TAUN0, TAUP0)
                ####################### END of HOLE Continuty Solver ###########
                # End of WHILE Loop for Poisson's eqn solver
                # print('inside while loop')
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
            Jtotal = Jelec + Jhole
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
    time3 = time.time()  # timing audit
    if not (config.messagesoff):
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
    return results


################################################
def Poisson_Schrodinger_DD_test(result, model):
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
    BETAN = model.BETAN
    BETAP = model.BETAP
    VSATN = model.VSATN
    VSATP = model.VSATP
    if comp_scheme in (4, 5, 6):
        logger.error(
            """aestimo_eh doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """aestimo_eh doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        exit()
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
    if not (config.messagesoff):
        # print "Ntotal ",Ntotal,"m**-3"
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
        ni[i] = sqrt(
            Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T))
        )  # Intrinsic carrier concentration [1/m^3] kb*T/q
        # print("%.3E" % (ni[i]*1e-6))
        # print(fi_e[i]-fi_h[i])
        if dop[i] == 1:
            dop[i] *= ni[i]
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * abs(dop[i])))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ni[i]))
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
        fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
    n = result.nf_result / ni
    p = result.pf_result / ni
    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )
        exit()
    # STARTING SELF CONSISTENT LOOP
    time2 = time.time()  # timing audit
    iteration = 1  # iteration counter
    previousE0 = 0  # (meV) energy of zeroth state for previous iteration(for testing convergence)
    previousfi0 = 0  # (meV) energy of  for previous iteration(for testing convergence)
    fitot = fi_h  # + Vapp #For initial iteration sum bandstructure and applied field
    fitotc = fi_e  # + Vapp
    Va_max = vmax  # 1.8#input()0.625
    # Va_max=0.625#input()0.625
    dVa = Each_Step  # *Vt#input()0.01
    dVa = dVa / Vt
    Each_Step = dVa
    Total_Steps = int(((Va_max - vmin) / Vt) / (Each_Step))
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
    fi_stat = fi
    if Va_max == 0:
        print("Va_max=0")
    else:
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
                fi[0] += Each_Step
            flag_conv_2 = True  # Convergence of the Poisson loop
            #% Initialize the First and Last Node for Poisson's eqn

            Va_t[vindex] = Va
            print("Va_t[", vindex, "]=", Va_t[vindex] * Vt)
            print("vindex=", vindex)
            # previousE0= 2   #(meV) energy of zeroth state for previous iteration(for testing convergence)
            while flag_conv_2:
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
                )
                ########### END of FIELD Dependant Mobility Calculation ###########
                n, p = Continuity3(
                    n, p, mun, mup, fi, fi_n, fi_p, Vt, Ldi, n_max, dx, TAUN0, TAUP0
                )
                # if config.quantum_effect:
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
            Jtotal = Jelec + Jhole
            """                
            pl.plot(xaxis*1e6,fitotc)
            pl.xlabel('Position (m)')
            pl.ylabel('Energy (meV)')
            pl.grid(True)
            """
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
    if not (config.messagesoff):
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
    return results


def Poisson_Schrodinger_DD_test_2(result, model):
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
    Cn0 = model.Cn0
    Cp0 = model.Cp0
    BETAN = model.BETAN
    BETAP = model.BETAP
    VSATN = model.VSATN
    VSATP = model.VSATP

    if comp_scheme in (4, 5, 6):
        logger.error(
            """aestimo_eh doesn't currently include exchange interactions
        in its valence band calculations."""
        )
        exit()
    if comp_scheme in (1, 3, 6):
        logger.error(
            """aestimo_eh doesn't currently include nonparabolicity effects in 
        its valence band calculations."""
        )
        exit()
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
    if not (config.messagesoff):
        # print "Ntotal ",Ntotal,"m**-3"
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
        ni[i] = sqrt(
            Nc[i] * Nv[i] * exp(-(fi_e[i] - fi_h[i]) / (kb * T))
        )  # Intrinsic carrier concentration [1/m^3] kb*T/q
        Ld_n_p[i] = sqrt(eps[i] * Vt / (q * abs(dop[i])))
        Ldi[i] = sqrt(eps[i] * Vt / (q * ni[i]))
        if dop[i] == 1:
            dop[i] *= ni[i]
        Half_Eg[i] = (fi_e[i] - fi_h[i]) / 2
        fi_e[i] = Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
        fi_h[i] = -Half_Eg[i] - kb * T * log(Nv[i] / Nc[i]) / 2
    n = result.nf_result / ni
    p = result.pf_result / ni
    if dx > min(Ld_n_p[:]) and 1 == 2:
        logger.error(
            """You are setting the grid size %g nm greater than the extrinsic Debye lengths %g nm""",
            dx * 1e9,
            min(Ld_n_p[:]) * 1e9,
        )
        # exit()
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
    ###############################################################
    len_ = xaxis[n_max - 1]

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
    # J=np.zeros((Total_Steps,n_max-1))
    lambda2 = np.zeros((Total_Steps, n_max))
    DV = np.zeros(Total_Steps)
    Emax = np.zeros(Total_Steps)

    nn, pp, fi_out = equi_np_fi(iteration, dop, Ppz_Psp, n_max, ni, model, Vt, surface)
    # xn = xm+1e-7
    # xp = xm-1e-7
    ## Scaling coefficients
    xs = len_
    ns1 = np.linalg.norm(dop, np.inf)

    ns2 = np.linalg.norm(Ppz_Psp_tmp, np.inf)

    ns = max(ns1, ns2)
        
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

    Vs = Vt
    us = max(max(mun0), max(mup0))
    Js = xs / (us * Vs * q * ns)
    xbar = len_  # [m]
    Vbar = Vt  # [V]
    mubar = max(max(mun0), max(mup0))  # [m^2 V^{-1} s^{-1}]
    tbar = xbar ** 2 / (mubar * Vbar)  # [s]
    Rbar = ns / tbar
    # [m^{-3} s^{-1}]
    CAubar = Rbar / ns ** 3  # [m^6 s^{-1}]
    idata.Cn = Cn0 / CAubar
    idata.Cp = Cp0 / CAubar
    ###############################################################
    if vmax == 0:
        print("Va_max=0")
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
            n_[vindex, :] = nn * ni
            p_[vindex, :] = pp * ni
            #
            Fn = Va * (xaxis <= xm)
            Fp = Fn
            #
            V_[vindex, :] = Fn - Vt * np.log(p_[vindex, :] / ni)
            #
            ## Scaling

            Fn_[vindex, :] = Fn - Vs * np.log(ni / ns)
            Fp_[vindex, :] = Fp + Vs * np.log(ni / ns)
            #
            idata.l2 = (Vs * eps[0 : n_max - 1]) / (q * ns * xs ** 2)
            idata.nis = ni / ns
            idata.dop = dop / ns
            idata.Ppz_Psp = Ppz_Psp / ns
            # mun,mup=Mobility2(mun0,mup0,fi,Vt,Ldi,VSATN,VSATP,BETAN,BETAP,n_max,dx)
            idata.mun = mun0 / us
            idata.mup = mup0 / us

            # sinodes = np.arange(len(xaxis))
            idata.TAUN0 = TAUN0 / tbar  # np.inf
            idata.TAUP0 = TAUP0 / tbar  # np.inf
            idata.theta = ni / ns

            idata.n = n_[vindex, :] / ns
            idata.p = p_[vindex, :] / ns
            idata.V = V_[vindex, :] / Vs
            idata.Fn = Fn_[vindex, :] / Vs
            idata.Fp = Fp_[vindex, :] / Vs
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
            toll = 1e-3
            maxit = 10
            ptoll = 1e-10
            pmaxit = 30
            verbose = 0

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

            [odata, it, res] = DDNnewtonmap(
                ni, fi_e, fi_h, xin, odata, toll, maxit, verbose, model, Vt
            )

            n_[vindex, :] = odata.n
            p_[vindex, :] = odata.p
            V_[vindex, :] = odata.V

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

            Bp = Ubernoulli(
                (V_[vindex, 1:n_max] - V_[vindex, 0 : n_max - 1])
                + (fi_n[1:n_max] - fi_n[0 : n_max - 1]),
                1,
            )
            Bm = Ubernoulli(
                (V_[vindex, 1:n_max] - V_[vindex, 0 : n_max - 1])
                + (fi_p[1:n_max] - fi_p[0 : n_max - 1]),
                0,
            )
            Jn[vindex, 0 : n_max - 1] = (
                -odata.mun[0 : n_max - 1]
                * (n_[vindex, 1:n_max] * Bp - n_[vindex, 0 : n_max - 1] * Bm)
                / (xin[1:n_max] - xin[0 : n_max - 1])
            )
            Jp[vindex, 0 : n_max - 1] = (
                odata.mup[0 : n_max - 1]
                * (p_[vindex, 1:n_max] * Bm - p_[vindex, 0 : n_max - 1] * Bp)
                / (xin[1:n_max] - xin[0 : n_max - 1])
            )
        ## Descaling
        n_ = n_ * ns
        p_ = p_ * ns
        V_ = V_ * Vs
        # J = abs (Jp+Jn)*Js
        Jtotal = abs(Jp + Jn) * us * q * ns
        Jtotal[:, n_max - 1] = Jtotal[:, n_max - 2]
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
        av_curr = Jtotal[:, n_max - 1]
        for i in range(1, n_max - 1):
            Ec_result[i] = fi_e[i] / q - V_[vindex, i]  # Values from the second Node%
            Ev_result[i] = fi_h[i] / q - V_[vindex, i]  # Values from the second Node%
            Ei_result[i] = Ec_result[i] - ((fi_e[i] - fi_h[i]) / (2 * q))
            ro_result[i] = -q * (n_[vindex, i] - p_[vindex, i] - ns * dop[i])
            el_field1_result[i] = -(V_[vindex, i + 1] - V_[vindex, i]) / (dx)
            el_field2_result[i] = -(V_[vindex, i + 1] - V_[vindex, i - 1]) / (2 * dx)
            #Efn_result[i] = Ec_result[i] + Vt * log(n_[vindex, i] / Nc[i])
            #Efp_result[i] = Ev_result[i] - Vt * log(p_[vindex, i] / Nv[i])
        # Efn_result=Fn_[vindex,:]
        # Efp_result=Fp_[vindex,:]
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
        if model.N_wells_virtual - 2 != 0 and 1 == 2:
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
    if not (config.messagesoff):
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
    # results.meff_state = meff_state
    results.E_statec = E_statec
    results.N_statec = N_statec
    # results.meff_statec = meff_statec
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
    return results


def save_and_plot2(result, model):
    xaxis = result.xaxis
    output_directory = config.output_directory + "_eh"

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    def saveoutput(fname, datatuple, header=""):
        fname2 = os.path.join(output_directory, fname)
        np.savetxt(
            fname2, np.column_stack(datatuple), fmt="%.6e", delimiter=" ", header=header
        )

    # Plotting results
    # if config.Drift_Diffusion_out:
    # saveoutput("av_curr.dat",(result.Va_t*Vt,result.av_curr*1e-4))
    for jjj in range(result.Total_Steps - 1, result.Total_Steps):
        vtt = result.Va_t[jjj]
        vt = vtt * Vt
        if config.Drift_Diffusion_out:
            if config.sigma_out:
                saveoutput("sigma_eh_%.2f.dat" % vt, (xaxis, result.ro_result))
            if config.electricfield_out:
                saveoutput(
                    "efield_eh_%.2f.dat" % vt,
                    (xaxis, result.el_field1_result, result.el_field2_result),
                )
            if config.potential_out:
                saveoutput(
                    "potn_eh_%.2f.dat" % vt,
                    (xaxis * 1e2, result.Ec_result, result.Ev_result),
                )
                saveoutput(
                    "np_data0_%.2f.dat" % vt,
                    (xaxis * 1e2, result.nf_result * 1e-6, result.pf_result * 1e-6),
                )
            if config.states_out and 1 == 2:
                for j in range(1, result.N_wells_virtual - 1):
                    rel_meff_state = [
                        meff / m_e for meff in result.meff_state_general[j]
                    ]  # going to report relative effective mass.
                    columns = (
                        range(model.subnumber_h),
                        result.E_state_general[j],
                        result.N_state_general[j],
                        rel_meff_state,
                    )
                    # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
                    header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
                    saveoutput(
                        "states_h_QWR%d_%.2f.dat" % (j, vt), columns, header=header
                    )
                    if config.probability_out:
                        saveoutput(
                            "wavefunctions_h_QWR%d_%.2f.dat" % (j, vt),
                            (xaxis, result.wfh_general[j].transpose()),
                        )
            if config.states_out and 1 == 2:
                for j in range(1, result.N_wells_virtual - 1):
                    rel_meff_statec = [
                        meff / m_e for meff in result.meff_statec_general[j]
                    ]  # going to report relative effective mass.
                    columns = (
                        range(model.subnumber_e),
                        result.E_statec_general[j],
                        result.N_statec_general[j],
                        rel_meff_statec,
                    )
                    # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
                    header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
                    saveoutput(
                        "states_e_QWR%d_%.2f.dat" % (j, vt), columns, header=header
                    )
                    if config.probability_out:
                        saveoutput(
                            "wavefunctions_e_QWR%d_%.2f.dat" % (j, vt),
                            (xaxis, result.wfe_general[j].transpose()),
                        )
    if config.resultviewer:
        span = np.ones(100000000)
        fig1 = pl.figure(figsize=(10, 8))
        pl.suptitle("Aestimo Results")
        pl.subplot(1, 1, 1)
        pl.plot(
            xaxis * 1e6,
            result.Ec_result,
            xaxis * 1e6,
            result.Ev_result,
            xaxis * 1e6,
            result.Ei_result,
            xaxis * 1e6,
            result.Efn_result,
            "r",
            xaxis * 1e6,
            result.Efp_result,
            "b",
        )
        if model.N_wells_virtual - 2 != 0:
            for j in range(1, result.N_wells_virtual - 1):
                I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
                i1 = I1 - I1
                i2 = I2 - I1
                for levelc, statec in zip(
                    result.E_statec_general[j, :], result.wfe_general[j, :, :]
                ):
                    # pl.axhline(levelc,0.1,0.9, hold=True,color='g',ls='--')
                    pl.plot(
                        xaxis[I1:I2] * 1e6,
                        statec[i1:i2] * config.wavefunction_scalefactor * 1e-3
                        + levelc * 1e-3,
                        "b",
                    )
                    pl.plot(
                        xaxis[I1:I2] * 1e6, levelc * span[I1:I2] * 1e-3, "g", ls="--"
                    )
                for level, state in zip(
                    result.E_state_general[j, :], result.wfh_general[j, :, :]
                ):
                    # pl.axhline(level,0.1,0.9,color='g',ls='--')
                    pl.plot(
                        xaxis[I1:I2] * 1e6,
                        state[i1:i2] * config.wavefunction_scalefactor * 1e-3
                        + level * 1e-3,
                        "b",
                    )
                    pl.plot(
                        xaxis[I1:I2] * 1e6, level * span[I1:I2] * 1e-3, "g", ls="--"
                    )
        pl.xlabel("x [um]")
        pl.ylabel("Energy [eV]")
        pl.title("Quasi Fermi Levels (Efn (red) & Efp (bleu)) vs Position", fontsize=12)
        pl.legend(("Ec", "Ev", "Ei", "Efn", "Efp"), loc="best", fontsize=12)
        pl.grid(True)

        fig2 = pl.figure(figsize=(10, 8))
        pl.suptitle(
            "1D Drift Diffusion Model for pn Diodes Results - at Applied Bias (%.2f)"
            % vt,
            fontsize=12,
        )
        pl.subplots_adjust(hspace=0.4, wspace=0.4)

        pl.subplot(2, 2, 1)
        pl.plot(xaxis * 1e6, result.ro_result * 1e-6)
        pl.xlabel("x [um]")
        pl.ylabel("Total Charge Density [C/cm^3]")
        pl.title("Total Charge Density vs Position ", fontsize=12)
        pl.legend(("Total Charge"), loc="best", fontsize=12)
        pl.grid(True)
        # Plotting Efield
        # figure(1)
        pl.subplot(2, 2, 2)
        pl.plot(
            xaxis * 1e6,
            result.el_field1_result * 1e-8,
            "r",
            xaxis * 1e6,
            result.el_field2_result * 1e-8,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Electric Field 1(red) & 2 (bleu) [MV/cm]")
        pl.title("Field Profile vs Position ", fontsize=12)
        pl.legend(("Electric Field 1", "Electric Field 2"), loc="best", fontsize=12)
        pl.grid(True)
        # Plotting Potential
        # figure(2)
        pl.subplot(2, 2, 3)
        pl.plot(xaxis * 1e6, result.Ec_result)
        pl.xlabel("x [um]")
        pl.ylabel("Conduction Band Energy (eV)")
        pl.title("Conduction Band vs Position ", fontsize=12)
        pl.legend(("Conduction Band"), loc="best", fontsize=12)
        pl.grid(True)
        # Plotting State(s)
        # figure(3)
        pl.subplot(2, 2, 4)
        pl.plot(result.Va_t , result.av_curr * 1e-4)
        pl.xlabel("Va [V]")
        pl.ylabel("Total Current Density [Amp/cm^2]")
        pl.title("I vs V Plot", fontsize=12)
        pl.legend(("Total Current"), loc="best", fontsize=12)
        pl.grid(True)
        pl.show()

        fig3 = pl.figure(figsize=(10, 8))
        pl.suptitle(
            "1D Drift Diffusion Model for pn Diodes Results - at Applied Bias (%.2f)"
            % vt,
            fontsize=12,
        )
        pl.subplots_adjust(hspace=0.4, wspace=0.4)
        pl.subplot(2, 2, 1)
        pl.plot(
            xaxis * 1e6,
            result.nf_result * 1e-6,
            "r",
            xaxis * 1e6,
            result.pf_result * 1e-6,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Electron  & Hole  Densities [1/cm^3]")
        pl.title("Electron (red) & Hole (bleu) Densities vs Position ", fontsize=12)
        pl.legend(("Electron", "Hole"), loc="best", fontsize=12)
        pl.grid(True)

        pl.subplot(2, 2, 2)
        pl.plot(
            xaxis * 1e6,
            result.Ec_result,
            xaxis * 1e6,
            result.Ev_result,
            xaxis * 1e6,
            result.Ei_result,
            xaxis * 1e6,
            result.Efn_result,
            "r",
            xaxis * 1e6,
            result.Efp_result,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Energy [eV]")
        pl.title("Quasi Fermi Levels (Efn (red) & Efp (bleu)) vs Position", fontsize=12)
        pl.legend(("Ec", "Ev", "Ei", "Efn", "Efp"), loc="best", fontsize=12)
        pl.grid(True)

        pl.subplot(2, 2, 3)
        pl.plot(xaxis * 1e6, Vt * result.fi_result)
        pl.xlabel("x [um]")
        pl.ylabel("Potential [V]")
        pl.title("Potential vs Position", fontsize=12)
        pl.legend(("fi"), loc="best", fontsize=12)
        pl.grid(True)

        pl.subplot(2, 2, 4)
        pl.plot(
            xaxis * 1e6, result.Efn_result, "r", xaxis * 1e6, result.Efp_result, "b"
        )
        pl.xlabel("x [um]")
        pl.ylabel("Energy [eV]")
        pl.title("Quasi Fermi Levels (Efn (red) & Efp (bleu)) vs Position", fontsize=12)
        pl.legend(("Efn", "Efp"), loc="best", fontsize=12)
        pl.grid(True)
        pl.show()
    return [fig1, fig2, fig3]


def save_and_plot(result, model):
    xaxis = result.xaxis
    
    output_directory = config.output_directory
    output_directory = os.path.join(examplesdir, output_directory)
    
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    def saveoutput(fname, datatuple, header=""):
        fname2 = os.path.join(output_directory, fname)
        np.savetxt(
            fname2, np.column_stack(datatuple), fmt="%.6e", delimiter=" ", header=header
        )

    if config.sigma_out:
        saveoutput("sigma_eh.dat", (xaxis, result.ro_result))
    if config.electricfield_out:
        saveoutput(
            "efield_eh.dat", (xaxis, result.el_field1_result, result.el_field2_result)
        )
    if config.potential_out:
        saveoutput("potn_eh.dat", (xaxis * 1e2, result.fitotc / q, result.fitot / q))
        saveoutput(
            "np_data0.dat",
            (xaxis * 1e2, result.nf_result * 1e-6, result.pf_result * 1e-6),
        )
    if config.states_out:
        for j in range(1, result.N_wells_virtual - 1):
            rel_meff_state = [
                meff / m_e for meff in result.meff_state_general[j]
            ]  # going to report relative effective mass.
            columns = (
                range(model.subnumber_h),
                result.E_state_general[j],
                result.N_state_general[j],
                rel_meff_state,
            )
            # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
            header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
            saveoutput("states_h_QWR%d.dat" % j, columns, header=header)
            if config.probability_out:
                saveoutput(
                    "wavefunctions_h_QWR%d.dat" % j,
                    (xaxis, result.wfh_general[j].transpose()),
                )
    if config.states_out:
        for j in range(1, result.N_wells_virtual - 1):
            rel_meff_statec = [
                meff / m_e for meff in result.meff_statec_general[j]
            ]  # going to report relative effective mass.
            columns = (
                range(model.subnumber_e),
                result.E_statec_general[j],
                result.N_statec_general[j],
                rel_meff_statec,
            )
            # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
            header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
            saveoutput("states_e_QWR%d.dat" % j, columns, header=header)
            if config.probability_out:
                saveoutput(
                    "wavefunctions_e_QWR%d.dat" % j,
                    (xaxis, result.wfe_general[j].transpose()),
                )
    # Resultviewer
    if config.resultviewer:
        span = np.ones(100000000)
        fig1 = pl.figure(figsize=(10, 8))
        pl.suptitle("Aestimo Results")
        pl.subplot(1, 1, 1)
        pl.plot(xaxis, result.fitot * J2meV, "k", xaxis, result.fitotc * J2meV, "k")
        for j in range(1, result.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
            i1 = I1 - I1
            i2 = I2 - I1
            for levelc, statec in zip(
                result.E_statec_general[j, :], result.wfe_general[j, :, :]
            ):
                # pl.axhline(levelc,0.1,0.9, hold=True,color='g',ls='--')
                pl.plot(
                    xaxis[I1:I2],
                    statec[i1:i2] * config.wavefunction_scalefactor + levelc,
                    "b",
                )
                pl.plot(xaxis[I1:I2], levelc * span[I1:I2], "g", ls="--")
            for level, state in zip(
                result.E_state_general[j, :], result.wfh_general[j, :, :]
            ):
                # pl.axhline(level,0.1,0.9,color='g',ls='--')
                pl.plot(
                    xaxis[I1:I2],
                    state[i1:i2] * config.wavefunction_scalefactor + level,
                    "b",
                )
                pl.plot(xaxis[I1:I2], level * span[I1:I2], "g", ls="--")
            # pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
        pl.plot(xaxis, result.EF * span[0 : model.n_max], "r", ls="--")
        # pl.axhline(result.E_F,0.1,0.9,color='r',ls='--')
        pl.xlabel("Position (m)")
        pl.ylabel("Energy (meV)")
        pl.grid(True)

        fig2 = pl.figure(figsize=(10, 8))
        pl.suptitle(
            "1D Drift Diffusion Model for pn Diodes Results - at Equilibrium Condition ",
            fontsize=12,
        )
        pl.subplots_adjust(hspace=0.4, wspace=0.4)

        # Plotting Sigma
        # figure(0)
        pl.subplot(2, 2, 1)
        pl.plot(xaxis * 1e6, result.ro_result * 1e-6)
        pl.xlabel("x [um]")
        pl.ylabel("Total Charge Density [C/cm^3]")
        pl.title("Total Charge Density vs Position ", fontsize=10)
        pl.grid(True)

        # Plotting Efield
        # figure(1)
        pl.subplot(2, 2, 2)
        pl.plot(
            xaxis * 1e6,
            result.el_field1_result * 1e-8,
            "r",
            xaxis * 1e6,
            result.el_field2_result * 1e-8,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Electric Field  [MV/cm]")
        pl.title("Field Profile 1(red) & 2 (bleu) vs Position ", fontsize=10)
        pl.grid(True)

        # Plotting Potential
        # figure(2)
        pl.subplot(2, 2, 3)
        pl.plot(xaxis * 1e6, Vt * result.fi_result)
        pl.xlabel("x [um]")
        pl.ylabel("Potential [V]")
        pl.title("Potential vs Position", fontsize=12)
        pl.legend(("fi"), loc="best", fontsize=12)
        pl.grid(True)

        pl.subplot(2, 2, 4)
        pl.plot(
            xaxis * 1e6,
            result.nf_result * 1e-6,
            "r",
            xaxis * 1e6,
            result.pf_result * 1e-6,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Electron  & Hole  Densities [1/cm^3]")
        pl.title("Electron (red)& Hole (bleu) Densities vs Position ", fontsize=10)
        pl.grid(True)
        pl.show()
    return [fig1, fig2]


def run_aestimo(input_obj):
    """A utility function that performs the standard simulation run
    for 'normal' input files. Input_obj can be a dict, class, named tuple or 
    module with the attributes needed to create the StructureFrom class, see 
    the class implementation or some of the sample-*.py files for details."""
    if not (config.messagesoff):
        logger.info("Aestimo_eh is starting...")
    # Initialise structure class
    model = StructureFrom(input_obj, database)

    # Perform the calculation
    result = Poisson_Schrodinger(model)
    if model.comp_scheme == 7:
        result_dd = Poisson_Schrodinger_DD(result, model)
    if model.comp_scheme == 8:
        result_dd = Poisson_Schrodinger_DD_test(result, model)
    if model.comp_scheme == 9:
        result_dd = Poisson_Schrodinger_DD_test_2(result, model)
    time4 = time.time()  # timing audit
    if not (config.messagesoff):

        logger.info("total running time (inc. loading libraries) %g s", (time4 - time0))
        logger.info("total running time (exc. loading libraries) %g s", (time4 - time1))
    # Write the simulation results in files
    if model.comp_scheme == 7 or model.comp_scheme == 8 or model.comp_scheme == 2:
        save_and_plot(result, model)
    if model.comp_scheme == 7 or model.comp_scheme == 8 or model.comp_scheme == 9:
        save_and_plot2(result_dd, model)
    if not (config.messagesoff):
        logger.info(
            """Simulation is finished. All files are closed. Please control the related files.
                    -----------------------------------------------------------------"""
        )
    return input_obj, model, result


if __name__ == "__main__":
    import optparse

    parser = optparse.OptionParser()
    parser.add_option(
        "-i",
        "--inputfile",
        action="store",
        dest="inputfile",
        default=config.inputfilename,
        help="chose input file to override default in config.py",
    )
    (options, args) = parser.parse_args()

    # Import from config file
    inputfile = __import__(options.inputfile)
    if not (config.messagesoff):
        logger.info("inputfile is %s", options.inputfile)
    run_aestimo(inputfile)
