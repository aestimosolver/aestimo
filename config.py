#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains aestimo's global configuration settings for aestimo.py
and aestimo_eh.py. It contains parameters for controlling the algorithms that are 
used to calculate the bandstructures. 

The 'inputfilename' variable defines the default input file used when aestimo.py
or main.py is run directly as a script. There are also parameters that define the
defaults for saving and presenting results; as well as for logging messages.
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
q = 1.602176e-19 #C
meV2J=1e-3*q #meV to Joules

# CONFIGURATION 

# Default Input File(s)
# -------------
#inputfilename = "sample_qw_barrierdope"
#inputfilename = "sample_qw_qwdope"
#inputfilename = "sample_moddop"
#inputfilename = "sample_qw_HarrisonCh3_3"
#inputfilename = "sample_qw_barrierdope_p"
#inputfilename = "sample_multi_qw_barrierdope_p"
#inputfilename = "sample_double_qw"
#inputfilename = "sample_qw_barrierdope_p_ingan"
#inputfilename = "sample_qw_barrierdope_ingaas"
#inputfilename = "sample_mqw_barrierdope_p_ingan"
#inputfilename = "sample_qw_barrierdope_p_cdzno"
#inputfilename = "sample_multi_qw_barrierdope_p_ingan"
#inputfilename = "sample_qw_wide_isbt"
#inputfilename = "sample_qw_barrierdope_p_InGaAsP"
#inputfilename = "sample_qw_barrierdope_p_AlGaInN"
#inputfilename = "sample_qw_barrierdope_p_AlGaInN_2"
inputfilename = "sample_pn"
#inputfilename = "sample_2qw_barrierdope_ingaas"
# Calculation
# -----------
# Aestimo
use_cython = True #provides a speed up for aestimo
# Shooting method parameters for Schrödinger Equation solution
delta_E = 0.5*meV2J #Energy step (Joules) for initial search. Initial delta_E is 1 meV. 
d_E = 1e-5*meV2J #Energy step (Joules) within Newton-Raphson method when improving the precision of the energy of a found level.
E_start = 0.0    #Energy to start shooting method from (if E_start = 0.0 uses minimum of energy of bandstructure)
Estate_convergence_test = 1e-9*meV2J
# FermiDirac
FD_d_E = 1e-9 #Initial and minimum Energy step (meV) for derivative calculation for Newton-Raphson method to find E_F
FD_convergence_test = 1e-6 #meV
np_d_E = 1.0 # Energy step (meV) for dispersion calculations
# Poisson Loop
"""damping:An adjustable parameter  (0 < damping < 1) is typically set to 0.5 at low carrier densities. With increasing
carrier densities, a smaller value of it is needed for rapid convergence."""
damping = 0.5    #averaging factor between iterations to smooth convergence.
Stern_damping=True#the extrapolated-convergence-factor method instead of the fixed-convergence-factor method
max_iterations=80 #maximum number of iterations.
convergence_test=1e-4 #convergence is reached when the ground state energy (meV) is stable to within this number between iterations.

# Aestimo_numpy_h
predic_correc=True#predictor corrector method
anti_crossing_length=0.0001 # the lower lenght limit to consider anti-crossing (nm), works with old versions
amort_wave_0=1.5#ratio of half well's width for wavefunction  to penetration into the the left adjacent barrier
amort_wave_1=1.5#ratio of half well's width for wavefunction to penetration into the the right adjacent barrier
strain =True # for aestimo_numpy_eh
piezo=False # directly calculationg the induced electric field,for old poisson solver, works with old versions
piezo1=True #indirectly using interface charges.
quantum_effect=True#temporary
#--------------
parameters=False

# Output Files
# ------------
output_directory = "output_"+inputfilename
parameters = True
electricfield_out = True
potential_out = True
sigma_out = True
probability_out = True
states_out = True
Drift_Diffusion_out=True

# Result Viewer
# -------------
resultviewer = True
wavefunction_scalefactor = 400 # scales wavefunctions when plotting QW diagrams
# Messages
# --------
messagesoff = True
logfile = "aestimo.log"

