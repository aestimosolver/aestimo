#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    
    INFORMATION:
    Aestimo's database module. Contains a materialproperty dict containing 
    sub-dicts of values for each material and similar alloyproperty dict for the
    alloys of the materials. See the source for details on the required keys for
    each material or alloy.

    REFERENCES:
  - GaAs,AlAs parameters:
    Properties of Semiconductor Alloys: Group-IV, III-V, and II-VI Semiconductors Sadao AdAchi?2009 John Wiley & Sons, Ltd.
    Basic Semiconductor Physics Second Edition, Prof. Chihiro Hamaguchi 2010 Springer
    Physics of Optoelectronic Devices, S-L.CHUANG,1995 by John Wiley & Sons. Inc
  
"""

# MATERIAL PROPERTIES
# materialproperties| Material : m_e | m_hh | epsilonStatic | Eg | Bowing_param | m_e_alpha |  Luttinger Parameters γ1,2 & 3 |Elastic constants C11,12|Lattice constant a0| Deformation potentials ac,av & b| delta splitt off|
materialproperty = {
'GaAs':{
'm_e':0.067, #conduction band effective mass (relative to electron mass)
'm_hh':0.45, #heavy hole band effective mass (used by aestimo_numpy_h)
'm_lh':0.087, #light hole band effective mass (used by aetsimo_numpy_h)
'epsilonStatic':12.90, #dielectric constant
'Eg':1.4223,#1.42 # (ev) band gap
'Ep':28.8, # (eV) k.p matrix element (used for non-parabolicity calculation (Vurgaftman2001)
'F':-1.94, # Kane parameter (used for non-parabolicity calculation (Vurgaftman2001)
'Band_offset':0.65, # conduction band/valence band offset ratio for GaAs - AlGaAs heterojunctions
'm_e_alpha':5.3782e18, # conduction band non-parabolicity variable for linear relation (Nelson approach)
# Valence band constants 
'delta':0.28, # (eV) Spin split-off energy gap
# below used by aestimo_numpy_h
'GA1':6.8, #luttinger parameter
'GA2':1.9, #luttinger parameter
'GA3':2.73, #luttinger parameter
'C11':11.879, # (GPa) Elastic Constants
'C12':5.376, # (GPa) Elastic Constants
'a0':5.6533, # (A)Lattice constant
'a0_sub':5.6533,
'Ac':-7.17, # (eV) deformation potentials (Van de Walle formalism)
'Av':1.16, # (eV) deformation potentials (Van de Walle formalism)
'B':-1.7, # (eV) shear deformation potential (Van de Walle formalism)
'TAUN0':0.1E-7,# Electron SRH lifetime
'TAUP0':0.1E-7,# Hole SRH lifetime
'mun0':0.1,# Electron Mobility in m2/V-s
'mup0':0.02,# Electron Mobility in m2/V-s
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,# Parameter in the calculation of the Field Dependant Mobility
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-6.92#Average Valence Band Energy or the absolute energy level
},
'AlAs':{
'm_e':0.15,
'm_hh':0.51,
'm_lh':0.18,
'epsilonStatic':10.06,
'Eg':3.0,#2.980,
'Ep':21.1,
'F':-0.48,
'Band_offset':0.53,
'm_e_alpha':0.0,
'GA1':3.45,
'GA2':0.68,
'GA3':1.29, 
'C11':11.879,
'C12':5.376,
'a0':5.66,
'a0_sub':5.6533,
'Ac':-5.64,
'Av':2.47,
'B':-1.5,
'delta':0.28,
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-7.49#Average Valence Band Energy or the absolute energy level
},
'InAs':{
'm_e':0.4,
'm_hh':0.26,
'm_lh':0.027,
'epsilonStatic':15.15,
'Eg':0.4,
'Ep':21.5,
'F':-2.90,
'Band_offset':0.63,
'm_e_alpha':0.0,
'GA1':20.4,
'GA2':8.3,
'GA3':9.1,
'C11':8.329,
'C12':4.526,
'a0':6.0583,
'a0_sub':5.6533,
'Ac':-5.08,
'Av':1.0,
'B':-1.8,
'delta':0.38,
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-6.67#Average Valence Band Energy or the absolute energy level
},
'InP':{
'm_e':0.073,
'm_hh':0.46,
'm_lh':0.12,
'epsilonStatic':12.50,
'Eg':1.35,
'Ep':20.7,
'F':-1.31,
'Band_offset':0.38,
'm_e_alpha':0.0,
'GA1':5.33,
'GA2':1.57,
'GA3':2.11,
'C11':8.329,
'C12':4.526,
'a0':5.8688,
'a0_sub':5.6533,
'Ac':-5.04,
'Av':1.27,
'B':-1.7,
'delta':0.11,
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-7.04#Average Valence Band Energy or the absolute energy level
},
'GaP':{
'm_e':0.82,
'm_hh':0.6,
'm_lh':0.6,
'epsilonStatic':11.1,
'Eg':2.261,
'Ep':31.4,
'F':-2.04,
'Band_offset':0.55,
'm_e_alpha':0.0,
'GA1':4.04,
'GA2':0.53,
'GA3':1.26,
'C11':14.05,
'C12':6.203,
'a0':5.4505,
'a0_sub':5.6533,
'Ac':-7.14,
'Av':1.70,
'B':-1.8,
'delta':0.08,
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-7.4#Average Valence Band Energy or the absolute energy level
},
'AlP':{
'm_e':0.22,
'm_hh':0.63,
'm_lh':0.2,
'epsilonStatic':10.464,
'Eg':2.48,
'Ep':17.7,
'F':-0.65,
'Band_offset':0.55,
'm_e_alpha':0.0,
'GA1':3.47,
'GA2':0.06,
'GA3':1.15,
'C11':15.0,
'C12':6.42,
'a0':5.4635,
'a0_sub':5.6533,
'Ac':-5.54,
'Av':3.15,
'B':-1.5,
'delta':0.04,
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-8.09#Average Valence Band Energy or the absolute energy level
},
'GaN':{ 
'm_e':0.21,#☻.22
'm_e_alpha':5.3782e18,
'm_hh':1.87,#1.1 1
'm_lh':0.14,#1.1
'm_so':0.5,#0.15
'epsilonStatic':9.5,
'Eg':3.5,# 3.44 Madelung (1991) 3.437
'Bowing_param':3.2,
'Band_offset':0.63,
'A1':-6.4,#-6.56 -0.91 5.65 -2.83 -3.13 -4.86
'A2':-0.5,
'A3':5.9,
'A4':-2.95,
'A5':-2.56,
'A6':-3.06,
'D1':-1.7,#-3.7 4.5 8.2 -4.1 -4.0 -5.5
'D2':6.3,
'D3':8.2,
'D4':-4.1,
'D5':-4,
'D6':-5.65,
'Ac':-4.60,
'a0_wz':3.189,#3.189 5.185
'C11':37.,
'C12':14.5,
'C13':10.3,
'C33':39.0,#40.4
'C44':9.5,
'C66':12.3,
'D15':-1.579e-12,#-1.579e-12      -1.7e-12
'D31':-1.253e-12,#-1.253e-12      -1.6e-12 
'D33':2.291e-12,# 2.291e-12     2.291e-12
'Psp':-0.029,#-0.029 -0.034 
'delta_so':0.015,#15.5mev
'delta_cr':0.022,#72.9mev
'a0_sub':3.189,
'TAUN0':0.1E-6,# Electron SRH lifetime
'TAUP0':0.1E-6,# Hole SRH lifetime
'mun0':0.1,# Electron Mobility in m2/V-s 0.1 0.2
'mup0':0.02,# Electron Mobility in m2/V-s 0.03 0.02
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2,# Parameter in the calculation of the Field Dependant Mobility
'BETAP':2,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,#3e5 Saturation Velocity of Electrons 1.5e5 2.5e5 2.7e5 2e5
'VSATP': 6e5, # Saturation Velocity of Holes
'AVb_E':-2.64#Average Valence Band Energy or the absolute energy level
},
'InN':{#now
'm_e':0.07,#0.11
'm_e_alpha':5.3782e18,
'm_hh':1.61,#1.61
'm_lh':0.11,#1.61
'm_so':0.5,#0.1
'epsilonStatic':15.3,
'Eg':0.608,#0.76 0.608
'Bowing_param':3.2,
'Band_offset':0.63,#
'A1':-9.09,#-9.28 -0.60 8.68 -4.34 -4.32 -6.08
'A2':-0.63,
'A3':8.46,
'A4':-4.23,
'A5':-4.36,
'A6':-6.34,
'D1':-1.76,#-3.7 4.5 8.2 -4.1 -4.0 -5.5
'D2':3.43,
'D3':5.19,
'D4':-2.595,
'D5':-2.33,
'D6':-5.5,
'Ac':-1.4,
'a0_wz':3.53,#3.548
'C11':22.3,
'C12':11.5,
'C13':9.2,
'C33':20.5,
'C44':9.5,
'C66':7.4,
'D15':-2.292e-12,#-2.292e-12 -1.1e-12
'D31':-3.147e-12,#-3.147e-12 -3.5e-12
'D33':6.201e-12,#6.201e-12
'Psp':-0.032,#-0.032 -0.042 
'delta_so':0.001,#5mev
'delta_cr':0.041,#40mev
'a0_sub':3.189,
'TAUN0':5.1e-8,# Electron SRH lifetime
'TAUP0':5.1e-8,# Hole SRH lifetime
'mun0':0.1,# Electron Mobility in m2/V-s 0.1 0.2
'mup0':0.02,# Electron Mobility in m2/V-s 0.03 0.02
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2,# Parameter in the calculation of the Field Dependant Mobility
'BETAP':1,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,#3e5 Saturation Velocity of Electrons 1.5e5
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-1.59#Average Valence Band Energy or the absolute energy level
},
'AlN':{
'm_e':0.22,
'm_e_alpha':5.3782e18,
'm_hh':1.2,
'm_lh':0.8,
'm_so':0.5,
'epsilonStatic':10.464,
'Eg':6.28,
'Bowing_param':0.0,
'Band_offset':0.55,
'A1':-3.95,
'A2':-0.27,
'A3':3.68,
'A4':-1.84,
'A5':-1.92,
'A6':-2.91,
'D1':-0.89,
'D2':4.27,
'D3':5.18,
'D4':-2.59,
'D5':-4,
'D6':3.4,
'Ac':-7.17,
'a0_wz':3.112,
'C11':39.6,
'C12':14,
'C13':10.8,
'C33':37.7,
'C44':11.6,
'D15':-2.0e-12,#−2.069e-12      −4.08e-12
'D31':-2.650e-12,#−2.298e-12    −2.65e-12
'D33':5.53e-12,#5.352e-12        5.53e-12
'Psp':-0.081,#-0.081 −0.09
'delta_so':0.019,#20.4mev
'delta_cr':-0.164,#−58.5mev
'a0_sub':3.189,
'TAUN0':0.1E-7,# Electron SRH lifetime
'TAUP0':0.1E-7,# Hole SRH lifetime
'mun0':0.1,# Electron Mobility in m2/V-s
'mup0':0.02,# Electron Mobility in m2/V-s
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':0.45,# Parameter in the calculation of the Field Dependant Mobility
'BETAP':0.45,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons 1.5e5
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-3.44#Average Valence Band Energy or the absolute energy level 1.5e5
},
'CdO':{
'm_e':0.12,
'm_e_alpha':5.3782e18,
'm_hh':1.2,
'm_lh':0.8,
'm_so':0.5,
'epsilonStatic':10.464,
'Eg':1.89,
'Bowing_param':0.0,
'Band_offset':0.65,
'A1':-3.78,
'A2':-0.44,
'A3':3.45,
'A4':-1.63,
'A5':-1.68,
'A6':-2.23,
'D1':-3.90,
'D2':-4.13,
'D3':-1.15,
'D4':1.22,
'D5':-1.53,
'D6':2.83,
'Ac':-1.4,
'a0_wz':3.45,#3.66
'C11':20.97,
'C12':12.11,
'C13':10.51,
'C33':21.09,
'C44':4.247,
'D15':-1.1e-12,
'D31':-1.1e-12,
'D33':0.22e-12,
'Psp':-0.099,
'delta_so':0.0126,
'delta_cr':0.0305,
'a0_sub':3.250,
'TAUN0':0.1E-7,# Electron SRH lifetime
'TAUP0':0.1E-7,# Hole SRH lifetime
'mun0':0.1,# Electron Mobility in m2/V-s
'mup0':0.02,# Electron Mobility in m2/V-s
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,# Parameter in the calculation of the Field Dependant Mobility
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-7.52#Average Valence Band Energy or the absolute energy level
},
'MgO':{
'm_e':0.24,
'm_e_alpha':5.3782e18,
'm_hh':1.2,
'm_lh':0.8,
'm_so':0.5,
'epsilonStatic':9.6,
'Eg':5.289,
'Bowing_param':0.1,
'Band_offset':0.65,
'A1':-3.78,
'A2':-0.44,
'A3':3.45,
'A4':-1.63,
'A5':-1.684,
'A6':-2.23,
'D1':-3.90,
'D2':-4.13,
'D3':1.15,
'D4':-1.22,
'D5':-1.53,
'D6':2.83,
'Ac':-7.17,
'a0_wz':3.199,
'C11':22.0,
'C12':9.42,
'C13':5.4,
'C33':21.6,
'C44':10.5,
'D15':-1.7e-12,
'D31':-7.887e-12,
'D33':0.34e-12,
'Psp':-0.068,
'delta_so':0.032,
'delta_cr':0.3172,
'a0_sub':3.250,
'TAUN0':0.1E-7,# Electron SRH lifetime
'TAUP0':0.1E-7,# Hole SRH lifetime
'mun0':0.1,# Electron Mobility in m2/V-s
'mup0':0.02,# Electron Mobility in m2/V-s
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,# Parameter in the calculation of the Field Dependant Mobility
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-7.52#Average Valence Band Energy or the absolute energy level
},
'ZnO':{
'm_e':0.24,
'm_e_alpha':5.3782e18,
'm_hh':1.2,
'm_lh':0.8,
'm_so':0.5,
'epsilonStatic':8.1,
'Eg':3.37,
'Bowing_param':0.0,
'Band_offset':0.65,
'A1':-3.78,
'A2':-0.44,
'A3':3.45,
'A4':-4.32,
'A5':-3.13,
'A6':-2.23,
'D1':-3.90,
'D2':-4.13,
'D3':-1.15,
'D4':1.22,
'D5':-1.53,
'D6':2.83,
'Ac':-6.05,
'a0_wz':3.250,
'C11':20.97,
'C12':12.11,
'C13':10.51,
'C33':21.09,
'C44':4.247,
'D15':-5e-12,
'D31':-5e-12,
'D33':1e-12,
'Psp':-0.05,
'delta_so':0.0126,#16mev
'delta_cr':0.0305,#43mev
'a0_sub':3.250,
'TAUN0':0.1E-7,# Electron SRH lifetime
'TAUP0':0.1E-7,# Hole SRH lifetime
'mun0':0.1,# Electron Mobility in m2/V-s
'mup0':0.02,# Electron Mobility in m2/V-s
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,# Parameter in the calculation of the Field Dependant Mobility
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-7.52#Average Valence Band Energy or the absolute energy level
},
'Si':{#http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
'm_e':1.18,
'm_hh':0.49,
'm_lh':0.16,
'epsilonStatic':11.7,
'Eg':1.12,
'Ep':21.1,#?
'F':-0.48,#?
'Band_offset':0.53,
'm_e_alpha':0.0,
'GA1':3.45,
'GA2':0.68,
'GA3':1.29, 
'C11':16.6,
'C12':6.4,
'C44':7.96,
'a0':5.431,
'a0_sub':5.4315,  # Just added as near value of a0 
'Ac':-5.64,
'Av':2.47,
'B':-1.5,
'delta':0.28,
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-7.03#Average Valence Band Energy or the absolute energy level
},
'Ge':{
'm_e':1.59,
'm_hh':0.33,
'm_lh':0.043,
'epsilonStatic':10.06,
'Eg':0.661,
'Ep':21.1,
'F':-0.48,
'Band_offset':0.53,
'm_e_alpha':0.0,
'GA1':3.45,
'GA2':0.68,
'GA3':1.29, 
'C11':12.60,
'C12':4.4,
'C44':6.77,
'a0':5.658,
'a0_sub':5.6585, # Just added as near value of a0 
'Ac':-5.64,
'Av':2.47,
'B':-1.5,
'delta':0.28,
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-6.35#Average Valence Band Energy or the absolute energy level
}
}

# ALLOY PROPERTIES
# alloyproperties| Alloy : m_e_x=0 | m_e_b  | eps_x=0 | eps_b | Eg | Bowing_param | m_e_alpha
alloyproperty = {
'AlGaAs':{
'Bowing_param':0.37,
'Band_offset':0.65,
'm_e_alpha':5.3782e18,
'delta_bowing_param':0.0,
'a0_sub':5.6533,
'Material1':'AlAs',
'Material2':'GaAs',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'InGaAs':{
'Bowing_param':0.58,
'Band_offset':0.63,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'a0_sub':5.6533,
'Material1':'InAs',
'Material2':'GaAs',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'InAsP':{
'Bowing_param':0.58,
'Band_offset':0.63,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'a0_sub':5.6533,
'Material1':'InAs',
'Material2':'InP',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'GaAsP':{
'Bowing_param':0.58,
'Band_offset':0.63,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'a0_sub':5.6533,
'Material1':'GaP',
'Material2':'GaAs',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'InGaP':{
'Bowing_param':0.65,
'Band_offset':0.33,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'a0_sub':5.6533,
'Material1':'InP',
'Material2':'GaP',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'AlInP':{
'Bowing_param':0.13,
'Band_offset':0.52,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'a0_sub':5.6533,
'Material1':'AlP',
'Material2':'InP',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'AlGaN':{
'Bowing_param':0.7,#0.7 1.3 1 1.3
'Band_offset':0.63,
'm_e_alpha':5.3782e18,
'a0_sub':3.189,
'c0_sub':4.982,
'Material1':'AlN',
'Material2':'GaN',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,#0.1 0.015
'mup0':0.1,#0.03 0.1
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'InGaN':{
'Bowing_param':1.4,#1.4 3.8 3
'Band_offset':0.63,#0.677 0.63
'm_e_alpha':0.0,
'a0_sub':3.189,
'c0_sub':4.982,
'Material1':'InN',
'Material2':'GaN',
'TAUN0':5.2e-8,
'TAUP0':5.2e-8,
'mun0':0.1,#0.1 0.015
'mup0':0.02,#0.03 0.1
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# 3e5Saturation Velocity of Electrons 1.5e5
'VSATP':6e5 #6e5 Saturation Velocity of Holes 1.5e5
},
'AlInN':{
'Bowing_param':3.6,#2.5
'Band_offset':0.677,
'm_e_alpha':0.0,
'a0_sub':3.189,
'c0_sub':4.982,
'Material1':'AlN',
'Material2':'InN',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'MgZnO':{
'Bowing_param':0.87,
'Band_offset':0.65,
'm_e_alpha':5.3782e18,
'a0_sub':3.250,
'c0_sub':4.982,
'Material1':'MgO',
'Material2':'ZnO',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'CdZnO':{
'Bowing_param':3.8,
'Band_offset':0.65,
'm_e_alpha':5.3782e18,
'a0_sub':3.250,
'c0_sub':4.982,
'Material1':'CdO',
'Material2':'ZnO',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
}
}
# ALLOY PROPERTIES
# alloyproperties| Alloy : m_e_x=0 | m_e_b  | eps_x=0 | eps_b | Eg | Bowing_param | m_e_alpha
alloyproperty4 = {
#AxB(1-x)CyD(1-y)=InxGa(1-x)AsyP(1-y)
'InGaAsP':{
'Bowing_param_ABC':0.58,
'Bowing_param_ABD':0.65,
'Bowing_param_ACD':0.58,
'Bowing_param_BCD':0.58,
'Band_offset':0.4,
'm_e_alpha':5.3782e18,
'delta_bowing_param':0.0,
'a0_sub':5.8688,
'Material1':'InAs',
'Material2':'GaAs',
'Material3':'InP',
'Material4':'GaP',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
},
'AlGaInN':{#BxCyD1-x-yA=AlxInyGa1-x-yN
'Bowing_param_ABC':3.6,#AlInN
'Bowing_param_ABD':1.3,#AlGaN
'Bowing_param_ACD':3.2,#InGaN
'Band_offset':0.7,
'm_e_alpha':5.3782e18,
'delta_bowing_param':0.0,
'a0_sub':3.189,
'c0_sub':4.982,
'Material1':'GaN',
'Material2':'InN',
'Material3':'AlN',
'TAUN0':0.1E-6,
'TAUP0':0.1E-6,
'mun0':0.15,
'mup0':0.1,
'Cn0':2.8e-31,# generation recombination model parameters [cm**6/s]
'Cp0':2.8e-32,# generation recombination model parameters [cm**6/s]
'BETAN':2.0,
'BETAP':1.0,# Parameter in the calculation of the Field Dependant Mobility
'VSATN':3e5,# Saturation Velocity of Electrons
'VSATP':6e5, # Saturation Velocity of Holes
'AVb_E':-2.1#Average Valence Band Energy or the absolute energy level
}
}



