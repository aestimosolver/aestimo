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

 Description:  Database file. Using lists for database entries.
                Absolutely, it is not using a parser subroutine.
                Quick and dirty solution for the code.
 References:
  - GaAs,AlAs parameters:
    Properties of Semiconductor Alloys: Group-IV, III-V and II-VI Semiconductors Sadao AdAchi?2009 John Wiley & Sons, Ltd.
    Basic Semiconductor Physics Second Edition,Prof. Chihiro Hamaguchi 2010 Springer
    Physics of Optoelectronic Devices ,S-L.CHUANG ,1995 by John Wiley & Sons. Inc
  
"""

# MATERIAL PROPERTIES
# materialproperties| Material : m_e | m_hh | epsilonStatic | Eg | Bowing_param | m_e_alpha |  Luttinger Parameters γ1,2 & 3 |Elastic constants C11,12|Lattice constant a0| Deformation potentials ac,av & b| delta splitt off|
materialproperty = {
'GaAs':{
'm_e':0.067, #conduction band effective mass (relative to electron mass)
'm_hh':0.45, #heavy hole band effective mass (used by aestimo_numpy_h)
'm_lh':0.087, #light hole band effective mass (used by aetsimo_numpy_h)
'epsilonStatic':12.90, #dielectric constant
'Eg':1.519,#1.42 # (ev) band gap
'Ep':28.8, # (eV) k.p matrix element (used for non-parabolicity calculation (Vurgaftman2001)
'F':-1.94, # Kane parameter (used for non-parabolicity calculation (Vurgaftman2001)
'Band_offset':0.65, # conduction band/valence band offset ratio for GaAs - AlGaAs heterojunctions
'm_e_alpha':5.3782e18, # conduction band non-parabolicity variable for linear relation (Nelson approach)
# Valence band constants 
'delta':0.341, # (eV) Spin split-off energy gap
# below used by aestimo_numpy_h
'GA1':6.8, #luttinger parameter
'GA2':1.9, #luttinger parameter
'GA3':2.73, #luttinger parameter
'C11':11.879, # (GPa) Elastic Constants
'C12':5.376, # (GPa) Elastic Constants
'a0':5.6533, # (A)Lattice constant
'Ac':-7.17, # (eV) deformation potentials (Van de Walle formalism)
'Av':1.16, # (eV) deformation potentials (Van de Walle formalism)
'B':-1.7, # (eV) shear deformation potential (Van de Walle formalism)
},
'AlAs':{
'm_e':0.15,
'm_hh':0.51,
'm_lh':0.18,
'epsilonStatic':10.06,
'Eg':3.099,#2.980,
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
'Ac':-5.64,
'Av':2.47,
'B':-1.5,
'delta':0.28
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
'a0':6.0584,
'Ac':-5.08,
'Av':1.0,
'B':-1.8,
'delta':0.38
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
'Ac':-5.04,
'Av':1.27,
'B':-1.7,
'delta':0.11
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
'Ac':-7.14,
'Av':1.70,
'B':-1.8,
'delta':0.08
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
'Ac':-5.54,
'Av':3.15,
'B':-1.5,
'delta':0.04
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
'Material1':'AlAs',
'Material2':'GaAs'
},
'InGaAs':{
'Bowing_param':0.58,
'Band_offset':0.63,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'Material1':'InAs',
'Material2':'GaAs'
},
'InGaP':{
'Bowing_param':0.65,
'Band_offset':0.33,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'Material1':'InP',
'Material2':'GaP'
},
'AlInP':{
'Bowing_param':0.13,
'Band_offset':0.52,
'm_e_alpha':0.0,
'delta_bowing_param':0.0,
'Material1':'AlP',
'Material2':'InP'
}
}



