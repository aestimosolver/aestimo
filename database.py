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

    You should hAve received a copy of the GNU General Public License
    along with this program. See ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt .

    For the list of contributors, see ~/AUTHORS

 Description:  Database file. Using lists for database entries.
                Absolutely, it is not using a parser subroutine.
                Quick and dirty solution for the code.
 References:
  - GaAs,AlAs parameters:
    Properties of Semiconductor Alloys: Group-IV, IIIV and IIVI Semiconductors Sadao AdAchi?2009 John Wiley & Sons, Ltd.
    Basic Semiconductor Physics Second Edition,Prof. Chihiro Hamaguchi 2010 Springer 
  
"""

# MATERIAL PROPERTIES
# materialproperties| Material : m_e | m_hh | epsilonStatic | Eg | Bowing_param | m_e_alpha |  Luttinger Parameters γ1,2 & 3 |Elastic constants C11,12|Lattice constant a0| Deformation potentials ac,av & b| delta splitt off|
materialproperty = {
'GaAs':{
'm_e':0.067,
'm_hh':0.45,
'm_lh':0.087,
'epsilonStatic':12.90,
'Eg':1.42,
'Bowing_param':0.0,
'Band_offset':0.65,
'm_e_alpha':5.3782e18,
'GA1':6.8,
'GA2':1.9,
'GA3':2.73, 
'C11':11.879,
'C12':5.376,
'a0':5.6533,
'Ac':-7.17,
'Av':1.16,
'B':-1.7,
'delta':0.34
},
'AlAs':{
'm_e':0.1,
'm_hh':0.51,
'm_lh':0.18,
'epsilonStatic':10.06,
'Eg':2.980,
'Bowing_param':0.0,
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
'Bowing_param':0.0,
'Band_offset':0.63,
'm_e_alpha':0.0,
'GA1':0.0,
'GA2':0.0,
'GA3':0.0, 
'C11':0.0,
'C12':0.0,
'a0':0.0,
'Ac':0.0,
'Av':0.0,
'B':0.0,
'delta':0.0
},
'InP':{
'm_e':0.073,
'm_hh':0.46,
'm_lh':0.12,
'epsilonStatic':12.50,
'Eg':1.35,
'Bowing_param':0.0,
'Band_offset':0.38,
'm_e_alpha':0.0,
'GA1':0.0,
'GA2':0.0,
'GA3':0.0, 
'C11':0.0,
'C12':0.0,
'a0':0.0,
'Ac':0.0,
'Av':0.0,
'B':0.0,
'delta':0.0
},
'GaP':{
'm_e':0.82,
'm_hh':0.6,
'm_lh':0.6,
'epsilonStatic':11.1,
'Eg':2.261,
'Bowing_param':0.0,
'Band_offset':0.55,
'm_e_alpha':0.0,
'GA1':0.0,
'GA2':0.0,
'GA3':0.0, 
'C11':0.0,
'C12':0.0,
'a0':0.0,
'Ac':0.0,
'Av':0.0,
'B':0.0,
'delta':0.0
},
'AlP':{
'm_e':0.22,
'm_hh':0.63,
'm_lh':0.2,
'epsilonStatic':10.464,
'Eg':2.48,
'Bowing_param':0.0,
'Band_offset':0.55,
'm_e_alpha':0.0,
'GA1':0.0,
'GA2':0.0,
'GA3':0.0, 
'C11':0.0,
'C12':0.0,
'a0':0.0,
'Ac':0.0,
'Av':0.0,
'B':0.0,
'delta':0.0
}
}

# ALLOY PROPERTIES
# alloyproperties| Alloy : m_e_x=0 | m_e_b  | eps_x=0 | eps_b | Eg | Bowing_param | m_e_alpha
alloyproperty = {
'AlGaAs':{
'Bowing_param':0.37,
'Band_offset':0.67,
'm_e_alpha':5.3782e18,
'Material1':'AlAs',
'Material2':'GaAs'
},
'InGaAs':{
'Bowing_param':0.58,
'Band_offset':0.63,
'm_e_alpha':0.0,
'Material1':'InAs',
'Material2':'GaAs'
},
'InGaP':{
'Bowing_param':0.65,
'Band_offset':0.33,
'm_e_alpha':0.0,
'Material1':'InP',
'Material2':'GaP'
},
'AlInP':{
'Bowing_param':0.13,
'Band_offset':0.52,
'm_e_alpha':0.0,
'Material1':'AlP',
'Material2':'InP'
}
}



