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
    Vurgaftman et al., J. Appl. Phys. 89 (11), 5815 (2001)
  - Si parameters:
  
"""

# MATERIAL PROPERTIES
# materialproperties| Material : cb_mass | vb_mass | epsilonStatic | Eg-bagil | V_CB | cb_mass_alpha
materialproperty = {'Si':   [0.156, 0.537, 11.7, 0.0, 0.0, 0.0],
                    'GaAs': [0.067, 0.500, 12.90, 0.0, 0.67, 5.3782e18],
                    'AlAs':  [0.15, 0.500, 10.06, 1.247, 0.67, 0.0]
                    }

# ALLOY PROPERTIES
# alloyproperties| Alloy : cb_mass_x=0 | cb_mass_b  | eps_x=0 | eps_b | Eg-bagil | V_CB | cb_mass_alpha
alloyproperty = {'AlGaAs':  [0.067, 0.083, 12.90, -2.84, 1.247, 0.67, 5.3782e18]
                }



