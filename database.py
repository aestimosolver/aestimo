#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
  Program:      Aestimo 1D Schrodinger-Poisson Solver
  Description:  Database file. Using lists for database entries.
                Absolutely, it is not using a parser subroutine.
                Quick and dirty solution for the code.
  References:
  - GaAs,AlAs parameters:
    Vurgaftman et al., J. Appl. Phys. 89 (11), 5815 (2001)
  - Si parameters:
  
"""

# MATERIAL PROPERTIES
# materialproperties| Material : cb_mass | vb_mass | epsilonStatic | Eg-bagil | V_CB |
materialproperty = {'Si':   [0.156, 0.537, 11.7, 0.0, 0.0],
                    'GaAs': [0.067, 0.500, 12.83, 0.0, 0.67],
                    'AlAs':  [0.15, 0.500, 12.83, 1.247, 0.67]
                    }

# ALLOY PROPERTIES
# alloyproperties| Alloy : cb_mass_x=0 | cb_mass_b  | eps_x=0 | eps_b | Eg-bagil | V_CB |
alloyproperty = {'AlGaAs':  [0.067, 0.083, 12.90, -2.84, 1.247, 0.67]
                }



