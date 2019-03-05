#!/bin/env python
# -*- coding: utf-8 -*-
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Copyright (C) 2013-2018 Sefer Bora Lisesivdin and Aestimo group

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
This module compiles the cythonised version of the psi_at_inf function used
by the aestimo.py shooting method (but not by aestimo_eh.py). See 
psi_at_inf_cython.pyx

Compile the cythonised function with the command:
   python setup_cython.py build_ext --inplace
or on windows:
   python setup_cython.py build_ext --inplace --compiler=mingw32

aestimo.py will then automatically use this faster version as long as the 
config.py module contains `use_cython = True`.
"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("psi_at_inf_cython", ["psi_at_inf_cython.pyx"])]

setup(
  name = 'psi_at_inf',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
  include_dirs=[numpy.get_include()]
)
