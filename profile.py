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
Simple script for profiling aestimo using the cProfile module.
"""

import cProfile
command = """import main"""
cProfile.runctx( command, globals(),locals(),filename="aestimo_numpy-0.8.3.profile")

#command = """import aestimo"""
#cProfile.runctx( command, globals(), locals(), filename="aestimo_t8.profile" )

#command = """import main"""
#cProfile.runctx( command, globals(), locals(), filename="aestimo_numpy.profile" )
