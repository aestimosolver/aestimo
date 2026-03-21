#!/usr/bin/env python
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
This file is one method of running aestimo. Simply define the input file in 
the the related line and run this script. We could also run aestimo.py directly
to achieve the same effect. 

Alternatively, many of the example input files show how we can transform them
into scripts that can be run directly to perform the simulations. That approach
also allows us to tailor each simulation even more to our needs.
"""
#import matplotlib.pyplot as pl
#import numpy as np
#import time
#import sys 

import config

import aestimo
    

# Import from config file
inputfile = ""
aestimo.logger.info("inputfile is %s" %inputfile)

if __name__=="__main__":
    aestimo.run_aestimo(inputfile)


