#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Aestimo 1D Schrodinger-Poisson Solver
 Version v.0.8
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

 Description:  This is the main file.
"""
import matplotlib.pyplot as pl
import numpy as np
#import sys 

import config

if False:
    import aestimo_numpy as aestimo
    import database
    
    # Import from config file
    inputfile = __import__(config.inputfilename)
    
    # Initialise structure class
    model = aestimo.StructureFrom(inputfile,database)
    
    # Perform the calculation
    result= aestimo.Poisson_Schrodinger(model)
    
    # Write the simulation results in files
    aestimo.save_and_plot(result,model)

else:
    import aestimo

    
print "Simulation is finished. All files are closed."
print "Please control the related files."
