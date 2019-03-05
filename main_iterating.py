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
This script shows how we can simulate a design several times while varying a
parameter over a range of values.
"""
import matplotlib.pyplot as pl
import numpy as np
import os

# aestimo modules
import aestimo
import config
inputfile = __import__(config.inputfilename) 
import database

### Example: Parameter to loop over.
#thickness of first layer of structure
thicknesses = [8,12,14,16,20] #nm
        
# Initialise structure class
model = aestimo.StructureFrom(inputfile,database)

# Looping over the parameter
"""In order to write the code correctly, it is necessary to have 
understood the Structure class and it's attributes. Some of the 
input file variables and the class's attribute use different names
but they are normally easy to match up. Alternatively, we could vary the 
(runtime values of the) variables within inputfile module's namespace
and then create a fresh Structure class instance.

Equally, we can vary the values in the database module if we want to."""

results = []
output_directory = config.output_directory+'-numpy' # will be our local copy of the original value

for thickness in thicknesses:
    model.material[0][0] = thickness
    #other examples -
    #model.Fapp = ... #equivalent to Fapplied
    #model.subnumber_e = ...
    #model.dx = gridfactor*1e-9
    #database.materialproperty['GaAs']['epsilonStatic']= ...
    
    model.create_structure_arrays() # update the instance's internals
    
    # Perform the calculation
    result= aestimo.Poisson_Schrodinger(model)
    
    results.append(result) #all the results can be stored for further analysis. 
    
    # Set output directory 
    # aestimo_numpy reads the output directory from the config module, so
    config.output_directory = os.path.join(output_directory,'dz0_%dnm' %thickness)

    # Write the simulation results in files
    aestimo.save_and_plot(result,model)
    
    #Plot QW representation
    #aestimo.QWplot(result)#,figno=None) # an alternative to save_and_plot function
                                            # which only plots the QW diagram and doesn't
                                            # save anything.

print("Simulation is finished. All files are closed.")
print("Please control the related files.")
