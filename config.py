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

 Description:  It is a simple config file, however there will
                be multi-input file execution feature in the future.
  
"""
# CONFIGURATION 

# Input File(s)
# -------------
inputfilename = "sample-qw-barrierdope"
#inputfilename = "sample-qw-qwdope"
#inputfilename = "sample-moddop"
#inputfilename = "sample-qw-HarrisonCh3_3"

# Output Files
# ------------
electricfield_out = True
potential_out = True
sigma_out = True
probability_out = True
states_out = True

# Result Viewer
# -------------
resultviewer = True

# Messages
# --------
messagesoff = False
logfile = 'aestimo.log'
