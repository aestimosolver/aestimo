=========================
Download and Installation
=========================

The latest version of the program is available in zipped form from the website: https://github.com/aestimosolver/aestimo/.

Prerequisites
=============

You will need to have a recent version of Python installed on your computer. For this, please refer to Python Website, where binary packages for most platforms can be found. Additionally, you need libraries called numpy and pylab. Both can be obtained from the Scipy Website.

For Macintosh, Python is preinstalled and related libraries can be found at Pythonmac Directory.

Running the Code
================

Most of the code is written in Python, and thus is platform independent. After extracting the aestimo_x.y.zip file to a folder, user may point the files that are written below in the related folder. Here x.y is the version number.

    main.py - The file that you need to run. 
    config.py - A simple configuration file. You must enter the input filename into this configuration file. database.py - A database for materials properties. 
    aestimo.py - Main program for conduction band calculations and gamma valley electrons. Its code is simple to understand. 
    aestimo_numpy.py - Main program which uses the Numpy library. Use this one for your conduction band calculations and gamma valley electrons. 
    aestimo_numpy_h.py - Calculator for valence band calculations and holes. 
    VBHM.py - A class file for 3x3 k.p method. 
    sample-X.py - Some samples files (X) are included in the package with prefix "sample-". 
    main_iterating.py - A script for simulating a design several times while varying a parameter over a range of values. 
    README.md - A readme file as you noticed. 
    COPYING.md - License of the software. 
    AUTHORS.md - List of the committers. 

First of all, user must prepare or use an input file. This file must specified in config.py file. There are other options in config.py file like necessary output files and on/off options for result viewer and in-run messages. After specifiying an input file in config.py, user can run the aestimo easily with executing the command

    ./aestimo.py
or

    ./aestimo_numpy.py

for conduction band calculations. For valence band calculations, aestimo uses a 3×3 k.p model which includes strain. After editing config.py for input file, execute the command

    ./aestimo_numpy_eh.py

For simulating a design several times while varying a parameter over a range of values, edit the main_iterating.py file for your needs, and then execute it as

    ./main_iterating.py

If the output file options are true in config.py file, results can be found in the outputs folder as:

    sigma.dat - Areal charge density 
    efield.dat - Electric field 
    potn.dat - Electric potential 
    states.dat - Values of states in the related quantum region 
    firststate.dat - Probability distribution of first quantum state
