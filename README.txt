AESTIMO 1-D SELF-CONSISTENT SCHRÖDINGER-POISSON SOLVER
Version 0.6
++++++++++++++++++++++++++++++++++++++++++++++++++++++

Overview
--------

Aestimo 1-D Self-consistent Schrödinger-Poisson Solver (simply aestimo) is a simple 1-dimensional (1-D) simulator for semiconductor heterostructures. Aestimo is started as a hobby at the beginning of 2012, and become a usable tool which can be used as a co-tool in an educational and scientific work.

Hope that it also works for you. Please do not hesitate to contact us in case of any bugs found.

The documentation material on this wiki is copyrighted (c) 2013. Reuse of the material on this wiki is permitted under either the GNU Free Documentation License or the CC-BY-SA License.

Current features
----------------

    Material and alloys: Si, GaAs, AlAs and AlGaAs,
    Band structure for electrons,
    Carrier concentrations for electrons,
    Electric field distribution,
    Electron wavefunctions.

Getting Started
---------------

See the examples subdirectory of the distribution. Also, detailed information can be found in "Using the Code" part of this document. Subscription to the aestimo-users mailing list is highly recommended for further support. For developers and people interested in latest development issues, there is an aestimo-devel mailing list.

License
-------

The software is released under the modified BSD license, which means that everyone is free to download, use, and modify the code without charge.

Download and Installation
-------------------------

The latest version of the program is available in zipped form from the website: https://bitbucket.org/sblisesivdin/aestimo/.

Prerequisites
-------------

You will need to have a recent version of Python installed on your computer. For this, please refer to Python Website, where binary packages for most platforms can be found. Additionally, you need libraries called numpy and pylab. Both can be obtained from the Scipy Website.

For Macintosh, Python is preinstalled and related libraries can be found at Pythonmac Directory.

Running the Code
----------------
Most of the code is written in Python, and thus is platform independent. After extracting the aestimo_x.y.zip file to a folder, user may point the files that are written below in the related folder. Here x.y is the version number.

main.py - The file that you need to run.
config.py - A simple configuration file. You must enter the input filename into this configuration file.
database.py - A database for materials properties.
aestimo.py - Main program.
sample-hetero.py - A sample input file.
readme.txt - A read me file as you noticed.

First of all, user must prepare or use an input file. This file must specified in config.py file. There are other options in config.py file like necessary output files and on/off options for result viewer and in-run messages. After specifiying an input file in config.py, user can run the aestimo easily with executing the command

./aestimo.py

or

./main.py

If the output file options are true in config.py file, results can be found in the outputs folder as:

sigma.dat - Areal charge density
efield.dat - Electric field
potn.dat - Electric potential
states.dat - Values of states in the related quantum region
firststate.dat - Probability distribution of first quantum state 
