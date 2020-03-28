AESTIMO 1-D SELF-CONSISTENT SCHRÖDINGER-POISSON SOLVER
======================================================
Version 2.0
-------------

Overview
--------

Aestimo 1-D Self-consistent Schrödinger-Poisson Solver (simply aestimo) is a simple 1-dimensional (1-D) simulator for semiconductor heterostructures. Aestimo was started as a hobby at the beginning of 2012, and become a usable tool which can be used as a co-tool in an educational and scientific work.

Hope that it also works for you. Please do not hesitate to contact us in case of any bugs found.

Current features
----------------

  * Material and alloys: GaAs, AlAs , InAs, InP, AlP, GaP, GaN, AlN, InN, CdO, MgO, ZnO, AlGaAs, InGaAs, InGaP, AlInP, InGaN, AlGaN, AlInN, MgZnO, CdZnO, InGaAsP, AlGaInN
  * Band structure for gamma electrons and heavy, light and split-off holes,
  * Effective-mass method for electrons and 3x3 k.p method for holes,
  * Carrier concentrations for gamma electrons and heavy, light and split-off holes,
  * Electric field distribution,
  * Electron wavefunctions,
  * Non-parabolicity,
  * External electric field,
  * Strain for valance band calculations,

Getting Started
---------------

See the examples subdirectory of the distribution. Also, detailed information can be found in "Using the Code" part of this document. Subscription to the aestimo-users mailing list is highly recommended for further support. For developers and people interested in latest development issues, there is an aestimo-devel mailing list.

License
-------

Aestimo is Copyrighted by (C) 2013-2020 AestimoSolver group. This software is distributed under the terms of the GNU General Public License v3, see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt . This means that everyone is free to use, change, share and share the changes.

Sefer Bora Lisesivdin is the initiator of the project, large contributions have since been made by Robert J. Steed and Hamza Hebal. For the full list of contributors, see ~/AUTHORS.

Get Help
--------
Before asking any question, please visit http://www.aestimosolver.org to read many tutorials which includes many important examples. Same tutorials are included in your /doc folder.

To ask a question to other possible users please send your question to email address: aestimo-users@googlegroups.com

Download and Installation
-------------------------

The latest version of the program is available in zipped form from the website: https://github.com/aestimosolver/aestimo.

Alternatively, aestimo can now be installed from PyPI via the command `pip install aestimo`.

(Note that if numpy is not installed before you install aestimo then there may be an compilation error for the cython extension but it seems that the extension gets compiled anyway, so that the error can be ignored.)

Prerequisites
-------------

You will need to have a version of Python 3 installed on your computer. For this, please refer to Python Website, where binary packages for most platforms can be found or search your distribution's package management system. Additionally, you need the following python libraries - numpy, scipy and matplotlib.

For Macintosh, Python is preinstalled and related libraries can be found at Pythonmac Directory.

Running the Code
----------------
Most of the code is written in Python, and thus is platform independent. After extracting the aestimo_x.y.zip file to a folder, user may point the files that are written below in the related folder. Here x.y is the version number.

  * main.py - The file that you need to run.
  * config.py - A simple configuration file. You must enter the input filename into this configuration file.
  * database.py - A database for materials properties.
  * aestimo.py - Main program which uses the Numpy library. Use this one for your conduction band calculations and gamma valley electrons.
  * aestimo_eh.py - Calculator for valence band calculations and holes.
  * VBHM.py - A class file for 3x3 k.p method.
  * sample-X.py - Some samples files (X) are included in the package with prefix "sample-".
  * main_iterating.py - A script for simulating a design several times while varying a parameter over a range of values.
  * README - A readme file as you noticed.
  * README_OUTPUTS - A readme about the structure of output files.
  * COPYING - License of the software.
  * AUTHORS - List of the committers.
  * /outputs - Output folder.
  * /outputs_eh - Output folder for aestimo_eh.

First of all, user must prepare or use an input file. This file must specified in `config.py` file. There are other options in `config.py` file like necessary output files and on/off options for result viewer and in-run messages. After specifiying an input file in `config.py`, user can run the aestimo easily with executing the command

    ./aestimo.py

for conduction band calculations. For valence band calculations, aestimo uses a 3x3 k.p model which includes strain. After editing `config.py` for input file, execute the command

    ./aestimo_eh.py

For simulating a design several times while varying a parameter over a range of values, edit the `main_iterating.py` file for your needs, and then execute it as

    ./main_iterating.py

If the output file options are true in ``config.py`` file, results can be found in the outputs folder. 

## Outputs Folder

In this directory, you can find 5 different files after a successful simulation. Each file have some data rows which are explained below.
Creation of the files can be controlled in config.py file.

efield.dat
----------

Row 1: Position (m), Row 2: Electric Field strength (V/m)

potn.dat
--------

Row 1: Position (m), Row 2:[V_cb + V_p] (J)

sigma.dat
---------

Row 1: Position (m), Row 2: Sigma (e/m^2)

states.dat
----------

Row 1: Number of state (number), Row 2: Carrier density of state (m^-2), Row 3: Energy of state (meV), Row 4: Effective mass of state (kg)

wavefunctions.dat
-----------------

Row 1: Position (m), Row 2: Psi^2
