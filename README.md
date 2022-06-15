Aestimo 1D - One dimensional Self-Consistent Schrödinger-Poisson Solver
======================================================
Version 3.0
-------------

Overview
--------

Aestimo 1D Self-consistent Schrödinger-Poisson Solver (simply Aestimo1D) is a simple 1-dimensional (1-D) simulator for semiconductor heterostructures. Aestimo1D was started as a hobby at the beginning of 2012, and become a usable tool which can be used as a co-tool in an educational and scientific work.

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

Aestimo is Copyrighted by (C) 2013-2022 AestimoSolver group. This software is distributed under the terms of the GNU General Public License v3, see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt . This means that everyone is free to use, change, share and share the changes.

Sefer Bora Lisesivdin is the initiator of the project, large contributions have since been made by Robert J. Steed and Hamza Hebal. For the full list of contributors, visit [related webpage.](https://www.aestimosolver.org/authors.html)

Get Help
--------
Before asking any question, please visit http://www.aestimosolver.org to read many tutorials which includes many important examples. Same tutorials are included in your /tutorials folder.

To ask a question to other possible users please send your question as an [new issue](https://github.com/aestimosolver/aestimo/issues/new/choose) on GitHub. Also, please look [previous issues](https://github.com/aestimosolver/aestimo/issues) before creating a new one.

Download and Installation
-------------------------

The latest version of the program is available in zipped form from the website: https://github.com/aestimosolver/aestimo.

Alternatively, aestimo can now be installed from PyPI via the command `pip install aestimo`.

Prerequisites
-------------

You will need to have a version of Python 3 installed on your computer. For this, please refer to Python Website, where binary packages for most platforms can be found or search your distribution's package management system. Additionally, you need the following python libraries - numpy, scipy and matplotlib.

For Macintosh, Python is preinstalled and related libraries can be found at Pythonmac Directory.

Running the Code
----------------
Most of the code is written in Python, and thus is platform independent. After extracting the aestimo_x.y.zip file to a folder, or installing with pip command, user may need add the PATH of aestimo.py file to system's PATH variable. On Linux systems this can be done by adding the following line at the end of the ~/.bashrc file.

    export PATH=/home/PATHTOAESTIMO:$PATH
    
Here, user must know the real path instead of `/home/PATHTOAESTIMO`.

After adding the PATH information, user can run `aestimo.py` in any folder. Please visit examples folder and run

    aestimo.py -i sample_1qw_barrierdope_ingaas.py

If the use want to see the results drawn after the end of simulation, -d argument can be used:

    aestimo.py -d -i sample_1qw_barrierdope_ingaas.py

In addition to running, `config.py`file also has important parameters for the simulation.

## Outputs

The results will be saved to a new folder with the name of input file. In this directory, you can find 5 different files after a successful simulation. Each file have some data rows which are explained below.

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
