Aestimo 1D - One dimensional Self-Consistent Schrödinger-Poisson Solver
======================================================
Version 3.0
-------------

Overview
--------

Aestimo 1D Self-consistent Schrödinger-Poisson Solver (simply Aestimo1D) is a simple one-dimensional (1-D) simulator for semiconductor heterostructures. It started as a hobby at the beginning of 2012 and has become a usable tool that can be used as a co-tool in educational and scientific work.

Hope that it also works for you. Please do not hesitate to contact us if you find any bugs.

Current features
----------------

  * Material and alloys: GaAs, AlAs , InAs, InP, AlP, GaP, GaN, AlN, InN, CdO, MgO, ZnO, AlGaAs, InGaAs, InGaP, AlInP, InGaN, AlGaN, AlInN, MgZnO, CdZnO, InGaAsP, AlGaInN
  * Band structure for gamma electrons and heavy, light, and split-off holes,
  * Effective-mass method for electrons and 3x3 k.p method for holes,
  * Carrier concentrations for gamma electrons and heavy, light, and split-off holes,
  * Electric field distribution,
  * Electron wavefunctions,
  * Non-parabolicity,
  * External electric field,
  * Strain for valance band calculations,

Getting Started
---------------

See the examples subdirectory of the distribution. For detailed information, see the "Using the Code" part of this document. For further support, it is highly recommended that you subscribe to the aestimo-users mailing list.

License
-------

Aestimo is Copyrighted by the (C) 2013-2022 AestimoSolver group. This software is distributed under the terms of the GNU General Public License v3 (see ~/the COPYING file or http://www.gnu.org/copyleft/gpl.txt). This means that everyone is free to use, change, and share the changes.

Sefer Bora Lisesivdin is the project's initiator, and large contributions have since been made by Robert J. Steed and Hamza Hebal. For the full list of contributors, visit [related webpage.](https://aestimosolver.github.io/authors.html)

Get Help
--------
Before asking any question, please visit http://www.aestimosolver.org to read many tutorials, which include many important examples. The same tutorials are included in your /tutorials folder.

To ask a question to other possible users, please send your question as a [new issue](https://github.com/aestimosolver/aestimo/issues/new/choose) on GitHub. Also, please look at [previous issues](https://github.com/aestimosolver/aestimo/issues) before creating a new one.

Download and Installation
-------------------------

The program's latest version is available in zipped form from the website: https://github.com/aestimosolver/aestimo.

Alternatively, aestimo can now be installed from PyPI via the command `pip install aestimo.`

Prerequisites
-------------

You will need to have a version of Python 3 installed on your computer. For this, please refer to Python Website, where binary packages for most platforms can be found or search your distribution's package management system. Additionally, you need the following Python libraries - numpy, scipy, and matplotlib.

For Macintosh, Python is preinstalled, and related libraries can be found in the Pythonmac Directory.

Running the Code
----------------
Most of the code is written in Python and thus is platform-independent. After extracting the aestimo_x.y.zip file to a folder or installing it with the pip command, the user may need to add the aestimo.py file's PATH to the system's PATH variable. On Linux systems, this can be done by adding the following line at the end of the ~/.bashrc file.

    export PATH=/home/PATHTOAESTIMO:$PATH
    
Here, the user must know the real path instead of `/home/PATHTOAESTIMO`

After adding the PATH information, the user can run `aestimo.py` in any folder. Please visit the examples folder and run.

    aestimo.py -i sample_1qw_barrierdope_ingaas.py

If the user wants to see the results drawn after the end of the simulation, the -d argument can be used:

    aestimo.py -d -i sample_1qw_barrierdope_ingaas.py

In addition to running, `config.py`file also has important parameters for the simulation.

## Outputs

The results will be saved to a new folder with the name of the input file. After a successful simulation, you can find 5 different files in this directory. Each file has some data rows, which are explained below.

efield.dat
----------

Row 1: Position (m), Row 2: Electric Field strength (V/m)

potn.dat
--------

Row 1: Position (m), Row 2:[V_cb + V_p] (eV)

sigma.dat
---------

Row 1: Position (m), Row 2: Sigma (C/m^3)

states.dat
----------

Row 1: Number of states (number), Row 2: Carrier density of state (m^-2), Row 3: Energy of state (meV), Row 4: Effective mass of state (kg)

wavefunctions.dat
-----------------

Row 1: Position (m), Row 2: Psi^2
