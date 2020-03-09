# Aestimo 1D
Aestimo is a one-dimensional (1D) self-consistent Schrödinger-Poisson solver for semiconductor heterostructures. Aestimo is started as a hobby at the beginning of 2012, and become an usable tool which can be used as a co-tool in an educational and/or scientific work.

Aestimo is a [GPL](http://www.gnu.org/copyleft/gpl.txt) software, written by [several authors](AUTHORS.md). The code is hosted on [Github](http://github.com).

Please,
* Report bugs or [browse opened issues](https://github.com/aestimosolver/aestimo/issues).
* [Browse the code](https://github.com/aestimosolver/aestimo).
* [Contribute!](DEVSGUIDE.md)
* or socialize with the community using [FB](https://www.facebook.com/aestimosolver).

# Latest Aestimo

Latest version of Aestimo is v.1.2 (November 6th, 2017).

<a href="https://github.com/sblisesivdin/aestimo/releases/download/v.1.2.0/aestimo_v.1.2.0.zip"><img src="files/download.png"/></a>

Please read [User's Guide](USERSGUIDE.md) before downloading or visit [release notes](CHANGELOG.md) page for more information.

## Cite Information
We will appreciate if any scientific work done using Aestimo v.1.2.0 will contain an acknowledgment and the following reference:

* **R. Steed, H. Hebal, and S. B. Lisesivdin. (2017, November 6). sblisesivdin/aestimo: Version v.1.2 (Version v.1.2.0). Zenodo. http://doi.org/10.5281/zenodo.1042657**

# News

## Version v.1.2.0 is released!

*2017/11/06 19:37 · Dr. S.Bora Lisesivdin*

Nearly one year passed and the new version (v.1.2.0) of Aestimo is ready! Thanks to Robert Steed and Hamza Hebal for their efforts to make Aestimo better. We hope that you enjoy.

Some new features of Aestimo:
* Quaternary alloys (type A_{x}B_{1-x}C_{y}D_{1-y})
* Added an improved model for modeling conduction band intersubband transitions
* Added a periodic boundary condition for the Electric field for modeling repeating structures.
* Small changes were made like to get Aestimo to work on python3.4, example files were renamed to be more pythonic etc ...
* Code is moved to GitHub.

Here are the [release notes](CHANGELOG.md).

Version 1.2 can be downloaded from this [link](https://github.com/sblisesivdin/aestimo/releases/download/v.1.2.0/aestimo_v.1.2.0.zip). Aestimo is GPL licensed software. Any improvements, bug reports or feedback are welcome.

## Moving to GitHub!

*2017/04/04 22:47 - Dr. S.Bora Lisesivdin*

Today, we moved all our code and issue management from BitBucket to GitHub. You can find the code at github.com/sblisesivdin/aestimo. Now, with this move, we want to open our code to the large user base of GitHub and expect new committers as many as we can find. This move also makes easier for everyone to contribute to aestimo.

Not a single line of code is lost, just moved to another place. This is how we try to make this change for everyone as comfortable as possible.

## Version 1.1 is released!

*2016/11/11 09:50 - Dr. Robert Steed* 

The Aestimo Team is proud to release version 1.1.0 of the Aestimo 1D Self-consistent Schrödinger-Poisson Solver. This release includes some improvements to the intersubband_optical_transitions module. Figure instances are now returned from plotting functions to enable custom alteration and saving of figures. Docstrings for the modules have also been improved for those that explore the codebase using an interactive terminal like ipython. Here are the release notes.

Version 1.1 can be downloaded from https://github.com/sblisesivdin/aestimo/archive/master.zip. Aestimo is GPL licensed software. Any improvements, bug reports or feedback are welcome.

## Aestimo is on PyPI!

*2016/11/11 09:50 - Dr. Robert Steed*

Aestimo has been uploaded to PyPI and can now be installed using a simple `pip install aestimo`. This should also compile the cython extension automatically. You can find it at https://pypi.python.org/pypi/aestimo.

Imagine - you can now distribute your aestimo input script and simply tell someone to do `pip install aestimo` before running the script! Particularly since Python is now a popular language for science, this is no longer something esoteric to suggest someone to do. Well, one can dream!

## New domain name aestimosolver.org

*2015/04/12 10:52 - Dr. S. Bora Lisesivdin*

We are happy to announce our new domain name http://www.aestimosolver.org/. Old domain http://aestimo.ndct.org will be changed to a forward page in a near future.

## New tutorial on optical absorption

*2014/11/22 15:30 - Dr. Robert Steed*

We have uploaded a new tutorial about modelling the optical absorption of intersubband transitions in quantum wells. There is a module in aestimo specifically for these calculations (intersubband_optical_transitions.py) since the intersubband transitions can be surprisingly shifted from the frequencies of their underlying energy levels due to interactions between the electrons occupying the quantum well.

The tutorial is released as an ipython notebook (tutorial5) and can also be found in the doc directory of the aestimo package.

## Version 1.0 and the "elementary" version are released!

*2014/08/29 16:12 - Dr. S. Bora Lisesivdin*

Aestimo Team is proud to release the version 1.0 of Aestimo 1D Self-consistent Schrödinger-Poisson Solver. This version includes many bugfixes, new organization of the structure, wurtzite material using, lots of change in intersubband optical transitions module. Dielectric constants are now handled more accurately. Many tutorials are added in the form of an ipython notebooks. Code is heavily modified and stabilized with more than 60 commits. This version is the most feature rich and stable version so ever. For more please visit release notes. Also another aestimo version, which is called “elementary” is released. This version is easy to understand in terms of coding, and therefore it can be used in undergraduate /graduate courses and forking new solvers. There will be no any other release for the “elementary” version.

The version 1.0 and “elementary” can be downloaded from https://github.com/sblisesivdin/aestimo/archive/master.zip Aestimo is GPL licensed software. We are always open to your contributions. Please use it, issue a bug, help to write use a better user's manual or commit to source!

## New Tutorial

*2014/06/09 19:04 - Dr. S. Bora Lisesivdin*

Developer Robert Steed uploaded a new tutorial about the latest state of aestimo solver and general usage. It includes information about all solvers available, all main files, modelling a structure. It includes many figures and it will answer most of your questions.

## Aestimo v.0.9 is released

*2013/11/10 11:54 - Dr. S. Bora Lisesivdin*

Aestimo Team is proud to release the version 0.9 of Aestimo 1D Self-consistent Schrödinger-Poisson Solver. This version includes many bugfixes, speed improvements, cython code additions, rewritten VBMAT-V part to use numpy better, merging conduction and valance band calculations and more. Code is heavily modified and stabilized.

The version 0.9 can be downloaded from https://github.com/sblisesivdin/aestimo/archive/master.zip Aestimo is GPL licensed software. We are always open to your contributions. Please use it, issue a bug, help to write use a better user's manual or commit to source!

## Version 0.8 is released

*2013/07/07 11:35 - Dr. S. Bora Lisesivdin*
Aestimo Team is proud to release the version 0.8 of Aestimo 1D Self-consistent Schrödinger-Poisson Solver. This version includes many new features, bugfixes and small corrections. The most important feature, which is added to Aestimo recently, is the implementation of strain included valence band calculation with 3×3 k.p model. Also, Numpy version is restructured, input file structure and sample inputs are changed and non-parabolicity of conduction band is implemented (Numpy version only). In addition, database is changed to a more clear-understable structure and exchange interaction potential is implemented (Numpy version only). Logging with timers, some customizations in config and a possibility of looping the simulation over a parameter are also the added to new version. Aestimo can work now with new materials InAs, InP, AlP, GaP and new alloys InGaAs, InGaP, AlInP in addition to GaAs and AlGaAs.

The version 0.8 can be downloaded from https://github.com/sblisesivdin/aestimo/archive/master.zip Aestimo is GPL licensed software. We are always open to your contributions. Please use it, issue a bug, help to write use a better user's manual or commit to source!

## New site for Aestimo

*2013/06/28 17:21 - Dr. S. Bora Lisesivdin*

Because of the Bitbucket wiki is suffers from lack of features, now we are moving our wiki to a new domain http://aestimo.ndct.org/. The hosting and subdomain is maintained by one of our committer. With the new site, we will be able to include figures and equations in our website, and can make this site usable for all documentation needs. This site will have everything a user need to learn, download and use the aestimo.

The documentation material on this wiki is copyrighted © 2013. Reuse of the material on this wiki is permitted under GNU Free Documentation License 1.3
