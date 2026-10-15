# Release notes

## v.4.0.0

- Add the desktop GUI and device project workflows, including quantum-well, diode, solar-cell, LED and laser presets.
- Bundle GUI project presets and reference data in pip distributions and provide a writable installed workspace without replacing user-edited files.
- Provide `aestimo`, `aestimo-gui` and the Windows GUI launcher. Report the local version without a network request.
- Avoid creating output/log files during module imports, CLI help or version reporting. Route calculation logs to the selected output directory without accumulating duplicate handlers, and correctly pass the CLI `-d` plot option.
- Add a user guide configured for Read the Docs and consolidate package metadata in `pyproject.toml`.
- Update the user guide, project citation, contributor credits and release history; correct obsolete installation, input and convergence instructions.
- Expand the publications list with verified 2023–2026 papers citing Aestimo.
- Restore the project's red styling and original banner. Put pip installation on the home page and expose authors, release notes and development directly in the navigation; keep internal integration records in the repository.
- Expand the GUI guide and first-calculation walkthrough, documenting saved result history, workspace paths and GUI/CLI input differences. Correct state export mass headers to identify the existing relative-mass values.
- Size plotting helper arrays to the actual grid instead of allocating 100 million float64 entries (about 800 MB) for each plot call.
- Report unavailable solver selections and grid-limit failures as catchable errors. CLI calculation failures return nonzero status, and malformed inputs receive concise usage errors instead of reporting success or stopping a worker.
- Add automated installed-package checks on Linux/Python 3.9, 3.12 and 3.13 and Windows/Python 3.12. Correct a solar plot label that prevented GUI imports on Python versions before 3.12.
- Distinguish numerical comparisons, synthetic baselines and unverified source claims in reference-data reporting.
- Preserve existing default GUI projects during automatic setup and use UTF-8 for project and user-database JSON files.
- Reject unconverged Mode 10 transport states and require spatial total-current conservation; stabilize near-equilibrium edge-flux evaluation.

## v.3.0.0 (July 13, 2022)
This release simplified and cleaned up the code and its usage, following the 2.0 release by approximately two years. 
- New argument parser. Users can use -h for help, -v for version and latest package information.
- The old electron-only executable was removed; the combined electron/hole solver formerly named aestimo_eh.py became aestimo.py.
- StructureFrom() is redesigned.
- New argument -d. Without using of this arg, aestimo will only calculate and write results to the output folder. No figures are shown.
- Save log and all output at the same folder named input file
- Logging is always on, removed from config.
- Input files can be executed from any folder. Therefore, it is removed from the config file.
- Config file is simplified.
- Cython support is removed. Aestimo is a 1D solver. Therefore, simulations took only few minutes to finish. We do not need to run it with cython.

## v.2.0 (March 28, 2020)
This is a major release nearly after 2 years. We would like to thank Hamza Hebal for his efforts and contributions. There are many new features added in this release:
- New DD solver base on Gummel and newton map.
- Updating structures files with new options and settings.
- Updating database with variables for new DD solver.
- Adding the stern damping method as an option.
- Adding a contact setting as a boundary for the Poisson equation.
- Improve piezo- and spontaneous- charges handling independently of wb and bw variables.
- Improve quantum well counting in the structure list.
- Including the option of the free quantum region.
- Correction in the position of interface charges for piezo option.
- The Quantum density of carrier is working now, according to the Modified QDD model.
- Adding new structure files.
- Setting for predictor-corrector method.-Setting for QW calculation limits around the well region.
- Organizing existing functions for Schrödinger and strain setting in separate functions.
- Improved method for locating energy levels for each QW.
- Quaternary alloys (type B_{x}C_{y}D_{1-x-y}A).
- Version 2.0 works with Python 3.7.6 and previous python3 versions.

## v.1.2.1 (July 22nd, 2018)

This is a minor release that fixes some issues with using aestimo as a system package in python3 and makes some changes to the intersubband_optical_transitions.py. This version will also be available on PyPI for installation via pip.

Changes include:

- intersubband_optical_transitions module now has model handling wavelength dependent dielectric constants - reliability of the results is unproven though.
- intersubband_optical_transitions module has made some small changes to improve handling of absorbing background dielectric constant - reliability of the results is unproven though.
- Fixed issues related to using aestimo as a package in python3.6
- Removing figsize from aestimo and aestimo_eh plots so that it can be configured via matplotlib configuration.

## v.1.2 (November 6th, 2017)

*	Quaternary alloys (type A_{x}B_{1-x}C_{y}D_{1-y})
*	Added an improved model for modelling conduction band intersubband transitions
*	Added a periodic boundary condition for the Electric field for modelling repeating structures.
*	Small changes were made like to get aestimo to work on python3.4, example files were renamed to be more pythonic etc ...
* 	Code is moved to GitHub.

## v.1.1 (November 8th, 2016)

*    Many small changes and improvements to the intersubband\_optical\_transitions module.
*    Plotting functions now return figure instances; helpful for altering and then saving versions of the plots.
*    Module docstrings have been improved.

## v.1.0 (August 29th, 2014)

*    Many bugfixes,
*    MQW structures ready to be used after the addition of anti-crossing barrier length,
*    A big reorganisation of the Structure and StructureFrom classes. This should make things more self-explanatory for anyone who wants to do something special; normal users shouldn't be affected,
*    Adding a new input parameter 'meff_method' which selects the effective mass functions to use when considering non-parabolicity. This also allows us to not require the material properties in the database to contain unnecessary attributes if we don't want to consider non-parabolicity,
*    Adding functionality to use valence band input files as scripts,
*    Changing file permissions and adding a missing parameter to some of the sample files,
*    Added a utility function to each simulator which carries out the default steps to run a simulation, altered some of the sample files and main.py to use the new function, 
*    Removing aestimo.py cases in main.py and main_interating.py,
*    Get help information is added. Information about obsolete aestimo.py is removed,
*    Adding many tutorials in the form of an ipython notebooks,
*    Allowing aestimo to be used as a package,
*    Moving 'aestimo is starting' message to main section to make using solvers as modules more consistent,
*    A big change to the code to permit the multi-quantum wells calculation,
*    Wurtzite materials are now used in aestimo with,
*    Lots of changes to intersubband optical transitions module. Dielectric constants are now handled more accurately. Although frequency dependent dielectric constants are still not implemented for the matrix model (classical model should be able to handle them),
*    Adding module for calculating dielectric constant and optical absorption of conduction band intersubband transitions,
*    More changes to saving and loading of simulation results,
*    Adding function for loading output data,
*    Adding calculation of fermilevel for non-parabolic subbands,
*    Working on k.p non-parabolicity (Vurgaftman),
*    Adding wavefunction_scalefactor variable to config, controls scaling of wavefunctions on QW diagrams,
*    Many small fixes.

## v.0.9 (November 10th, 2013)

*    Many bugfixes,
*    Speed improvements,
*    Cython code additions,
*    Rewritten VBMAT-V part to use numpy better,
*    Merged conduction and valence band calculations.
*    This is a mostly bugfix version. Code is heavily modified and stabilized with 28 commits.

## v.0.8 (July 7th, 2013)

*    Numpy version is restructured,
*    Input file structure and sample inputs are changed,
*    Non-parabolicity of conduction band is implemented (Numpy version only),
*    Database is changed to a more clearer structure,
*    Exchange Interaction Potential is implemented (Numpy version only),
*    Logging with timers and some customizations in config,
*    A sample input to show methods of looping the simulation over a parameter is added,
*    New materials InAs, InP, AlP, GaP and new alloys InGaAs, InGaP, AlInP,
*    Strain included valence band calculation with 3×3 k.p model.
*    Many bugfixes and small corrections,

## v.0.7 (April 27th, 2013)

*    New feature: External electric field,
*    Code is optimized heavily.
*    Numpy version of Aestimo.

## v.0.6 (February, 2013)

*    New configuration properties for outputs and resultviewer.

## v.0.5 (January, 2013)

*    Minor bug fix for calculating electric field.
*    Database support for alloys.

## v.0.4 (November, 2012)

*    Major enhancements.
*    Heterojunction support.
*    Implementing shooting method for Schrödinger eq. (from P. Harrison)

## v.0.3 (October, 2012)

*    Input file support and database support for materials are added.
*    AlAs, GaAs and AlGaAs are added.

## v.0.2.1 (September, 2012)

*    Few minor corrections.
*    Plotting.
*    First documentation draft.

## v.0.2 (June, 2012)

*    Simple Poisson solver for silicon p-n homojunction.

## v.0.1 (March, 2012)

*    Simple electric field calculator for silicon.

The project began in early 2012. Its website moved from the Bitbucket wiki to
a separate domain in 2013, and the code/issue management moved to GitHub in 2017.
The historical elementary solver is now the separate
[Aestimo EDU project](https://github.com/aestimosolver/aestimo-edu).
