# About Aestimo 1D

Aestimo is a one-dimensional semiconductor heterostructure simulator intended
for education and research. Sefer Bora Lisesivdin started the project in early
2012; Robert Steed, Hamza Hebal and other contributors developed it further.

A structure is represented by layers and spatially varying material parameters
along the growth direction. Schrödinger–Poisson calculations connect confined
carrier states with the electrostatic potential. Drift–diffusion routes extend
the calculation to carrier transport. The desktop GUI supplies project editing
and device workflows alongside the Python and command-line interfaces.

## What the calculation represents

The growth-direction model resolves a one-dimensional layer profile. It does
not resolve arbitrary lateral device geometry. Input material parameters,
boundary conditions, mesh spacing and the selected solver determine what a
result means. Matching a device preset name is insufficient to establish
agreement with a measurement.

The historical site lists conduction-band electron and valence-band hole states,
carrier densities, electric fields, strain and quantum confinement. Availability
depends on the current solver route: see [Models and materials](models.md) for
the distinction between supported routes and historical option labels.

## Who can use it?

Students can start from a complete example, change a layer width or composition,
and inspect how the band profile and confined states respond. Researchers can
adapt the source and material database, keeping track of units, convergence and
model assumptions. Basic Python knowledge is useful for modifying input scripts.

- [Install Aestimo](installation.md)
- [Run a first calculation](quickstart.md)
- [Use the desktop GUI](gui.md)
- [Read the tutorials](tutorials.md)
- [Authors and contributors](authors.md)

(citing-aestimo)=
## Citing Aestimo

Please cite the project paper in scientific work using Aestimo:

H. Hebal, Z. Koziol, S. B. Lisesivdin and R. Steed,
“General-purpose open-source 1D self-consistent Schrödinger–Poisson Solver:
Aestimo 1D”, *Computational Materials Science* **186**, 110015 (2021).
[DOI: 10.1016/j.commatsci.2020.110015](https://doi.org/10.1016/j.commatsci.2020.110015).

Also record the actual Aestimo version or commit, input structure, material
parameters and numerical controls used in the calculation. The project paper
predates the current GUI and does not describe every later device workflow.

See the [historical publication list](publications.md) for examples collected
on the old site.

## License and documentation credits

Aestimo's source code is licensed under the GNU GPL v3 or later; see
[COPYING.md](https://github.com/aestimosolver/aestimo/blob/master/COPYING.md).
The guide incorporates the original project website's manual and historical
records. The original manual notice is **Copyright (C) 2013–2020 Aestimo group**,
with reuse declared under the [GNU Free Documentation License 1.3](https://www.gnu.org/licenses/fdl-1.3.html).
Source attribution and migration details are retained in the repository's
[migration record](https://github.com/aestimosolver/aestimo/blob/master/docs/website-migration.md).
