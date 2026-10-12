# Development

## Project aims

Aestimo is developed as an educational and research tool. Its source should be
readable enough to follow a calculation, while numerical methods should be
checked for accuracy and convergence. Python input files let users describe
layers, material parameters and numerical controls with ordinary Python
statements; the desktop interface offers JSON projects for the same broader
simulation work.

Contributions to the solver, material database, examples, tutorials and
interface are welcome. The [source and issue tracker](https://github.com/aestimosolver/aestimo)
are hosted on GitHub. Development discussion takes place on the
[Aestimo development group](https://groups.google.com/g/aestimo-devel).

## Contribute a change

1. Open an issue to discuss a substantial feature or scientific-model change.
   For a bug, include the input, Aestimo version, full error and relevant output.
2. Fork the repository and make a focused change. Explain the resulting behavior
   and update examples or documentation where needed.
3. Run the relevant checks. For solver changes, compare representative numerical
   outputs and include the physical assumptions, units and convergence evidence.
4. Submit a pull request to `aestimosolver/aestimo`, including the commands used
   to verify it. Documentation corrections can be small standalone contributions.

## Set up a checkout

Clone your fork or the upstream repository, then install it in an activated
virtual environment:

```bash
git clone https://github.com/aestimosolver/aestimo.git
cd aestimo
python -m pip install -e .
```

See [Installation](installation.md) for Python and desktop requirements.
Editable installation uses the checkout's code. Package metadata, dependencies,
data-file declarations and launch commands are defined in `pyproject.toml`.

## Run software checks

From the checkout root on Linux/macOS:

```bash
MPLBACKEND=Agg python -m unittest discover -s examples -p 'test_*.py'
```

On Windows PowerShell:

```powershell
$env:MPLBACKEND = "Agg"
python -m unittest discover -s examples -p 'test_*.py'
```

The runtime workflow checks installed packages on Linux/Python 3.9, 3.12 and
3.13 and Windows/Python 3.12. Linux jobs also run the example regression suite.
The installed-package checker verifies launch commands, a representative CLI
calculation, plotting and bundled resources without opening a GUI window.
See [Validation and reference data](validation.md) for interpreting numerical
and reference comparisons.

## Build the documentation

```bash
python -m pip install -r docs/requirements.txt
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html`. The build uses Python 3.12 and renders
Markdown with MyST. It does not import the solver or run simulations.
Read the Docs uses the configuration in `.readthedocs.yaml` at the repository
root and the pinned dependencies in `docs/requirements.txt`.

The public guide covers installation, calculations and project information.
Release checklists and historical integration records are maintained separately
in the repository for contributors.

## Build distributions

```bash
python -m pip install build twine
python -m build
python -m twine check dist/*
```

Install the resulting wheel in a fresh environment outside the checkout and
run `python tools/check_release_install.py --source .` with that environment's
interpreter. A local build does not upload the package or create a release.
The [release preparation checklist](https://github.com/aestimosolver/aestimo/blob/master/docs/releasing.md)
records the remaining publication steps.

## Background reading

Background references on quantum confinement, computational physics and
semiconductor-device simulation:

- Paul Harrison, *Quantum Wells, Wires and Dots*.
- Morten Hjorth-Jensen, *Computational Physics*, University of Oslo lecture notes
  (2007 and later notes).
- Jaan Kiusalaas, *Numerical Methods in Engineering with Python*, Cambridge (2010).
- Jonathan Edward Moussa, *The Schrodinger-Poisson Self Consistency in Layered
  Quantum Semiconductor Structures*, MSc thesis (2003).
- Nuno Sucena Almeida, *Development of a Multi-physics Simulation Framework for
  Semiconductor Materials and Devices*, thesis (2005).
- Jeffery Wayne Allen, *A Simulation Tool Suite for the Modelling and Optimization
  of Multiple Quantum Well Structures*, dissertation (2005).
