# Development and documentation

## Editable installation

```bash
python -m pip install -e .
```

Package metadata, dependencies, data-file declarations and launch commands live
in `pyproject.toml`. `setup.py` is only a compatibility entry point.

## Build distributions

```bash
python -m pip install build twine
python -m build
python -m twine check dist/*
```

Check a wheel outside the source checkout in a fresh virtual environment. Verify
`pip check`, CLI help, GUI imports and bundled preset/reference resources.
A successful local build does not publish anything on PyPI.

## Build the documentation

From the repository root:

```bash
python -m pip install -r docs/requirements.txt
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` to read the result. Sphinx does not import
the solver, start Tk or execute simulations. Existing Markdown review documents
are rendered through MyST rather than copied into separate manuals.

Read the Docs uses `.readthedocs.yaml`, Python 3.12 and the pinned documentation
requirements. In the Read the Docs dashboard, import or connect
`aestimosolver/aestimo`, choose `master` for the development version and enable
builds. The configuration is ready for that service, but merging this PR alone
does not create or activate a hosted project. Set the public documentation link
in the README and package metadata after the actual project URL is confirmed.

## Example regressions

```bash
MPLBACKEND=Agg python -m unittest discover -s examples -p 'test_*.py'
```

Detailed live GUI testing is outside this documentation/packaging increment.
Solver refactoring should be a separate change, with representative numerical
outputs compared before and after changes.

## Project aims and contribution workflow

The original website describes Aestimo as an educational and research tool
whose source should remain readable. Preserve that aim when adding features:
keep interfaces understandable and distinguish numerical correctness from
agreement with an experimental device. Python input scripts remain supported
alongside the GUI's JSON projects.

Develop in a fork, document the changed behavior and run checks appropriate to
the change. Submit a focused pull request to the upstream repository when the
work is ready. Include enough input and numerical detail to reproduce a solver
fix. Scientific-model changes need their own evidence; passing an installation
check is not a substitute for a physical validation argument.

Development discussion also takes place on the
[Aestimo development group](https://groups.google.com/g/aestimo-devel).
For bugs and reproducible support requests, use
[GitHub issues](https://github.com/aestimosolver/aestimo/issues).

## Reading list retained from the original site

These references formed part of the original development guidance. They are
background reading, not a claim that every current solver follows every method
in these works.

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

The old site's empty contribution placeholders have been replaced by the
workflow above. See the [migration record](website-migration.md) for sources.
