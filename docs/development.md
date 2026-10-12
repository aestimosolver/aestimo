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
