# Aestimo 1D

A one-dimensional semiconductor heterostructure simulator with
Schrödinger–Poisson and drift–diffusion solvers, a command-line interface and a
desktop GUI. Aestimo supports educational and scientific calculations of
quantum wells and semiconductor devices.

## Install

Install the published package in an activated Python environment:

```bash
python -m pip install aestimo
```

The desktop GUI described below is part of the upcoming 4.0.0 release.
Until it is published on PyPI, install the preparation source:

```bash
python -m pip install "git+https://github.com/aestimosolver/aestimo.git@master"
aestimo --help
aestimo-gui
```

Python 3.9 or newer is required. The GUI also needs Tk and a graphical desktop.
From a downloaded or cloned checkout, use `python -m pip install .` instead.
Installed launch commands require no manual source-directory `PATH` changes.

## Run a calculation

From a source checkout:

```bash
aestimo -i examples/sample_1qw_barrierdope_ingaas.py
```

Add `-d` to display plots. The CLI accepts Python inputs and GUI JSON projects;
results are written to `<input-stem>_output/` in the current working directory.

Installed GUI presets and reference CSVs are copied to `~/Aestimo/examples`
on first launch without overwriting existing user files. Set `AESTIMO_WORKSPACE`
to choose a different workspace root. A checkout uses its own `examples/`.

## Documentation

Read the [Aestimo documentation](https://aestimo.readthedocs.io/)
for installation, calculations and project information. The guide source is
maintained in `docs/`.

- [About and citation](docs/about.md)
- [Models and materials](docs/models.md)
- [Installation](docs/installation.md)
- [Your first calculation](docs/first-calculation.md)
- [CLI quick start](docs/quickstart.md)
- [GUI workflow](docs/gui.md)
- [Inputs and units](docs/inputs.md)
- [Outputs](docs/outputs.md)
- [Tutorials](docs/tutorials.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Authors and contributors](docs/authors.md)
- [Release notes](RELEASE_NOTES.md)
- [Publications](docs/publications.md)
- [Development](docs/development.md)
- [Validation and reference-data status](docs/validation.md)

The source version is recorded in `pyproject.toml`. Preset names, plot labels
and passing regression tests do not independently certify experimental accuracy.

## Support and license

Report bugs or ask questions through [GitHub issues](https://github.com/aestimosolver/aestimo/issues).

Aestimo is distributed under the GNU GPL v3 or later; see [COPYING.md](COPYING.md).
The project was initiated by Sefer Bora Lisesivdin, with major contributions
from Robert J. Steed, Hamza Hebal and other
[contributors](https://github.com/aestimosolver/aestimo/graphs/contributors).
