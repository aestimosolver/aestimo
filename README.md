# Aestimo 1D

A one-dimensional semiconductor heterostructure simulator with
Schrödinger–Poisson and drift–diffusion solvers, a command-line interface and a
desktop GUI. Aestimo supports educational and scientific calculations of
quantum wells and semiconductor devices.

## Install the current source

In an activated Python virtual environment:

```bash
python -m pip install "git+https://github.com/aestimosolver/aestimo.git@master"
aestimo --help
aestimo-gui
```

Python 3.9 or newer is required. The GUI also needs Tk and a graphical desktop.
From a downloaded or cloned checkout, use `python -m pip install .` instead.
Installed launch commands require no manual source-directory `PATH` changes.

`python -m pip install aestimo` installs the release available on PyPI, which
may differ from the current source. The merged GUI work has not been published
as a new PyPI release by this change.

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

The [user guide](docs/index.rst) is maintained in `docs/` and configured for
Read the Docs. A hosted URL will be linked here after the project is connected.

- [Installation](docs/installation.md)
- [CLI quick start](docs/quickstart.md)
- [GUI workflow](docs/gui.md)
- [Inputs and units](docs/inputs.md)
- [Outputs](docs/outputs.md)
- [Tutorials](docs/tutorials.md)
- [Development and documentation builds](docs/development.md)
- [Validation and reference-data status](docs/validation.md)

The source version is recorded in `pyproject.toml`. Preset names, plot labels
and passing regression tests do not independently certify experimental accuracy.

## Support and license

Report bugs or ask questions through [GitHub issues](https://github.com/aestimosolver/aestimo/issues).

Aestimo is distributed under the GNU GPL v3 or later; see [COPYING.md](COPYING.md).
The project was initiated by Sefer Bora Lisesivdin, with major contributions
from Robert J. Steed, Hamza Hebal and other
[contributors](https://github.com/aestimosolver/aestimo/graphs/contributors).
