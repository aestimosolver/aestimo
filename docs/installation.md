# Installation

## Install the current source

Use Python 3.9 or newer. The integration was checked with a clean Python 3.12
installation; other declared Python versions are not all verified.
Create a virtual environment before installing:

```bash
python -m venv .venv
```

Activate it on Linux/macOS:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the current upstream development version:

```bash
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/aestimosolver/aestimo.git@master"
```

This command requires Git. Alternatively, download and extract the repository,
then run `python -m pip install .` from its root. Use a commit SHA instead of
`master` in the Git URL when you need a reproducible installation.

Installation creates `aestimo` and `aestimo-gui` commands in the active environment.
No manual addition of the source directory to `PATH` is needed.

## Install a published release

```bash
python -m pip install aestimo
```

This installs the version currently published on PyPI. It may differ from
`master` and may not contain the GUI documented here. This documentation change
does not publish a new release. Inspect the installed distribution with:

```bash
python -m pip show aestimo
python -m pip check
```

## Desktop requirements

The GUI needs Tk and a graphical desktop in addition to its Python dependencies.
Check Tk with `python -m tkinter`. If Tk is missing on Debian/Ubuntu, install
the system package `python3-tk` for your Python installation. Python distributions
on other platforms may provide Tk through their own installers.

For command-line calculations on a machine without a display, use a
non-interactive Matplotlib backend, for example `MPLBACKEND=Agg` on Linux.

## Installed examples

The wheel includes JSON GUI presets and their reference CSV files. The GUI copies
missing resources to `~/Aestimo/examples` on first launch and preserves existing
user files. Set `AESTIMO_WORKSPACE` to choose a different workspace root:

```bash
export AESTIMO_WORKSPACE="$HOME/my-aestimo-work"
aestimo-gui
```

The examples directory will be `$AESTIMO_WORKSPACE/examples`. In a source checkout,
the GUI uses the checkout's existing `examples/` directory instead. Python input
scripts and historical notebooks are available in the repository; they are not
all distributed as installed GUI presets.
