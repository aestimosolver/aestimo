# Installation

## Install with pip

In an activated Python environment, install the published release with:

```bash
python -m pip install aestimo
```

[Aestimo 4.0.0 is available on PyPI](https://pypi.org/project/aestimo/4.0.0/)
and includes the desktop GUI and command-line launchers. To update an existing
installation:

```bash
python -m pip install --upgrade aestimo
```

For a reproducible installation of this release, use `aestimo==4.0.0` instead
of `aestimo`.

Check the installed package and its dependencies with:

```bash
python -m pip show aestimo
python -m pip check
```

## Create a Python environment

Aestimo 4.0.0 requires Python 3.9 or newer. Installation and CLI checks
cover Python 3.9, 3.12 and 3.13 on Linux, and Python 3.12 on Windows. These
checks import the GUI without opening a window; other Python/platform
combinations have not all been verified.

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

(install-source)=
## Install from source

To install the 4.0.0 release directly from the main Aestimo repository:

```bash
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/aestimosolver/aestimo.git@v4.0.0"
```

This command requires Git. Alternatively, download and extract the tagged release,
then run `python -m pip install .` from its root. To install the development
version, replace `v4.0.0` with `master` in the Git URL.

Installation creates `aestimo` and `aestimo-gui` commands in the active environment.
No manual addition of the source directory to `PATH` is needed.

## Using an existing conda installation

You can use conda instead of `venv`. With conda already installed:

```bash
conda create --name aes python=3.12 pip
conda activate aes
python -m pip install .
```

Run the last command from a downloaded or cloned source checkout. The active
Python's pip installs Aestimo's declared dependencies. `pylab` is not a separate
required package, and no manual `chmod` or `site-packages` path edits are needed.

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
