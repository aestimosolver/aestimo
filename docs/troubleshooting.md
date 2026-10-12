# Troubleshooting

## The command is not found

Activate the environment where Aestimo was installed and run
`python -m pip show aestimo`. With a matching interpreter,
`python -m aestimo --help` provides another way to invoke the CLI.
A current pip installation creates launch commands; editing the source directory
into `PATH` or making files inside `site-packages` executable is unnecessary.

## Tk is unavailable or there is no display

The GUI requires Tk and a graphical desktop. Check Tk with
`python -m tkinter`; install it through the Python distributor or system package
manager if missing. WSL is an optional environment, not a requirement for all
Windows installations. A WSL GUI also needs display support.

For headless CLI calculations on Linux, use
`MPLBACKEND=Agg aestimo -i /path/to/input.py` without `-d`.

## The input file cannot be loaded

Check the actual path, filename extension, Python syntax and field names.
Python input files execute in order and can depend on previously defined
variables. Start from a complete example for the current source.
For JSON, save a project through the GUI and check its structure rather than
renaming a Python input to `.json`.

## A calculation does not converge

First inspect `aestimo.log` for the solver, failing bias and failure reason.
Check layer thicknesses, alloy fractions, doping units, grid spacing, contacts
and whether the chosen solver supports the requested options.

Vary the mesh and compare the observable of interest at successive refinements.
Smaller spacing increases grid size and computational cost; there is no universal
`gridfactor < 1 nm` rule. For a bias sweep, a smaller step may help diagnose
where the calculation fails.

Keep the intended physical structure when investigating numerical convergence.
Reducing doping changes the device; it is not a general numerical fix.
Tightening or loosening a tolerance alone does not prove mesh convergence.
For Mode 10, inspect both residual and spatial total-current checks in the
[current-conservation report](mode10-current-conservation.md).

## Results are missing or an older result was overwritten

The CLI uses `<input-stem>_output/` under the current working directory.
Two inputs with the same stem can target the same directory.
Use separate working directories and keep the input with its outputs.
Output names also depend on solver and configuration; see [Outputs](outputs.md).

Installed GUI examples are writable copies. Existing files are retained on an
upgrade, so an edited or older preset may remain in the workspace. Compare it
with the bundled preset when reproducing a result from another version.

## Reporting a problem

Include the Aestimo version or commit, Python/platform, command, input,
material edits and the relevant log excerpt in a
[GitHub issue](https://github.com/aestimosolver/aestimo/issues).
