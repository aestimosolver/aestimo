# Command-line quick start

## Run an example

From a source checkout, install the package and run a Python input:

```bash
python -m pip install .
aestimo --help
aestimo -i examples/sample_1qw_barrierdope_ingaas.py
```

To display the plots after a calculation:

```bash
aestimo -d -i examples/sample_1qw_barrierdope_ingaas.py
```

The CLI also accepts JSON projects saved by the GUI:

```bash
aestimo -i /path/to/project.json
```

Python inputs execute as Python modules. JSON projects use a separate loader
that translates GUI field names into solver fields; these formats are not simply
interchangeable text representations. See [Inputs](inputs.md).

## Command reference

| Option | Purpose |
| --- | --- |
| `-h`, `--help` | Show available arguments and exit. |
| `-i`, `--input` | Load a Python input file or JSON project. |
| `-d`, `--drawfigures` | Display figures at the end of the calculation. |
| `-v`, `--version` | Print the local Aestimo version and exit. |

Help/version reporting and module imports do not create calculation output
files. Exit status 0 indicates normal CLI completion, 1 indicates a rejected
calculation (such as an unsupported solver or exceeded grid limit), and 2
indicates a usage/input-loading error. Missing paths, malformed Python/JSON and
a JSON root that is not a project object are input-loading errors.

For a rejected calculation, inspect `aestimo.log` in its output directory.
The error message states the rejected option or limit; no successful-completion
message or numerical results should be inferred from that log's existence.

Results are written below the current working directory as
`<input-stem>_output/`. For example, `project.json` produces `project_output/`.
Use separate working directories for inputs with the same filename stem.
Read `aestimo.log` there for convergence messages and errors.

For a non-interactive shell on Linux:

```bash
MPLBACKEND=Agg aestimo -i /path/to/project.json
```

Successful completion is a numerical result for the chosen model. Check the
mesh, convergence and model assumptions before interpreting physical observables.
