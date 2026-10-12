# Preparing Aestimo 4.0.0

## Scope

Prepare 4.0.0 from the integrated GUI and documentation source. Until release,
restrict work to installation, reproducibility, documentation and necessary bug
fixes. Defer new features and broad package/solver refactoring.

Work on `sblisesivdin/aestimo:master`. Agree separately when to send the final
release fixes upstream. The canonical release should reference the reviewed
upstream source, with the fork synchronized to it.

## Checks before publication

1. Confirm that `pyproject.toml`, `aestimo.__version__`, wheel metadata and the
   proposed `v4.0.0` tag refer to the same version/source.
2. Build the wheel and sdist, run `twine check`, and rebuild the documentation
   with warnings treated as errors.
3. Install the wheel in a fresh environment, run `pip check`, and run the
   installed-package smoke check below. Rebuild a wheel from the sdist when
   making packaging changes.
4. Run the example regression suite. Retain the tested Python/platform and
   commands in the release record; do not equate passing tests with experimental
   validation.
5. Finalize the changelog and record known numerical and platform limitations.
6. Connect the documentation project to Read the Docs, confirm a successful
   hosted build, and add its actual URL to the README and package metadata.
7. Confirm access to the existing PyPI `aestimo` project and the chosen upload
   route. The GitHub app connection does not establish PyPI publishing access.

From the checkout root:

```bash
python -m pip install -r docs/requirements.txt build twine
python -m build
python -m twine check dist/*
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
```

Create a separate environment and install the built wheel, using its interpreter
for the check. On Linux/macOS:

```bash
python -m venv /tmp/aestimo-release-env
/tmp/aestimo-release-env/bin/python -m pip install dist/aestimo-4.0.0-py3-none-any.whl
/tmp/aestimo-release-env/bin/python -m pip check
/tmp/aestimo-release-env/bin/python tools/check_release_install.py --source "$PWD"
MPLBACKEND=Agg /tmp/aestimo-release-env/bin/python -m unittest discover -s examples -p 'test_*.py'
```

On Windows, use an equivalent temporary environment and its `Scripts/python.exe`.
The smoke check runs a real installed CLI calculation in a temporary working
directory, verifies launcher metadata, imports the GUI without opening Tk, checks
preset/reference bytes and confirms preservation of edited presets.

## Publication and verification

After the agreed source is upstream, create the `v4.0.0` tag and GitHub release
from that exact commit, publish the wheel and sdist to PyPI, and select the
matching documentation version. Verify `python -m pip install aestimo==4.0.0`
from PyPI in a fresh environment. Update installation text that currently says
the GUI source is not yet published, using the actual publication status.

Do not create a new tag or upload under a version that already exists. If final
verification fails after publication, prepare a new patch version rather than
replacing release artifacts.

## Current evidence

- GUI integration and documentation were merged upstream as PRs #61 and #62.
- PR #62's documentation/packaging workflow passed on the merge commit.
- The documentation increment passed 81 example regressions, strict Sphinx
  builds, clean Python 3.12 installation and wheel/sdist resource checks.
- A real installed CLI run of `sample_1qw_barrierdope_ingaas.py` completed on
  Python 3.12/Linux during release preparation.

Hosted Read the Docs publication, the final release source, final artifact
checks and PyPI publication are still pending. This document records preparation,
not a declaration that the release has shipped.
