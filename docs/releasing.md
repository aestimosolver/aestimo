# Aestimo 4.0.0 release record

## Published release

Aestimo 4.0.0 was published on October 20, 2026:

- Source commit: `62ded8c5308a749dbdb5277f517b51dda10d402f`.
- [GitHub release and tag](https://github.com/aestimosolver/aestimo/releases/tag/v4.0.0).
- [PyPI release](https://pypi.org/project/aestimo/4.0.0/), containing
  `aestimo-4.0.0-py3-none-any.whl` and `aestimo-4.0.0.tar.gz`.
- All jobs in the [final source compatibility run](https://github.com/aestimosolver/aestimo/actions/runs/39430673632)
  passed. Linux jobs ran 98 regressions; Windows coverage used the installed-package
  smoke check. Python 3.9 uses Ubuntu 24.04 because setup-python has no Python 3.9
  build for Ubuntu 26.04.
- The project owner built and checked the distributions, verified a clean wheel
  installation, uploaded both files to PyPI, then repeated the dependency and
  installed-package checks in a fresh environment using `aestimo==4.0.0` from PyPI.
  All checks passed.

The sections below retain the preparation scope, checklist and earlier evidence.

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
5. Finalize `RELEASE_NOTES.md` and record known numerical and platform limitations.
6. Verify the latest hosted build at https://aestimo.readthedocs.io/ matches the
   final release source, including the banner, navigation and release notes.
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
The smoke check runs a real installed CLI calculation with and without plots
using a non-interactive backend and verifies that plotting leaves numerical
exports unchanged. It uses a temporary working directory, verifies launcher
metadata, imports the GUI without opening Tk, checks preset/reference bytes and
confirms preservation of edited presets.

## Publication and verification

After the agreed source is upstream, create the `v4.0.0` tag and GitHub release
from that exact commit, publish the wheel and sdist to PyPI, and select the
matching documentation version. Verify `python -m pip install aestimo==4.0.0`
from PyPI in a fresh environment. Revise installation text using the actual publication status.

Do not create a new tag or upload under a version that already exists. If final
verification fails after publication, prepare a new patch version rather than
replacing release artifacts.

## Current evidence

- GUI integration and documentation were merged upstream as PRs #61 and #62.
- PR #62's documentation/packaging workflow passed on the merge commit.
- The documentation increment passed 81 example regressions, strict Sphinx
  builds, clean Python 3.12 installation and wheel/sdist resource checks.
- Release startup fixes passed 85 regressions, including import/help/version
  file-creation checks, missing-input errors, plot-flag forwarding and log routing.
  A representative CLI example retained identical numerical arrays across all
  eight exported data files before and after the startup changes.
- A real installed CLI run of `sample_1qw_barrierdope_ingaas.py` completed on
  Python 3.12/Linux during release preparation.
- The `Runtime compatibility` workflow builds the wheel from the sdist, installs
  it with dependencies and runs `pip check` and the installed-package smoke
  check on Linux/Python 3.9, 3.12 and 3.13 and Windows/Python 3.12. Each Linux
  job also runs all 94 example regressions. Windows coverage is limited to the
  smoke check; no job opens a live GUI window.
- The [compatibility run for commit e8c65f0](https://github.com/sblisesivdin/aestimo/actions/runs/38251470347)
  records these release-preparation results. The workflow publishes no artifacts
  to PyPI and creates no release or tag.

- PRs #63 and #64 were merged upstream. Documentation/packaging and all four
  runtime compatibility jobs passed on merge commit
  `6c58c73d8bfec09e296793ffbc68fa2df2c1c8bf`.
- The project owner confirmed https://aestimo.readthedocs.io/ as the documentation
  address. README and package metadata link to it; the source-installation
  instructions now use the upstream repository.

The final source checks, local artifact checks and PyPI publication are complete
as recorded above. Verification of the latest hosted Read the Docs build remains
a separate deployment check.
