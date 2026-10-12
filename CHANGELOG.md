# Changelog

## 4.0.0 — release preparation

This entry describes the current source. A 4.0.0 release has not yet been
published by this preparation commit.

- Add the desktop GUI and device project workflows, including quantum-well,
  diode, solar-cell, LED and laser presets.
- Bundle GUI project presets and reference data in pip distributions and provide
  a writable installed workspace without replacing user-edited files.
- Provide `aestimo`, `aestimo-gui` and the Windows GUI launcher. Report the local
  version without a network request.
- Avoid creating output/log files during module imports, CLI help or version
  reporting. Route calculation logs to the selected output directory without
  accumulating duplicate handlers, and correctly pass the CLI `-d` plot option.
- Add a user guide configured for Read the Docs and consolidate package metadata
  in `pyproject.toml`.
- Distinguish numerical comparisons, synthetic baselines and unverified source
  claims in reference-data reporting.
- Reject unconverged Mode 10 transport states and require spatial total-current
  conservation; stabilize near-equilibrium edge-flux evaluation.

### Compatibility and limitations

Python inputs, GUI JSON projects and the existing solver interfaces remain in
use. Output filenames differ from the older five-file documentation; consult
`docs/outputs.md`. Tk and a graphical desktop are required for the GUI.

Clean installation and automated regression checks have been performed on
Python 3.12/Linux. Detailed live GUI testing was waived for the integration;
this is not a claim of complete cross-platform testing. Source metadata declares
Python 3.9 or newer, but every declared version/platform has not been checked.

Reference provenance and physical accuracy still require independent review.
Low-current mesh accuracy and platform-dependent floating-point precision remain
limitations; numerical acceptance does not certify experimental agreement.
See `docs/validation.md` and the linked audit reports.
