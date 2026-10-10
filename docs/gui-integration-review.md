# GUI v4 integration review

This branch integrates the upstream GUI work for review. It is not approval to
release version 4.0.0 or evidence that all device models are experimentally validated.

## Source and merge decisions

- Base: `sblisesivdin/aestimo` master, `19674cf5293c5806b27556b3855ffc08dad7ca38`.
- Contribution: `aestimosolver/aestimo` feat/gui,
  `b32582a375b9d42f26a0c65882981331140bd016`.
- The two repositories' master tips were identical when this review began.
- The contribution and master diverged from `fa0b911`; 36 files conflicted.
  Master's three commits squash the early GUI changes, whose original commits
  are also in the contribution's history. Those changes are therefore not
  three independently missing fixes.
- Conflicts use the evolved contribution versions of the GUI, plotting, core,
  packaging, ignore rules and JSON presets. The existing master-only root example
  `sample_1qw_barrierdope_ingaas.json` is retained. Both parent histories are kept.

## Reproducible tests

In a source checkout with the declared dependencies installed:

```bash
MPLBACKEND=Agg python -m unittest discover -s examples -p 'test_*.py' -v
```

The original 45-test suite reproduced as 39 passes and six failures. All six
GUI failures required pre-existing LED/laser solver output directories that
were absent from the repository. The amended suite has 47 passing tests on
Python 3.12, including two new checks against unrelated working-directory data.

- LED/laser rendering tests create deterministic synthetic `av_curr.dat`
  sweeps in temporary directories. They check figures, current-density units,
  terminal-voltage series resistance and finite laser metrics. They do not
  claim experimental agreement or validate the drift-diffusion solver.
- These fixtures explicitly carry `MODEL-BASED / NOT EXPERIMENTALLY VALIDATED`.
  Previously asserted threshold values from unavailable calculation artifacts
  are replaced with range checks appropriate to rendering tests. Existing
  analytical laser characterization tests remain intact.
- LED/laser figure builders only read their requested output directory or its
  `sim_output` subdirectory; they no longer read unrelated CWD output folders.
- QW benchmark functions accept an optional `output_dir`. Tests use temporary
  directories rather than overwriting tracked reports and images. The
  developer-specific Windows artifact-copy paths are removed.
- Publication export tests also write into temporary directories.

This command covers the example unittest suite. It does not exercise a live
desktop session, an installed wheel, every legacy calculation mode or a full
Mode 10 device solve.

## Mode 10 acceptance and failure handling

The next increment replaces the potential-correction shortcut with a requirement
that the complete residual infinity norm is below tolerance. The default
normalized tolerance remains `0.02`; this change does not establish that this
threshold is adequate for every device. Equilibrium is checked before any
factorization, and the state after the last permitted update is checked too.

- Non-finite states, residuals, Jacobians and Newton corrections, and singular
  linear systems fail the voltage step. Current arrays from a previous accepted
  step are cleared; failed states are not used to calculate/export currents.
- `last_diagnostics` records voltage, iteration count, residual norm, component
  norms for Poisson/electron/hole equations when available, and the stop reason.
- The core raises `NewtonConvergenceError` before storing a failed voltage step.
  The error includes voltage and diagnostics. Existing GUI workers catch it and
  take their error path. Output from an earlier run in an existing directory is
  not deleted by this change; a failed run must not be treated as a fresh result.
- `dd_max_iterations` is now read from the input into the model and honored by
  the Newton solver (default 25, non-negative integer). A zero budget checks
  only the initial residual and cannot bypass the core failure guard.
- A single-bias Mode 10 calculation now executes the solver, instead of skipping
  it as an equilibrium-only legacy branch. Modes 7–9 retain their existing
  routing; their convergence behavior has not been audited in this increment.

The full example unittest suite now has **61 passing tests**, including 14 new
Mode 10 tests. Central differences check the analytic Jacobian on a small
nonequilibrium system. Equilibrium and a small-bias solve check zero current and
total-current conservation. Failure tests cover stalled continuity, iteration
limits, singular/non-finite systems, invalid controls, stale currents, and core
routing that refuses to export a partially failed sweep.

A real 20-node Si p-n calculation completes a two-point low-bias sweep. As an
additional smoke check, the unmodified `examples/sample_pn.py` preset completes
41 steps from 0 to 0.8 V. The final residual infinity norm is approximately
`8.95e-6` (default acceptance tolerance `0.02`). This smoke result verifies
execution and numerical acceptance, not agreement with experiment or device
accuracy under grid refinement. A live GUI session was not exercised.

## Reference provenance and numerical agreement

The next increment separates numerical agreement from reference provenance.
The QW report requires both NRMSE ≤ 5% and residual-based R² ≥ 0.95 for
`METRICS PASSED`. These are provisional comparison thresholds, not a universal
physical acceptance rule. Pearson r² is displayed separately and never decides
acceptance; a failed comparison is not automatically called calibrated.
Invalid, constant or unpaired data are marked `NOT ASSESSABLE`/`NOT ASSESSED`.
Curve interpolation is restricted to the simulation domain without extrapolation.

QW plot statistics are computed from the actual plotted curve pairs rather than
literal/default R² and RMSE values. Regenerated Dingle and Miller Markdown
reports now mark their comparisons `METRICS FAILED`. Previously cached QW PNGs
with obsolete statistics are removed; regenerate them with the benchmark
scripts after checkout. Numerical data and solver model parameters are unchanged.

GUI fallback and legacy trace-record labels no longer promote a preset name,
a requested legacy label or a high score to experimental validation. Unvalidated
labels also no longer receive the green badge through a substring match.
The preset audit uses CSV reference origins instead of certifying devices by
filename. Calibration is identified only when an explicit fitting claim exists.

`docs/reference-data-audit.json` inventories all 23 reference CSVs with SHA-256
fingerprints and inherited source/extraction claims: 3 synthetic/model references,
1 literature-based unverified reference, and 19 claimed measurements or
digitizations awaiting source review. **No original measurement, paper figure,
table or extraction project has been independently matched in this increment.**
The source citations are preserved as claims, not newly verified facts.
Si I–V is model-defined in its header; Si C–V and InGaAs I–V explicitly state
synthetic origins. None can establish experimental agreement.

Rebuild the inventories using:

```bash
python examples/audit_reference_data.py
python examples/audit_device_examples.py
```

The example unittest suite now has **69 passing tests**. New regressions cover
perfect correlation with wrong offsets/signs, invalid and constant data, honest
no-overlap comparisons, plot-score recomputation and source-status distinctions.
The actual Miller benchmark tests require both inaccurate comparisons to fail.
The generated reports retain RMSE/MAPE and both types of R² for review.

## Mode 10 mesh and tolerance investigation

`examples/mode10_mesh_audit.py` runs nine dark Si pn cases (three meshes,
three tolerances) with isolated output. See `docs/mode10-mesh-audit.md` and
its JSON diagnostics. Eight cases complete 21 biases; the finest mesh at
1e-8 fails at 0.02 V and is retained as a failure.

The inherited 0.02 tolerance accepts a spatial total-current span of about
277 mA/cm² at 0.4 V on the 5 nm mesh, despite a median current of only
0.303 mA/cm². A small difference between median currents across meshes would
hide this conservation defect. At 1e-7 the span falls to about 8.90e-6 mA/cm²,
but low-bias currents remain comparable to or below numerical variations.
This investigation does not certify grid convergence or physical accuracy.

Mode 10 accepts the positive finite input `dd_residual_tolerance` (default
0.02 retained during review). It enables explicit tolerance studies without
patching solver code. The suite now has **72 passing tests**, including input
control and comparisons that cannot substitute a failed finest reference.
Current conservation, near-equilibrium cancellation and equation scaling are
release blockers, alongside the existing scientific and installation checks.

## Current conservation and low-current precision

Mode 10 now requires the total-current span over all mesh edges to be within
`dd_current_atol + dd_current_rtol * abs(interior_median_current)` in addition
to the residual criterion. Defaults are 1e-8 mA/cm² absolute and 1e-3 relative,
provisional numerical controls requiring physical review. Failed conservation
continues Newton updates or raises before export; currents are never flattened.

Continuity and export share a stable SG flux calculation; near-equilibrium
subtractions use expm1/log1p. Potential/carrier updates use longdouble while
sparse factorization stays float64. The tested Linux platform has 64 significant
state bits; Windows and other platforms may offer no extra precision and have
not been checked. Equations and material parameters are unchanged.

`docs/mode10-current-conservation.md` and `docs/mode10-current-audit.json`
record the current result; the previous mesh JSON is a historical baseline.
All nine standard mesh/tolerance cases now complete 21 biases. Three stricter
absolute-tolerance cases fail honestly near equilibrium. Low-bias grid errors
and the platform-dependent precision floor remain unresolved. The full example
suite has **77 passing tests**, including five new flux/conservation regressions.

## GUI example resources in installed packages

This increment stays within the GUI integration scope. The previous wheel
contained zero example JSON files and zero reference CSVs, leaving preset
selection and reference loading unavailable after installation.

The wheel and sdist now include `aestimo_examples`: 53 project presets, the
existing audit JSON, 23 CSVs and the reference README. The examples directory
is mapped to that package without duplicating the source data. Both setup
metadata paths declare the same package/data mapping.

`aeslibs.gui_resources` locates checkout or installed resources. Installed GUI
launches copy missing presets/references into `~/Aestimo/examples` (or
`$AESTIMO_WORKSPACE/examples`) so the existing example-relative save/output
paths are writable. The copy step preserves existing user files and leaves
package resources unchanged. Source checkout behavior remains unchanged.
The GUI's manual and async diode comparison paths now resolve preset-relative
references through this directory rather than assuming an adjacent source tree.

Verification on Python 3.12:

- **81 example unittest tests pass**, including four resource/path regressions.
- A wheel was built and installed with `--no-deps` in a temporary target.
  Checks ran outside the checkout, imported the installed GUI, resolved every
  preset reference, and verified exact reference bytes and preserved user files.
  Dependencies came from the test runtime; this was not a fresh dependency install.
- A wheel rebuilt from the sdist contains matching presets, CSVs and GUI resource
  code. `aestimo --help` works from the installed target.
- Both GUI entry points load/call `main` with a mocked window. **No live Tk
  window, display interaction or target-platform GUI behavior was tested.**

No solver equations, acceptance settings, material parameters or CLI behavior
are changed in this packaging increment. No PyPI release is published.

## Clean installation and owner-approved GUI scope

The owner confirmed that the GUI was working and requested skipping detailed
GUI checks on 2026-10-05. Live Tk interaction and thread/shutdown testing are
therefore **out of scope for this integration review**, not a gate that keeps
this task open. No additional solver development is introduced in this increment.

A new Python 3.12 virtual environment installed the wheel and all declared
runtime dependencies from scratch, without system-site packages or PYTHONPATH.
`pip check` found no broken requirements. Installed GUI imports and preset
resource initialization succeeded outside the source checkout. The source's
**81 unittest tests also passed using this clean environment**. Resolved
versions were NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, CustomTkinter 6.0.0,
Pillow 12.3.0, darkdetect 0.8.0 and packaging 26.3. These are observations of
this test environment, not a claim that every supported platform/version is tested.

The declared build minimum setuptools 61.0.0 fails before metadata generation
on Python 3.12 (`pkgutil.ImpImporter` is absent). Building with setuptools
77.0.3 and wheel 0.43.0 succeeds with the current license metadata and package
mapping. The build requirement now specifies `setuptools>=77.0.3`, matching
the tested build backend. Runtime dependencies and supported-Python metadata
are unchanged. No new version or PyPI release is published.

## Separate follow-up work before a formal release

1. Extend the Mode 10 audit to mesh convergence, heterogeneous devices and
   illumination/recombination regimes. The residual acceptance and failure
   propagation fixes have been completed, but the new Si mesh audit exposes
   inadequate current conservation at the inherited residual-only tolerance.
   A current-conservation guard and stable flux evaluation now address this;
   low-current grid accuracy and platform precision still need review.
2. Establish observable-specific acceptance criteria using measurement uncertainty
   and grid/model convergence. The previous correlation-only rule has been
   replaced; the provisional QW thresholds still require scientific review.
3. Obtain and review original figure/page and extraction projects for all CSVs
   claiming measurement/digitization origins. Reference origins are now
   distinguished, but independent experimental validation remains pending.
4. Review calibration of physical material parameters, including the Si
   bandgap adjustment described in `VALIDATION_RESULTS.md`.
5. Compare representative legacy Schrödinger-Poisson outputs with master.
6. Other supported Python/platform combinations remain unverified. Wheel/sdist
   resources, installed lookup, entry-point wiring and clean dependency installation
   have been checked on Python 3.12. Detailed GUI testing is waived by the owner
   for this PR and is not a remaining integration task.

The agreed GUI integration and installation review is complete: 81 tests pass,
clean Python 3.12 installation is checked, and detailed GUI testing is waived
by the owner. This PR can move to ready for review; merging remains pending
collaborator feedback. The items above are separate follow-up work, not an
invitation to expand this GUI PR. Issues are disabled on this fork, so this
numbered list is the current follow-up tracker. Do not infer formal release
readiness from the inherited version number or validation labels.
