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

## Remaining work before release

1. Audit Mode 10: finite-difference verification of the analytic Jacobian,
   equilibrium, current conservation, mesh convergence and failure propagation.
   A small potential correction currently permits success even when residuals
   remain large; callers can retain unconverged voltage-step results.
2. Replace misleading acceptance rules and validation labels. The Miller QCSE
   report currently passes a heavy-hole comparison with NRMSE 36.56% and MAPE
   112% because Pearson correlation is high. Correlation alone is insufficient.
3. Separate synthetic references, calibrated models and sourced experimental
   measurements. InGaAs I-V and Si C-V reference files explicitly contain
   synthetic data. Request figure/page and digitization provenance for CSVs
   claiming experimental origins.
4. Review calibration of physical material parameters, including the Si
   bandgap adjustment described in `VALIDATION_RESULTS.md`.
5. Compare representative legacy Schrödinger-Poisson outputs with master.
6. Verify wheel data resources, supported Python versions, CLI entry points,
   and live GUI thread/shutdown behavior on target desktop systems.

Keep this contribution in a draft integration PR until these scientific and
installation checks have been completed. Do not infer release readiness from
the version number or validation strings inherited from the contribution.
