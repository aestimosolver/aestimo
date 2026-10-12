# Validation and reference data

A simulation result, a reference comparison and experimental validation are
separate kinds of evidence. A device preset provides a starting configuration;
its filename or the word "Validation" in a menu is not an acceptance result.

## Read the status as well as the score

| Status used by the current reporting code | Meaning |
| --- | --- |
| `MODEL-BASED / NOT EXPERIMENTALLY VALIDATED` | A model result without an independently established experimental comparison. |
| `SYNTHETIC REFERENCE / NOT EXPERIMENTAL` | Comparison against an analytical or generated reference rather than measurements. |
| `REFERENCE COMPARISON / PROVENANCE UNVERIFIED` | A reference is supplied, but its claimed measurement/digitization origin still needs source matching. |
| `CALIBRATED MODEL / PROVENANCE UNVERIFIED` | A fitting claim exists; it does not independently establish the reference's origin or predictive accuracy. |

For curve agreement, review errors, offsets, amplitude, sign and overlap of the
comparison domain. A large correlation score can coexist with a wrong curve.
A metric marked not assessable should not be treated as a passing comparison.

## What to check for your calculation

- Identify the reference's actual origin from its CSV header and audit record.
- Confirm the observable, units, temperature, geometry and comparison domain.
- Check numerical convergence for the observable you want to report.
- Record any fitted material or contact parameters and separate fitted data from
  independent comparisons.

Synthetic diode baselines are useful for numerical regression. Agreement with
them does not turn them into wafer measurements. Likewise, a literature citation
in a traceability record is a source claim until the original figure/table and
extraction details are matched.

## Scientific and numerical reports

- [Current-conservation audit](mode10-current-conservation.md): Mode 10 residual
  and spatial total-current acceptance, with limitations near equilibrium.
- [Earlier mesh audit](mode10-mesh-audit.md): a historical tolerance investigation;
  read it alongside the subsequent conservation work.
- [Reference provenance dossier](reference-data-provenance.md): contributor-supplied
  source and extraction claims.
- [GUI integration review](gui-integration-review.md): historical integration and
  installation evidence, plus separately tracked scientific follow-up work.

Automated regression and package-installation checks assess software behavior.
They do not certify every device model or physical observable.
