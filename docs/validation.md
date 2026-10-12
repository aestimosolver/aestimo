# Validation and reference data

The GUI integration review checked 81 example tests and a clean Python 3.12
package installation. Those checks do not independently establish experimental
accuracy for all devices or observables.

The repository distinguishes numerical regression, comparisons with analytical
models, calibration and comparisons against claimed measurement/digitization
sources. Correlation alone is insufficient to demonstrate agreement: an offset
or wrong amplitude can retain high correlation.

## Reading the reference material

- [Current-conservation audit](mode10-current-conservation.md): Mode 10 numerical
  acceptance and limitations near equilibrium.
- [Mesh audit](mode10-mesh-audit.md): the earlier mesh/tolerance investigation;
  read it alongside the subsequent conservation work.
- [Provenance dossier](reference-data-provenance.md): contributor-supplied citations,
  claimed figure/table origins and extraction details. These claims still require
  independent source matching; inclusion here does not certify them.
- [GUI integration review](gui-integration-review.md): historical review steps,
  installation evidence and separately tracked scientific follow-up work.

Synthetic references support numerical checks, not experimental validation.
Check the reference CSV headers and audit records before citing a comparison.
