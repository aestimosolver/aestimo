# Outputs and interpretation

## Where results are written

| Interface / workflow | Location |
| --- | --- |
| CLI | `<input-stem>_output/` under the current working directory. |
| Main GUI workers | `<project-name>_output/` beneath the GUI examples directory. |
| GUI Solar Study | Separate dark, illuminated and temperature output directories under the working directory. |

Saving a GUI JSON to another directory does not relocate all calculation outputs.
Check **Console** and `aestimo.log` for the selected workflow. Reusing a filename
stem can reuse its output directory.

## Schrödinger–Poisson files

The standard exporter writes whitespace-separated numerical columns:

| File | Columns / meaning |
| --- | --- |
| `aestimo.log` | Calculation messages, iterations and diagnostics. |
| `potn_eh_equi_cond.dat` | Position (m), conduction-band profile (eV), valence-band profile (eV). |
| `np_data0_equi_cond.dat` | Position (m), electron density (cm⁻³), hole density (cm⁻³). |
| `efield_eh_equi_cond.dat` | Position and the two exported electric-field arrays. |
| `sigma_eh_equi_cond.dat` | Position and the exported charge-profile array. |
| `states_e_QWR*_equi_cond.dat` | State index, energy (meV), subband population (m⁻²), effective-mass ratio `m*/m_e`. |
| `states_h_QWR*_equi_cond.dat` | State index, energy (meV), subband population (m⁻²), effective-mass ratio `m*/m_e`. |
| `wavefunctions_e_QWR*_equi_cond.dat` | Position and electron wavefunction amplitudes for the region's states. |
| `wavefunctions_h_QWR*_equi_cond.dat` | Position and hole wavefunction amplitudes for the region's states. |

`QWR` identifies a quantum region. The wavefunction exports contain amplitudes,
not precomputed `psi²`. Normalization and grid spacing matter when deriving a
probability density from those amplitudes.

Band bending describes the spatial energy profile; confined-state energies and
carrier populations answer different questions. Keep their units and reference
energy distinct when comparing a band plot with a state table.

## Transport and device workflows

Drift–diffusion routes produce `av_curr.dat` and bias-specific files, such as
`potn_eh_0.40.dat`, `np_data0_0.40.dat`, and `states_e_QWR1_0.40.dat`.
The state columns use the same energy, population and relative-mass conventions
as the equilibrium state files.

`av_curr.dat` pairs bias with the route's average-current values. Conversion
between current density and total device current requires the device area and
the route's unit convention; do not apply an area factor twice or infer units
solely from an axis label.

QW/LED/laser characterization can also write `traceability_record.json` and
`validation_report.md`. These record the workflow's comparison information;
their presence alone is not proof of experimental validation. Some figures are
held in the GUI until saved through the plot toolbar.

## Which files should be present?

Availability depends on the solver, the quantum regions and `config.py` flags
including `potential_out`, `states_out`, `probability_out`, `sigma_out`,
`electricfield_out` and `Drift_Diffusion_out`. There is no universal fixed number
of output files. A failed transport sweep can intentionally stop before export.

For a reproducible record, retain the input, material database changes,
Aestimo version or commit, solver/grid controls, relevant logs and numerical
outputs together with the exported figures.
