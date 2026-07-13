# R6 Terminal-Current Observable And Unit Audit

## Scope

This follow-up reuses the completed reverse-transport validation outputs. It does not introduce new device sweeps or change the authoritative JSON structure. Its purpose is to determine whether the saved solver state contains a physically usable reverse terminal-current observable.

## Unit Audit

- The active scheme-7 current routine computes `Jelec` and `Jhole` in SI `A/m^2`, then the scheme-7 path multiplies them by `0.1` before `Jtotal`, `av_curr`, and the saved result arrays are formed.
- Therefore the native saved current array is best interpreted as `mA/cm^2`; converting it to `A/cm^2` requires multiplying by `1e-3`.
- The earlier reverse-validation export used the historical `A/m^2 -> A/cm^2` factor of `1e-4`, so the previously reported nonzero current magnitudes are low by a factor of 10. Exact-zero and step-like behavior are unchanged.

## Corrected Representative Values At -1 V

- `R1_selective_bc0p0`: official `0.000e+00` A/cm^2, whole median `0.000e+00` A/cm^2, abs-max internal proxy `9.882e+01` A/cm^2.
- `R2_ohmic_bc0p6`: official `-2.555e-06` A/cm^2, whole median `0.000e+00` A/cm^2, abs-max internal proxy `1.273e-03` A/cm^2.
- `R4_tat_disabled`: official `1.145e-08` A/cm^2, whole median `1.525e-09` A/cm^2, abs-max internal proxy `2.295e-05` A/cm^2.
- `R4_tat_strong`: official `0.000e+00` A/cm^2, whole median `0.000e+00` A/cm^2, abs-max internal proxy `2.647e+02` A/cm^2.
- `R5_smoothed_interfaces`: official `0.000e+00` A/cm^2, whole median `0.000e+00` A/cm^2, abs-max internal proxy `9.882e+01` A/cm^2.

## Current-Conservation Diagnosis

- In `R1_selective_bc0p0`, the -1 V spatial profile has median current `0.000e+00` A/cm^2 and abs-max current `9.882e+01` A/cm^2, giving a conservation ratio of `9.882e+31`.
- In `R4_tat_disabled`, the -1 V spatial profile has median current `1.555e-09` A/cm^2 and abs-max current `2.295e-05` A/cm^2, giving a conservation ratio of `3.098e+03`.
- A steady one-dimensional terminal current should be nearly position independent. Ratios this large mean a contact or edge probe cannot be promoted to a physical terminal current without first fixing current continuity and boundary treatment.

## Decision

A corrected unit conversion increases the magnitude of nonzero currents, but it does not recover a gradual, conserved reverse leakage branch. The best current signal remains spatially localized and contact-sensitive rather than terminal and conserved.

Conclusion B remains supported: the present solver/contact framework fundamentally limits quantitative high-indium reverse-bias leakage simulation in this workflow. The scientifically defensible path is to report forward/calibrated behavior with the existing solver only after validation, and to treat reverse leakage as a methodology limitation unless a new high-field/contact transport implementation is added.

## Generated Artifacts

- Summary table: `results\reverse_transport\r6_terminal_observable\terminal_observable_summary.csv`
- Representative reverse observables plot: `results\reverse_transport\r6_terminal_observable\representative_reverse_observables.png`
- Representative -1 V spatial current plot: `results\reverse_transport\r6_terminal_observable\representative_spatial_current_profiles_m1v.png`
- Corrected per-case current tables: `results/reverse_transport/r6_terminal_observable/corrected_current_tables/`
