# R4 High-Field Transport Validation

This stage checks whether the implemented Hurkx-like transport enhancements meaningfully change the reverse branch and documents which high-field mechanisms are missing from the active solver path.

| Case | Official J(-1 V) (A/cm^2) | Whole J(-1 V) (A/cm^2) | Official Onset (V) | Whole Onset (V) | Timeouts | All Biases Reached |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| [R4_tat_baseline](R4_tat_baseline/numerical_stability_notes.md) | 0.000e+00 | 0.000e+00 | N/A | N/A | 0 | True |
| [R4_tat_disabled](R4_tat_disabled/numerical_stability_notes.md) | 1.145e-09 | 1.525e-10 | N/A | N/A | 0 | True |
| [R4_tat_strong](R4_tat_strong/numerical_stability_notes.md) | 0.000e+00 | 0.000e+00 | N/A | N/A | 0 | True |
| [R4_trap_density_50](R4_trap_density_50/numerical_stability_notes.md) | 0.000e+00 | 0.000e+00 | N/A | N/A | 0 | True |
