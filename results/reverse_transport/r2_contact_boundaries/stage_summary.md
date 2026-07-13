# R2 Contact Boundary Analysis

This stage compares the actual contact-boundary modes available in the active solver path and checks whether reverse injection is artificially suppressed at the contacts.

| Case | Official J(-1 V) (A/cm^2) | Whole J(-1 V) (A/cm^2) | Official Onset (V) | Whole Onset (V) | Timeouts | All Biases Reached |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| [R2_selective_bc0p6](R2_selective_bc0p6/numerical_stability_notes.md) | 0.000e+00 | 0.000e+00 | N/A | N/A | 0 | True |
| [R2_ohmic_bc0p6](R2_ohmic_bc0p6/numerical_stability_notes.md) | -2.555e-07 | 0.000e+00 | N/A | N/A | 0 | True |
| [R2_selective_bc0p0](R2_selective_bc0p0/numerical_stability_notes.md) | 0.000e+00 | 0.000e+00 | N/A | N/A | 0 | True |
| [R2_ohmic_bc0p0](R2_ohmic_bc0p0/numerical_stability_notes.md) | 0.000e+00 | 0.000e+00 | N/A | N/A | 0 | True |
