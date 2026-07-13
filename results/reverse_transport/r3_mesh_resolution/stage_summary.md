# R3 Mesh And Field Resolution

This stage tests mesh sensitivity with the least blocked contact offset while keeping the device stack fixed. Native local nonuniform refinement is not exposed by this workflow, so uniform grid refinement is used instead.

| Case | Official J(-1 V) (A/cm^2) | Whole J(-1 V) (A/cm^2) | Official Onset (V) | Whole Onset (V) | Timeouts | All Biases Reached |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| [R3_mesh_1p0nm](R3_mesh_1p0nm/numerical_stability_notes.md) | 0.000e+00 | 0.000e+00 | N/A | N/A | 0 | True |
| [R3_mesh_0p5nm](R3_mesh_0p5nm/numerical_stability_notes.md) | -1.559e-11 | 7.134e-18 | N/A | N/A | 105 | False |
| [R3_mesh_0p2nm](R3_mesh_0p2nm/numerical_stability_notes.md) | -1.553e-11 | 1.011e-17 | N/A | N/A | 105 | False |
