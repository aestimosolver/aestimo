# Stage A Contact Calibration

This stage inspects the active contact model used by the solver and sweeps the effective contact barrier proxy (`bc_right`) plus external series resistance.

| Case | Leakage @ -1 V (A/cm^2) | n | Turn-on @ 0.5 A/cm^2 (V) | Curvature std(n) | Score |
| --- | ---: | ---: | ---: | ---: | ---: |
| [A01_barrier_0p60](A01_barrier_0p60/metrics.md) | 2.000e-09 | 1.007 | >1.00 | 0.460 | 17.621 |
| [A02_barrier_0p30](A02_barrier_0p30/metrics.md) | 2.000e-09 | 2.397 | >1.00 | 276385399196283291389701128192.000 | 16.482 |
| [A03_barrier_0p00](A03_barrier_0p00/metrics.md) | 2.000e-09 | 2.991 | 0.905 | 0.453 | 13.648 |
| [A04_rs_50](A04_rs_50/metrics.md) | 2.000e-09 | 2.991 | 0.915 | 0.459 | 13.650 |
| [A05_rs_200](A05_rs_200/metrics.md) | 2.000e-09 | 2.991 | >0.78 | 0.597 | 16.623 |

Best stage candidate: `A03_barrier_0p00`
