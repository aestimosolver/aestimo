# Mode 10 dark Si pn mesh and tolerance audit

This is a numerical investigation, not experimental validation or a grid-convergence certificate.
It uses `examples/sample_pn.py`: two 1000 nm Si layers, p/n doping 1e17 cm⁻³,
300 K, dark conditions, voltage 0–0.4 V in 0.02 V steps. Each case starts
independently from equilibrium. Mesh spacings are 10, 5 and 2.5 nm (200, 400
and 800 nodes). The Newton iteration budget is 100 updates per bias.

Reproduce in a source checkout with its dependencies installed:

```bash
MPLBACKEND=Agg python examples/mode10_mesh_audit.py --output docs/mode10-mesh-audit.json
```

The script isolates solver outputs in temporary directories, records each
accepted/failed bias and compares matched biases with the finest mesh at each
tolerance. It never substitutes a coarser reference when the finest case fails.
Raw diagnostics are in `mode10-mesh-audit.json`. Current units are mA/cm²;
current span is max(Jn+Jp)−min(Jn+Jp) over every mesh edge. In steady state,
total current should be spatially constant even when recombination transfers
current between electrons and holes.

## Observations at 0.4 V

| Residual tolerance | Mesh (nm) | Median current (mA/cm²) | Spatial current span (mA/cm²) | Residual infinity norm |
|---|---:|---:|---:|---:|
| 0.02 | 10 | 0.303336 | 272.568 | 6.13e-3 |
| 0.02 | 5 | 0.303137 | 276.764 | 1.51e-3 |
| 0.02 | 2.5 | 0.302841 | 282.206 | 4.13e-4 |
| 1e-7 | 10 | 0.303303 | 7.06e-6 | 7.19e-10 |
| 1e-7 | 5 | 0.303107 | 8.90e-6 | 3.47e-9 |
| 1e-7 | 2.5 | 0.302812 | 1.64e-5 | 1.63e-8 |

All three meshes complete all 21 biases at 0.02 and 1e-7. At 1e-8, the
10 and 5 nm cases complete, but the 2.5 nm case stops at 0.02 V after 100
updates with residual 1.79e-8. It is recorded as a failed case, not passed.

**The inherited default tolerance 0.02 is insufficient for spatial current
conservation in this example.** Agreement of median currents at 0.4 V would
conceal this defect: the 5 nm median differs from the finest mesh by only
0.098%, while its current span is about 277 mA/cm².

At 1e-7, the 5 nm median differs from the 2.5 nm result by about 0.097%
at 0.4 V and 0.691% at 0.2 V; the 10 nm differences are 0.162% and 7.16%.
This limited agreement does not certify the whole sweep. At 0.02 V the
median currents are only 3.16e-7, 3.99e-7 and 4.09e-8 mA/cm², respectively,
and their spatial spans exceed the median magnitudes. Relative errors there
are dominated by small currents and numerical cancellation. No universal
relative-only pass criterion is applied. The 1e-8 failure on the finest mesh
also shows that lowering a scalar residual threshold alone is not sufficient.

## Input control and remaining work

Mode 10 now accepts `dd_residual_tolerance`, a positive finite normalized
residual infinity-norm threshold passed directly to `solve_step`. The default
remains 0.02 to preserve existing behavior during this draft investigation;
**default acceptance must not be interpreted as physical current accuracy**.
`dd_max_iterations` controls the update budget. These controls apply to the
fully coupled Mode 10 path; they do not change legacy modes 7–9.

For Python inputs, for example:

```python
dd_residual_tolerance = 1e-7
dd_max_iterations = 100
```

The 1e-7 value is a useful investigated setting here, not a recommended
universal default. Before release, investigate equation scaling, roundoff and
current calculation near equilibrium; establish absolute and relative current
conservation criteria and a guarded acceptance policy. Extend the audit to
higher bias, additional mesh levels, heterojunctions and illumination, then
compare physical observables and model parameters with independently reviewed
references. No solver-equation or material-parameter change is made here.
