# Mode 10 current conservation and roundoff

## Behavior

An accepted Mode 10 state must satisfy both the normalized residual tolerance
and total-current conservation over **all mesh edges**:

`max(Jtotal) - min(Jtotal) <= dd_current_atol + dd_current_rtol * abs(median_current)`

The median uses the existing interior 20–80% interval. Defaults are
`dd_current_atol = 1e-8` mA/cm² and `dd_current_rtol = 1e-3` (dimensionless).
The absolute tolerance must be positive and finite; the relative tolerance
must be finite and non-negative. These are provisional numerical defaults,
not measurement uncertainty or universal physical acceptance thresholds.
The inherited `dd_residual_tolerance = 0.02` remains an independent necessary
condition. A small residual alone can no longer accept a nonconserved current.

A conservation failure continues Newton updates within `dd_max_iterations`.
If the budget is exhausted, diagnostics carry the observed span, limit and
failure reason. Current arrays are cleared and the existing core failure guard
prevents exporting a partial sweep. Neither the spatial span nor the median
is artificially flattened to make a state pass.

## Numerical changes

Continuity residuals and exported currents now use one shared
Scharfetter–Gummel edge-flux routine. Near electrochemical equilibrium it
factors the subtraction of large opposing drift/diffusion terms into an
`expm1` of the departure from equilibrium. Relative carrier logarithms use
`log1p` for small departures; larger departures retain the usual SG expression
to avoid exponential overflow. Exact equilibrium produces zero flux even
with large equilibrium carrier contrasts.

Potential and carrier states retain Newton updates in NumPy `longdouble`,
while the sparse Jacobian, factorization and residual vector remain float64.
Potential differences are formed from the departure from equilibrium in both
Poisson and continuity equations. The analytic Jacobian still passes central
difference checks. The equations, recombination model, mobility model and
material parameters are unchanged by this refactoring.

`longdouble` provides 64 significant bits on the tested Linux platform. Its
precision depends on the platform and can equal float64, particularly on
Windows. No desktop/platform equivalence is claimed. Conservation still gates
acceptance on platforms where extra precision is unavailable; a case can fail
rather than export an inaccurate result.

## Reproduced cases

```bash
MPLBACKEND=Agg python examples/mode10_mesh_audit.py --output docs/mode10-current-audit.json
```

The earlier `mode10-mesh-audit.json` is retained as historical evidence from
commit `4d3bf8d`, before the conservation guard and precision changes.
The new JSON contains nine cases with default conservation controls (three
meshes × three residual tolerances), plus three stricter absolute-tolerance
cases. Every solver run uses a temporary output directory.

All nine standard cases complete 21 biases from 0 to 0.4 V. The previously
failed finest-mesh case at residual tolerance 1e-8 now completes. At the
inherited residual tolerance 0.02:

| Mesh (nm) | Current at 0.02 V (mA/cm²) | Spatial span at 0.02 V | Current at 0.4 V (mA/cm²) | Spatial span at 0.4 V |
|---:|---:|---:|---:|---:|
| 10 | 3.28702e-7 | 1.73e-9 | 0.303303 | 7.05e-6 |
| 5 | 3.98699e-7 | 3.58e-9 | 0.303107 | 4.74e-6 |
| 2.5 | 4.03605e-7 | 7.82e-9 | 0.302812 | 5.03e-6 |

Before this change, the 5 nm default-residual case at 0.4 V had a spatial
span of about 277 mA/cm². That case now needs another Newton update and
passes with a span of about 4.74e-6 mA/cm². At 0.02 V the accepted 5 nm
state needs four updates, compared with two under the old residual-only rule.

The three stricter cases set `dd_current_atol=1e-9` mA/cm². All are rejected
at 0.02 V after 100 updates: spans remain about 1.73e-9, 3.58e-9 and 7.82e-9,
above their absolute-plus-relative limits. Those failures are retained in the
JSON. More Newton iterations do not remove this numerical floor.

The 5 nm median at 0.02 V is about 1.22% below the 2.5 nm result; the 10 nm
result is about 18.6% below it. Thus current conservation is improved, but
**grid convergence and low-current physical accuracy are not certified**.
Passing an absolute tolerance also does not guarantee a fixed relative
accuracy for currents near or below that tolerance. Additional meshes,
platform checks, physical acceptance criteria, heterogeneous devices and
illumination remain pending.

The unmodified `sample_pn.py` preset also completes all 41 biases through
0.8 V with the new default guards and its original 25-update budget. Its last
residual is 8.95e-6; total-current span is 5.92 mA/cm² against a relative-plus-
absolute allowance of about 1209.84 mA/cm². This is numerical acceptance only,
not validation of the large predicted high-bias current or of device physics.

## Tests

The example unittest suite has 77 passing tests. New regressions check zero
equilibrium flux across large carrier contrasts, retention of tiny potential
fluxes, agreement with the direct SG expression away from equilibrium,
rejection of residual-accepted but nonconserved states, and invalid current
controls. Existing finite-difference Jacobian, current conservation, core
failure propagation and GUI/report tests still pass.
