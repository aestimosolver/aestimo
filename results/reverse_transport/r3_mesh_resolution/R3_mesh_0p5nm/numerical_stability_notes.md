# Numerical Stability Notes

- Case: `R3_mesh_0p5nm`
- Official current at -1.0 V: -1.559e-11 A/cm^2
- Whole-device probe current at -1.0 V: 7.134e-18 A/cm^2
- Abs-max current proxy at -1.0 V: 1.427e+01 A/cm^2
- All target biases reached: False
- Timeout events: 105
- Max Gummel iterations accumulated at a bias point: 10000
- Minimum attempted sub-step: 0.0001953125

## Notes

- Uniform mesh refinement to test field and current sensitivity to spatial resolution.
- The official `av_curr` series contains exact zero-like points, so a numerical floor or extraction dead-zone is present at least in the reported terminal current.
- The whole-device current probe does not match the official right-boundary median extraction at -1 V, indicating extraction-location sensitivity.
- At least one target reverse bias was not fully reached before the adaptive stepping logic hit its minimum step size.
