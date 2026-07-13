# Numerical Stability Notes

- Case: `R4_tat_disabled`
- Official current at -1.0 V: 1.145e-09 A/cm^2
- Whole-device probe current at -1.0 V: 1.525e-10 A/cm^2
- Abs-max current proxy at -1.0 V: 2.295e-06 A/cm^2
- All target biases reached: True
- Timeout events: 0
- Max Gummel iterations accumulated at a bias point: 1124
- Minimum attempted sub-step: 0.04999999999999993

## Notes

- Disables the Hurkx activation path by setting `tat_field` above the solver threshold.
- The official `av_curr` series contains exact zero-like points, so a numerical floor or extraction dead-zone is present at least in the reported terminal current.
- The whole-device current probe does not match the official right-boundary median extraction at -1 V, indicating extraction-location sensitivity.
