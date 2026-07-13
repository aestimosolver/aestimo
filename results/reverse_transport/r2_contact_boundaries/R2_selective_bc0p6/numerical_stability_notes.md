# Numerical Stability Notes

- Case: `R2_selective_bc0p6`
- Official current at -1.0 V: 0.000e+00 A/cm^2
- Whole-device probe current at -1.0 V: 0.000e+00 A/cm^2
- Abs-max current proxy at -1.0 V: 1.226e-04 A/cm^2
- All target biases reached: True
- Timeout events: 0
- Max Gummel iterations accumulated at a bias point: 1488
- Minimum attempted sub-step: 0.04999999999999993

## Notes

- Selective photovoltaic contact boundary conditions from the baseline JSON.
- The official `av_curr` series contains exact zero-like points, so a numerical floor or extraction dead-zone is present at least in the reported terminal current.
