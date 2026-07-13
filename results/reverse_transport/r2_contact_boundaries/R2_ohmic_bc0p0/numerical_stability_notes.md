# Numerical Stability Notes

- Case: `R2_ohmic_bc0p0`
- Official current at -1.0 V: 0.000e+00 A/cm^2
- Whole-device probe current at -1.0 V: 0.000e+00 A/cm^2
- Abs-max current proxy at -1.0 V: 1.022e+01 A/cm^2
- All target biases reached: True
- Timeout events: 0
- Max Gummel iterations accumulated at a bias point: 1107
- Minimum attempted sub-step: 0.04999999999999993

## Notes

- Combines the least restrictive electrostatic boundary with the least restrictive contact-carrier boundary type available in the current solver family.
- The official `av_curr` series contains exact zero-like points, so a numerical floor or extraction dead-zone is present at least in the reported terminal current.
