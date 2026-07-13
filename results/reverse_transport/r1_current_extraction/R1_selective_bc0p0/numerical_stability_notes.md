# Numerical Stability Notes

- Case: `R1_selective_bc0p0`
- Official current at -1.0 V: 0.000e+00 A/cm^2
- Whole-device probe current at -1.0 V: 0.000e+00 A/cm^2
- Abs-max current proxy at -1.0 V: 9.882e+00 A/cm^2
- All target biases reached: True
- Timeout events: 0
- Max Gummel iterations accumulated at a bias point: 1997
- Minimum attempted sub-step: 0.04999999999999993

## Notes

- Reduces only the electrostatic right-contact offset so the current-extraction artifact can be observed in a less contact-blocked reverse branch.
- The official `av_curr` series contains exact zero-like points, so a numerical floor or extraction dead-zone is present at least in the reported terminal current.
