# Numerical Stability Notes

- Case: `R4_trap_density_50`
- Official current at -1.0 V: 0.000e+00 A/cm^2
- Whole-device probe current at -1.0 V: 0.000e+00 A/cm^2
- Abs-max current proxy at -1.0 V: 1.224e+01 A/cm^2
- All target biases reached: True
- Timeout events: 0
- Max Gummel iterations accumulated at a bias point: 1864
- Minimum attempted sub-step: 0.04999999999999993

## Notes

- Raises only the phenomenological trap-center density that multiplies the active Hurkx-like path.
- The official `av_curr` series contains exact zero-like points, so a numerical floor or extraction dead-zone is present at least in the reported terminal current.
