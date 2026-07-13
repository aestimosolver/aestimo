# Reverse Transport Validation Summary

This study isolates reverse-bias methodology limitations before any further device sweeps.

## High-Level Outcome

- Contact-unmasked reference case: `R2_selective_bc0p0`
- Official `J(-1 V)`: 0.000e+00 A/cm^2
- Whole-device probe `J(-1 V)`: 0.000e+00 A/cm^2
- Abs-max probe `J(-1 V)`: 9.882e+00 A/cm^2

## Conclusion

The reverse branch behavior is dominated by methodology limits rather than by a simple lack of defect strength.

A realistic gradual high-indium leakage tail was not recovered. The study instead points to a framework limitation that combines:

- extraction sensitivity, because `av_curr` samples a near-contact median rather than a robust terminal current observable
- selective-contact suppression of minority-carrier injection under reverse bias
- strong bias-step / damping sensitivity near reverse onset
- and the absence of several physically relevant high-field leakage channels beyond Hurkx-like TAT

Under these conditions, conclusion **B** is currently better supported: the present solver/contact framework fundamentally limits accurate high-indium reverse-bias simulation in its current form.

A next-priority solver task would be to replace or augment the reverse-bias contact/current methodology before resuming calibration or any design sweeps.
