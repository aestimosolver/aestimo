# Calibration Assessment

The staged calibration campaign completed through contact, defect, mobility, polarization, and candidate-search passes using `examples/untitled_project.json` as the fixed baseline structure.

## Outcome

The stop condition was **not** met.

- Best overall case: `A03_barrier_0p00`
- Reverse leakage at `-1.0 V`: `2.0e-9 A/cm^2`
- Ideality factor: `2.991`
- Turn-on at `|J| = 0.5 A/cm^2`: `0.905 V`

This is still outside the requested high-indium InGaN target window:

- leakage target: `1e-5 to 1e-2 A/cm^2`
- ideality target: `2.0 to 3.0`
- turn-on target: `0.5 to 0.8 V`

## Stage Findings

### Stage A Contact Calibration

- The active contact-barrier proxy in solver 7 is the electrostatic boundary `bc_right`.
- `work_function_left/right` are not an active transport boundary in this solver path.
- External `Rs` can be post-processed for measured I-V comparison, but it does not repair the intrinsic forward bottleneck.
- Reducing `bc_right` from `0.6 V` to `0.0 V` improved turn-on from `>1.0 V` to `0.905 V` and moved ideality into the requested range.
- Reverse leakage remained orders of magnitude too low.

### Stage B Defect Calibration

- Shorter SRH lifetimes, stronger TAT, larger trap-density scaling, and trap-energy offsets all increased internal SRH/TAT dominance.
- None of these changes lifted the extracted reverse leakage at `-1.0 V` out of the `~1e-9 A/cm^2` range.
- Several cases produced unphysical forward-fit behavior, including negative ideality factors or near-zero extracted turn-on.

### Stage C Mobility Calibration

- Aggressive mobility reduction lowered forward current further and did not improve reverse leakage.
- Mobility degradation alone therefore cannot bridge the gap to the literature-calibrated device behavior.

### Stage D Polarization Calibration

- Partial polarization relaxation and reduced piezoelectric contribution had negligible impact on the terminal metrics relative to the best defect-limited chain.
- In the current solver path, polarization scaling is not the dominant missing ingredient behind the leakage shortfall.

## Important Reverse-Bias Observation

The reverse branch shows a numerical/physical inconsistency in multiple cases:

- `A03_barrier_0p00` raw reverse current stays at `0` from `-1.05 V` through `-0.80 V`, then jumps abruptly to `2.53 A/cm^2` by `-0.75 V`.
- `B03_tatfield_1e6` raw reverse current stays at `0` through `-0.30 V`, then jumps abruptly above `1 A/cm^2` at `-0.25 V`.

That shape is not a realistic high-indium leakage tail. It suggests the current solver/boundary formulation is not producing a trustworthy reverse-leakage observable for this structure, even when the internal recombination profile is strongly defect dominated.

## Conclusion

The present calibration pass shows that the device can be made more non-ideal in forward bias, but not in the experimentally required way. Under the current solver-7 boundary formulation, the limiting issue is no longer "insufficient defect strength" alone. The stronger evidence is that reverse leakage extraction itself is not behaving physically.

## Recommended Next Step

Do **not** proceed to thickness, composition, illuminated, or design-map sweeps yet.

Instead, the next technical task should be to correct the reverse-bias/contact methodology before resuming calibration:

- verify the reverse-bias sweep protocol and current extraction
- inspect whether solver 7 is appropriate for reverse-leakage calibration of this homojunction
- implement or expose a physically grounded contact-injection / tunneling boundary model if available
- then rerun the staged calibration starting from the `A03_barrier_0p00` contact setting
