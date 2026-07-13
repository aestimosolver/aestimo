# Reverse Transport Discussion

The reverse-bias validation study was completed without changing the baseline layer stack from `examples/untitled_project.json`. The question was not whether a new parameter set could be found, but whether the present reverse-transport methodology can produce a physically credible high-indium InGaN leakage branch at all.

## R1 Current Extraction Validation

- In the baseline JSON case, the reported `av_curr` at `-1.0 V` is exactly `0.0 A/cm^2`.
- The whole-device median probe is also `0.0 A/cm^2`.
- The internal absolute-maximum current proxy at the same bias is already `1.23e-4 A/cm^2`.
- In the contact-unmasked `bc_right = 0.0 V` case, `av_curr(-1.0 V)` is still `0.0 A/cm^2`, while the internal absolute-maximum current proxy rises to `9.88 A/cm^2`.

This means the flat near-zero reverse branch is not explained by output clipping in `av_curr.dat`. The file writer is passing through the solver result directly. The more important issue is that the solver’s official current observable is sampled from a right-boundary median region that can remain zero even when large internal current density exists elsewhere in the device.

## R2 Contact Boundary Analysis

- Selective photovoltaic contacts with the original `bc_right = 0.6 V` keep both the official and whole-device reverse leakage at `0.0 A/cm^2` at `-1.0 V`.
- Standard ohmic contacts with `bc_right = 0.6 V` produce a small nonzero official current, `-2.55e-7 A/cm^2`, but the whole-device median remains `0.0 A/cm^2`.
- Removing the electrostatic boundary offset (`bc_right = 0.0 V`) increases internal current dramatically for both selective and ohmic contact treatments, yet the official terminal current at `-1.0 V` still remains `0.0 A/cm^2`.

The boundary-condition study therefore points to two coupled limitations:

- the continuity-equation boundary conditions strongly shape whether reverse carrier supply is allowed,
- and the terminal current observable is not robust enough to report that transport consistently.

## R3 Mesh And Field Resolution

- Uniform refinement from `1.0 nm` to `0.5 nm` and `0.2 nm` did not recover realistic reverse leakage.
- The refined meshes produced only tiny official leakage values near `-1.55e-11 A/cm^2` at `-1.0 V`.
- Both refined cases accumulated `105` timeout events and required minimum sub-steps of about `1.95e-4 V`.

So mesh refinement does change the numerical behavior, but mainly by exposing solver fragility rather than converging toward a stable physical leakage tail. This is important: the reverse-transport problem is not simply an under-resolved field spike.

## R4 High-Field Transport Validation

- The active solver family clearly contains a Hurkx-like field-enhanced TAT factor controlled by `tat_field`.
- Disabling that path (`tat_field = 1e10 V/m`) produces a small nonzero official reverse current, `1.15e-9 A/cm^2`, and a much smaller internal current proxy, `2.29e-6 A/cm^2`.
- Keeping the baseline TAT path or strengthening it drives the internal absolute-maximum current proxy upward (`9.88 A/cm^2`, `26.47 A/cm^2`, `12.24 A/cm^2` for the tested cases), while the official terminal current at `-1.0 V` remains pinned at `0.0 A/cm^2`.

This is a critical result. The implemented high-field enhancement changes internal transport strongly, but it does not translate into a stable, gradual, terminal reverse-leakage branch. That is consistent with a methodology failure, not just a tuning failure.

Also, no explicit implementation was found for:

- Poole-Frenkel emission
- band-to-band tunneling
- Fowler-Nordheim or field-emission transport
- Schottky barrier lowering

Those missing mechanisms matter for high-indium reverse leakage and are likely part of the remaining physics gap.

## R5 Polarization-Charge Stability

- Full polarization, partial polarization scaling, and smoothed polarization interfaces all gave essentially the same reported reverse leakage at `-1.0 V`.
- The internal current proxy also remained qualitatively unchanged.

So polarization treatment is not the dominant source of the reverse-branch failure in the present solver path. It may still influence the electrostatics, but it is not what is fundamentally preventing realistic leakage extraction here.

## Final Interpretation

The combined evidence favors conclusion **B**:

the present solver/contact/current-extraction framework fundamentally limits accurate reverse-bias simulation of high-indium InGaN homojunctions in its current form.

The strongest reasons are:

- terminal current is extracted from a contact-local median region that can remain identically zero while internal current densities become large,
- selective reverse-bias carrier supply is strongly restricted by the implemented boundary conditions,
- refined meshes reveal convergence sensitivity rather than stable physical recovery,
- and the available high-field model set is narrower than what realistic high-indium leakage likely requires.

## Practical Next Step

Before any further calibration or parameter sweeps, the reverse-bias methodology itself should be revised. The next technically meaningful improvement would be:

- replace or supplement the current extraction method with a true terminal-current or integral-flux observable,
- rework the reverse-bias contact model so minority-carrier injection and extraction are physically represented,
- then reassess whether Hurkx-like TAT plus defect recombination is sufficient, or whether additional high-field leakage mechanisms must be implemented.
