# Reverse-Transport Methodology Recommendation

The reverse-bias validation and the R6 terminal-current audit indicate that the current scheme-7 contact/current framework should not be used for quantitative high-indium InGaN reverse-leakage prediction.

Key reasons:

- The saved current units in the active scheme-7 path are best interpreted as `mA/cm^2`, so previous nonzero current magnitudes exported by the validation script are low by a factor of 10; this does not change the zero-current or step-onset failure mode.
- The official `av_curr` observable can remain exactly zero at -1 V while internal current proxies are large.
- Spatial current profiles at -1 V are not conserved well enough to define a robust terminal current from an alternate probe.
- Contact boundary changes and TAT strength change internal currents, but they do not produce a stable, gradual, terminal reverse branch.

Detailed evidence is in `results\reverse_transport\r6_terminal_observable\terminal_current_methodology_report.md`.

Recommended paper treatment:

- Present reverse-bias results as a solver-methodology validation rather than as calibrated device predictions.
- State that quantitative reverse leakage for high-indium InGaN requires a contact/high-field transport extension, such as physically parameterized thermionic-field emission, Poole-Frenkel or field-enhanced SRH, band-to-band tunneling where applicable, and a terminal-current formulation with current-conservation checks.
- Keep later design sweeps restricted to observables that pass validation, or clearly label reverse-leakage trends as qualitative until the transport framework is extended.
