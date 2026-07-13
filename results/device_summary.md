# Device Summary

Source JSON: `examples/untitled_project.json`
Solver: `7: SP-Drift Diffusion (Sequential)`
Material system: `Wurtzite`
Temperature: 300.0 K
Total thickness: 100.0 nm
Mesh step: 1.00 nm

## Layer Stack

| Layer | Material | Composition x | Thickness (nm) | Doping Type | Doping (cm^-3) | Type |
| --- | --- | ---: | ---: | --- | ---: | --- |
| 1 | InGaN | 0.57 | 40.0 | p | 1.000e+16 | barrier |
| 2 | InGaN | 0.57 | 60.0 | n | 2.000e+17 | barrier |

## Active Models Enabled

- Poisson-Schrodinger core: `True`
- Drift-diffusion transport: `True`
- Polarization: `True`
- Quantum regions: `False`
- SRH recombination: `True`
- Auger recombination: `True`
- TAT enhancement active: `True`
- Optical generation in project JSON: `True`

## Boundary Conditions

- Left surface boundary: 0.000 V
- Right surface boundary: 0.600 V
- Project voltage window: 0.00 V to 1.60 V in 0.05 V steps

## Baseline Simulation Adjustments

- Dark baseline only: `G_optical` forced from 1.000e+21 to `0.0 cm^-3 s^-1`
- Reverse leakage diagnostic only: voltage window extended to `-1.00 V to 0.00 V` with the device unchanged
