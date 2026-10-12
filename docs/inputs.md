# Inputs and units

## Python inputs

A Python input defines attributes consumed by `StructureFrom`. The source
repository's `examples/` directory contains complete inputs, including
`sample_1qw_barrierdope_ingaas.py`.

Important fields include:

| Field | Meaning / input unit |
| --- | --- |
| `T` | Temperature, K. |
| `computation_scheme` | Solver selection; `comp_scheme`, if provided, takes precedence. |
| `gridfactor` | Spatial grid spacing, nm. |
| `maxgridpoints` | Upper limit for grid point count. |
| `mat_type` | Crystal system, e.g. `Zincblende`. |
| `subnumber_e`, `subnumber_h` | Requested electron and hole subband counts. |
| `Fapplied` | Applied electric field, V/m. |
| `vmin`, `vmax`, `Each_Step` | Bias sweep start, end and step, V. |
| `device_area` | Device area, cm². |
| `G_optical` | Optical generation rate, cm⁻³ s⁻¹. |

The modern seven-field `material` layer format is:

```python
# thickness_nm, material, mole, mole_y, doping_cm3, doping_type, layer_type
material = [
    [10.0, "GaAs", 0.0, 0.0, 0.0, "n", "w"],
    [20.0, "AlGaAs", 0.3, 0.0, 1e17, "n", "b"],
]
```

This fragment only illustrates the layer format; start a calculation from a
complete input. `w` denotes a well and `b` a barrier. Alloy fractions and the
second composition field depend on the chosen material. Older notebooks can
use earlier formats and require adaptation.

The solver converts several input quantities to SI internally. Do not change
input units to match internal arrays without checking the loader and model.

## JSON projects

Use **Save Project** in the GUI to obtain the current JSON structure. The CLI
maps `layers` to `material` and translates these common fields:

| JSON field | Solver field |
| --- | --- |
| `temp` | `T` |
| `grid_step` | `gridfactor` |
| `max_pts` | `maxgridpoints` |
| `mat_sys` | `mat_type` |
| `vstep` | `Each_Step` |
| `area` | `device_area` |
| `solver` | Numeric solver index extracted from the selection string. |

Other device-specific options should follow a preset for the corresponding
workflow. These tables describe common fields, not an exhaustive schema.

## Numerical controls

The historical `config.py` module contains shared convergence, plotting and
output settings. Input attributes also provide controls for some device solvers;
not every control is read from the same location.

For Mode 10, `dd_residual_tolerance`, `dd_current_atol` and `dd_current_rtol`
control residual and spatial total-current acceptance. Read the
[current-conservation audit](mode10-current-conservation.md) before changing
these settings. Tightening a tolerance can cause a calculation to fail rather
than establish mesh convergence.
