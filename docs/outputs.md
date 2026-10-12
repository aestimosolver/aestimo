# Outputs

The CLI writes to `<input-stem>_output/` under the current working directory.
GUI workflows may select their own directories; check the console and log.
The number and naming of files depend on the solver and output controls.

## Common files in the current exporter

| File | Content |
| --- | --- |
| `aestimo.log` | Calculation messages and diagnostics. |
| `potn_eh_equi_cond.dat` | Position in m; conduction and valence energy profiles in eV. |
| `np_data0_equi_cond.dat` | Position in m; electron and hole densities in cm⁻³. |
| `efield_eh_equi_cond.dat` | Position and electric-field arrays. |
| `sigma_eh_equi_cond.dat` | Position and charge-profile array. |
| `states_e_QWR*_equi_cond.dat` | Electron states for a quantum region; columns described by the file header. |
| `states_h_QWR*_equi_cond.dat` | Hole states for a quantum region; columns described by the file header. |
| `av_curr.dat` | Bias sweep and average current values for drift–diffusion routes. |

Bias-dependent exports also use names such as `potn_eh_0.40.dat` and
`states_e_QWR1_0.40.dat`. Available files depend on the calculation route and
`config.py` flags including `potential_out`, `states_out` and
`Drift_Diffusion_out`.

The old documentation's five-file list (`efield.dat`, `potn.dat`, `sigma.dat`,
`states.dat`, `wavefunctions.dat`) describes earlier output conventions.
Use the current exporter and each file's header when analysing a new result.
Current units and normalization must be checked for the selected device route;
do not infer a unit solely from a filename or a GUI axis label.

When reporting results, retain the input project, material changes, package
version or commit, solver/grid controls and relevant convergence logs.
