# Reverse-Bias Code Findings

## Current Extraction

- Scheme 7 computes the reported `av_curr` from the median of the last 10% of `Jtotal` nodes near the right boundary.
- `av_curr.dat` is written directly from `result.Va_t` and `result.av_curr`; there is no file-output clipping or floor in the plotting layer.
- In the active scheme-7 path, `Jelec` and `Jhole` are converted to `mA/cm^2` before `Jtotal` and `av_curr` are formed. Earlier reverse-validation exports used the historical `A/m^2 -> A/cm^2` factor, so nonzero current magnitudes in the original R1-R5 tables are low by a factor of 10; the exact-zero and step-like behavior is unchanged.

## Boundary Conditions

- In photovoltaic mode, the continuity solver uses selective contacts: majority carrier Dirichlet and minority-carrier zero-flux Neumann conditions.
- In non-photovoltaic mode, the continuity solver applies standard Dirichlet conditions to both carriers at both contacts.
- `bc_left` and `bc_right` enter as electrostatic surface offsets in the equilibrium/contact potential initialization rather than as a separate Schottky injection model.

## High-Field Physics

- Implemented reverse-field enhancement in the active solver family is a Hurkx-like TAT multiplier controlled by `tat_field`, `trap_density_scale`, and `trap_energy_offset_ev`.
- No explicit Poole-Frenkel, band-to-band tunneling, Fowler-Nordheim field emission, or Schottky-barrier-lowering implementation was found in the active drift-diffusion solver path.

## Mesh

- The project path exposed through the JSON and solver input supports a uniform `grid_step`; no native nonuniform local junction-refinement control was found in this workflow.
