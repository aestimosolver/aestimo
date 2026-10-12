# Models and materials

## Physical scope

Aestimo represents a layered heterostructure along one spatial direction.
The effective-mass electron model and valence-band k·p routines calculate
confined states from the band profile; the electrostatic solution couples
charge density and potential. Transport routes add carrier continuity and
drift–diffusion equations.

Changing a well width changes confinement; changing composition changes the
material parameters and band offsets. Doping, contacts and an applied field
can reshape the electrostatic potential. Interpret trends within the selected
model and check numerical convergence before comparing closely spaced energies
or small currents.

## Solver selection in the current source

The table below describes the standard core dispatcher in
`aestimo.run_aestimo`. GUI device workers can add their own characterization
and studies; [GUI workflow selection](gui.md) explains that distinction.

| Scheme | Current route |
| --- | --- |
| 0 | Schrödinger route without the self-consistent Poisson iteration used by scheme 2. |
| 2 | Self-consistent Schrödinger–Poisson route. |
| 7 | Schrödinger–Poisson followed by the standard drift–diffusion route. |
| 8 | Alternative coupled drift–diffusion route. |
| 9 | Gummel/Newton-map drift–diffusion route. |
| 10 | Fully coupled Newton transport selection, with residual and total-current acceptance checks. |
| 11 | Separate multiband k·p route; its exports differ from the standard routes. |

Historical labels 1/3/6 describe nonparabolicity, and 4/5/6 describe exchange.
The current combined electron/hole `Poisson_Schrodinger` routine explicitly
rejects those selections. Their presence in an older example comment is not
evidence that the route works in this source. Use a complete current preset
for the desired workflow rather than changing the scheme number alone.

## Database entries

The current source database contains these keys:

- Materials: AlAs, AlN, AlP, CdO, GaAs, GaN, GaP, Ge, InAs, InN, InP, MgO, Si, ZnO.
- Ternary alloy entries: AlGaAs, AlGaN, AlInN, AlInP, CdZnO, GaAsP, InAsP,
  InGaAs, InGaN, InGaP, MgZnO.
- Quaternary alloy entries: AlGaInN, InGaAsP.

A database entry does not establish that every solver route handles the material
equally or that its parameters have been calibrated for a particular device.
Use the correct crystal system and composition convention, and record any
database edits. See [Inputs and units](inputs.md) and
[Validation](validation.md) for numerical and provenance limitations.
