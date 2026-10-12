# Desktop GUI

Launch the application in the environment where Aestimo is installed:

```bash
aestimo-gui
```

On Windows, `aestimo-gui-win` starts the same application without a separate
console window. Both launchers require Tk and a graphical desktop.

## Start from a project

Use **Example Presets** in the left sidebar to load a bundled configuration,
or **Load Project** to open a saved JSON project. Loading a project sets the
layer table, physics settings and solver controls together. The initial window
loads the GaAs Tobin preset when that file is available; select the example you
intend to study before running it.

Use **Save Project** to save your configuration under a new name before editing
it. This saves settings, not the material database or a complete results archive.
The name of the saved JSON file becomes the project name used for subsequent
outputs. See [Your first calculation](first-calculation.md) for a guided start.

## Tabs and controls

| Tab | What to review |
| --- | --- |
| **Structure** | Layer order, thickness, material, alloy composition, doping and well/barrier role. The layer diagram and inspector help check the stack along the growth direction. |
| **Physics & Environment** | Temperature, applied field, boundary values, bias sweep, device mode and reference/comparison settings. |
| **Solver & Grid** | Solver Model, Grid Step, Max Grid Points, Material System, subband counts and optional quantum controls. |
| **Plotting & Results** | Figures, available device metrics and the Simulation Run selector. |
| **Console** | Calculation messages, warnings and errors; use Copy All or Save Log... to retain them. |
| **Database** | Material/alloy properties and separate database save/load controls. |

In **Structure**, **+ Add Layer** adds a layer, **+ Substrate/Buffer** inserts a
buffer-style layer, and **Browse Materials** opens the material browser.
**Clear All** removes the current layer stack. Review the displayed units:
layer thickness is in nm and doping concentration is in cm⁻³.

The GUI's **Applied Electric Field** is labelled **kV/cm**. The Python input's
`Fapplied` is in **V/m**. Do not transfer a number between these fields without
converting units: 1 kV/cm = 100,000 V/m.

## Device mode and solver are separate choices

**Simulation Device / Mode** selects a workflow; **Solver Model** selects a
numerical route within applicable workflows. Keep both selections from the
preset for a first calculation.

| Workflow | Typical result |
| --- | --- |
| Generic Diode / LED | Standard layered-structure calculation and available transport outputs. |
| Diode Simulation | Bias-dependent diode curves and reference comparisons. |
| Solar Cell / Photodetector | Illuminated-device results and solar figures of merit. |
| Solar Study | A sequence of dark, illuminated and temperature calculations. |
| LED Simulation / Characterization | LED characteristics, spectra and available efficiency metrics. |
| Laser Diode Simulation / Characterization | Laser characteristics, threshold and available spectral metrics. |
| QW characterization presets | Confined states, transitions, overlaps and quantum-well comparison plots. |

The QW presets load their dedicated workflow even though that choice is not
listed in the ordinary device dropdown. Project names also select some
specialized comparison suites: keep `dingle` or `miller` in a copied benchmark
project's filename when you intend to retain its corresponding width/field study.

The **Enable Advanced QW Solver** checkbox and **Enable Quantum Corrected DD**
control different calculations. Their settings should follow the selected
workflow and a suitable example; enabling both is not a general accuracy upgrade.
[Models and materials](models.md) explains solver restrictions.

## Run and inspect

Click **Run Simulation**. The run button is disabled while a calculation is
active; progress and messages appear in the sidebar and **Console**. Review the
completion status before interpreting plots.

**Plotting & Results → Simulation Run** selects a run from the current session
history. Loading a project may also display previously saved results found on
disk. A plot appearing immediately after loading does not mean a new calculation
was performed. Check the selected run and its project/time information.

The plot toolbar can save a displayed figure. The structure diagram has a
separate export action supporting PNG, SVG and PDF. Preserve the project and
numerical files with exported figures; a figure alone is not the full record.

## Output location

The main GUI workers use `<project-name>_output` beneath the examples directory.
Solar Study also uses separate working-directory outputs for its dark, light
and temperature runs. Loading a JSON from another directory does not by itself
move the GUI workspace there. Check **Console** for the actual paths and see
[Outputs](outputs.md).

## Material database and project files

In **Database**, **Update In-Memory** applies edited properties to the running
application. **Save to File** exports the database, **Load from File** imports
one, and **Reset Defaults** restores the startup values. These are separate from
**Save Project** and **Load Project**. Retain the database file alongside any
project that depends on modified material parameters.

## Installed presets and references

Installed resources are copied to `~/Aestimo/examples`, or to the workspace
chosen with `AESTIMO_WORKSPACE`. Existing files are preserved on upgrades; a
user-edited preset is not automatically replaced by a new bundled preset.
A source checkout uses its own `examples/` directory.

The comparison section is labelled **Experimental Validation** in the GUI,
but the origin of the loaded reference determines what can be concluded.
Read [Validation and reference data](validation.md) for the status labels and
for the difference between numerical agreement, calibration and measurement.
