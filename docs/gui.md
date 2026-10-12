# Desktop GUI

Start the installed desktop application with:

```bash
aestimo-gui
```

Windows installations also provide `aestimo-gui-win`, which uses the GUI launcher
without a console window. Tk and a graphical desktop are required.

## Basic workflow

1. Load an existing project with **Load Project**, or edit the layer table.
2. In **Structure**, review layer thicknesses, materials, alloy fractions and doping.
3. In **Physics & Environment**, review temperature, bias and device settings.
4. In **Solver & Grid**, choose the solver and spatial grid for the calculation.
5. Save a JSON project with **Save Project**, then run the calculation.
6. Inspect the results tab and **Console** messages. Keep the saved project with
   the outputs when recording a calculation.

The **Database** tab allows inspection and editing of material properties and
saving/loading database files. A project preset and a material database are
separate inputs; record any database changes used in a study.

## Presets and workspace

Installed presets are copied into a writable workspace as described in
[Installation](installation.md). Existing user files are preserved rather than
replaced on an upgrade. Compare your saved presets with updated bundled presets
when reproducing a calculation from another version.

The GUI includes quantum-well, diode, solar-cell, LED and laser workflows.
A preset is a starting configuration; its name or a plot label does not certify
experimental agreement. Consult [Validation and reference data](validation.md)
for the status of numerical comparisons and source claims.

CLI output paths and GUI workflow output paths can differ. Consult the GUI
console for the actual output directory of the selected workflow.
