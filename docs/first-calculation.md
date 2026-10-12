# Your first calculation

Choose the desktop workflow to explore a preset visually, or the CLI workflow
for a small, reproducible Schrödinger–Poisson calculation.

## Desktop: inspect a quantum-well preset

1. Start `aestimo-gui` in your installed environment.
2. Select **Dingle 1975 Confinement (Validation)** from **Example Presets**.
3. Inspect **Structure** and **Solver & Grid**. Keep the preset's settings for
   the first run.
4. Use **Save Project** and name the copy `my_dingle_qw.json`. Keeping `dingle`
   in the filename retains the named comparison suite.
5. Click **Run Simulation** and follow **Console** until it finishes.
6. Open **Plotting & Results** and check the selected **Simulation Run**.

The QW workflow displays band edges and confined-state wavefunctions, transition
information, overlaps and carrier profiles. The Dingle comparison additionally
runs a predefined well-width series. Its comparison curve is a study across
multiple structures; it is not simply the band diagram of the one displayed
layer stack.

As a physical trend, stronger confinement in a narrower well generally raises
electron confinement energies relative to the well edge. In a finite,
self-consistent heterostructure the result also depends on barriers, masses and
charge redistribution. Inspect the actual computed trend and convergence rather
than assigning a universal numerical shift.

The historical word "Validation" remains in this preset's name. Check the
reference status and comparison metrics described in [Validation](validation.md)
before using the result as evidence of experimental agreement.

## CLI: run the barrier-doped single-well example

This workflow uses a source checkout, because the complete Python example is
provided there. From the checkout root, with Aestimo installed:

```bash
aestimo -i examples/sample_1qw_barrierdope_ingaas.py
```

Despite its historical filename, the current layer list contains a **GaAs well
with AlGaAs barriers**. The input selects scheme 2, temperature 300 K and a
1 nm grid. Read the actual material list when identifying the structure.

A successful run creates:

```text
sample_1qw_barrierdope_ingaas_output/
```

Open `aestimo.log` there. Inspect the band profiles in
`potn_eh_equi_cond.dat`, the carrier profiles in `np_data0_equi_cond.dat`, and
available `states_e_QWR*_equi_cond.dat` / `states_h_QWR*_equi_cond.dat` files.
State-file energies are in meV and effective masses are ratios to the free
electron mass. [Outputs](outputs.md) explains the columns.

To display the figures after calculating:

```bash
aestimo -d -i examples/sample_1qw_barrierdope_ingaas.py
```

## Make one controlled change

Keep the original input. Copy it to a new filename, change one physical or
numerical setting, and run the copy. For example, compare `gridfactor = 1.0`
with `gridfactor = 0.5`, keeping the layer structure and doping unchanged.
The smaller spacing increases the number of grid points; compare energies or
another selected observable to assess sensitivity to the mesh.

When changing a layer thickness, also review expressions derived from the total
thickness and grid size. Python inputs execute in order. Do not change the
physical structure and mesh simultaneously if you want to identify which change
caused a difference in the result.

Use a new input stem or separate working directory for each retained run to
avoid replacing earlier output. Keep the input, any database changes, log and
Aestimo version/commit together.
