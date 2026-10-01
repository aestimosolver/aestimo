# Autonomous Agent Directives for Aestimo 1D

## 1. Autonomous Execution & Authority

The AI agent is granted full autonomous execution permissions for this repository. The agent must:
- **Proactively Run Commands**: Directly execute shell and Python commands to test, profile, verify, benchmark, or inspect the codebase without pausing to ask for permission.
- **Proactively Create & Edit Files**: Create new files, refactor existing code, update configuration files, and generate benchmark scripts as needed to accomplish user goals.
- **Autonomous Error Resolution**: When errors, test failures, or numerical divergences occur, immediately investigate the root cause, formulate the fix, apply modifications, and verify the resolution end-to-end before concluding the task.
- **No Unnecessary Prompts**: Do not interrupt the workflow or ask the user for confirmation on routine operations (e.g., executing scripts, installing packages, editing files, running tests, checking logs).

---

## 2. Project Architecture & Technical Context

**Aestimo 1D** is a high-performance 1D Schrödinger-Poisson and Drift-Diffusion semiconductor device simulator with a modern CustomTkinter GUI.

### Key Components:
- `aestimo.py`: Central simulation coordinator, `StructureFrom` model builder, layer grid setup, material assignment, and solver routing.
- `aestimo_gui.py`: CustomTkinter GUI frontend featuring async thread-safe simulation workers (`sim_queue`), interactive layer editor, 5-figure plotting suite, and real-time solar metrics panel.
- `aeslibs/newton_raphson.py`: Mode 10 (`comp_scheme = 10`) Fully-Coupled Newton-Raphson drift-diffusion solver using an exact analytic block-tridiagonal Jacobian and sparse direct LU decomposition (`scipy.sparse.linalg.spsolve`).
- `aeslibs/aestimo_poisson1d.py`: Core 1D Poisson-Schrödinger routines, equilibrium Poisson potential calculators, and Fermi-Dirac integrals.
- `characterize_solar.py`: Pure intrinsic photovoltaic characterization extracting $J_{\text{sc}}$, $V_{\text{oc}}$, $\text{FF}$, $\eta$, $P_{\text{max}}$, $V_{\text{mpp}}$, $J_{\text{mpp}}$ directly from numerical current density profiles.
- `database.py`: Material property parameters (energy gaps, band offsets, mobilities, dielectric constants, effective masses) for zincblende and wurtzite semiconductors.
- `examples/`: Standalone `.json` project configurations and matching `.py` benchmark scripts (e.g., `gaas_tobin1990_benchmark.json` / `.py`).

---

## 3. Core Physics & Development Principles

1. **Strict Physical Intrinsic Integrity**:
   - Never inject analytical superposition shortcuts or unphysical empirical approximations into drift-diffusion solutions.
   - All diode and solar figures of merit must be computed intrinsically from the converged numerical solution.
2. **Boundary Conditions & References**:
   - The equilibrium Poisson surface reference (`surface`) must remain `[0.0, 0.0]` for bulk/ohmic contacts.
   - Contact work functions (e.g., $\Phi_{\text{left}} = 5.2\text{ eV}$, $\Phi_{\text{right}} = 4.1\text{ eV}$) must be handled strictly through `work_function_left` and `work_function_right` in drift-diffusion boundary routines.
3. **Crystal Structure Specifics**:
   - Zincblende materials (GaAs, AlGaAs, Si, InP, InGaAs) do not possess wurtzite spontaneous polarization. Only calculate polarization when `mat_type == 'Wurtzite'`.
4. **GUI Thread Safety**:
   - All heavy simulation runs must execute in background daemon threads and post results back to the main thread via `self.sim_queue`.
   - Always force the non-interactive Matplotlib backend (`matplotlib.use('Agg', force=True)`) inside background worker threads to avoid Tkinter thread deadlocks.
5. **Windows & Character Encoding**:
   - Always specify `encoding='utf-8'` when reading or writing files in Python scripts to ensure UTF-8 compatibility across Windows environments.

---

## 4. Verification Workflow

Before concluding any change:
1. Run headless Python tests or script simulations to verify numerical convergence and error margins.
2. Ensure all bias steps converge cleanly without incomplete convergence warnings.
3. Validate that figures and metric dictionaries are properly formatted and rendered without exceptions.
