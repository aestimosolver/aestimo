"""
gui_headless_runner.py
======================
Replicates *exactly* what the Aestimo GUI does when a user:
  1. Clicks "Load Project" -> selects test_solar_cell.json
  2. Clicks "RUN SIMULATION"
  3. Views the Results tab

Then performs the same Eff-vs-T sweep that the GUI's "RUN FULL SOLAR STUDY"
button would trigger, and saves a clean table + plot.

Run with:
    python gui_headless_runner.py
"""

import os, sys, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── ensure project root is importable ─────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aestimo
from aestimo import run_aestimo

# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 1 – Load the example JSON (same as GUI's load_project)
# ═══════════════════════════════════════════════════════════════════════════════
JSON_FILE = os.path.join(os.path.dirname(__file__), "test_solar_cell.json")
print(f"\n{'='*60}")
print(f"  GUI Headless Runner")
print(f"  Loading: {JSON_FILE}")
print(f"{'='*60}\n")

with open(JSON_FILE) as f:
    config = json.load(f)

# ═══════════════════════════════════════════════════════════════════════════════
#  Build helpers that mirror the GUI's run_simulation_worker exactly
# ═══════════════════════════════════════════════════════════════════════════════

def build_material_list(config):
    """Mirror GUI's layer extraction (lines 869-885)."""
    material_list = []
    for l in config["layers"]:
        th    = float(l["thickness"])
        mat   = l["material"]
        x     = float(l["mole"])
        y     = float(l.get("mole_y", 0.0))
        dop   = float(l["doping"])
        dtype = l["doping_type"]
        ltype = l["type"][0]
        if dtype == "i":
            dtype = "n"
        material_list.append([th, mat, x, y, dop, dtype, ltype])
    return material_list


def build_dop_profile(material_list, dx_m):
    """Mirror GUI's doping array builder (lines 1001-1012)."""
    tot_thick = sum(r[0] for r in material_list) * 1e-9
    n_max = int(tot_thick / dx_m)
    dop_arr = np.zeros(n_max)
    curr = 0
    for row in material_list:
        th_m  = row[0] * 1e-9
        val   = row[4]
        dtype = row[5]
        if dtype == 'p':
            val = -val
        steps = int(th_m / dx_m)
        end   = min(curr + steps, n_max)
        dop_arr[curr:end] = val * 1e6   # cm⁻³ → m⁻³
        curr = end
    return dop_arr


def make_input_object(config, T_override=None, G_override=None, label="GUI_RUN"):
    """
    Mirror GUI's InputObject (lines 909-1015).
    Adds photovoltaic_mode = True (needed by Scheme 9).
    """
    material_list = build_material_list(config)

    scheme_id  = int(config.get("solver", "9: Coupled Newton (Robust)").split(":")[0])
    grid_step  = float(config.get("grid_step", 5.0))
    max_pts    = int(config.get("max_pts", 200000))
    sub_e      = int(config.get("sub_e", 5))
    sub_h      = int(config.get("sub_h", 5))
    mat_sys    = config.get("mat_type", config.get("mat_sys", "Wurtzite"))

    T          = T_override if T_override is not None else float(config.get("temp", 300))
    F_app      = float(config.get("field", 0.0)) * 1e5

    vmin       = 0.0          # Always start at 0 V for solar I-V
    vmax       = 2.2          # Sweep forward past Voc
    vstep      = 0.1          # 0.1 V steps (coarse – fast)
    bc_left    = float(config.get("bc_left", 0.0))
    bc_right   = float(config.get("bc_right", 0.0))
    G_opt      = G_override if G_override is not None else float(config.get("G_optical", 5e18))
    tat_field  = float(config.get("tat_field", 1e12))
    rs_val     = float(config.get("rs", 0.0))
    area_cm2   = float(config.get("area", 1e-4))

    dx_m = grid_step * 1e-9

    class InputObject:
        pass

    obj = InputObject()
    obj.T_val              = T
    obj.T                  = T
    obj.computation_scheme = scheme_id
    obj.comp_scheme        = scheme_id
    obj.subnumber_h        = sub_h
    obj.subnumber_e        = sub_e
    obj.gridfactor         = grid_step
    obj.maxgridpoints      = max_pts
    obj.mat_type           = mat_sys
    obj.dx                 = dx_m
    obj.material           = material_list
    obj.Fapplied           = F_app
    obj.tat_field          = tat_field
    obj.G_optical          = G_opt
    obj.vmax               = vmax
    obj.vmin               = vmin
    obj.Each_Step          = vstep
    obj.surface            = np.array([bc_left, bc_right])
    obj.Quantum_Regions    = config.get("Quantum_Regions", False)
    obj.Quantum_Regions_boundary = np.zeros((1, 2))
    obj.Rs                 = rs_val
    obj.device_area_m2     = area_cm2 * 1e-4
    obj.device_area        = area_cm2
    obj.photovoltaic_mode  = True
    obj.enable_polarization = config.get("enable_polarization", False)
    obj.surface_recomb     = (0, 0)
    obj.dop_profile        = build_dop_profile(material_list, dx_m)
    obj.__file__           = os.path.abspath(f"GUI_HEADLESS_{label}.py")
    obj.inputfilename      = f"GUI_HEADLESS_{label}"
    return obj


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 2 – Single simulation at 300 K (mirrors "RUN SIMULATION" button)
# ═══════════════════════════════════════════════════════════════════════════════
OUT_DIR = os.path.join(os.path.dirname(__file__), "gui_headless_output")
os.makedirs(OUT_DIR, exist_ok=True)
aestimo.output_directory = OUT_DIR

print("Running single simulation at 300 K (mirrors GUI 'RUN SIMULATION') …")
obj300 = make_input_object(config, T_override=300, label="300K")
_, _, result300, _ = run_aestimo(obj300, drawFigures=False, show=False)

# Load I-V
iv_file = os.path.join(OUT_DIR, "av_curr.dat")
try:
    iv_data = np.loadtxt(iv_file)
    print(f"\n  I-V loaded: {len(iv_data)} voltage points")
    print(f"  V range   : {iv_data[:,0].min():.2f} → {iv_data[:,0].max():.2f} V")
    print(f"  J(0 V)    : {iv_data[0,1]:.4e} A/m²")
except Exception as e:
    print(f"  [WARN] Could not load av_curr.dat: {e}")

# ═══════════════════════════════════════════════════════════════════════════════
#  Helper: extract PV parameters from a raw av_curr.dat array
# ═══════════════════════════════════════════════════════════════════════════════

def extract_pv_params(iv_data, Pin_mW_cm2=100.0):
    """
    iv_data: Nx2  col0=V(V)  col1=J(A/m²)
    Works with a forward-bias sweep (V increases 0 -> Voc+).
    Jsc = |J(V=0)| in mA/cm²
    Voc = V where J=0 in forward bias
    """
    V  = iv_data[:, 0]
    J  = iv_data[:, 1]         # Aestimo now outputs mA/cm² directly in results/files

    # Sort by V ascending
    order = np.argsort(V)
    V = V[order]
    J = J[order]

    # Jsc: magnitude of J at V closest to 0
    idx0 = np.argmin(np.abs(V))
    Jsc  = abs(J[idx0])

    # Voc: first zero crossing in forward bias (V>0 region)
    Voc = 0.0
    fwd = V >= 0
    Vf  = V[fwd]
    Jf  = J[fwd]
    for i in range(len(Jf)-1):
        if Jf[i] * Jf[i+1] <= 0 and Jf[i] != Jf[i+1]:
            Voc = Vf[i] - Jf[i] * (Vf[i+1]-Vf[i]) / (Jf[i+1]-Jf[i])
            break

    if Voc <= 0 or Jsc <= 0:
        return Jsc, Voc, 0.0, 0.0

    # Power  P = V * J  (J is +ve photocurrent in our convention)
    P = Vf * np.abs(Jf)
    mask = (Vf >= 0) & (Vf <= Voc) & np.isfinite(P)
    if mask.sum() < 2:
        return Jsc, Voc, 0.0, 0.0

    Pmax = P[mask].max()
    FF   = Pmax / (Jsc * Voc) * 100
    Eff  = Pmax / Pin_mW_cm2 * 100
    return Jsc, Voc, FF, Eff


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 3 – Dark Simulation & Temperature sweep (mirrors GUI's "Solar Study")
# ═══════════════════════════════════════════════════════════════════════════════
TEMPS = np.arange(200, 525, 25)   # 200 K … 500 K in 25 K steps

print(f"\n{'='*60}")
print(f"  Starting Solar Study: {TEMPS[0]} K -> {TEMPS[-1]} K")
print(f"  (mirrors GUI 'RUN FULL SOLAR STUDY')")
print(f"{'='*60}\n")

# -- 3.1 Dark Case at 300K --
print("  Running DARK Case (300 K) ...", end=" ", flush=True)
out_dark = os.path.join(OUT_DIR, "Dark_300")
os.makedirs(out_dark, exist_ok=True)
aestimo.output_directory = out_dark
try:
    obj_dark = make_input_object(config, T_override=300.0, G_override=0.0, label="Dark_300")
    run_aestimo(obj_dark, drawFigures=False, show=False)
    print("DONE")
except Exception as e:
    print(f"FAILED ({e})")

# -- 3.2 Temperature Sweep (Illuminated) --
sweep_results = []   # [(T, Jsc, Voc, FF, Eff), …]
for T in TEMPS:
    out_t = os.path.join(OUT_DIR, f"T_{int(T)}")
    os.makedirs(out_t, exist_ok=True)
    aestimo.output_directory = out_t

    print(f"  T = {int(T)} K …", end=" ", flush=True)
    try:
        obj = make_input_object(config, T_override=float(T), label=f"{int(T)}K")
        run_aestimo(obj, drawFigures=False, show=False)

        iv_path = os.path.join(out_t, "av_curr.dat")
        iv = np.loadtxt(iv_path)
        Jsc, Voc, FF, Eff = extract_pv_params(iv)
        sweep_results.append((T, Jsc, Voc, FF, Eff))
        print(f"Jsc={Jsc:.2f} mA/cm², Voc={Voc:.3f} V, FF={FF:.1f}%, Eff={Eff:.2f}%")
    except Exception as e:
        print(f"FAILED ({e})")
        sweep_results.append((T, np.nan, np.nan, np.nan, np.nan))

# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 4 – Print clean table
# ═══════════════════════════════════════════════════════════════════════════════
res = np.array(sweep_results)
T_arr, Jsc_arr, Voc_arr, FF_arr, Eff_arr = res.T

print()
print("=" * 65)
print(" GUI Headless Runner — InGaN Solar Cell — Eff vs Temperature")
print(" (test_solar_cell.json, Scheme 9, AM1.5G G_opt)")
print("=" * 65)
print(f"{'Temp (K)':>9} {'Jsc (mA/cm²)':>13} {'Voc (V)':>9} {'FF (%)':>8} {'Eff (%)':>9}")
print("-" * 55)
for row in sweep_results:
    T, Jsc, Voc, FF, Eff = row
    print(f"{T:>9.0f} {Jsc:>13.3f} {Voc:>9.4f} {FF:>8.2f} {Eff:>9.3f}")
print()

# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 5 – Generate the plot  (same style as eff_vs_temp_plot.py)
# ═══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 9))
fig.patch.set_facecolor('#1a1a2e')
gs = GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.4)

def styled_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor('#16213e')
    ax.tick_params(colors='#e0e0e0', labelsize=9)
    for s in ['bottom', 'left']:
        ax.spines[s].set_color('#555577')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title, color='#aab4d4', fontsize=10, fontweight='bold', pad=8)
    ax.set_xlabel(xlabel, color='#c0c8e0', fontsize=9)
    ax.set_ylabel(ylabel, color='#c0c8e0', fontsize=9)
    ax.grid(True, color='#2d3561', linewidth=0.6, linestyle='--')

ax1 = fig.add_subplot(gs[0, :2])
ax1.plot(T_arr, Eff_arr, color='#f7b731', lw=2.5, marker='o', markersize=5, label='Efficiency (%)')
ax1.fill_between(T_arr, Eff_arr, alpha=0.12, color='#f7b731')
ax1.axvline(300, color='#aaaaaa', ls=':', lw=1.3, label='300 K')
idx300 = np.argmin(np.abs(T_arr - 300))
ax1.annotate(f" {Eff_arr[idx300]:.2f}%",
             xy=(300, Eff_arr[idx300]), color='#f7b731', fontsize=10, fontweight='bold')
styled_ax(ax1, 'Power Conversion Efficiency vs Temperature', 'Temperature (K)', 'Efficiency (%)')
ax1.legend(facecolor='#16213e', edgecolor='#555577', labelcolor='#e0e0e0', fontsize=9)

ax2 = fig.add_subplot(gs[0, 2])
ax2.plot(T_arr, Voc_arr, color='#26de81', lw=2.0, marker='s', markersize=4)
styled_ax(ax2, 'Open-Circuit Voltage', 'Temperature (K)', 'Voc (V)')

ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(T_arr, Jsc_arr, color='#45aaf2', lw=2.0, marker='^', markersize=4)
styled_ax(ax3, 'Short-Circuit Current', 'Temperature (K)', 'Jsc (mA/cm²)')

ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(T_arr, FF_arr, color='#fd9644', lw=2.0, marker='D', markersize=4)
styled_ax(ax4, 'Fill Factor', 'Temperature (K)', 'FF (%)')

ax5 = fig.add_subplot(gs[1, 2])
# I-V curve at 300 K
try:
    iv300 = np.loadtxt(os.path.join(OUT_DIR, "T_300", "av_curr.dat"))
    J300 = iv300[:, 1] * 0.1   # mA/cm²
    ax5.plot(iv300[:, 0], -J300, color='#fc5c65', lw=2.0)
    ax5.axhline(0, color='#888', lw=0.8, ls='--')
    styled_ax(ax5, 'I-V Curve at 300 K', 'Voltage (V)', 'J (mA/cm²)')
except Exception:
    styled_ax(ax5, 'I-V Curve at 300 K', 'Voltage (V)', 'J (mA/cm²)')

fig.suptitle(
    'InGaN Solar Cell (test_solar_cell.json) — GUI Headless Run — Eff vs Temperature',
    color='#e8ecff', fontsize=12, fontweight='bold', y=0.98
)

out_png = os.path.join(os.path.dirname(__file__), 'gui_eff_vs_temp.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight', facecolor='#1a1a2e')
print(f"Plot saved to: {out_png}")
print("\nDone — GUI headless run complete.\n")
