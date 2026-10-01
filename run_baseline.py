"""
run_baseline.py  –  Phase 1: Baseline validation for In0.57GaN p-n solar cell
Uses solver 9 (DDG Gummel map) which is confirmed stable for this structure.
Dark I-V sweep: 0 V → 1.6 V in 0.05 V steps.
"""

import json
import os
import sys
import shutil
import numpy as np

sys.path.insert(0, os.path.abspath('.'))
import aestimo
from aeslibs.experimental_validation import (
    load_current_from_avcurr,
    apply_parasitic_resistances,
    calculate_ideality_factor,
)

RESULTS = "results"
OUTPUT_DIR = "baseline_dark_output"

def setup_directories():
    for d in ["band_diagrams","iv_curves","polarization","thickness_sweep",
               "defect_sweep","high_indium","illuminated","design_maps"]:
        os.makedirs(os.path.join(RESULTS, d), exist_ok=True)

def build_config(json_path):
    with open(json_path) as f:
        raw = json.load(f)

    # ── Layer list → material matrix ──────────────────────────────────────────
    mat_list = []
    for l in raw.get("layers", []):
        mat_list.append([
            float(l["thickness"]),
            l["material"],
            float(l.get("mole",   0.0)),
            float(l.get("mole_y", 0.0)),
            float(l["doping"]),
            l["doping_type"],
            l["type"],
        ])
    raw["material"] = mat_list
    raw.pop("layers", None)

    # ── Type-cast all numeric strings ─────────────────────────────────────────
    for k in ["temp","field","bc_left","bc_right","vmin","vmax","vstep",
              "grid_step","area","rs","rsh","G_optical","diffusion_len","tat_field",
              "taun0","taup0"]:
        if k in raw:
            try: raw[k] = float(raw[k])
            except: pass

    # ── Standard aliases ──────────────────────────────────────────────────────
    raw["Fapplied"]    = raw.get("field",  0.0)
    raw["Each_Step"]   = raw.get("vstep",  0.05)
    raw["device_area"] = raw.get("area",   1e-4)
    if "mat_sys" in raw: raw["mat_type"] = raw["mat_sys"]
    for alias, key in [("sub_e","subnumber_e"), ("sub_h","subnumber_h"),
                        ("max_pts","maxgridpoints")]:
        if alias in raw:
            raw[key] = int(raw[alias])

    # ── Solver: use scheme 9 (DDG Gummel map) – confirmed convergent ──────────
    raw["comp_scheme"] = 9

    # ── Device mode ───────────────────────────────────────────────────────────
    raw["photovoltaic_mode"] = "Solar Cell" in str(raw.get("device_type",""))

    # ── Disable experimental validation (no exp file in baseline) ────────────
    raw["enable_experimental_validation"] = False

    return raw

# ─────────────────────────────────────────────────────────────────────────────
def print_device_summary(raw_json):
    print("\n" + "="*60)
    print("DEVICE SUMMARY  –  untitled_project.json")
    print("="*60)
    for i, l in enumerate(raw_json.get("layers", []), 1):
        print("  Layer %d: %s In=%s  %s nm  %s cm-3 %s  type=%s" % (
            i, l['material'], l['mole'], l['thickness'],
            l['doping'], l['doping_type'].upper(), l['type']))
    print("\n  Solver   : %s" % raw_json.get('solver'))
    print("  V sweep  : %s -> %s V (step %s V)" % (
          raw_json.get('vmin'), raw_json.get('vmax'), raw_json.get('vstep')))
    print("  Device   : %s" % raw_json.get('device_type'))
    print("  G_optical: %s cm-3 s-1  (dark = 0)" % raw_json.get('G_optical'))
    print("="*60 + "\n")

# ─────────────────────────────────────────────────────────────────────────────
def main():
    setup_directories()

    # ── Read raw JSON for summary (before layer conversion) ───────────────────
    with open("examples/untitled_project.json") as f:
        raw_json = json.load(f)
    print_device_summary(raw_json)

    # ── Build solver config ───────────────────────────────────────────────────
    config = build_config("examples/untitled_project.json")
    config["G_optical"]          = 0.0           # dark simulation
    config["enable_polarization"] = True          # keep physics correct
    config["comp_scheme"]        = 9
    config["__file__"]           = "baseline_dark"

    # ── Reverse Transport Recovery Fixes ──────────────────────────────────────
    # Activate Trap-Assisted Tunneling (TAT) to remove 'flat' reverse current
    config["tat_field"]          = 1e8            # V/m (Realistic for InGaN)
    # Refine mesh to resolve high-field gradients and prevent 'jumps'
    config["grid_step"]          = 0.2            # nm
    # Ensure stability for high-indium structures
    config["damping"]            = 0.2

    print("Running baseline DARK I-V  (solver 9, 0 -> 1.6 V, step 0.05 V) ...")
    print("  Each voltage step ~5 min  ->  expected ~165 min total\n")

    result = aestimo.run_aestimo(config, drawFigures=False, show=False)

    # ── Post-process ──────────────────────────────────────────────────────────
    area = config["device_area"]
    rs   = config.get("rs",  11.7)
    rsh  = config.get("rsh", 5000.0)

    try:
        sim_v0, sim_i0 = load_current_from_avcurr(OUTPUT_DIR, device_area_cm2=area)
        sim_v, sim_i   = apply_parasitic_resistances(sim_v0, sim_i0, Rs=rs, Rsh=rsh)

        # Ideality factor
        v_mid, n_arr = calculate_ideality_factor(sim_v, sim_i, temperature=300.0)
        n_avg = float(np.mean(n_arr)) if len(n_arr) else float("nan")

        # Save CSV
        csv_path = os.path.join(RESULTS, "iv_curves", "baseline_dark_iv.csv")
        np.savetxt(csv_path,
                   np.column_stack((sim_v, sim_i)),
                   delimiter=",", header="Voltage(V),Current(A)", comments="")

        print("\n" + "-"*60)
        print("  Baseline Calibration Report")
        print("-"*60)
        print("  Data points collected : %d" % len(sim_v))
        print("  V range               : %.3f -> %.3f V" % (sim_v.min(), sim_v.max()))
        print("  Average ideality n    : %.2f" % n_avg)
        print("  CSV saved to          : %s" % csv_path)
        print("-"*60 + "\n")

        # Write calibration report
        rpt = os.path.join(RESULTS, "baseline_calibration_report.txt")
        with open(rpt, "w") as f:
            f.write("Baseline Calibration Report\n")
            f.write("="*40 + "\n")
            f.write(f"Structure   : In0.57GaN p-n homojunction\n")
            f.write(f"p-layer     : 40 nm  1e16 cm-3\n")
            f.write(f"n-layer     : 60 nm  2e17 cm-3\n")
            f.write(f"Solver      : 9 (DDG Gummel map)\n")
            f.write(f"Illumination: DARK (G=0)\n")
            f.write(f"Polarization: enabled\n")
            f.write(f"Average η   : {n_avg:.2f}\n")
            f.write(f"V points    : {len(sim_v)}\n")
        print(f"  Report saved to {rpt}")

    except Exception as e:
        print(f"Post-processing error: {e}")
        import traceback; traceback.print_exc()

    # Copy band diagram if present
    for png in [f"cb_vb_0.0.png", f"cb_vb_0.png"]:
        src = os.path.join(OUTPUT_DIR, png)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(RESULTS, "band_diagrams",
                                          "baseline_equilibrium_band_diagram.png"))
            print(f"  Band diagram saved.")
            break

    print("\n✓ Baseline simulation COMPLETE.")

if __name__ == "__main__":
    main()
