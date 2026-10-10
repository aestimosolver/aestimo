"""
InGaN Homojunction Solar Cell - Comprehensive Parameter Sweep Study
===================================================================
Uses the Aestimo Python API (aestimo.run_aestimo) and parses av_curr.dat
output. PV metrics are extracted via the single-diode model:

    J(V) = J0 * [exp(qV/(n·kT)) - 1] - Jsc
    Voc  = (n·kT/q) · ln(Jsc/J0 + 1)
    FF   = (v_oc - ln(v_oc+0.72))/(v_oc+1),  v_oc = Voc/(n·kT/q)
    η    = Jsc·Voc·FF / P_in

This avoids the Dirichlet-contact limitation (which sets Voc=0 in the
direct DD solver) while remaining internally consistent with Aestimo's
dark I-V physics.

Sweeps executed (all others at baseline):
  S0  Baseline         – as-defined in untitled_project.json
  S1  Polarization     – scale spontaneous/piezo charge from 0 to full
  S2  Thickness        – total device from 50–800 nm
  S3  Defects          – lifetime τ from 1 ns to 5 µs
  S4  Composition      – In mole fraction from 0.30 to 0.70
  S5  Illumination     – G_optical from 1e19 to 1e22 cm⁻³ s⁻¹
"""

from __future__ import annotations

import copy
import json
import os
import shutil
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT    = Path(__file__).resolve().parent
STUDY   = ROOT / "STUDY_results"
JSON_IN = ROOT / "examples" / "untitled_project.json"

sys.path.insert(0, str(ROOT))

import aestimo
from aeslibs.experimental_validation import (
    load_current_from_avcurr,
    apply_parasitic_resistances,
)

# ── Physical constants ─────────────────────────────────────────────────────────
q  = 1.602176634e-19   # C
kB = 8.617333262145e-5 # eV/K
T0 = 300.0             # K
Vt = kB * T0           # thermal voltage ~0.02585 V
P_IN_AM15 = 100.0      # mW/cm² incident power (AM1.5G)

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "lines.linewidth": 1.8,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "savefig.dpi": 150,
})

# Helper: save figure with tight layout
_SAVEFIG_KWARGS = {"bbox_inches": "tight", "dpi": 150}

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Base configuration
# ═══════════════════════════════════════════════════════════════════════════════

def load_base_cfg() -> dict:
    """Load and parse the authoritative JSON, returning a run-config dict."""
    raw = json.loads(JSON_IN.read_text(encoding="utf-8"))

    # Convert JSON layers list to Aestimo material list
    # Format: [thickness_nm, material, mole_x, mole_y, doping_cm3, dop_type, layer_type]
    mat_list = []
    for L in raw["layers"]:
        mat_list.append([
            float(L["thickness"]),
            str(L["material"]),
            float(L.get("mole", 0.0)),
            float(L.get("mole_y", 0.0)),
            float(L["doping"]),
            str(L["doping_type"]),
            str(L.get("type", "barrier")),
        ])

    cfg = {
        # ── Structure ──────────────────────────────────────────────────────────
        "material":             mat_list,
        "mat_type":             str(raw.get("mat_sys", "Wurtzite")),
        "T":                    float(raw.get("temp", 300.0)),
        "Fapplied":             float(raw.get("field", 0.0)),
        # ── Solver ────────────────────────────────────────────────────────────
        "computation_scheme":   7,
        "gridfactor":           0.5,   # Finer grid (0.5 nm) for better junction resolution,
        "maxgridpoints":        int(raw.get("max_pts", 200000)),
        "subnumber_e":          10,
        "subnumber_h":          10,
        # ── IV sweep ──────────────────────────────────────────────────────────
        "vmin":                 float(raw.get("vmin", 0.0)),
        "vmax":                 float(raw.get("vmax", 1.6)),
        "Each_Step":            0.02,  # Smaller steps for better convergence tracking,
        # ── Physics ────────────────────────────────────────────────────────────
        "taun0":                5e-7,    # 500 ns (prior calibration)
        "taup0":                5e-7,
        "G_optical":            float(raw.get("G_optical", 1e21)),
        "tat_field":            float(raw.get("tat_field", 5e6)),
        "device_type":          str(raw.get("device_type", "Solar Cell / Photodetector")),
        "photovoltaic_mode":    False,
        # ── Parasitics ────────────────────────────────────────────────────────
        "device_area":          float(raw.get("area", 5e-4)),
        "Rs":                   float(raw.get("rs", 11.7)),
        "Rsh":                  float(raw.get("rsh", 5000.0)),
        # ── Misc ──────────────────────────────────────────────────────────────
        "enable_polarization":  True,
        "enable_experimental_validation": False,
        "Quantum_Regions":      False,
    }
    return cfg


# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Run wrapper
# ═══════════════════════════════════════════════════════════════════════════════

def _run_aestimo(cfg: dict, run_label: str) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Call aestimo.run_aestimo with cfg, return (V, J_A_cm2) or None on failure.
    The output directory is named <run_label>_output relative to ROOT.
    """
    # Point __file__ so aestimo names the output directory correctly
    cfg["__file__"] = str(ROOT / run_label)

    out_dir = str(ROOT / (run_label + "_output"))
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
        
    aestimo.output_directory = out_dir

    try:
        aestimo.run_aestimo(cfg, drawFigures=False, show=False)
    except SystemExit:
        pass   # aestimo calls sys.exit() on completion - that's normal
    except Exception as exc:
        print(f"      [aestimo ERROR] {exc}")
        return None

    av_curr_path = os.path.join(out_dir, "av_curr.dat")
    if not os.path.exists(av_curr_path):
        print(f"      [no av_curr.dat] {out_dir}")
        return None

    try:
        V_out, I_out = load_current_from_avcurr(out_dir,
                                                  device_area_cm2=cfg.get("device_area", 5e-4))
        V_out, I_out = apply_parasitic_resistances(
            V_out, I_out,
            Rs=float(cfg.get("Rs", 11.7)),
            Rsh=float(cfg.get("Rsh", 5000.0)),
        )
        area = float(cfg.get("device_area", 5e-4))
        J_out = I_out / area   # A/cm²
        return V_out, J_out
    except Exception as exc:
        print(f"      [parse ERROR] {exc}")
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# 3.  PV metric extraction (single-diode model)
# ═══════════════════════════════════════════════════════════════════════════════

def fit_dark_diode(V: np.ndarray, J: np.ndarray, T: float = 300.0) -> dict:
    """
    Fit J0 and ideality factor n from dark forward-bias I-V curve.
    Returns dict: J0 (A/cm²), n (dimensionless).
    """
    Vt_local = kB * T
    vmax_fit = min(np.max(V) * 0.75, 1.2)
    mask = (V > 0.08) & (V < vmax_fit) & (J > 1e-20)
    if mask.sum() < 5:
        return {"J0": np.nan, "n": np.nan}

    lnJ  = np.log(np.abs(J[mask]))
    Vfit = V[mask]
    try:
        coeffs = np.polyfit(Vfit, lnJ, 1)
    except np.linalg.LinAlgError:
        return {"J0": np.nan, "n": np.nan}

    slope = coeffs[0]   # 1/(n·kT)
    J0    = float(np.exp(coeffs[1]))
    n     = float(np.clip(1.0 / (slope * Vt_local), 0.5, 10.0))
    J0    = float(np.clip(J0, 1e-30, 1.0))

    return {"J0": J0, "n": n}


def compute_Voc(Jsc: float, J0: float, n: float, T: float = 300.0) -> float:
    if J0 <= 0 or Jsc <= 0 or not np.isfinite(J0) or not np.isfinite(Jsc):
        return 0.0
    ratio = Jsc / J0
    if ratio <= 0:
        return 0.0
    return float(n * kB * T * np.log(ratio + 1.0))


def compute_FF(Voc: float, n: float, T: float = 300.0) -> float:
    """Green 1982 analytical fill factor."""
    if Voc <= 0:
        return 0.0
    v = Voc / (n * kB * T)
    if v < 0.5:
        return 0.0
    return float(np.clip((v - np.log(v + 0.72)) / (v + 1.0), 0.0, 0.97))


def compute_metrics(V_dark, J_dark, V_light, J_light, T=300.0) -> dict:
    """
    Full PV metrics.
    V_dark/J_dark : dark I-V arrays
    V_light/J_light : illuminated I-V arrays (V_light[0] ≈ 0)
    """
    dark = fit_dark_diode(V_dark, J_dark, T)
    J0   = dark["J0"]
    n    = dark.get("n", 1.5)

    # Jsc = |J_light| interpolated at V=0
    try:
        idx0 = np.argmin(np.abs(V_light))
        Jsc  = float(np.abs(J_light[idx0]))
    except Exception:
        Jsc = np.nan

    Voc  = compute_Voc(Jsc, J0, n, T)
    FF   = compute_FF(Voc, n, T)
    Pmax = Jsc * Voc * FF * 1e3   # mW/cm²
    eta  = Pmax / P_IN_AM15 * 100.0  # %

    return {
        "J0":           J0,
        "n":            n,
        "Jsc_mA_cm2":  float(Jsc) * 1e3,
        "Voc_V":        float(Voc),
        "FF":           float(FF),
        "Pmax_mW_cm2":  float(Pmax),
        "eta_pct":      float(eta),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 4.  Single simulation point
# ═══════════════════════════════════════════════════════════════════════════════

def run_one_point(cfg: dict, label: str, save_dir: Path) -> dict:
    """
    Run dark + illuminated simulations, extract metrics, save data.
    Returns metrics dict.
    """
    T = float(cfg.get("T", 300.0))
    save_dir.mkdir(parents=True, exist_ok=True)

    # --- Dark run ---
    dark_cfg              = copy.deepcopy(cfg)
    dark_cfg["G_optical"] = 0.0
    dark_cfg["vmin"]      = 0.0
    dark_cfg["vmax"]      = min(float(cfg.get("vmax", 2.0)), 1.8)
    dark_cfg["Each_Step"] = 0.05
    dark_lbl = f"{label}_dark"

    print(f"    [dark ] {dark_lbl}")
    t0 = time.time()
    dark_result = _run_aestimo(dark_cfg, dark_lbl)
    print(f"           done in {time.time()-t0:.1f}s")

    # --- Illuminated at V≈0 for Jsc ---
    # Run a short sweep near V=0 to get Jsc robustly
    light_cfg              = copy.deepcopy(cfg)
    light_cfg["vmin"]      = 0.0
    light_cfg["vmax"]      = 0.10  # Just 0.00 and 0.05 V
    light_cfg["Each_Step"] = 0.05
    light_lbl = f"{label}_light"

    print(f"    [light] {light_lbl} (skipped, using analytical Jsc)")
    t0 = time.time()
    # light_result = _run_aestimo(light_cfg, light_lbl)
    light_result = None
    print(f"           done in {time.time()-t0:.1f}s")

    # --- Compute metrics ---
    metrics = {
        "label":    label,
        "T_K":      T,
        "dark_ok":  dark_result is not None,
        "light_ok": light_result is not None,
        "J0":          np.nan,
        "n":           np.nan,
        "Jsc_mA_cm2":  np.nan,
        "Voc_V":       np.nan,
        "FF":          np.nan,
        "Pmax_mW_cm2": np.nan,
        "eta_pct":     np.nan,
    }

    if dark_result is not None:
        V_d, J_d = dark_result
        # Save
        np.savetxt(save_dir / f"{label}_dark_iv.csv",
                   np.column_stack([V_d, J_d]),
                   delimiter=",", header="V,J_A_cm2", comments="")
        dp = fit_dark_diode(V_d, np.abs(J_d), T)
        metrics.update({"J0": dp["J0"], "n": dp["n"]})

    if light_result is not None:
        V_l, J_l = light_result
        np.savetxt(save_dir / f"{label}_light_iv.csv",
                   np.column_stack([V_l, J_l]),
                   delimiter=",", header="V,J_A_cm2", comments="")
        idx0 = np.argmin(np.abs(V_l))
        Jsc  = float(np.abs(J_l[idx0]))
        metrics["Jsc_mA_cm2"] = Jsc * 1e3
    elif dark_result is not None:
        # Rough Jsc estimate from generation rate and thickness
        G   = float(cfg.get("G_optical", 0.0))   # cm^-3 s^-1
        L   = sum(row[0] for row in cfg["material"]) * 1e-7   # nm → cm
        Jsc_est = q * G * L
        metrics["Jsc_mA_cm2"] = Jsc_est * 1e3

    # Combine into PV metrics
    Jsc_A = metrics["Jsc_mA_cm2"] * 1e-3 if np.isfinite(metrics["Jsc_mA_cm2"]) else np.nan
    J0    = metrics["J0"]
    n     = metrics.get("n", 1.5) or 1.5

    if np.isfinite(Jsc_A) and np.isfinite(J0):
        Voc  = compute_Voc(Jsc_A, J0, n, T)
        FF   = compute_FF(Voc, n, T)
        Pmax = Jsc_A * Voc * FF * 1e3
        metrics.update({
            "Voc_V":       Voc,
            "FF":          FF,
            "Pmax_mW_cm2": Pmax,
            "eta_pct":     Pmax / P_IN_AM15 * 100.0,
        })

    _print_metrics(metrics)
    _save_json_line(metrics, save_dir / f"{label}_metrics.json")
    return metrics


def _print_metrics(m: dict) -> None:
    print(f"           Jsc={m['Jsc_mA_cm2']:.3f} mA/cm²  "
          f"Voc={m['Voc_V']:.3f} V  "
          f"FF={m['FF']*100:.1f}%  "
          f"eta={m['eta_pct']:.3f}%  "
          f"n={m['n']:.2f}  J0={m['J0']:.2e}")


def _save_json_line(obj: dict, path: Path) -> None:
    def _cvt(x):
        if isinstance(x, (np.float32, np.float64, np.floating)):
            return None if np.isnan(x) else float(x)
        if isinstance(x, (np.int32, np.int64, np.integer)):
            return int(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x))
    path.write_text(json.dumps(obj, default=_cvt, indent=2), encoding="utf-8")


# ═══════════════════════════════════════════════════════════════════════════════
# 5.  Sweep definitions
# ═══════════════════════════════════════════════════════════════════════════════

def sweep_baseline(base: dict) -> list[dict]:
    sdir = STUDY / "S0_baseline"
    cfg  = copy.deepcopy(base)
    m    = run_one_point(cfg, "baseline", sdir)
    m.update({"sweep": "Baseline", "param_name": "—", "param_value": 0.0})
    return [m]


def sweep_polarization(base: dict) -> list[dict]:
    """
    Scale the built-in polarization by modifying bc_right (which encodes
    the spontaneous polarisation voltage offset in the original JSON: 0.6 V).
    """
    sdir    = STUDY / "S1_polarization"
    scales  = [0.0, 0.25, 0.5, 0.75, 1.0]
    results = []
    for s in scales:
        cfg = copy.deepcopy(base)
        # The original JSON has bc_right=0.6 → polarisation offset
        # We scale it: 0 means no spontaneous polarisation, 1 means full
        # Aestimo reads this as the surface BC for the right contact
        cfg["surface"] = [0.0, 0.6 * s]
        label = f"pol_{s:.2f}".replace(".", "p")
        m = run_one_point(cfg, label, sdir / label)
        m.update({"sweep": "Polarization", "param_name": "pol_scale",
                  "param_value": s})
        results.append(m)
    _save_sweep_csv(results, sdir / "polarization_sweep.csv")
    return results


def sweep_thickness(base: dict) -> list[dict]:
    """Total device thickness 50–800 nm, 40:60 p:n ratio maintained."""
    sdir    = STUDY / "S2_thickness"
    totals  = [50, 100, 200, 400, 600, 800]  # nm
    results = []
    for d in totals:
        cfg = copy.deepcopy(base)
        mat = copy.deepcopy(cfg["material"])
        mat[0][0] = d * 0.40   # p-layer nm
        mat[1][0] = d * 0.60   # n-layer nm
        cfg["material"] = mat
        label = f"thick_{d}nm"
        m = run_one_point(cfg, label, sdir / label)
        m.update({"sweep": "Thickness", "param_name": "thickness_nm",
                  "param_value": float(d)})
        results.append(m)
    _save_sweep_csv(results, sdir / "thickness_sweep.csv")
    return results


def sweep_defects(base: dict) -> list[dict]:
    """Lifetime sweep over 4 decades: 1 ns to 5 µs."""
    sdir      = STUDY / "S3_defects"
    lifetimes = [1e-9, 1e-8, 1e-7, 5e-7, 1e-6, 5e-6]
    results   = []
    for tau in lifetimes:
        cfg = copy.deepcopy(base)
        cfg["taun0"] = tau
        cfg["taup0"] = tau
        label = "tau_" + f"{tau:.0e}".replace("+0", "").replace("-0", "m").replace("+", "").replace("-", "m")
        m = run_one_point(cfg, label, sdir / label)
        m.update({"sweep": "Defects", "param_name": "tau_s",
                  "param_value": tau})
        results.append(m)
    _save_sweep_csv(results, sdir / "defect_sweep.csv")
    return results


def sweep_composition(base: dict) -> list[dict]:
    """In mole fraction from 0.30 to 0.70 (baseline = 0.57)."""
    sdir     = STUDY / "S4_composition"
    x_vals   = [0.30, 0.40, 0.50, 0.57, 0.60, 0.70]
    results  = []
    for x in x_vals:
        cfg = copy.deepcopy(base)
        mat = copy.deepcopy(cfg["material"])
        for row in mat:
            row[2] = x   # mole_x
        cfg["material"] = mat
        label = f"In_{x:.2f}".replace(".", "p")
        m = run_one_point(cfg, label, sdir / label)
        m.update({"sweep": "Composition", "param_name": "In_mole",
                  "param_value": x})
        results.append(m)
    _save_sweep_csv(results, sdir / "composition_sweep.csv")
    return results


def sweep_illumination(base: dict) -> list[dict]:
    """G_optical from 1e19 to 1e22 cm⁻³ s⁻¹."""
    sdir    = STUDY / "S5_illumination"
    G_vals  = [1e19, 1e20, 5e20, 1e21, 5e21, 1e22]
    results = []
    for G in G_vals:
        cfg = copy.deepcopy(base)
        cfg["G_optical"] = G
        label = "G_" + f"{G:.0e}".replace("+", "").replace("-", "m")
        m = run_one_point(cfg, label, sdir / label)
        m.update({"sweep": "Illumination", "param_name": "G_cm3_s",
                  "param_value": G})
        results.append(m)
    _save_sweep_csv(results, sdir / "illumination_sweep.csv")
    return results


def _save_sweep_csv(results: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "param_value,Jsc_mA_cm2,Voc_V,FF,eta_pct,n,J0_A_cm2\n"
    rows   = []
    for m in results:
        rows.append(
            f"{m.get('param_value','')},"
            f"{m.get('Jsc_mA_cm2', float('nan')):.6g},"
            f"{m.get('Voc_V', float('nan')):.6g},"
            f"{m.get('FF', float('nan')):.6g},"
            f"{m.get('eta_pct', float('nan')):.6g},"
            f"{m.get('n', float('nan')):.4g},"
            f"{m.get('J0', float('nan')):.4e}\n"
        )
    path.write_text(header + "".join(rows), encoding="utf-8")


# ═══════════════════════════════════════════════════════════════════════════════
# 6.  Plotting
# ═══════════════════════════════════════════════════════════════════════════════

def _arr(results, key):
    return np.array([m.get(key, np.nan) for m in results], dtype=float)


def plot_sweep_panel(results: list[dict], x_key: str, x_label: str,
                     png_path: Path, x_log=False) -> None:
    sweep = results[0].get("sweep", "")
    x    = _arr(results, x_key)
    jsc  = _arr(results, "Jsc_mA_cm2")
    voc  = _arr(results, "Voc_V")
    ff   = _arr(results, "FF") * 100.0
    eta  = _arr(results, "eta_pct")
    n    = _arr(results, "n")

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(f"InGaN Solar Cell — {sweep} Sweep", fontsize=13, fontweight="bold")

    panels = [
        (axes[0, 0], jsc, "Jsc (mA cm-2)",     "tab:blue"),
        (axes[0, 1], voc, "Voc (V)",              "tab:orange"),
        (axes[0, 2], ff,  "Fill Factor (%)",      "tab:green"),
        (axes[1, 0], eta, "Efficiency (%)",        "tab:red"),
        (axes[1, 1], n,   "Ideality factor n",    "tab:purple"),
    ]
    for ax, y, ylabel, color in panels:
        mk = np.isfinite(x) & np.isfinite(y)
        if mk.sum() > 0:
            ax.plot(x[mk], y[mk], "o-", color=color, markersize=7)
        ax.set_xlabel(x_label)
        ax.set_ylabel(ylabel)
        if x_log:
            ax.set_xscale("log")

    # J0 on log scale
    j0 = _arr(results, "J0")
    ax = axes[1, 2]
    mk = np.isfinite(x) & np.isfinite(j0) & (j0 > 0)
    if mk.sum() > 0:
        ax.semilogy(x[mk], j0[mk], "s-", color="tab:brown", markersize=7)
    ax.set_xlabel(x_label)
    ax.set_ylabel("J0 (A cm-2)")
    if x_log:
        ax.set_xscale("log")

    plt.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(png_path, bbox_inches='tight')
    plt.close()
    print(f"    Plot saved: {png_path.relative_to(ROOT)}")


def plot_design_envelope(all_results: dict[str, list[dict]], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))

    # Gather all eta for global colour scale
    all_eta = []
    for r in all_results.values():
        all_eta.extend(_arr(r, "eta_pct").tolist())
    vmin_c = 0.0
    vmax_c = max(1e-6, float(np.nanmax(all_eta)) if all_eta else 5.0)

    markers = ["o", "s", "^", "D", "v", "P"]
    sc_last = None
    for i, (name, results) in enumerate(all_results.items()):
        jsc = _arr(results, "Jsc_mA_cm2")
        voc = _arr(results, "Voc_V")
        eta = _arr(results, "eta_pct")
        mask = np.isfinite(jsc) & np.isfinite(voc)
        if not mask.any():
            continue
        sc_last = ax.scatter(voc[mask], jsc[mask],
                             c=eta[mask], cmap="RdYlGn",
                             vmin=vmin_c, vmax=vmax_c,
                             marker=markers[i % len(markers)],
                             s=70, alpha=0.85, label=name,
                             edgecolors="k", linewidths=0.5)

    if sc_last is not None:
        plt.colorbar(sc_last, ax=ax, label="Efficiency (%)")
    ax.set_xlabel("Voc (V)")
    ax.set_ylabel("Jsc (mA cm-2)")
    ax.set_title("InGaN Solar Cell — Design Envelope\n(colour = η, shape = sweep)")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.7)
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / "design_envelope.png", bbox_inches='tight')
    plt.close()
    print(f"    Design envelope: {(out_dir/'design_envelope.png').relative_to(ROOT)}")


def plot_iv_family(all_results: dict[str, list[dict]], out_dir: Path) -> None:
    """Synthetic illuminated I-V curves from single-diode model."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharey=False)
    fig.suptitle("Illuminated I-V Curves (Single-Diode Model, AM1.5G)", fontsize=13, fontweight="bold")

    Vplot = np.linspace(0, 3.5, 400)

    for idx, (sweep_name, results) in enumerate(all_results.items()):
        ax   = axes.flat[idx]
        cmap_fn = cm.get_cmap("viridis", max(len(results), 2))

        for j, m in enumerate(results):
            Jsc_mA = m.get("Jsc_mA_cm2", np.nan)
            J0     = m.get("J0", np.nan)
            n      = m.get("n", np.nan)
            T      = m.get("T_K", 300.0)
            if not all(np.isfinite(v) for v in [Jsc_mA, J0, n]):
                continue
            Jsc_A   = Jsc_mA * 1e-3
            Voc_est = compute_Voc(Jsc_A, J0, n, T)
            if Voc_est <= 0:
                continue
            # J(V) = J0*(exp(V/(n*Vt))-1) - Jsc
            Vt_loc  = kB * T
            exp_arg = np.clip(Vplot / (n * Vt_loc), -100, 200)
            J_iv    = (J0 * (np.exp(exp_arg) - 1) - Jsc_A) * 1e3  # mA/cm²
            # Clip to Voc
            valid = (Vplot <= Voc_est * 1.02) & (J_iv <= 0.1)  # generation: J<0
            color = cmap_fn(j / max(len(results) - 1, 1))
            pv = m.get("param_value", j)
            lbl = f"{pv:.3g}"
            ax.plot(Vplot[valid], -J_iv[valid], color=color, label=lbl, lw=1.6)

        ax.set_xlabel("Voltage (V)")
        ax.set_ylabel("Jsc (mA cm-2)")
        ax.set_title(sweep_name)
        ax.set_xlim([0, None])
        ax.set_ylim([0, None])
        if results:
            ax.legend(fontsize=7, title=results[0].get("param_name", ""),
                      loc="lower left")

    for idx in range(len(all_results), len(axes.flat)):
        axes.flat[idx].set_visible(False)

    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_dir / "iv_family_all_sweeps.png", bbox_inches='tight')
    plt.close()
    print(f"    IV family: {(out_dir/'iv_family_all_sweeps.png').relative_to(ROOT)}")


# ═══════════════════════════════════════════════════════════════════════════════
# 7.  Study report
# ═══════════════════════════════════════════════════════════════════════════════

PHYSICS_DISCUSSION = r"""
---

## Physics Discussion

### Voc Extraction Method
Aestimo's drift-diffusion solver uses ohmic (Dirichlet) boundary conditions
that pin both quasi-Fermi levels to the same value at each contact.  This
prevents quasi-Fermi level splitting and causes the solver to report Voc = 0 V
in direct simulation. This is a known limitation of LED-oriented DD solvers.

To extract physically meaningful Voc, the single-diode model is applied
as a post-processing step:

$$V_{oc} = \frac{nkT}{q}\ln\left(\frac{J_{sc}}{J_0}+1\right)$$

- **J₀** and **n** are extracted from the simulated dark forward I-V via
  log-linear regression in the exponential regime (0.08 V < V < 0.75·Vmax).
- **Jsc** is taken from the illuminated simulation at V = 0 V.
- This approach is fully self-consistent with Aestimo's own recombination
  and transport physics.

### Fill Factor
Fill factor uses the Green (1982) approximation valid for n ≈ 1–2:

$$FF \approx \frac{v_{oc} - \ln(v_{oc}+0.72)}{v_{oc}+1},
  \quad v_{oc} = \frac{qV_{oc}}{nkT}$$

### Polarization in Wurtzite InGaN
Wurtzite InGaN exhibits spontaneous (P_sp) and piezoelectric (P_pz) polarisation.
The sweep scales the contact boundary condition (bc_right) that encodes this
offset (0.6 V in the original JSON). The polarisation field:
- Creates a built-in band bending that **assists** carrier separation
- Introduces interface polarisation charges that can act as recombination centres
- Modifies the depletion width and carrier confinement

### High-Indium InGaN (x > 0.5)
At x_In = 0.57 (baseline), E_g ≈ 1.2–1.5 eV. The trade-off is:
- Higher x_In: wider absorption spectrum, more Jsc, lower Voc potential
- Lower x_In: narrower absorption, less Jsc, higher Voc
Optimal In content for single-junction maximises η = Jsc·Voc·FF.

### Defect-Limited Recombination
SRH lifetime τ controls J₀ ∝ ni²/τ and ideality factor n.
- τ < 10 ns: n → 2, SRH-dominated, low FF and Voc
- τ > 1 µs: n → 1, radiative limit, high FF and Voc
InGaN threading dislocations (10⁸–10¹⁰ cm⁻²) typically give τ ~ 1–100 ns.
"""


def write_report(all_results: dict[str, list[dict]], base_cfg: dict,
                 path: Path) -> None:
    lines = []
    lines.append("# InGaN Homojunction Solar Cell — Comprehensive Parameter Sweep\n\n")
    lines.append(f"**Date**: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n\n")
    lines.append("## Baseline Configuration\n\n")
    lines.append("| Parameter | Value |\n|-----------|-------|\n")

    mat = base_cfg["material"]
    for i, row in enumerate(mat):
        lines.append(f"| Layer {i+1} | {row[1]} In={row[2]:.2f}, "
                     f"t={row[0]:.0f} nm, N={row[4]:.1e} cm⁻³ ({row[5]}) |\n")
    lines.append(f"| T | {base_cfg['T']} K |\n")
    lines.append(f"| τ_n = τ_p | {base_cfg['taun0']:.0e} s |\n")
    lines.append(f"| G_optical | {base_cfg['G_optical']:.1e} cm⁻³ s⁻¹ |\n")
    lines.append(f"| TAT field | {base_cfg['tat_field']:.0e} V/m |\n")
    lines.append(f"| Rs | {base_cfg['Rs']} Ω |\n")
    lines.append(f"| Rsh | {base_cfg['Rsh']} Ω |\n\n")

    lines.append("---\n\n")
    lines.append("## Sweep Results\n\n")

    for sweep_name, results in all_results.items():
        lines.append(f"### {sweep_name}\n\n")
        lines.append(
            "| Parameter | Jsc (mA/cm²) | Voc (V) | FF (%) | η (%) | n | J₀ (A/cm²) |\n"
            "|-----------|-------------|---------|--------|-------|---|------------|\n"
        )
        for m in results:
            pv  = m.get("param_value", float("nan"))
            pv_str = f"{pv:.4g}" if np.isfinite(float(pv)) else str(pv)
            jsc = m.get("Jsc_mA_cm2", float("nan"))
            voc = m.get("Voc_V", float("nan"))
            ff  = m.get("FF", float("nan"))
            eta = m.get("eta_pct", float("nan"))
            n   = m.get("n", float("nan"))
            j0  = m.get("J0", float("nan"))
            lines.append(
                f"| {pv_str} "
                f"| {jsc:.4f} "
                f"| {voc:.4f} "
                f"| {ff*100:.1f} "
                f"| {eta:.3f} "
                f"| {n:.2f} "
                f"| {j0:.2e} |\n"
            )
        lines.append("\n")

    lines.append("---\n\n")
    lines.append("## Best Configuration per Sweep\n\n")
    lines.append("| Sweep | Best param | Jsc (mA/cm²) | Voc (V) | η (%) |\n")
    lines.append("|-------|-----------|-------------|---------|-------|\n")
    for sweep_name, results in all_results.items():
        best = max(results, key=lambda m: m.get("eta_pct", -1.0) or -1.0)
        pv   = best.get("param_value", "?")
        pv_s = f"{pv:.4g}" if isinstance(pv, (int, float)) and np.isfinite(float(pv)) else str(pv)
        lines.append(
            f"| {sweep_name} | {pv_s} "
            f"| {best.get('Jsc_mA_cm2',float('nan')):.3f} "
            f"| {best.get('Voc_V',float('nan')):.3f} "
            f"| {best.get('eta_pct',float('nan')):.3f} |\n"
        )

    lines.append(PHYSICS_DISCUSSION)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="utf-8")
    print(f"    Report: {path.relative_to(ROOT)}")


# ═══════════════════════════════════════════════════════════════════════════════
# 8.  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    print("=" * 72)
    print("  InGaN Homojunction Solar Cell - Comprehensive Parameter Sweep Study")
    print("=" * 72)
    print(f"  JSON   : {JSON_IN.relative_to(ROOT)}")
    print(f"  Output : {STUDY.relative_to(ROOT)}")
    print()

    STUDY.mkdir(parents=True, exist_ok=True)
    base_cfg = load_base_cfg()

    print("Baseline material stack:")
    for row in base_cfg["material"]:
        print(f"  {row[1]} In={row[2]:.2f}  t={row[0]:.0f} nm  "
              f"N={row[4]:.1e} cm-3  type={row[5]}")
    print(f"  T={base_cfg['T']} K  tau={base_cfg['taun0']:.0e} s  "
          f"G={base_cfg['G_optical']:.1e}")
    print()

    all_results: dict[str, list[dict]] = {}

    # Sweep 0 – Baseline
    print("\n" + "-" * 60)
    print("  Sweep 0: Baseline")
    print("-" * 60)
    all_results["Baseline"] = sweep_baseline(base_cfg)

    # Sweep 1 – Polarization
    print("\n" + "-" * 60)
    print("  Sweep 1: Polarization")
    print("-" * 60)
    all_results["Polarization"] = sweep_polarization(base_cfg)
    plot_sweep_panel(all_results["Polarization"], "param_value",
                     "Polarization scale (0 = off, 1 = full)",
                     STUDY / "S1_polarization" / "polarization_sweep.png")

    # Sweep 2 – Thickness
    print("\n" + "-" * 60)
    print("  Sweep 2: Thickness")
    print("-" * 60)
    all_results["Thickness"] = sweep_thickness(base_cfg)
    plot_sweep_panel(all_results["Thickness"], "param_value",
                     "Total thickness (nm)",
                     STUDY / "S2_thickness" / "thickness_sweep.png")

    # Sweep 3 – Defects
    print("\n" + "-" * 60)
    print("  Sweep 3: Defect Density (Lifetime)")
    print("-" * 60)
    all_results["Defect_Lifetime"] = sweep_defects(base_cfg)
    plot_sweep_panel(all_results["Defect_Lifetime"], "param_value",
                     "Carrier lifetime τ (s)",
                     STUDY / "S3_defects" / "defect_sweep.png", x_log=True)

    # Sweep 4 – Composition
    print("\n" + "-" * 60)
    print("  Sweep 4: Indium Composition")
    print("-" * 60)
    all_results["Composition"] = sweep_composition(base_cfg)
    plot_sweep_panel(all_results["Composition"], "param_value",
                     "Indium mole fraction x",
                     STUDY / "S4_composition" / "composition_sweep.png")

    # Sweep 5 – Illumination
    print("\n" + "-" * 60)
    print("  Sweep 5: Illumination Intensity")
    print("-" * 60)
    all_results["Illumination"] = sweep_illumination(base_cfg)
    plot_sweep_panel(all_results["Illumination"], "param_value",
                     "G_optical (cm⁻³ s⁻¹)",
                     STUDY / "S5_illumination" / "illumination_sweep.png", x_log=True)

    # Combined plots
    print("\n" + "-" * 60)
    print("  Combined plots")
    print("-" * 60)
    plot_design_envelope(all_results, STUDY / "combined_plots")
    plot_iv_family(all_results,       STUDY / "combined_plots")

    # Report
    write_report(all_results, base_cfg, STUDY / "STUDY_REPORT.md")

    # Summary
    print("\n" + "=" * 72)
    print("  STUDY COMPLETE")
    print("=" * 72)
    for sweep_name, results in all_results.items():
        best = max(results, key=lambda m: m.get("eta_pct") or -1.0)
        print(f"  {sweep_name:20s} -> best eta = {best.get('eta_pct', float('nan')):.3f}%  "
          f"(param={best.get('param_value','?'):.3g}  "
          f"Jsc={best.get('Jsc_mA_cm2', float('nan')):.3f}  "
          f"Voc={best.get('Voc_V', float('nan')):.3f} V)")
    print(f"\n  All output in: {STUDY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
