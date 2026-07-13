from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
RESULTS_ROOT = ROOT / "results"
JSON_PATH = ROOT / "examples" / "untitled_project.json"
LITERATURE_REF_PATH = ROOT / "examples" / "experimental_data" / "ingan_pn_experimental_iv.csv"
FORWARD_RAW_DIR = RESULTS_ROOT / "_tmp_dark_test_output"
REVERSE_RAW_DIR = RESULTS_ROOT / "_tmp_dark_reverse_test_output"
AV_CURR_TO_A_CM2 = 1e-3

REQUIRED_FOLDERS = [
    "band_diagrams",
    "iv_curves",
    "polarization",
    "thickness_sweep",
    "defect_sweep",
    "high_indium",
    "illuminated",
    "design_maps",
]

THERMAL_VOLTAGE_300K = 8.617333262145e-5 * 300.0


@dataclass
class Curve:
    voltage_v: np.ndarray
    current_density_a_cm2: np.ndarray


def ensure_results_layout() -> None:
    RESULTS_ROOT.mkdir(exist_ok=True)
    for name in REQUIRED_FOLDERS:
        (RESULTS_ROOT / name).mkdir(parents=True, exist_ok=True)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required baseline artifact not found: {path}")


def load_authoritative_config(path: Path) -> dict:
    raw = json.loads(path.read_text(encoding="utf-8"))
    layers = []
    for idx, layer in enumerate(raw["layers"], start=1):
        layers.append(
            {
                "index": idx,
                "material": str(layer["material"]),
                "mole": float(layer.get("mole", 0.0)),
                "mole_y": float(layer.get("mole_y", 0.0)),
                "thickness_nm": float(layer["thickness"]),
                "doping_cm3": float(layer["doping"]),
                "doping_type": str(layer["doping_type"]),
                "layer_type": str(layer["type"]),
            }
        )

    return {
        "source_json": str(path.relative_to(ROOT)).replace("\\", "/"),
        "layers": layers,
        "temp_k": float(raw["temp"]),
        "field_kv_cm": float(raw["field"]),
        "bc_left_v": float(raw["bc_left"]),
        "bc_right_v": float(raw["bc_right"]),
        "vmin_v": float(raw["vmin"]),
        "vmax_v": float(raw["vmax"]),
        "vstep_v": float(raw["vstep"]),
        "solver_label": str(raw["solver"]),
        "solver_id": int(str(raw["solver"]).split(":")[0]),
        "grid_step_nm": float(raw["grid_step"]),
        "max_pts": int(raw["max_pts"]),
        "mat_sys": str(raw["mat_sys"]),
        "sub_e": int(raw["sub_e"]),
        "sub_h": int(raw["sub_h"]),
        "quantum_regions": bool(raw.get("Quantum_Regions", False)),
        "quantum_regions_boundary": raw.get("Quantum_Regions_boundary", [[0.0, 0.0]]),
        "tat_field_v_m": float(raw["tat_field"]),
        "project_exp_file": str(raw.get("exp_file", "")),
        "area_cm2": float(raw["area"]),
        "rs_ohm": float(raw["rs"]),
        "rs_mode": str(raw.get("rs_mode", "")),
        "rsh_ohm": float(raw["rsh"]),
        "device_type": str(raw.get("device_type", "")),
        "g_optical_cm3_s": float(raw.get("G_optical", 0.0)),
        "graded_junction": bool(raw.get("graded_junc", False)),
        "diffusion_len_nm": float(raw.get("diffusion_len", 0.0)),
        "enable_polarization": bool(raw.get("enable_polarization", raw.get("polarization", True))),
        "active_models": {
            "poisson_schrodinger": True,
            "drift_diffusion": True,
            "quantum_regions": bool(raw.get("Quantum_Regions", False)),
            "polarization": bool(raw.get("enable_polarization", raw.get("polarization", True))),
            "optical_generation_project_setting": float(raw.get("G_optical", 0.0)) > 0.0,
            "srh_recombination": True,
            "auger_recombination": True,
            "hurkx_like_tat_enhancement": float(raw["tat_field"]) < 1.0e9,
            "field_dependent_mobility": True,
        },
    }


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_csv(path: Path, header: list[str], data: np.ndarray) -> None:
    np.savetxt(path, data, delimiter=",", header=",".join(header), comments="", fmt="%.9e")


def format_sci(value: float) -> str:
    return f"{value:.3e}"


def summarize_device(config: dict) -> dict:
    total_thickness_nm = sum(layer["thickness_nm"] for layer in config["layers"])
    p_layers = [layer for layer in config["layers"] if layer["doping_type"].lower() == "p"]
    n_layers = [layer for layer in config["layers"] if layer["doping_type"].lower() == "n"]
    compositions = sorted({layer["mole"] for layer in config["layers"]})
    return {
        "source_json": config["source_json"],
        "total_layers": len(config["layers"]),
        "total_thickness_nm": total_thickness_nm,
        "composition_x_values": compositions,
        "p_layers": p_layers,
        "n_layers": n_layers,
        "layers": config["layers"],
        "mesh": {
            "grid_step_nm": config["grid_step_nm"],
            "max_pts": config["max_pts"],
            "estimated_grid_points": int(round(total_thickness_nm / config["grid_step_nm"])),
        },
        "contacts": {
            "surface_boundary_left_v": config["bc_left_v"],
            "surface_boundary_right_v": config["bc_right_v"],
            "device_area_cm2": config["area_cm2"],
            "series_resistance_ohm": config["rs_ohm"],
            "shunt_resistance_ohm": config["rsh_ohm"],
        },
        "solver": {
            "solver_label": config["solver_label"],
            "solver_id": config["solver_id"],
            "temperature_k": config["temp_k"],
            "voltage_window_project_v": [config["vmin_v"], config["vmax_v"]],
            "voltage_step_v": config["vstep_v"],
            "material_system": config["mat_sys"],
            "subbands_electron": config["sub_e"],
            "subbands_hole": config["sub_h"],
        },
        "active_models": config["active_models"],
    }


def write_device_summary(config: dict, summary: dict) -> None:
    summary_json = RESULTS_ROOT / "device_summary.json"
    summary_md = RESULTS_ROOT / "device_summary.md"
    write_json(summary_json, summary)

    lines = [
        "# Device Summary",
        "",
        f"Source JSON: `{config['source_json']}`",
        f"Solver: `{config['solver_label']}`",
        f"Material system: `{config['mat_sys']}`",
        f"Temperature: {config['temp_k']:.1f} K",
        f"Total thickness: {summary['total_thickness_nm']:.1f} nm",
        f"Mesh step: {config['grid_step_nm']:.2f} nm",
        "",
        "## Layer Stack",
        "",
        "| Layer | Material | Composition x | Thickness (nm) | Doping Type | Doping (cm^-3) | Type |",
        "| --- | --- | ---: | ---: | --- | ---: | --- |",
    ]
    for layer in config["layers"]:
        lines.append(
            "| {index} | {material} | {mole:.2f} | {thickness_nm:.1f} | {doping_type} | {doping} | {layer_type} |".format(
                index=layer["index"],
                material=layer["material"],
                mole=layer["mole"],
                thickness_nm=layer["thickness_nm"],
                doping_type=layer["doping_type"],
                doping=format_sci(layer["doping_cm3"]),
                layer_type=layer["layer_type"],
            )
        )

    lines.extend(
        [
            "",
            "## Active Models Enabled",
            "",
            f"- Poisson-Schrodinger core: `{summary['active_models']['poisson_schrodinger']}`",
            f"- Drift-diffusion transport: `{summary['active_models']['drift_diffusion']}`",
            f"- Polarization: `{summary['active_models']['polarization']}`",
            f"- Quantum regions: `{summary['active_models']['quantum_regions']}`",
            f"- SRH recombination: `{summary['active_models']['srh_recombination']}`",
            f"- Auger recombination: `{summary['active_models']['auger_recombination']}`",
            f"- TAT enhancement active: `{summary['active_models']['hurkx_like_tat_enhancement']}`",
            f"- Optical generation in project JSON: `{summary['active_models']['optical_generation_project_setting']}`",
            "",
            "## Boundary Conditions",
            "",
            f"- Left surface boundary: {config['bc_left_v']:.3f} V",
            f"- Right surface boundary: {config['bc_right_v']:.3f} V",
            f"- Project voltage window: {config['vmin_v']:.2f} V to {config['vmax_v']:.2f} V in {config['vstep_v']:.2f} V steps",
            "",
            "## Baseline Simulation Adjustments",
            "",
            f"- Dark baseline only: `G_optical` forced from {config['g_optical_cm3_s']:.3e} to `0.0 cm^-3 s^-1`",
            "- Reverse leakage diagnostic only: voltage window extended to `-1.00 V to 0.00 V` with the device unchanged",
        ]
    )
    summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_band_bundle(raw_dir: Path) -> dict[str, np.ndarray]:
    require(raw_dir / "potn_eh_equi_cond.dat")
    require(raw_dir / "np_data0_equi_cond.dat")
    require(raw_dir / "sigma_eh_equi_cond.dat")
    require(raw_dir / "efield_eh_equi_cond.dat")
    potn = np.loadtxt(raw_dir / "potn_eh_equi_cond.dat")
    carriers = np.loadtxt(raw_dir / "np_data0_equi_cond.dat")
    sigma = np.loadtxt(raw_dir / "sigma_eh_equi_cond.dat")
    efield = np.loadtxt(raw_dir / "efield_eh_equi_cond.dat")
    return {
        "x_m": potn[:, 0],
        "ec_eV": potn[:, 1],
        "ev_eV": potn[:, 2],
        "n_cm3": carriers[:, 1],
        "p_cm3": carriers[:, 2],
        "sigma_c_m3": sigma[:, 1],
        "efield_1_v_m": efield[:, 1],
        "efield_2_v_m": efield[:, 2],
    }


def save_band_outputs(bundle: dict[str, np.ndarray]) -> dict:
    out_dir = RESULTS_ROOT / "band_diagrams"
    x_nm = bundle["x_m"] * 1e9
    band_csv = np.column_stack(
        [
            x_nm,
            bundle["ec_eV"],
            bundle["ev_eV"],
            bundle["n_cm3"],
            bundle["p_cm3"],
            bundle["sigma_c_m3"] * 1e-6,
            bundle["efield_1_v_m"] * 1e-5,
            bundle["efield_2_v_m"] * 1e-5,
        ]
    )
    write_csv(
        out_dir / "baseline_equilibrium_band_diagram.csv",
        [
            "x_nm",
            "Ec_eV",
            "Ev_eV",
            "electron_density_cm3",
            "hole_density_cm3",
            "charge_density_C_cm3",
            "electric_field_1_V_cm",
            "electric_field_2_V_cm",
        ],
        band_csv,
    )

    metrics = {
        "device_length_nm": float(x_nm[-1] - x_nm[0]),
        "ec_drop_eV": float(bundle["ec_eV"][0] - bundle["ec_eV"][-1]),
        "ev_drop_eV": float(bundle["ev_eV"][0] - bundle["ev_eV"][-1]),
        "peak_abs_field_mv_cm": float(
            max(np.max(np.abs(bundle["efield_1_v_m"])), np.max(np.abs(bundle["efield_2_v_m"]))) * 1e-5 * 1e3
        ),
        "peak_abs_charge_density_C_cm3": float(np.max(np.abs(bundle["sigma_c_m3"])) * 1e-6),
        "max_electron_density_cm3": float(np.max(bundle["n_cm3"])),
        "max_hole_density_cm3": float(np.max(bundle["p_cm3"])),
    }
    write_json(out_dir / "baseline_equilibrium_band_diagram_metrics.json", metrics)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    axes[0].plot(x_nm, bundle["ec_eV"], label="Ec", color="#006d77", linewidth=2.0)
    axes[0].plot(x_nm, bundle["ev_eV"], label="Ev", color="#ae2012", linewidth=2.0)
    axes[0].set_ylabel("Energy (eV)")
    axes[0].set_title("Baseline Equilibrium Band Diagram")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].semilogy(x_nm, np.clip(bundle["n_cm3"], 1.0, None), label="n", color="#1d3557", linewidth=1.8)
    axes[1].semilogy(x_nm, np.clip(bundle["p_cm3"], 1.0, None), label="p", color="#d62828", linewidth=1.8)
    axes[1].set_xlabel("Position (nm)")
    axes[1].set_ylabel("Carrier Density (cm$^{-3}$)")
    axes[1].grid(alpha=0.25, which="both")
    axes[1].legend()

    fig.savefig(out_dir / "baseline_equilibrium_band_diagram.png", dpi=220)
    plt.close(fig)
    return metrics


def load_curve(path: Path, area_cm2: float, sign_from_voltage: bool = False) -> Curve:
    require(path)
    raw = np.loadtxt(path)
    voltage_v = raw[:, 0]
    current_density_a_cm2 = raw[:, 1] * AV_CURR_TO_A_CM2
    if sign_from_voltage:
        sign = np.where(voltage_v < 0.0, -1.0, 1.0)
        current_density_a_cm2 = sign * np.abs(current_density_a_cm2)
    return Curve(voltage_v=voltage_v, current_density_a_cm2=current_density_a_cm2)


def load_literature_reference() -> Curve:
    voltage = []
    current_density = []
    for line in LITERATURE_REF_PATH.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        v_str, i_str = line.split(",")
        voltage.append(float(v_str))
        current_density.append(float(i_str) / 5.0e-4)
    return Curve(voltage_v=np.asarray(voltage), current_density_a_cm2=np.asarray(current_density))


def interpolate_voltage_at_current(curve: Curve, target_a_cm2: float) -> float | None:
    j = np.abs(curve.current_density_a_cm2)
    if not np.any(j >= target_a_cm2):
        return None
    idx = int(np.argmax(j >= target_a_cm2))
    if idx == 0:
        return float(curve.voltage_v[0])
    v0, v1 = curve.voltage_v[idx - 1], curve.voltage_v[idx]
    j0, j1 = j[idx - 1], j[idx]
    if math.isclose(j0, j1):
        return float(v1)
    alpha = (target_a_cm2 - j0) / (j1 - j0)
    return float(v0 + alpha * (v1 - v0))


def interpolate_current_at_voltage(curve: Curve, voltage_v: float) -> float:
    return float(np.interp(voltage_v, curve.voltage_v, curve.current_density_a_cm2))


def fit_ideality(curve: Curve, vmin: float, vmax: float) -> dict:
    mask = (curve.voltage_v >= vmin) & (curve.voltage_v <= vmax) & (curve.current_density_a_cm2 > 0.0)
    voltage = curve.voltage_v[mask]
    current = curve.current_density_a_cm2[mask]
    if len(voltage) < 3:
        return {"window_v": [vmin, vmax], "points": int(len(voltage)), "ideality_factor": None}
    slope, intercept = np.polyfit(voltage, np.log(current), 1)
    ideality = 1.0 / (slope * THERMAL_VOLTAGE_300K)
    return {
        "window_v": [vmin, vmax],
        "points": int(len(voltage)),
        "ideality_factor": float(ideality),
        "ln_j_slope_per_v": float(slope),
        "ln_j_intercept": float(intercept),
    }


def local_ideality_trace(curve: Curve) -> np.ndarray:
    positive = np.abs(curve.current_density_a_cm2) > 0.0
    voltage = curve.voltage_v[positive]
    current = np.abs(curve.current_density_a_cm2[positive])
    if len(voltage) < 3:
        return np.zeros((0, 2))
    dlnj_dv = np.gradient(np.log(current), voltage)
    ideality = 1.0 / np.maximum(dlnj_dv * THERMAL_VOLTAGE_300K, 1.0e-30)
    return np.column_stack([voltage, ideality])


def write_curve_csv(path: Path, curve: Curve) -> None:
    write_csv(
        path,
        ["voltage_v", "current_density_a_cm2", "current_density_abs_a_cm2"],
        np.column_stack([curve.voltage_v, curve.current_density_a_cm2, np.abs(curve.current_density_a_cm2)]),
    )


def evaluate_baseline(forward: Curve, reverse: Curve, reference: Curve) -> dict:
    sim_ideality = fit_ideality(forward, 0.10, 0.30)
    ref_ideality = fit_ideality(
        Curve(
            voltage_v=reference.voltage_v[reference.voltage_v >= 0.0],
            current_density_a_cm2=np.abs(reference.current_density_a_cm2[reference.voltage_v >= 0.0]),
        ),
        0.10,
        0.30,
    )

    sim_turn_on = interpolate_voltage_at_current(forward, 5.0e-1)
    ref_turn_on = interpolate_voltage_at_current(
        Curve(
            voltage_v=reference.voltage_v[reference.voltage_v >= 0.0],
            current_density_a_cm2=np.abs(reference.current_density_a_cm2[reference.voltage_v >= 0.0]),
        ),
        5.0e-1,
    )
    sim_leakage = abs(interpolate_current_at_voltage(reverse, -0.95))
    ref_leakage = abs(interpolate_current_at_voltage(reference, -0.95))

    leakage_orders = None
    if sim_leakage > 0.0 and ref_leakage > 0.0:
        leakage_orders = float(math.log10(sim_leakage / ref_leakage))

    turn_on_delta = None if sim_turn_on is None or ref_turn_on is None else float(sim_turn_on - ref_turn_on)
    ideality_delta = None
    if sim_ideality["ideality_factor"] is not None and ref_ideality["ideality_factor"] is not None:
        ideality_delta = float(sim_ideality["ideality_factor"] - ref_ideality["ideality_factor"])

    physically_reasonable = (
        sim_turn_on is not None
        and ref_turn_on is not None
        and turn_on_delta is not None
        and abs(turn_on_delta) <= 0.30
        and ideality_delta is not None
        and abs(ideality_delta) <= 0.50
        and leakage_orders is not None
        and abs(leakage_orders) <= 2.0
    )

    return {
        "method_notes": {
            "turn_on_definition": "Voltage where |J| first reaches 0.5 A/cm^2 on the forward dark curve.",
            "ideality_fit_window": "Linear fit of ln(J) vs V from 0.10 V to 0.30 V.",
            "leakage_definition": "Absolute reverse current density at -0.95 V; -1.00 V point excluded because the first raw reverse step initializes to zero in av_curr.dat.",
            "reasonableness_heuristics": {
                "max_turn_on_delta_v": 0.30,
                "max_ideality_delta": 0.50,
                "max_leakage_gap_decades": 2.0,
            },
        },
        "simulation": {
            "forward_current_density_at_0p60V_a_cm2": interpolate_current_at_voltage(forward, 0.60),
            "forward_current_density_at_0p70V_a_cm2": interpolate_current_at_voltage(forward, 0.70),
            "turn_on_voltage_at_0p5A_cm2_v": sim_turn_on,
            "turn_on_voltage_at_0p1A_cm2_v": interpolate_voltage_at_current(forward, 1.0e-1),
            "turn_on_voltage_at_0p01A_cm2_v": interpolate_voltage_at_current(forward, 1.0e-2),
            "reverse_leakage_density_at_minus_0p95V_a_cm2": sim_leakage,
            "ideality_fit_0p10_to_0p30_v": sim_ideality,
        },
        "literature_reference": {
            "source": str(LITERATURE_REF_PATH.relative_to(ROOT)).replace("\\", "/"),
            "forward_current_density_at_0p60V_a_cm2": interpolate_current_at_voltage(reference, 0.60),
            "forward_current_density_at_0p70V_a_cm2": interpolate_current_at_voltage(reference, 0.70),
            "turn_on_voltage_at_0p5A_cm2_v": ref_turn_on,
            "turn_on_voltage_at_0p1A_cm2_v": interpolate_voltage_at_current(reference, 1.0e-1),
            "turn_on_voltage_at_0p01A_cm2_v": interpolate_voltage_at_current(reference, 1.0e-2),
            "reverse_leakage_density_at_minus_0p95V_a_cm2": ref_leakage,
            "ideality_fit_0p10_to_0p30_v": ref_ideality,
        },
        "comparison": {
            "turn_on_delta_v": turn_on_delta,
            "ideality_delta": ideality_delta,
            "leakage_gap_decades_log10_sim_over_ref": leakage_orders,
            "baseline_physically_reasonable": physically_reasonable,
        },
    }


def save_iv_outputs(forward: Curve, reverse: Curve, reference: Curve, metrics: dict) -> None:
    out_dir = RESULTS_ROOT / "iv_curves"
    write_curve_csv(out_dir / "baseline_dark_iv_forward.csv", forward)
    write_curve_csv(out_dir / "baseline_dark_iv_reverse.csv", reverse)

    combined_voltage = np.concatenate([reverse.voltage_v[:-1], forward.voltage_v])
    combined_current = np.concatenate([reverse.current_density_a_cm2[:-1], forward.current_density_a_cm2])
    combined_curve = Curve(combined_voltage, combined_current)
    write_curve_csv(out_dir / "baseline_dark_iv_combined.csv", combined_curve)

    ideality_trace = local_ideality_trace(forward)
    write_csv(
        out_dir / "baseline_dark_iv_ideality_trace.csv",
        ["voltage_v", "local_ideality_factor"],
        ideality_trace,
    )

    reference_csv = np.column_stack(
        [
            reference.voltage_v,
            reference.current_density_a_cm2,
            np.abs(reference.current_density_a_cm2),
        ]
    )
    write_csv(
        out_dir / "baseline_dark_iv_literature_reference.csv",
        ["voltage_v", "current_density_a_cm2", "current_density_abs_a_cm2"],
        reference_csv,
    )

    write_json(out_dir / "baseline_dark_iv_metrics.json", metrics)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
    axes[0].plot(forward.voltage_v, forward.current_density_a_cm2, color="#005f73", linewidth=2.0)
    axes[0].set_title("Baseline Dark I-V (Forward Sweep)")
    axes[0].set_xlabel("Voltage (V)")
    axes[0].set_ylabel("Current Density (A/cm$^2$)")
    axes[0].grid(alpha=0.25)

    axes[1].semilogy(reverse.voltage_v, np.clip(np.abs(reverse.current_density_a_cm2), 1.0e-15, None), label="Reverse", color="#9b2226", linewidth=1.8)
    axes[1].semilogy(forward.voltage_v, np.clip(np.abs(forward.current_density_a_cm2), 1.0e-15, None), label="Forward", color="#0a9396", linewidth=1.8)
    axes[1].set_title("Baseline Dark I-V Magnitude")
    axes[1].set_xlabel("Voltage (V)")
    axes[1].set_ylabel("|J| (A/cm$^2$)")
    axes[1].grid(alpha=0.25, which="both")
    axes[1].legend()
    fig.savefig(out_dir / "baseline_dark_iv.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
    axes[0].semilogy(
        reference.voltage_v,
        np.clip(np.abs(reference.current_density_a_cm2), 1.0e-15, None),
        label="Literature-based reference",
        color="#bb3e03",
        linewidth=2.0,
    )
    axes[0].semilogy(
        combined_curve.voltage_v,
        np.clip(np.abs(combined_curve.current_density_a_cm2), 1.0e-15, None),
        label="Simulation",
        color="#005f73",
        linewidth=1.8,
    )
    axes[0].set_title("Baseline Dark I-V vs High-Indium Literature Trend")
    axes[0].set_xlabel("Voltage (V)")
    axes[0].set_ylabel("|J| (A/cm$^2$)")
    axes[0].grid(alpha=0.25, which="both")
    axes[0].legend()

    if len(ideality_trace) > 0:
        axes[1].plot(ideality_trace[:, 0], ideality_trace[:, 1], color="#005f73", linewidth=1.8, label="Simulation")
    sim_n = metrics["simulation"]["ideality_fit_0p10_to_0p30_v"]["ideality_factor"]
    ref_n = metrics["literature_reference"]["ideality_fit_0p10_to_0p30_v"]["ideality_factor"]
    if sim_n is not None:
        axes[1].axhline(sim_n, color="#005f73", linestyle="--", linewidth=1.2, alpha=0.7)
    if ref_n is not None:
        axes[1].axhline(ref_n, color="#bb3e03", linestyle="--", linewidth=1.2, alpha=0.8, label="Reference fit")
    axes[1].set_title("Forward Ideality Behavior")
    axes[1].set_xlabel("Voltage (V)")
    axes[1].set_ylabel("n")
    axes[1].set_ylim(1.0, 4.5)
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.savefig(out_dir / "baseline_dark_iv_vs_literature.png", dpi=220)
    plt.close(fig)


def write_calibration_report(config: dict, summary: dict, metrics: dict) -> None:
    report_path = RESULTS_ROOT / "baseline_calibration_report.md"
    sim = metrics["simulation"]
    ref = metrics["literature_reference"]
    comp = metrics["comparison"]

    lines = [
        "# Baseline Calibration Report",
        "",
        f"Authoritative source: `{config['source_json']}`",
        "",
        "## Baseline Structure",
        "",
        f"- Layer stack: {len(summary['layers'])} layer high-indium InGaN homojunction",
        f"- Composition: InGaN with x = {summary['composition_x_values'][0]:.2f}",
        f"- Thickness: {summary['total_thickness_nm']:.1f} nm total",
        f"- Doping: p = {format_sci(summary['p_layers'][0]['doping_cm3'])} cm^-3, n = {format_sci(summary['n_layers'][0]['doping_cm3'])} cm^-3",
        f"- Contacts: left = {config['bc_left_v']:.3f} V, right = {config['bc_right_v']:.3f} V",
        "",
        "## Modified Parameters",
        "",
        f"- Dark baseline: `G_optical` changed from {config['g_optical_cm3_s']:.3e} to `0.0 cm^-3 s^-1`",
        "- Leakage diagnostic only: voltage window changed from `0.00 to 1.60 V` to `-1.00 to 0.00 V`",
        "",
        "## Comparison Against High-Indium Literature Trend",
        "",
        "| Metric | Simulation | Literature-Based Reference | Assessment |",
        "| --- | ---: | ---: | --- |",
        "| Reverse leakage at -0.95 V (A/cm^2) | {sim_leak:.3e} | {ref_leak:.3e} | {leak_assess} |".format(
            sim_leak=sim["reverse_leakage_density_at_minus_0p95V_a_cm2"],
            ref_leak=ref["reverse_leakage_density_at_minus_0p95V_a_cm2"],
            leak_assess="Too low by {:.2f} decades".format(abs(comp["leakage_gap_decades_log10_sim_over_ref"]))
            if comp["leakage_gap_decades_log10_sim_over_ref"] is not None and comp["leakage_gap_decades_log10_sim_over_ref"] < 0.0
            else "Within heuristic envelope",
        ),
        "| Turn-on at |J| = 0.5 A/cm^2 (V) | {sim_turn} | {ref_turn:.3f} | {turn_assess} |".format(
            sim_turn=f'>{config["vmax_v"]:.2f}' if sim["turn_on_voltage_at_0p5A_cm2_v"] is None else f'{sim["turn_on_voltage_at_0p5A_cm2_v"]:.3f}',
            ref_turn=ref["turn_on_voltage_at_0p5A_cm2_v"],
            turn_assess="Too late / too soft" if sim["turn_on_voltage_at_0p5A_cm2_v"] is None or (comp["turn_on_delta_v"] is not None and comp["turn_on_delta_v"] > 0.30) else "Within heuristic envelope",
        ),
        "| Ideality factor fit (0.10 V to 0.30 V) | {sim_n:.3f} | {ref_n:.3f} | {n_assess} |".format(
            sim_n=sim["ideality_fit_0p10_to_0p30_v"]["ideality_factor"],
            ref_n=ref["ideality_fit_0p10_to_0p30_v"]["ideality_factor"],
            n_assess="Too ideal; not defect-dominated enough"
            if comp["ideality_delta"] is not None and comp["ideality_delta"] < -0.50
            else "Within heuristic envelope",
        ),
        "",
        "## Interpretation",
        "",
        "- The equilibrium band diagram is internally consistent, but the dark transport does not match the defect-limited behavior expected for high-indium InGaN homojunction devices.",
        "- Reverse leakage is many orders of magnitude lower than the literature-based experimental trend, which points to insufficient defect-assisted leakage and tunneling in the baseline.",
        "- The forward knee is too delayed, and the fitted ideality factor is closer to diffusion-limited behavior than the SRH-dominated regime typically reported for high-indium material.",
        "",
        "## Validation Decision",
        "",
        f"- Baseline physically reasonable: `{comp['baseline_physically_reasonable']}`",
        "- Parameter sweeps should remain blocked until the baseline is recalibrated against the high-indium reference regime.",
        "",
        "## Reference Basis",
        "",
        f"- Numerical comparison file: `{ref['source']}`",
        "- This repository reference file encodes a literature-based high-indium InGaN p-n homojunction trend with high defect density, elevated leakage, and ideality above 2.",
    ]

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_results_layout()
    require(JSON_PATH)
    require(FORWARD_RAW_DIR / "av_curr.dat")
    require(REVERSE_RAW_DIR / "av_curr.dat")
    require(LITERATURE_REF_PATH)

    config = load_authoritative_config(JSON_PATH)
    summary = summarize_device(config)
    write_device_summary(config, summary)

    band_bundle = load_band_bundle(FORWARD_RAW_DIR)
    save_band_outputs(band_bundle)

    forward = load_curve(FORWARD_RAW_DIR / "av_curr.dat", config["area_cm2"], sign_from_voltage=False)
    reverse = load_curve(REVERSE_RAW_DIR / "av_curr.dat", config["area_cm2"], sign_from_voltage=True)
    reference = load_literature_reference()
    metrics = evaluate_baseline(forward, reverse, reference)
    save_iv_outputs(forward, reverse, reference, metrics)
    write_calibration_report(config, summary, metrics)

    print("Baseline packaging complete.")
    print(f"Device summary: {RESULTS_ROOT / 'device_summary.md'}")
    print(f"Band diagram: {RESULTS_ROOT / 'band_diagrams' / 'baseline_equilibrium_band_diagram.png'}")
    print(f"Dark I-V: {RESULTS_ROOT / 'iv_curves' / 'baseline_dark_iv.png'}")
    print(f"Calibration report: {RESULTS_ROOT / 'baseline_calibration_report.md'}")


if __name__ == "__main__":
    main()
