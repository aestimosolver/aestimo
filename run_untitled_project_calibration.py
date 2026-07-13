from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aestimo
import config as aestimo_config
from aestimo import run_aestimo
from aeslibs.experimental_validation import apply_parasitic_resistances

from run_untitled_project_baseline import load_authoritative_config


ROOT = Path(__file__).resolve().parent
JSON_PATH = ROOT / "examples" / "untitled_project.json"
CAL_ROOT = ROOT / "results" / "calibration"
LITERATURE_REF_PATH = ROOT / "examples" / "experimental_data" / "ingan_pn_experimental_iv.csv"
AV_CURR_TO_A_CM2 = 1e-3

FORWARD_WINDOW = {"vmin": 0.0, "vmax": 0.80, "vstep": 0.05}
REVERSE_WINDOW = {"vmin": -1.05, "vmax": 0.0, "vstep": 0.05}
PROFILE_FORWARD_BIAS_V = 0.60
PROFILE_REVERSE_BIAS_V = -1.00
THERMAL_VOLTAGE_300K = 8.617333262145e-5 * 300.0


@dataclass
class Curve:
    voltage_v: np.ndarray
    current_density_a_cm2: np.ndarray
    current_a: np.ndarray


@dataclass
class SimulationRecord:
    stage: str
    case_id: str
    label: str
    notes: str
    overrides: dict[str, Any]
    postprocess: dict[str, Any]
    case_dir: Path
    forward_curve: Curve
    reverse_curve: Curve
    metrics: dict[str, Any]


class CalibrationInput:
    def __init__(self, config_dict: dict[str, Any], input_name: str):
        for key, value in config_dict.items():
            setattr(self, key, value)

        self.material = []
        for layer in config_dict["layers"]:
            self.material.append(
                [
                    float(layer["thickness"]),
                    layer["material"],
                    float(layer.get("mole", 0.0)),
                    float(layer.get("mole_y", 0.0)),
                    float(layer.get("doping", 0.0)),
                    layer.get("doping_type", "n"),
                    layer.get("type", "barrier")[:1],
                ]
            )

        self.T = float(config_dict.get("temp", 300.0))
        self.gridfactor = float(config_dict.get("grid_step", 1.0))
        self.maxgridpoints = int(config_dict.get("max_pts", 200000))
        self.mat_type = config_dict.get("mat_sys", "Wurtzite")
        self.subnumber_e = int(config_dict.get("sub_e", 1))
        self.subnumber_h = int(config_dict.get("sub_h", 1))
        self.computation_scheme = int(str(config_dict.get("solver", "7")).split(":")[0])
        self.comp_scheme = self.computation_scheme
        self.Fapplied = float(config_dict.get("field", 0.0)) * 1e5
        self.surface = [
            float(config_dict.get("bc_left", 0.0)),
            float(config_dict.get("bc_right", 0.0)),
        ]
        self.work_function_left = float(config_dict.get("bc_left", 0.0))
        self.work_function_right = float(config_dict.get("bc_right", 0.0))
        self.vmin = float(config_dict.get("vmin", 0.0))
        self.vmax = float(config_dict.get("vmax", 0.0))
        self.Each_Step = float(config_dict.get("vstep", 0.05))
        self.device_area = float(config_dict.get("area", 1.0e-4))
        self.tat_field = float(config_dict.get("tat_field", 1.0e10))
        self.G_optical = float(config_dict.get("G_optical", 0.0))
        self.photovoltaic_mode = bool(config_dict["photovoltaic_mode"]) if "photovoltaic_mode" in config_dict else "Solar Cell" in str(config_dict.get("device_type", ""))
        self.enable_polarization = bool(
            config_dict.get("enable_polarization", config_dict.get("polarization", True))
        )
        self.polarization_smoothing_nm = float(config_dict.get("polarization_smoothing_nm", 0.0))
        self.Quantum_Regions = bool(config_dict.get("Quantum_Regions", False))
        self.Quantum_Regions_boundary = config_dict.get("Quantum_Regions_boundary", [[0.0, 0.0]])
        self.__file__ = str(JSON_PATH.resolve())
        self.inputfilename = input_name


def ensure_layout() -> None:
    CAL_ROOT.mkdir(parents=True, exist_ok=True)


def sanitize_name(value: str) -> str:
    return value.lower().replace(" ", "_").replace("/", "_").replace(".", "p").replace("-", "m")


def deep_copy_config(config: dict[str, Any]) -> dict[str, Any]:
    return json.loads(json.dumps(config))


def merge_case_config(base_config: dict[str, Any], overrides: dict[str, Any], voltage_window: dict[str, float]) -> dict[str, Any]:
    merged = deep_copy_config(base_config)
    merged["G_optical"] = 0.0
    merged["vmin"] = voltage_window["vmin"]
    merged["vmax"] = voltage_window["vmax"]
    merged["vstep"] = voltage_window["vstep"]

    layer_overrides = overrides.get("layer_overrides", {})
    if layer_overrides:
        for layer_idx, layer_patch in layer_overrides.items():
            for key, value in layer_patch.items():
                merged["layers"][layer_idx][key] = value

    for key, value in overrides.items():
        if key == "layer_overrides":
            continue
        merged[key] = value
    return merged


def run_dark_window(
    base_config: dict[str, Any],
    case_id: str,
    case_dir: Path,
    overrides: dict[str, Any],
    voltage_window: dict[str, float],
) -> tuple[Any, Any]:
    sim_config = merge_case_config(base_config, overrides, voltage_window)
    out_dir = case_dir / f"raw_{sanitize_name(case_id)}_{sanitize_name(f'{voltage_window['vmin']}_{voltage_window['vmax']}')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    aestimo.output_directory = str(out_dir)
    aestimo_config.Drift_Diffusion_out = True
    aestimo_config.potential_out = False
    aestimo_config.electricfield_out = False
    aestimo_config.sigma_out = False
    aestimo_config.states_out = False
    aestimo_config.probability_out = False

    input_obj = CalibrationInput(sim_config, input_name=case_id)
    _, model, result, _ = run_aestimo(input_obj, drawFigures=False, show=False)
    return model, result


def curve_from_result(result: Any, area_cm2: float, postprocess: dict[str, Any]) -> Curve:
    voltage_v = np.asarray(result.Va_t, dtype=float)
    raw_current_density_a_cm2 = np.asarray(result.av_curr, dtype=float) * AV_CURR_TO_A_CM2
    current_a = raw_current_density_a_cm2 * area_cm2

    external_rs = float(postprocess.get("external_rs_ohm", 0.0))
    if external_rs > 0.0:
        voltage_v, current_a = apply_parasitic_resistances(voltage_v, current_a, Rs=external_rs, Rsh=1.0e12)

    current_density_a_cm2 = current_a / area_cm2
    return Curve(voltage_v=voltage_v, current_density_a_cm2=current_density_a_cm2, current_a=current_a)


def load_reference_curve(area_cm2: float) -> Curve:
    voltage = []
    current = []
    for line in LITERATURE_REF_PATH.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        v_str, i_str = line.split(",")
        voltage.append(float(v_str))
        current.append(float(i_str))
    current = np.asarray(current, dtype=float)
    return Curve(
        voltage_v=np.asarray(voltage, dtype=float),
        current_density_a_cm2=current / area_cm2,
        current_a=current,
    )


def interpolate_current(curve: Curve, voltage_v: float) -> float:
    return float(np.interp(voltage_v, curve.voltage_v, curve.current_density_a_cm2))


def interpolate_voltage_at_abs_current(curve: Curve, target_a_cm2: float) -> float | None:
    magnitude = np.abs(curve.current_density_a_cm2)
    if not np.any(magnitude >= target_a_cm2):
        return None
    idx = int(np.argmax(magnitude >= target_a_cm2))
    if idx == 0:
        return float(curve.voltage_v[0])
    v0, v1 = curve.voltage_v[idx - 1], curve.voltage_v[idx]
    j0, j1 = magnitude[idx - 1], magnitude[idx]
    if math.isclose(j0, j1):
        return float(v1)
    alpha = (target_a_cm2 - j0) / (j1 - j0)
    return float(v0 + alpha * (v1 - v0))


def fit_ideality(curve: Curve, vmin: float = 0.10, vmax: float = 0.30) -> dict[str, Any]:
    mask = (curve.voltage_v >= vmin) & (curve.voltage_v <= vmax) & (curve.current_density_a_cm2 > 0.0)
    vv = curve.voltage_v[mask]
    jj = curve.current_density_a_cm2[mask]
    if len(vv) < 3:
        return {"window_v": [vmin, vmax], "points": int(len(vv)), "ideality_factor": None}
    slope, intercept = np.polyfit(vv, np.log(jj), 1)
    ideality = 1.0 / (slope * THERMAL_VOLTAGE_300K)
    return {
        "window_v": [vmin, vmax],
        "points": int(len(vv)),
        "ideality_factor": float(ideality),
        "ln_j_slope_per_v": float(slope),
        "ln_j_intercept": float(intercept),
    }


def local_ideality_trace(curve: Curve) -> np.ndarray:
    mask = (curve.voltage_v >= 0.05) & (curve.current_density_a_cm2 > 0.0)
    v = curve.voltage_v[mask]
    j = curve.current_density_a_cm2[mask]
    if len(v) < 3:
        return np.zeros((0, 2))
    dlnj_dv = np.gradient(np.log(j), v)
    nloc = 1.0 / np.maximum(dlnj_dv * THERMAL_VOLTAGE_300K, 1e-30)
    return np.column_stack([v, nloc])


def curvature_metric(local_ideality: np.ndarray) -> float | None:
    if len(local_ideality) < 5:
        return None
    mask = (local_ideality[:, 0] >= 0.2) & (local_ideality[:, 0] <= 0.8)
    if np.count_nonzero(mask) < 4:
        return None
    return float(np.std(local_ideality[mask, 1]))


def score_metrics(metrics: dict[str, Any]) -> float:
    leakage = metrics["reverse_leakage_at_minus_1V_a_cm2"]
    ideality = metrics["ideality_fit"]["ideality_factor"]
    turn_on = metrics["turn_on_0p5Acm2_v"]

    log_leak = math.log10(max(leakage, 1e-30))
    leakage_penalty = 0.0
    if log_leak < -5.0:
        leakage_penalty = (-5.0 - log_leak) ** 2
    elif log_leak > -2.0:
        leakage_penalty = (log_leak + 2.0) ** 2

    ideality_penalty = 10.0
    if ideality is not None:
        if ideality < 2.0:
            ideality_penalty = (2.0 - ideality) ** 2
        elif ideality > 3.0:
            ideality_penalty = (ideality - 3.0) ** 2
        else:
            ideality_penalty = 0.0

    turn_on_penalty = 5.0
    if turn_on is not None:
        if turn_on < 0.5:
            turn_on_penalty = (0.5 - turn_on) ** 2
        elif turn_on > 0.8:
            turn_on_penalty = (turn_on - 0.8) ** 2
        else:
            turn_on_penalty = 0.0
    else:
        turn_on_penalty = 3.0

    curvature = metrics["forward_curvature_std_n"]
    curvature_bonus = 0.0 if curvature is None else -min(curvature, 2.0) * 0.1
    return leakage_penalty + ideality_penalty + turn_on_penalty + curvature_bonus


def result_array(result: Any, base_name: str) -> np.ndarray:
    value = getattr(result, f"{base_name}_", None)
    if value is None:
        value = getattr(result, base_name)
    return np.asarray(value, dtype=float)


def result_bias_slice(result: Any, base_name: str, bias_index: int) -> np.ndarray:
    data = result_array(result, base_name)
    if data.ndim == 1:
        return data
    return np.asarray(data[bias_index, :], dtype=float)


def compute_recombination_components(
    model: Any,
    result: Any,
    bias_index: int,
) -> dict[str, np.ndarray]:
    n_m3 = result_bias_slice(result, "nf_result", bias_index)
    p_m3 = result_bias_slice(result, "pf_result", bias_index)
    ni_m3 = np.asarray(getattr(model, "ni_phys"), dtype=float)
    taun0_s = np.asarray(model.TAUN0, dtype=float)
    taup0_s = np.asarray(model.TAUP0, dtype=float)
    cn0_m6_s = np.asarray(model.Cn0, dtype=float)
    cp0_m6_s = np.asarray(model.Cp0, dtype=float)
    efield_v_m = np.abs(result_bias_slice(result, "el_field1_result", bias_index))

    trap_density_scale = max(getattr(model, "trap_density_scale", 1.0), 1e-12)
    trap_energy_offset_ev = getattr(model, "trap_energy_offset_ev", 0.0)
    trap_arg = np.clip(trap_energy_offset_ev / THERMAL_VOLTAGE_300K, -40.0, 40.0)
    n1 = ni_m3 * np.exp(trap_arg)
    p1 = ni_m3 * np.exp(-trap_arg)

    gamma = np.zeros_like(efield_v_m)
    tat_field = float(getattr(model, "tat_field", 1.0e10))
    if tat_field < 1.0e9:
        mask = efield_v_m > 1.0e4
        ratio = efield_v_m[mask] / tat_field
        gamma[mask] = 2.0 * np.sqrt(3.0 * np.pi) * ratio * np.exp(np.clip(ratio**2, 0.0, 20.0))
    gamma *= trap_density_scale

    np_minus_ni2 = n_m3 * p_m3 - ni_m3**2
    denom = taup0_s * (n_m3 + n1) + taun0_s * (p_m3 + p1)
    denom = np.maximum(denom, 1e-30)
    srh_tat_m3_s = np_minus_ni2 * ((1.0 + gamma) * trap_density_scale) / denom
    auger_m3_s = (cn0_m6_s * n_m3 + cp0_m6_s * p_m3) * np_minus_ni2
    total_m3_s = srh_tat_m3_s + auger_m3_s

    return {
        "srh_tat_cm3_s": srh_tat_m3_s * 1e-6,
        "auger_cm3_s": auger_m3_s * 1e-6,
        "total_cm3_s": total_m3_s * 1e-6,
        "n_cm3": n_m3 * 1e-6,
        "p_cm3": p_m3 * 1e-6,
    }


def find_bias_index(result: Any, target_bias_v: float) -> int:
    voltage = np.asarray(result.Va_t, dtype=float)
    return int(np.argmin(np.abs(voltage - target_bias_v)))


def extract_profile_tables(model: Any, result: Any) -> dict[str, dict[str, np.ndarray]]:
    x_nm = np.asarray(result.xaxis, dtype=float) * 1e9
    profiles: dict[str, dict[str, np.ndarray]] = {}
    for label, target_v in [("forward_0p60V", PROFILE_FORWARD_BIAS_V), ("reverse_m1p00V", PROFILE_REVERSE_BIAS_V)]:
        idx = find_bias_index(result, target_v)
        recomb = compute_recombination_components(model, result, idx)
        profiles[label] = {
            "x_nm": x_nm,
            "electric_field_v_cm": result_bias_slice(result, "el_field1_result", idx) * 1e-2,
            "charge_density_c_cm3": result_bias_slice(result, "ro_result", idx) * 1e-6,
            "srh_tat_cm3_s": recomb["srh_tat_cm3_s"],
            "auger_cm3_s": recomb["auger_cm3_s"],
            "total_recomb_cm3_s": recomb["total_cm3_s"],
            "n_cm3": recomb["n_cm3"],
            "p_cm3": recomb["p_cm3"],
            "bias_v_actual": np.asarray(result.Va_t, dtype=float)[idx],
        }
    return profiles


def extract_combined_profile_tables(
    model_fwd: Any,
    result_fwd: Any,
    model_rev: Any,
    result_rev: Any,
) -> dict[str, dict[str, np.ndarray]]:
    forward_profiles = extract_profile_tables(model_fwd, result_fwd)
    reverse_profiles = extract_profile_tables(model_rev, result_rev)
    return {
        "forward_0p60V": forward_profiles["forward_0p60V"],
        "reverse_m1p00V": reverse_profiles["reverse_m1p00V"],
    }


def evaluate_case_metrics(
    case_id: str,
    stage: str,
    base_config: dict[str, Any],
    forward_curve: Curve,
    reverse_curve: Curve,
    reference_curve: Curve,
    profile_data: dict[str, dict[str, np.ndarray]],
) -> dict[str, Any]:
    local_n = local_ideality_trace(forward_curve)
    ideality = fit_ideality(forward_curve)
    leakage = abs(interpolate_current(reverse_curve, -1.0))
    turn_on = interpolate_voltage_at_abs_current(forward_curve, 0.5)
    turn_on_0p1 = interpolate_voltage_at_abs_current(forward_curve, 0.1)

    contact_forward_profile = profile_data["forward_0p60V"]
    left_layer = base_config["layers"][0]
    right_layer = base_config["layers"][-1]
    left_doping = float(left_layer["doping"])
    right_doping = float(right_layer["doping"])
    left_majority = (
        contact_forward_profile["p_cm3"][0]
        if left_layer["doping_type"].lower() == "p"
        else contact_forward_profile["n_cm3"][0]
    )
    right_majority = (
        contact_forward_profile["p_cm3"][-1]
        if right_layer["doping_type"].lower() == "p"
        else contact_forward_profile["n_cm3"][-1]
    )

    integrated_srh_tat = float(
        np.trapz(np.abs(contact_forward_profile["srh_tat_cm3_s"]), contact_forward_profile["x_nm"] * 1e-7)
    )
    integrated_auger = float(
        np.trapz(np.abs(contact_forward_profile["auger_cm3_s"]), contact_forward_profile["x_nm"] * 1e-7)
    )

    metrics = {
        "stage": stage,
        "case_id": case_id,
        "reverse_leakage_at_minus_1V_a_cm2": leakage,
        "ideality_fit": ideality,
        "turn_on_0p5Acm2_v": turn_on,
        "turn_on_0p1Acm2_v": turn_on_0p1,
        "current_density_at_0p60V_a_cm2": interpolate_current(forward_curve, 0.60),
        "current_density_at_0p70V_a_cm2": interpolate_current(forward_curve, 0.70),
        "forward_curvature_std_n": curvature_metric(local_n),
        "local_ideality_trace_points": int(len(local_n)),
        "left_contact_majority_injection_ratio_at_0p60V": float(left_majority / max(left_doping, 1.0)),
        "right_contact_majority_injection_ratio_at_0p60V": float(right_majority / max(right_doping, 1.0)),
        "integrated_srh_tat_forward_0p60V": integrated_srh_tat,
        "integrated_auger_forward_0p60V": integrated_auger,
        "defect_dominated_forward_0p60V": bool(integrated_srh_tat > 5.0 * max(integrated_auger, 1e-30)),
        "reference": {
            "reverse_leakage_at_minus_1V_a_cm2": abs(interpolate_current(reference_curve, -1.0)),
            "turn_on_0p5Acm2_v": interpolate_voltage_at_abs_current(reference_curve, 0.5),
            "ideality_fit": fit_ideality(
                Curve(
                    voltage_v=reference_curve.voltage_v[reference_curve.voltage_v >= 0.0],
                    current_density_a_cm2=np.abs(
                        reference_curve.current_density_a_cm2[reference_curve.voltage_v >= 0.0]
                    ),
                    current_a=np.abs(reference_curve.current_a[reference_curve.voltage_v >= 0.0]),
                )
            ),
        },
    }
    metrics["score"] = score_metrics(metrics)
    metrics["meets_stop_condition"] = bool(
        1.0e-5 <= leakage <= 1.0e-2
        and ideality["ideality_factor"] is not None
        and 2.0 <= ideality["ideality_factor"] <= 3.0
        and turn_on is not None
        and turn_on < 0.8
        and (metrics["forward_curvature_std_n"] or 0.0) >= 0.25
    )
    return metrics


def save_curve_csv(path: Path, curve: Curve) -> None:
    data = np.column_stack([curve.voltage_v, curve.current_density_a_cm2, curve.current_a])
    np.savetxt(
        path,
        data,
        delimiter=",",
        header="voltage_v,current_density_a_cm2,current_a",
        comments="",
        fmt="%.9e",
    )


def load_curve_csv(path: Path) -> Curve:
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return Curve(
        voltage_v=np.asarray(data[:, 0], dtype=float),
        current_density_a_cm2=np.asarray(data[:, 1], dtype=float),
        current_a=np.asarray(data[:, 2], dtype=float),
    )


def save_profile_csv(path: Path, table: dict[str, np.ndarray]) -> None:
    header = [
        "x_nm",
        "electric_field_v_cm",
        "charge_density_c_cm3",
        "n_cm3",
        "p_cm3",
        "srh_tat_cm3_s",
        "auger_cm3_s",
        "total_recomb_cm3_s",
    ]
    data = np.column_stack([table[key] for key in header])
    np.savetxt(path, data, delimiter=",", header=",".join(header), comments="", fmt="%.9e")


def save_case_outputs(
    record: SimulationRecord,
    model_fwd: Any,
    result_fwd: Any,
    result_rev: Any,
    reference_curve: Curve,
    profile_data: dict[str, dict[str, np.ndarray]],
) -> None:
    record.case_dir.mkdir(parents=True, exist_ok=True)

    effective_case = {
        "source_json": str(JSON_PATH.relative_to(ROOT)).replace("\\", "/"),
        "stage": record.stage,
        "case_id": record.case_id,
        "label": record.label,
        "overrides": record.overrides,
        "postprocess": record.postprocess,
    }
    (record.case_dir / "effective_case.json").write_text(json.dumps(effective_case, indent=2), encoding="utf-8")
    (record.case_dir / "metrics.json").write_text(json.dumps(record.metrics, indent=2), encoding="utf-8")
    save_curve_csv(record.case_dir / "dark_iv_forward.csv", record.forward_curve)
    save_curve_csv(record.case_dir / "dark_iv_reverse.csv", record.reverse_curve)
    save_curve_csv(record.case_dir / "literature_reference.csv", reference_curve)

    for label, table in profile_data.items():
        save_profile_csv(record.case_dir / f"{label}_profile.csv", table)

    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    ax.plot(record.forward_curve.voltage_v, record.forward_curve.current_density_a_cm2, color="#005f73", linewidth=2.0)
    ax.set_title(f"{record.label} Dark I-V")
    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Current Density (A/cm$^2$)")
    ax.grid(alpha=0.25)
    fig.savefig(record.case_dir / "dark_iv_linear.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    ax.semilogy(
        reference_curve.voltage_v,
        np.clip(np.abs(reference_curve.current_density_a_cm2), 1.0e-15, None),
        label="Literature reference",
        color="#bb3e03",
        linewidth=2.0,
    )
    combined_v = np.concatenate([record.reverse_curve.voltage_v[:-1], record.forward_curve.voltage_v])
    combined_j = np.concatenate([record.reverse_curve.current_density_a_cm2[:-1], record.forward_curve.current_density_a_cm2])
    ax.semilogy(combined_v, np.clip(np.abs(combined_j), 1.0e-15, None), label="Simulation", color="#005f73", linewidth=1.8)
    ax.set_title(f"{record.label} Dark I-V (Semi-Log)")
    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("|J| (A/cm$^2$)")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    fig.savefig(record.case_dir / "dark_iv_semilog.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.5), sharex=True, constrained_layout=True)
    for idx, (profile_name, table) in enumerate(profile_data.items()):
        axes[idx].plot(table["x_nm"], table["electric_field_v_cm"], color="#1d3557", linewidth=1.8)
        axes[idx].set_ylabel("Field (V/cm)")
        axes[idx].set_title(f"{profile_name} at {table['bias_v_actual']:.2f} V")
        axes[idx].grid(alpha=0.25)
    axes[-1].set_xlabel("Position (nm)")
    fig.suptitle(f"{record.label} Electric Field Profiles")
    fig.savefig(record.case_dir / "electric_field_profiles.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.5), sharex=True, constrained_layout=True)
    for idx, (profile_name, table) in enumerate(profile_data.items()):
        axes[idx].plot(table["x_nm"], table["srh_tat_cm3_s"], label="SRH+TAT", color="#9b2226", linewidth=1.7)
        axes[idx].plot(table["x_nm"], table["auger_cm3_s"], label="Auger", color="#0a9396", linewidth=1.5)
        axes[idx].plot(table["x_nm"], table["total_recomb_cm3_s"], label="Total", color="#6a4c93", linewidth=1.5, alpha=0.85)
        axes[idx].axhline(0.0, color="black", linewidth=0.8, alpha=0.4)
        axes[idx].set_ylabel("Rate (cm$^{-3}$ s$^{-1}$)")
        axes[idx].set_title(f"{profile_name} at {table['bias_v_actual']:.2f} V")
        axes[idx].grid(alpha=0.25)
        axes[idx].legend()
    axes[-1].set_xlabel("Position (nm)")
    fig.suptitle(f"{record.label} Recombination Profiles")
    fig.savefig(record.case_dir / "recombination_profiles.png", dpi=220)
    plt.close(fig)

    metrics_md = [
        "# Metrics",
        "",
        f"- Reverse leakage at -1.0 V: {record.metrics['reverse_leakage_at_minus_1V_a_cm2']:.3e} A/cm^2",
        f"- Ideality factor (0.10-0.30 V): {record.metrics['ideality_fit']['ideality_factor'] if record.metrics['ideality_fit']['ideality_factor'] is not None else 'N/A'}",
        f"- Turn-on at |J| = 0.5 A/cm^2: {record.metrics['turn_on_0p5Acm2_v'] if record.metrics['turn_on_0p5Acm2_v'] is not None else f'>{np.max(record.forward_curve.voltage_v):.2f}'} V",
        f"- Turn-on at |J| = 0.1 A/cm^2: {record.metrics['turn_on_0p1Acm2_v'] if record.metrics['turn_on_0p1Acm2_v'] is not None else f'>{np.max(record.forward_curve.voltage_v):.2f}'} V",
        f"- Forward curvature std(n): {record.metrics['forward_curvature_std_n']}",
        f"- Defect-dominated forward profile: {record.metrics['defect_dominated_forward_0p60V']}",
        f"- Composite score: {record.metrics['score']:.4f}",
        f"- Meets stop condition: {record.metrics['meets_stop_condition']}",
    ]
    (record.case_dir / "metrics.md").write_text("\n".join(metrics_md) + "\n", encoding="utf-8")
    (record.case_dir / "calibration_notes.md").write_text(record.notes + "\n", encoding="utf-8")


def run_case(
    base_config: dict[str, Any],
    reference_curve: Curve,
    stage: str,
    case_id: str,
    label: str,
    overrides: dict[str, Any],
    notes: str,
    postprocess: dict[str, Any] | None = None,
) -> SimulationRecord:
    if postprocess is None:
        postprocess = {}

    case_dir = CAL_ROOT / stage / case_id
    metrics_path = case_dir / "metrics.json"
    forward_csv_path = case_dir / "dark_iv_forward.csv"
    reverse_csv_path = case_dir / "dark_iv_reverse.csv"
    if metrics_path.exists() and forward_csv_path.exists() and reverse_csv_path.exists():
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        return SimulationRecord(
            stage=stage,
            case_id=case_id,
            label=label,
            notes=notes,
            overrides=overrides,
            postprocess=postprocess,
            case_dir=case_dir,
            forward_curve=load_curve_csv(forward_csv_path),
            reverse_curve=load_curve_csv(reverse_csv_path),
            metrics=metrics,
        )

    model_fwd, result_fwd = run_dark_window(base_config, case_id, case_dir, overrides, FORWARD_WINDOW)
    model_rev, result_rev = run_dark_window(base_config, case_id + "_rev", case_dir, overrides, REVERSE_WINDOW)

    area_cm2 = float(base_config["area_cm2"])
    forward_curve = curve_from_result(result_fwd, area_cm2, postprocess)
    reverse_curve = curve_from_result(result_rev, area_cm2, postprocess)
    profile_data = extract_combined_profile_tables(model_fwd, result_fwd, model_rev, result_rev)
    metrics = evaluate_case_metrics(case_id, stage, base_config, forward_curve, reverse_curve, reference_curve, profile_data)

    record = SimulationRecord(
        stage=stage,
        case_id=case_id,
        label=label,
        notes=notes,
        overrides=overrides,
        postprocess=postprocess,
        case_dir=case_dir,
        forward_curve=forward_curve,
        reverse_curve=reverse_curve,
        metrics=metrics,
    )
    save_case_outputs(record, model_fwd, result_fwd, result_rev, reference_curve, profile_data)
    return record


def choose_best(records: list[SimulationRecord]) -> SimulationRecord:
    return min(records, key=lambda item: item.metrics["score"])


def write_stage_summary(stage: str, title: str, records: list[SimulationRecord], intro: str) -> None:
    stage_dir = CAL_ROOT / stage
    lines = [f"# {title}", "", intro, "", "| Case | Leakage @ -1 V (A/cm^2) | n | Turn-on @ 0.5 A/cm^2 (V) | Curvature std(n) | Score |", "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for record in records:
        n_val = record.metrics["ideality_fit"]["ideality_factor"]
        turn_on = record.metrics["turn_on_0p5Acm2_v"]
        lines.append(
            "| [{case}]({case}/metrics.md) | {leak:.3e} | {nval} | {ton} | {curv} | {score:.3f} |".format(
                case=record.case_id,
                leak=record.metrics["reverse_leakage_at_minus_1V_a_cm2"],
                nval="N/A" if n_val is None else f"{n_val:.3f}",
                ton=f">{np.max(record.forward_curve.voltage_v):.2f}" if turn_on is None else f"{turn_on:.3f}",
                curv="N/A" if record.metrics["forward_curvature_std_n"] is None else f"{record.metrics['forward_curvature_std_n']:.3f}",
                score=record.metrics["score"],
            )
        )
    best = choose_best(records)
    lines.extend(["", f"Best stage candidate: `{best.case_id}`"])
    (stage_dir / "stage_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_overall_summary(records: list[SimulationRecord], best_record: SimulationRecord) -> None:
    lines = [
        "# Calibration Summary",
        "",
        f"Cases run: {len(records)}",
        f"Best candidate: `{best_record.case_id}`",
        "",
        "| Stage | Case | Leakage @ -1 V (A/cm^2) | n | Turn-on @ 0.5 A/cm^2 (V) | Score | Stop Condition |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for record in records:
        n_val = record.metrics["ideality_fit"]["ideality_factor"]
        turn_on = record.metrics["turn_on_0p5Acm2_v"]
        lines.append(
            "| {stage} | [{case}]({stage}/{case}/metrics.md) | {leak:.3e} | {nval} | {ton} | {score:.3f} | {stop} |".format(
                stage=record.stage,
                case=record.case_id,
                leak=record.metrics["reverse_leakage_at_minus_1V_a_cm2"],
                nval="N/A" if n_val is None else f"{n_val:.3f}",
                ton=f">{np.max(record.forward_curve.voltage_v):.2f}" if turn_on is None else f"{turn_on:.3f}",
                score=record.metrics["score"],
                stop=record.metrics["meets_stop_condition"],
            )
        )
    (CAL_ROOT / "calibration_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_layout()
    authoritative = load_authoritative_config(JSON_PATH)
    base_project = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    base_project["area_cm2"] = float(base_project["area"])
    reference_curve = load_reference_curve(authoritative["area_cm2"])
    all_records: list[SimulationRecord] = []

    stage_a_records: list[SimulationRecord] = []
    stage_a_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_a_contact",
            "A01_barrier_0p60",
            "Stage A Contact: bc_right = 0.60 V",
            {"bc_right": 0.60},
            "Contact inspection: `work_function_left/right` are placeholders in the current solver. The active contact-barrier proxy in this calibration is the electrostatic surface boundary `bc_right`, while contact resistance is treated as an external series parasitic for solver-7 compatibility.",
            {"external_rs_ohm": 11.7},
        )
    )
    stage_a_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_a_contact",
            "A02_barrier_0p30",
            "Stage A Contact: bc_right = 0.30 V",
            {"bc_right": 0.30},
            "Reduced the effective right-contact barrier proxy by halving `bc_right` while keeping the semiconductor stack fixed.",
            {"external_rs_ohm": 11.7},
        )
    )
    stage_a_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_a_contact",
            "A03_barrier_0p00",
            "Stage A Contact: bc_right = 0.00 V",
            {"bc_right": 0.00},
            "Removed the right-contact electrostatic offset to test whether the baseline turn-on was dominated by an overly restrictive contact boundary.",
            {"external_rs_ohm": 11.7},
        )
    )
    best_barrier = choose_best(stage_a_records)
    stage_a_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_a_contact",
            "A04_rs_50",
            "Stage A Contact: external Rs = 50 Ohm",
            dict(best_barrier.overrides),
            "Series-resistance probe applied externally because solver 7 does not solve contact resistivity self-consistently.",
            {"external_rs_ohm": 50.0},
        )
    )
    stage_a_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_a_contact",
            "A05_rs_200",
            "Stage A Contact: external Rs = 200 Ohm",
            dict(best_barrier.overrides),
            "Aggressive series-resistance probe to verify whether contact resistivity can explain the forward soft turn-on. This is expected to worsen, not improve, the knee if contacts are already resistive.",
            {"external_rs_ohm": 200.0},
        )
    )
    write_stage_summary(
        "stage_a_contact",
        "Stage A Contact Calibration",
        stage_a_records,
        "This stage inspects the active contact model used by the solver and sweeps the effective contact barrier proxy (`bc_right`) plus external series resistance.",
    )
    all_records.extend(stage_a_records)
    best_after_a = choose_best(stage_a_records)

    stage_b_records: list[SimulationRecord] = []
    stage_b_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_b_defects",
            "B01_tau_5e-10",
            "Stage B Defects: tau = 5e-10 s",
            {**best_after_a.overrides, "taun0": 5.0e-10, "taup0": 5.0e-10},
            "Shortened both SRH lifetimes by about two orders of magnitude from the InGaN alloy default to emulate a defect-richer high-indium absorber.",
            best_after_a.postprocess,
        )
    )
    stage_b_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_b_defects",
            "B02_tau_1e-10",
            "Stage B Defects: tau = 1e-10 s",
            {**best_after_a.overrides, "taun0": 1.0e-10, "taup0": 1.0e-10},
            "Pushed the SRH lifetime deeper into the defect-limited regime while leaving contact settings unchanged.",
            best_after_a.postprocess,
        )
    )
    best_tau = choose_best(stage_b_records)
    stage_b_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_b_defects",
            "B03_tatfield_1e6",
            "Stage B Defects: tat_field = 1e6 V/m",
            {**best_tau.overrides, "tat_field": 1.0e6},
            "Strengthened Hurkx-like trap-assisted tunneling by lowering the characteristic field threshold by a factor of five.",
            best_tau.postprocess,
        )
    )
    best_tat = choose_best(stage_b_records)
    stage_b_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_b_defects",
            "B04_trap_density_50",
            "Stage B Defects: trap_density_scale = 50",
            {**best_tat.overrides, "trap_density_scale": 50.0},
            "Introduced a phenomenological deep-defect density scaling that boosts SRH/TAT center density without changing the layer stack.",
            best_tat.postprocess,
        )
    )
    best_density = choose_best(stage_b_records)
    stage_b_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_b_defects",
            "B05_trap_energy_0p08",
            "Stage B Defects: trap_energy_offset = +0.08 eV",
            {**best_density.overrides, "trap_energy_offset_ev": 0.08},
            "Shifted the effective trap level slightly toward the conduction band to test whether asymmetry in deep centers improves the leakage/ideality tradeoff.",
            best_density.postprocess,
        )
    )
    write_stage_summary(
        "stage_b_defects",
        "Stage B Defect-Assisted Leakage Calibration",
        stage_b_records,
        "This stage strengthens SRH generation, TAT enhancement, and deep-defect occupancy one knob at a time while preserving the contact settings selected from Stage A.",
    )
    all_records.extend(stage_b_records)
    best_after_b = choose_best(stage_b_records)

    stage_c_records: list[SimulationRecord] = []
    stage_c_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_c_mobility",
            "C01_mu_e0p7_h0p3",
            "Stage C Mobility: mun x0.7, mup x0.3",
            {**best_after_b.overrides, "mun0_scale": 0.7, "mup0_scale": 0.3},
            "Moderately reduced electron mobility and aggressively reduced hole mobility to mimic defect- and alloy-scattering in high-indium InGaN.",
            best_after_b.postprocess,
        )
    )
    stage_c_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_c_mobility",
            "C02_mu_e0p5_h0p1",
            "Stage C Mobility: mun x0.5, mup x0.1",
            {**best_after_b.overrides, "mun0_scale": 0.5, "mup0_scale": 0.1},
            "More aggressive mobility degradation to test whether forward curvature can be made visibly non-ideal without destroying the turn-on target.",
            best_after_b.postprocess,
        )
    )
    write_stage_summary(
        "stage_c_mobility",
        "Stage C Mobility Calibration",
        stage_c_records,
        "This stage reduces low-field carrier mobility to emulate composition- and defect-limited transport while keeping the calibrated contact and defect settings fixed.",
    )
    all_records.extend(stage_c_records)
    best_after_c = choose_best(stage_c_records)

    stage_d_records: list[SimulationRecord] = []
    stage_d_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_d_polarization",
            "D01_polscale_0p60",
            "Stage D Polarization: overall scale 0.60",
            {**best_after_c.overrides, "polarization_scale": 0.60},
            "Applied partial overall polarization relaxation to test whether overestimated built-in fields are suppressing forward conduction.",
            best_after_c.postprocess,
        )
    )
    stage_d_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_d_polarization",
            "D02_piezo_0p40",
            "Stage D Polarization: piezo scale 0.40",
            {**best_after_c.overrides, "piezoelectric_scale": 0.40, "spontaneous_scale": 1.0},
            "Reduced only the piezoelectric contribution while leaving spontaneous polarization intact to mimic partial strain relaxation and defect screening.",
            best_after_c.postprocess,
        )
    )
    write_stage_summary(
        "stage_d_polarization",
        "Stage D Polarization Calibration",
        stage_d_records,
        "This stage checks whether the high-indium forward bottleneck is driven by overestimated polarization-induced band bending.",
    )
    all_records.extend(stage_d_records)
    best_after_d = choose_best(stage_d_records)

    candidate_records: list[SimulationRecord] = []
    candidate_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_e_candidates",
            "E01_best_chain",
            "Stage E Candidate: best chain",
            dict(best_after_d.overrides),
            "Carries forward the best-scoring chain from Stages A-D without introducing any new physics knob.",
            best_after_d.postprocess,
        )
    )
    candidate_records.append(
        run_case(
            base_project,
            reference_curve,
            "stage_e_candidates",
            "E02_aggressive_defect_transport",
            "Stage E Candidate: aggressive defect + moderated transport",
            {
                **best_after_a.overrides,
                "taun0": 1.0e-10,
                "taup0": 1.0e-10,
                "tat_field": 1.0e6,
                "trap_density_scale": 50.0,
                "trap_energy_offset_ev": 0.08,
                "mun0_scale": 0.7,
                "mup0_scale": 0.3,
                "piezoelectric_scale": 0.40,
            },
            "Explicitly combines the best defect-assisted leakage knobs with moderated mobility degradation and reduced piezoelectric field to search for a defect-dominated practical operating point.",
            best_after_a.postprocess,
        )
    )
    write_stage_summary(
        "stage_e_candidates",
        "Stage E Candidate Search",
        candidate_records,
        "This stage stores candidate calibrated structures that best balance leakage, ideality, turn-on, and visible non-ideal curvature.",
    )
    all_records.extend(candidate_records)

    best_record = choose_best(all_records)
    build_overall_summary(all_records, best_record)
    print(f"Calibration complete. Best case: {best_record.case_id} ({best_record.stage})")
    print(f"Summary: {CAL_ROOT / 'calibration_summary.md'}")


if __name__ == "__main__":
    main()
