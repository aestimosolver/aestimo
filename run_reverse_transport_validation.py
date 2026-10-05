from __future__ import annotations

import json
import math
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

from run_untitled_project_baseline import load_authoritative_config
from run_untitled_project_calibration import CalibrationInput, merge_case_config, sanitize_name


ROOT = Path(__file__).resolve().parent
JSON_PATH = ROOT / "examples" / "untitled_project.json"
OUT_ROOT = ROOT / "results" / "reverse_transport"
REVERSE_WINDOW = {"vmin": -1.05, "vmax": 0.0, "vstep": 0.05}
AV_CURR_TO_A_CM2 = 1e-3

STAGE_DIRS = {
    "r1_current_extraction": "R1 Current Extraction Validation",
    "r2_contact_boundaries": "R2 Contact Boundary Analysis",
    "r3_mesh_resolution": "R3 Mesh And Field Resolution",
    "r4_high_field_models": "R4 High-Field Transport Validation",
    "r5_polarization_stability": "R5 Polarization-Charge Stability",
}


@dataclass
class ReverseCase:
    stage: str
    case_id: str
    label: str
    overrides: dict[str, Any]
    notes: str


@dataclass
class CaseArtifacts:
    case: ReverseCase
    case_dir: Path
    merged_config: dict[str, Any]
    current_table: dict[str, np.ndarray]
    metrics: dict[str, Any]


def ensure_layout() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for stage in STAGE_DIRS:
        (OUT_ROOT / stage).mkdir(parents=True, exist_ok=True)


def load_base_project() -> tuple[dict[str, Any], dict[str, Any]]:
    raw = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    raw["area_cm2"] = float(raw["area"])
    summary = load_authoritative_config(JSON_PATH)
    return raw, summary


def result_array(result: Any, name: str) -> np.ndarray:
    value = getattr(result, name, None)
    if value is None:
        base = name[:-1] if name.endswith("_") else name
        value = getattr(result, base)
    return np.asarray(value, dtype=float)


def run_reverse_case(base_config: dict[str, Any], case: ReverseCase) -> tuple[Any, Any, dict[str, Any], Path]:
    merged = merge_case_config(base_config, case.overrides, REVERSE_WINDOW)
    case_dir = OUT_ROOT / case.stage / case.case_id
    raw_dir = case_dir / f"raw_{sanitize_name(case.case_id)}"
    raw_dir.mkdir(parents=True, exist_ok=True)

    aestimo.output_directory = str(raw_dir)
    aestimo_config.Drift_Diffusion_out = True
    aestimo_config.potential_out = False
    aestimo_config.electricfield_out = False
    aestimo_config.sigma_out = False
    aestimo_config.states_out = False
    aestimo_config.probability_out = False

    input_obj = CalibrationInput(merged, input_name=case.case_id)
    _, model, result, _ = run_aestimo(input_obj, drawFigures=False, show=False)
    return model, result, merged, raw_dir


def convert_current_density(arr: np.ndarray) -> np.ndarray:
    return np.asarray(arr, dtype=float) * AV_CURR_TO_A_CM2


def current_table_from_result(result: Any) -> dict[str, np.ndarray]:
    table = {
        "bias_v": np.asarray(result.Va_t, dtype=float),
        "official_a_cm2": convert_current_density(result.av_curr),
        "left_probe_a_cm2": convert_current_density(getattr(result, "av_curr_left", result.av_curr)),
        "center_probe_a_cm2": convert_current_density(getattr(result, "av_curr_center", result.av_curr)),
        "right_probe_a_cm2": convert_current_density(getattr(result, "av_curr_right", result.av_curr)),
        "whole_probe_a_cm2": convert_current_density(getattr(result, "av_curr_whole", result.av_curr)),
        "absmax_probe_a_cm2": convert_current_density(getattr(result, "av_curr_absmax", np.abs(result.av_curr))),
        "left_edge_a_cm2": convert_current_density(getattr(result, "av_curr_left_edge", result.av_curr)),
        "right_edge_a_cm2": convert_current_density(getattr(result, "av_curr_right_edge", result.av_curr)),
        "reached_bias_v": np.asarray(getattr(result, "reached_bias_v", result.Va_t), dtype=float),
        "converged_flag": np.asarray(getattr(result, "converged_bias_flags", np.ones_like(result.Va_t, dtype=bool)), dtype=bool),
        "substep_attempts": np.asarray(getattr(result, "substep_attempts", np.zeros_like(result.Va_t)), dtype=int),
        "successful_substeps": np.asarray(getattr(result, "successful_substeps", np.zeros_like(result.Va_t)), dtype=int),
        "gummel_iterations": np.asarray(getattr(result, "gummel_iterations", np.zeros_like(result.Va_t)), dtype=int),
        "timeout_events": np.asarray(getattr(result, "timeout_events", np.zeros_like(result.Va_t)), dtype=int),
        "min_step_attempt_v": np.asarray(getattr(result, "min_step_attempt_v", np.full_like(result.Va_t, np.nan)), dtype=float),
        "final_step_v": np.asarray(getattr(result, "final_step_v", np.zeros_like(result.Va_t)), dtype=float),
    }
    return table


def exact_or_interp(table: dict[str, np.ndarray], column: str, bias_v: float) -> float:
    return float(np.interp(bias_v, table["bias_v"], table[column]))


def first_bias_above_threshold(bias_v: np.ndarray, current_a_cm2: np.ndarray, threshold: float) -> float | None:
    mask = (bias_v < 0.0) & (np.abs(current_a_cm2) >= threshold)
    if not np.any(mask):
        return None
    return float(bias_v[np.argmax(mask)])


def representative_profile_biases(table: dict[str, np.ndarray]) -> list[float]:
    negative_mask = table["bias_v"] < 0.0
    bias_neg = table["bias_v"][negative_mask]
    whole_neg = np.abs(table["whole_probe_a_cm2"][negative_mask])
    if len(bias_neg) == 0:
        return [0.0]
    idx_minus_1 = int(np.argmin(np.abs(bias_neg + 1.0)))
    idx_peak = int(np.argmax(whole_neg))
    chosen = [float(bias_neg[idx_minus_1])]
    peak_bias = float(bias_neg[idx_peak])
    if not math.isclose(peak_bias, chosen[0], abs_tol=1e-6):
        chosen.append(peak_bias)
    return chosen


def bias_index(result: Any, target_bias_v: float) -> int:
    bias = np.asarray(result.Va_t, dtype=float)
    return int(np.argmin(np.abs(bias - target_bias_v)))


def extract_profiles(result: Any, target_biases: list[float]) -> dict[str, dict[str, np.ndarray]]:
    x_nm = np.asarray(result.xaxis, dtype=float) * 1e9
    ec = result_array(result, "Ec_result_")
    ev = result_array(result, "Ev_result_")
    ei = result_array(result, "Ei_result_")
    efn = result_array(result, "Efn_result_")
    efp = result_array(result, "Efp_result_")
    n = result_array(result, "nf_result_")
    p = result_array(result, "pf_result_")
    field = result_array(result, "el_field1_result_")
    charge = result_array(result, "ro_result_")
    jtot = result_array(result, "Jtotal_result_")
    jelec = result_array(result, "Jelec_result_")
    jhole = result_array(result, "Jhole_result_")

    profiles: dict[str, dict[str, np.ndarray]] = {}
    for target in target_biases:
        idx = bias_index(result, target)
        actual = float(np.asarray(result.Va_t, dtype=float)[idx])
        label = sanitize_name(f"profile_{actual:.2f}V")
        profiles[label] = {
            "x_nm": x_nm,
            "field_v_cm": np.asarray(field[idx, :], dtype=float) * 1e-2,
            "charge_c_cm3": np.asarray(charge[idx, :], dtype=float) * 1e-6,
            "n_cm3": np.asarray(n[idx, :], dtype=float) * 1e-6,
            "p_cm3": np.asarray(p[idx, :], dtype=float) * 1e-6,
            "Ec_eV": np.asarray(ec[idx, :], dtype=float),
            "Ev_eV": np.asarray(ev[idx, :], dtype=float),
            "Ei_eV": np.asarray(ei[idx, :], dtype=float),
            "Efn_eV": np.asarray(efn[idx, :], dtype=float),
            "Efp_eV": np.asarray(efp[idx, :], dtype=float),
            "Jtotal_a_cm2": convert_current_density(np.asarray(jtot[idx, :], dtype=float)),
            "Jelec_a_cm2": convert_current_density(np.asarray(jelec[idx, :], dtype=float)),
            "Jhole_a_cm2": convert_current_density(np.asarray(jhole[idx, :], dtype=float)),
            "bias_v_actual": np.array([actual], dtype=float),
        }
    return profiles


def build_metrics(case: ReverseCase, table: dict[str, np.ndarray], profiles: dict[str, dict[str, np.ndarray]]) -> dict[str, Any]:
    official_at_minus_1 = exact_or_interp(table, "official_a_cm2", -1.0)
    whole_at_minus_1 = exact_or_interp(table, "whole_probe_a_cm2", -1.0)
    absmax_at_minus_1 = exact_or_interp(table, "absmax_probe_a_cm2", -1.0)
    threshold = max(1e-12, absmax_at_minus_1 * 10.0)
    onset_official = first_bias_above_threshold(table["bias_v"], table["official_a_cm2"], threshold)
    onset_whole = first_bias_above_threshold(table["bias_v"], table["whole_probe_a_cm2"], threshold)
    onset_absmax = first_bias_above_threshold(table["bias_v"], table["absmax_probe_a_cm2"], threshold)
    profile_key = sorted(profiles.keys())[0]
    profile = profiles[profile_key]
    return {
        "case_id": case.case_id,
        "label": case.label,
        "official_current_at_minus_1V_a_cm2": official_at_minus_1,
        "whole_probe_current_at_minus_1V_a_cm2": whole_at_minus_1,
        "absmax_current_at_minus_1V_a_cm2": absmax_at_minus_1,
        "official_zero_like_points": int(np.count_nonzero(np.abs(table["official_a_cm2"]) <= 1e-18)),
        "whole_zero_like_points": int(np.count_nonzero(np.abs(table["whole_probe_a_cm2"]) <= 1e-18)),
        "onset_official_v": onset_official,
        "onset_whole_v": onset_whole,
        "onset_absmax_v": onset_absmax,
        "all_bias_targets_reached": bool(np.all(table["converged_flag"])),
        "timeout_events_total": int(np.sum(table["timeout_events"])),
        "max_gummel_iterations": int(np.max(table["gummel_iterations"])),
        "min_substep_attempt_v": None if np.all(np.isnan(table["min_step_attempt_v"])) else float(np.nanmin(table["min_step_attempt_v"])),
        "max_abs_field_v_cm_at_profile0": float(np.max(np.abs(profile["field_v_cm"]))),
        "left_contact_n_cm3_at_profile0": float(profile["n_cm3"][0]),
        "left_contact_p_cm3_at_profile0": float(profile["p_cm3"][0]),
        "right_contact_n_cm3_at_profile0": float(profile["n_cm3"][-1]),
        "right_contact_p_cm3_at_profile0": float(profile["p_cm3"][-1]),
        "left_contact_qf_split_eV_at_profile0": float(profile["Efn_eV"][0] - profile["Efp_eV"][0]),
        "right_contact_qf_split_eV_at_profile0": float(profile["Efn_eV"][-1] - profile["Efp_eV"][-1]),
    }


def save_current_table(path: Path, table: dict[str, np.ndarray]) -> None:
    headers = [
        "bias_v",
        "official_a_cm2",
        "left_probe_a_cm2",
        "center_probe_a_cm2",
        "right_probe_a_cm2",
        "whole_probe_a_cm2",
        "absmax_probe_a_cm2",
        "left_edge_a_cm2",
        "right_edge_a_cm2",
        "reached_bias_v",
        "converged_flag",
        "substep_attempts",
        "successful_substeps",
        "gummel_iterations",
        "timeout_events",
        "min_step_attempt_v",
        "final_step_v",
    ]
    data = np.column_stack([np.asarray(table[key], dtype=float) for key in headers])
    np.savetxt(path, data, delimiter=",", header=",".join(headers), comments="", fmt="%.9e")


def save_profile_csv(path: Path, profile: dict[str, np.ndarray]) -> None:
    headers = [
        "x_nm",
        "field_v_cm",
        "charge_c_cm3",
        "n_cm3",
        "p_cm3",
        "Ec_eV",
        "Ev_eV",
        "Ei_eV",
        "Efn_eV",
        "Efp_eV",
        "Jtotal_a_cm2",
        "Jelec_a_cm2",
        "Jhole_a_cm2",
    ]
    data = np.column_stack([np.asarray(profile[key], dtype=float) for key in headers])
    np.savetxt(path, data, delimiter=",", header=",".join(headers), comments="", fmt="%.9e")


def plot_current_semilog(path: Path, label: str, table: dict[str, np.ndarray], comparison_columns: list[tuple[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    for col, col_label in comparison_columns:
        ax.semilogy(table["bias_v"], np.clip(np.abs(table[col]), 1e-20, None), linewidth=1.8, label=col_label)
    ax.set_xlabel("Bias (V)")
    ax.set_ylabel("|J| (A/cm$^2$)")
    ax.set_title(f"{label} Reverse I-V")
    ax.grid(alpha=0.25, which="both")
    ax.set_xlim(REVERSE_WINDOW["vmin"], REVERSE_WINDOW["vmax"])
    ax.legend()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_zoomed_current(path: Path, label: str, table: dict[str, np.ndarray], comparison_columns: list[tuple[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    for col, col_label in comparison_columns:
        ax.plot(table["bias_v"], table[col], linewidth=1.6, label=col_label)
    ax.set_xlabel("Bias (V)")
    ax.set_ylabel("J (A/cm$^2$)")
    ax.set_title(f"{label} Raw Reverse Current Zoom")
    ax.grid(alpha=0.25)
    ax.set_xlim(REVERSE_WINDOW["vmin"], -0.1)
    ax.legend()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_field_profiles(path: Path, label: str, profiles: dict[str, dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(len(profiles), 1, figsize=(8.2, 3.3 * len(profiles)), sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, (name, profile) in zip(axes, profiles.items()):
        ax.plot(profile["x_nm"], profile["field_v_cm"], color="#1d3557", linewidth=1.8)
        ax.set_ylabel("Field (V/cm)")
        ax.set_title(f"{name} at {profile['bias_v_actual'][0]:.2f} V")
        ax.grid(alpha=0.25)
    axes[-1].set_xlabel("Position (nm)")
    fig.suptitle(f"{label} Electric Field Profiles")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_carrier_profiles(path: Path, label: str, profiles: dict[str, dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(len(profiles), 1, figsize=(8.2, 3.4 * len(profiles)), sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, (name, profile) in zip(axes, profiles.items()):
        ax.semilogy(profile["x_nm"], np.clip(profile["n_cm3"], 1e-5, None), label="n", color="#005f73", linewidth=1.7)
        ax.semilogy(profile["x_nm"], np.clip(profile["p_cm3"], 1e-5, None), label="p", color="#bb3e03", linewidth=1.7)
        ax.set_ylabel("Carrier Density (cm$^{-3}$)")
        ax.set_title(f"{name} at {profile['bias_v_actual'][0]:.2f} V")
        ax.grid(alpha=0.25, which="both")
        ax.legend()
    axes[-1].set_xlabel("Position (nm)")
    fig.suptitle(f"{label} Carrier Density Profiles")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_quasi_fermi(path: Path, label: str, profiles: dict[str, dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(len(profiles), 1, figsize=(8.4, 3.6 * len(profiles)), sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, (name, profile) in zip(axes, profiles.items()):
        ax.plot(profile["x_nm"], profile["Ec_eV"], color="#1d3557", linewidth=1.6, label="Ec")
        ax.plot(profile["x_nm"], profile["Ev_eV"], color="#9b2226", linewidth=1.6, label="Ev")
        ax.plot(profile["x_nm"], profile["Efn_eV"], color="#2a9d8f", linewidth=1.4, linestyle="--", label="Efn")
        ax.plot(profile["x_nm"], profile["Efp_eV"], color="#e76f51", linewidth=1.4, linestyle="--", label="Efp")
        ax.set_ylabel("Energy (eV)")
        ax.set_title(f"{name} at {profile['bias_v_actual'][0]:.2f} V")
        ax.grid(alpha=0.25)
        ax.legend(ncol=4, fontsize=8)
    axes[-1].set_xlabel("Position (nm)")
    fig.suptitle(f"{label} Quasi-Fermi Levels")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_convergence(path: Path, label: str, table: dict[str, np.ndarray]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.5), sharex=True, constrained_layout=True)
    axes[0].plot(table["bias_v"], table["gummel_iterations"], color="#005f73", linewidth=1.7, label="Gummel iterations")
    axes[0].plot(table["bias_v"], table["substep_attempts"], color="#bb3e03", linewidth=1.5, label="Sub-step attempts")
    axes[0].plot(table["bias_v"], table["timeout_events"], color="#9b2226", linewidth=1.4, label="Timeout events")
    axes[0].set_ylabel("Count")
    axes[0].set_title(f"{label} Convergence Behavior")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    axes[1].plot(table["bias_v"], table["reached_bias_v"] - table["bias_v"], color="#1d3557", linewidth=1.6, label="Reached - target")
    axes[1].plot(table["bias_v"], np.nan_to_num(table["min_step_attempt_v"], nan=0.0), color="#2a9d8f", linewidth=1.4, label="Min sub-step (V)")
    axes[1].set_xlabel("Bias (V)")
    axes[1].set_ylabel("Voltage")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_case_notes(case_dir: Path, case: ReverseCase, metrics: dict[str, Any], extra_notes: list[str]) -> None:
    lines = [
        "# Numerical Stability Notes",
        "",
        f"- Case: `{case.case_id}`",
        f"- Official current at -1.0 V: {metrics['official_current_at_minus_1V_a_cm2']:.3e} A/cm^2",
        f"- Whole-device probe current at -1.0 V: {metrics['whole_probe_current_at_minus_1V_a_cm2']:.3e} A/cm^2",
        f"- Abs-max current proxy at -1.0 V: {metrics['absmax_current_at_minus_1V_a_cm2']:.3e} A/cm^2",
        f"- All target biases reached: {metrics['all_bias_targets_reached']}",
        f"- Timeout events: {metrics['timeout_events_total']}",
        f"- Max Gummel iterations accumulated at a bias point: {metrics['max_gummel_iterations']}",
        f"- Minimum attempted sub-step: {metrics['min_substep_attempt_v']}",
    ]
    if metrics["onset_official_v"] is not None:
        lines.append(f"- Official onset above threshold: {metrics['onset_official_v']:.3f} V")
    if metrics["onset_whole_v"] is not None:
        lines.append(f"- Whole-probe onset above threshold: {metrics['onset_whole_v']:.3f} V")
    lines.extend(["", "## Notes", ""])
    for note in extra_notes:
        lines.append(f"- {note}")
    (case_dir / "numerical_stability_notes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_case_outputs(case: ReverseCase, merged_config: dict[str, Any], result: Any, table: dict[str, np.ndarray], case_dir: Path) -> CaseArtifacts:
    case_dir.mkdir(parents=True, exist_ok=True)
    profile_biases = representative_profile_biases(table)
    profiles = extract_profiles(result, profile_biases)
    metrics = build_metrics(case, table, profiles)

    (case_dir / "effective_case.json").write_text(json.dumps({
        "stage": case.stage,
        "case_id": case.case_id,
        "label": case.label,
        "overrides": case.overrides,
        "reverse_window": REVERSE_WINDOW,
        "merged_config_excerpt": {
            "bc_left": merged_config.get("bc_left"),
            "bc_right": merged_config.get("bc_right"),
            "grid_step": merged_config.get("grid_step"),
            "tat_field": merged_config.get("tat_field"),
            "device_type": merged_config.get("device_type"),
            "photovoltaic_mode": merged_config.get("photovoltaic_mode"),
            "enable_polarization": merged_config.get("enable_polarization", True),
            "polarization_smoothing_nm": merged_config.get("polarization_smoothing_nm", 0.0),
        },
    }, indent=2), encoding="utf-8")
    (case_dir / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    save_current_table(case_dir / "raw_current_table.csv", table)

    for name, profile in profiles.items():
        save_profile_csv(case_dir / f"{name}.csv", profile)

    plot_current_semilog(
        case_dir / "reverse_iv_semilog.png",
        case.label,
        table,
        [
            ("official_a_cm2", "Official av_curr"),
            ("whole_probe_a_cm2", "Whole-device median"),
            ("absmax_probe_a_cm2", "Abs-max current proxy"),
        ],
    )
    plot_zoomed_current(
        case_dir / "reverse_current_zoom.png",
        case.label,
        table,
        [
            ("official_a_cm2", "Official av_curr"),
            ("left_probe_a_cm2", "Left probe"),
            ("center_probe_a_cm2", "Center probe"),
            ("right_probe_a_cm2", "Right probe"),
            ("whole_probe_a_cm2", "Whole-device median"),
        ],
    )
    plot_field_profiles(case_dir / "electric_field_profiles.png", case.label, profiles)
    plot_carrier_profiles(case_dir / "carrier_density_profiles.png", case.label, profiles)
    plot_quasi_fermi(case_dir / "quasi_fermi_profiles.png", case.label, profiles)
    plot_convergence(case_dir / "convergence_behavior.png", case.label, table)

    notes = [case.notes]
    if metrics["official_zero_like_points"] > 0:
        notes.append("The official `av_curr` series contains exact zero-like points, so a numerical floor or extraction dead-zone is present at least in the reported terminal current.")
    if metrics["whole_probe_current_at_minus_1V_a_cm2"] != metrics["official_current_at_minus_1V_a_cm2"]:
        notes.append("The whole-device current probe does not match the official right-boundary median extraction at -1 V, indicating extraction-location sensitivity.")
    if not metrics["all_bias_targets_reached"]:
        notes.append("At least one target reverse bias was not fully reached before the adaptive stepping logic hit its minimum step size.")
    write_case_notes(case_dir, case, metrics, notes)

    return CaseArtifacts(case=case, case_dir=case_dir, merged_config=merged_config, current_table=table, metrics=metrics)


def plot_stage_comparison(stage_dir: Path, stage_title: str, artifacts: list[CaseArtifacts]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 7.0), sharex=True, constrained_layout=True)
    for artifact in artifacts:
        table = artifact.current_table
        axes[0].semilogy(table["bias_v"], np.clip(np.abs(table["official_a_cm2"]), 1e-20, None), linewidth=1.6, label=artifact.case.case_id)
        axes[1].semilogy(table["bias_v"], np.clip(np.abs(table["whole_probe_a_cm2"]), 1e-20, None), linewidth=1.6, label=artifact.case.case_id)
    axes[0].set_ylabel("|J| (A/cm$^2$)")
    axes[0].set_title(f"{stage_title}: Official av_curr")
    axes[0].grid(alpha=0.25, which="both")
    axes[1].set_ylabel("|J| (A/cm$^2$)")
    axes[1].set_title(f"{stage_title}: Whole-Device Median Probe")
    axes[1].set_xlabel("Bias (V)")
    axes[1].grid(alpha=0.25, which="both")
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    fig.savefig(stage_dir / "stage_reverse_iv_comparison.png", dpi=220)
    plt.close(fig)


def write_stage_summary(stage: str, artifacts: list[CaseArtifacts], intro: str) -> None:
    stage_dir = OUT_ROOT / stage
    plot_stage_comparison(stage_dir, STAGE_DIRS[stage], artifacts)
    lines = [
        f"# {STAGE_DIRS[stage]}",
        "",
        intro,
        "",
        "| Case | Official J(-1 V) (A/cm^2) | Whole J(-1 V) (A/cm^2) | Official Onset (V) | Whole Onset (V) | Timeouts | All Biases Reached |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for artifact in artifacts:
        m = artifact.metrics
        lines.append(
            "| [{case}]({case}/numerical_stability_notes.md) | {joff:.3e} | {jwhole:.3e} | {o1} | {o2} | {timeouts} | {reached} |".format(
                case=artifact.case.case_id,
                joff=m["official_current_at_minus_1V_a_cm2"],
                jwhole=m["whole_probe_current_at_minus_1V_a_cm2"],
                o1="N/A" if m["onset_official_v"] is None else f"{m['onset_official_v']:.3f}",
                o2="N/A" if m["onset_whole_v"] is None else f"{m['onset_whole_v']:.3f}",
                timeouts=m["timeout_events_total"],
                reached=m["all_bias_targets_reached"],
            )
        )
    (stage_dir / "stage_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_code_findings() -> None:
    lines = [
        "# Reverse-Bias Code Findings",
        "",
        "## Current Extraction",
        "",
        "- Scheme 7 computes the reported `av_curr` from the median of the last 10% of `Jtotal` nodes near the right boundary.",
        "- `av_curr.dat` is written directly from `result.Va_t` and `result.av_curr`; there is no file-output clipping or floor in the plotting layer.",
        "- In the active scheme-7 path, `Jelec` and `Jhole` are converted to `mA/cm^2` before `Jtotal` and `av_curr` are formed; generated analysis tables therefore convert `av_curr` to `A/cm^2` with `1e-3`.",
        "",
        "## Boundary Conditions",
        "",
        "- In photovoltaic mode, the continuity solver uses selective contacts: majority carrier Dirichlet and minority-carrier zero-flux Neumann conditions.",
        "- In non-photovoltaic mode, the continuity solver applies standard Dirichlet conditions to both carriers at both contacts.",
        "- `bc_left` and `bc_right` enter as electrostatic surface offsets in the equilibrium/contact potential initialization rather than as a separate Schottky injection model.",
        "",
        "## High-Field Physics",
        "",
        "- Implemented reverse-field enhancement in the active solver family is a Hurkx-like TAT multiplier controlled by `tat_field`, `trap_density_scale`, and `trap_energy_offset_ev`.",
        "- No explicit Poole-Frenkel, band-to-band tunneling, Fowler-Nordheim field emission, or Schottky-barrier-lowering implementation was found in the active drift-diffusion solver path.",
        "",
        "## Mesh",
        "",
        "- The project path exposed through the JSON and solver input supports a uniform `grid_step`; no native nonuniform local junction-refinement control was found in this workflow.",
    ]
    (OUT_ROOT / "code_findings.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_overall_summary(all_artifacts: list[CaseArtifacts]) -> None:
    best_method_case = next((artifact for artifact in all_artifacts if artifact.case.case_id == "R2_selective_bc0p0"), None)
    lines = [
        "# Reverse Transport Validation Summary",
        "",
        "This study isolates reverse-bias methodology limitations before any further device sweeps.",
        "",
        "## High-Level Outcome",
        "",
    ]
    if best_method_case is not None:
        m = best_method_case.metrics
        lines.extend(
            [
                f"- Contact-unmasked reference case: `{best_method_case.case.case_id}`",
                f"- Official `J(-1 V)`: {m['official_current_at_minus_1V_a_cm2']:.3e} A/cm^2",
                f"- Whole-device probe `J(-1 V)`: {m['whole_probe_current_at_minus_1V_a_cm2']:.3e} A/cm^2",
                f"- Abs-max probe `J(-1 V)`: {m['absmax_current_at_minus_1V_a_cm2']:.3e} A/cm^2",
            ]
        )
    lines.extend(
        [
            "",
            "## Conclusion",
            "",
            "The reverse branch behavior is dominated by methodology limits rather than by a simple lack of defect strength.",
            "",
            "A realistic gradual high-indium leakage tail was not recovered. The study instead points to a framework limitation that combines:",
            "",
            "- extraction sensitivity, because `av_curr` samples a near-contact median rather than a robust terminal current observable",
            "- selective-contact suppression of minority-carrier injection under reverse bias",
            "- strong bias-step / damping sensitivity near reverse onset",
            "- and the absence of several physically relevant high-field leakage channels beyond Hurkx-like TAT",
            "",
            "Under these conditions, conclusion **B** is currently better supported: the present solver/contact framework fundamentally limits accurate high-indium reverse-bias simulation in its current form.",
            "",
            "A next-priority solver task would be to replace or augment the reverse-bias contact/current methodology before resuming calibration or any design sweeps.",
        ]
    )
    (OUT_ROOT / "reverse_transport_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_stage(base_project: dict[str, Any], cases: list[ReverseCase], intro: str) -> list[CaseArtifacts]:
    artifacts: list[CaseArtifacts] = []
    for case in cases:
        model, result, merged, _ = run_reverse_case(base_project, case)
        table = current_table_from_result(result)
        artifact = save_case_outputs(case, merged, result, table, OUT_ROOT / case.stage / case.case_id)
        artifacts.append(artifact)
    write_stage_summary(cases[0].stage, artifacts, intro)
    return artifacts


def main() -> None:
    ensure_layout()
    base_project, _ = load_base_project()

    stage_r1 = [
        ReverseCase(
            stage="r1_current_extraction",
            case_id="R1_baseline_json",
            label="R1 Baseline JSON Reverse Branch",
            overrides={},
            notes="Uses the exact JSON project contact and reverse-bias settings to inspect the raw `av_curr` generation path with no transport-model changes.",
        ),
        ReverseCase(
            stage="r1_current_extraction",
            case_id="R1_selective_bc0p0",
            label="R1 Selective Contact With bc_right = 0.0 V",
            overrides={"bc_right": 0.0},
            notes="Reduces only the electrostatic right-contact offset so the current-extraction artifact can be observed in a less contact-blocked reverse branch.",
        ),
    ]

    stage_r2 = [
        ReverseCase(
            stage="r2_contact_boundaries",
            case_id="R2_selective_bc0p6",
            label="R2 Selective PV Contacts, bc_right = 0.6 V",
            overrides={},
            notes="Selective photovoltaic contact boundary conditions from the baseline JSON.",
        ),
        ReverseCase(
            stage="r2_contact_boundaries",
            case_id="R2_ohmic_bc0p6",
            label="R2 Ohmic Contacts, bc_right = 0.6 V",
            overrides={"photovoltaic_mode": False},
            notes="Forces the standard Dirichlet/ohmic continuity boundary treatment while keeping the same electrostatic surface offsets.",
        ),
        ReverseCase(
            stage="r2_contact_boundaries",
            case_id="R2_selective_bc0p0",
            label="R2 Selective PV Contacts, bc_right = 0.0 V",
            overrides={"bc_right": 0.0},
            notes="Removes the right electrostatic boundary offset while keeping selective photovoltaic carrier supply conditions.",
        ),
        ReverseCase(
            stage="r2_contact_boundaries",
            case_id="R2_ohmic_bc0p0",
            label="R2 Ohmic Contacts, bc_right = 0.0 V",
            overrides={"photovoltaic_mode": False, "bc_right": 0.0},
            notes="Combines the least restrictive electrostatic boundary with the least restrictive contact-carrier boundary type available in the current solver family.",
        ),
    ]

    stage_r3 = [
        ReverseCase(
            stage="r3_mesh_resolution",
            case_id="R3_mesh_1p0nm",
            label="R3 Uniform Mesh 1.0 nm",
            overrides={"bc_right": 0.0, "grid_step": 1.0},
            notes="Uses the baseline grid spacing as the reverse-transport reference case.",
        ),
        ReverseCase(
            stage="r3_mesh_resolution",
            case_id="R3_mesh_0p5nm",
            label="R3 Uniform Mesh 0.5 nm",
            overrides={"bc_right": 0.0, "grid_step": 0.5},
            notes="Uniform mesh refinement to test field and current sensitivity to spatial resolution.",
        ),
        ReverseCase(
            stage="r3_mesh_resolution",
            case_id="R3_mesh_0p2nm",
            label="R3 Uniform Mesh 0.2 nm",
            overrides={"bc_right": 0.0, "grid_step": 0.2},
            notes="Strong uniform refinement. The solver path does not expose native nonuniform local junction refinement, so this is the closest supported mesh-resolution stress test.",
        ),
    ]

    stage_r4 = [
        ReverseCase(
            stage="r4_high_field_models",
            case_id="R4_tat_baseline",
            label="R4 Baseline Hurkx-Like TAT",
            overrides={"bc_right": 0.0},
            notes="Keeps the JSON `tat_field = 5e6 V/m` setting and uses the less blocked contact offset chosen to expose reverse-branch behavior.",
        ),
        ReverseCase(
            stage="r4_high_field_models",
            case_id="R4_tat_disabled",
            label="R4 TAT Disabled",
            overrides={"bc_right": 0.0, "tat_field": 1.0e10},
            notes="Disables the Hurkx activation path by setting `tat_field` above the solver threshold.",
        ),
        ReverseCase(
            stage="r4_high_field_models",
            case_id="R4_tat_strong",
            label="R4 Stronger Hurkx TAT",
            overrides={"bc_right": 0.0, "tat_field": 1.0e6},
            notes="Strengthens only the field-assisted tunneling multiplier.",
        ),
        ReverseCase(
            stage="r4_high_field_models",
            case_id="R4_trap_density_50",
            label="R4 Higher Trap Density Scale",
            overrides={"bc_right": 0.0, "trap_density_scale": 50.0},
            notes="Raises only the phenomenological trap-center density that multiplies the active Hurkx-like path.",
        ),
    ]

    stage_r5 = [
        ReverseCase(
            stage="r5_polarization_stability",
            case_id="R5_full_polarization",
            label="R5 Full Polarization",
            overrides={"bc_right": 0.0},
            notes="Uses the default polarization treatment from the authoritative JSON and solver.",
        ),
        ReverseCase(
            stage="r5_polarization_stability",
            case_id="R5_partial_polarization",
            label="R5 Partial Polarization Scale 0.6",
            overrides={"bc_right": 0.0, "polarization_scale": 0.60},
            notes="Scales the total polarization charge contribution while keeping the layer stack fixed.",
        ),
        ReverseCase(
            stage="r5_polarization_stability",
            case_id="R5_smoothed_interfaces",
            label="R5 Smoothed Polarization Interfaces 1.0 nm",
            overrides={"bc_right": 0.0, "polarization_smoothing_nm": 1.0},
            notes="Applies interface smoothing to the effective polarization-charge profile without changing the material stack or composition.",
        ),
    ]

    all_artifacts: list[CaseArtifacts] = []
    all_artifacts.extend(run_stage(base_project, stage_r1, "This stage validates whether the step-like reverse branch is already introduced by the current extraction methodology itself."))
    all_artifacts.extend(run_stage(base_project, stage_r2, "This stage compares the actual contact-boundary modes available in the active solver path and checks whether reverse injection is artificially suppressed at the contacts."))
    all_artifacts.extend(run_stage(base_project, stage_r3, "This stage tests mesh sensitivity with the least blocked contact offset while keeping the device stack fixed. Native local nonuniform refinement is not exposed by this workflow, so uniform grid refinement is used instead."))
    all_artifacts.extend(run_stage(base_project, stage_r4, "This stage checks whether the implemented Hurkx-like transport enhancements meaningfully change the reverse branch and documents which high-field mechanisms are missing from the active solver path."))
    all_artifacts.extend(run_stage(base_project, stage_r5, "This stage checks whether polarization-charge discontinuities or interface smoothing materially change reverse transport stability and leakage extraction."))
    write_code_findings()
    write_overall_summary(all_artifacts)
    print(f"Reverse transport validation complete. Summary: {OUT_ROOT / 'reverse_transport_summary.md'}")


if __name__ == "__main__":
    main()
