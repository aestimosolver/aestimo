from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
REVERSE_ROOT = ROOT / "results" / "reverse_transport"
OUT_ROOT = REVERSE_ROOT / "r6_terminal_observable"

# The reverse-transport validation tables were exported using the historical
# A/m^2 -> A/cm^2 conversion. The active scheme-7 path multiplies the DD current
# by 0.1 before saving, so those arrays are best interpreted as mA/cm^2.
CORRECTED_SCALE = 10.0
CURRENT_COLUMNS = (
    "official_a_cm2",
    "left_probe_a_cm2",
    "center_probe_a_cm2",
    "right_probe_a_cm2",
    "whole_probe_a_cm2",
    "absmax_probe_a_cm2",
    "left_edge_a_cm2",
    "right_edge_a_cm2",
)
PROFILE_CURRENT_COLUMNS = ("Jtotal_a_cm2", "Jelec_a_cm2", "Jhole_a_cm2")
REPRESENTATIVE_CASES = (
    "R1_selective_bc0p0",
    "R2_ohmic_bc0p6",
    "R4_tat_disabled",
    "R4_tat_strong",
    "R5_smoothed_interfaces",
)


@dataclass
class CaseAnalysis:
    stage: str
    case_id: str
    case_dir: Path
    table: np.ndarray
    table_headers: list[str]
    metrics: dict[str, float | int | str | bool | None]


def ensure_out() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)


def iter_case_dirs() -> Iterable[Path]:
    for stage_dir in sorted(REVERSE_ROOT.iterdir()):
        if not stage_dir.is_dir() or not stage_dir.name.startswith("r"):
            continue
        if stage_dir.name == OUT_ROOT.name:
            continue
        for case_dir in sorted(stage_dir.iterdir()):
            if (case_dir / "raw_current_table.csv").exists():
                yield case_dir


def load_csv_table(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        headers = next(reader)
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    return headers, data


def corrected_current(data: np.ndarray, column: str) -> np.ndarray:
    return np.asarray(data[column], dtype=float) * CORRECTED_SCALE


def interp_at(data: np.ndarray, column: str, bias_v: float) -> float:
    return float(np.interp(bias_v, np.asarray(data["bias_v"], dtype=float), corrected_current(data, column)))


def profile_paths(case_dir: Path) -> list[Path]:
    return sorted(case_dir.glob("profile_*.csv"))


def load_profile(path: Path) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",", names=True, dtype=float)


def profile_bias_from_name(path: Path) -> float | None:
    stem = path.stem.removeprefix("profile_")
    sign = -1.0 if stem.startswith("m") else 1.0
    stem = stem[1:] if stem[:1] in {"m", "p"} else stem
    stem = stem.removesuffix("v").replace("p", ".")
    try:
        return sign * float(stem)
    except ValueError:
        return None


def conservation_metrics(profile: np.ndarray) -> dict[str, float | int]:
    j = np.asarray(profile["Jtotal_a_cm2"], dtype=float) * CORRECTED_SCALE
    abs_j = np.abs(j)
    finite = np.isfinite(j)
    if not np.any(finite):
        return {
            "profile_points": int(j.size),
            "j_median_a_cm2": float("nan"),
            "j_absmax_a_cm2": float("nan"),
            "j_p95_a_cm2": float("nan"),
            "j_p05_a_cm2": float("nan"),
            "zero_like_fraction": float("nan"),
            "nonnegative_fraction": float("nan"),
            "conservation_ratio": float("nan"),
        }
    j_f = j[finite]
    abs_f = abs_j[finite]
    absmax = float(np.max(abs_f))
    absmed = float(np.median(abs_f))
    return {
        "profile_points": int(j_f.size),
        "j_median_a_cm2": float(np.median(j_f)),
        "j_absmax_a_cm2": absmax,
        "j_p95_a_cm2": float(np.percentile(j_f, 95)),
        "j_p05_a_cm2": float(np.percentile(j_f, 5)),
        "zero_like_fraction": float(np.mean(abs_f <= 1e-18)),
        "nonnegative_fraction": float(np.mean(j_f >= 0.0)),
        "conservation_ratio": float(absmax / max(absmed, 1e-30)),
    }


def analyze_case(case_dir: Path) -> CaseAnalysis:
    headers, data = load_csv_table(case_dir / "raw_current_table.csv")
    metrics_path = case_dir / "metrics.json"
    prior_metrics = json.loads(metrics_path.read_text(encoding="utf-8")) if metrics_path.exists() else {}
    stage = case_dir.parent.name
    case_id = case_dir.name

    metrics: dict[str, float | int | str | bool | None] = {
        "stage": stage,
        "case_id": case_id,
        "official_j_m1_v_a_cm2": interp_at(data, "official_a_cm2", -1.0),
        "whole_j_m1_v_a_cm2": interp_at(data, "whole_probe_a_cm2", -1.0),
        "left_j_m1_v_a_cm2": interp_at(data, "left_probe_a_cm2", -1.0),
        "right_j_m1_v_a_cm2": interp_at(data, "right_probe_a_cm2", -1.0),
        "left_edge_j_m1_v_a_cm2": interp_at(data, "left_edge_a_cm2", -1.0),
        "right_edge_j_m1_v_a_cm2": interp_at(data, "right_edge_a_cm2", -1.0),
        "absmax_j_m1_v_a_cm2": interp_at(data, "absmax_probe_a_cm2", -1.0),
        "official_zero_like_points": int(np.count_nonzero(np.abs(corrected_current(data, "official_a_cm2")) <= 1e-18)),
        "whole_zero_like_points": int(np.count_nonzero(np.abs(corrected_current(data, "whole_probe_a_cm2")) <= 1e-18)),
        "timeout_events_total": int(prior_metrics.get("timeout_events_total", 0)),
        "max_gummel_iterations": int(prior_metrics.get("max_gummel_iterations", 0)),
        "all_bias_targets_reached": bool(prior_metrics.get("all_bias_targets_reached", False)),
    }

    closest_profile_path: Path | None = None
    closest_distance = float("inf")
    for path in profile_paths(case_dir):
        bias = profile_bias_from_name(path)
        if bias is None:
            continue
        distance = abs(bias + 1.0)
        if distance < closest_distance:
            closest_distance = distance
            closest_profile_path = path
    if closest_profile_path is not None:
        profile = load_profile(closest_profile_path)
        metrics.update({f"profile_m1_{key}": value for key, value in conservation_metrics(profile).items()})
        metrics["profile_m1_source"] = closest_profile_path.name
    else:
        metrics["profile_m1_source"] = None

    return CaseAnalysis(stage=stage, case_id=case_id, case_dir=case_dir, table=data, table_headers=headers, metrics=metrics)


def write_summary_csv(cases: list[CaseAnalysis]) -> Path:
    path = OUT_ROOT / "terminal_observable_summary.csv"
    headers = [
        "stage",
        "case_id",
        "official_j_m1_v_a_cm2",
        "whole_j_m1_v_a_cm2",
        "left_j_m1_v_a_cm2",
        "right_j_m1_v_a_cm2",
        "left_edge_j_m1_v_a_cm2",
        "right_edge_j_m1_v_a_cm2",
        "absmax_j_m1_v_a_cm2",
        "official_zero_like_points",
        "whole_zero_like_points",
        "profile_m1_source",
        "profile_m1_j_median_a_cm2",
        "profile_m1_j_absmax_a_cm2",
        "profile_m1_j_p95_a_cm2",
        "profile_m1_j_p05_a_cm2",
        "profile_m1_zero_like_fraction",
        "profile_m1_conservation_ratio",
        "timeout_events_total",
        "max_gummel_iterations",
        "all_bias_targets_reached",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for case in cases:
            writer.writerow({key: case.metrics.get(key, "") for key in headers})
    return path


def save_corrected_current_tables(cases: list[CaseAnalysis]) -> None:
    for case in cases:
        out_dir = OUT_ROOT / "corrected_current_tables" / case.stage
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"{case.case_id}_corrected_current_table.csv"
        rows = []
        for row_index in range(case.table.shape[0]):
            row = {}
            for header in case.table_headers:
                value = case.table[header][row_index]
                row[header] = value * CORRECTED_SCALE if header in CURRENT_COLUMNS else value
            rows.append(row)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=case.table_headers)
            writer.writeheader()
            writer.writerows(rows)


def plot_representative_reverse_iv(cases: list[CaseAnalysis]) -> Path:
    path = OUT_ROOT / "representative_reverse_observables.png"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    selected = [case for case in cases if case.case_id in REPRESENTATIVE_CASES]
    for case in selected:
        bias = np.asarray(case.table["bias_v"], dtype=float)
        official = np.abs(corrected_current(case.table, "official_a_cm2"))
        right_edge = np.abs(corrected_current(case.table, "right_edge_a_cm2"))
        absmax = np.abs(corrected_current(case.table, "absmax_probe_a_cm2"))
        axes[0].semilogy(bias, np.maximum(official, 1e-30), marker="o", linewidth=1.2, label=case.case_id)
        axes[1].semilogy(bias, np.maximum(right_edge, 1e-30), marker="o", linewidth=1.2, label=f"{case.case_id} edge")
        axes[1].semilogy(bias, np.maximum(absmax, 1e-30), linestyle="--", linewidth=1.1, alpha=0.75, label=f"{case.case_id} absmax")
    axes[0].set_title("Official terminal observable")
    axes[1].set_title("Edge and internal-current proxies")
    for axis in axes:
        axis.set_xlabel("Bias (V)")
        axis.set_ylabel("|J| (A/cm$^2$)")
        axis.grid(True, which="both", alpha=0.25)
        axis.set_xlim(-1.05, 0.0)
    axes[1].legend(fontsize=7, loc="best")
    axes[0].legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_spatial_current_profiles(cases: list[CaseAnalysis]) -> Path:
    path = OUT_ROOT / "representative_spatial_current_profiles_m1v.png"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for case in cases:
        if case.case_id not in REPRESENTATIVE_CASES:
            continue
        profile_path = case.case_dir / "profile_m1p00v.csv"
        if not profile_path.exists():
            continue
        profile = load_profile(profile_path)
        x_nm = np.asarray(profile["x_nm"], dtype=float)
        j = np.asarray(profile["Jtotal_a_cm2"], dtype=float) * CORRECTED_SCALE
        axes[0].plot(x_nm, j, linewidth=1.3, label=case.case_id)
        axes[1].semilogy(x_nm, np.maximum(np.abs(j), 1e-30), linewidth=1.3, label=case.case_id)
    axes[0].set_title("Signed spatial current")
    axes[1].set_title("Absolute spatial current")
    for axis in axes:
        axis.set_xlabel("Position (nm)")
        axis.set_ylabel("Jtotal (A/cm$^2$)")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def write_report(cases: list[CaseAnalysis], summary_csv: Path, iv_plot: Path, profile_plot: Path) -> Path:
    def find(case_id: str) -> CaseAnalysis:
        return next(case for case in cases if case.case_id == case_id)

    r1 = find("R1_selective_bc0p0")
    r2 = find("R2_ohmic_bc0p6")
    r4_off = find("R4_tat_disabled")
    r4_on = find("R4_tat_strong")
    r5 = find("R5_smoothed_interfaces")

    lines = [
        "# R6 Terminal-Current Observable And Unit Audit",
        "",
        "## Scope",
        "",
        "This follow-up reuses the completed reverse-transport validation outputs. It does not introduce new device sweeps or change the authoritative JSON structure. Its purpose is to determine whether the saved solver state contains a physically usable reverse terminal-current observable.",
        "",
        "## Unit Audit",
        "",
        "- The active scheme-7 current routine computes `Jelec` and `Jhole` in SI `A/m^2`, then the scheme-7 path multiplies them by `0.1` before `Jtotal`, `av_curr`, and the saved result arrays are formed.",
        "- Therefore the native saved current array is best interpreted as `mA/cm^2`; converting it to `A/cm^2` requires multiplying by `1e-3`.",
        "- The earlier reverse-validation export used the historical `A/m^2 -> A/cm^2` factor of `1e-4`, so the previously reported nonzero current magnitudes are low by a factor of 10. Exact-zero and step-like behavior are unchanged.",
        "",
        "## Corrected Representative Values At -1 V",
        "",
        f"- `R1_selective_bc0p0`: official `{r1.metrics['official_j_m1_v_a_cm2']:.3e}` A/cm^2, whole median `{r1.metrics['whole_j_m1_v_a_cm2']:.3e}` A/cm^2, abs-max internal proxy `{r1.metrics['absmax_j_m1_v_a_cm2']:.3e}` A/cm^2.",
        f"- `R2_ohmic_bc0p6`: official `{r2.metrics['official_j_m1_v_a_cm2']:.3e}` A/cm^2, whole median `{r2.metrics['whole_j_m1_v_a_cm2']:.3e}` A/cm^2, abs-max internal proxy `{r2.metrics['absmax_j_m1_v_a_cm2']:.3e}` A/cm^2.",
        f"- `R4_tat_disabled`: official `{r4_off.metrics['official_j_m1_v_a_cm2']:.3e}` A/cm^2, whole median `{r4_off.metrics['whole_j_m1_v_a_cm2']:.3e}` A/cm^2, abs-max internal proxy `{r4_off.metrics['absmax_j_m1_v_a_cm2']:.3e}` A/cm^2.",
        f"- `R4_tat_strong`: official `{r4_on.metrics['official_j_m1_v_a_cm2']:.3e}` A/cm^2, whole median `{r4_on.metrics['whole_j_m1_v_a_cm2']:.3e}` A/cm^2, abs-max internal proxy `{r4_on.metrics['absmax_j_m1_v_a_cm2']:.3e}` A/cm^2.",
        f"- `R5_smoothed_interfaces`: official `{r5.metrics['official_j_m1_v_a_cm2']:.3e}` A/cm^2, whole median `{r5.metrics['whole_j_m1_v_a_cm2']:.3e}` A/cm^2, abs-max internal proxy `{r5.metrics['absmax_j_m1_v_a_cm2']:.3e}` A/cm^2.",
        "",
        "## Current-Conservation Diagnosis",
        "",
        f"- In `R1_selective_bc0p0`, the -1 V spatial profile has median current `{r1.metrics['profile_m1_j_median_a_cm2']:.3e}` A/cm^2 and abs-max current `{r1.metrics['profile_m1_j_absmax_a_cm2']:.3e}` A/cm^2, giving a conservation ratio of `{r1.metrics['profile_m1_conservation_ratio']:.3e}`.",
        f"- In `R4_tat_disabled`, the -1 V spatial profile has median current `{r4_off.metrics['profile_m1_j_median_a_cm2']:.3e}` A/cm^2 and abs-max current `{r4_off.metrics['profile_m1_j_absmax_a_cm2']:.3e}` A/cm^2, giving a conservation ratio of `{r4_off.metrics['profile_m1_conservation_ratio']:.3e}`.",
        "- A steady one-dimensional terminal current should be nearly position independent. Ratios this large mean a contact or edge probe cannot be promoted to a physical terminal current without first fixing current continuity and boundary treatment.",
        "",
        "## Decision",
        "",
        "A corrected unit conversion increases the magnitude of nonzero currents, but it does not recover a gradual, conserved reverse leakage branch. The best current signal remains spatially localized and contact-sensitive rather than terminal and conserved.",
        "",
        "Conclusion B remains supported: the present solver/contact framework fundamentally limits quantitative high-indium reverse-bias leakage simulation in this workflow. The scientifically defensible path is to report forward/calibrated behavior with the existing solver only after validation, and to treat reverse leakage as a methodology limitation unless a new high-field/contact transport implementation is added.",
        "",
        "## Generated Artifacts",
        "",
        f"- Summary table: `{summary_csv.relative_to(ROOT)}`",
        f"- Representative reverse observables plot: `{iv_plot.relative_to(ROOT)}`",
        f"- Representative -1 V spatial current plot: `{profile_plot.relative_to(ROOT)}`",
        "- Corrected per-case current tables: `results/reverse_transport/r6_terminal_observable/corrected_current_tables/`",
    ]
    path = OUT_ROOT / "terminal_current_methodology_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_top_level_update(report: Path) -> Path:
    path = REVERSE_ROOT / "reverse_transport_methodology_recommendation.md"
    lines = [
        "# Reverse-Transport Methodology Recommendation",
        "",
        "The reverse-bias validation and the R6 terminal-current audit indicate that the current scheme-7 contact/current framework should not be used for quantitative high-indium InGaN reverse-leakage prediction.",
        "",
        "Key reasons:",
        "",
        "- The saved current units in the active scheme-7 path are best interpreted as `mA/cm^2`, so previous nonzero current magnitudes exported by the validation script are low by a factor of 10; this does not change the zero-current or step-onset failure mode.",
        "- The official `av_curr` observable can remain exactly zero at -1 V while internal current proxies are large.",
        "- Spatial current profiles at -1 V are not conserved well enough to define a robust terminal current from an alternate probe.",
        "- Contact boundary changes and TAT strength change internal currents, but they do not produce a stable, gradual, terminal reverse branch.",
        "",
        f"Detailed evidence is in `{report.relative_to(ROOT)}`.",
        "",
        "Recommended paper treatment:",
        "",
        "- Present reverse-bias results as a solver-methodology validation rather than as calibrated device predictions.",
        "- State that quantitative reverse leakage for high-indium InGaN requires a contact/high-field transport extension, such as physically parameterized thermionic-field emission, Poole-Frenkel or field-enhanced SRH, band-to-band tunneling where applicable, and a terminal-current formulation with current-conservation checks.",
        "- Keep later design sweeps restricted to observables that pass validation, or clearly label reverse-leakage trends as qualitative until the transport framework is extended.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    ensure_out()
    cases = [analyze_case(case_dir) for case_dir in iter_case_dirs()]
    cases = [case for case in cases if case.case_id != "SMOKE_baseline_json"]
    if not cases:
        raise RuntimeError(f"No reverse-transport cases found under {REVERSE_ROOT}")
    summary_csv = write_summary_csv(cases)
    save_corrected_current_tables(cases)
    iv_plot = plot_representative_reverse_iv(cases)
    profile_plot = plot_spatial_current_profiles(cases)
    report = write_report(cases, summary_csv, iv_plot, profile_plot)
    recommendation = write_top_level_update(report)
    print(f"Wrote {summary_csv}")
    print(f"Wrote {report}")
    print(f"Wrote {recommendation}")


if __name__ == "__main__":
    main()
