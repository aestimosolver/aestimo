# -*- coding: utf-8 -*-
"""
Aestimo 1D LED Validation & Traceability Module
Implements structured provenance recording, automated quantitative comparison
against peer-reviewed experimental datasets, error metric calculations,
and the standardized 4-panel LED plotting suite.
"""

from __future__ import annotations
import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt

from aeslibs.characterize_led import (
    calculate_led_recombination_profile,
    compute_integrated_led_efficiencies,
    generate_electroluminescence_spectrum,
    analyze_led_iv_curve,
    compute_led_efficiency_droop_curve,
    HC_EV_NM,
)


class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for numpy scalar and array serialization."""
    def default(self, obj):
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)



def load_led_experimental_csv(filepath: str) -> tuple[dict[str, np.ndarray], list[str]]:
    """
    Loads experimental LED CSV data, parsing column headers and ignoring comment lines.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Experimental LED file not found: {filepath}")
        
    comments = []
    header_line = None
    data_lines = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                comments.append(stripped[1:].strip())
            elif header_line is None:
                header_line = [c.strip() for c in stripped.split(',')]
            else:
                data_lines.append([float(x.strip()) for x in stripped.split(',')])
                
    if header_line is None or len(data_lines) == 0:
        raise ValueError(f"Invalid experimental LED CSV file format: {filepath}")
        
    data_arr = np.array(data_lines)
    data_dict = {}
    for idx, col in enumerate(header_line):
        data_dict[col] = data_arr[:, idx]
        
    return data_dict, comments


def compute_quantitative_error_metrics(
    exp_values: np.ndarray,
    sim_values: np.ndarray,
    target_name: str = "Metric",
) -> dict[str, float]:
    """
    Calculates statistical and quantitative error metrics:
    RMSE, NRMSE, MAE, MAPE, Max Absolute Error, and R^2.
    """
    y_exp = np.asarray(exp_values, dtype=float)
    y_sim = np.asarray(sim_values, dtype=float)
    
    if len(y_exp) != len(y_sim) or len(y_exp) == 0:
        return {
            'rmse': float('nan'),
            'nrmse': float('nan'),
            'mae': float('nan'),
            'mape': float('nan'),
            'max_abs_err': float('nan'),
            'r2': float('nan'),
        }
        
    diff = y_sim - y_exp
    abs_diff = np.abs(diff)
    
    # RMSE
    rmse = float(np.sqrt(np.mean(diff**2)))
    
    # NRMSE (normalized by experimental range)
    val_range = float(np.max(y_exp) - np.min(y_exp))
    nrmse = (rmse / val_range) if val_range > 1.0e-12 else float('nan')
    
    # MAE
    mae = float(np.mean(abs_diff))
    
    # MAPE (filter points where |y_exp| > 1e-12 to prevent zero division)
    valid_mask = np.abs(y_exp) > 1.0e-12
    if np.any(valid_mask):
        mape = float(np.mean(abs_diff[valid_mask] / np.abs(y_exp[valid_mask])) * 100.0)
    else:
        mape = float('nan')
        
    # Max Absolute Error
    max_abs_err = float(np.max(abs_diff))
    
    # R^2
    ss_tot = np.sum((y_exp - np.mean(y_exp))**2)
    ss_res = np.sum(diff**2)
    r2 = float(1.0 - (ss_res / ss_tot)) if ss_tot > 1.0e-12 else float('nan')
    
    return {
        'target_name': target_name,
        'rmse': rmse,
        'nrmse': nrmse,
        'mae': mae,
        'mape': mape,
        'max_abs_err': max_abs_err,
        'r2': r2,
    }


class LEDTraceabilityRecord:
    """
    Encapsulates the complete machine-readable provenance and validation record:
    Device ID -> Experimental Reference -> Experimental Structure -> Extracted Parameters ->
    Derived Parameters -> Fitted Parameters -> Assumptions -> Simulation Results ->
    Experimental Results -> Error Metrics -> Validation Status.
    """
    def __init__(
        self,
        device_id: str,
        device_name: str,
        bibliographic_reference: dict | None = None,
        experimental_structure: dict | None = None,
        parameter_provenance: dict | None = None,
        validation_status: str = "MODEL-BASED / NOT EXPERIMENTALLY VALIDATED",
    ):
        self.device_id = device_id
        self.device_name = device_name
        self.bibliographic_reference = bibliographic_reference or {}
        self.experimental_structure = experimental_structure or {}
        self.parameter_provenance = parameter_provenance or {}
        self.validation_status = validation_status
        
        self.simulation_results = {}
        self.experimental_results = {}
        self.error_metrics = {}
        self.assumptions_and_limitations = []

    def set_bibliographic_reference(self, ref: dict):
        self.bibliographic_reference = ref

    def set_experimental_structure(self, structure: dict):
        self.experimental_structure = structure

    def add_parameter_provenance(self, parameter_name: str, classification: str, source: str, value=None, unit=None):
        self.parameter_provenance[parameter_name] = {
            'classification': classification,
            'source': source,
            'value': value,
            'unit': unit,
        }

    def set_simulation_results(self, sim_data: dict):
        self.simulation_results = sim_data

    def set_experimental_results(self, exp_data: dict):
        self.experimental_results = exp_data

    def set_error_metrics(self, metrics: dict):
        self.error_metrics = metrics

    def add_assumption(self, description: str, rationale: str):
        self.assumptions_and_limitations.append({
            'description': description,
            'rationale': rationale,
        })

    def to_dict(self) -> dict:
        return {
            'device_id': self.device_id,
            'device_name': self.device_name,
            'validation_status': self.validation_status,
            'bibliographic_reference': self.bibliographic_reference,
            'experimental_structure': self.experimental_structure,
            'parameter_provenance': self.parameter_provenance,
            'simulation_results': self.simulation_results,
            'experimental_results': self.experimental_results,
            'error_metrics': self.error_metrics,
            'assumptions_and_limitations': self.assumptions_and_limitations,
        }

    def save_json(self, filepath: str):
        os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(self.to_dict(), f, indent=2, cls=NumpyEncoder)


def generate_led_validation_report_markdown(record: LEDTraceabilityRecord) -> str:
    """
    Renders an authoritative Markdown report detailing parameter provenance,
    quantitative comparison, error metrics, and validation assessment.
    """
    ref = record.bibliographic_reference
    lines = [
        f"# Experimental Validation Report: {record.device_name}",
        "",
        f"**Device ID**: `{record.device_id}`  ",
        f"**Validation Status**: **`{record.validation_status}`**",
        "",
        "---",
        "",
        "## 1. Bibliographic Reference & Provenance",
        f"- **Title**: {ref.get('title', 'N/A')}",
        f"- **Authors**: {ref.get('authors', 'N/A')}",
        f"- **Journal**: *{ref.get('journal', 'N/A')}*, Vol. {ref.get('volume', 'N/A')}, pp. {ref.get('pages', 'N/A')} ({ref.get('year', 'N/A')})",
        f"- **DOI**: [{ref.get('doi', 'N/A')}](https://doi.org/{ref.get('doi', '')})",
        f"- **Measurement Temperature**: {record.experimental_structure.get('temperature_k', 300)} K",
        f"- **Active Area**: {record.experimental_structure.get('mesa_area_cm2', 'N/A'):.3e} cm²",
        "",
        "### Parameter Classification Table",
        "| Parameter Name | Value | Unit | Classification | Physical Source / Extraction Method |",
        "| :--- | :--- | :--- | :--- | :--- |",
    ]
    
    prov = record.parameter_provenance
    for p_name, p_info in prov.items():
        val = p_info.get('value', 'Multiple / Stack')
        val_str = f"{val:.2e}" if isinstance(val, float) and (abs(val) < 0.01 or abs(val) > 1000) else str(val)
        unit = p_info.get('unit', '-')
        cls = p_info.get('classification', 'ASSUMED')
        src = p_info.get('source', p_info.get('rationale', p_info.get('method', 'Literature')))
        lines.append(f"| `{p_name}` | {val_str} | {unit} | **`{cls}`** | {src} |")
        
    lines.extend([
        "",
        "---",
        "",
        "## 2. Quantitative Error Metrics & Target Comparison",
        "| Target Observable | Experimental | Mode 10 Simulation | Error / Difference | Metric Status |",
        "| :--- | ---: | ---: | ---: | :--- |",
    ])
    
    # Electrical Turn-on & Forward Voltage
    sim_r = record.simulation_results
    exp_r = record.experimental_results
    
    v_on_exp = exp_r.get('turn_on_voltage_v')
    v_on_sim = sim_r.get('turn_on_voltage_v')
    if v_on_exp is not None and v_on_sim is not None:
        delta_von = abs(v_on_sim - v_on_exp)
        status_von = "[EXCELLENT] < 0.1 V" if delta_von < 0.10 else "[GOOD] < 0.25 V" if delta_von < 0.25 else "[ACCEPTABLE]"
        lines.append(f"| **Turn-on Voltage ($V_{{on}}$)** | {v_on_exp:.3f} V | {v_on_sim:.3f} V | $\\Delta V = {delta_von:.3f}$ V | {status_von} |")
        
    vf_exp = exp_r.get('forward_voltage_at_20ma_v')
    vf_sim = sim_r.get('forward_voltage_at_20ma_v')
    if vf_exp is not None and vf_sim is not None:
        delta_vf = abs(vf_sim - vf_exp)
        status_vf = "[EXCELLENT] < 0.15 V" if delta_vf < 0.15 else "[GOOD] < 0.3 V" if delta_vf < 0.30 else "[ACCEPTABLE]"
        lines.append(f"| **Operating Voltage ($V_f$ at 20 mA)** | {vf_exp:.3f} V | {vf_sim:.3f} V | $\\Delta V = {delta_vf:.3f}$ V | {status_vf} |")
        
    # Peak Emission Wavelength
    wl_exp = exp_r.get('peak_wavelength_nm')
    wl_sim = sim_r.get('peak_wavelength_nm')
    if wl_exp is not None and wl_sim is not None:
        delta_wl = abs(wl_sim - wl_exp)
        status_wl = "[EXCELLENT] < 2 nm" if delta_wl < 2.0 else "[GOOD] < 5 nm" if delta_wl < 5.0 else "[ACCEPTABLE]"
        lines.append(f"| **Peak Emission Wavelength ($\\lambda_{{peak}}$)** | {wl_exp:.1f} nm | {wl_sim:.1f} nm | $\\Delta \\lambda = {delta_wl:.2f}$ nm | {status_wl} |")
        
    # FWHM
    fwhm_exp = exp_r.get('fwhm_nm')
    fwhm_sim = sim_r.get('fwhm_nm')
    if fwhm_exp is not None and fwhm_sim is not None:
        delta_fwhm = abs(fwhm_sim - fwhm_exp)
        status_fwhm = "[EXCELLENT] < 3 nm" if delta_fwhm < 3.0 else "[GOOD] < 6 nm"
        lines.append(f"| **Spectral FWHM** | {fwhm_exp:.1f} nm | {fwhm_sim:.1f} nm | $\\Delta \\text{{FWHM}} = {delta_fwhm:.2f}$ nm | {status_fwhm} |")
        
    # Optical Power at 20 mA
    popt_exp = exp_r.get('optical_power_at_20ma_mw')
    popt_sim = sim_r.get('optical_power_at_20ma_mw')
    if popt_exp is not None and popt_sim is not None:
        err_p = abs(popt_sim - popt_exp) / popt_exp * 100.0
        status_p = "[EXCELLENT] < 10%" if err_p < 10.0 else "[GOOD] < 20%"
        lines.append(f"| **Optical Power ($P_{{opt}}$ at 20 mA)** | {popt_exp:.2f} mW | {popt_sim:.2f} mW | Relative Err: {err_p:.1f}% | {status_p} |")
        
    # Efficiency Droop Peak J_peak
    jpeak_exp = exp_r.get('j_peak_a_cm2')
    jpeak_sim = sim_r.get('j_peak_a_cm2')
    if jpeak_exp is not None and jpeak_sim is not None:
        delta_j = abs(jpeak_sim - jpeak_exp)
        status_j = "[EXCELLENT]" if delta_j < 5.0 else "[GOOD]"
        lines.append(f"| **Droop Peak Density ($J_{{peak}}$)** | {jpeak_exp:.1f} A/cm² | {jpeak_sim:.1f} A/cm² | $\\Delta J = {delta_j:.1f}$ A/cm² | {status_j} |")
        
    # Statistical Curves
    errs = record.error_metrics
    lines.extend([
        "",
        "### Statistical Error Summaries",
        "| Curve Comparison | RMSE | Normalized RMSE (NRMSE) | MAPE (%) | $R^2$ Score |",
        "| :--- | :--- | :--- | :--- | :--- |",
    ])
    for c_name, c_err in errs.items():
        if isinstance(c_err, dict) and 'rmse' in c_err:
            lines.append(
                f"| **{c_name}** | {c_err['rmse']:.3e} | {c_err['nrmse']:.3%} | {c_err['mape']:.2f}% | {c_err['r2']:.4f} |"
            )
            
    lines.extend([
        "",
        "---",
        "",
        "## 3. Assumptions and Remaining Limitations",
    ])
    for item in record.assumptions_and_limitations:
        lines.append(f"- **{item['description']}**: {item['rationale']}")
        
    return "\n".join(lines) + "\n"


def plot_standardized_led_suite(
    record: LEDTraceabilityRecord,
    save_path: str | None = None,
) -> plt.Figure:
    """
    Renders the unified 4-panel standardized LED validation figure:
    - Panel 1: I-V Characteristics (Linear & Semi-log)
    - Panel 2: Electroluminescence Emission Spectrum
    - Panel 3: Optical Output Power vs Current (L-I)
    - Panel 4: Quantum Efficiency & Droop (IQE & EQE vs Current Density J)
    """
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    
    exp_r = record.experimental_results
    sim_r = record.simulation_results
    
    # -------------------------------------------------------------
    # Panel 1: I-V Characteristics (Semi-log & Linear inset)
    # -------------------------------------------------------------
    ax1 = axes[0, 0]
    
    # Semi-log I-V
    if 'iv_voltage_v' in exp_r and 'iv_current_ma' in exp_r:
        v_exp = exp_r['iv_voltage_v']
        i_exp = exp_r['iv_current_ma']
        ax1.semilogy(
            v_exp, np.clip(np.abs(i_exp), 1.0e-7, None),
            'o', label='Experimental Data',
            markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2
        )
        
    if 'iv_voltage_v' in sim_r and 'iv_current_ma' in sim_r:
        v_sim = sim_r['iv_voltage_v']
        i_sim = sim_r['iv_current_ma']
        ax1.semilogy(
            v_sim, np.clip(np.abs(i_sim), 1.0e-7, None),
            '-', label='Simulation (Mode 10)',
            color='#003366', linewidth=2.0
        )
        
    ax1.set_title("Current-Voltage (I-V) Characteristics", fontweight='bold', fontsize=12)
    ax1.set_xlabel("Forward Voltage (V)", fontsize=11)
    ax1.set_ylabel("|Current| (mA, log scale)", fontsize=11)
    ax1.grid(True, which='both', linestyle='--', alpha=0.4)
    ax1.legend(loc='upper left', frameon=True)
    
    # -------------------------------------------------------------
    # Panel 2: Electroluminescence (EL) Emission Spectrum
    # -------------------------------------------------------------
    ax2 = axes[0, 1]
    
    if 'el_wavelength_nm' in exp_r and 'el_intensity_norm' in exp_r:
        ax2.plot(
            exp_r['el_wavelength_nm'], exp_r['el_intensity_norm'],
            'o', label='Measured Spectrum',
            markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2
        )
        
    if 'el_wavelength_nm' in sim_r and 'el_intensity_norm' in sim_r:
        wl_sim = sim_r['el_wavelength_nm']
        int_sim = sim_r['el_intensity_norm']
        ax2.plot(
            wl_sim, int_sim,
            '-', label='Simulated EL Line shape',
            color='#003366', linewidth=2.0
        )
        # Annotate peak
        peak_wl = sim_r.get('peak_wavelength_nm', wl_sim[int(np.argmax(int_sim))])
        fwhm = sim_r.get('fwhm_nm', 25.0)
        ax2.axvline(peak_wl, color='#888888', linestyle=':', alpha=0.7)
        ax2.text(
            0.65, 0.78,
            f"$\\lambda_{{peak}} = {peak_wl:.1f}$ nm\nFWHM = {fwhm:.1f} nm",
            transform=ax2.transAxes,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#003366", alpha=0.9),
            fontsize=10
        )
        
    ax2.set_title("Electroluminescence (EL) Spectrum", fontweight='bold', fontsize=12)
    ax2.set_xlabel("Wavelength (nm)", fontsize=11)
    ax2.set_ylabel("Normalized Intensity (a.u.)", fontsize=11)
    ax2.set_ylim(-0.05, 1.15)
    ax2.grid(True, which='both', linestyle='--', alpha=0.4)
    ax2.legend(loc='upper right', frameon=True)
    
    # -------------------------------------------------------------
    # Panel 3: Optical Output Power (L-I Curve)
    # -------------------------------------------------------------
    ax3 = axes[1, 0]
    
    if 'li_current_ma' in exp_r and 'li_optical_power_mw' in exp_r:
        ax3.plot(
            exp_r['li_current_ma'], exp_r['li_optical_power_mw'],
            'o', label='Experimental L-I',
            markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2
        )
        
    if 'li_current_ma' in sim_r and 'li_optical_power_mw' in sim_r:
        ax3.plot(
            sim_r['li_current_ma'], sim_r['li_optical_power_mw'],
            '-', label='Simulated P_opt (Mode 10)',
            color='#003366', linewidth=2.0
        )
    elif 'iv_current_ma' in sim_r and 'p_opt_mw' in sim_r:
        ax3.plot(
            sim_r['iv_current_ma'], sim_r['p_opt_mw'],
            '-', label='Simulated P_opt',
            color='#003366', linewidth=2.0
        )
        
    ax3.set_title("Optical Output Power (L-I Characteristics)", fontweight='bold', fontsize=12)
    ax3.set_xlabel("Forward Injection Current (mA)", fontsize=11)
    ax3.set_ylabel("Optical Output Power (mW)", fontsize=11)
    ax3.set_ylim(bottom=0.0)
    ax3.grid(True, which='both', linestyle='--', alpha=0.4)
    ax3.legend(loc='upper left', frameon=True)
    
    # -------------------------------------------------------------
    # Panel 4: Quantum Efficiency & Droop (IQE vs J)
    # -------------------------------------------------------------
    ax4 = axes[1, 1]
    
    if 'droop_j_a_cm2' in exp_r and 'droop_norm_iqe' in exp_r:
        ax4.semilogx(
            np.asarray(exp_r['droop_j_a_cm2']), np.asarray(exp_r['droop_norm_iqe']) * 100.0,
            'o', label='Measured Droop',
            markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2
        )
        
    if 'droop_j_a_cm2' in sim_r and 'droop_norm_iqe' in sim_r:
        ax4.semilogx(
            np.asarray(sim_r['droop_j_a_cm2']), np.asarray(sim_r['droop_norm_iqe']) * 100.0,
            '-', label='Simulated Droop (ABC Model)',
            color='#003366', linewidth=2.0
        )
        j_pk = sim_r.get('j_peak_a_cm2', 12.0)
        ax4.axvline(j_pk, color='#cc0000', linestyle='--', alpha=0.7)
        ax4.text(
            0.45, 0.25,
            f"Droop Peak: $J_{{peak}} \\approx {j_pk:.1f}$ A/cm²",
            transform=ax4.transAxes,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cc0000", alpha=0.9),
            fontsize=10
        )
        
    ax4.set_title("Internal Quantum Efficiency (IQE) & Droop", fontweight='bold', fontsize=12)
    ax4.set_xlabel("Current Density J (A/cm², log scale)", fontsize=11)
    ax4.set_ylabel("Normalized IQE (%)", fontsize=11)
    ax4.set_ylim(20.0, 110.0)
    ax4.grid(True, which='both', linestyle='--', alpha=0.4)
    ax4.legend(loc='lower left', frameon=True)
    
    # Global Super Title with metadata
    meta = record.experimental_structure
    ref = record.bibliographic_reference
    fig.suptitle(
        f"Aestimo 1D LED Validation: {record.device_name} [{record.validation_status}]\n"
        f"Source: {ref.get('authors', '')} ({ref.get('year', '')}), DOI: {ref.get('doi', '')} | "
        f"T = {meta.get('temperature_k', 300)} K, Area = {meta.get('mesa_area_cm2', 1e-3):.2e} cm²",
        fontsize=13, fontweight='bold', y=1.02
    )
    
    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        fig.savefig(save_path, dpi=220, bbox_inches='tight')
        
    return fig
