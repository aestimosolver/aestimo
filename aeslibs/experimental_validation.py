#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Experimental Data Validation Module for Aestimo

This module provides functionality to validate drift-diffusion simulation results
against experimental I-V and C-V measurements.

Functions:
- load_experimental_data: Load experimental data from CSV files
- calculate_current_from_simulation: Compute I-V from simulation carrier densities
- calculate_capacitance_from_simulation: Compute C-V from simulation results
- compute_error_metrics: Calculate validation metrics (RMSE, MAE, R²)
- plot_iv_comparison: Generate I-V comparison plots
- plot_cv_comparison: Generate C-V comparison plots
- generate_validation_report: Create comprehensive validation report
- calculate_ideality_factor: Extract local ideality factor n(V)
- apply_parasitic_resistances: Apply Rs and Rsh to simulated data
"""

import numpy as np
import os
from scipy import integrate, constants

# Physical constants
q = constants.e  # Elementary charge (C)
kb = constants.k  # Boltzmann constant (J/K)


def load_experimental_data(filepath):
    """
    Load experimental data from CSV file.
    
    Parameters:
    -----------
    filepath : str
        Path to CSV file containing experimental data
        Format: voltage (V), current (A) or capacitance (F)
        Lines starting with # are treated as comments
    
    Returns:
    --------
    voltage : np.ndarray
        Applied voltage values (V)
    data : np.ndarray
        Measured data values (A for I-V, F for C-V)
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Experimental data file not found: {filepath}")
    
    # Load data, skipping comment lines
    raw_data = np.loadtxt(filepath, delimiter=',', comments='#')
    
    voltage = raw_data[:, 0]
    data = raw_data[:, 1]
    
    return voltage, data

def load_experimental_metadata(filepath):
    """
    Extract metadata (Doping, Area) from experimental file header.
    
    Returns:
    --------
    metadata : dict
        Dictionary with keys 'area_cm2', 'doping_p', 'doping_n' if found.
    """
    metadata = {}
    if not os.path.exists(filepath):
        return metadata
        
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    break # Stop at first non-comment line
                
                # Parse Area
                if 'Area:' in line:
                    # Example: # Area: 1e-4 cm^2
                    parts = line.split('Area:')[1].split('cm^2')[0]
                    try:
                        metadata['area_cm2'] = float(parts.strip())
                    except:
                        pass
                
                # Parse Doping
                if 'Doping:' in line:
                    # Example: # Doping: p-side: 1e18 cm^-3, n-side: 1e18 cm^-3
                    # Simplified parsing logic
                    try:
                        line_content = line.split('Doping:')[1]
                        if 'p-side:' in line_content:
                            p_part = line_content.split('p-side:')[1].split(',')[0].replace('cm^-3','').strip()
                            metadata['doping_p'] = float(p_part)
                        if 'n-side:' in line_content:
                            n_part = line_content.split('n-side:')[1].replace('cm^-3','').strip()
                            metadata['doping_n'] = float(n_part)
                    except:
                        pass
    except Exception as e:
        print(f"Warning: Could not parse metadata from {filepath}: {e}")
        
    return metadata


def get_silicon_mobility(doping, carrier_type='n'):
    """
    Caughey-Thomas mobility model for Silicon at 300K.
    
    Parameters:
    -----------
    doping : np.ndarray
        Local doping concentration (cm^-3)
    carrier_type : str
        'n' for electrons, 'p' for holes
    """
    abs_doping = np.abs(doping)
    if carrier_type == 'n':
        # Parameters for electrons in Si
        mu_min, mu_max, N_ref, alpha = 92.0, 1360.0, 1.3e17, 0.91
    else:
        # Parameters for holes in Si
        mu_min, mu_max, N_ref, alpha = 47.7, 460.0, 6.3e16, 0.76
        
    return mu_min + (mu_max - mu_min) / (1.0 + (abs_doping / N_ref)**alpha)


def load_current_from_avcurr(output_dir, device_area_cm2=1e-4):
    """
    Load current directly from simulation output (av_curr.dat).
    
    This is more accurate than recalculating from carrier profiles because
    the simulation already computes current using the same discretization
    scheme internally.
    
    Parameters:
    -----------
    output_dir : str
        Path to simulation output directory
    device_area_cm2 : float
        Device cross-sectional area in cm² (default: 1e-4 cm² = 100x100 µm)
    
    Returns:
    --------
    voltages : np.ndarray
        Applied voltages (V)
    currents : np.ndarray
        Device current (A)
    """
    av_curr_file = os.path.join(output_dir, "av_curr.dat")
    
    if not os.path.exists(av_curr_file):
        raise FileNotFoundError(f"av_curr.dat not found in {output_dir}")
    
    data = np.loadtxt(av_curr_file)
    voltages = data[:, 0]
    J_density = data[:, 1]  # Active scheme-7 output is stored as mA/cm^2.
    
    # Convert mA/cm^2 to A/cm^2.
    J_cm2 = J_density * 1e-3
    
    # Convert to current: I = J × Area
    currents = J_cm2 * device_area_cm2
    
    # (Removed offset correction as it zeroes out the physical Jsc in light simulations)
    
    return voltages, currents


def apply_parasitic_resistances(voltages, currents, Rs=0.0, Rsh=1e12):
    """
    Apply series and shunt resistances to simulated I-V data.
    
    V_ext = V_int + I * Rs
    I_ext = I_int + V_ext / Rsh
    
    Parameters:
    -----------
    voltages : np.ndarray
        Internal simulation voltages (V)
    currents : np.ndarray
        Internal simulation currents (A)
    Rs : float
        Series resistance (Ohms)
    Rsh : float
        Shunt resistance (Ohms)
    
    Returns:
    --------
    v_ext : np.ndarray
        External (terminal) voltages (V)
    i_ext : np.ndarray
        External (terminal) currents (A)
    """
    if Rs == 0.0 and Rsh > 1e11:
        return voltages, currents
        
    i_ext = currents + voltages / Rsh
    v_ext = voltages + i_ext * Rs
    
    return v_ext, i_ext


def calculate_ideality_factor(voltages, currents, temperature=300):
    """
    Calculate the local ideality factor n(V) from I-V data.
    
    n = (q / (kb * T)) * (d(ln I) / dV)^-1
    
    Parameters:
    -----------
    voltages : np.ndarray
        Voltage values (V)
    currents : np.ndarray
        Current values (A)
    temperature : float
        Temperature (K)
        
    Returns:
    --------
    v_mid : np.ndarray
        Midpoint voltages (V)
    n : np.ndarray
        Local ideality factor
    """
    # Filter for forward bias and positive current
    mask = (voltages > 0.1) & (currents > 1e-12)
    v_f = voltages[mask]
    i_f = currents[mask]
    
    if len(v_f) < 3:
        return np.array([]), np.array([])
        
    ln_i = np.log(i_f)
    dv = np.diff(v_f)
    dlni = np.diff(ln_i)
    
    # Avoid division by zero
    valid = (dv > 0) & (np.abs(dlni) > 1e-15)
    dv = dv[valid]
    dlni = dlni[valid]
    v_mid = (v_f[:-1][valid] + v_f[1:][valid]) / 2.0
    
    thermal_volt = kb * temperature / q
    n = (1.0 / thermal_volt) * (dv / (dlni + 1e-20))
    
    # Clip n to reasonable values to avoid noise artifacts
    n = np.clip(n, 0.5, 5.0)
    
    return v_mid, n


def calculate_current_from_simulation(sim_results, device_area=1e-4, temperature=300):
    """
    Calculate current using Scharfetter-Gummel discretization.
    
    Parameters:
    -----------
    sim_results : dict
        Dictionary containing simulation results
    device_area : float
        Device cross-sectional area (cm²)
    temperature : float
        Temperature (K)
    """
    from aeslibs.func_lib import Ubern
    
    Vt = kb * temperature / q  # Thermal voltage
    
    voltages = sim_results['voltages']
    n_profiles = sim_results['n_data']
    p_profiles = sim_results['p_data']
    E_profiles = sim_results['electric_field'] # V/cm
    x = sim_results['position'] # m
    x_cm = x * 100.0
    
    # We use carrier concentration as a proxy for doping in quasi-neutral regions
    # for mobility calculation, or better, we can assume it's intrinsic if we don't have dop Profile.
    # But for this PN junction example, we know doping is 1e18.
    # Let's use the actual carrier concentration to get a local effective mobility.
    
    currents = []
    
    for i, V in enumerate(voltages):
        n = n_profiles[i]  # cm^-3
        p = p_profiles[i]  # cm^-3
        E = E_profiles[i]  # V/cm
        
        # Local mobility (cm^2/Vs)
        mu_n = get_silicon_mobility(n, 'n')
        mu_p = get_silicon_mobility(p, 'p')
        
        # Scharfetter-Gummel for Jn and Jp
        # J_i+1/2 = (q*mu*Vt/dx) * [n_i+1 * B(dV/Vt) - n_i * B(-dV/Vt)]
        # dV = V_i+1 - V_i = -E_avg * dx
        
        dx = np.diff(x_cm)
        E_mid = (E[:-1] + E[1:]) / 2.0
        dV_norm = -E_mid * dx / Vt
        
        Jn_mid = []
        Jp_mid = []
        
        for j in range(len(dx)):
            bp, bn = Ubern(dV_norm[j])
            
            # Midpoint mobility for better accuracy
            mu_n_mid = (mu_n[j] + mu_n[j+1]) / 2.0
            mu_p_mid = (mu_p[j] + mu_p[j+1]) / 2.0
            
            # Electron current
            jn = (q * mu_n_mid * Vt / dx[j]) * (n[j+1] * bp - n[j] * bn)
            # Hole current
            jp = (q * mu_p_mid * Vt / dx[j]) * (p[j] * bp - p[j+1] * bn)
            
            Jn_mid.append(jn)
            Jp_mid.append(jp)
            
        J_total = np.array(Jn_mid) + np.array(Jp_mid)
        
        # Use median instead of mean to be robust against numerical noise in junction
        J_avg = np.median(J_total)
        
        I = J_avg * device_area
        currents.append(I)
        
        if i % 10 == 0:
            print(f"V={V:.2f}V, I={I:.4e} A, J_avg={J_avg:.4e} A/cm2")
            
    currents = np.array(currents)
    
    # Offset correction (equilibrium noise)
    if len(currents) > 0:
        offset = currents[0]
        currents = currents - offset
        
    return voltages, currents


def calculate_capacitance_from_simulation(sim_results, device_area=1e-4):
    """
    Calculate capacitance from simulation charge distribution.
    
    C = dQ/dV where Q is the total charge in depletion region
    
    Parameters:
    -----------
    sim_results : dict
        Dictionary containing simulation results
    device_area : float
        Device cross-sectional area (cm²)
    
    Returns:
    --------
    voltages : np.ndarray
        Applied voltages (V)
    capacitances : np.ndarray
        Capacitance values (F)
    """
    voltages = sim_results['voltages']
    
    # This is a simplified implementation
    # For proper implementation, compute depletion width from charge profile
    # and calculate C = epsilon * A / W
    
    # Placeholder - return zeros for now
    # Real implementation would analyze the charge distribution
    capacitances = np.zeros_like(voltages)
    
    return voltages, capacitances


def compute_error_metrics(experimental, simulated):
    """
    Compute error metrics between experimental and simulated data.
    
    Parameters:
    -----------
    experimental : np.ndarray
        Experimental data
    simulated : np.ndarray
        Simulated data (interpolated to match experimental points)
    
    Returns:
    --------
    metrics : dict
        Dictionary containing:
        - 'rmse': Root mean square error
        - 'mae': Mean absolute error
        - 'mape': Mean absolute percentage error
        - 'r2': R² coefficient of determination
    """
    # Ensure arrays are same length
    if len(experimental) != len(simulated):
        raise ValueError("Experimental and simulated data must have same length")
    
    # Remove any NaN or Inf values
    mask = np.isfinite(experimental) & np.isfinite(simulated)
    exp = experimental[mask]
    sim = simulated[mask]
    
    if len(exp) == 0:
        return {'rmse': np.nan, 'mae': np.nan, 'mape': np.nan, 'r2': np.nan}
    
    # RMSE: Root Mean Square Error
    rmse = np.sqrt(np.mean((exp - sim)**2))
    
    # MAE: Mean Absolute Error
    mae = np.mean(np.abs(exp - sim))
    
    # MAPE: Mean Absolute Percentage Error
    # Filter out very small experimental values to avoid artifacts at 0V
    valid_mask = np.abs(exp) > 1e-12
    if np.any(valid_mask):
        mape = np.mean(np.abs((exp[valid_mask] - sim[valid_mask]) / exp[valid_mask])) * 100
    else:
        mape = 0.0
    
    # R²: Coefficient of Determination
    ss_res = np.sum((exp - sim)**2)
    ss_tot = np.sum((exp - np.mean(exp))**2)
    epsilon = 1e-20
    r2 = 1 - (ss_res / (ss_tot + epsilon))
    
    # Log-RMSE: Better for data covering many orders of magnitude
    # Filter for positive values
    pos_mask = (exp > 1e-15) & (sim > 1e-15)
    if np.any(pos_mask):
        log_rmse = np.sqrt(np.mean((np.log10(exp[pos_mask]) - np.log10(sim[pos_mask]))**2))
    else:
        log_rmse = np.nan

    return {
        'rmse': rmse,
        'mae': mae,
        'mape': mape,
        'r2': r2,
        'log_rmse': log_rmse
    }


def plot_iv_comparison(exp_voltage, exp_current, sim_voltage, sim_current, 
                       output_path=None, show=True):
    """
    Generate I-V comparison plot between experimental and simulation data.
    
    Parameters:
    -----------
    exp_voltage : np.ndarray
        Experimental voltage values (V)
    exp_current : np.ndarray
        Experimental current values (A)
    sim_voltage : np.ndarray
        Simulated voltage values (V)
    sim_current : np.ndarray
        Simulated current values (A)
    output_path : str, optional
        Path to save the plot
    show : bool
        Whether to display the plot
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib not available. Skipping plot generation.")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Linear scale plot (Current in mA)
    ax1.plot(exp_voltage, exp_current * 1e3, 'o', label='Experimental', 
             markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2)
    ax1.plot(sim_voltage, sim_current * 1e3, '-', label='Simulation (Mode 10)', 
             color='#003366', linewidth=2.0)
    ax1.set_xlabel('Voltage (V)', fontsize=12)
    ax1.set_ylabel('Current (mA)', fontsize=12)
    ax1.set_title('I-V Characteristics (Linear Scale)', fontsize=13, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.4)
    ax1.legend(fontsize=10, loc='best')
    
    # Restrict X-axis to experimental range for clarity
    v_min, v_max = exp_voltage.min(), exp_voltage.max()
    padding = (v_max - v_min) * 0.05
    ax1.set_xlim(v_min - padding, v_max + padding)

    # Y-axis fitting (Linear in mA)
    i_min, i_max = exp_current.min() * 1e3, exp_current.max() * 1e3
    padding_i = (i_max - i_min) * 0.1
    ax1.set_ylim(i_min - padding_i, i_max + padding_i)
    
    # Semi-log plot (for absolute current magnitude in mA)
    exp_i_abs = np.abs(exp_current) * 1e3
    sim_i_abs = np.abs(sim_current) * 1e3
    
    ax2.semilogy(exp_voltage, exp_i_abs, 'o', label='Experimental',
                 markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2)
    ax2.semilogy(sim_voltage, sim_i_abs, '-', label='Simulation (Mode 10)', 
                 color='#003366', linewidth=2.0)
    
    ax2.set_xlabel('Voltage (V)', fontsize=12)
    ax2.set_ylabel('|Current| (mA, log scale)', fontsize=12)
    ax2.set_title('I-V Characteristics (Semi-log)', fontsize=13, fontweight='bold')
    ax2.grid(True, which='both', linestyle='--', alpha=0.4)
    ax2.legend(fontsize=10, loc='best')
    ax2.set_xlim(v_min - padding, v_max + padding)
    
    # Y-axis fitting (Log in mA)
    pos_exp = exp_i_abs[exp_i_abs > 1e-12]
    if len(pos_exp) > 0:
        y_min_log = max(pos_exp.min() / 3, 1e-10)
        y_max_log = pos_exp.max() * 3
        ax2.set_ylim(y_min_log, y_max_log)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"I-V comparison plot saved to: {output_path}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_cv_comparison(exp_voltage, exp_capacitance, sim_voltage, sim_capacitance,
                       output_path=None, show=True):
    """
    Generate C-V comparison plot between experimental and simulation data.
    
    Parameters:
    -----------
    exp_voltage : np.ndarray
        Experimental voltage values (V)
    exp_capacitance : np.ndarray
        Experimental capacitance values (F)
    sim_voltage : np.ndarray
        Simulated voltage values (V)
    sim_capacitance : np.ndarray
        Simulated capacitance values (F)
    output_path : str, optional
        Path to save the plot
    show : bool
        Whether to display the plot
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib not available. Skipping plot generation.")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # C-V plot
    ax1.plot(exp_voltage, exp_capacitance * 1e12, 'o', label='Experimental',
             markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2)
    ax1.plot(sim_voltage, sim_capacitance * 1e12, '-', label='Simulation',
             color='#003366', linewidth=2.0)
    ax1.set_xlabel('Voltage (V)', fontsize=12)
    ax1.set_ylabel('Capacitance (pF)', fontsize=12)
    ax1.set_title('C-V Characteristics', fontsize=13, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.4)
    ax1.legend(fontsize=10, loc='best')
    
    # 1/C² plot (Mott-Schottky)
    C_sq_inv_exp = 1 / (exp_capacitance**2)
    C_sq_inv_sim = 1 / (sim_capacitance**2 + 1e-30)  # avoid division by zero
    
    ax2.plot(exp_voltage, C_sq_inv_exp * 1e-20, 'o', label='Experimental',
             markersize=5, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.2)
    ax2.plot(sim_voltage, C_sq_inv_sim * 1e-20, '-', label='Simulation',
             color='#003366', linewidth=2.0)
    ax2.set_xlabel('Voltage (V)', fontsize=12)
    ax2.set_ylabel('1/C² (×10²⁰ F⁻²)', fontsize=12)
    ax2.set_title('Mott-Schottky Plot', fontsize=13, fontweight='bold')
    ax2.grid(True, which='both', linestyle='--', alpha=0.4)
    ax2.legend(fontsize=10, loc='best')
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"C-V comparison plot saved to: {output_path}")
    
    if show:
        plt.show()
    else:
        plt.close()


def generate_validation_report(metrics, output_path=None):
    """
    Generate a validation report with error metrics.
    
    Parameters:
    -----------
    metrics : dict
        Dictionary of error metrics from compute_error_metrics()
    output_path : str, optional
        Path to save the report text file
    
    Returns:
    --------
    report : str
        Formatted validation report
    """
    report = "="*60 + "\n"
    report += "REFERENCE COMPARISON REPORT\n"
    report += "Numerical agreement only; source provenance and independent validation are not established.\n"
    report += "="*60 + "\n\n"
    
    report += "Error Metrics:\n"
    report += "-" * 40 + "\n"
    report += f"RMSE (Root Mean Square Error):  {metrics['rmse']:.6e}\n"
    report += f"MAE (Mean Absolute Error):      {metrics['mae']:.6e}\n"
    report += f"MAPE (Mean Absolute % Error):   {metrics['mape']:.2f}%\n"
    report += f"Log-RMSE (Logarithmic Error):   {metrics.get('log_rmse', 0.0):.4f} decades\n"
    report += f"R² (Coefficient of Determination): {metrics['r2']:.4f}\n"
    report += "\n"
    
    # Interpretation
    report += "Interpretation:\n"
    report += "-" * 40 + "\n"
    if metrics['r2'] > 0.95:
        report += "[EXCELLENT] Excellent agreement (R2 > 0.95)\n"
    elif metrics['r2'] > 0.90:
        report += "[GOOD] Good agreement (R2 > 0.90)\n"
    elif metrics['r2'] > 0.80:
        report += "[OK] Acceptable agreement (R2 > 0.80)\n"
    else:
        report += "[POOR] Poor agreement (R2 < 0.80)\n"
    
    if metrics['mape'] < 10:
        report += "[EXCELLENT] Low percentage error (MAPE < 10%)\n"
    elif metrics['mape'] < 20:
        report += "[OK] Moderate percentage error (MAPE < 20%)\n"
    else:
        report += "[POOR] High percentage error (MAPE > 20%)\n"
    
    if 'avg_n' in metrics:
        report += f"Extracted Avg Ideality Factor: {metrics['avg_n']:.2f}\n"

    report += "="*60 + "\n"
    
    print(report)
    
    if output_path:
        with open(output_path, 'w') as f:
            f.write(report)
        print(f"Validation report saved to: {output_path}")
    
    return report


def run_validation_workflow(input_obj, output_dir):
    """
    Standard workflow to perform validation after a simulation.
    Used as a core hook in aestimo.py and individual scripts.
    """
    def get_val(obj, key, default=None):
        if isinstance(obj, dict):
            return obj.get(key, default)
        return getattr(obj, key, default)

    # 1. Identify parameters from input_obj
    enable = get_val(input_obj, 'enable_experimental_validation', False)
    if not enable:
        return
        
    print("\n" + "="*60)
    print("EXPERIMENTAL VALIDATION HOOK")
    print("="*60)
    
    exp_file = get_val(input_obj, 'experimental_iv_file', None)
    if not exp_file:
        print("No experimental data file specified for validation. Skipping.")
        return
    if not os.path.exists(exp_file):
        # Fallback search path logic for specified file
        fname = os.path.basename(exp_file)
        paths_to_try = [
            os.path.join("examples", exp_file),
            os.path.join("examples", "experimental_data", fname),
            os.path.join("..", "examples", "experimental_data", fname)
        ]
        found = False
        for p in paths_to_try:
            if os.path.exists(p):
                exp_file = p
                found = True
                break
        if not found:
            print(f"Warning: Experimental file '{exp_file}' not found. Skipping validation.")
            return

    area = get_val(input_obj, 'device_area', 1e-4)
    # Handle both names used in different script versions
    rs = get_val(input_obj, 'Rs_ext', 0.0) 
    rsh = get_val(input_obj, 'Rsh_ext', 1e12)
    rs_mode = get_val(input_obj, 'rs_mode', "External (Fast)")

    # 2. Load Data
    try:
        exp_v, exp_i = load_experimental_data(exp_file)
        sim_v_int, sim_i_int = load_current_from_avcurr(output_dir, device_area_cm2=area)
    except Exception as e:
        print(f"Validation Error: Could not load data files: {e}")
        return
    
    # 3. Apply Parasitics (avoid double counting if Internal)
    if rs_mode == "Internal (Self-Consistent)":
        # Rs is already in sim_v_int
        sim_v, sim_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=0.0, Rsh=rsh)
    else:
        sim_v, sim_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=rs, Rsh=rsh)
    
    # 4. Metrics & Report
    # Use global interp to align with experimental points
    sim_i_interp = np.interp(exp_v, sim_v, sim_i)
    metrics = compute_error_metrics(exp_i, sim_i_interp)
    
    # Optional: Ideality Factor
    temp = get_val(input_obj, 'T_val', get_val(input_obj, 'T', 300.0))
    v_mid, n_sim = calculate_ideality_factor(sim_v, sim_i, temperature=temp)
    if len(n_sim) > 0:
        metrics['avg_n'] = np.mean(n_sim)

    report_path = os.path.join(output_dir, "validation_report.txt")
    generate_validation_report(metrics, output_path=report_path)
    
    # 5. Plot
    plot_path = os.path.join(output_dir, "iv_comparison.png")
    print(f"Generating I-V comparison plot: {plot_path}")
    plot_iv_comparison(exp_v, exp_i, sim_v, sim_i, output_path=plot_path, show=False)
    
    print("\nValidation success! Files saved to output directory.")
    print("="*60 + "\n")


if __name__ == "__main__":
    print("Experimental Validation Module for Aestimo")
    print("This module provides tools for validating simulation results")
    print("against experimental I-V and C-V data.")
