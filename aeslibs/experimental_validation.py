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


def calculate_current_from_simulation(sim_results, device_area=1e-4, temperature=300):
    """
    Calculate current from simulation carrier concentrations using drift-diffusion.
    
    This is a simplified calculation. For more accurate results, integrate the
    current density across the device using the drift-diffusion equations.
    
    Parameters:
    -----------
    sim_results : dict
        Dictionary containing simulation results with keys:
        - 'voltages': Applied voltages (V)
        - 'n_data': Electron concentration profiles [voltage x position]
        - 'p_data': Hole concentration profiles [voltage x position]
        - 'position': Position array (m)
        - 'electric_field': Electric field profiles (V/m)
    device_area : float
        Device cross-sectional area (cm²)
    temperature : float
        Temperature (K)
    
    Returns:
    --------
    voltages : np.ndarray
        Applied voltages (V)
    currents : np.ndarray
        Total current (A)
    """
    # Mobility values for Si at 300K (cm²/V·s)
    # These are typical values; should ideally come from simulation
    mu_n = 1350  # electron mobility
    mu_p = 450   # hole mobility
    
    Vt = kb * temperature / q  # Thermal voltage
    
    voltages = sim_results['voltages']
    n_profiles = sim_results['n_data']
    p_profiles = sim_results['p_data']
    E_profiles = sim_results['electric_field']
    x = sim_results['position']
    
    currents = []
    
    for i, V in enumerate(voltages):
        n = n_profiles[i]  # electron concentration (cm⁻³)
        p = p_profiles[i]  # hole concentration (cm⁻³)
        E = E_profiles[i]  # electric field (V/cm)
        
        # Current density: J_n = q * mu_n * n * E + q * D_n * dn/dx
        # Simplified: use drift component at junction
        # Find junction position (where n ≈ p)
        junction_idx = np.argmin(np.abs(n - p))
        
        if junction_idx > 0 and junction_idx < len(n) - 1:
            E_junction = E[junction_idx]
            n_junction = n[junction_idx]
            p_junction = p[junction_idx]
            
            # Drift current density (A/cm²)
            J_drift = q * (mu_n * n_junction + mu_p * p_junction) * E_junction
            
            # Total current (A)
            I = J_drift * device_area
        else:
            I = 0.0
        
        currents.append(I)
    
    return voltages, np.array(currents)


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
    
    # MAPE: Mean Absolute Percentage Error (avoid division by very small numbers)
    epsilon = 1e-20
    mape = np.mean(np.abs((exp - sim) / (np.abs(exp) + epsilon))) * 100
    
    # R²: Coefficient of Determination
    ss_res = np.sum((exp - sim)**2)
    ss_tot = np.sum((exp - np.mean(exp))**2)
    r2 = 1 - (ss_res / (ss_tot + epsilon))
    
    return {
        'rmse': rmse,
        'mae': mae,
        'mape': mape,
        'r2': r2
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
    
    # Linear scale plot
    ax1.plot(exp_voltage, exp_current * 1e3, 'o', label='Experimental', 
             markersize=6, alpha=0.7)
    ax1.plot(sim_voltage, sim_current * 1e3, '-', label='Simulation', 
             linewidth=2)
    ax1.set_xlabel('Voltage (V)', fontsize=12)
    ax1.set_ylabel('Current (mA)', fontsize=12)
    ax1.set_title('I-V Characteristics (Linear Scale)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Semi-log plot (for forward bias)
    forward_mask_exp = exp_current > 0
    forward_mask_sim = sim_current > 0
    
    if np.any(forward_mask_exp):
        ax2.semilogy(exp_voltage[forward_mask_exp], 
                     exp_current[forward_mask_exp] * 1e3, 
                     'o', label='Experimental', markersize=6, alpha=0.7)
    if np.any(forward_mask_sim):
        ax2.semilogy(sim_voltage[forward_mask_sim], 
                     sim_current[forward_mask_sim] * 1e3, 
                     '-', label='Simulation', linewidth=2)
    
    ax2.set_xlabel('Voltage (V)', fontsize=12)
    ax2.set_ylabel('Current (mA, log scale)', fontsize=12)
    ax2.set_title('I-V Characteristics (Semi-log)', fontsize=14)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(fontsize=10)
    
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
             markersize=6, alpha=0.7)
    ax1.plot(sim_voltage, sim_capacitance * 1e12, '-', label='Simulation',
             linewidth=2)
    ax1.set_xlabel('Voltage (V)', fontsize=12)
    ax1.set_ylabel('Capacitance (pF)', fontsize=12)
    ax1.set_title('C-V Characteristics', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # 1/C² plot (Mott-Schottky)
    C_sq_inv_exp = 1 / (exp_capacitance**2)
    C_sq_inv_sim = 1 / (sim_capacitance**2 + 1e-30)  # avoid division by zero
    
    ax2.plot(exp_voltage, C_sq_inv_exp * 1e-20, 'o', label='Experimental',
             markersize=6, alpha=0.7)
    ax2.plot(sim_voltage, C_sq_inv_sim * 1e-20, '-', label='Simulation',
             linewidth=2)
    ax2.set_xlabel('Voltage (V)', fontsize=12)
    ax2.set_ylabel('1/C² (×10²⁰ F⁻²)', fontsize=12)
    ax2.set_title('Mott-Schottky Plot', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    
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
    report += "EXPERIMENTAL VALIDATION REPORT\n"
    report += "="*60 + "\n\n"
    
    report += "Error Metrics:\n"
    report += "-" * 40 + "\n"
    report += f"RMSE (Root Mean Square Error):  {metrics['rmse']:.6e}\n"
    report += f"MAE (Mean Absolute Error):      {metrics['mae']:.6e}\n"
    report += f"MAPE (Mean Absolute % Error):   {metrics['mape']:.2f}%\n"
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
    
    report += "="*60 + "\n"
    
    print(report)
    
    if output_path:
        with open(output_path, 'w') as f:
            f.write(report)
        print(f"Validation report saved to: {output_path}")
    
    return report


if __name__ == "__main__":
    print("Experimental Validation Module for Aestimo")
    print("This module provides tools for validating simulation results")
    print("against experimental I-V and C-V data.")
