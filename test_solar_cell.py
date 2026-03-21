#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Test script for Solar Cell photogeneration
import numpy as np
import os
import sys

# Setup paths
root_dir = os.getcwd()
sys.path.append(root_dir)

from aestimo import run_aestimo

from characterize_solar import analyze_iv_curve
import matplotlib.pyplot as plt

def test_solar_cell():
    print("Running Improved Solar Cell Verification...")
    
    # Create a simple class to mimic the object structure anticipated by Aestimo
    class InputObject:
        def __init__(self, dictionary):
            for key, value in dictionary.items():
                setattr(self, key, value)
            self.__file__ = "test_solar_cell.py"
            self.inputfilename = "test_solar_cell"
                
    # Base InGaN Structure
    config_dict = {
        'material': [
            [250.0, 'InGaN', 0.57, 0.0, 1e16, 'p', 'b'],
            [250.0, 'InGaN', 0.57, 0.0, 2e17, 'n', 'b']
        ],
        'T': 300.0,
        'gridfactor': 5.0, # Improved stability for DD solver
        'maxgridpoints': 200000,
        'mat_type': 'Wurtzite',
        'computation_scheme': 9, # Coupled Newton (more stable for high injection)
        'vmin': 0.0,
        'vmax': -1.8, # Negative bias to lower barrier (Extended to capture full Voc)
        'Each_Step': -0.1,
        'fval': 0.0, 
        'tat_field': 1e12, # Disabled for logical baseline
        'enable_polarization': False, # CRITICAL: Disable polarization to prevent p-type inversion in InGaN
        'G_optical': 5e18, # cm^-3 s^-1 (→ 5e24 m^-3 s^-1, calibrated for Jsc ~20-30 mA/cm²)
        'device_area': 1.0e-4 
    }
    
    # Setup Doping Profile
    gf = config_dict['gridfactor']
    N_p = int(250.0 / gf)
    N_n = int(250.0 / gf)
    dop_p = -1e16 * 1e6
    dop_n =  2e17 * 1e6
    dop_arr = np.zeros(N_p + N_n)
    dop_arr[:N_p] = dop_p
    dop_arr[N_p:] = dop_n
    config_dict['dop_profile'] = dop_arr
    
    # Setup output directory
    out_dir = os.path.abspath(os.path.join(os.getcwd(), "examples", "test_solar_cell_output"))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    import aestimo
    aestimo.output_directory = out_dir

    # 1. Run Simulation
    print("Executing simulation...")
    config_obj = InputObject(config_dict)
    _, model, result_ps, figures = run_aestimo(config_obj, drawFigures=False, show=False)
    # run_aestimo with scheme 8 already calls Poisson_Schrodinger_DD_test and saves files
    
    # 2. Extract and Analyze Data
    iv_path = os.path.join(out_dir, 'av_curr.dat')
    if not os.path.exists(iv_path):
        print(f"[ERROR] {iv_path} not found!")
        return

    data = np.loadtxt(iv_path)
    voltages = data[:, 0]
    # av_curr.dat contains current density in A/m^2
    j_density_am2 = data[:, 1]
    
    # Direct conversion
    # av_curr.dat is in A/m²
    j_raw = data[:, 1]
    v_raw = data[:, 0]
    
    # 1. Handle Voltage Sign (Forward Bias was Negative Voltage in simulation)
    # Physical Forward Voltage = -v_raw 
    # (Assuming simulation sweep was 0 -> -1.5)
    
    # Sort by voltage (increasing)
    v_phys = -v_raw
    sort_idx = np.argsort(v_phys)
    v = v_phys[sort_idx]
    j_density_am2 = j_raw[sort_idx]
    
    # Convert to mA/cm²
    # 1 A/m² = 0.1 mA/cm²
    j = j_density_am2 / 10.0
    
    # 2. Calculate Metrics
    # Jsc (current at V=0). Note: Current flows N->P (positive in code convention)
    # Standard solar cell convention: J is negative (photocurrent).
    # But here J is positive.
    # We can stick to physical values or invert for standard plotting.
    # Let's keep J as is, but P = V * J.
    # If J > 0 and V > 0, P > 0 (Generating? No, Consuming).
    # Solar Cell Generates power. V > 0, I < 0 (standard).
    # Or V > 0, I > 0 (load convention).
    # Here I > 0 at V=0.
    # At V=Voc, I=0.
    # Power P = V * I. 
    # Efficiency is Max(V*I).
    
    jsc = np.interp(0.0, v, j) # Should be ~ +11.7
    
    # Voc (where J crosses 0)
    # Search for crossing
    voc = 0.0
    for k in range(len(v)-1):
        if j[k] * j[k+1] <= 0:
            # Linear interpolation
            slope = (v[k+1] - v[k]) / (j[k+1] - j[k])
            voc = v[k] - j[k] * slope
            break
            
    # If no crossing found, maybe huge leakage or breakdown?
    if voc == 0.0 and np.min(j) < 0:
         # It crossed somewhere
         voc = np.interp(0.0, j, v) # only works if monotonic decreasing
    elif voc == 0.0 and len(v) > 0 and j[-1] > 0:
         # Never crossed 0?
         voc = v[-1] # Approximation
    
    # Power
    p_curve = v * j # mW/cm²
    mpp_idx = np.argmax(p_curve)
    pmpp = p_curve[mpp_idx]
    vmpp = v[mpp_idx]
    jmpp = j[mpp_idx]
    
    # Theoretical inputs
    pin = config_dict.get('device_area', 1e-4) * 100.0 # mW ? No.
    # Pin = 100 mW/cm^2. Efficiency = Pmpp / 100.
    eta = pmpp # if Pmpp is in mW/cm^2 and Pin=100
    
    # Fill Factor
    ff = (pmpp / (jsc * voc)) * 100.0 if (jsc*voc) != 0 else 0.0
    
    print("\n----------------------------------------")
    print("SOLAR CELL VERIFICATION REPORT")
    print("----------------------------------------")
    print(f"Jsc:  {jsc:.4f} mA/cm² (Positive = N->P flow)")
    print(f"Voc:  {voc:.4f} V")
    print(f"Pmax: {pmpp:.4f} mW/cm²")
    print(f"FF:   {ff:.2f} %")
    print(f"Efficiency (eta): {eta:.2f} %")
    print("----------------------------------------")

    # 3. Visualization
    plt.figure(figsize=(10, 6))
    ax1 = plt.gca()
    
    # Plot J-V
    ax1.plot(v, j, 'b-o', label='J-V Curve')
    ax1.set_xlabel('Voltage (V)')
    ax1.set_ylabel('Current Density (mA/cm²)')
    ax1.set_ylim(bottom=min(-5, np.min(j)*1.1), top=max(15, np.max(j)*1.1))
    ax1.grid(True)
    ax1.legend(loc='upper right')
    
    # Plot Power on twin axis
    ax2 = ax1.twinx()
    ax2.plot(v, p_curve, 'r--', label='Power (mW/cm²)')
    ax2.set_ylabel('Power Density (mW/cm²)')
    ax2.legend(loc='lower left')
    
    plt.title(f"InGaN Solar Cell IV Characteristics\n(Polarization Disabled, V_scan Inverted)")
    output_plot = os.path.join(out_dir, 'iv_solar_test_final.png')
    plt.savefig(output_plot)
    print(f"Final plot saved to: {output_plot}")
    ax2.plot(v, p_curve, 'g--', label='P-V Curve')
    
    # Mark MPP
    ax1.plot(vmpp, jmpp, 'ro', label='MPP')
    
    ax1.set_xlabel('Voltage (V)')
    ax1.set_ylabel('Current Density (mA/cm²)', color='b')
    ax2.set_ylabel('Power Density (mW/cm²)', color='g')
    plt.title(f"InGaN Solar Cell IV Characteristics\nFF={ff:.1f}%, Efficiency={eta:.2f}%")
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0, color='black', lw=1)
    ax1.axvline(0, color='black', lw=1)
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    plot_path = os.path.join(out_dir, 'iv_solar_test.png')
    plt.savefig(plot_path)
    print(f"Plot saved to: {plot_path}")
    
    if voc > 0.1 and ff > 10:
        print("\n[SUCCESS] PV characterization complete and verified.")
    else:
        print("\n[WARNING] Poor PV performance detected. Check parameters.")

if __name__ == "__main__":
    test_solar_cell()
