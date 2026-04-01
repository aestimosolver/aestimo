# Solar Characterization and Automation
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import shutil
from pathlib import Path

# Setup paths to import aestimo
root_dir = os.getcwd()
sys.path.append(root_dir)
import aestimo
from aestimo import run_aestimo

def analyze_iv_curve(voltage, current, area_cm2=1.0, pin_mw_cm2=100.0):
    """
    Analyzes an I-V curve to extract solar cell metrics.
    """
    sort_idx = np.argsort(voltage)
    v = voltage[sort_idx]
    i = current[sort_idx]
    
    j = i # i is already in mA/cm² (from aestimo.py fix)

    # Spike Filtering (Robustness against numerical artifacts)
    # Solar cells should have J that is fairly smooth and monotonically increasing 
    # (becoming less negative) in the power quadrant. Massive negative jumps are artifacts.
    j_filtered = j.copy()
    valid_limit = len(j)
    jsc0 = abs(j[0])
    for idx in range(1, len(j)):
        # If current suddenly jumps more than 50% of Jsc in a single step (0.05V), 
        # or becomes unphysically huge (>100x Jsc), it's likely a numerical divergence.
        if v[idx] > 0.05 and (j[idx] < j[idx-1] - 0.5 * max(jsc0, 1e-3) or abs(j[idx]) > 1e4):
            valid_limit = idx
            break
            
    if len(v) < 2:
        # Return empty metrics if simulation diverged immediately
        return {
            'jsc': 0.0, 'voc': 0.0, 'pmpp': 0.0, 'vmpp': 0.0, 'jmpp': 0.0,
            'ff': 0.0, 'eta': 0.0, 'rs': 0.0, 'rsh': 0.0, 'v': v, 'j': j, 'p': np.zeros_like(j)
        }
    
    v = v[:valid_limit]
    j = j[:valid_limit]
    
    # 1. Jsc (at V=0)
    jsc = -np.interp(0.0, v, j) 
    
    # 2. Voc (at J=0)
    if np.min(j) < 0 and np.max(j) > 0:
        voc = np.interp(0.0, j, v)
    else:
        voc = 0.0
        
    # 3. Power Density
    p_density = v * (-j) # mW/cm^2
    
    # Find MPP
    valid_mask = (v >= 0) & (v <= (voc if voc > 0 else np.max(v)))
    if np.any(valid_mask):
        p_valid = p_density[valid_mask]
        v_valid = v[valid_mask]
        j_valid = j[valid_mask]
        mpp_idx = np.argmax(p_valid)
        pmpp = p_valid[mpp_idx]
        vmpp = v_valid[mpp_idx]
        jmpp = -j_valid[mpp_idx]
    else:
        pmpp, vmpp, jmpp = 0.0, 0.0, 0.0

    # 4. Fill Factor
    ff = (pmpp / (voc * jsc) * 100.0) if (voc > 0 and jsc > 0) else 0.0
        
    # 5. Efficiency
    eta = (pmpp / pin_mw_cm2) * 100.0 if pin_mw_cm2 > 0 else 0.0
    
    # 6. Resistances (Requires at least 2 points)
    if len(v) >= 2:
        dv = np.gradient(v)
        dj = np.gradient(j)
        slope = dj / dv # mA/cm^2 / V
    else:
        slope = np.zeros_like(v)
    
    # Rsh at V=0
    idx_sc = np.argmin(np.abs(v - 0.0))
    slope_sc = slope[idx_sc]
    rsh = 1.0 / (slope_sc * 1e-3) if slope_sc > 0 else np.inf
        
    # Rs at Voc
    if voc > 0:
        idx_oc = np.argmin(np.abs(v - voc))
        slope_oc = slope[idx_oc]
        rs = 1.0 / (slope_oc * 1e-3) if slope_oc > 0 else 0.0
    else:
        rs = 0.0
        
    return {
        'jsc': jsc, 'voc': voc, 'pmpp': pmpp, 'vmpp': vmpp, 'jmpp': jmpp,
        'ff': ff, 'eta': eta, 'rs': rs, 'rsh': rsh, 'v': v, 'j': j, 'p': p_density
    }

def run_simulation(config_dict, label):
    """Utility to run Aestimo with a specific configuration and return metrics."""
    print(f"\n>>> Running Simulation: {label}...")
    
    # Force scheme 8 (Gummel) for solar
    config_dict['comp_scheme'] = 8
    config_dict['computation_scheme'] = 8 

    class InputObject:
        def __init__(self, d):
            for k,v in d.items(): setattr(self, k, v)
            self.__file__ = f"sim_{label}.py"
            self.inputfilename = f"sim_{label}"

    # Setup output directory
    out_dir = os.path.abspath(os.path.join(root_dir, "examples", f"solar_study_{label}_output"))
    if os.path.exists(out_dir): shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    
    # IMPORTANT: We must set the GLOBAL output_directory in aestimo
    aestimo.output_directory = out_dir
    print(f"DEBUG: Setting aestimo.output_directory = {out_dir}")
    
    # Build doping profile manually
    N_p = int(250.0 / 5.0)
    N_n = int(250.0 / 5.0)
    dop_arr = np.zeros(N_p + N_n)
    dop_arr[:N_p] = -config_dict['material'][0][4] * 1e6
    dop_arr[N_p:] =  config_dict['material'][1][4] * 1e6
    config_dict['dop_profile'] = dop_arr

    # Run
    import config
    config.Drift_Diffusion_out = True
    config.potential_out = True
    config.sigma_out = True
    
    print(f"DEBUG: vmin={config_dict.get('vmin')}, vmax={config_dict.get('vmax')}, Each_Step={config_dict.get('Each_Step')}")
    
    run_aestimo(InputObject(config_dict), drawFigures=False, show=False)
    
    # run_aestimo might have changed aestimo.output_directory!
    actual_out = aestimo.output_directory
    print(f"DEBUG: Solver finished. Final aestimo.output_directory = {actual_out}")
    
    # Check if file exists
    iv_path = os.path.join(actual_out, 'av_curr.dat')
    if not os.path.exists(iv_path):
        print(f"CRITICAL ERROR: {iv_path} not found after simulation!")
        if os.path.exists(actual_out):
            print(f"Contents of {actual_out}: {os.listdir(actual_out)}")
        else:
            print(f"Directory {actual_out} does not exist!")
    
    # Analyze
    data = np.loadtxt(iv_path)
    metrics = analyze_iv_curve(data[:,0], data[:,1], area_cm2=config_dict['device_area'])
    return metrics

def main():
    # Base configuration for InGaN p-n
    base_config = {
        'material': [
            [250.0, 'InGaN', 0.57, 0.0, 1e16, 'p', 'b'],
            [250.0, 'InGaN', 0.57, 0.0, 2e17, 'n', 'b']
        ],
        'T': 300.0,
        'gridfactor': 5.0,
        'maxgridpoints': 200000,
        'mat_type': 'Wurtzite',
        'comp_scheme': 8,
        'vmin': -0.1,
        'vmax': 1.0, 
        'Each_Step': 0.05,
        'tat_field': 1e10, 
        'G_optical': 0.0,
        'device_area': 1.0e-4 
    }

    # 1. Dark vs Light Comparison
    dark_metrics = run_simulation(base_config.copy(), "Dark")
    
    light_config = base_config.copy()
    light_config['G_optical'] = 1e18
    light_metrics = run_simulation(light_config, "Light_1e18")

    # 2. Temperature Sweep
    temp_results = []
    temps = [250, 300, 350, 400]
    for T in temps:
        t_config = light_config.copy()
        t_config['T'] = float(T)
        m = run_simulation(t_config, f"Temp_{T}K")
        temp_results.append((T, m))

    # --- PLOTTING ---
    plt.figure(figsize=(15, 10))

    # Subplot 1: Dark vs Light J-V
    plt.subplot(2, 2, 1)
    plt.plot(dark_metrics['v'], dark_metrics['j'], 'k--', label='Dark')
    plt.plot(light_metrics['v'], light_metrics['j'], 'r-', label='Light (1e18)')
    plt.axhline(0, color='gray', lw=0.5)
    plt.axvline(0, color='gray', lw=0.5)
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current Density (mA/cm²)')
    plt.title('Dark vs Illuminated J-V')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Subplot 2: Power-Voltage Curve (Generation Quadrant)
    plt.subplot(2, 2, 2)
    plt.plot(light_metrics['v'], light_metrics['p'], 'g-', label='Power (Light)')
    plt.axvline(light_metrics['vmpp'], color='orange', ls=':', label=f'MPP: {light_metrics['pmpp']:.2f} mW/cm²')
    plt.xlabel('Voltage (V)')
    plt.ylabel('Power Density (mW/cm²)')
    plt.title('P-V Characteristic')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Subplot 3: Temperature Dependence of Voc
    plt.subplot(2, 2, 3)
    t_vals = [r[0] for r in temp_results]
    voc_vals = [r[1]['voc'] for r in temp_results]
    plt.plot(t_vals, voc_vals, 'bo-')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Voc (V)')
    plt.title('Voc vs Temperature')
    plt.grid(True, alpha=0.3)

    # Subplot 4: Temperature Dependence of Efficiency
    plt.subplot(2, 2, 4)
    eta_vals = [r[1]['eta'] for r in temp_results]
    plt.plot(t_vals, eta_vals, 'ro-')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Efficiency (%)')
    plt.title('Efficiency vs Temperature')
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('solar_characterization_summary.png')
    print("\nSummary plot saved as 'solar_characterization_summary.png'")

    # Report for Light 1e18
    print("\n" + "="*40)
    print("FINAL CHARACTERIZATION (Light 1e18, 300K)")
    print("="*40)
    for k, v in light_metrics.items():
        if isinstance(v, (float, int, np.float64)):
            print(f"{k+':':<10} {v:>15.4g}")

if __name__ == "__main__":
    main()

