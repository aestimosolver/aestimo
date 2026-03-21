import sys, os, json
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '.')
import aestimo

def main():
    json_path = 'examples/ingan_solar_high_perf.json'
    print(f"Loading {json_path}...")
    with open(json_path, 'r') as f:
        config = json.load(f)
    
    # Mocking input object like in run_user_json.py
    class InputObject: pass
    in_obj = InputObject()
    for k, v in config.items():
        setattr(in_obj, k, v)
    
    # Required properties for aestimo
    in_obj.T_val = float(config["temp"])
    in_obj.mat_type = config["mat_sys"]
    in_obj.gridfactor = 1.0 # 1nm grid for speed
    in_obj.maxgridpoints = 10000
    in_obj.dx = in_obj.gridfactor * 1e-9
    in_obj.vmax = 2.2 # Sweep slightly past Voc avoids overflow
    in_obj.vmin = 0.0
    in_obj.Each_Step = 0.05
    in_obj.surface = np.array([float(config["bc_left"]), float(config["bc_right"])])
    in_obj.subnumber_h = int(config["sub_h"])
    in_obj.subnumber_e = int(config["sub_e"])
    # Convert field from mV/cm to V/m
    in_obj.Fapplied = float(config["field"]) * 1e5
    
    # Avoid overriding work functions to let the solver determine the built-in potential naturally
    in_obj.surface_recomb = (0, 0)
    in_obj.enable_polarization = True
    in_obj.Quantum_Regions = False
    in_obj.Quantum_Regions_boundary = np.zeros((1, 2))
    in_obj.__file__ = os.path.abspath(json_path)
    in_obj.Rs = 0.0
    in_obj.computation_scheme = 9 
    in_obj.tat_field = 1e12 # Disable TAT 
    in_obj.G_optical = 2e21 * 1e6 # Approx 1-sun level for InGaN 20% (convert cm^-3 to m^-3 s^-1)
    in_obj.tau = float(config.get("tau", 2e-6))
    in_obj.photovoltaic_mode = True

    # Map layers and build doping profile (Sync with GUI)
    material_list = []
    tot_m = sum(float(l["thickness"]) for l in config["layers"]) * 1e-9
    n_max = int(tot_m / in_obj.dx)
    dop_arr = np.zeros(n_max)
    curr_idx = 0
    
    for l in config["layers"]:
        th = float(l["thickness"])
        dop_val = float(l["doping"])
        material_list.append([th, l["material"], float(l["mole"]), 0.0, dop_val, l["doping_type"], l["type"][0]])
        val = dop_val
        if l["doping_type"] == 'p': val = -val
        steps = int((th * 1e-9) / in_obj.dx)
        end_idx = min(curr_idx + steps, n_max)
        dop_arr[curr_idx:end_idx] = val * 1e6 # cm-3 to m-3
        curr_idx = end_idx
        
    in_obj.material = material_list
    in_obj.dop_profile = dop_arr
    in_obj.n_max = n_max
    
    aestimo.output_directory = 'final_plots_output'
    if not os.path.exists(aestimo.output_directory):
        os.makedirs(aestimo.output_directory)
    
    print("Running simulation...")
    # Capture results for diagnostics
    _, _, res, _ = aestimo.run_aestimo(in_obj, drawFigures=False, show=False)
    
    # Load and process IV data
    iv_data = np.loadtxt(os.path.join(aestimo.output_directory, 'av_curr.dat'))
    v_raw = iv_data[:, 0]
    j_raw = iv_data[:, 1] * 0.1 # mA/cm2
    
    # Find Voc
    voc = 0.0
    for i in range(len(v_raw) - 1):
        if j_raw[i] < 0 and j_raw[i+1] >= 0:
            voc = v_raw[i] + (0 - j_raw[i]) * (v_raw[i+1] - v_raw[i]) / (j_raw[i+1] - j_raw[i])
    
    # Filter spikes
    jsc = abs(j_raw[0])
    limit_idx = len(v_raw)
    for i in range(1, len(v_raw)):
        if v_raw[i] > 0.5 and (j_raw[i] < j_raw[i-1] - 0.5): 
             limit_idx = i
             break
    
    v_iv = v_raw[:limit_idx]
    j_iv = j_raw[:limit_idx]
    p_iv = v_iv * (-j_iv)
    
    pv_mask = (v_iv >= 0) & (v_iv <= (voc + 0.1 if voc > 0 else 2.3))
    if np.any(pv_mask):
        p_max = np.max(p_iv[pv_mask])
        v_mpp = v_iv[pv_mask][np.argmax(p_iv[pv_mask])]
    else:
        p_max, v_mpp = 0, 0
    
    ff = (p_max / (voc * jsc) * 100.0) if (voc > 0 and jsc > 0) else 0.0
    eff = (p_max / 100.0) * 100.0
    
    # --- DIAGNOSTIC: BAND DIAGRAM AT V_MAX ---
    plt.figure(figsize=(10, 6))
    x_um = res.xaxis * 1e6
    plt.plot(x_um, res.Ec_result, 'k-', lw=2, label='Ec')
    plt.plot(x_um, res.Ev_result, 'k-', lw=2, label='Ev')
    plt.plot(x_um, res.Efn_result, 'r--', label='Efn')
    plt.plot(x_um, res.Efp_result, 'b--', label='Efp')
    plt.xlabel('Position (µm)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Diagnostic Band Diagram at V = {v_raw[res.Total_Steps-1]:.2f} V')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('band_diagnostic.png')
    plt.close()

    # --- FINAL SUMMARY PLOT ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
    ax1.plot(v_iv, j_iv, color='#1f77b4', lw=4, label='Simulation')
    ax1.axhline(0, color='k', lw=1.5, ls='--')
    ax1.axvline(0, color='k', lw=1.5, ls='--')
    ax1.set_xlabel('Voltage (V)', fontsize=13)
    ax1.set_ylabel('Current (mA/cm²)', fontsize=13)
    ax1.set_title('InGaN Solar Cell: J-V Characteristic', fontsize=15, fontweight='bold')
    ax1.set_xlim([-0.1, 2.3])
    ax1.set_ylim([-jsc*1.4, jsc*0.4])
    ax1.grid(True, linestyle=':', alpha=0.7)
    ax1.legend()
    
    ax2.plot(v_iv[v_iv >= 0], p_iv[v_iv >= 0], color='#d62728', lw=4)
    ax2.set_xlabel('Voltage (V)', fontsize=13)
    ax2.set_ylabel('Power (mW/cm²)', fontsize=13)
    ax2.set_title('Power-Voltage Curve', fontsize=15, fontweight='bold')
    ax2.set_xlim([0, 2.3])
    ax2.set_ylim([0, p_max*1.5 if p_max > 0 else 1.0])
    ax2.grid(True, linestyle=':', alpha=0.7)
    
    if p_max > 0:
        ax2.scatter(v_mpp, p_max, color='black', s=120, zorder=5)
        ax2.annotate(f'MPP: {p_max:.3f} mW/cm²', xy=(v_mpp, p_max), xytext=(v_mpp+0.1, p_max*1.2),
                    arrowprops=dict(facecolor='black', shrink=0.08, width=1.5), fontsize=12, fontweight='bold')

    summary_box_text = (
        f"SOLAR CELL PERFORMANCE SUMMARY\n"
        f"------------------------------\n"
        f"Jsc: {jsc:.3f} mA/cm²\n"
        f"Voc: {voc:.3f} V\n"
        f"Pmax: {p_max:.3f} mW/cm²\n"
        f"FF: {ff:.2f}%\n"
        f"Eff: {eff:.3f}%"
    )
    plt.figtext(0.5, 0.05, summary_box_text, ha="center", fontsize=14, family='monospace',
                bbox={"facecolor":"#f9f9f9", "edgecolor":"#1f77b4", "boxstyle":"round,pad=1.2", "alpha":0.95})
    
    plt.tight_layout(rect=[0, 0.2, 1, 1])
    plt.savefig('final_results_summary.png', dpi=150)
    plt.show()

    print(f"Results: Jsc={jsc:.4f} mA/cm2, Voc={voc:.4f} V, Pmax={p_max:.4f} mW/cm2, FF={ff:.2f}%, Eff={eff:.4f}%")

if __name__ == "__main__":
    main()
