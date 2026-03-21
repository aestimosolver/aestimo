import sys, os, json
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '.')
import aestimo

def run_sim(tat_field_val):
    with open('examples/ingan_solar_high_perf.json', 'r') as f:
        config = json.load(f)
    
    # Create mock input object
    class InputObject: pass
    in_obj = InputObject()
    
    in_obj.T = float(config["temp"])
    in_obj.T_val = in_obj.T
    in_obj.mat_type = config["mat_sys"]
    in_obj.gridfactor = float(config["grid_step"])
    in_obj.maxgridpoints = int(config["max_pts"])
    in_obj.dx = in_obj.gridfactor * 1e-9
    
    in_obj.vmax = float(config["vmax"])
    in_obj.vmin = float(config["vmin"])
    in_obj.Each_Step = float(config["vstep"])
    in_obj.surface = np.array([float(config["bc_left"]), float(config["bc_right"])])
    in_obj.photovoltaic_mode = config.get("photovoltaic_mode", False)
    in_obj.G_optical = float(config.get("G_optical", 0.0))
    in_obj.tau = float(config.get("tau", 1e-6))
    
    in_obj.subnumber_h = int(config["sub_h"])
    in_obj.subnumber_e = int(config["sub_e"])
    in_obj.Fapplied = float(config["field"]) * 1e5
    in_obj.tat_field = tat_field_val # V/m
    in_obj.work_function_left = 4.2
    in_obj.work_function_right = 5.2
    in_obj.surface_recomb = (1e3, 1e3)
    in_obj.enable_polarization = True
    in_obj.Quantum_Regions = False
    in_obj.Quantum_Regions_boundary = np.zeros((1, 2))
    in_obj.__file__ = os.path.abspath('sweep_tat.py')
    in_obj.Rs = 0.0
    in_obj.computation_scheme = 7 # Poisson_Schrodinger_DD

    # Map layers
    material_list = []
    for l in config["layers"]:
        material_list.append([
            float(l["thickness"]), 
            l["material"], 
            float(l["mole"]), 
            0.0, 
            float(l["doping"]), 
            l["doping_type"], 
            l["type"][0]
        ])
    in_obj.material = material_list

    out_subdir = f"tat_{tat_field_val:.0e}"
    aestimo.output_directory = os.path.join('sweep_output/tat_sweep', out_subdir)
    if not os.path.exists(aestimo.output_directory):
        os.makedirs(aestimo.output_directory)
    
    aestimo.run_aestimo(in_obj, drawFigures=False, show=False)
    
    res_path = os.path.join(aestimo.output_directory, 'av_curr.dat')
    data = np.loadtxt(res_path)
    return data # V, J_m2

def main():
    tat_vals = [1e12, 5e7, 1e7, 1e6] # 1e12 is effectively disabled
    labels = ["Disabled (1e12)", "TAT (5e7)", "TAT (1e7)", "TAT (1e6)"]
    
    plt.figure(figsize=(10, 6))
    
    for val, label in zip(tat_vals, labels):
        print(f"Running simulation with tat_field = {val:.2e}...")
        data = run_sim(val)
        v = data[:, 0]
        j = data[:, 1] / 10.0 # mA/cm2
        plt.plot(v, j, label=label)
        
        # Calculate Jsc
        jsc = abs(j[0])
        print(f"  -> Jsc: {jsc:.4f} mA/cm2")

    plt.axhline(0, color='black', lw=1)
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current Density (mA/cm²)')
    plt.title('Impact of Trap-Assisted Tunneling (TAT) on InGaN Solar Cell')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim([-0.6, 0.2])
    plt.xlim([0, 2.4])
    plt.savefig('tat_impact_iv.png', dpi=150)
    print("TAT impact plot saved to tat_impact_iv.png")

if __name__ == "__main__":
    main()
