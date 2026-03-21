import sys, os, json
import numpy as np
sys.path.insert(0, '.')
import aestimo

def run_verified_json(json_path):
    with open(json_path, 'r') as f:
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
    in_obj.tat_field = 1e12
    in_obj.work_function_left = 4.2
    in_obj.work_function_right = 5.2
    in_obj.surface_recomb = (1e3, 1e3)
    in_obj.enable_polarization = True
    in_obj.Quantum_Regions = False
    in_obj.Quantum_Regions_boundary = np.zeros((1, 2))
    in_obj.__file__ = os.path.abspath(json_path)
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

    aestimo.output_directory = 'sweep_output/user_verification'
    if not os.path.exists(aestimo.output_directory):
        os.makedirs(aestimo.output_directory)
    
    print(f"Running simulation for {config['name']}...")
    aestimo.run_aestimo(in_obj, drawFigures=False, show=False)
    
    # Load and check results
    res_path = os.path.join(aestimo.output_directory, 'av_curr.dat')
    if os.path.exists(res_path):
        data = np.loadtxt(res_path)
        jsc = np.abs(data[0, 1]) / 10.0 # Convert A/m2 to mA/cm2
        # Voc search
        voc = 0.0
        for i in range(len(data)-1):
            if data[i, 1] < 0 and data[i+1, 1] > 0:
                # Linear interpolation for better accuracy
                v1, j1 = data[i]
                v2, j2 = data[i+1]
                voc = v1 - j1 * (v2 - v1) / (j2 - j1)
                break
        print(f"RESULTS:")
        print(f"Jsc: {jsc:.4f} mA/cm2")
        print(f"Voc: {voc:.4f} V")
    else:
        print("Error: av_curr.dat not found.")

if __name__ == "__main__":
    run_verified_json('examples/ingan_solar_high_perf.json')
