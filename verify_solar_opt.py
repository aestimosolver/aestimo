import sys, os
sys.path.append(os.getcwd())
import aestimo
import numpy as np
import matplotlib.pyplot as plt
import json

def verify_solar():
    # Load the JSON config for InGaN solar cell
    json_path = os.path.join("examples", "ingan_solar_optimized.json")
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found.")
        return

    with open(json_path, 'r') as f:
        config = json.load(f)

    # Reconstruct InputObject logic from aestimo_gui.py
    
    # Layers
    material_list = []
    for l in config["layers"]:
        th = float(l["thickness"])
        mat = l["material"]
        x = float(l["mole"])
        y = float(l.get("mole_y", 0.0))
        dop = float(l["doping"])
        dtype = l["doping_type"] 
        ltype = l["type"][0]
        
        if dtype == "i": dtype = "n"
        
        material_list.append([th, mat, x, y, dop, dtype, ltype])

    if not material_list:
        raise ValueError("Structure is empty.")
    
    # Physics / Solver
    if ":" in str(config["solver"]):
        scheme_id = int(config["solver"].split(":")[0])
    else:
        scheme_id = int(config["solver"])
        
    grid_step = float(config["grid_step"])
    max_pts = int(config["max_pts"])
    sub_e = int(config["sub_e"])
    sub_h = int(config["sub_h"])
    mat_sys = config["mat_sys"]
    
    T = float(config["temp"])
    F_app = float(config["field"]) * 1e5
    
    val_vmin = float(config["vmin"])
    val_vmax = float(config["vmax"])
    val_vstep = float(config["vstep"])
    
    bc_left = float(config["bc_left"])
    bc_right = float(config["bc_right"])

    # Doping Profile Calculation
    tot_thick = sum(row[0] for row in material_list) * 1e-9
    dx_m = grid_step * 1e-9
    n_max = int(tot_thick / dx_m)
    
    dop_arr = np.zeros(n_max)
    curr = 0
    for row in material_list:
        th_m = row[0] * 1e-9
        val = row[4]
        dtype = row[5]
        if dtype == 'p': val = -val
        
        steps = int(th_m / dx_m)
        end = min(curr + steps, n_max)
        dop_arr[curr:end] = val * 1e6 # cm-3 to m-3
        curr = end

    # Construct object instance directly instead of class definition tricky scope
    class InputObject:
        pass
        
    input_obj = InputObject()
    input_obj.T = T
    input_obj.T_val = T
    input_obj.computation_scheme = scheme_id
    input_obj.subnumber_h = sub_h
    input_obj.subnumber_e = sub_e
    input_obj.gridfactor = grid_step
    input_obj.maxgridpoints = max_pts
    input_obj.mat_type = mat_sys
    input_obj.dx = grid_step * 1e-9 # m
    
    input_obj.material = material_list
    
    input_obj.Fapplied = F_app
    
    # TAT Field
    input_obj.tat_field = float(config.get("tat_field", 1e10))
    
    # Solar Cell Optical Generation (Ensure units are m^-3)
    input_obj.G_optical = float(config.get("G_optical", 0.0)) * 1e6
    
    input_obj.vmax = val_vmax
    input_obj.vmin = val_vmin
    input_obj.Each_Step = val_vstep
    
    input_obj.surface = np.array([bc_left, bc_right])
    
    # Photovoltaic Mode Settings
    input_obj.photovoltaic_mode = True # Always enable for solar cell verification
    input_obj.work_function_left = float(config.get("wf_left", 4.2))
    input_obj.work_function_right = float(config.get("wf_right", 5.2))
    input_obj.surface_recomb = (float(config.get("s_left", 1e7)), float(config.get("s_right", 1e7)))
    
    input_obj.enable_polarization = config.get("enable_polarization", True)
    
    input_obj.Quantum_Regions = config.get("Quantum_Regions", False)
    qr_b_str = config.get("Quantum_Regions_boundary", "[[0.0, 0.0]]")
    try:
        if isinstance(qr_b_str, str):
            Quantum_Regions_boundary = np.array(json.loads(qr_b_str))
        else:
            Quantum_Regions_boundary = np.array(qr_b_str)
    except:
        Quantum_Regions_boundary = np.zeros((1,2))
    input_obj.Quantum_Regions_boundary = Quantum_Regions_boundary # Aestimo checks this attribute potentially? Or uses config.

    input_obj.dop_profile = dop_arr
    
    rs_val = float(config.get("rs", 0.0))
    rs_mode = config.get("rs_mode", "External (Fast)")
    
    if rs_mode == "Internal (Self-Consistent)":
        input_obj.Rs = rs_val
    else:
        input_obj.Rs = 0.0
    
    input_obj.device_area_m2 = float(config.get("area", 1e-4)) * 1e-4

    # Fake __file__ for output directory
    input_obj.__file__ = os.path.abspath("examples/ingan_solar_cell_script.py")

    print(f"Running simulation with G_optical={input_obj.G_optical:.2e} m^-3 s^-1")
    
    # Run
    # Force run_aestimo to not show plots
    aestimo.output_directory = os.path.join("examples", "ingan_solar_cell_output")
    if not os.path.isdir(aestimo.output_directory):
        os.makedirs(aestimo.output_directory)
        
    # run_aestimo expects class or object? It checks attributes. Instance is fine if StructureFrom handles it.
    # StructureFrom takes (input_obj, database). 
    # aestimo.run_aestimo takes (input_obj, drawFigures, show)
    # Inside run_aestimo: model = StructureFrom(input_obj, database)
    # StructureFrom accesses input_obj.material etc.
    # So instance should work fine.
    
    aestimo.run_aestimo(input_obj, drawFigures=False, show=False)
    
    # Load results
    out_dir = aestimo.output_directory
    data_path = os.path.join(out_dir, 'av_curr.dat')
    if not os.path.exists(data_path):
        print("Error: av_curr.dat not found. Simulation failed?")
        return

    data = np.loadtxt(data_path)
    
    # Process Results
    v = data[:,0]
    j = data[:,1] / 10.0 # A/m2 -> mA/cm2
    
    plt.figure()
    plt.plot(v, j, '-o')
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current Density (mA/cm2)')
    plt.title('InGaN Solar Cell Verification')
    plt.grid(True)
    plt.savefig('verify_solar.png')
    
    print("Raw I-V Data:")
    for v_val, j_val in zip(v, j):
        print(f"V: {v_val:.3f} V, J: {j_val:.4e} mA/cm2")
    
    # Calculate Metrics
    # Jsc = J at V=0
    try:
        Jsc = np.interp(0, v, j)
    except:
        Jsc = 0.0
    
    # Voc = V where J=0
    Voc = 0.0
    # Find zero crossing
    for i in range(len(v)-1):
        if j[i] * j[i+1] <= 0:
            slope = (j[i+1] - j[i]) / (v[i+1] - v[i])
            if slope != 0:
                Voc = v[i] - j[i] / slope
            break
            
    # Power
    P = v * j # mW/cm2
    
    # Solar cell generates power when I and V have opposite signs (4th quadrant usually)
    # In our simulation, J is light generated current (negative?) plus diode current.
    # At V=0, J = Jsc (positive or negative?). 
    # Let's check Jsc sign. 
    # If Jsc is positive, then we generate power when V>0 and J>0? No, usually Jsc is defined as current at V=0 short circuit.
    # Photodiode convention: I_light flows opposite to forward bias. So I is negative.
    # So power generation is V > 0 and I < 0.
    # If Jsc printed above is 0.037, it might be magnitude.
    # Let's look at the plot... or raw data.
    
    # We will assume generation quadrant is where V*J < 0
    # But effectively we want the max power point.
    
    # Indices where power is generated (V and J have opposite signs)
    gen_idx = np.where(v * j < 0)[0]
    
    if len(gen_idx) > 0:
        P_gen = P[gen_idx]
        Pmax = np.max(np.abs(P_gen))
    else:
        Pmax = 0.0
    
    # Fill Factor
    FF = 0.0
    if Jsc != 0 and Voc != 0:
        FF = Pmax / (np.abs(Jsc) * np.abs(Voc)) * 100
        
    print("-" * 30)
    print("InGaN Solar Cell Simulation Results")
    print("-" * 30)
    print(f"{'Parameter':<10} | {'Value':<10} | {'Unit':<5}")
    print("-" * 30)
    print(f"{'Jsc':<10} | {abs(Jsc):<10.4f} | mA/cm2")
    print(f"{'Voc':<10} | {Voc:<10.4f} | V")
    print(f"{'FF':<10} | {FF:<10.2f} | %")
    print(f"{'Pmax':<10} | {Pmax:<10.4f} | mW/cm2")
    print("-" * 30)

if __name__ == "__main__":
    verify_solar()
