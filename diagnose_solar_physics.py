import sys, os
sys.path.append(os.getcwd())
import aestimo
import numpy as np
import matplotlib.pyplot as plt
import json

def diagnose():
    # Load config
    json_path = os.path.join("examples", "ingan_solar_cell.json")
    with open(json_path, 'r') as f:
        config = json.load(f)

    # Simplified InputObject construction
    class InputObject:
        pass
    
    input_obj = InputObject()
    # Copy essential params
    input_obj.computation_scheme = int(config["solver"].split(":")[0] if ":" in str(config["solver"]) else config["solver"])
    input_obj.gridfactor = float(config["grid_step"])
    input_obj.maxgridpoints = int(config["max_pts"])
    input_obj.mat_type = config["mat_sys"]
    input_obj.dx = input_obj.gridfactor * 1e-9
    input_obj.T = float(config["temp"])
    input_obj.subnumber_e = int(config["sub_e"])
    input_obj.subnumber_h = int(config["sub_h"])
    
    # Construct material list
    mat_list = []
    for l in config["layers"]:
        dop = float(l["doping"])
        if l["doping_type"] == 'p': dop = -dop
        if l["doping_type"] == 'i': dop = 0
        mat_list.append([
            float(l["thickness"]), 
            l["material"], 
            float(l["mole"]), 
            float(l.get("mole_y",0)), 
            dop, 
            l["doping_type"], 
            l["type"][0]
        ])
    input_obj.material = mat_list
    
    # Physics flags
    input_obj.Fapplied = 0.0 # Force zero for equilibrium diagnostic
    input_obj.tat_field = 1e10
    input_obj.G_optical = 0.0 # Dark/Equilibrium
    input_obj.vmax = 0.0
    input_obj.vmin = 0.0
    input_obj.Each_Step = 0.01
    input_obj.surface = np.array([0.0, 0.0]) # Flat band / Ohmic?
    input_obj.Quantum_Regions = False
    input_obj.Quantum_Regions_boundary = np.zeros((1,2))
    input_obj.Rs = 0.0
    input_obj.device_area_m2 = 1e-8
    
    # Calculate doping profile
    tot_thick = sum(row[0] for row in mat_list) * 1e-9
    dx_m = input_obj.gridfactor * 1e-9
    n_max = int(tot_thick / dx_m)
    
    dop_arr = np.zeros(n_max)
    curr = 0
    for row in mat_list:
        th_m = row[0] * 1e-9
        val = row[4]  # Already signed
        steps = int(th_m / dx_m)
        end = min(curr + steps, n_max)
        dop_arr[curr:end] = val * 1e6  # cm-3 to m-3
        curr = end
    
    input_obj.dop_profile = dop_arr
    
    # Fake file path
    input_obj.__file__ = os.path.abspath("examples/diagnostic.py")
    
    print("Running Aestimo simulation at equilibrium (0V)...")
    aestimo.output_directory = "diagnostic_output"
    if not os.path.exists(aestimo.output_directory): 
        os.makedirs(aestimo.output_directory)
    
    # Run simulation
    _, model, res, _ = aestimo.run_aestimo(input_obj, drawFigures=False, show=False)
    
    print(f"\n--- Diagnostic Report ---")
    print(f"Temperature: {input_obj.T} K")
    print(f"Material System: {input_obj.mat_type}")
    print(f"Layers: {len(input_obj.material)}")
    
    # Print layer structure
    print(f"\n--- Layer Structure ---")
    for i, layer in enumerate(input_obj.material):
        thick, mat, x, y, dop, dtype, ltype = layer
        print(f"Layer {i+1}: {thick:.1f} nm {mat} (x={x:.2f}), doping={abs(dop):.2e} cm^-3 ({dtype}-type)")
    
    # Load and analyze band diagram
    try:
        data = np.loadtxt(os.path.join(aestimo.output_directory, 'band_grad.dat'))
        x = data[:,0]  # nm
        Ec = data[:,1]  # eV
        Ev = data[:,2]  # eV
        Efn = data[:,3]  # eV
        Efp = data[:,4]  # eV
        
        # Calculate bandgap
        Eg = Ec - Ev
        Eg_avg = np.mean(Eg)
        
        # Plot
        plt.figure(figsize=(12,7))
        plt.subplot(2,1,1)
        plt.plot(x, Ec, 'b-', label='Ec', linewidth=2)
        plt.plot(x, Ev, 'r-', label='Ev', linewidth=2)
        plt.plot(x, Efn, 'g--', label='Efn', linewidth=1.5)
        plt.plot(x, Efp, 'm--', label='Efp', linewidth=1.5)
        plt.legend()
        plt.ylabel('Energy (eV)')
        plt.title('Equilibrium Band Diagram (0V, Dark)')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(2,1,2)
        plt.plot(x, Eg, 'k-', linewidth=2)
        plt.ylabel('Bandgap (eV)')
        plt.xlabel('Position (nm)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('diag_equi.png', dpi=150)
        print(f"\nBand diagram saved to diag_equi.png")
        
        # Calculate built-in potential
        # Vbi ≈ difference in conduction band edge between p and n sides
        # At equilibrium, Fermi level is flat, so Ec difference gives Vbi
        Vbi = abs(Ec[0] - Ec[-1])
        
        print(f"\n--- Band Structure Analysis ---")
        print(f"Bandgap (average): {Eg_avg:.4f} eV")
        print(f"Bandgap (left edge): {Eg[0]:.4f} eV")
        print(f"Bandgap (right edge): {Eg[-1]:.4f} eV")
        print(f"\nEc(left, p-side): {Ec[0]:.4f} eV")
        print(f"Ec(right, n-side): {Ec[-1]:.4f} eV")
        print(f"Built-in Potential (Vbi): {Vbi:.4f} V")
        
        # Theoretical max Voc
        Voc_max = min(Vbi, Eg_avg - 0.3)  # Upper limit (accounting for recombination losses)
        print(f"\nTheoretical max Voc ≈ {Voc_max:.4f} V")
        print(f"Target Voc range: 1.5-3.0 V")
        
        # Diagnosis
        print(f"\n--- Diagnosis ---")
        if Vbi < 1.5:
            print(f"⚠️  WARNING: Vbi ({Vbi:.4f} V) is below target range!")
            print(f"\nPossible causes:")
            print(f"  1. Bandgap too narrow: {Eg_avg:.4f} eV")
            print(f"     → For Voc > 1.5V, need Eg > 1.8 eV")
            print(f"  2. Doping too low (check layer structure above)")
            print(f"\nRecommendations:")
            if Eg_avg < 1.8:
                current_In = input_obj.material[0][2]
                target_In = max(0.1, current_In - 0.1)
                print(f"  → Reduce In content from {current_In:.2f} to ~{target_In:.2f}")
                print(f"     (InGaN bandgap: In=0.3 → 2.0eV, In=0.5 → 1.5eV, In=0.7 → 1.0eV)")
            print(f"  → Increase doping to 1e18 - 1e19 cm^-3")
        else:
            print(f"✓ Vbi ({Vbi:.4f} V) is in acceptable range")
            print(f"  Low Voc in actual simulation may be due to:")
            print(f"  - High recombination (check carrier lifetimes)")
            print(f"  - Poor carrier collection (check mobilities)")
            print(f"  - Incorrect current scaling")
            
    except Exception as e:
        print(f"\nError reading band diagram: {e}")
        import traceback
        traceback.print_exc()
    
    print("\nDiagnostic complete.")

if __name__ == "__main__":
    diagnose()
