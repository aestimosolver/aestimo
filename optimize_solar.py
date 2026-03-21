import sys, os
sys.path.append(os.getcwd())
import aestimo
import numpy as np
import matplotlib.pyplot as plt
import json
from itertools import product

def run_solar_simulation(In_fraction, p_doping, n_doping, p_thick, n_thick, G_opt=1e21):
    """Run a single solar cell simulation with given parameters."""
    
    # Build config
    config = {
        "layers": [
            {
                "material": "InGaN",
                "mole": str(In_fraction),
                "mole_y": "0.0",
                "thickness": str(p_thick),
                "doping": str(p_doping),
                "doping_type": "p",
                "type": "barrier"
            },
            {
                "material": "InGaN",
                "mole": str(In_fraction),
                "mole_y": "0.0",
                "thickness": str(n_thick),
                "doping": str(n_doping),
                "doping_type": "n",
                "type": "barrier"
            }
        ],
        "temp": "300.0",
        "field": "0.0",
        "bc_left": "0.0",
        "bc_right": "0.0",
        "vmin": "0.0",
        "vmax": "1.6",
        "vstep": "0.04",
        "solver": "7",
        "grid_step": "1.0",
        "max_pts": "200000",
        "mat_sys": "Wurtzite",
        "sub_e": "2",
        "sub_h": "2",
        "area": "5.0e-4",
        "rs": "0.0",
        "rsh": "1e6",
        "tat_field": "5.0e6",
        "taun0": "5.0e-7",
        "taup0": "5.0e-7",
        "G_optical": str(G_opt),
        "enable_polarization": False
    }
    
    # Build InputObject
    class InputObject:
        pass
    
    input_obj = InputObject()
    input_obj.computation_scheme = 7
    input_obj.gridfactor = 1.0
    input_obj.maxgridpoints = 200000
    input_obj.mat_type = "Wurtzite"
    input_obj.dx = 1e-9
    input_obj.T = 300.0
    input_obj.subnumber_e = 2
    input_obj.subnumber_h = 2
    
    # Material list
    mat_list = []
    for l in config["layers"]:
        dop = float(l["doping"])
        if l["doping_type"] == 'p': dop = -dop
        mat_list.append([
            float(l["thickness"]),
            l["material"],
            float(l["mole"]),
            float(l.get("mole_y", 0)),
            dop,
            l["doping_type"],
            l["type"][0]
        ])
    input_obj.material = mat_list
    
    input_obj.Fapplied = 0.0
    input_obj.tat_field = 5e6
    input_obj.G_optical = float(G_opt) * 1e6  # cm-3 to m-3
    input_obj.vmax = 1.6
    input_obj.vmin = 0.0
    input_obj.Each_Step = 0.04
    input_obj.surface = np.array([0.0, 0.0])
    input_obj.Quantum_Regions = False
    input_obj.Quantum_Regions_boundary = np.zeros((1,2))
    input_obj.Rs = 0.0
    input_obj.device_area_m2 = 5e-8
    
    # Doping profile
    tot_thick = sum(row[0] for row in mat_list) * 1e-9
    dx_m = 1e-9
    n_max = int(tot_thick / dx_m)
    
    dop_arr = np.zeros(n_max)
    curr = 0
    for row in mat_list:
        th_m = row[0] * 1e-9
        val = row[4]
        steps = int(th_m / dx_m)
        end = min(curr + steps, n_max)
        dop_arr[curr:end] = val * 1e6
        curr = end
    
    input_obj.dop_profile = dop_arr
    input_obj.__file__ = os.path.abspath("optimize_temp.py")
    
    # Run simulation
    aestimo.output_directory = "optimize_output"
    if not os.path.exists(aestimo.output_directory):
        os.makedirs(aestimo.output_directory)
    
    try:
        _, model, res, _ = aestimo.run_aestimo(input_obj, drawFigures=False, show=False)
        
        # Load IV data
        data = np.loadtxt(os.path.join(aestimo.output_directory, 'av_curr.dat'))
        v = data[:,0]
        j = data[:,1] / 10.0  # A/m2 -> mA/cm2
        
        # Calculate metrics
        Jsc = np.interp(0, v, j)
        
        # Voc
        Voc = 0.0
        for i in range(len(v)-1):
            if j[i] * j[i+1] <= 0:
                slope = (j[i+1] - j[i]) / (v[i+1] - v[i])
                if slope != 0:
                    Voc = v[i] - j[i] / slope
                break
        
        # Power
        P = v * j
        gen_idx = np.where(v * j < 0)[0]
        if len(gen_idx) > 0:
            Pmax = np.max(np.abs(P[gen_idx]))
        else:
            Pmax = 0.0
        
        # FF
        FF = 0.0
        if abs(Jsc) > 0 and abs(Voc) > 0:
            FF = Pmax / (abs(Jsc) * abs(Voc)) * 100
        
        return {
            'Jsc': abs(Jsc),
            'Voc': Voc,
            'FF': FF,
            'Pmax': Pmax,
            'success': True
        }
        
    except Exception as e:
        print(f"Simulation failed: {e}")
        return {
            'Jsc': 0,
            'Voc': 0,
            'FF': 0,
            'Pmax': 0,
            'success': False
        }

def optimize():
    """Automated parameter optimization."""
    
    print("=" * 60)
    print("InGaN Solar Cell Automated Optimization")
    print("=" * 60)
    print("\nTarget Performance:")
    print("  Voc: 1.5 - 3.0 V")
    print("  Jsc: 5 - 20 mA/cm²")
    print("  FF: > 70%")
    print("\n" + "=" * 60)
    
    # Parameter ranges
    In_fractions = [0.20, 0.25, 0.30, 0.35]  # Lower In = higher Eg
    p_dopings = [1e17, 5e17, 1e18]  # cm-3
    n_dopings = [1e18, 5e18, 1e19]  # cm-3
    
    # Fixed for now
    p_thick = 40  # nm
    n_thick = 60  # nm
    G_opt = 1e21  # cm-3 s-1
    
    results = []
    best_result = None
    best_score = -np.inf
    
    total_runs = len(In_fractions) * len(p_dopings) * len(n_dopings)
    run_count = 0
    
    print(f"\nRunning {total_runs} simulations...")
    print("-" * 60)
    
    for In, p_dop, n_dop in product(In_fractions, p_dopings, n_dopings):
        run_count += 1
        print(f"\n[{run_count}/{total_runs}] In={In:.2f}, p={p_dop:.1e}, n={n_dop:.1e} cm-3")
        
        result = run_solar_simulation(In, p_dop, n_dop, p_thick, n_thick, G_opt)
        
        if result['success']:
            print(f"  -> Jsc={result['Jsc']:.4f} mA/cm2, Voc={result['Voc']:.4f} V, FF={result['FF']:.2f}%")
            
            # Score function (prioritize Voc and Jsc in target range)
            score = 0
            if 1.5 <= result['Voc'] <= 3.0:
                score += 100
            if 5 <= result['Jsc'] <= 20:
                score += 100
            if result['FF'] > 70:
                score += 50
            
            # Bonus for being closer to middle of ranges
            score += (1 - abs(result['Voc'] - 2.25) / 0.75) * 20  # Voc closer to 2.25V
            score += (1 - abs(result['Jsc'] - 12.5) / 7.5) * 20   # Jsc closer to 12.5 mA/cm2
            
            result['score'] = score
            result['In'] = In
            result['p_dop'] = p_dop
            result['n_dop'] = n_dop
            results.append(result)
            
            if score > best_score:
                best_score = score
                best_result = result
                print(f"  * NEW BEST (score={score:.1f})")
        else:
            print(f"  X Failed")
    
    print("\n" + "=" * 60)
    print("OPTIMIZATION COMPLETE")
    print("=" * 60)
    
    if best_result:
        print(f"\n* BEST CONFIGURATION:")
        print(f"  In fraction: {best_result['In']:.2f}")
        print(f"  p-doping: {best_result['p_dop']:.2e} cm-3")
        print(f"  n-doping: {best_result['n_dop']:.2e} cm-3")
        print(f"\n  Performance:")
        print(f"  Jsc: {best_result['Jsc']:.4f} mA/cm2")
        print(f"  Voc: {best_result['Voc']:.4f} V")
        print(f"  FF: {best_result['FF']:.2f} %")
        print(f"  Pmax: {best_result['Pmax']:.4f} mW/cm2")
        print(f"  Score: {best_result['score']:.1f}")
        
        # Save best config
        best_config = {
            "layers": [
                {
                    "material": "InGaN",
                    "mole": str(best_result['In']),
                    "mole_y": "0.0",
                    "thickness": str(p_thick),
                    "doping": str(best_result['p_dop']),
                    "doping_type": "p",
                    "type": "barrier"
                },
                {
                    "material": "InGaN",
                    "mole": str(best_result['In']),
                    "mole_y": "0.0",
                    "thickness": str(n_thick),
                    "doping": str(best_result['n_dop']),
                    "doping_type": "n",
                    "type": "barrier"
                }
            ],
            "temp": "300.0",
            "field": "0.0",
            "bc_left": "0.0",
            "bc_right": "0.0",
            "vmin": "0.0",
            "vmax": "1.6",
            "vstep": "0.04",
            "solver": "7: Standard DD",
            "grid_step": "1.0",
            "max_pts": "200000",
            "mat_sys": "Wurtzite",
            "sub_e": "2",
            "sub_h": "2",
            "area": "5.0e-4",
            "rs": "0.0",
            "rsh": "1000000.0",
            "tat_field": "5.0e6",
            "taun0": "5.0e-7",
            "taup0": "5.0e-7",
            "G_optical": str(G_opt),
            "enable_polarization": False,
            "exp_file": ""
        }
        
        with open("examples/ingan_solar_optimized.json", 'w') as f:
            json.dump(best_config, f, indent=4)
        
        print(f"\nSaved to examples/ingan_solar_optimized.json")
        
        # Plot top 5 results
        sorted_results = sorted([r for r in results if r['success']], key=lambda x: x['score'], reverse=True)[:5]
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        for i, r in enumerate(sorted_results):
            label = f"In={r['In']:.2f}, p={r['p_dop']:.1e}, n={r['n_dop']:.1e}"
            axes[0].scatter(r['Voc'], r['Jsc'], s=100, label=label if i < 3 else "")
        axes[0].axhspan(5, 20, alpha=0.2, color='green', label='Target Jsc')
        axes[0].axvspan(1.5, 3.0, alpha=0.2, color='blue', label='Target Voc')
        axes[0].set_xlabel('Voc (V)')
        axes[0].set_ylabel('Jsc (mA/cm2)')
        axes[0].legend(fontsize=8)
        axes[0].grid(True, alpha=0.3)
        axes[0].set_title('Top 5 Configurations')
        
        # In vs Voc
        In_vals = [r['In'] for r in results if r['success']]
        Voc_vals = [r['Voc'] for r in results if r['success']]
        axes[1].scatter(In_vals, Voc_vals, alpha=0.6)
        axes[1].set_xlabel('In Fraction')
        axes[1].set_ylabel('Voc (V)')
        axes[1].axhspan(1.5, 3.0, alpha=0.2, color='green')
        axes[1].grid(True, alpha=0.3)
        axes[1].set_title('In Composition vs Voc')
        
        # Doping vs Jsc
        n_dop_vals = [r['n_dop'] for r in results if r['success']]
        Jsc_vals = [r['Jsc'] for r in results if r['success']]
        axes[2].scatter(n_dop_vals, Jsc_vals, alpha=0.6)
        axes[2].set_xlabel('n-doping (cm-3)')
        axes[2].set_ylabel('Jsc (mA/cm2)')
        axes[2].set_xscale('log')
        axes[2].axhspan(5, 20, alpha=0.2, color='green')
        axes[2].grid(True, alpha=0.3)
        axes[2].set_title('n-Doping vs Jsc')
        
        plt.tight_layout()
        plt.savefig('optimization_results.png', dpi=150)
        print(f"Saved optimization_results.png")
        
    else:
        print("\nX No successful simulations")
    
    print("\n" + "=" * 60)

if __name__ == "__main__":
    optimize()
