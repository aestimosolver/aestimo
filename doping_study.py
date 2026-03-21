
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# Add current directory to path
sys.path.append(os.getcwd())

from aestimo import run_aestimo
from characterize_solar import analyze_iv_curve

def run_test(doping, disable_auger=False):
    print(f"\n--- Testing Doping: {doping:.1e} cm^-3 (Auger={not disable_auger}) ---")
    
    # Monkey-patch database if needed
    import database
    if disable_auger:
        orig_cn0 = database.materialproperty['GaN']['Cn0']
        orig_cp0 = database.materialproperty['GaN']['Cp0']
        database.materialproperty['GaN']['Cn0'] = 0.0
        database.materialproperty['GaN']['Cp0'] = 0.0
        database.materialproperty['InN']['Cn0'] = 0.0
        database.materialproperty['InN']['Cp0'] = 0.0
        # For alloys, they use material property values in create_structure_arrays
    
    # Configuration
    config = {
        "layers": [
            {"thickness": 30, "material": "InGaN", "mole": 0.2, "doping": doping, "doping_type": "p", "type": "barrier"},
            {"thickness": 450, "material": "InGaN", "mole": 0.2, "doping": doping, "doping_type": "n", "type": "barrier"}
        ],
        "T": 300,
        "vmax": 2.4,
        "vmin": 0.0,
        "Each_Step": 0.02,
        "grid_step": 0.2,
        "comp_scheme": 7,
        "enable_polarization": False,
        "tau": 2e-6,
        "G_optical": 5e20,
        "work_function_left": 7.0,
        "work_function_right": 4.0,
        "mat_type": "Wurtzite"
    }

    class InputObject:
        def __init__(self, d):
            for k, v in d.items():
                setattr(self, k, v)
            material_list = []
            for l in d["layers"]:
                th = float(l["thickness"])
                material_list.append([th, l["material"], float(l["mole"]), 0.0, float(l["doping"]), l["doping_type"], l["type"][0]])
            self.material = material_list
            self.mat_sys = "Wurtzite"

    in_obj = InputObject(config)
    _, model, res, _ = run_aestimo(in_obj)
    
    if disable_auger:
        database.materialproperty['GaN']['Cn0'] = orig_cn0
        database.materialproperty['GaN']['Cp0'] = orig_cp0
        # InN reset omitted for brevity but should be done if reused

    # Metrics
    metrics = analyze_iv_curve(res.Va_t, res.av_curr)
    print(f"Results: Jsc={metrics['jsc']:.2f}, Voc={metrics['voc']:.2f}, FF={metrics['ff']:.2f}%, Pmax={metrics['pmpp']:.2f}")
    return metrics

dopings = [1e17, 1e16]
results = []
# Test 1: Standard Auger
results.append(run_test(1e17, disable_auger=False))
# Test 2: Lower doping
results.append(run_test(1e16, disable_auger=False))
# Test 3: Standard doping, No Auger
results.append(run_test(1e17, disable_auger=True))

# Quick Plot
plt.figure(figsize=(8,5))
plt.plot([1, 2, 3], [r['ff'] for r in results], 'o-', label='FF (%)')
plt.plot([1, 2, 3], [r['jsc'] * 10 for r in results], 's-', label='Jsc (mA/cm2) x 10')
plt.xticks([1, 2, 3], ['1e17', '1e16', '1e17 (No Auger)'])
plt.xscale('log')
plt.xlabel('Doping (cm-3)')
plt.ylabel('Value')
plt.title('Doping Impact on InGaN Solar Cell (No Polarization)')
plt.legend()
plt.grid(True)
plt.savefig('doping_study.png')
print("\nStudy complete. saved to doping_study.png")
