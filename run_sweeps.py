import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath('.'))

import aestimo
from aeslibs.experimental_validation import (
    load_current_from_avcurr,
    apply_parasitic_resistances,
    calculate_ideality_factor
)

def get_base_config():
    with open('examples/untitled_project.json', 'r') as f:
        config = json.load(f)
        
    config['G_optical'] = 1e21
    config['enable_polarization'] = True
    config['comp_scheme'] = 9 # Use more stable Gummel map solver
    
    # Format layers
    mat_list = []
    for l in config.get('layers', []):
        mat_list.append([
            float(l.get('thickness', 0.0)),
            l.get('material', ''),
            float(l.get('mole', 0.0)),
            float(l.get('mole_y', 0.0)),
            float(l.get('doping', 0.0)),
            l.get('doping_type', 'n'),
            l.get('type', 'b')
        ])
    config['material'] = mat_list
    if 'layers' in config: del config['layers']
    
    config['Fapplied'] = float(config.get('field', 0.0))
    config['Each_Step'] = float(config.get('vstep', 0.05))
    config['device_area'] = float(config.get('area', 1e-4))
    if 'mat_sys' in config: config['mat_type'] = config['mat_sys']
    if 'tat_field' in config: config['tat_field'] = float(config['tat_field'])
    
    config['enable_experimental_validation'] = False
    
    return config

def extract_solar_metrics(v, j):
    # j is in A/cm^2, v in Volts
    # Convert j to mA/cm^2 for standard reporting
    j_ma = j * 1000.0
    
    if len(v) < 2 or len(j) < 2:
        return 0, 0, 0, 0
        
    # Jsc is J at V=0
    jsc = np.interp(0.0, v, j_ma)
    
    # Voc is V at J=0
    # Find zero crossing
    try:
        # np.interp needs x to be monotonically increasing, and we want to interpolate V(J)
        # However, J goes from negative (Jsc) to positive
        v_interp = np.interp(0.0, j_ma, v)
        voc = v_interp
    except:
        voc = 0.0
        
    # Power
    p = v * j_ma
    pmax = np.min(p) # Power is negative in generation convention
    
    # FF
    if jsc * voc != 0:
        ff = abs(pmax / (jsc * voc))
    else:
        ff = 0.0
        
    # Efficiency (G_optical = 1e21, we can approximate incident power, but for now just relative eff)
    eff = abs(pmax) # arbitrary units since we don't know incident spectrum power perfectly
    
    return abs(jsc), voc, ff * 100.0, eff

def run_thickness_sweep():
    print("Running Thickness Sweep...")
    thicknesses = np.arange(20, 101, 20)
    results = []
    
    for t in thicknesses:
        config = get_base_config()
        # Layer 1 is barrier p-type (40nm). Layer 2 is barrier n-type (was 60nm).
        # We sweep layer 2 thickness
        config['material'][1][0] = float(t)
        config['__file__'] = f'thickness_{t}nm'
        
        print(f" Simulating thickness {t} nm...")
        try:
            aestimo.run_aestimo(config, drawFigures=False, show=False)
            
            output_dir = f'thickness_{t}nm_output'
            area = config['device_area']
            sim_v_int, sim_i_int = load_current_from_avcurr(output_dir, device_area_cm2=area)
            sim_v, sim_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=float(config.get('rs', 11.7)), Rsh=float(config.get('rsh', 5000)))
            
            j = sim_i / area
            jsc, voc, ff, eff = extract_solar_metrics(sim_v, j)
            results.append((t, jsc, voc, ff, eff))
            
        except Exception as e:
            print(f"  Failed for thickness {t}: {e}")
            
    # Save results
    with open("results/thickness_sweep/thickness_results.csv", "w") as f:
        f.write("Thickness(nm),Jsc(mA/cm2),Voc(V),FF(%),Pmax\n")
        for r in results:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}\n")
    print("Thickness sweep completed.")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "thickness":
        run_thickness_sweep()
    else:
        print("Please specify a sweep: thickness")
