import os
import sys
import json
import numpy as np
import time

sys.path.insert(0, os.path.abspath('.'))
import aestimo
from aeslibs.experimental_validation import load_current_from_avcurr, apply_parasitic_resistances

def get_base_config():
    with open('examples/untitled_project.json', 'r') as f:
        config = json.load(f)
        
    config['enable_polarization'] = True
    config['comp_scheme'] = 9 # Use Gummel map
    
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

def extract_metrics(v, j):
    j_ma = j * 1000.0
    if len(v) < 2 or len(j) < 2: return 0,0,0,0
    jsc = np.interp(0.0, v, j_ma)
    try:
        voc = np.interp(0.0, j_ma, v)
    except:
        voc = 0.0
    p = v * j_ma
    pmax = np.min(p)
    ff = abs(pmax / (jsc * voc)) if jsc * voc != 0 else 0.0
    return abs(jsc), voc, ff * 100.0, abs(pmax)

def run_defect_sweep():
    print("Running Defect Sweep...")
    # Sweep taun0 and taup0 logarithmically from 1e-9 to 1e-6 (4 points)
    tau_vals = np.logspace(-9, -6, 4)
    results = []
    for tau in tau_vals:
        config = get_base_config()
        config['taun0'] = float(tau)
        config['taup0'] = float(tau)
        config['__file__'] = f'defect_{tau:.1e}'
        config['G_optical'] = 1e21
        print(f" Simulating tau={tau:.1e}...")
        try:
            aestimo.run_aestimo(config, drawFigures=False, show=False)
            out_dir = f"defect_{tau:.1e}_output"
            sim_v_int, sim_i_int = load_current_from_avcurr(out_dir, device_area_cm2=config['device_area'])
            sim_v, sim_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=float(config.get('rs', 11.7)), Rsh=float(config.get('rsh', 5000)))
            jsc, voc, ff, pmax = extract_metrics(sim_v, sim_i / config['device_area'])
            results.append((tau, jsc, voc, ff, pmax))
        except Exception as e:
            print(f"  Failed for tau={tau}: {e}")
            
    with open("results/defect_sweep/defect_results.csv", "w") as f:
        f.write("Tau(s),Jsc(mA/cm2),Voc(V),FF(%),Pmax\n")
        for r in results: f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}\n")
    print("Defect sweep completed.")

def run_indium_sweep():
    print("Running Indium Sweep...")
    # Sweep mole of both layers from 0.40 to 0.70 (4 points)
    moles = np.linspace(0.40, 0.70, 4)
    results = []
    for m in moles:
        config = get_base_config()
        config['material'][0][2] = float(m)
        config['material'][1][2] = float(m)
        config['__file__'] = f'indium_{m:.2f}'
        config['G_optical'] = 1e21
        print(f" Simulating indium={m:.2f}...")
        try:
            aestimo.run_aestimo(config, drawFigures=False, show=False)
            out_dir = f"indium_{m:.2f}_output"
            sim_v_int, sim_i_int = load_current_from_avcurr(out_dir, device_area_cm2=config['device_area'])
            sim_v, sim_i = apply_parasitic_resistances(sim_v_int, sim_i_int, Rs=float(config.get('rs', 11.7)), Rsh=float(config.get('rsh', 5000)))
            jsc, voc, ff, pmax = extract_metrics(sim_v, sim_i / config['device_area'])
            results.append((m, jsc, voc, ff, pmax))
        except Exception as e:
            print(f"  Failed for indium={m}: {e}")
            
    with open("results/high_indium/indium_results.csv", "w") as f:
        f.write("Mole,Jsc(mA/cm2),Voc(V),FF(%),Pmax\n")
        for r in results: f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}\n")
    print("Indium sweep completed.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "defect": run_defect_sweep()
        elif sys.argv[1] == "indium": run_indium_sweep()
    else:
        print("Please specify sweep: defect or indium")
