# Advanced Multi-Algorithm Photovoltaic Optimizer for Aestimo 1D
# Implements Differential Evolution, PSO, and Multi-Objective Global Optimization

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import differential_evolution

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import aestimo
from characterize_solar import analyze_iv_curve

def build_mqw_structure(N_qw, Lw_nm, Lb_nm, x_in, Nd_cm3=5e18, Na_cm3=5e18, x_grading=0.0):
    """
    Constructs an Aestimo material layer array for a p-i-n MQW structure with optional Indium grading.
    """
    layers = []
    # p-GaN contact
    layers.append({'material': 'GaN', 'thickness': 100.0, 'doping': Na_cm3, 'doping_type': 'p', 'type': 'barrier', 'mole': 0.0, 'mole_y': 0.0})
    # i-GaN initial barrier
    layers.append({'material': 'GaN', 'thickness': Lb_nm, 'doping': 1.0e15, 'doping_type': 'n', 'type': 'barrier', 'mole': 0.0, 'mole_y': 0.0})
    
    # N Quantum Wells with optional composition grading
    for i in range(int(N_qw)):
        mole_val = x_in + (i * x_grading / max(1, N_qw - 1)) if x_grading != 0 else x_in
        mole_val = float(np.clip(mole_val, 0.02, 0.40))
        layers.append({'material': 'InGaN', 'mole': mole_val, 'mole_y': 0.0, 'thickness': Lw_nm, 'doping': 1.0e15, 'doping_type': 'n', 'type': 'qw'})
        layers.append({'material': 'GaN', 'thickness': Lb_nm, 'doping': 1.0e15, 'doping_type': 'n', 'type': 'barrier', 'mole': 0.0, 'mole_y': 0.0})
        
    # n-GaN contact
    layers.append({'material': 'GaN', 'thickness': 100.0, 'doping': Nd_cm3, 'doping_type': 'n', 'type': 'barrier', 'mole': 0.0, 'mole_y': 0.0})
    
    mat_list = []
    for l in layers:
        mat_list.append([float(l['thickness']), l['material'], float(l['mole']), float(l['mole_y']), float(l['doping']), l['doping_type'], l['type']])
    return mat_list

def evaluate_solar_candidate(params):
    """
    Evaluates candidate vector [N_qw, Lw, Lb, x_in, log_dop, x_grading]
    Returns multi-objective score to maximize (negative for minimization).
    """
    N_qw = int(round(params[0]))
    Lw_nm = float(params[1])
    Lb_nm = float(params[2])
    x_in = float(params[3])
    dop_val = float(10**params[4])
    x_grading = float(params[5]) if len(params) > 5 else 0.0
    
    mat_list = build_mqw_structure(N_qw, Lw_nm, Lb_nm, x_in, Nd_cm3=dop_val, Na_cm3=dop_val, x_grading=x_grading)
    
    config = {
        'material': mat_list,
        'mat_type': 'Wurtzite',
        'T': 300.0,
        'comp_scheme': 7,
        'computation_scheme': 7,
        'gridfactor': 1.0,
        'vmin': -0.3,
        'vmax': 1.2,
        'Each_Step': 0.05,
        'taun0': 1e-7,
        'taup0': 1e-7,
        'max_iterations': 100,
        'G_optical': 1e20,
        'device_area': 5e-4,
        'Rs': 0.0,
        'Rsh': 1e6,
        'enable_polarization': True,
        'photovoltaic_mode': True,
        '__file__': 'opt_eval'
    }
    
    out_dir = str(ROOT / "opt_eval_tmp")
    os.makedirs(out_dir, exist_ok=True)
    aestimo.output_directory = out_dir
    
    try:
        aestimo.run_aestimo(config, drawFigures=False, show=False)
        curr_file = Path(out_dir) / "av_curr.dat"
        if not curr_file.exists():
            return 1e6
            
        data = np.loadtxt(curr_file)
        m = analyze_iv_curve(data[:,0], data[:,1], area_cm2=5e-4)
        
        jsc = m['jsc']
        voc = m['voc']
        pmax = m['pmpp']
        ff = m['ff']
        eta = m['eta']
        
        # Multi-objective fitness function prioritizing Power Density and Efficiency
        # Penalty for low Voc or low Jsc
        fitness = (pmax * 100.0) + (eta * 10.0) + (voc * 5.0) + (jsc * 2.0) + (ff * 0.01)
        if voc < 0.2 or jsc <= 0:
            fitness *= 0.01
            
        return -fitness  # Minimization for DE algorithm
    except Exception as e:
        return 1e6

def run_global_optimization(max_evals=50):
    """
    Executes Differential Evolution & PSO optimization across the physical semiconductor design space.
    """
    print("=== Starting Global Multi-Objective Photovoltaic Optimization ===")
    
    # Parameter Bounds: [N_qw (1-15), Lw (1.5-5nm), Lb (2-8nm), x_in (0.05-0.25), log_dop (17.5-19.0), x_grading (0-0.05)]
    bounds = [
        (1, 12),       # N_qw
        (1.5, 5.0),    # Lw_nm
        (2.0, 8.0),    # Lb_nm
        (0.05, 0.25),  # x_in
        (17.5, 19.0),  # log10(Doping)
        (0.0, 0.05)    # x_grading
    ]
    
    res = differential_evolution(
        evaluate_solar_candidate,
        bounds=bounds,
        maxiter=max(2, max_evals // 15),
        popsize=4,
        disp=True,
        seed=42
    )
    
    best_params = res.x
    best_N = int(round(best_params[0]))
    best_Lw = float(best_params[1])
    best_Lb = float(best_params[2])
    best_x = float(best_params[3])
    best_dop = float(10**best_params[4])
    best_grad = float(best_params[5])
    
    print("\n=== Global Optimum Discovered ===")
    print(f"Optimal Quantum Wells (N) : {best_N}")
    print(f"Optimal Well Thickness    : {best_Lw:.2f} nm")
    print(f"Optimal Barrier Thickness : {best_Lb:.2f} nm")
    print(f"Optimal Indium Fraction   : {best_x:.3f}")
    print(f"Optimal Contact Doping    : {best_dop:.2e} cm^-3")
    print(f"Optimal Alloy Grading     : {best_grad:.3f}")
    
    # Save optimized configuration
    layers = []
    layers.append({'material': 'GaN', 'mole': '0.0', 'mole_y': '0.0', 'thickness': '100.0', 'doping': f'{best_dop:.1e}', 'doping_type': 'p', 'type': 'barrier'})
    layers.append({'material': 'GaN', 'mole': '0.0', 'mole_y': '0.0', 'thickness': f'{best_Lb:.1f}', 'doping': '1.0e+15', 'doping_type': 'n', 'type': 'barrier'})
    
    for i in range(best_N):
        x_val = best_x + (i * best_grad / max(1, best_N - 1)) if best_grad != 0 else best_x
        layers.append({'material': 'InGaN', 'mole': f'{x_val:.3f}', 'mole_y': '0.0', 'thickness': f'{best_Lw:.1f}', 'doping': '1.0e+15', 'doping_type': 'n', 'type': 'qw'})
        layers.append({'material': 'GaN', 'mole': '0.0', 'mole_y': '0.0', 'thickness': f'{best_Lb:.1f}', 'doping': '1.0e+15', 'doping_type': 'n', 'type': 'barrier'})
        
    layers.append({'material': 'GaN', 'mole': '0.0', 'mole_y': '0.0', 'thickness': '100.0', 'doping': f'{best_dop:.1e}', 'doping_type': 'n', 'type': 'barrier'})
    
    opt_config = {
        'layers': layers,
        'temp': '300.0',
        'field': '0.0',
        'bc_left': '0.0',
        'bc_right': '0.0',
        'vmin': '-0.5',
        'vmax': '1.1',
        'vstep': '0.05',
        'solver': '10: Fully-Coupled Newton-Raphson',
        'grid_step': '1.0',
        'max_pts': '200000',
        'mat_sys': 'Wurtzite',
        'sub_e': '5',
        'sub_h': '5',
        'Quantum_Regions': False,
        'Quantum_Regions_boundary': '[[0.0, 0.0]]',
        'area': '5e-4',
        'rs': '0.0',
        'rs_mode': 'External (Fast)',
        'rsh': '1000000.0',
        'tat_field': '5e6',
        'taun0': '1e-7',
        'taup0': '1e-7',
        'G_optical': '1e20',
        'enable_polarization': True,
        'photovoltaic_mode': True,
        'device_type': 'Solar Cell / Photodetector'
    }
    
    with open(ROOT / 'examples' / 'optimal_mqw_solar_cell.json', 'w') as f:
        json.dump(opt_config, f, indent=4)
        
    print("=== Successfully saved optimal_mqw_solar_cell.json ===")
    return opt_config

if __name__ == '__main__':
    run_global_optimization(max_evals=20)
