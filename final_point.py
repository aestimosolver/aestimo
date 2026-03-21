#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import sys

# Setup paths
root_dir = os.getcwd()
sys.path.append(root_dir)

from aestimo import run_aestimo
from characterize_solar import analyze_iv_curve

def run_temp_point(T_kelvin):
    class InputObject:
        def __init__(self, dictionary):
            for key, value in dictionary.items():
                setattr(self, key, value)
            self.__file__ = "study_temp_efficiency.py"
            self.inputfilename = f"test_solar_T_{int(T_kelvin)}"
                
    config_dict = {
        'material': [
            [250.0, 'InGaN', 0.57, 0.0, 1e16, 'p', 'b'],
            [250.0, 'InGaN', 0.57, 0.0, 2e17, 'n', 'b']
        ],
        'T': T_kelvin,
        'gridfactor': 1.0,
        'maxgridpoints': 200000,
        'mat_type': 'Wurtzite',
        'computation_scheme': 8, 
        'vmin': 0.0,
        'vmax': 1.0, 
        'Each_Step': 0.05,
        'fval': 0.0, 
        'tat_field': 1e12,
        'G_optical': 1e22, 
        'device_area': 1.0e-4 
    }
    
    gf = 1.0
    N_p = int(250.0 / gf)
    N_n = int(250.0 / gf)
    dop_p = -1e16 * 1e6
    dop_n =  2e17 * 1e6
    dop_arr = np.zeros(N_p + N_n)
    dop_arr[:N_p] = dop_p
    dop_arr[N_p:] = dop_n
    config_dict['dop_profile'] = dop_arr
    
    config_obj = InputObject(config_dict)
    _, model, result_ps, figures = run_aestimo(config_obj, drawFigures=False, show=False)
    
    import aestimo
    out_dir = aestimo.output_directory
    data = np.loadtxt(os.path.join(out_dir, 'av_curr.dat'))
    voltages = data[:, 0]
    j_density_am2 = data[:, 1]
    
    j_density_acm2 = j_density_am2 * 1e-4
    total_current_a = j_density_acm2 * config_dict['device_area']
    
    metrics = analyze_iv_curve(voltages, total_current_a, area_cm2=config_dict['device_area'], pin_mw_cm2=100.0)
    return metrics

m = run_temp_point(250.0)
print(f"{250:<10} | {m['jsc']:<15.4f} | {m['voc']:<10.4f} | {m['ff']:<10.2f} | {m['eta']:<15.4f}")
