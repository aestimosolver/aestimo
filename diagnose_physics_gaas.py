import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add current directory to path
sys.path.append(os.getcwd())
import aestimo
from aestimo import run_aestimo

# Configuration closely matching test_solar_cell.py but for a single point
# Using GaAs (Zincblende)
G_OPTICAL = 5e18 

class InputObject:
    def __init__(self, d):
        for key, val in d.items():
            setattr(self, key, val)

def run_diagnostic():
    # Define structure: GaAs p-n junction
    gridfactor = 5.0 # nm
    thickness_p = 200
    thickness_n = 500
    
    N_p = int(thickness_p / gridfactor)
    N_n = int(thickness_n / gridfactor)
    
    # Material structure
    # GaAs
    material = []
    material.append( (thickness_p, 'GaAs', 0.0, 0.0) )
    material.append( (thickness_n, 'GaAs', 0.0, 0.0) )
    
    # Doping profile (m^-3)
    # Start with 1e18 cm^-3 (standard)
    dop_p = -1e18 * 1e6
    dop_n = 1e18 * 1e6
    
    dop_profile = np.concatenate([
        np.full(N_p, dop_p),
        np.full(N_n, dop_n)
    ])
    
    alloy_profile = np.full(len(material), 0.0)
    alloy_profile_y = np.full(len(material), 0.0)
    
    input_config = {
        'material': material,
        'alloy_profile': alloy_profile,
        'alloy_profile_y': alloy_profile_y,
        'substrate': 'GaAs',
        'computation_scheme': 9, # Coupled Newton
        'gridfactor': gridfactor,
        'maxgridpoints': 200000,
        'mat_type': 'Zincblende',
        'mat_crys_strc': 'Zincblende',
        'subnumber_e': 5,
        'subnumber_h': 5,
        'T': 300.0,
        'Field': 0.0,
        'vmax': 0.1, 
        'vmin': 0.0,
        'Each_Step': 0.1,
        'fval': 0.0,
        'tat_field': 1e12,
        'G_optical': G_OPTICAL, 
        'device_area': 1.0e-4,
        'dop_profile': dop_profile
    }
    
    print("Running GaAs diagnostic simulation at 0V...")
    aestimo.output_directory = os.path.join(os.getcwd(), "diagnostic_gaas_output")
    if not os.path.exists(aestimo.output_directory):
        os.makedirs(aestimo.output_directory)
        
    input_obj = InputObject(input_config)
    _, model, result, _ = run_aestimo(input_obj, drawFigures=False, show=False)
    
    print("\nExtracting physical parameters...")
    
    # Extract data from FILES
    import glob
    
    # Potential file: potn_eh_0.00.dat
    pot_file = os.path.join(aestimo.output_directory, "potn_eh_0.00.dat")
    if os.path.exists(pot_file):
        print(f"Reading potential from {pot_file}")
        pot_data = np.loadtxt(pot_file)
        if pot_data.shape[1] >= 2:
            x_plot = pot_data[:, 0] # nm
            fi_plot = pot_data[:, 1] # V
            fi = fi_plot
            print(f"Nodes: {len(fi)}")
            print(f"Left (p-side) potential: {fi[0]:.4f} V")
            print(f"Right (n-side) potential: {fi[-1]:.4f} V")
        else:
            print(f"Error: pot file has unexpected shape {pot_data.shape}")
            return
    else:
        print(f"Error: Potential file {pot_file} not found")
        return

    # Carrier density file
    np_file = os.path.join(aestimo.output_directory, "np_data0_0.00.dat")
    if os.path.exists(np_file):
        print(f"Reading carrier densities from {np_file}")
        np_data = np.loadtxt(np_file)
        if np_data.shape[1] >= 3:
            n = np_data[:, 1]
            p = np_data[:, 2]
            print(f"Left n: {n[0]:.4e} m^-3")
            print(f"Left p: {p[0]:.4e} m^-3")
            print(f"Right n: {n[-1]:.4e} m^-3")
            print(f"Right p: {p[-1]:.4e} m^-3")
        else:
             print(f"Error: np file has unexpected shape {np_data.shape}")
             return
    else:
        print(f"Error: Carrier density file {np_file} not found")
        return

if __name__ == "__main__":
    try:
        run_diagnostic()
    except Exception as e:
        print(f"Diagnostic failed: {e}")
        import traceback
        traceback.print_exc()
