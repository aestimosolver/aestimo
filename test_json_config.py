"""
Quick test to verify test_solar_cell.json works with the simulation backend
"""
import json
import sys
import os

# Load the JSON config
with open('test_solar_cell.json', 'r') as f:
    config = json.load(f)

print("=" * 60)
print("TESTING JSON CONFIGURATION")
print("=" * 60)
print(f"\nLoaded configuration from test_solar_cell.json")
print(f"Solver: {config['solver']}")
print(f"Material: {config['layers'][0]['material']} (In content: {config['layers'][0]['mole']})")
print(f"G_optical: {config['G_optical']:.2e} cm^-3 s^-1")
print(f"Device area: {config['area']} cm^2")
print(f"Voltage range: {config['vmin']} to {config['vmax']} V (step: {config['vstep']})")

# Parse solver scheme
scheme_id = int(config["solver"].split(":")[0])
print(f"\nParsed computation scheme: {scheme_id}")

# Setup similar to GUI worker
import numpy as np
from aestimo import run_aestimo

class InputObject:
    def __init__(self, d):
        for key, val in d.items():
            setattr(self, key, val)

# Build doping profile from layers
gridfactor = config['grid_step']
layers = config['layers']
dop_arr = []
for layer in layers:
    thickness_nm = layer['thickness']
    N_layer = int(thickness_nm / gridfactor)
    doping_val = layer['doping']
    if layer['doping_type'] == 'p':
        doping_val = -abs(doping_val)
    else:
        doping_val = abs(doping_val)
    # Convert cm^-3 to m^-3
    doping_val_m3 = doping_val * 1e6
    dop_arr.extend([doping_val_m3] * N_layer)

dop_arr = np.array(dop_arr)
print(f"\nDoping profile: {len(dop_arr)} points")
print(f"  p-region: {dop_arr[0]:.2e} m^-3")
print(f"  n-region: {dop_arr[-1]:.2e} m^-3")

# Build material structure in the format StructureFrom expects
# Each element should be: (thickness_nm, material_name, x_alloy, y_alloy)
material_structure = []

for layer in layers:
    thickness_nm = layer['thickness']
    material_name = layer['material']
    x_alloy = layer['mole']
    y_alloy = layer.get('mole_y', 0.0)
    material_structure.append((thickness_nm, material_name, x_alloy, y_alloy))

# Create input configuration
input_config = {
    'material': material_structure,
    'substrate': 'GaN',
    'computation_scheme': scheme_id,
    'gridfactor': gridfactor,
    'maxgridpoints': config.get('max_pts', 200000),
    'maxgridpoints': config.get('max_pts', 200000),
    'mat_type': config.get('mat_type', 'Wurtzite'),
    'subnumber_e': config.get('sub_e', 5),
    'subnumber_h': config.get('sub_h', 5),
    'T': config.get('temp', 300.0),
    'Field': config.get('field', 0.0),
    'vmax': config['vmax'],
    'vmin': config['vmin'],
    'Each_Step': config.get('vstep', 0.05),
    'fval': config.get('bc_left', 0.0),
    'tat_field': config.get('tat_field', 1e12),
    'G_optical': config.get('G_optical', 0.0),
    'device_area': config.get('area', 1e-4),
    'dop_profile': dop_arr
}

print(f"\nStarting simulation with scheme {scheme_id}...")
print("=" * 60)

try:
    import aestimo
    aestimo.output_directory = os.path.join(os.getcwd(), "examples", "test_json_output")
    os.makedirs(aestimo.output_directory, exist_ok=True)
    
    input_obj = InputObject(input_config)
    print(f"DEBUG: InputObject.mat_crys_strc = '{getattr(input_obj, 'mat_crys_strc', 'NOT SET')}'")
    _, model, result, figures = run_aestimo(input_obj, drawFigures=False, show=False)
    
    print("\n" + "=" * 60)
    print("[SUCCESS] SIMULATION COMPLETED")
    print("=" * 60)
    print(f"\nOutput directory: {aestimo.output_directory}")
    
    # Check for IV data
    iv_file = os.path.join(aestimo.output_directory, 'av_curr.dat')
    if os.path.exists(iv_file):
        data = np.loadtxt(iv_file)
        print(f"\nGenerated {len(data)} voltage points")
        print(f"Voltage range: {data[0,0]:.3f} to {data[-1,0]:.3f} V")
        print(f"Current range: {data[:,1].min():.3e} to {data[:,1].max():.3e} A/m²")
    
    print("\n[SUCCESS] JSON configuration is valid and simulation runs successfully!")
    
except Exception as e:
    print("\n" + "=" * 60)
    print("[ERROR] SIMULATION FAILED")
    print("=" * 60)
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
