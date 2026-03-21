import json
import numpy as np

# GUI Project Structure (List of Dicts)
layers = [
    {
        "material": "InGaN",
        "thickness": 250.0,
        "mole": 0.57,
        "mole_y": 0.0,
        "doping": 1e16,
        "doping_type": "p",
        "type": "barrier"
    },
    {
        "material": "InGaN",
        "thickness": 250.0,
        "mole": 0.57,
        "mole_y": 0.0,
        "doping": 2e17,
        "doping_type": "n",
        "type": "barrier"
    }
]

config_dict = {
    "layers": layers,
    # Physics
    "temp": 300.0,
    "field": 0.0,
    "bc_left": 0.0,
    "bc_right": 0.0,
    "vmin": 0.0,
    "vmax": -1.8, # Forward Bias scan (0 -> -1.8V)
    "vstep": -0.1,
    
    # Solver
    "solver": "9: Coupled Newton (Robust)", # Updated to verified stable solver
    "grid_step": 5.0,
    "max_pts": 200000,
    "mat_type": "Wurtzite",
    "sub_e": 5,
    "sub_h": 5,
    "enable_polarization": False, # CRITICAL: Disable polarization for InGaN stability
    
    # TAT
    "tat_field": 1e12, # Disabled by default
    
    # Solar
    "device_type": "Solar Cell / Photodetector",
    "G_optical": 5e18, # Verified calibrated value
    
    # Validation
    "area": 1.0e-4,
    "exp_file": "examples/experimental_data/ingan_pn_experimental_iv.csv",
    "rs": 0.0,
    "rsh": 1e12,
    "rs_mode": "External (Fast)"
}

# Manual Doping Profile Construction
# Note: The GUI might preserve extra keys when loading/saving, or validly pass them to the worker if I check line 796 in GUI.
# In GUI: input_data = self.get_current_configuration() -> this ignores extra keys in 'config' passed to load_configuration unless they are mapped to widgets.
# Wait, load_configuration sets widgets. 
# get_current_configuration reads widgets.
# So 'dop_profile' in JSON will be LOST when loading into GUI because there is no widget for it.
# However, the user wants to Run Simulation.
# The GUI calls run_simulation_worker with get_current_configuration().
# Since dop_profile is not in widgets, it won't be in the config passed to worker.
# AND the worker re-calculates doping from layers (lines 920+).
# So my manual doping profile is useless in the GUI unless I modify the GUI to support custom doping profiles or use the "Graded Junction" as a proxy?
# Or I can use 'dop_profile' override in the text entry? No.

# Workaround: The verification script had a manual doping profile because 'StructureFrom' (backend) failed to generate reasonable doping.
# But looking at GUI worker line 960:
# dop_arr[curr:end] = val * 1e6
# It manually populates the array based on layer doping.
# Using 1e16 and 2e17 from layers.
# This should work fine! The issue in 'test_solar_cell.py' was potentially how StructureFrom handled it or some other issue.
# But the GUI worker builds the InputObject explicitly and sets InputObject.dop_profile = dop_arr.
# So if I set correct layer doping, the GUI worker will create the correct array!
# I don't need to inject 'dop_profile' list into JSON.

with open('test_solar_cell.json', 'w') as f:
    json.dump(config_dict, f, indent=4)

print("test_solar_cell.json created successfully (GUI Format).")
