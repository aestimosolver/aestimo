import json
import sys
from pathlib import Path

# Add project root to path so we can import aestimo
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import aestimo

def main():
    print("Starting minimal debug run...")
    
    # Load base config from JSON
    json_path = ROOT / "examples" / "untitled_project.json"
    with open(json_path, 'r') as f:
        raw = json.load(f)
    
    # Convert layers to material matrix as required by aestimo
    mat_list = []
    for l in raw.get('layers', []):
        mat_list.append([
            float(l['thickness']), 
            l['material'], 
            float(l.get('mole', 0.0)), 
            float(l.get('mole_y', 0.0)), 
            float(l['doping']), 
            l['doping_type'], 
            l['type']
        ])

    # Create a minimal config for a quick test
    config = {
        # Structure
        "material": mat_list,
        "mat_type": str(raw.get("mat_sys", "Wurtzite")),
        "T": float(raw.get("temp", 300.0)),
        
        # Solver (minimal, fast settings)
        "comp_scheme": 7,  # SP-Drift Diffusion (Sequential) - as in JSON
        "gridfactor": 2.0, # Very coarse grid for speed
        "subnumber_e": 3,
        "subnumber_h": 3,
        
        # IV sweep (very short)
        "vmin": 0.0,
        "vmax": 0.1,  # Only go to 0.1V
        "Each_Step": 0.1,
        
        # Physics
        "taun0": 5e-7,
        "taup0": 5e-7,
        "G_optical": 0.0,  # Dark simulation for simplicity
        
        # Parasitics
        "device_area": float(raw.get("area", 5e-4)),
        "Rs": float(raw.get("rs", 11.7)),
        "Rsh": float(raw.get("rsh", 5000.0)),
        
        # Misc
        "enable_polarization": True,
        "Quantum_Regions": False,
        "photovoltaic_mode": True,
        "__file__": "debug_simple_run"  # Output directory name
    }
    
    print("Configuration loaded. Running aestimo...")
    try:
        # This should run quickly if everything is working
        aestimo.run_aestimo(config, drawFigures=False, show=False)
        print("SUCCESS: Simulation completed without stalling.")
    except Exception as e:
        print(f"ERROR: Simulation failed with exception: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()