
import numpy as np
import os
import sys

# Ensure local modules are found
sys.path.append(os.getcwd())

from aeslibs.experimental_validation import (
    load_experimental_data, plot_iv_comparison
)

def regenerate():
    print("Regenerating IV comparison plot...")

    # 1. Load Experimental Data
    exp_file = "examples/experimental_data/si_pn_experimental_iv.csv"
    if not os.path.exists(exp_file):
        print(f"Error: Experimental file not found: {exp_file}")
        return
    
    exp_voltage, exp_current = load_experimental_data(exp_file)
    print(f"Loaded experimental data: {len(exp_voltage)} points")

    # 2. Load Simulation Data
    # Priority: av_curr_verify.dat (from verify_fix.py), then output folder
    sim_files = [
        "av_curr_verify.dat",
        "sample_pn_with_experimental_validation_output/av_curr.dat"
    ]
    
    sim_file = None
    for f in sim_files:
        if os.path.exists(f):
            # Check if file has content
            if os.path.getsize(f) > 0:
                sim_file = f
                break
    
    if not sim_file:
        print("Error: No simulation current file found (av_curr.dat). Run simulation first.")
        return

    print(f"Using simulation data from: {sim_file}")
    
    try:
        sim_data = np.loadtxt(sim_file)
        # Handle case where file might be empty or 1 line
        if sim_data.ndim == 1:
             sim_data = sim_data.reshape(1, -1)
             
        sim_voltage = sim_data[:, 0]
        sim_current = sim_data[:, 1]
    except Exception as e:
        print(f"Error loading simulation data: {e}")
        return

    # 3. Plot
    output_dir = "sample_pn_with_experimental_validation_output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output_path = os.path.join(output_dir, "iv_comparison.png")
    
    plot_iv_comparison(exp_voltage, exp_current, sim_voltage, sim_current, 
                       output_path=output_path, show=False)
    
    print("Done!")

if __name__ == "__main__":
    regenerate()
