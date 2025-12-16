import sys
import os
import shutil
import aestimo
from pathlib import Path

# We want to verify that run_aestimo changes the output directory 
# from default 'output' to '<filename>_output'

def test_output_naming():
    print("Starting Output Naming Verification...")
    
    # 1. Setup Mock Input
    # We need an input object that looks like a module or has __file__
    import numpy as np
    # 2. Reset global output_directory in aestimo to generic 'output'
    # This simulates the default state when aestimo.py is imported or run
    default_out = os.path.join(os.getcwd(), 'output')
    
    class MockInput:
        __file__ = os.path.abspath("test_project_mysim.py")
        T = 300
        computation_scheme = 2 # Schrodinger-Poisson
        subnumber_h = 1
        subnumber_e = 1
        Fapplied = 0
        vmax=0 # No sweep to be fast
        vmin=0
        Each_Step=1
        gridfactor=10 # Coarse grid for speed
        maxgridpoints=100
        mat_type='Zincblende'
        # A simple structure: GaAs well
        material = [[10, 'GaAs', 0, 0, 0, 'n', 'w']]
        
        # Physics params usually needed
        TAUN0 = 1e-9
        TAUP0 = 1e-9
        mun0 = 1000
        mup0 = 100
        Cn0 = 1e-30
        Cp0 = 1e-30
        BETAN = 1
        BETAP = 1
        VSATN = 1e7
        VSATP = 1e7
        
        # Doping profiles
        n_max = int(100) # matches maxgridpoints ideally? 
        # Actually logic inside will calc n_max from gridfactor/material width.
        # But if we provide dop_profile directly, it might use it?
        # Let's provide basic arrays.
        dop_n = np.zeros(100)
        dop_p = np.zeros(100)
        dop_profile = np.zeros(100)
        xaxis = np.arange(100)

        surface = [0, 0] # Left/Right BC
    
    # 2. Reset global output_directory in aestimo to generic 'output'
    # This simulates the default state when aestimo.py is imported or run
    default_out = os.path.join(os.getcwd(), 'output')
    aestimo.output_directory = default_out
    
    # PROBLEM: aestimo.py initializes logger on import (in the else block at end of file)
    # This opens aestimo.log in default_out, locking it.
    # We must close handlers to delete the directory.
    import logging
    logger = logging.getLogger('aestimo')
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)
    
    # Ensure expectation path doesn't exist yet or is clean
    expected_name = "test_project_mysim_output"
    expected_path = os.path.join(os.getcwd(), expected_name)
    
    # Now we can safely remove
    if os.path.exists(expected_path):
        try:
            shutil.rmtree(expected_path)
        except PermissionError:
            print(f"Warning: Could not remove {expected_path}, checking if we can write to it.")
            
    if os.path.exists(default_out):
        try:
            shutil.rmtree(default_out)
        except PermissionError:
             print(f"Warning: Could not remove {default_out}. Log file might still be locked.")

    # Re-init logger if needed, or let run_aestimo handle it?
    # run_aestimo doesn't init logger, it assumes it's there. 
    # But wait, run_aestimo doesn't check logger.
    # We should probably re-add a handler to stdout so we don't crash if it logs?
    # Actually, logging without handlers just outputs to stderr or nowhere (warnings).
    # Ideally, we want the logger to point to the NEW directory if it changes.
    # But currently aestimo.py doesn't re-init logger when output_directory changes in run_aestimo?
    # Let's check that logic later. For now, let's verify directory creation.
        
    print(f"Initial global output: {aestimo.output_directory}")
    print(f"Expected output after run: {expected_path}")
    
    # 3. Run Aestimo
    # We expect it to run, fail/succeed simulation doesn't matter much as long as it reaches the directory logic
    # The directory logic is AFTER simulation in run_aestimo (lines 4635+), so simulation must complete.
    # We provided minimal valid material, hopefully it runs fast.
    
    try:
        aestimo.run_aestimo(MockInput, drawFigures=False, show=False)
    except Exception as e:
        print(f"Simulation crashed: {e}")
        # Even if it crashed, check if it changed directory if the crash was late? 
        # Actually logic is at end, so it must succeed.
        # If it crashes, we need to fix MockInput to be valid enough.
    
    # 4. Verify
    print(f"Final global output: {aestimo.output_directory}")
    
    if aestimo.output_directory == expected_path:
        print("SUCCESS: Output directory was correctly renamed.")
        # Check if created
        if os.path.isdir(expected_path):
             print(f"SUCCESS: Directory '{expected_name}' was created.")
        else:
             print(f"WARNING: Directory path updated but folder not found on disk.")
    else:
        print(f"FAILURE: Output directory is {aestimo.output_directory}, expected {expected_path}")

if __name__ == "__main__":
    test_output_naming()
