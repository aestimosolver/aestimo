import numpy as np
import os
import sys

# Add current dir to path
sys.path.append(os.getcwd())

import aestimo
from aeslibs.experimental_validation import (
    load_experimental_data,
    calculate_current_from_simulation,
    compute_error_metrics
)

# Load the example script as a module
import importlib.util
spec = importlib.util.spec_from_file_location("sample_pn", "sample_pn_with_experimental_validation.py")
sample_pn = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sample_pn)

# LIGTHWEIGHT OVERRIDES for FAST EXIT
sample_pn.gridfactor = 50   # Broad grid for speed
sample_pn.vmin = 0.0
sample_pn.vmax = 0.2
sample_pn.Each_Step = 0.1   # Only 3 steps (0.0, 0.1, 0.2)
sample_pn.inputfilename = "fast_verify_output"

input_obj = vars(sample_pn)

# Run simulation
print("Running LIGHTWEIGHT simulation to verify sparse refactor...")
results = aestimo.run_aestimo(input_obj, drawFigures=False, show=False)

print("\nFAST EXIT SUCCESSFUL: Simulation completed in seconds with new sparse solvers.")
