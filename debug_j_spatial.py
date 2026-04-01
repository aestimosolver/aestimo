import numpy as np
import os
import sys

# Setup paths
root_dir = os.getcwd()
sys.path.append(root_dir)
import aestimo
from aestimo import run_aestimo

# Define a minimal solar cell input
class TinyInput:
    def __init__(self):
        self.T = 300.0
        self.computation_scheme = 9
        self.comp_scheme = 9
        self.subnumber_h = 5
        self.subnumber_e = 5
        self.gridfactor = 1.0
        self.dx = 1e-9
        self.maxgridpoints = 100
        self.mat_type = "Wurtzite"
        self.material = [
            [250.0, "InGaN", 0.57, 0.0, 1e16, "p", "barrier"],
            [250.0, "InGaN", 0.57, 0.0, 2e17, "n", "barrier"]
        ]
        self.vmin = 0.0
        self.vmax = 0.1
        self.Each_Step = 0.1
        self.G_optical = 5e18
        self.tat_field = 1e12
        self.device_area_m2 = 1e-8
        self.device_area = 0.0001
        self.inputfilename = "TINY_TEST"
        self.__file__ = "TINY_TEST.py"
        self.photovoltaic_mode = True
        self.enable_polarization = False
        self.work_function_left = 7.0
        self.work_function_right = 4.0
        self.surface_recomb = (0, 0)
        self.Quantum_Regions = False
        self.Quantum_Regions_boundary = np.zeros((1,2))
        self.Rs = 0.0
        
        # Build doping profile manually
        self.dx = 1e-9
        self.maxgridpoints = 500
        dop_arr = np.zeros(500)
        dop_arr[0:250] = -1e16 * 1e6 # m-3
        dop_arr[250:500] = 2e17 * 1e6 # m-3
        self.dop_profile = dop_arr

# Run simulation
print("Running tiny simulation...")
aestimo.output_directory = os.path.join(os.getcwd(), "TINY_TEST_output")
inp = TinyInput()
# We need to monkey-patch Current2 to print spatial data
import aestimo
orig_current2 = aestimo.Current2

def mock_current2(*args, **kwargs):
    results = orig_current2(*args, **kwargs)
    Jelec = results[2]
    Jhole = results[5]
    Jtotal = Jelec + Jhole
    vindex = args[0]
    print(f"DEBUG SPATIAL J at vindex {vindex}:")
    print(f"  Node 0: {Jtotal[vindex, 0]} A/m^2")
    print(f"  Node 100: {Jtotal[vindex, 100]} A/m^2")
    print(f"  Node 250 (junction): {Jtotal[vindex, 250]} A/m^2")
    print(f"  Node 400: {Jtotal[vindex, 400]} A/m^2")
    print(f"  Node 499: {Jtotal[vindex, 499]} A/m^2")
    return results

aestimo.Current2 = mock_current2

try:
    run_aestimo(inp, drawFigures=False, show=False)
except Exception as e:
    print(f"FAILED: {e}")
