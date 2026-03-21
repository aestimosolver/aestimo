import sys
import os
import numpy as np
import inspect

# Mocking parts of GUI for verification
class MockWidget:
    def __init__(self):
        self.children = []
    def destroy(self): pass
    def winfo_children(self): return self.children
    def pack(self, **kwargs): pass
    def grid(self, **kwargs): pass
    def grid_forget(self): pass
    def configure(self, **kwargs): pass

try:
    import aestimo_gui
    from aestimo_gui import AestimoGUI
    print("SUCCESS: aestimo_gui.py imported successfully.")
except Exception as e:
    print(f"FAILURE: aestimo_gui.py has errors: {e}")
    sys.exit(1)

# Test the analyze_iv_curve import
from characterize_solar import analyze_iv_curve
print("SUCCESS: analyze_iv_curve imported from characterize_solar.")

# Minimal verification of update_solar_metrics logic
# We won't run the full App (requires Display), but we can check if it initializes
# By actually running it in a sub-process or just checking the function exists.
if hasattr(AestimoGUI, 'update_solar_metrics'):
    print("SUCCESS: AestimoGUI has update_solar_metrics method.")
else:
    print("FAILURE: AestimoGUI missing update_solar_metrics method.")

# Check for Solar Study logic
if hasattr(AestimoGUI, 'start_solar_study_thread'):
    print("SUCCESS: AestimoGUI has start_solar_study_thread method.")
else:
    print("FAILURE: AestimoGUI missing start_solar_study_thread method.")

# Test Config Logic (Dry Run)
try:
    # We need a mock instance or similar. 
    # Since AestimoGUI.__init__ starts a lot of things, we might just check the methods exist and variable names in source.
    src_get = inspect.getsource(AestimoGUI.get_current_configuration)
    src_load = inspect.getsource(AestimoGUI.load_configuration)
    
    if "self.g_opt_entry.get()" in src_get:
        print("SUCCESS: get_current_configuration uses self.g_opt_entry.get()")
    else:
        print("FAILURE: get_current_configuration missing self.g_opt_entry.get() or using wrong name.")
        
    if "self.set_entry(self.g_opt_entry" in src_load:
        print("SUCCESS: load_configuration uses self.g_opt_entry")
    else:
        print("FAILURE: load_configuration missing self.g_opt_entry or using wrong name.")

    # Check for correct setup method
    src_setup_phys = inspect.getsource(AestimoGUI.setup_physics_tab)
    if 'setattr(self, "g_opt_entry", entry)' in inspect.getsource(AestimoGUI.create_input_row) or 'var_name' in inspect.getsource(AestimoGUI.create_input_row):
        # We know create_input_row does setattr. 
        # Let's check if setup_physics_tab sets the correct var_name
        if '"g_opt_entry"' in src_setup_phys:
             print("SUCCESS: setup_physics_tab initializes g_opt_entry")
        else:
             print("FAILURE: setup_physics_tab does not initialize g_opt_entry with correct name")
except Exception as e:
    print(f"VERIFICATION ERROR: {e}")
