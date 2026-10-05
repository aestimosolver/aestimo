import numpy as np
import matplotlib.pyplot as plt
import config
from aestimo import solve_device # Assuming solve_device is the main function in aestimo.py

def validate_tat_effect():
    # Test parameters
    test_voltage = -0.5  # Reverse bias point
    tat_fields = [1e10, 1e9, 1e8, 1e7] # From disabled to strong enhancement
    results = []

    print(f"Starting Reverse Transport Validation at V = {test_voltage}V")
    print(f"{'TAT Field (V/m)':<20} | {'Current (A/m^2)':<20}")
    print("-" * 42)

    for tf in tat_fields:
        # Override the model parameter
        # Since the model is often passed as an object or read from config,
        # we ensure the solver uses this specific tat_field.
        
        # This is a conceptual wrapper. In the actual aestimo.py, 
        # we need to ensure the 'model' object has the attribute.
        try:
            # We simulate a run for a single voltage point
            # result = solve_device(voltage=test_voltage, tat_field=tf)
            # current = result.av_curr
            
            # For the purpose of this validation script, we will 
            # call the solver and capture the output.
            # Since I don't have the exact function signature of solve_device,
            # I will use a placeholder that the user can adapt or I will 
            # refine once I read the full aestimo.py.
            
            # Placeholder for actual solver call
            current = 0.0 # Replace with actual solver call
            results.append((tf, current))
            print(f"{tf:<20.1e} | {current:<20.5e}")
        except Exception as e:
            print(f"Error testing tf={tf}: {e}")

    # Plotting the sensitivity
    plt.figure(figsize=(8, 6))
    plt.semilogy([r[0] for r in results], [r[1] for r in results], 'o-')
    plt.xscale('log')
    plt.xlabel('TAT Field Threshold (V/m)')
    plt.ylabel('Reverse Current Density (A/m^2)')
    plt.title(f'Reverse Leakage Sensitivity to TAT Field at {test_voltage}V')
    plt.grid(True, which="both", ls="-")
    plt.savefig('reverse_tat_sensitivity.png')
    print("\nSensitivity plot saved as 'reverse_tat_sensitivity.png'")

if __name__ == "__main__":
    validate_tat_effect()
```

Actually, to make this truly useful, I need to read the **entire** `aestimo.py` to see exactly how to call the solver for a single point. I will read `aestimo.py` now.

```tool
TOOL_NAME: read_file
BEGIN_ARG: filepath
"aestimo.py"