import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add current directory to path
sys.path.append(os.getcwd())
import aestimo
from aestimo import run_aestimo

# Configuration closely matching test_solar_cell.py but for a single point
# Using the corrected G_optical
G_OPTICAL = 5e18 # matches the current test_solar_cell.py

class InputObject:
    def __init__(self, d):
        for key, val in d.items():
            setattr(self, key, val)

def run_diagnostic():
    # Define structure: InGaN p-n junction
    # p-type (left) 200nm, n-type (right) 500nm
    gridfactor = 5.0 # nm
    thickness_p = 200
    thickness_n = 500
    
    N_p = int(thickness_p / gridfactor)
    N_n = int(thickness_n / gridfactor)
    
    # Material structure
    # (thickness, material, x, y)
    material = []
    material.append( (thickness_p, 'InGaN', 0.57, 0.0) )
    material.append( (thickness_n, 'InGaN', 0.57, 0.0) )
    
    # Doping profile (m^-3)
    # p-type: 1e19 cm^-3 = 1e25 m^-3 (negative)
    # n-type: 1e19 cm^-3 = 1e25 m^-3 (positive)
    dop_p = -1e19 * 1e6
    dop_n = 1e19 * 1e6
    
    dop_profile = np.concatenate([
        np.full(N_p, dop_p),
        np.full(N_n, dop_n)
    ])
    
    alloy_profile = np.full(len(material), 0.57)
    alloy_profile_y = np.full(len(material), 0.0)
    
    input_config = {
        'material': material,
        'alloy_profile': alloy_profile,
        'alloy_profile_y': alloy_profile_y,
        'substrate': 'GaN',
        'computation_scheme': 9, # Coupled Newton
        'gridfactor': gridfactor,
        'maxgridpoints': 200000,
        'enable_polarization': False, # Verify if this fixes potential inversion
        'mat_type': 'Wurtzite', # Correct parameter name expected by aestimo/InputObject?
        # Actually structure_from checks getattr(input_obj, "mat_crys_strc", "Zincblende")
        # But let's check what test_solar_cell does. 
        # test_solar_cell sets 'mat_type'.
        # If aestimo defaults to Zincblende, we need to be careful.
        # Let's try matching test_solar_cell exactly.
        'mat_crys_strc': 'Wurtzite',
        'subnumber_e': 5,
        'subnumber_h': 5,
        'T': 300.0,
        'Field': 0.0,
        'vmax': 0.1, # Run small step to force loop to enter
        'vmin': 0.0,
        'Each_Step': 0.1,
        'fval': 0.0,
        'tat_field': 1e12,
        'G_optical': G_OPTICAL, 
        'device_area': 1.0e-4,
        'dop_profile': dop_profile
    }
    
    print("Running diagnostic simulation at 0V...")
    aestimo.output_directory = os.path.join(os.getcwd(), "diagnostic_output")
    if not os.path.exists(aestimo.output_directory):
        os.makedirs(aestimo.output_directory)
        
    input_obj = InputObject(input_config)
    _, model, result, _ = run_aestimo(input_obj, drawFigures=False, show=False)
    
    print("\nExtracting physical parameters...")
    
    # Extract data from result object first, then model
    # For scheme 9, result is result_dd (AttrDict usually)
    
    # Extract data from FILES because result object is empty
    import glob
    
    # Potential file: potn_eh_0.00.dat
    pot_file = os.path.join(aestimo.output_directory, "potn_eh_0.00.dat")
    if os.path.exists(pot_file):
        print(f"Reading potential from {pot_file}")
        pot_data = np.loadtxt(pot_file)
        # Assuming format: x(nm)  fi(V)  (or similar)
        # Let's check shape
        if pot_data.shape[1] >= 2:
            x_plot = pot_data[:, 0] # nm
            fi_plot = pot_data[:, 1] # V
            # If there's a 3rd column, it might be something else, but index 1 is usually potential
            fi = fi_plot
        else:
            print(f"Error: pot file has unexpected shape {pot_data.shape}")
            return
    else:
        print(f"Error: Potential file {pot_file} not found")
        return

    # Carrier density file: np_data0_0.00.dat
    np_file = os.path.join(aestimo.output_directory, "np_data0_0.00.dat")
    if os.path.exists(np_file):
        print(f"Reading carrier densities from {np_file}")
        np_data = np.loadtxt(np_file)
        # Assuming format: x(nm)  n(1/m^3)  p(1/m^3)
        if np_data.shape[1] >= 3:
            n = np_data[:, 1]
            p = np_data[:, 2]
        else:
             print(f"Error: np file has unexpected shape {np_data.shape}")
             return
    else:
        print(f"Error: Carrier density file {np_file} not found")
        return

    # Set x axis for later use
    x = x_plot

    n = getattr(model, 'n', None)
    p = getattr(model, 'p', None)
    Ec = getattr(model, 'Ec', None)
    Ev = getattr(model, 'Ev', None)
    Fn = getattr(model, 'Fn', None) # Quasi-Fermi levels
    p = getattr(result, 'pf_result', None)
    
    # Check ni if available
    ni_arr = getattr(model, 'ni', None)
    if ni_arr is not None:
        print(f"Mean ni: {np.mean(ni_arr):.4e} m^-3")
    else:
        print("ni array not found in model")

    # Fallback to model if needed (sometimes stored there during computation)
    # If Ec/Ev are not stored, calculate them
    # Ec = -q*fi - chi + delta (?)
    # Aestimo usually calculates Ec_conduction_band in the plotting routines.
    # Let's use the potential profile to diagnose.
    
    x = np.arange(len(fi)) * gridfactor
    
    print(f"Nodes: {len(fi)}")
    print(f"Left (p-side) potential: {fi[0]:.4f} V")
    print(f"Right (n-side) potential: {fi[-1]:.4f} V")
    print(f"Built-in Potential Vbi: {fi[-1] - fi[0]:.4f} V")
    
    print(f"Left n: {n[0]:.4e} m^-3")
    print(f"Left p: {p[0]:.4e} m^-3")
    print(f"Right n: {n[-1]:.4e} m^-3")
    print(f"Right p: {p[-1]:.4e} m^-3")
    
    # Verify p-n junction logic
    # p-side: p >> n. fi should be lower.
    # n-side: n >> p. fi should be higher.
    
    p_side_check = p[0] > n[0]
    n_side_check = n[-1] > p[-1]
    potential_check = fi[-1] > fi[0]
    
    print(f"P-side majority carrier correct (p>n): {p_side_check}")
    print(f"N-side majority carrier correct (n>p): {n_side_check}")
    print(f"Potential slope correct (n > p): {potential_check}")
    
    # Check charge neutrality at edges
    # rho = p - n + dop
    dop_p_val = dop_profile[0]
    dop_n_val = dop_profile[-1]
    
    rho_left = p[0] - n[0] + dop_p_val
    rho_right = p[-1] - n[-1] + dop_n_val
    
    print(f"Charge neutrality left (should be ~0): {rho_left:.4e}")
    print(f"Charge neutrality right (should be ~0): {rho_right:.4e}")
    
    # Plotting
    plt.figure(figsize=(10,12))
    
    plt.subplot(3,1,1)
    plt.plot(x, fi, label='Potential (V)')
    plt.title('Electrostatic Potential')
    plt.grid(True)
    plt.legend()
    
    plt.subplot(3,1,2)
    plt.plot(x, np.log10(np.abs(n)+1e-30), 'b-', label='n')
    plt.plot(x, np.log10(np.abs(p)+1e-30), 'r-', label='p')
    plt.plot(x, np.log10(np.abs(dop_profile)+1e-30), 'g--', label='doping') 
    # Note: doping is signed, abs for log
    plt.title('Carrier Densities (log10 m^-3)')
    plt.grid(True)
    plt.legend()
    
    plt.subplot(3,1,3)
    # Band Diagram estimate
    # Ec = -fi - affinity (assume const affinity for now)
    # Ev = Ec - Eg
    # For qualitative check, just use -fi
    plt.plot(x, -fi, 'k-', label='-Potential (~Band Edge)')
    plt.title('Band Bending Qualitative')
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('diagnostic_plot.png')
    print("Diagnostic plot saved to diagnostic_plot.png")

if __name__ == "__main__":
    try:
        run_diagnostic()
    except Exception as e:
        print(f"Diagnostic failed: {e}")
        import traceback
        traceback.print_exc()

