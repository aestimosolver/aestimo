import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

# --- Experimental Data from CALIBRATION_STATUS.md ---
# Voltage (V), Current (A)
EXPERIMENTAL_DATA = {
    0.2: 2.2e-9,
    0.4: 49.4e-9,
    0.6: 1.29e-6,
    0.8: 38.6e-6
}
RS_VALUE = 5.0  # Series Resistance in Ohms

def load_latest_results():
    """
    Searches for the most recent current table in the results directory.
    Looking for files named 'raw_current_table.csv' or similar.
    """
    search_pattern = os.path.join("results", "**", "raw_current_table.csv")
    files = glob.glob(search_pattern, recursive=True)
    if not files:
        # Fallback: check for any .dat or .csv in results that might contain IV data
        search_pattern = os.path.join("results", "**", "*.csv")
        files = glob.glob(search_pattern, recursive=True)
        
    if not files:
        raise FileNotFoundError("No simulation result files found in results/ directory.")
    
    # Get the most recently modified file
    latest_file = max(files, key=os.path.getmtime)
    print(f"Loading latest results from: {latest_file}")
    return pd.read_csv(latest_file)

def apply_rs_correction(voltage, current, rs):
    """
    Calculates the effective voltage at the device junction.
    V_junction = V_applied - I * Rs
    """
    return voltage - (current * rs)

def analyze():
    try:
        df = load_latest_results()
    except FileNotFoundError as e:
        print(e)
        return

    # Assuming the CSV has columns 'Voltage' and 'Current' (or similar)
    # We will normalize column names to lowercase for robustness
    df.columns = [c.lower() for c in df.columns]
    v_col = 'voltage' if 'voltage' in df.columns else df.columns[0]
    i_col = 'current' if 'current' in df.columns else df.columns[1]

    # 1. Extract simulated values at experimental voltages
    sim_results = []
    exp_v = sorted(EXPERIMENTAL_DATA.keys())
    exp_i = [EXPERIMENTAL_DATA[v] for v in exp_v]

    print("\n--- Comparison Analysis ---")
    print(f"{'Voltage (V)':<12} | {'Exp Current (A)':<18} | {'Sim Current (A)':<18} | {'Ratio (Sim/Exp)':<15}")
    print("-" * 65)

    for v in exp_v:
        # Find the closest simulated voltage point
        idx = (df[v_col] - v).abs().idxmin()
        sim_i = df[i_col].iloc[idx]
        ratio = sim_i / EXPERIMENTAL_DATA[v]
        sim_results.append(sim_i)
        print(f"{v:<12.2f} | {EXPERIMENTAL_DATA[v]:<18.2e} | {sim_i:<18.2e} | {ratio:<15.2f}")

    # 2. Series Resistance Correction
    # We calculate what the current would be if we shifted the voltage by I*Rs
    # Since we have a full IV curve, we can interpolate the current at (V_applied - I*Rs)
    corrected_results = []
    for v in exp_v:
        # Iterative approach to find I such that I = f(V_applied - I*Rs)
        # Start with the raw simulated current
        idx = (df[v_col] - v).abs().idxmin()
        current_guess = df[i_col].iloc[idx]
        
        for _ in range(5): # 5 iterations usually converge for small Rs
            v_eff = v - (current_guess * RS_VALUE)
            # Interpolate current from the simulated curve at v_eff
            current_guess = np.interp(v_eff, df[v_col], df[i_col])
        
        corrected_results.append(current_guess)

    print("\n--- With Rs = 5 Ohm Correction ---")
    print(f"{'Voltage (V)':<12} | {'Exp Current (A)':<18} | {'Corr Current (A)':<18} | {'Ratio (Corr/Exp)':<15}")
    print("-" * 65)
    for i, v in enumerate(exp_v):
        ratio = corrected_results[i] / EXPERIMENTAL_DATA[v]
        print(f"{v:<12.2f} | {EXPERIMENTAL_DATA[v]:<18.2e} | {corrected_results[i]:<18.2e} | {ratio:<15.2f}")

    # 3. Plotting
    plt.figure(figsize=(10, 6))
    plt.semilogy(exp_v, exp_i, 'ko', label='Experimental', markersize=8)
    plt.semilogy(exp_v, sim_results, 'b-o', label='Simulated (Baseline)', linewidth=2)
    plt.semilogy(exp_v, corrected_results, 'r-o', label=f'Simulated + Rs({RS_VALUE}Ω)', linewidth=2)
    
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current (A)')
    plt.title('IV Curve Calibration: Baseline vs Experimental')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    plt.savefig('results/calibration_comparison_plot.png')
    print("\nPlot saved to results/calibration_comparison_plot.png")

if __name__ == "__main__":
    analyze()