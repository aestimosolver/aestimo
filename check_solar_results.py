import numpy as np
import matplotlib.pyplot as plt
import os

def check_results():
    out_dir = "test_solar_cell_output"
    iv_file = os.path.join(out_dir, "av_curr.dat")
    
    if not os.path.exists(iv_file):
        print(f"Error: {iv_file} not found.")
        return

    data = np.loadtxt(iv_file)
    v = data[:, 0]
    j_raw = data[:, 1] # A/m^2
    j = j_raw * 0.1 # mA/cm^2
    
    plt.figure(figsize=(8, 6))
    plt.plot(v, j, 'o-')
    plt.axhline(0, color='k', lw=1)
    plt.axvline(0, color='k', lw=1)
    plt.xlabel("Voltage (V)")
    plt.ylabel("Current (mA/cm²)")
    plt.title("J-V Curve - InGaN Solar Cell (High Injection)")
    plt.grid(True)
    plt.savefig("diagnostic_iv_plot.png")
    print("Plot saved as diagnostic_iv_plot.png")

    # Metrics
    jsc = -np.interp(0, v, j)
    voc_candidates = []
    if np.min(j) < 0 and np.max(j) > 0:
        voc = np.interp(0, j, v)
        print(f"Jsc: {jsc:.4f} mA/cm²")
        print(f"Voc: {voc:.4f} V")
    else:
        print("Warning: No VOC found in the simulated range.")
        print(f"J range: {np.min(j):.4g} to {np.max(j):.4g}")

if __name__ == "__main__":
    check_results()
