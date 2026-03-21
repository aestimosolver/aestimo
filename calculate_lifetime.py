"""
Calculate required SRH lifetime to match experimental Js

For a Si p-n junction with Na = Nd = 1e18 cm^-3:
- Experimental Js ≈ 1e-12 A/cm²
- Theoretical: Js = q * (Dp*pn0/Lp + Dn*np0/Ln)
- Where Lp = sqrt(Dp * tau_p), Ln = sqrt(Dn * tau_n)

Current simulation gives ~143x too high current at 0.4V
This suggests we need shorter lifetimes to increase recombination

Strategy:
1. Reduce TAUN0 and TAUP0 by factor of ~10-100
2. This will reduce diffusion lengths and thus Js
"""

import numpy as np

# Constants
q = 1.602e-19  # C
kb = 1.38e-23  # J/K
T = 300  # K
Vt = kb * T / q  # V

# Si parameters at 300K
ni = 1e10  # cm^-3 (intrinsic carrier concentration)
Na = Nd = 1e18  # cm^-3 (doping)

# Minority carrier concentrations at equilibrium
pn0 = ni**2 / Nd  # holes in n-side
np0 = ni**2 / Na  # electrons in p-side

print(f"Minority carriers: pn0 = {pn0:.2e} cm^-3, np0 = {np0:.2e} cm^-3")

# Mobility from Caughey-Thomas for 1e18 cm^-3
# (These match what we implemented)
mu_n = 92 + (1360 - 92) / (1 + (1e18 / 1.3e17)**0.91)  # cm²/Vs
mu_p = 47.7 + (460 - 47.7) / (1 + (1e18 / 6.3e16)**0.76)  # cm²/Vs

print(f"Mobilities: mu_n = {mu_n:.1f} cm²/Vs, mu_p = {mu_p:.1f} cm²/Vs")

# Diffusion coefficients (Einstein relation)
Dn = mu_n * Vt  # cm²/s
Dp = mu_p * Vt  # cm²/s

print(f"Diffusion: Dn = {Dn:.2f} cm²/s, Dp = {Dp:.2f} cm²/s")

# Target Js
Js_target = 1e-12  # A/cm²

# Current lifetime (simulation default)
tau_current = 1e-7  # s

# Calculate Ln and Lp with current lifetime
Ln_current = np.sqrt(Dn * tau_current)
Lp_current = np.sqrt(Dp * tau_current)

print(f"\nCurrent lifetime: {tau_current:.2e} s")
print(f"Current diffusion lengths: Ln = {Ln_current*1e4:.1f} µm, Lp = {Lp_current*1e4:.1f} µm")

# Calculate current Js
Js_current = q * (Dp * pn0 / Lp_current + Dn * np0 / Ln_current)
print(f"Current Js (theoretical): {Js_current:.2e} A/cm²")

# Required reduction factor
reduction_factor = Js_current / Js_target
print(f"\nReduction factor needed: {reduction_factor:.1f}x")

# Since Js ∝ 1/sqrt(tau), we need tau_new = tau_current / reduction_factor²
tau_new = tau_current / (reduction_factor**2)
print(f"Required lifetime: {tau_new:.2e} s")

# Verify
Ln_new = np.sqrt(Dn * tau_new)
Lp_new = np.sqrt(Dp * tau_new)
Js_new = q * (Dp * pn0 / Lp_new + Dn * np0 / Ln_new)

print(f"\nWith new lifetime:")
print(f"  Diffusion lengths: Ln = {Ln_new*1e4:.1f} µm, Lp = {Lp_new*1e4:.1f} µm")
print(f"  Js = {Js_new:.2e} A/cm²")
print(f"  Match factor: {Js_new/Js_target:.2f}x")

print(f"\n=== RECOMMENDED VALUES ===")
print(f"TAUN0 = {tau_new:.2e} s  (currently: {tau_current:.2e} s)")
print(f"TAUP0 = {tau_new:.2e} s  (currently: {tau_current:.2e} s)")
