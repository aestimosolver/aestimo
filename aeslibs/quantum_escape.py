# Quantum Well Carrier Escape Engine for Aestimo 1D
# Implements Thermionic Emission & WKB Quantum Mechanical Tunneling

import numpy as np

# Constants
H_BAR = 1.054571817e-34     # J s
Q_ELEM = 1.60217663e-19    # C
M_ELEM = 9.1093837015e-31  # kg
K_BOLTZ = 1.380649e-23     # J/K

def calculate_thermionic_emission_rate(barrier_height_ev, lw_m, eff_mass_rel, temp_k=300.0):
    """
    Calculates thermionic emission rate (s^-1) over heterojunction energy barrier.
    """
    if barrier_height_ev <= 0:
        return 1.0e13 # Instantaneous escape
        
    m_eff = eff_mass_rel * M_ELEM
    vt_j = K_BOLTZ * temp_k
    barrier_j = barrier_height_ev * Q_ELEM
    
    # Thermionic attempt frequency: nu_th = sqrt(kB T / (2 pi m* Lw^2))
    nu_th = np.sqrt(vt_j / (2.0 * np.pi * m_eff * (lw_m**2)))
    rate_th = nu_th * np.exp(-barrier_j / vt_j)
    
    return np.clip(rate_th, 1.0e3, 1.0e13)

def calculate_wkb_tunneling_rate(barrier_height_ev, lb_m, e_field_v_m, eff_mass_rel):
    """
    Calculates field-assisted WKB quantum mechanical tunneling rate (s^-1) through barrier.
    """
    if barrier_height_ev <= 0 or lb_m <= 0:
        return 1.0e13
        
    m_eff = eff_mass_rel * M_ELEM
    barrier_j = barrier_height_ev * Q_ELEM
    field_abs = max(abs(e_field_v_m), 1.0e4)
    
    # WKB transmission probability T_wkb = exp(-4 sqrt(2 m*) (Delta E)^(3/2) / (3 q hbar F))
    exponent = (4.0 * np.sqrt(2.0 * m_eff) * (barrier_j**1.5)) / (3.0 * Q_ELEM * H_BAR * field_abs)
    t_wkb = np.exp(-np.clip(exponent, 0.0, 80.0))
    
    # Attempt frequency nu_0 = hbar pi / (2 m* Lw^2)
    nu_0 = (H_BAR * np.pi) / (2.0 * m_eff * (3.0e-9**2))
    rate_tun = nu_0 * t_wkb
    
    return np.clip(rate_tun, 1.0e3, 1.0e13)

def compute_effective_escape_lifetime(barrier_height_ev, lw_m, lb_m, e_field_v_m, eff_mass_rel, temp_k=300.0):
    """
    Calculates total effective carrier escape lifetime tau_esc (s) from quantum well.
    1/tau_esc = 1/tau_th + 1/tau_tun
    """
    r_th = calculate_thermionic_emission_rate(barrier_height_ev, lw_m, eff_mass_rel, temp_k)
    r_tun = calculate_wkb_tunneling_rate(barrier_height_ev, lb_m, e_field_v_m, eff_mass_rel)
    
    r_total = r_th + r_tun
    tau_esc = 1.0 / max(r_total, 1.0e3)
    return tau_esc
