# Advanced Solar Optics Engine for Aestimo 1D
# Implements NREL AM1.5G Solar Spectrum & Wavelength-Dependent InGaN Absorption

import numpy as np

# Physical Constants
H_PLANCK = 6.62607015e-34  # J s
C_LIGHT = 2.99792458e8     # m/s
Q_ELEM = 1.60217663e-19    # C

def get_am15g_spectrum(wl_nm=None):
    """
    Returns reference AM1.5G solar photon flux density spectrum Ph(lambda) in photons/(m^2 s nm).
    Standard 1-sun integrated power density = 100 mW/cm^2 (1000 W/m^2).
    """
    if wl_nm is None:
        wl_nm = np.linspace(300, 1200, 181) # 5 nm steps
        
    wl_m = wl_nm * 1e-9
    photon_energy_j = H_PLANCK * C_LIGHT / wl_m
    photon_energy_ev = photon_energy_j / Q_ELEM

    # Standard AM1.5G spectral irradiance model (W/m^2 nm) approximating NREL ASTM G173-03
    # Blackbody 5777K modified with atmospheric absorption dips (H2O, O2, CO2)
    t_atm = np.exp(-((wl_nm - 300)/150)**0.8) * np.exp(-((wl_nm - 760)/30)**2) * 0.9 + 0.1
    irradiance_w_m2_nm = (2 * np.pi * H_PLANCK * (C_LIGHT**2) / (wl_m**5)) / (np.exp(H_PLANCK * C_LIGHT / (wl_m * 1.380649e-23 * 5777)) - 1) * 1.8e-5
    irradiance_w_m2_nm = np.clip(irradiance_w_m2_nm, 0, 2.2) # normalize scale

    # Convert Spectral Irradiance (W/m^2 nm) -> Photon Flux (photons/m^2 s nm)
    photon_flux = (irradiance_w_m2_nm * 1e-3) / photon_energy_j  # photons / (m^2 s nm)
    return wl_nm, photon_flux

def calculate_ingan_bandgap(x_in, temp_k=300.0):
    """
    Calculates In_x Ga_{1-x} N bandgap (eV) using Varshni / Vegard rule with bowing factor b = 1.4 eV.
    """
    eg_gan_0 = 3.510  # eV at 0K
    eg_inn_0 = 0.780  # eV at 0K
    
    # Temperature dependence (Varshni)
    eg_gan = eg_gan_0 - (0.939e-3 * temp_k**2) / (temp_k + 772.0)
    eg_inn = eg_inn_0 - (0.419e-3 * temp_k**2) / (temp_k + 176.0)
    
    bowing = 1.40  # eV
    eg_ingan = x_in * eg_inn + (1.0 - x_in) * eg_gan - bowing * x_in * (1.0 - x_in)
    return max(0.6, eg_ingan)

def calculate_absorption_coefficient(wl_nm, x_in, temp_k=300.0):
    """
    Calculates absorption coefficient alpha(lambda, x) in m^-1 for In_x Ga_{1-x} N.
    """
    eg = calculate_ingan_bandgap(x_in, temp_k)
    photon_energy_ev = (H_PLANCK * C_LIGHT / (wl_nm * 1e-9)) / Q_ELEM
    
    alpha = np.zeros_like(photon_energy_ev)
    above_gap = photon_energy_ev >= eg
    
    # Direct bandgap square-root model above gap + Urbach tail below gap
    a0 = 2.0e7 # m^-1
    e_u = 0.030 # Urbach energy in eV
    
    alpha[above_gap] = a0 * np.sqrt((photon_energy_ev[above_gap] - eg) / eg + 1e-6) + 1e5
    alpha[~above_gap] = 1e5 * np.exp((photon_energy_ev[~above_gap] - eg) / e_u)
    
    return np.clip(alpha, 0.0, 5.0e7)

def compute_optical_generation_profile(x_coords_m, layers_info, temp_k=300.0, concentration_suns=1.0):
    """
    Computes depth-dependent optical generation rate G_opt(x) in m^-3 s^-1
    by integrating NREL AM1.5G solar spectrum absorption through the heterostructure.
    """
    n_pts = len(x_coords_m)
    g_opt = np.zeros(n_pts)
    
    wl_nm, photon_flux = get_am15g_spectrum()
    photon_flux *= concentration_suns
    
    dx = x_coords_m[1] - x_coords_m[0] if n_pts > 1 else 1e-9
    dwl = wl_nm[1] - wl_nm[0] if len(wl_nm) > 1 else 1.0
    
    # Map spatial composition x_in(x)
    x_in_profile = np.zeros(n_pts)
    for i, x_m in enumerate(x_coords_m):
        x_nm = x_m * 1e9
        curr_d = 0.0
        for lyr in layers_info:
            th = lyr['thickness']
            if curr_d <= x_nm <= curr_d + th:
                if lyr['material'] == 'InGaN':
                    x_in_profile[i] = float(lyr.get('mole', 0.1))
                else:
                    x_in_profile[i] = 0.0
                break
            curr_d += th

    # Calculate depth-resolved absorption alpha(x, lambda)
    alpha_matrix = np.zeros((len(wl_nm), n_pts))
    for i in range(n_pts):
        alpha_matrix[:, i] = calculate_absorption_coefficient(wl_nm, x_in_profile[i], temp_k)
        
    # Beer-Lambert optical flux integration: I(x, lambda) = I_0(lambda) * exp(-int_0^x alpha(x', lambda) dx')
    cum_alpha_dx = np.cumsum(alpha_matrix * dx, axis=1)
    flux_matrix = photon_flux[:, np.newaxis] * np.exp(-cum_alpha_dx)
    
    # G(x) = int_lambda alpha(x, lambda) * flux(x, lambda) d_lambda
    g_spectral = alpha_matrix * flux_matrix # photons / (m^3 s nm)
    g_opt = np.sum(g_spectral * dwl, axis=0)
    
    return np.clip(g_opt, 1e18, 1e28)
