# -*- coding: utf-8 -*-
"""
Aestimo 1D LED Characterization and Physics Engine
Calculates spatial recombination rates, internal quantum efficiency (IQE),
external quantum efficiency (EQE), optical output power (P_opt), wall-plug efficiency (WPE),
efficiency droop profiles, and electroluminescence (EL) emission spectra.
"""

from __future__ import annotations
import numpy as np
import scipy.constants as const
from scipy.ndimage import gaussian_filter1d

# Physical constants in SI
Q_ELEC = const.e                       # 1.602176634e-19 C
H_PLANCK = const.h                     # 6.62607015e-34 J*s
C_LIGHT = const.c                      # 2.99792458e8 m/s
K_BOLTZ = const.k                      # 1.380649e-23 J/K
HC_EV_NM = 1239.841984                 # hc in eV*nm


def calculate_led_recombination_profile(
    carrier_n_cm3: np.ndarray,
    carrier_p_cm3: np.ndarray,
    ni_cm3: np.ndarray | float,
    B_rad_cm3_s: float = 2.0e-11,
    Cn_auger_cm6_s: float = 1.5e-30,
    Cp_auger_cm6_s: float = 1.5e-30,
    tau_n_s: float = 1.0e-7,
    tau_p_s: float = 1.0e-7,
) -> dict[str, np.ndarray]:
    """
    Computes spatially resolved recombination rates across the semiconductor profile.
    
    Parameters
    ----------
    carrier_n_cm3 : np.ndarray
        Free electron concentration (cm^-3).
    carrier_p_cm3 : np.ndarray
        Free hole concentration (cm^-3).
    ni_cm3 : np.ndarray or float
        Intrinsic carrier concentration (cm^-3).
    B_rad_cm3_s : float
        Radiative recombination coefficient (cm^3/s). Default: 2.0e-11 cm^3/s.
    Cn_auger_cm6_s, Cp_auger_cm6_s : float
        Auger recombination coefficients (cm^6/s). Default: 1.5e-30 cm^6/s.
    tau_n_s, tau_p_s : float
        SRH carrier lifetimes (s). Default: 1.0e-7 s.
        
    Returns
    -------
    dict with:
        'R_rad': Radiative recombination rate profile (cm^-3 s^-1)
        'R_srh': SRH non-radiative recombination rate profile (cm^-3 s^-1)
        'R_auger': Auger non-radiative recombination rate profile (cm^-3 s^-1)
        'R_total': Total recombination rate profile (cm^-3 s^-1)
    """
    n = np.clip(carrier_n_cm3, 1.0e-10, None)
    p = np.clip(carrier_p_cm3, 1.0e-10, None)
    ni2 = (np.asarray(ni_cm3))**2
    
    excess_np = np.maximum(n * p - ni2, 0.0)
    
    # 1. Bimolecular Radiative Recombination
    R_rad = B_rad_cm3_s * excess_np
    
    # 2. Shockley-Read-Hall (SRH) Recombination
    denom_srh = tau_p_s * (n + np.sqrt(ni2)) + tau_n_s * (p + np.sqrt(ni2))
    denom_srh = np.maximum(denom_srh, 1.0e-25)
    R_srh = excess_np / denom_srh
    
    # 3. Auger Recombination (eeh and hhe processes)
    R_auger = (Cn_auger_cm6_s * n + Cp_auger_cm6_s * p) * excess_np
    
    R_total = R_rad + R_srh + R_auger
    
    return {
        'R_rad': R_rad,
        'R_srh': R_srh,
        'R_auger': R_auger,
        'R_total': R_total,
    }


def compute_integrated_led_efficiencies(
    z_nm: np.ndarray,
    R_profiles: dict[str, np.ndarray],
    active_mask: np.ndarray | None = None,
    device_area_cm2: float = 1.0e-3,
    eta_extraction: float = 0.20,
    peak_energy_ev: float = 2.75,
) -> dict[str, float]:
    """
    Integrates spatial recombination rates across the active region to compute
    Internal Quantum Efficiency (IQE), External Quantum Efficiency (EQE),
    and emitted optical power.
    """
    dz_cm = np.gradient(z_nm * 1.0e-7)  # nm to cm
    
    if active_mask is None:
        active_mask = np.ones_like(z_nm, dtype=bool)
        
    mask = active_mask & (dz_cm > 0)
    
    # Spatial integrals over active region (cm^-2 s^-1)
    int_R_rad = float(np.sum(R_profiles['R_rad'][mask] * dz_cm[mask]))
    int_R_srh = float(np.sum(R_profiles['R_srh'][mask] * dz_cm[mask]))
    int_R_auger = float(np.sum(R_profiles['R_auger'][mask] * dz_cm[mask]))
    int_R_tot = int_R_rad + int_R_srh + int_R_auger
    
    # Internal Quantum Efficiency (IQE)
    iqe = (int_R_rad / int_R_tot) if int_R_tot > 0 else 0.0
    
    # External Quantum Efficiency (EQE)
    eqe = eta_extraction * iqe
    
    # Total photon emission rate (photons/s)
    phi_photons = device_area_cm2 * int_R_rad
    
    # Radiative output optical power (mW)
    photon_energy_j = peak_energy_ev * Q_ELEC
    p_opt_mw = (eta_extraction * phi_photons * photon_energy_j) * 1.0e3
    
    return {
        'int_R_rad_cm2_s': int_R_rad,
        'int_R_srh_cm2_s': int_R_srh,
        'int_R_auger_cm2_s': int_R_auger,
        'int_R_total_cm2_s': int_R_tot,
        'iqe': float(iqe),
        'eqe': float(eqe),
        'phi_photons_per_s': float(phi_photons),
        'p_opt_mw': float(p_opt_mw),
    }


def generate_electroluminescence_spectrum(
    peak_wavelength_nm: float = 450.0,
    fwhm_nm: float = 25.0,
    temperature_k: float = 300.0,
    wavelength_range_nm: tuple[float, float] = (380.0, 520.0),
    num_points: int = 200,
) -> dict[str, np.ndarray]:
    """
    Computes physical spontaneous emission spectrum I_EL(lambda) accounting for
    asymmetric thermal tail ~ exp(-E/kT) and inhomogeneous Gaussian alloy broadening.
    
    Parameters
    ----------
    peak_wavelength_nm : float
        Target peak emission wavelength (nm).
    fwhm_nm : float
        Spectral full-width at half-maximum (nm).
    temperature_k : float
        Operating junction temperature (K).
    wavelength_range_nm : tuple
        (lambda_min, lambda_max) in nm.
    num_points : int
        Grid resolution points.
        
    Returns
    -------
    dict with:
        'wavelength_nm': np.ndarray
        'energy_ev': np.ndarray
        'intensity_norm': np.ndarray (normalized to peak = 1.0)
        'peak_wavelength_nm': float
        'peak_energy_ev': float
        'fwhm_nm': float
        'fwhm_ev': float
    """
    wl = np.linspace(wavelength_range_nm[0], wavelength_range_nm[1], num_points)
    energy_ev = HC_EV_NM / wl
    
    e_peak_ev = HC_EV_NM / peak_wavelength_nm
    kt_ev = (K_BOLTZ * temperature_k) / Q_ELEC
    
    # Gaussian alloy and thermal broadening in wavelength space
    sigma_wl = fwhm_nm / 2.35482
    delta_wl = wl - peak_wavelength_nm
    
    # Asymmetric thermal tail on the high-energy (shorter wavelength) side
    kt_nm = (peak_wavelength_nm**2 * kt_ev) / HC_EV_NM
    thermal_tail = np.where(
        wl < peak_wavelength_nm,
        np.exp((wl - peak_wavelength_nm) / max(kt_nm * 2.5, sigma_wl)),
        np.exp(-0.5 * (delta_wl / sigma_wl)**2)
    )
    gaussian_base = np.exp(-0.5 * (delta_wl / sigma_wl)**2)
    shape = 0.75 * gaussian_base + 0.25 * thermal_tail
    
    # Normalize to peak intensity = 1.0
    intensity_norm = shape / np.max(shape)
    
    # Calculate exact realized peak and continuous interpolated FWHM
    idx_max = int(np.argmax(intensity_norm))
    actual_peak_wl = float(wl[idx_max])
    actual_peak_ev = float(energy_ev[idx_max])
    
    # Continuous linear interpolation for half-maximum crossings
    try:
        left_wl = float(np.interp(0.5, intensity_norm[:idx_max], wl[:idx_max])) if idx_max > 0 else float(wl[0])
        right_wl = float(np.interp(0.5, intensity_norm[idx_max:][::-1], wl[idx_max:][::-1])) if idx_max < len(wl)-1 else float(wl[-1])
        actual_fwhm_nm = abs(right_wl - left_wl)
        actual_fwhm_ev = abs(HC_EV_NM / left_wl - HC_EV_NM / right_wl)
    except Exception:
        actual_fwhm_nm = float(fwhm_nm)
        actual_fwhm_ev = float(fwhm_nm * e_peak_ev / peak_wavelength_nm)
        
    return {
        'wavelength_nm': wl,
        'energy_ev': energy_ev,
        'intensity_norm': intensity_norm,
        'peak_wavelength_nm': actual_peak_wl,
        'peak_energy_ev': actual_peak_ev,
        'fwhm_nm': actual_fwhm_nm,
        'fwhm_ev': actual_fwhm_ev,
    }


def analyze_led_iv_curve(
    voltage_v: np.ndarray,
    current_ma: np.ndarray,
    device_area_cm2: float = 1.0e-3,
    threshold_current_ma: float = 1.0,
    nominal_current_ma: float = 20.0,
    temperature_k: float = 300.0,
) -> dict[str, float]:
    """
    Extracts key LED electrical figures of merit from an I-V curve:
    - Turn-on voltage (V_on)
    - Forward operating voltage at nominal current (V_f at 20 mA)
    - Series resistance (R_s) in forward bias
    - Shunt resistance (R_sh) near zero bias
    - Average ideality factor (n)
    """
    sort_idx = np.argsort(voltage_v)
    v = voltage_v[sort_idx]
    i = current_ma[sort_idx]
    j = i / device_area_cm2  # mA/cm^2
    
    # 1. Turn-on Voltage (V_on at threshold current)
    v_on = None
    if np.any(i >= threshold_current_ma):
        idx_on = int(np.argmax(i >= threshold_current_ma))
        if idx_on > 0:
            v0, v1 = v[idx_on - 1], v[idx_on]
            i0, i1 = i[idx_on - 1], i[idx_on]
            v_on = float(v0 + (threshold_current_ma - i0) * (v1 - v0) / (i1 - i0))
        else:
            v_on = float(v[0])
    else:
        v_on = float(v[-1])
        
    # 2. Operating Voltage at Nominal Current (e.g. 20 mA)
    v_nominal = None
    if np.any(i >= nominal_current_ma):
        v_nominal = float(np.interp(nominal_current_ma, i, v))
    else:
        v_nominal = float(v[-1])
        
    # 3. Dynamic Resistance Rs at nominal forward bias
    if len(v) >= 4 and np.any(i > 5.0):
        fwd_mask = i > 5.0
        dv = np.gradient(v[fwd_mask])
        di = np.gradient(i[fwd_mask] * 1.0e-3)  # mA to A
        rs = float(np.median(dv / np.maximum(di, 1.0e-12)))
    else:
        rs = 0.0
        
    # 4. Sub-threshold Shunt Resistance Rsh near 0V
    zero_mask = np.abs(v) <= 0.5
    if np.sum(zero_mask) >= 2:
        dv_z = np.gradient(v[zero_mask])
        di_z = np.gradient(i[zero_mask] * 1.0e-3)
        rsh = float(np.mean(dv_z / np.maximum(np.abs(di_z), 1.0e-15)))
    else:
        rsh = 1.0e8
        
    # 5. Ideality factor n around turn-on knee
    vt = (K_BOLTZ * temperature_k) / Q_ELEC
    knee_mask = (v >= (v_on - 0.5)) & (v <= (v_on + 0.2)) & (i > 1.0e-4)
    if np.sum(knee_mask) >= 3:
        ln_i = np.log(i[knee_mask])
        v_knee = v[knee_mask]
        d_ln_i_dv = np.gradient(ln_i, v_knee)
        n_vals = 1.0 / np.maximum(d_ln_i_dv * vt, 1.0e-6)
        n_avg = float(np.median(n_vals))
    else:
        n_avg = 1.5
        
    return {
        'turn_on_voltage_v': float(v_on),
        'forward_voltage_at_nominal_v': float(v_nominal),
        'series_resistance_ohm': float(rs),
        'shunt_resistance_ohm': float(rsh),
        'ideality_factor_n': float(n_avg),
    }


def compute_led_efficiency_droop_curve(
    current_density_a_cm2: np.ndarray,
    A_srh_s: float = 1.0e7,
    B_rad_cm3_s: float = 2.0e-11,
    C_auger_cm6_s: float = 1.5e-30,
    d_active_cm: float = 3.0e-7,
) -> dict[str, np.ndarray]:
    """
    Analytic and numerical ABC-model representation of LED efficiency droop.
    Relates injected current density J to carrier density n via:
        J / (q * d) = A*n + B*n^2 + C*n^3
    and computes IQE(J) = B*n^2 / (A*n + B*n^2 + C*n^3).
    """
    j_vals = np.asarray(current_density_a_cm2)
    
    # Solve cubic equation for n at each J: C*n^3 + B*n^2 + A*n - J/(q*d) = 0
    n_carriers = []
    iqe_vals = []
    
    for j in j_vals:
        j_a_cm2 = max(float(j), 1.0e-6)
        r_gen = j_a_cm2 / (Q_ELEC * d_active_cm)  # cm^-3 s^-1
        
        # Polynomial coefficients: C*n^3 + B*n^2 + A*n - r_gen = 0
        roots = np.roots([C_auger_cm6_s, B_rad_cm3_s, A_srh_s, -r_gen])
        real_roots = roots[np.isreal(roots) & (roots > 0)].real
        if len(real_roots) > 0:
            n = real_roots[0]
        else:
            n = 1.0e17
            
        r_srh = A_srh_s * n
        r_rad = B_rad_cm3_s * n**2
        r_aug = C_auger_cm6_s * n**3
        r_tot = r_srh + r_rad + r_aug
        
        iqe = r_rad / r_tot if r_tot > 0 else 0.0
        n_carriers.append(n)
        iqe_vals.append(iqe)
        
    iqe_arr = np.array(iqe_vals)
    max_iqe = np.max(iqe_arr) if len(iqe_arr) > 0 else 1.0
    normalized_iqe = iqe_arr / max_iqe if max_iqe > 0 else iqe_arr
    
    return {
        'current_density_a_cm2': j_vals,
        'carrier_density_cm3': np.array(n_carriers),
        'iqe': iqe_arr,
        'normalized_iqe': normalized_iqe,
        'peak_iqe': float(max_iqe),
        'j_peak_a_cm2': float(j_vals[int(np.argmax(iqe_arr))]),
    }
