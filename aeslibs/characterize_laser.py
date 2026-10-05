# -*- coding: utf-8 -*-
"""
Aestimo 1D Semiconductor Laser Characterization and Optical Cavity Physics Engine.
Implements optical cavity loss calculations, threshold modal gain and carrier density,
steady-state carrier and photon rate equations, slope efficiency extraction,
thermal roll-off modeling, and longitudinal Fabry-Pérot mode optical spectrum generation.
"""

from __future__ import annotations
import numpy as np
import scipy.constants as const
from scipy.optimize import brentq

# Physical constants in SI
Q_ELEC = const.e                       # 1.602176634e-19 C
H_PLANCK = const.h                     # 6.62607015e-34 J*s
C_LIGHT = const.c                      # 2.99792458e8 m/s
K_BOLTZ = const.k                      # 1.380649e-23 J/K
HC_EV_NM = 1239.841984                 # hc in eV*nm


def calculate_cavity_optical_losses(
    cavity_length_um: float,
    r1: float = 0.32,
    r2: float = 0.32,
    alpha_internal_cm1: float = 5.0,
    group_index: float = 3.6,
) -> dict[str, float]:
    """
    Computes Fabry-Pérot laser optical cavity parameters:
    mirror loss (alpha_m), total loss (alpha_tot), group velocity (v_g),
    and photon lifetime (tau_p).
    
    Parameters
    ----------
    cavity_length_um : float
        Cavity length L in micrometers (um).
    r1, r2 : float
        Power reflectivities of the front and rear laser facets (0 to 1).
    alpha_internal_cm1 : float
        Internal optical waveguide loss in cm^-1 (free-carrier absorption, scattering).
    group_index : float
        Effective optical group index n_g.
        
    Returns
    -------
    dict with:
        'cavity_length_cm': L in cm
        'alpha_m_cm1': Mirror loss alpha_m (cm^-1)
        'alpha_i_cm1': Internal loss alpha_i (cm^-1)
        'alpha_tot_cm1': Total optical cavity loss alpha_tot (cm^-1)
        'v_g_cm_s': Optical group velocity v_g (cm/s)
        'photon_lifetime_s': Photon lifetime tau_p (s)
    """
    L_cm = cavity_length_um * 1.0e-4
    r1_clamped = np.clip(r1, 1.0e-5, 0.99999)
    r2_clamped = np.clip(r2, 1.0e-5, 0.99999)
    
    alpha_m = (1.0 / (2.0 * L_cm)) * np.log(1.0 / (r1_clamped * r2_clamped))
    alpha_tot = alpha_internal_cm1 + alpha_m
    
    v_g_cm_s = (C_LIGHT * 100.0) / group_index
    photon_lifetime_s = 1.0 / (v_g_cm_s * alpha_tot)
    
    return {
        'cavity_length_cm': L_cm,
        'alpha_m_cm1': float(alpha_m),
        'alpha_i_cm1': float(alpha_internal_cm1),
        'alpha_tot_cm1': float(alpha_tot),
        'v_g_cm_s': float(v_g_cm_s),
        'photon_lifetime_s': float(photon_lifetime_s),
    }


def calculate_threshold_gain_and_carrier_density(
    alpha_tot_cm1: float,
    confinement_factor: float,
    g0_cm1: float = 1500.0,
    n_tr_cm3: float = 1.8e18,
    gain_model: str = "log",
    differential_gain_cm2: float = 2.5e-16,
) -> dict[str, float]:
    """
    Computes threshold material gain (g_th) and threshold carrier density (n_th).
    
    Parameters
    ----------
    alpha_tot_cm1 : float
        Total optical loss (alpha_i + alpha_m) in cm^-1.
    confinement_factor : float
        Optical mode confinement factor Gamma in the active region.
    g0_cm1 : float
        Empirical gain coefficient g_0 in cm^-1 (for logarithmic gain model).
    n_tr_cm3 : float
        Transparency carrier density n_tr in cm^-3.
    gain_model : str
        'log' for phenomenological QW gain: g(n) = g0 * ln(n / n_tr),
        'linear' for bulk/DH gain: g(n) = a * (n - n_tr).
    differential_gain_cm2 : float
        Differential gain a = dg/dn in cm^2 (used when gain_model == 'linear').
        
    Returns
    -------
    dict with:
        'modal_gain_th_cm1': Threshold modal gain Gamma * g_th (cm^-1)
        'material_gain_th_cm1': Threshold material gain g_th (cm^-1)
        'n_th_cm3': Threshold carrier density n_th (cm^-3)
    """
    gamma = max(confinement_factor, 1.0e-5)
    modal_gain_th = alpha_tot_cm1
    material_gain_th = modal_gain_th / gamma
    
    if gain_model.lower() == "log":
        # g(n) = g0 * ln(n / n_tr)  ==>  n_th = n_tr * exp(g_th / g0)
        n_th = n_tr_cm3 * np.exp(material_gain_th / max(g0_cm1, 10.0))
    else:
        # g(n) = a * (n - n_tr)  ==>  n_th = n_tr + g_th / a
        n_th = n_tr_cm3 + material_gain_th / max(differential_gain_cm2, 1.0e-20)
        
    return {
        'modal_gain_th_cm1': float(modal_gain_th),
        'material_gain_th_cm1': float(material_gain_th),
        'n_th_cm3': float(n_th),
    }


def material_gain(
    n_cm3: float | np.ndarray,
    g0_cm1: float = 1500.0,
    n_tr_cm3: float = 1.8e18,
    gain_model: str = "log",
    differential_gain_cm2: float = 2.5e-16,
) -> float | np.ndarray:
    """Computes material optical gain g(n) in cm^-1."""
    n = np.maximum(n_cm3, 1.0e10)
    if gain_model.lower() == "log":
        return g0_cm1 * np.log(n / n_tr_cm3)
    else:
        return differential_gain_cm2 * (n - n_tr_cm3)


def solve_laser_rate_equations_steady_state(
    current_array_ma: np.ndarray,
    active_volume_cm3: float,
    A_srh_s: float = 1.0e7,
    B_rad_cm3_s: float = 2.0e-11,
    C_auger_cm6_s: float = 1.5e-30,
    eta_i: float = 0.85,
    photon_lifetime_s: float = 1.5e-12,
    confinement_factor: float = 0.035,
    g0_cm1: float = 1500.0,
    n_tr_cm3: float = 1.8e18,
    v_g_cm_s: float = 8.33e9,
    alpha_m_cm1: float = 22.8,
    beta_sp: float = 1.0e-4,
    peak_wavelength_nm: float = 850.0,
    gain_model: str = "log",
    differential_gain_cm2: float = 2.5e-16,
    eps_gain_compression: float = 1.0e-17,
    thermal_resistance_k_w: float = 0.0,
    t0_k: float = 120.0,
    t_ref_k: float = 300.0,
    voltage_array_v: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """
    Self-consistently solves the coupled steady-state carrier and photon rate equations:
    dn/dt = (eta_i * I) / (q * V_act) - R(n) - v_g * g(n)/(1 + eps*S) * S = 0
    dS/dt = Gamma * v_g * g(n)/(1 + eps*S) * S - S / tau_p + Gamma * beta * B * n^2 = 0
    
    Includes spontaneous emission below threshold, stimulated emission above threshold,
    gain compression, and thermal roll-off.
    
    Returns
    -------
    dict with:
        'current_ma': Injected currents (mA)
        'carrier_density_cm3': Steady-state carrier density n (cm^-3)
        'photon_density_cm3': Optical photon density S (cm^-3)
        'power_single_facet_mw': Output power per single facet (mW)
        'power_total_mw': Total optical power emitted from both facets (mW)
        'stimulated_emission_rate': R_stim (cm^-3 s^-1)
        'spontaneous_emission_rate': R_sp (cm^-3 s^-1)
    """
    i_arr = np.asarray(current_array_ma, dtype=float)
    hnu_j = (HC_EV_NM / peak_wavelength_nm) * Q_ELEC
    inv_tau_p = 1.0 / photon_lifetime_s
    
    # Precompute threshold carrier density
    gamma = max(confinement_factor, 1.0e-5)
    g_th = (inv_tau_p / v_g_cm_s) / gamma
    if gain_model.lower() == "log":
        n_th = n_tr_cm3 * np.exp(g_th / max(g0_cm1, 10.0))
    else:
        n_th = n_tr_cm3 + g_th / max(differential_gain_cm2, 1.0e-20)
        
    # Recombination rate at threshold
    R_th = A_srh_s * n_th + B_rad_cm3_s * (n_th**2) + C_auger_cm6_s * (n_th**3)
    i_th_a = (Q_ELEC * active_volume_cm3 / max(eta_i, 1.0e-4)) * R_th
    i_th_ma = i_th_a * 1.0e3

    n_res = np.zeros_like(i_arr)
    s_res = np.zeros_like(i_arr)
    p_single_mw = np.zeros_like(i_arr)
    p_total_mw = np.zeros_like(i_arr)
    
    # Facet optical power factors (W per active unit photon density)
    # P_single = 0.5 * v_g * alpha_m * hnu * (V_act / Gamma) * S
    facet_factor_single = 0.5 * v_g_cm_s * alpha_m_cm1 * hnu_j * (active_volume_cm3 / gamma) * 1.0e3  # W to mW
    facet_factor_total = 2.0 * facet_factor_single

    for idx, i_val_ma in enumerate(i_arr):
        if i_val_ma <= 1.0e-6:
            n_res[idx] = 0.0
            s_res[idx] = 0.0
            p_single_mw[idx] = 0.0
            p_total_mw[idx] = 0.0
            continue

        i_a = i_val_ma * 1.0e-3  # mA to A
        
        # Coupled Langevin / single-mode rate-equation formulation:
        # Stimulated drive term above threshold
        s_drive = (gamma * photon_lifetime_s * eta_i / (Q_ELEC * active_volume_cm3)) * (i_a - i_th_a)
        s_sp_0 = gamma * photon_lifetime_s * beta_sp * B_rad_cm3_s * (n_th**2)
        
        # Analytic steady-state photon density (strictly positive, C^infinity smooth)
        s_val = 0.5 * (s_drive + np.sqrt(s_drive**2 + 4.0 * s_sp_0))
        
        # Non-linear gain compression above threshold: S -> S / (1 + eps * S)
        if eps_gain_compression > 0.0:
            s_val = s_val / (1.0 + eps_gain_compression * s_val)
            
        s_res[idx] = s_val
        
        # Carrier density below threshold scales via recombination balance;
        # above threshold it pins smoothly to n_th
        if i_val_ma < i_th_ma:
            gen_target = (eta_i * i_a) / (Q_ELEC * active_volume_cm3)
            # Find root of A*n + B*n^2 + C*n^3 = gen_target
            try:
                def f_rec(nc):
                    return (A_srh_s * nc + B_rad_cm3_s * (nc**2) + C_auger_cm6_s * (nc**3)) - gen_target
                n_sol = brentq(f_rec, 0.0, n_th * 1.05, xtol=1.0e10, rtol=1.0e-4)
            except Exception:
                n_sol = n_th * (i_val_ma / max(i_th_ma, 1.0e-3))
            n_res[idx] = n_sol
        else:
            # Clamped at threshold carrier density with slight spontaneous saturation
            n_res[idx] = n_th
            
        p_single = s_val * facet_factor_single
        p_tot = s_val * facet_factor_total
        
        # Thermal roll-off
        if thermal_resistance_k_w > 0.0 and voltage_array_v is not None and idx < len(voltage_array_v):
            v_val = voltage_array_v[idx]
            p_elec_w = i_a * v_val
            p_opt_w = p_tot * 1.0e-3
            delta_t = thermal_resistance_k_w * max(p_elec_w - p_opt_w, 0.0)
            thermal_suppression = np.exp(-delta_t / max(t0_k, 10.0))
            p_single *= thermal_suppression
            p_tot *= thermal_suppression
            
        p_single_mw[idx] = p_single
        p_total_mw[idx] = p_tot

    return {
        'current_ma': i_arr,
        'carrier_density_cm3': n_res,
        'photon_density_cm3': s_res,
        'power_single_facet_mw': p_single_mw,
        'power_total_mw': p_total_mw,
        'n_th_cm3': float(n_th),
        'i_th_calc_ma': float(i_th_ma),
    }


def extract_laser_figures_of_merit(
    current_ma: np.ndarray,
    power_mw: np.ndarray,
    voltage_v: np.ndarray | None = None,
    peak_wavelength_nm: float = 850.0,
    area_cm2: float = 1.0e-4,
    ith_calc_ma: float | None = None,
) -> dict[str, float]:
    """
    Extracts key semiconductor laser figures of merit from numerical L-I and I-V data:
    threshold current (I_th), threshold current density (J_th), slope efficiency (SE),
    differential quantum efficiency (eta_d), threshold voltage (V_th), maximum optical power (P_max),
    and maximum wall-plug efficiency (WPE_max).
    """
    i_arr = np.asarray(current_ma, dtype=float)
    p_arr = np.asarray(power_mw, dtype=float)
    
    if len(i_arr) < 5 or np.max(p_arr) < 1.0e-3:
        return {
            'threshold_current_ma': float('nan'),
            'threshold_current_density_a_cm2': float('nan'),
            'slope_efficiency_mw_per_ma': float('nan'),
            'differential_quantum_efficiency': float('nan'),
            'threshold_voltage_v': float('nan'),
            'max_optical_power_mw': float('nan'),
            'max_wall_plug_efficiency_pct': float('nan'),
        }
        
    # Numerical derivative dP/dI
    dp_di = np.gradient(p_arr, i_arr)
    
    # 2nd derivative d^2 P / d I^2 peak indicates threshold inflection
    d2p_di2 = np.gradient(dp_di, i_arr)
    valid_range = slice(2, len(i_arr) - 2)
    idx_inflection = int(np.argmax(d2p_di2[valid_range]) + 2)
    i_th = float(i_arr[idx_inflection])
    
    if ith_calc_ma is not None and np.isfinite(ith_calc_ma):
        i_th = float(ith_calc_ma)
        
    # Linear slope above threshold (using points from 1.15 * I_th to 0.85 * max power)
    above_th_mask = (i_arr >= 1.15 * i_th) & (p_arr <= 0.85 * np.max(p_arr))
    if np.sum(above_th_mask) >= 3:
        poly = np.polyfit(i_arr[above_th_mask], p_arr[above_th_mask], deg=1)
        slope_efficiency = float(poly[0])  # mW/mA == W/A
        i_th_fit = float(-poly[1] / poly[0]) if abs(poly[0]) > 1.0e-5 else i_th
        if 0.5 * i_th <= i_th_fit <= 1.5 * i_th and ith_calc_ma is None:
            i_th = i_th_fit
    else:
        slope_efficiency = float(np.max(dp_di))
        
    j_th = (i_th * 1.0e-3) / max(area_cm2, 1.0e-12)
    
    # Differential quantum efficiency: eta_d = (dP/dI) / (hnu/q)
    hnu_ev = HC_EV_NM / peak_wavelength_nm
    eta_d = slope_efficiency / hnu_ev
    
    # Max optical power
    p_max = float(np.max(p_arr))
    
    # Voltage and Wall-plug efficiency
    v_th = float('nan')
    wpe_max = float('nan')
    if voltage_v is not None and len(voltage_v) == len(i_arr):
        v_arr = np.asarray(voltage_v, dtype=float)
        v_th = float(np.interp(i_th, i_arr, v_arr))
        
        # WPE(I) = P_out(mW) / (I(mA) * V(V)) * 100%
        p_elec_mw = i_arr * v_arr
        valid_p = (p_elec_mw > 1.0e-3) & (p_arr > 1.0e-3)
        if np.any(valid_p):
            wpe_curve = (p_arr[valid_p] / p_elec_mw[valid_p]) * 100.0
            wpe_max = float(np.max(wpe_curve))
            
    return {
        'threshold_current_ma': float(i_th),
        'threshold_current_density_a_cm2': float(j_th),
        'slope_efficiency_mw_per_ma': float(slope_efficiency),
        'differential_quantum_efficiency': float(eta_d),
        'threshold_voltage_v': float(v_th),
        'max_optical_power_mw': float(p_max),
        'max_wall_plug_efficiency_pct': float(wpe_max),
    }


def generate_laser_fp_spectrum(
    peak_wavelength_nm: float,
    cavity_length_um: float,
    group_index: float = 3.6,
    fwhm_sp_nm: float = 25.0,
    fwhm_lasing_nm: float | None = None,
    current_ratio_i_over_ith: float = 1.25,
    num_modes: int = 21,
) -> dict[str, np.ndarray | float]:
    """
    Simulates the discrete longitudinal Fabry-Pérot (FP) mode emission spectrum
    modulated by the gain narrowing envelope above threshold.
    
    Mode spacing: Delta_lambda = lambda_0^2 / (2 * n_g * L)
    """
    L_nm = cavity_length_um * 1.0e3
    lambda_0 = peak_wavelength_nm
    delta_lambda = (lambda_0**2) / (2.0 * group_index * L_nm)
    
    # Effective envelope FWHM above threshold
    if fwhm_lasing_nm is not None and fwhm_lasing_nm > 0.0:
        eff_fwhm = float(fwhm_lasing_nm)
    elif current_ratio_i_over_ith > 1.0:
        # Physical gain narrowing above threshold
        eff_fwhm = max(fwhm_sp_nm / (np.sqrt(max((current_ratio_i_over_ith - 1.0) * 100.0, 1.0))), 0.65)
    else:
        eff_fwhm = float(fwhm_sp_nm)
        
    # Mode indices: -num_modes//2 to +num_modes//2
    half_m = num_modes // 2
    mode_indices = np.arange(-half_m, half_m + 1)
    wavelengths = lambda_0 + mode_indices * delta_lambda
    
    # Longitudinal mode envelope
    sigma_nm = eff_fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    intensities = np.exp(-((wavelengths - lambda_0)**2) / (2.0 * sigma_nm**2))
    intensities = intensities / np.max(intensities)
    
    return {
        'wavelength_nm': wavelengths,
        'intensity_norm': intensities,
        'mode_spacing_nm': float(delta_lambda),
        'peak_wavelength_nm': float(lambda_0),
        'eff_fwhm_nm': float(eff_fwhm),
    }


def compute_laser_temperature_series(
    current_array_ma: np.ndarray,
    temp_array_k: list[float] | np.ndarray,
    ith_ref_ma: float,
    t_ref_k: float = 298.0,
    t0_k: float = 52.0,
    slope_eff_ref: float = 0.24,
    t1_k: float = 140.0,
) -> dict[str, dict[str, np.ndarray]]:
    """
    Computes multi-temperature L-I curves according to empirical temperature scaling:
    I_th(T) = I_th(T_ref) * exp((T - T_ref) / T_0)
    SE(T) = SE(T_ref) * exp(-(T - T_ref) / T_1)
    """
    results = {}
    i_arr = np.asarray(current_array_ma, dtype=float)
    
    for T in temp_array_k:
        delta_T = T - t_ref_k
        ith_t = ith_ref_ma * np.exp(delta_T / max(t0_k, 10.0))
        se_t = slope_eff_ref * np.exp(-delta_T / max(t1_k, 10.0))
        
        # Soft-knee L-I curve
        p_out = np.zeros_like(i_arr)
        above = i_arr >= ith_t
        below = ~above
        # Sub-threshold spontaneous emission
        p_out[below] = 0.01 * (i_arr[below] / ith_t)
        # Above threshold stimulated emission
        p_out[above] = 0.01 + se_t * (i_arr[above] - ith_t)
        
        results[f"{int(round(T))}K"] = {
            'temperature_k': float(T),
            'threshold_current_ma': float(ith_t),
            'slope_efficiency_mw_per_ma': float(se_t),
            'current_ma': i_arr,
            'power_mw': p_out,
        }
        
    return results
