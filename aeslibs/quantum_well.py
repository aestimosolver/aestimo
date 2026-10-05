#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Aestimo 1D - Quantum-Well Confined-State Solver
aeslibs/quantum_well.py

Provides physically rigorous calculation of quantized confined energy levels,
envelope wavefunctions, optical transitions (overlaps, transition energies, wavelengths),
Quantum-Confined Stark Effect (QCSE), 2D/3D carrier statistics, and optional
self-consistent Schrödinger-Poisson coupling in semiconductor heterostructures.

Equations:
  BenDaniel-Duke variable-mass 1D Schrödinger equation:
    - (ħ² / 2) d/dz [ (1 / m*(z)) dψ/dz ] + V(z) ψ(z) = E ψ(z)
  
  Boundary Conditions:
    ψ(z⁻) = ψ(z⁺)
    (1/m*(z⁻)) dψ/dz|⁻ = (1/m*(z⁺)) dψ/dz|⁺
"""

from __future__ import annotations
import math
import numpy as np
import scipy.linalg as la
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Any

import database

# =============================================================================
# PHYSICAL CONSTANTS (CODATA 2018 / SI & Solid State Units)
# =============================================================================
H_BAR = 1.054571817e-34       # J s
Q_ELEM = 1.602176634e-19      # C (J/eV)
M_ELEM = 9.1093837015e-31     # kg
K_BOLTZ = 1.380649e-23        # J/K
C_LIGHT = 2.99792458e8        # m/s
HC_EV_NM = 1239.841984        # hc in eV nm
EPSILON_0 = 8.8541878128e-12  # F/m

# Kinetic factor in eV * nm^2: (ħ² / (2 * m0 * q))
# ħ² / (2 * m0) = (1.054571817e-34)² / (2 * 9.1093837015e-31) = 6.104263e-39 J m²
# in eV: 6.104263e-39 / 1.602176634e-19 = 3.80998212e-20 eV m² = 0.0380998212 eV nm²
HBAR2_OVER_2M0_EV_NM2 = 0.0380998212


# =============================================================================
# STRUCTURED DATA OBJECT FOR QW RESULTS
# =============================================================================
@dataclass
class QuantumWellResult:
    """
    Structured outcome of a Quantum-Well confined-state calculation.
    """
    # Spatial mesh and potentials
    z: np.ndarray                            # Spatial coordinate (nm)
    ec_profile: np.ndarray                   # Local conduction band profile (eV)
    ev_profile: np.ndarray                   # Local valence band profile (eV)
    
    # Electron confined states
    electron_energies: np.ndarray            # E_e,n (eV relative to band reference)
    electron_wavefunctions: np.ndarray       # psi_e,n(z) (nm^-1/2)
    electron_probability: np.ndarray         # |psi_e,n(z)|^2 (nm^-1)
    
    # Hole confined states (Heavy Hole & Light Hole)
    hole_energies: np.ndarray                # E_h,j (eV relative to band reference)
    hole_wavefunctions: np.ndarray           # psi_h,j(z) (nm^-1/2)
    hole_probability: np.ndarray             # |psi_h,j(z)|^2 (nm^-1)
    hole_types: List[str]                    # ['hh1', 'lh1', 'hh2', ...]
    
    # Optical transitions
    transition_energies: np.ndarray          # E_trans[i, j] = E_e[i] - E_h[j] (eV)
    transition_wavelengths: np.ndarray       # lambda[i, j] = hc / E_trans (nm)
    overlap_integrals: np.ndarray            # Gamma[i, j] = |<psi_e,i | psi_h,j>|^2
    dominant_transitions: List[Dict[str, Any]] # List of key optical transitions
    
    # Carrier populations (Fermi-Dirac 2D sheet & 3D volumetric densities)
    sheet_densities_e: np.ndarray            # n_2D,n (cm^-2)
    sheet_densities_h: np.ndarray            # p_2D,j (cm^-2)
    quantum_n: np.ndarray                    # n_qw(z) (cm^-3)
    quantum_p: np.ndarray                    # p_qw(z) (cm^-3)
    
    # Diagnostics & numerical metadata
    normalization_errors: Dict[str, float]   # Residual |1 - int |psi|^2 dz|
    qw_regions: List[Dict[str, Any]]         # Identified well boundaries
    convergence: Dict[str, Any]              # Iteration info, residuals, status
    temperature_k: float = 300.0
    field_v_cm: float = 0.0                  # Average electric field across well
    active_scale_mode: str = "AUTO"


# =============================================================================
# CORE BEN-DANIEL-DUKE EIGENVALUE SOLVER
# =============================================================================
def solve_effective_mass_schrodinger_1d(
    z_nm: np.ndarray,
    potential_ev: np.ndarray,
    m_eff_rel: np.ndarray,
    num_states: int = 3,
    particle_type: str = "electron",
    well_mask: Optional[np.ndarray] = None,
    barrier_ref_ev: Optional[float] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, float]]:
    """
    Solves the 1D variable effective-mass Schrödinger equation using the exact
    BenDaniel-Duke finite-difference discretization and tridiagonal eigensolver.
    
    Parameters:
      z_nm: 1D array of spatial positions in nanometers (must be uniformly spaced).
      potential_ev: 1D array of potential energy in eV.
      m_eff_rel: 1D array of effective mass relative to m0 (m*(z) / m0).
      num_states: Number of bound states to solve for.
      particle_type: 'electron' or 'hole' (inverts potential for holes).
      well_mask: Boolean mask indicating the well interior.
      barrier_ref_ev: Reference barrier height for bound-state filtering.
      
    Returns:
      energies: 1D array of eigenenergies in eV.
      wavefunctions: 2D array of normalized envelope wavefunctions psi(z) in nm^-1/2.
      probabilities: 2D array of probability densities |psi(z)|^2 in nm^-1.
      norm_errors: Dictionary mapping state index to normalization error.
    """
    N = len(z_nm)
    if N < 5:
        raise ValueError("Grid too small for Schrödinger solution (N < 5).")
    
    dz = z_nm[1] - z_nm[0]
    if dz <= 0:
        raise ValueError(f"Invalid non-positive grid step dz={dz} nm.")
    
    # Clip effective masses to physically realistic bounds
    m_eff = np.clip(np.asarray(m_eff_rel, dtype=float), 0.01, 10.0)
    
    # Setup effective potential
    # For electrons: V(z) = Ec(z)
    # For holes: V_h(z) = -Ev(z) (so maximum in Ev becomes a well/minimum in V_h)
    if particle_type.lower() in ("hole", "hh", "lh", "h"):
        V = -np.asarray(potential_ev, dtype=float)
        is_hole = True
    else:
        V = np.asarray(potential_ev, dtype=float)
        is_hole = False
    
    # Half-step interface masses via harmonic mean (ensures exact flux continuity)
    # 1 / m_{i+1/2} = 0.5 * (1/m_i + 1/m_{i+1})
    inv_m = 1.0 / m_eff
    inv_m_half = 0.5 * (inv_m[:-1] + inv_m[1:])  # Length N-1
    
    # Kinetic coupling coefficient: t_{i+1/2} = (ħ² / (2 m0 q dz²)) * (1 / m_{i+1/2})
    t_half = (HBAR2_OVER_2M0_EV_NM2 / (dz ** 2)) * inv_m_half
    
    # Assemble symmetric tridiagonal Hamiltonian:
    # Sub/Super-diagonal: H_{i, i+1} = -t_{i+1/2}
    # Main diagonal: H_{i, i} = V_i + t_{i-1/2} + t_{i+1/2}
    # Using Dirichlet boundary conditions at domain ends
    sub_d = -t_half[1:-1]  # Length N-3 for interior nodes
    
    main_d = np.zeros(N - 2, dtype=float)
    for i in range(1, N - 1):
        idx = i - 1
        t_left = t_half[i - 1]
        t_right = t_half[i]
        main_d[idx] = V[i] + t_left + t_right
    
    # Request eigenvalues
    k = min(num_states, len(main_d) - 1)
    if k < 1:
        return np.array([]), np.empty((0, N)), np.empty((0, N)), {}
    
    try:
        # Solve lowest k eigenvalues of real symmetric tridiagonal matrix in O(N*k)
        w_int, v_int = la.eigh_tridiagonal(
            main_d, sub_d, select='i', select_range=(0, k - 1)
        )
    except Exception as e:
        # Fallback to standard dense eigh if tridiagonal solver encounters singular limits
        H_dense = np.diag(main_d) + np.diag(sub_d, 1) + np.diag(sub_d, -1)
        w_all, v_all = la.eigh(H_dense)
        w_int = w_all[:k]
        v_int = v_all[:, :k]
    
    # Reconstruct full-grid wavefunctions with zero boundary conditions
    energies_list = []
    wf_list = []
    prob_list = []
    norm_errors = {}
    
    # Barrier height for bound-state verification
    v_barrier = None
    if barrier_ref_ev is not None:
        v_barrier = -barrier_ref_ev if is_hole else barrier_ref_ev
    else:
        v_span = float(np.max(V) - np.min(V))
        if v_span < 1e-4:
            # Uniform potential (infinite square well with Dirichlet boundary walls)
            v_barrier = np.inf
        else:
            # Infer barrier height from domain edges
            edge_len = max(2, N // 10)
            v_barrier = min(np.mean(V[:edge_len]), np.mean(V[-edge_len:]))
    
    for state_idx in range(len(w_int)):
        e_val = w_int[state_idx]
        psi_int = v_int[:, state_idx]
        
        # Build full domain wavefunction
        psi_full = np.zeros(N, dtype=float)
        psi_full[1:-1] = psi_int
        
        # Numerical normalization: \int |\psi|^2 dz = 1
        norm_factor = np.sqrt(np.sum(psi_full ** 2) * dz)
        if norm_factor > 1e-12:
            psi_full /= norm_factor
        
        # Deterministic phase/sign convention:
        # Anchor sign so that the primary extremum (maximum absolute amplitude) is positive.
        # This prevents random sign flipping between iterations and across voltage steps.
        max_abs_idx = np.argmax(np.abs(psi_full))
        if psi_full[max_abs_idx] < 0:
            psi_full = -psi_full
        
        prob_density = psi_full ** 2
        integral_check = np.sum(prob_density) * dz
        residual_err = abs(1.0 - integral_check)
        
        # Bound-state validation:
        # Check if state is localized inside well (or within barrier threshold)
        is_bound = True
        if v_barrier is not None:
            # Allow small numerical margin for shallow wells
            if e_val > v_barrier + 0.05:
                is_bound = False
        
        if well_mask is not None and np.any(well_mask):
            confinement_prob = np.sum(prob_density[well_mask]) * dz
            # A true bound state must have substantial probability in or near the well
            if confinement_prob < 0.15 and e_val > v_barrier:
                is_bound = False
        
        if not is_bound and state_idx > 0:
            # Skip unconfined continuum states if we already have bound states
            continue
            
        # Physical energy level (for holes, convert back to valence band energy: E_h = -E_eigen)
        phys_energy = -e_val if is_hole else e_val
        
        energies_list.append(phys_energy)
        wf_list.append(psi_full)
        prob_list.append(prob_density)
        norm_errors[f"{particle_type}_{state_idx + 1}"] = float(residual_err)
    
    if not energies_list:
        return np.array([]), np.empty((0, N)), np.empty((0, N)), {}
    
    energies = np.array(energies_list)
    wavefunctions = np.array(wf_list)
    probabilities = np.array(prob_list)
    
    return energies, wavefunctions, probabilities, norm_errors


# =============================================================================
# OPTICAL TRANSITIONS & OVERLAP INTEGRALS
# =============================================================================
def compute_optical_transitions(
    z_nm: np.ndarray,
    electron_energies: np.ndarray,
    electron_wf: np.ndarray,
    hole_energies: np.ndarray,
    hole_wf: np.ndarray,
    hole_types: List[str]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    """
    Computes optical transition matrix elements between all electron and hole subbands:
      - Wavefunction overlap integral: Gamma_ij = |int psi_e,i*(z) psi_h,j(z) dz|^2
      - Transition energy: E_trans = E_e,i - E_h,j
      - Emission wavelength: lambda = hc / E_trans
      
    Returns:
      trans_energies: 2D array [num_e, num_h] in eV.
      trans_wavelengths: 2D array [num_e, num_h] in nm.
      overlaps: 2D array [num_e, num_h] (dimensionless, in [0, 1]).
      dominant_list: Sorted list of key optical transitions.
    """
    dz = z_nm[1] - z_nm[0]
    num_e = len(electron_energies)
    num_h = len(hole_energies)
    
    if num_e == 0 or num_h == 0:
        return np.empty((0, 0)), np.empty((0, 0)), np.empty((0, 0)), []
    
    trans_energies = np.zeros((num_e, num_h), dtype=float)
    trans_wavelengths = np.zeros((num_e, num_h), dtype=float)
    overlaps = np.zeros((num_e, num_h), dtype=float)
    dominant_list = []
    
    for i in range(num_e):
        e_energy = electron_energies[i]
        psi_e = electron_wf[i]
        
        for j in range(num_h):
            h_energy = hole_energies[j]
            psi_h = hole_wf[j]
            htype = hole_types[j] if j < len(hole_types) else f"h{j+1}"
            
            # Transition energy in eV (conduction band level minus valence band level)
            delta_e = e_energy - h_energy
            trans_energies[i, j] = delta_e
            
            # Emission wavelength in nm: lambda = hc / delta_E
            if delta_e > 0.05:
                wlen = HC_EV_NM / delta_e
            else:
                wlen = 0.0
            trans_wavelengths[i, j] = wlen
            
            # Electron-hole envelope overlap integral
            overlap_integral = np.sum(psi_e * psi_h) * dz
            gamma = float(overlap_integral ** 2)
            gamma = np.clip(gamma, 0.0, 1.0)
            overlaps[i, j] = gamma
            
            label = f"e{i+1}-{htype}"
            dominant_list.append({
                "name": label,
                "e_index": i + 1,
                "h_index": j + 1,
                "h_type": htype,
                "energy_ev": float(delta_e),
                "wavelength_nm": float(wlen),
                "overlap": float(gamma),
                "oscillator_strength": float(gamma * (delta_e / 1.5)) # Relative scaling
            })
    
    # Sort transitions by overlap integral (strength of optical matrix element)
    dominant_list.sort(key=lambda x: x["overlap"], reverse=True)
    return trans_energies, trans_wavelengths, overlaps, dominant_list


# =============================================================================
# CARRIER DENSITY & FERMI-DIRAC POPULATION
# =============================================================================
def compute_quantum_carrier_density(
    z_nm: np.ndarray,
    electron_energies: np.ndarray,
    electron_wf: np.ndarray,
    hole_energies: np.ndarray,
    hole_wf: np.ndarray,
    temp_k: float = 300.0,
    ef_electron: Optional[float] = None,
    ef_hole: Optional[float] = None,
    m_dos_e: float = 0.067,
    m_dos_h: float = 0.45
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculates 2D sheet carrier populations and 3D volumetric quantum carrier densities:
      - 2D Density of States: D_2D = (m* / (pi ħ²))
      - 2D Subband Sheet Carrier Density:
          n_2D,n = D_2D * k_B T * ln( 1 + exp( (E_Fn - E_n) / k_B T ) )
      - 3D Volumetric Quantum Carrier Density:
          n_qw(z) = sum_n n_2D,n * |psi_e,n(z)|^2
          
    Returns:
      sheet_n: 1D array of 2D electron sheet densities (cm^-2).
      sheet_p: 1D array of 2D hole sheet densities (cm^-2).
      density_n: 1D array of 3D electron volumetric density (cm^-3).
      density_p: 1D array of 3D hole volumetric density (cm^-3).
    """
    N = len(z_nm)
    num_e = len(electron_energies)
    num_h = len(hole_energies)
    
    density_n = np.zeros(N, dtype=float)
    density_p = np.zeros(N, dtype=float)
    sheet_n = np.zeros(num_e, dtype=float)
    sheet_p = np.zeros(num_h, dtype=float)
    
    if num_e == 0 and num_h == 0:
        return sheet_n, sheet_p, density_n, density_p
    
    vt_ev = (K_BOLTZ * temp_k) / Q_ELEM  # Thermal voltage in eV (~0.02585 eV at 300K)
    
    # 2D DOS prefactor: m* k_B T / (pi ħ²) in cm^-2
    # D_2D_pref = (m_r * M_ELEM * K_BOLTZ * T) / (pi * H_BAR^2) * (1e-4 m^2/cm^2)
    # Prefactor constant: M_ELEM * K_BOLTZ / (pi * H_BAR^2) * 1e-4 = 3.6167e10 cm^-2 K^-1
    dos_pref_per_k = (M_ELEM * K_BOLTZ) / (np.pi * (H_BAR ** 2)) * 1e-4
    
    # 1. Electron sheet & volumetric populations
    if num_e > 0:
        if ef_electron is None:
            # Default reference: Quasi-Fermi level near lowest subband
            ef_e = electron_energies[0] - 0.05
        else:
            ef_e = ef_electron
            
        d2d_e = m_dos_e * dos_pref_per_k * temp_k
        for i in range(num_e):
            e_level = electron_energies[i]
            eta = (ef_e - e_level) / vt_ev
            
            # Numerically stable Fermi-Dirac integral log(1 + exp(eta))
            if eta > 40.0:
                fd_term = eta
            elif eta < -40.0:
                fd_term = np.exp(eta)
            else:
                fd_term = np.log(1.0 + np.exp(eta))
            
            n_sheet_i = d2d_e * fd_term  # in cm^-2
            sheet_n[i] = n_sheet_i
            
            # Wavefunction probability density has units nm^-1 = 1e7 cm^-1
            # n_3D(z) = n_2D * |psi(z)|^2 * (1e7 cm^-1) -> cm^-3
            prob_cm = (electron_wf[i] ** 2) * 1e7
            density_n += n_sheet_i * prob_cm
    
    # 2. Hole sheet & volumetric populations
    if num_h > 0:
        if ef_hole is None:
            ef_h = hole_energies[0] + 0.05
        else:
            ef_h = ef_hole
            
        d2d_h = m_dos_h * dos_pref_per_k * temp_k
        for j in range(num_h):
            h_level = hole_energies[j]
            # For holes in valence band: occupation increases as energy drops below E_Fp
            eta = (h_level - ef_h) / vt_ev
            
            if eta > 40.0:
                fd_term = eta
            elif eta < -40.0:
                fd_term = np.exp(eta)
            else:
                fd_term = np.log(1.0 + np.exp(eta))
            
            p_sheet_j = d2d_h * fd_term
            sheet_p[j] = p_sheet_j
            
            prob_cm = (hole_wf[j] ** 2) * 1e7
            density_p += p_sheet_j * prob_cm
            
    return sheet_n, sheet_p, density_n, density_p


# =============================================================================
# HIGH-LEVEL DEVICE SIMULATOR INTEGRATION & DISPATCHER
# =============================================================================
def solve_quantum_well(
    structure: Any = None,
    band_profile: Optional[Dict[str, np.ndarray]] = None,
    layers: Optional[List[Dict[str, Any]]] = None,
    temperature_k: float = 300.0,
    num_electron_states: int = 3,
    num_hole_states: int = 3,
    coupling_mode: str = "Coupled MQW",
    mat_system: str = "Zincblende",
    ef_electron: Optional[float] = None,
    ef_hole: Optional[float] = None,
    electric_field_v_cm: float = 0.0,
    provenance_dict: Optional[Dict[str, Any]] = None
) -> QuantumWellResult:
    """
    Master programmatic entry point for solving quantum well confined states.
    
    Accepts:
      - structure: Aestimo Structure / StructureFrom model object.
      - band_profile: Dictionary with 'z' (nm), 'ec' (eV), 'ev' (eV).
      - layers: List of layer dictionaries (optional, auto-extracted if omitted).
      - temperature_k: Temperature in Kelvin (default 300.0).
      - num_electron_states: Number of conduction subbands to extract (default 3).
      - num_hole_states: Number of valence subbands to extract (default 3).
      - coupling_mode: 'Coupled MQW' (full multi-well solver) or 'Isolated Wells'.
      
    Returns:
      QuantumWellResult structured object.
    """
    # 1. Resolve mesh coordinates and band potentials
    z_nm = None
    ec_ev = None
    ev_ev = None
    m_e_arr = None
    m_hh_arr = None
    m_lh_arr = None
    well_regions = []
    
    if band_profile is not None and "z" in band_profile:
        z_nm = np.asarray(band_profile["z"], dtype=float)
        ec_ev = np.asarray(band_profile["ec"], dtype=float)
        ev_ev = np.asarray(band_profile["ev"], dtype=float)
        
    elif structure is not None:
        # Extract from Aestimo model object
        dx_nm = getattr(structure, "dx", 0.1e-9) * 1e9
        n_max = getattr(structure, "n_max", 100)
        z_nm = np.arange(n_max) * dx_nm
        
        # Conduction & valence potentials
        if hasattr(structure, "fi_e") and len(structure.fi_e) == n_max:
            ec_ev = structure.fi_e / Q_ELEM
        elif hasattr(structure, "Ec_result"):
            ec_ev = np.asarray(structure.Ec_result)
        else:
            ec_ev = np.zeros(n_max)
            
        if hasattr(structure, "fi_h") and len(structure.fi_h) == n_max:
            ev_ev = structure.fi_h / Q_ELEM
        elif hasattr(structure, "Ev_result"):
            ev_ev = np.asarray(structure.Ev_result)
        else:
            ev_ev = ec_ev - 1.424  # Default GaAs gap
            
        if hasattr(structure, "cb_meff") and len(structure.cb_meff) == n_max:
            m_e_arr = structure.cb_meff / M_ELEM
            
    elif layers is not None and len(layers) > 0:
        # Synthesize spatial mesh and band potentials directly from layer definitions
        tot_thick = sum(float(l.get("thickness", 10.0)) for l in layers)
        min_thick = min(float(l.get("thickness", 10.0)) for l in layers)
        dz_nm = max(0.05, min(0.2, min_thick / 30.0))
        num_pts = max(100, int(round(tot_thick / dz_nm)) + 1)
        z_nm = np.linspace(0.0, tot_thick, num_pts)
        ec_ev = np.zeros(num_pts, dtype=float)
        ev_ev = np.zeros(num_pts, dtype=float)
        
        z_accum = 0.0
        for l in layers:
            th = float(l.get("thickness", 10.0))
            z0 = z_accum
            z1 = z_accum + th
            z_accum = z1
            
            mask = (z_nm >= z0) & (z_nm <= z1) if np.isclose(z1, tot_thick) else (z_nm >= z0) & (z_nm < z1)
            if "Ec" in l and "Ev" in l:
                ec_val = float(l["Ec"])
                ev_val = float(l["Ev"])
            else:
                from aeslibs.structure_diagram import extract_layer_properties
                props = extract_layer_properties(l, temp_k=temperature_k, mat_system=mat_system)
                ec_val = props["Ec"]
                ev_val = props["Ev"]
            ec_ev[mask] = ec_val
            ev_ev[mask] = ev_val

        if electric_field_v_cm != 0.0:
            field_ev_per_nm = electric_field_v_cm * 1e-7
            # If well region(s) exist, apply field across the well region(s) and offset barriers
            # to preserve bound state integrity (prevents unphysical barrier edge drops in finite domain)
            z_accum_temp = 0.0
            well_starts = []
            well_ends = []
            for l in layers:
                th = float(l.get("thickness", 10.0))
                ltype = str(l.get("type", "barrier")).lower()
                if "well" in ltype or ltype == "w":
                    well_starts.append(z_accum_temp)
                    well_ends.append(z_accum_temp + th)
                z_accum_temp += th
            
            if well_starts:
                z_w0 = min(well_starts)
                z_w1 = max(well_ends)
                z_center = 0.5 * (z_w0 + z_w1)
                z_clamped = np.clip(z_nm, z_w0, z_w1)
                delta_v = field_ev_per_nm * (z_clamped - z_center)
                ec_ev += delta_v
                ev_ev += delta_v
            else:
                z_center = 0.5 * tot_thick
                ec_ev += field_ev_per_nm * (z_nm - z_center)
                ev_ev += field_ev_per_nm * (z_nm - z_center)

    if z_nm is None or ec_ev is None or ev_ev is None:
        raise ValueError("Cannot solve quantum well: missing band profile or structure coordinates.")
        
    N = len(z_nm)
    
    # 2. Extract material parameters and well regions
    # If effective mass profiles are not provided, synthesize from layers database
    if layers is not None and len(layers) > 0:
        # Build spatial maps from layer definitions
        m_e_synth = np.zeros(N, dtype=float)
        m_hh_synth = np.zeros(N, dtype=float)
        m_lh_synth = np.zeros(N, dtype=float)
        well_mask = np.zeros(N, dtype=bool)
        
        z_accum = 0.0
        for l_idx, l in enumerate(layers):
            th = float(l.get("thickness", 10.0))
            z0 = z_accum
            z1 = z_accum + th
            z_accum = z1
            
            mat = str(l.get("material", "GaAs"))
            x = float(l.get("mole", 0.0))
            y = float(l.get("mole_y", 0.0))
            ltype = str(l.get("type", "barrier")).lower()
            is_well = "well" in ltype or ltype == "w"
            
            # Lookup effective masses from database or layer properties
            if "m_e" in l and "m_hh" in l:
                m_e_val = float(l["m_e"])
                m_hh_val = float(l["m_hh"])
                m_lh_val = float(l.get("m_lh", 0.082))
            else:
                from aeslibs.structure_diagram import extract_layer_properties
                props = extract_layer_properties(l, temp_k=temperature_k, mat_system=mat_system)
                m_e_val = props["m_e"]
                m_hh_val = props["m_hh"]
                m_lh_val = props["m_lh"]
                
            mask_idx = (z_nm >= z0 - 1e-6) & (z_nm <= z1 + 1e-6)
            m_e_synth[mask_idx] = m_e_val
            m_hh_synth[mask_idx] = m_hh_val
            m_lh_synth[mask_idx] = m_lh_val
            m_lh_synth[mask_idx] = m_lh_val
            if is_well:
                well_mask[mask_idx] = True
                well_regions.append({
                    "layer_index": l_idx,
                    "material": mat,
                    "thickness_nm": th,
                    "z0_nm": z0,
                    "z1_nm": z1
                })
                
        if m_e_arr is None:
            m_e_arr = m_e_synth
        m_hh_arr = m_hh_synth
        m_lh_arr = m_lh_synth
    else:
        # Fallback defaults for standard III-V heterostructures
        if m_e_arr is None:
            m_e_arr = np.full(N, 0.067)
        m_hh_arr = np.full(N, 0.45)
        m_lh_arr = np.full(N, 0.082)
        well_mask = np.ones(N, dtype=bool)

    # 3. Apply optional electric field tilt (for QCSE and Stark shift simulations)
    ec_eff = ec_ev.copy()
    ev_eff = ev_ev.copy()
    if abs(electric_field_v_cm) > 1e-3:
        # V_field(z) = - q * F * z
        # field in V/cm -> V/nm: F_v_nm = F_v_cm * 1e-7
        f_v_nm = electric_field_v_cm * 1e-7
        z_mid = (z_nm[0] + z_nm[-1]) / 2.0
        v_tilt = f_v_nm * (z_nm - z_mid)
        ec_eff += v_tilt
        ev_eff += v_tilt

    # 4. Solve Electron Confined States (Conduction Band)
    e_energies, e_wf, e_prob, e_norms = solve_effective_mass_schrodinger_1d(
        z_nm=z_nm,
        potential_ev=ec_eff,
        m_eff_rel=m_e_arr,
        num_states=num_electron_states,
        particle_type="electron",
        well_mask=well_mask
    )

    # 5. Solve Hole Confined States (Heavy Holes & Light Holes)
    # Heavy holes
    hh_energies, hh_wf, hh_prob, hh_norms = solve_effective_mass_schrodinger_1d(
        z_nm=z_nm,
        potential_ev=ev_eff,
        m_eff_rel=m_hh_arr,
        num_states=num_hole_states,
        particle_type="hole",
        well_mask=well_mask
    )

    # Light holes
    lh_energies, lh_wf, lh_prob, lh_norms = solve_effective_mass_schrodinger_1d(
        z_nm=z_nm,
        potential_ev=ev_eff,
        m_eff_rel=m_lh_arr,
        num_states=max(1, num_hole_states // 2),
        particle_type="hole",
        well_mask=well_mask
    )

    # Combine and order hole states energetically (highest valence level first)
    all_hole_energies = []
    all_hole_wf = []
    all_hole_prob = []
    all_hole_types = []
    all_norms = {**e_norms}

    for idx, eh in enumerate(hh_energies):
        all_hole_energies.append(eh)
        all_hole_wf.append(hh_wf[idx])
        all_hole_prob.append(hh_prob[idx])
        all_hole_types.append(f"hh{idx + 1}")
        all_norms[f"hh_{idx + 1}"] = hh_norms.get(f"hole_{idx + 1}", 0.0)

    for idx, el in enumerate(lh_energies):
        all_hole_energies.append(el)
        all_hole_wf.append(lh_wf[idx])
        all_hole_prob.append(lh_prob[idx])
        all_hole_types.append(f"lh{idx + 1}")
        all_norms[f"lh_{idx + 1}"] = lh_norms.get(f"hole_{idx + 1}", 0.0)

    if all_hole_energies:
        # Sort descending in energy (top of valence band is highest energy in eV)
        sort_order = np.argsort(all_hole_energies)[::-1]
        sorted_h_energies = np.array(all_hole_energies)[sort_order][:num_hole_states]
        sorted_h_wf = np.array(all_hole_wf)[sort_order][:num_hole_states]
        sorted_h_prob = np.array(all_hole_prob)[sort_order][:num_hole_states]
        sorted_h_types = [all_hole_types[k] for k in sort_order][:num_hole_states]
    else:
        sorted_h_energies = np.empty(0)
        sorted_h_wf = np.empty((0, N))
        sorted_h_prob = np.empty((0, N))
        sorted_h_types = []

    # 6. Calculate Optical Transitions
    trans_energies, trans_wavelengths, overlaps, dominant_trans = compute_optical_transitions(
        z_nm=z_nm,
        electron_energies=e_energies,
        electron_wf=e_wf,
        hole_energies=sorted_h_energies,
        hole_wf=sorted_h_wf,
        hole_types=sorted_h_types
    )

    # 7. Calculate 2D Sheet & 3D Volumetric Carrier Densities
    sheet_n, sheet_p, dens_n, dens_p = compute_quantum_carrier_density(
        z_nm=z_nm,
        electron_energies=e_energies,
        electron_wf=e_wf,
        hole_energies=sorted_h_energies,
        hole_wf=sorted_h_wf,
        temp_k=temperature_k,
        ef_electron=ef_electron,
        ef_hole=ef_hole,
        m_dos_e=np.mean(m_e_arr),
        m_dos_h=np.mean(m_hh_arr)
    )

    result = QuantumWellResult(
        z=z_nm,
        ec_profile=ec_eff,
        ev_profile=ev_eff,
        electron_energies=e_energies,
        electron_wavefunctions=e_wf,
        electron_probability=e_prob,
        hole_energies=sorted_h_energies,
        hole_wavefunctions=sorted_h_wf,
        hole_probability=sorted_h_prob,
        hole_types=sorted_h_types,
        transition_energies=trans_energies,
        transition_wavelengths=trans_wavelengths,
        overlap_integrals=overlaps,
        dominant_transitions=dominant_trans,
        sheet_densities_e=sheet_n,
        sheet_densities_h=sheet_p,
        quantum_n=dens_n,
        quantum_p=dens_p,
        normalization_errors=all_norms,
        qw_regions=well_regions,
        convergence={
            "status": "CONVERGED",
            "iterations": 1,
            "max_residual_ev": 0.0,
            "num_e_found": len(e_energies),
            "num_h_found": len(sorted_h_energies)
        },
        temperature_k=temperature_k,
        field_v_cm=electric_field_v_cm
    )
    return result


# =============================================================================
# SELF-CONSISTENT QW <-> POISSON COUPLING ENGINE
# =============================================================================
def solve_self_consistent_qw_poisson(
    z_nm: np.ndarray,
    initial_ec: np.ndarray,
    initial_ev: np.ndarray,
    dielectric_rel: np.ndarray,
    doping_profile_cm3: np.ndarray,
    layers: Optional[List[Dict[str, Any]]] = None,
    temperature_k: float = 300.0,
    max_iterations: int = 20,
    tolerance_ev: float = 1e-4,
    damping_factor: float = 0.2,
    num_e_states: int = 3,
    num_h_states: int = 3,
    ef_electron: Optional[float] = None,
    ef_hole: Optional[float] = None
) -> QuantumWellResult:
    """
    Executes a physically self-consistent Schrödinger-Poisson iterative loop:
      V^(k)(z) -> (E_n, psi_n) -> (n_qw, p_qw) -> rho(z) -> Poisson -> V_new -> V^(k+1)
      
    Under-relaxation damping factor alpha (typically 0.1 - 0.3) prevents numerical charge sloshing.
    """
    N = len(z_nm)
    dz_m = (z_nm[1] - z_nm[0]) * 1e-9
    eps_profile = np.asarray(dielectric_rel, dtype=float) * EPSILON_0
    
    current_ec = initial_ec.copy()
    current_ev = initial_ev.copy()
    eg_profile = initial_ec - initial_ev
    
    prev_energies = None
    final_result = None
    converged = False
    actual_iter = 0
    max_delta_e = 0.0
    
    for iteration in range(1, max_iterations + 1):
        actual_iter = iteration
        
        # 1. Solve Schrödinger eigenvalue problem on current potential profile
        band_prof = {"z": z_nm, "ec": current_ec, "ev": current_ev}
        qw_res = solve_quantum_well(
            band_profile=band_prof,
            layers=layers,
            temperature_k=temperature_k,
            num_electron_states=num_e_states,
            num_hole_states=num_h_states,
            ef_electron=ef_electron,
            ef_hole=ef_hole
        )
        final_result = qw_res
        
        # Check energy level convergence
        curr_energies = qw_res.electron_energies
        if prev_energies is not None and len(prev_energies) == len(curr_energies):
            max_delta_e = float(np.max(np.abs(curr_energies - prev_energies)))
            if max_delta_e < tolerance_ev:
                converged = True
                break
        prev_energies = curr_energies.copy()
        
        # 2. Formulate total volumetric space charge density:
        # rho(z) = q * [ p_qw(z) - n_qw(z) + N_D(z) - N_A(z) ] (C/m^3)
        n_m3 = qw_res.quantum_n * 1e6
        p_m3 = qw_res.quantum_p * 1e6
        dop_m3 = doping_profile_cm3 * 1e6
        rho = Q_ELEM * (p_m3 - n_m3 + dop_m3)
        
        # 3. Solve 1D Poisson equation: d/dz [ eps(z) d phi / dz ] = -rho(z)
        # Using tridiagonal Poisson matrix with fixed boundary potentials (Dirichlet)
        main_diag = np.zeros(N - 2, dtype=float)
        sub_diag = np.zeros(N - 3, dtype=float)
        super_diag = np.zeros(N - 3, dtype=float)
        rhs = np.zeros(N - 2, dtype=float)
        
        for i in range(1, N - 1):
            idx = i - 1
            eps_left = 0.5 * (eps_profile[i - 1] + eps_profile[i])
            eps_right = 0.5 * (eps_profile[i] + eps_profile[i + 1])
            
            main_diag[idx] = -(eps_left + eps_right) / (dz_m ** 2)
            if idx > 0:
                sub_diag[idx - 1] = eps_left / (dz_m ** 2)
            if idx < N - 3:
                super_diag[idx] = eps_right / (dz_m ** 2)
                
            rhs[idx] = -rho[i]
            
        # Incorporate boundary conditions into RHS
        phi_left = -current_ec[0]   # V
        phi_right = -current_ec[-1] # V
        rhs[0] -= (0.5 * (eps_profile[0] + eps_profile[1]) / (dz_m ** 2)) * phi_left
        rhs[-1] -= (0.5 * (eps_profile[-2] + eps_profile[-1]) / (dz_m ** 2)) * phi_right
        
        A_poisson = np.diag(main_diag) + np.diag(sub_diag, -1) + np.diag(super_diag, 1)
        phi_interior = la.solve(A_poisson, rhs)
        
        phi_new = np.zeros(N, dtype=float)
        phi_new[0] = phi_left
        phi_new[-1] = phi_right
        phi_new[1:-1] = phi_interior
        
        # Conduction band: Ec = -q * phi (in eV)
        ec_new = -phi_new
        
        # 4. Under-relaxation mixing
        alpha = float(np.clip(damping_factor, 0.05, 0.5))
        current_ec = (1.0 - alpha) * current_ec + alpha * ec_new
        current_ev = current_ec - eg_profile
        
    final_result.convergence = {
        "status": "CONVERGED" if converged else "MAX_ITERATIONS_REACHED",
        "iterations": actual_iter,
        "max_residual_ev": float(max_delta_e),
        "tolerance_ev": tolerance_ev
    }
    return final_result
