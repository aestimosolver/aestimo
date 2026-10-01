#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Aestimo 1D - Publication-Quality Structure Diagram & Band Alignment Module
aeslibs/structure_diagram.py

Provides scientifically rigorous, peer-reviewed-standard cross-sectional layer-stack
schematics, dynamic thickness scaling (TRUE, SCHEMATIC, AUTO), coupled energy band
alignments, periodic MQW grouping, device-specific physics annotations (LED, Laser, Solar),
interactive layer inspection, model consistency validation, and publication-ready vector
and raster export (PDF, SVG, 300+ DPI PNG).
"""

import os
import re
import math
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.patheffects as path_effects

import database

# Physical Constants
Q_E = 1.602176634e-19       # C (J/eV)
EPS_0 = 8.8541878128e-12    # F/m
M_E0 = 9.1093837015e-31     # kg

# ---------------------------------------------------------------------------
# 1. Chemical Formula & Text Formatting
# ---------------------------------------------------------------------------

SUB_MAP = str.maketrans("0123456789.", "₀₁₂₃₄₅₆₇₈₉.")

def format_chemical_formula(material, mole=0.0, mole_y=0.0, use_mathtext=True):
    """
    Formats semiconductor chemical formulas into publication-grade typography.
    e.g., AlGaAs (x=0.6) -> Al₀.₆Ga₀.₄As or $\\mathrm{Al}_{0.60}\\mathrm{Ga}_{0.40}\\mathrm{As}$
          InGaN (x=0.15) -> In₀.₁₅Ga₀.₈₅N or $\\mathrm{In}_{0.15}\\mathrm{Ga}_{0.85}\\mathrm{N}$
          InGaAsP (x=0.7, y=0.82) -> $\\mathrm{In}_{0.70}\\mathrm{Ga}_{0.30}\\mathrm{As}_{0.82}\\mathrm{P}_{0.18}$
    """
    mat = str(material).strip()
    x = float(mole) if mole is not None else 0.0
    y = float(mole_y) if mole_y is not None else 0.0

    # Binary materials or x=0
    if mat in database.materialproperty or (x == 0.0 and y == 0.0):
        if use_mathtext:
            return f"$\\mathrm{{{mat}}}$"
        return mat

    # Ternary alloys
    if mat in database.alloyproperty:
        alloy_info = database.alloyproperty[mat]
        m1 = alloy_info.get("Material1", "")
        m2 = alloy_info.get("Material2", "")

        # Standard naming patterns
        if mat == "AlGaAs":
            if use_mathtext:
                return f"$\\mathrm{{Al}}_{{{x:.2f}}}\\mathrm{{Ga}}_{{{1.0-x:.2f}}}\\mathrm{{As}}$"
            return f"Al{x:.2f}Ga{1.0-x:.2f}As"
        elif mat == "InGaN":
            if use_mathtext:
                return f"$\\mathrm{{In}}_{{{x:.2f}}}\\mathrm{{Ga}}_{{{1.0-x:.2f}}}\\mathrm{{N}}$"
            return f"In{x:.2f}Ga{1.0-x:.2f}N"
        elif mat == "InGaAs":
            if use_mathtext:
                return f"$\\mathrm{{In}}_{{{x:.2f}}}\\mathrm{{Ga}}_{{{1.0-x:.2f}}}\\mathrm{{As}}$"
            return f"In{x:.2f}Ga{1.0-x:.2f}As"
        elif mat == "AlGaN":
            if use_mathtext:
                return f"$\\mathrm{{Al}}_{{{x:.2f}}}\\mathrm{{Ga}}_{{{1.0-x:.2f}}}\\mathrm{{N}}$"
            return f"Al{x:.2f}Ga{1.0-x:.2f}N"
        elif mat == "InGaP":
            if use_mathtext:
                return f"$\\mathrm{{In}}_{{{x:.2f}}}\\mathrm{{Ga}}_{{{1.0-x:.2f}}}\\mathrm{{P}}$"
            return f"In{x:.2f}Ga{1.0-x:.2f}P"
        elif mat == "GaAsP":
            if use_mathtext:
                return f"$\\mathrm{{GaAs}}_{{{1.0-x:.2f}}}\\mathrm{{P}}_{{{x:.2f}}}$"
            return f"GaAs{1.0-x:.2f}P{x:.2f}"
        else:
            if use_mathtext:
                return f"$\\mathrm{{{mat}}}\\ (x={x:.2f})$"
            return f"{mat} (x={x:.2f})"

    # Quaternary alloys
    if mat in database.alloyproperty4 or mat == "InGaAsP":
        if use_mathtext:
            return f"$\\mathrm{{In}}_{{{x:.2f}}}\\mathrm{{Ga}}_{{{1.0-x:.2f}}}\\mathrm{{As}}_{{{y:.2f}}}\\mathrm{{P}}_{{{1.0-y:.2f}}}$"
        return f"In{x:.2f}Ga{1.0-x:.2f}As{y:.2f}P{1.0-y:.2f}"

    return mat


def format_thickness_str(d_nm):
    """Formats thickness with appropriate unit (nm or um)."""
    if d_nm >= 1000.0:
        return f"{d_nm / 1000.0:.2f} um"
    elif d_nm >= 10.0:
        return f"{d_nm:.1f} nm"
    else:
        return f"{d_nm:.2f} nm"


def format_doping_str(doping_val, doping_type):
    """Formats doping type and concentration."""
    dtype = str(doping_type).strip().lower()
    val = float(doping_val) if doping_val is not None else 0.0

    if dtype in ("i", "undoped") or val <= 0:
        return "undoped (i)"
    
    # Scientific formatting e.g. 2.0e18 -> 2.0e18 cm^-3
    exp = int(math.floor(math.log10(val)))
    mantissa = val / (10 ** exp)
    if abs(mantissa - round(mantissa)) < 0.05:
        mant_str = f"{round(mantissa):d}"
    else:
        mant_str = f"{mantissa:.1f}"
    
    return f"{dtype} = {mant_str}e{exp} cm^-3"


# ---------------------------------------------------------------------------
# 2. Material & Electronic Property Calculation
# ---------------------------------------------------------------------------

def calculate_varshni_eg(matprops, temp_k=300.0):
    """Calculates temperature-dependent bandgap using Varshni equation."""
    eg0 = matprops.get("Eg", 1.42)
    alpha = matprops.get("alpha_varshni", 0.0)
    beta = matprops.get("beta_varshni", 0.0)
    if alpha > 0 and beta > 0:
        return eg0 - (alpha * (temp_k ** 2)) / (beta + temp_k)
    return eg0


_LAYER_PROPERTY_CACHE = {}

def extract_layer_properties(layer, temp_k=300.0, mat_system="Zincblende", provenance_dict=None):
    """
    Extracts complete electronic and optical properties for a layer from database.py.
    Returns dictionary with Eg, offsets, dielectric constants, effective masses,
    refractive index, and 4-tier provenance metadata.
    Uses in-memory cache for ultra-fast lookup.
    """
    mat = str(layer.get("material", "GaAs")).strip()
    x = float(layer.get("mole", 0.0))
    y = float(layer.get("mole_y", 0.0))
    thickness = float(layer.get("thickness", 10.0))
    doping = float(layer.get("doping", 0.0))
    doping_type = str(layer.get("doping_type", "n")).strip().lower()
    layer_type = str(layer.get("type", "barrier")).strip().lower()

    # Check cache for intrinsic material parameters
    cache_key = (mat, round(x, 4), round(y, 4), round(temp_k, 1), mat_system)
    cached = _LAYER_PROPERTY_CACHE.get(cache_key)

    if cached is not None:
        eg = cached["Eg"]
        band_offset = cached["Band_offset"]
        ec = cached["Ec"]
        ev = cached["Ev"]
        m_e = cached["m_e"]
        m_hh = cached["m_hh"]
        m_lh = cached["m_lh"]
        eps_r = cached["eps_r"]
        eps_inf = cached["eps_inf"]
        refractive_index = cached["n_refractive"]
    else:
        # Defaults
        eg = 1.424
        band_offset = 0.65
        m_e = 0.067
        m_hh = 0.45
        m_lh = 0.082
        eps_r = 12.9
        eps_inf = 10.89
        refractive_index = math.sqrt(eps_inf) if eps_inf > 0 else 3.5

        # 1. Binary material lookup
        if mat in database.materialproperty:
            props = database.materialproperty[mat]
            eg = calculate_varshni_eg(props, temp_k)
            band_offset = props.get("Band_offset", 0.65)
            m_e = props.get("m_e", 0.067)
            m_hh = props.get("m_hh", 0.45)
            m_lh = props.get("m_lh", 0.082)
            eps_r = props.get("epsilonStatic", 12.9)
            eps_inf = props.get("epsilonHigh", props.get("epsilon_inf", eps_r * 0.85))
            refractive_index = math.sqrt(eps_inf)

        # 2. Ternary alloy lookup
        elif mat in database.alloyproperty:
            alloy = database.alloyproperty[mat]
            m1_name = alloy.get("Material1", "")
            m2_name = alloy.get("Material2", "")
            m1 = database.materialproperty.get(m1_name, {})
            m2 = database.materialproperty.get(m2_name, {})

            eg1 = calculate_varshni_eg(m1, temp_k)
            eg2 = calculate_varshni_eg(m2, temp_k)
            bowing = alloy.get("Bowing_param", 0.0)
            eg = x * eg1 + (1.0 - x) * eg2 - bowing * x * (1.0 - x)
            band_offset = alloy.get("Band_offset", 0.65)

            m_e = x * m1.get("m_e", 0.067) + (1.0 - x) * m2.get("m_e", 0.067)
            m_hh = x * m1.get("m_hh", 0.45) + (1.0 - x) * m2.get("m_hh", 0.45)
            m_lh = x * m1.get("m_lh", 0.082) + (1.0 - x) * m2.get("m_lh", 0.082)
            eps_r = x * m1.get("epsilonStatic", 12.9) + (1.0 - x) * m2.get("epsilonStatic", 12.9)
            eps1_inf = m1.get("epsilonHigh", m1.get("epsilonStatic", 12.9) * 0.85)
            eps2_inf = m2.get("epsilonHigh", m2.get("epsilonStatic", 12.9) * 0.85)
            eps_inf = x * eps1_inf + (1.0 - x) * eps2_inf
            refractive_index = math.sqrt(eps_inf)

        # 3. Quaternary alloy lookup
        elif mat in database.alloyproperty4 or mat == "InGaAsP":
            alloy = database.alloyproperty4.get(mat, {})
            m1 = database.materialproperty.get(alloy.get("Material1", "InAs"), {})
            m2 = database.materialproperty.get(alloy.get("Material2", "GaAs"), {})
            m3 = database.materialproperty.get(alloy.get("Material3", "InP"), {})
            m4 = database.materialproperty.get(alloy.get("Material4", "GaP"), {})

            eg_abc = x * m1.get("Eg", 0.354) + (1.0 - x) * m2.get("Eg", 1.424) - alloy.get("Bowing_param_ABC", 0.475) * x * (1.0 - x)
            eg_abd = x * m3.get("Eg", 1.344) + (1.0 - x) * m4.get("Eg", 2.26) - alloy.get("Bowing_param_ABD", 0.0) * x * (1.0 - x)
            denom = (x * (1.0 - x) + y * (1.0 - y))
            if denom > 1e-6:
                eg = (x * (1.0 - x) * (y * eg_abc + (1.0 - y) * eg_abd)) / denom
            else:
                eg = 0.80

            band_offset = alloy.get("Band_offset", 0.40)
            m_e = x * 0.04 + (1.0 - x) * 0.067
            m_hh = 0.45
            m_lh = 0.08
            eps_r = 13.0
            eps_inf = 11.0
            refractive_index = 3.45

        # Band edges (relative to reference in eV)
        ec = band_offset * eg
        ev = -(1.0 - band_offset) * eg

        _LAYER_PROPERTY_CACHE[cache_key] = {
            "Eg": eg,
            "Band_offset": band_offset,
            "Ec": ec,
            "Ev": ev,
            "m_e": m_e,
            "m_hh": m_hh,
            "m_lh": m_lh,
            "eps_r": eps_r,
            "eps_inf": eps_inf,
            "n_refractive": refractive_index
        }


    # Provenance lookup
    prov = "ASSUMED"
    if provenance_dict and isinstance(provenance_dict, dict):
        layer_prov = provenance_dict.get("layers", {})
        if isinstance(layer_prov, dict):
            prov = layer_prov.get(mat, provenance_dict.get("general", "ASSUMED"))
        elif isinstance(layer_prov, str):
            prov = layer_prov
        elif "parameters" in provenance_dict:
            prov = provenance_dict.get("parameters", {}).get("thickness", "EXPERIMENTAL")

    return {
        "material": mat,
        "mole": x,
        "mole_y": y,
        "thickness": thickness,
        "doping": doping,
        "doping_type": doping_type,
        "type": layer_type,
        "Eg": eg,
        "Band_offset": band_offset,
        "Ec": ec,
        "Ev": ev,
        "m_e": m_e,
        "m_hh": m_hh,
        "m_lh": m_lh,
        "eps_r": eps_r,
        "eps_inf": eps_inf,
        "n_refractive": refractive_index,
        "provenance": prov
    }


# ---------------------------------------------------------------------------
# 3. Layer Role Inference & Structure Classification
# ---------------------------------------------------------------------------

def infer_layer_role(layer, idx, total_layers, all_layers=None):
    """
    Infers the structural role of a layer based on its position, thickness,
    doping, and type (barrier/well).
    Roles: Substrate, Buffer, n-Cladding, p-Cladding, SCH Waveguide,
           Quantum Well (QW), Barrier, Electron Blocking Layer (EBL),
           Contact Cap, Window/Absorber.
    """
    ltype = str(layer.get("type", "barrier")).lower()
    d = float(layer.get("thickness", 10.0))
    dtype = str(layer.get("doping_type", "n")).lower()
    dop = float(layer.get("doping", 0.0))
    mat = str(layer.get("material", ""))
    x = float(layer.get("mole", 0.0))

    if ltype in ("well", "w"):
        return "Quantum Well (Active)"

    # Bottom layer with large thickness or substrate keyword
    if idx == total_layers - 1 or idx == 0:
        if d >= 200.0 or "sub" in mat.lower():
            if idx == total_layers - 1:
                return "Substrate / Buffer"
            else:
                return "Contact Cap / Window"

    # Electron Blocking Layer (EBL) characteristic: thin p-AlGaN or high-bandgap p-layer adjacent to active
    if dtype == "p" and ("AlGaN" in mat or (mat == "AlGaAs" and x >= 0.20)) and (10.0 <= d <= 40.0):
        return "Electron Blocking Layer (EBL)"

    # Waveguide / SCH (Separate Confinement Heterostructure)
    if "barrier" in ltype and (30.0 <= d <= 250.0) and dtype in ("i", "undoped"):
        return "SCH Waveguide Guide"

    # Cladding layers: thick, doped
    if d >= 100.0 and dop >= 1e17:
        if dtype == "p":
            return "p-Cladding"
        elif dtype == "n":
            return "n-Cladding"

    # Contact caps: top highly doped layer
    if idx == 0 and dop >= 1e18:
        return f"{dtype}⁺ Contact Cap"

    if ltype in ("barrier", "b"):
        return "Barrier"

    return "Semiconductor Layer"


# ---------------------------------------------------------------------------
# 4. Multi-Quantum-Well (MQW) Periodic Sequence Detection
# ---------------------------------------------------------------------------

def detect_mqw_regions(layers):
    """
    Detects repeating periodic sequences of alternating quantum wells and barriers.
    Returns list of detected MQW regions with start_idx, end_idx, period_count,
    well_layer, barrier_layer, and total_thickness.
    """
    if len(layers) < 3:
        return []

    mqw_regions = []
    i = 0
    while i < len(layers) - 1:
        # Check if layer i is a well or barrier followed by opposite
        t_i = str(layers[i].get("type", "")).lower()
        if "well" in t_i or t_i == "w":
            # Well found, check next for barrier
            well_idx = i
            barrier_idx = i + 1
            if barrier_idx < len(layers):
                t_next = str(layers[barrier_idx].get("type", "")).lower()
                if "barrier" in t_next or t_next == "b":
                    # We have a candidate pair (Well + Barrier)
                    w_ref = layers[well_idx]
                    b_ref = layers[barrier_idx]
                    count = 1
                    curr = barrier_idx + 1

                    while curr + 1 < len(layers):
                        w_cand = layers[curr]
                        b_cand = layers[curr + 1]
                        # Check match
                        w_match = (
                            w_cand.get("material") == w_ref.get("material") and
                            abs(float(w_cand.get("thickness", 0)) - float(w_ref.get("thickness", 0))) < 0.5 and
                            abs(float(w_cand.get("mole", 0)) - float(w_ref.get("mole", 0))) < 0.05
                        )
                        b_match = (
                            b_cand.get("material") == b_ref.get("material") and
                            abs(float(b_cand.get("thickness", 0)) - float(b_ref.get("thickness", 0))) < 0.5 and
                            abs(float(b_cand.get("mole", 0)) - float(b_ref.get("mole", 0))) < 0.05
                        )
                        if w_match and b_match:
                            count += 1
                            curr += 2
                        else:
                            break

                    if count >= 2:
                        end_idx = curr - 1
                        total_th = sum(float(layers[k].get("thickness", 0.0)) for k in range(well_idx, end_idx + 1))
                        mqw_regions.append({
                            "start_idx": well_idx,
                            "end_idx": end_idx,
                            "period_count": count,
                            "well_material": w_ref.get("material"),
                            "barrier_material": b_ref.get("material"),
                            "well_thick": float(w_ref.get("thickness", 0)),
                            "barrier_thick": float(b_ref.get("thickness", 0)),
                            "total_thickness": total_th
                        })
                        i = end_idx + 1
                        continue
        i += 1

    return mqw_regions


# ---------------------------------------------------------------------------
# 5. Dynamic Thickness Scaling (TRUE, SCHEMATIC, AUTO)
# ---------------------------------------------------------------------------

def calculate_display_geometry(layers, scale_mode="AUTO", total_span=100.0, mqw_mode="Detailed"):
    """
    Computes spatial bounding intervals [z_start, z_end] for each layer.
    
    Modes:
      - 'TRUE SCALE': Visual width directly proportional to physical thickness in nm.
      - 'SCHEMATIC SCALE': Non-linear compression allowing thin QWs (3 nm) and thick
        substrates (2000 nm) to both be clearly readable.
      - 'AUTO': Selects SCHEMATIC if max(d)/min(d) > 15.0; otherwise TRUE.
    
    Returns:
      geometry list: list of dicts with physical [z0_phys, z1_phys],
                     display [z0_disp, z1_disp], display thickness, physical thickness,
                     and active scale indicator string.
    """
    if not layers:
        return {"layers": [], "active_scale": "Empty", "total_phys_nm": 0.0}

    phys_thicknesses = [max(0.1, float(l.get("thickness", 10.0))) for l in layers]
    total_phys_nm = sum(phys_thicknesses)
    min_d = min(phys_thicknesses)
    max_d = max(phys_thicknesses)
    ratio = max_d / max(0.1, min_d)

    # Determine active scale
    active_scale = scale_mode.upper()
    if active_scale == "AUTO":
        active_scale = "SCHEMATIC SCALE" if ratio > 15.0 else "TRUE SCALE"

    N = len(layers)
    geom_layers = []

    # Calculate physical coordinates (accumulated along z)
    phys_z = 0.0
    phys_bounds = []
    for d in phys_thicknesses:
        phys_bounds.append((phys_z, phys_z + d))
        phys_z += d

    # Calculate display coordinates
    if active_scale == "TRUE SCALE":
        scale_factor = total_span / max(1.0, total_phys_nm)
        disp_z = 0.0
        for i, (p0, p1) in enumerate(phys_bounds):
            d_phys = phys_thicknesses[i]
            d_disp = d_phys * scale_factor
            geom_layers.append({
                "layer_idx": i,
                "z0_phys": p0,
                "z1_phys": p1,
                "d_phys": d_phys,
                "z0_disp": disp_z,
                "z1_disp": disp_z + d_disp,
                "d_disp": d_disp,
            })
            disp_z += d_disp
    else:
        # SCHEMATIC SCALE: Non-linear power-law compression
        # w_i = d_min_disp + (d_i / d_ref)**0.32 * w_scale
        d_ref = max(1.0, np.median(phys_thicknesses))
        weights = []
        for d in phys_thicknesses:
            w = 1.0 + (d / d_ref) ** 0.35
            weights.append(w)
        total_weight = sum(weights)
        scale_factor = total_span / max(1.0, total_weight)

        disp_z = 0.0
        for i, (p0, p1) in enumerate(phys_bounds):
            d_phys = phys_thicknesses[i]
            d_disp = weights[i] * scale_factor
            geom_layers.append({
                "layer_idx": i,
                "z0_phys": p0,
                "z1_phys": p1,
                "d_phys": d_phys,
                "z0_disp": disp_z,
                "z1_disp": disp_z + d_disp,
                "d_disp": d_disp,
            })
            disp_z += d_disp

    return {
        "layers": geom_layers,
        "active_scale": active_scale,
        "total_phys_nm": total_phys_nm,
        "ratio": ratio
    }


# ---------------------------------------------------------------------------
# 6. Physical & Numerical Consistency Validation
# ---------------------------------------------------------------------------

def parse_and_validate_structure(layers, mat_system="Zincblende", temp_k=300.0, provenance_dict=None):
    """
    Validates user layer configuration against physical and database rules.
    Returns:
      valid: bool
      errors: list of error strings
      warnings: list of warning strings
      parsed_layers: list of enriched layer dictionaries
    """
    errors = []
    warnings = []
    parsed = []

    if not layers or len(layers) == 0:
        return False, ["Structure is empty. Please define at least one layer."], [], []

    total_th = 0.0
    has_well = False
    has_doped = False

    for idx, l in enumerate(layers):
        mat = str(l.get("material", "")).strip()
        try:
            th = float(l.get("thickness", 0.0))
        except (ValueError, TypeError):
            errors.append(f"Layer {idx + 1}: Invalid non-numeric thickness.")
            th = 0.0

        if th <= 0.0:
            errors.append(f"Layer {idx + 1} ({mat}): Thickness must be > 0 (got {th} nm).")

        total_th += th

        # Material check
        known = (
            mat in database.materialproperty or
            mat in database.alloyproperty or
            mat in database.alloyproperty4
        )
        if not known:
            errors.append(f"Layer {idx + 1}: Unknown material '{mat}' not found in database.py.")

        # Mole fraction check
        try:
            x = float(l.get("mole", 0.0))
            if x < 0.0 or x > 1.0:
                errors.append(f"Layer {idx + 1} ({mat}): Mole fraction x must be in [0, 1] (got {x}).")
        except (ValueError, TypeError):
            errors.append(f"Layer {idx + 1} ({mat}): Non-numeric mole fraction x.")
            x = 0.0

        try:
            y = float(l.get("mole_y", 0.0))
            if y < 0.0 or y > 1.0:
                errors.append(f"Layer {idx + 1} ({mat}): Mole fraction y must be in [0, 1] (got {y}).")
        except (ValueError, TypeError):
            errors.append(f"Layer {idx + 1} ({mat}): Non-numeric mole fraction y.")
            y = 0.0

        # Doping check
        try:
            dop = float(l.get("doping", 0.0))
            if dop < 0.0:
                errors.append(f"Layer {idx + 1} ({mat}): Doping cannot be negative (got {dop}).")
            if dop > 1e15:
                has_doped = True
        except (ValueError, TypeError):
            errors.append(f"Layer {idx + 1} ({mat}): Non-numeric doping.")
            dop = 0.0

        dtype = str(l.get("doping_type", "n")).lower()
        if dtype not in ("n", "p", "i", "undoped"):
            warnings.append(f"Layer {idx + 1} ({mat}): Unrecognized doping type '{dtype}', defaulting to 'n'.")

        ltype = str(l.get("type", "barrier")).lower()
        if "well" in ltype:
            has_well = True

        # Extract electronic properties
        props = extract_layer_properties(l, temp_k, mat_system, provenance_dict)
        props["role"] = infer_layer_role(l, idx, len(layers), layers)
        parsed.append(props)

    if total_th > 50000.0:
        warnings.append(f"Total structure thickness is very large ({total_th/1000.0:.1f} µm). Ensure grid resolution is sufficient.")

    valid = len(errors) == 0
    return valid, errors, warnings, parsed


# ---------------------------------------------------------------------------
# 7. Energy Band Profile Extraction
# ---------------------------------------------------------------------------

def calculate_flat_band_profile(layers, mat_system="Zincblende", temp_k=300.0, geom=None):
    """
    Computes flat-band conduction Ec(z) and valence Ev(z) band profiles
    along the growth direction coordinate, matching layer interfaces.
    """
    if not layers:
        return None
    if geom is None:
        geom = calculate_display_geometry(layers)

    geom_layers = geom["layers"]
    z_pts = []
    ec_pts = []
    ev_pts = []

    for i, gl in enumerate(geom_layers):
        l = layers[i]
        props = extract_layer_properties(l, temp_k, mat_system)
        z0 = gl["z0_disp"]
        z1 = gl["z1_disp"]
        ec = props["Ec"]
        ev = props["Ev"]

        # Insert points to form crisp rectangular step-heterojunctions
        z_pts.extend([z0, z1])
        ec_pts.extend([ec, ec])
        ev_pts.extend([ev, ev])

    return {
        "z": np.array(z_pts),
        "ec": np.array(ec_pts),
        "ev": np.array(ev_pts)
    }


# ---------------------------------------------------------------------------
# 8. Publication-Quality Rendering Engine
# ---------------------------------------------------------------------------

# Color Palettes
THEMES = {
    "Dark (GUI)": {
        "fig_bg": "#1E1E24",
        "ax_bg": "#252530",
        "text": "#EAECEE",
        "subtext": "#BDC3C7",
        "border": "#3A3D4D",
        "grid": "#333544",
        "contact": "#474B5E",
        "contact_edge": "#7F8C8D",
        "p_fill": "#78283B",
        "p_edge": "#E74C3C",
        "n_fill": "#1B4965",
        "n_edge": "#3498DB",
        "i_fill": "#1F4E43",
        "i_edge": "#2ECC71",
        "well_fill": "#C0392B",
        "well_edge": "#F39C12",
        "ec_line": "#E74C3C",
        "ev_line": "#3498DB",
        "ef_line": "#F1C40F",
        "badge_bg": "#14151B",
        "badge_border": "#474B5E",
        "highlight": "#00E5FF",
        "annotation_photons": "#F1C40F",
        "annotation_carrier_e": "#3498DB",
        "annotation_carrier_h": "#E74C3C",
        "annotation_field": "#9B59B6"
    },
    "Publication (Light)": {
        "fig_bg": "#FFFFFF",
        "ax_bg": "#FAFAFA",
        "text": "#1A1A1A",
        "subtext": "#4A4A4A",
        "border": "#2C3E50",
        "grid": "#E0E0E0",
        "contact": "#D5D8DC",
        "contact_edge": "#5D6D7E",
        "p_fill": "#FADBD8",
        "p_edge": "#900C3F",
        "n_fill": "#D4E6F1",
        "n_edge": "#1B4F72",
        "i_fill": "#D5F5E3",
        "i_edge": "#196F3D",
        "well_fill": "#FCF3CF",
        "well_edge": "#B7950B",
        "ec_line": "#C0392B",
        "ev_line": "#2980B9",
        "ef_line": "#D4AC0D",
        "badge_bg": "#F4F6F6",
        "badge_border": "#BDC3C7",
        "highlight": "#D35400",
        "annotation_photons": "#B7950B",
        "annotation_carrier_e": "#2980B9",
        "annotation_carrier_h": "#C0392B",
        "annotation_field": "#8E44AD"
    }
}


def render_device_diagram(
    fig,
    layers,
    config=None,
    view_mode="Structure + Bands",
    scale_mode="AUTO",
    mqw_mode="Detailed",
    show_annotations=True,
    theme="Dark (GUI)",
    selected_layer_idx=None,
    device_type="Generic Diode / LED",
    temp_k=300.0,
    provenance_dict=None,
    qw_result=None,
    show_qw_states=False,
    qw_display_mode="Wavefunction ψ",
    show_qw_transitions=False,
    show_electron_states=True,
    show_hole_states=True
):
    """
    Master rendering routine. Populates `fig` with publication-standard visualization.
    
    Modes:
      - 'Structure + Bands': Dual aligned subplots sharing growth z-axis with tie-lines.
      - 'Layer Structure': Dedicated full cross-sectional heterostructure stack.
      - 'Band Diagram': Flat-band conduction and valence band energy profile.
    """
    palette = THEMES.get(theme, THEMES["Dark (GUI)"])
    fig.patch.set_facecolor(palette["fig_bg"])

    if not layers:
        fig.clear()
        fig._struct_view_mode = None
        ax = fig.add_subplot(111)
        ax.set_facecolor(palette["ax_bg"])
        ax.text(0.5, 0.5, "No Semiconductor Layers Defined\nClick '+ Add Layer' or Load an Example",
                ha="center", va="center", color=palette["subtext"], fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])
        return None

    # Calculate geometry
    geom = calculate_display_geometry(layers, scale_mode=scale_mode, total_span=100.0, mqw_mode=mqw_mode)
    band_data = calculate_flat_band_profile(layers, temp_k=temp_k, geom=geom)

    expected_axes_count = 2 if view_mode == "Structure + Bands" else 1
    can_reuse_axes = (
        getattr(fig, "_struct_view_mode", None) == view_mode and
        len(fig.axes) == expected_axes_count
    )

    if can_reuse_axes:
        for ax in fig.axes:
            ax.cla()
    else:
        fig.clear()
        fig._struct_view_mode = view_mode

    # 1. View: Structure + Band Diagram (Canonical dual layout)
    if view_mode == "Structure + Bands":
        if can_reuse_axes:
            ax_top, ax_bot = fig.axes[0], fig.axes[1]
        else:
            gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1.0], hspace=0.25)
            ax_top = fig.add_subplot(gs[0])
            ax_bot = fig.add_subplot(gs[1], sharex=ax_top)

        _draw_layer_stack_horizontal(
            ax_top, layers, geom, palette, selected_layer_idx,
            show_annotations=show_annotations, device_type=device_type
        )
        _draw_band_profile_horizontal(
            ax_bot, band_data, geom, palette, selected_layer_idx,
            qw_result=qw_result, show_qw_states=show_qw_states,
            qw_display_mode=qw_display_mode, show_qw_transitions=show_qw_transitions,
            show_electron_states=show_electron_states, show_hole_states=show_hole_states,
            layers=layers, temp_k=temp_k
        )
        _draw_tie_lines(ax_top, ax_bot, geom, palette)

        # Scale and provenance badge on top panel
        _draw_scale_badge(ax_top, geom, palette)

    # 2. View: Layer Structure Only
    elif view_mode == "Layer Structure":
        ax = fig.axes[0] if can_reuse_axes else fig.add_subplot(111)
        _draw_layer_stack_detailed(
            ax, layers, geom, palette, selected_layer_idx,
            show_annotations=show_annotations, device_type=device_type
        )
        _draw_scale_badge(ax, geom, palette)

    # 3. View: Band Diagram Only
    elif view_mode == "Band Diagram":
        ax = fig.axes[0] if can_reuse_axes else fig.add_subplot(111)
        _draw_band_profile_horizontal(
            ax, band_data, geom, palette, selected_layer_idx, is_standalone=True,
            qw_result=qw_result, show_qw_states=show_qw_states,
            qw_display_mode=qw_display_mode, show_qw_transitions=show_qw_transitions,
            show_electron_states=show_electron_states, show_hole_states=show_hole_states,
            layers=layers, temp_k=temp_k
        )
        _draw_scale_badge(ax, geom, palette)

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        try:
            fig.tight_layout()
        except Exception:
            fig.subplots_adjust(left=0.12, right=0.96, top=0.94, bottom=0.12, hspace=0.32)
    return geom


def _draw_layer_stack_horizontal(ax, layers, geom, palette, selected_layer_idx=None, show_annotations=True, device_type="Generic Diode / LED"):
    """Draws horizontal cross-sectional layer stack aligned along growth coordinate z."""
    ax.set_facecolor(palette["ax_bg"])
    y_min, y_max = 0.0, 1.0

    geom_layers = geom["layers"]
    total_z = 100.0

    # Draw Contact Anode/Left
    contact_w = 4.0
    rect_left = patches.Rectangle((-contact_w, y_min), contact_w, y_max - y_min,
                                  facecolor=palette["contact"], edgecolor=palette["contact_edge"],
                                  hatch="///", linewidth=1.2, zorder=2)
    ax.add_patch(rect_left)
    ax.text(-contact_w / 2.0, 0.5, "Top Contact\n(Anode)", color=palette["text"],
            fontsize=8, ha="center", va="center", rotation=90, weight="bold")

    # Draw Semiconductor Layers
    for i, gl in enumerate(geom_layers):
        l = layers[i]
        z0 = gl["z0_disp"]
        z1 = gl["z1_disp"]
        w = gl["d_disp"]
        d_phys = gl["d_phys"]

        dtype = str(l.get("doping_type", "n")).lower()
        ltype = str(l.get("type", "barrier")).lower()
        mat = str(l.get("material", ""))

        # Palette selection
        if "well" in ltype:
            fc = palette["well_fill"]
            ec = palette["well_edge"]
            lw = 1.8
        elif dtype == "p":
            fc = palette["p_fill"]
            ec = palette["p_edge"]
            lw = 1.0
        elif dtype == "n":
            fc = palette["n_fill"]
            ec = palette["n_edge"]
            lw = 1.0
        else:
            fc = palette["i_fill"]
            ec = palette["i_edge"]
            lw = 1.0

        # Selected layer highlight
        is_sel = (selected_layer_idx is not None and selected_layer_idx == i)
        if is_sel:
            ec = palette["highlight"]
            lw = 2.5

        rect = patches.Rectangle((z0, y_min), w, y_max - y_min,
                                 facecolor=fc, edgecolor=ec, linewidth=lw,
                                 zorder=3, alpha=0.92)
        ax.add_patch(rect)

        # Label inside layer if width permits, or staggered
        formula = format_chemical_formula(mat, l.get("mole"), l.get("mole_y"), use_mathtext=True)
        th_str = format_thickness_str(d_phys)

        if w >= 12.0:
            ax.text(z0 + w / 2.0, 0.65, formula, color=palette["text"],
                    fontsize=8.5, ha="center", va="center", weight="bold")
            ax.text(z0 + w / 2.0, 0.35, th_str, color=palette["subtext"],
                    fontsize=7.5, ha="center", va="center")
        elif w >= 6.0:
            ax.text(z0 + w / 2.0, 0.5, formula, color=palette["text"],
                    fontsize=7.5, ha="center", va="center", rotation=90)
        else:
            # Thin layer (QW) indicator line above
            ax.text(z0 + w / 2.0, 1.08, f"QW\n{th_str}", color=palette["well_edge"],
                    fontsize=7.0, ha="center", va="bottom", weight="bold")

    # Draw Contact Cathode/Right
    rect_right = patches.Rectangle((total_z, y_min), contact_w, y_max - y_min,
                                   facecolor=palette["contact"], edgecolor=palette["contact_edge"],
                                   hatch="///", linewidth=1.2, zorder=2)
    ax.add_patch(rect_right)
    ax.text(total_z + contact_w / 2.0, 0.5, "Back Contact\n(Cathode)", color=palette["text"],
            fontsize=8, ha="center", va="center", rotation=270, weight="bold")

    # Annotations
    if show_annotations:
        _draw_annotations_horizontal(ax, layers, geom, palette, device_type)

    ax.set_xlim(-contact_w - 2.0, total_z + contact_w + 2.0)
    ax.set_ylim(-0.15, 1.35)
    ax.set_yticks([])
    ax.set_ylabel("Heterostructure Stack", color=palette["text"], fontsize=9, weight="bold")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_color(palette["border"])


def _draw_band_profile_horizontal(
    ax, band_data, geom, palette, selected_layer_idx=None, is_standalone=False,
    qw_result=None, show_qw_states=False, qw_display_mode="Wavefunction ψ",
    show_qw_transitions=False, show_electron_states=True, show_hole_states=True,
    layers=None, temp_k=300.0
):
    """Draws conduction and valence band profiles horizontally aligned with layer stack."""
    ax.set_facecolor(palette["ax_bg"])
    if not band_data:
        return

    z = band_data["z"]
    ec = band_data["ec"]
    ev = band_data["ev"]

    # Band curves
    ax.plot(z, ec, color=palette["ec_line"], linewidth=2.0, label="$E_c$ (Conduction Band)", zorder=4)
    ax.plot(z, ev, color=palette["ev_line"], linewidth=2.0, label="$E_v$ (Valence Band)", zorder=4)

    # Fermi reference line (E = 0)
    ax.axhline(0.0, color=palette["ef_line"], linestyle="--", linewidth=1.2, alpha=0.7, label="$E_F$ Reference", zorder=3)

    # Fill bandgap region lightly
    ax.fill_between(z, ev, ec, color=palette["text"], alpha=0.03, zorder=2)

    # Highlight selected layer
    if selected_layer_idx is not None and selected_layer_idx < len(geom["layers"]):
        gl = geom["layers"][selected_layer_idx]
        ax.axvspan(gl["z0_disp"], gl["z1_disp"], color=palette["highlight"], alpha=0.18, zorder=1)

    # Quantum Well Confined States & Optical Transitions Overlay
    if show_qw_states and layers:
        _draw_qw_confined_states(
            ax, band_data, geom, palette,
            qw_result=qw_result, show_qw_states=show_qw_states,
            qw_display_mode=qw_display_mode, show_qw_transitions=show_qw_transitions,
            show_electron_states=show_electron_states, show_hole_states=show_hole_states,
            layers=layers, temp_k=temp_k
        )

    ax.set_xlabel("Device Growth Axis Coordinate $z$ (Schematic / Proportional Display)", color=palette["text"], fontsize=9)
    ax.set_ylabel("Energy (eV)", color=palette["text"], fontsize=9, weight="bold")
    ax.tick_params(colors=palette["subtext"], labelsize=8)
    ax.grid(True, color=palette["grid"], linestyle=":", alpha=0.6)
    ax.legend(loc="lower left", fontsize=8, facecolor=palette["ax_bg"], edgecolor=palette["border"], labelcolor=palette["text"])

    for spine in ax.spines.values():
        spine.set_color(palette["border"])


def _draw_qw_confined_states(
    ax, band_data, geom, palette,
    qw_result=None, show_qw_states=False, qw_display_mode="Wavefunction ψ",
    show_qw_transitions=False, show_electron_states=True, show_hole_states=True,
    layers=None, temp_k=300.0
):
    """
    Renders quantized electron and hole eigenenergy levels, envelope wavefunctions
    or probability densities, and optical transition arrows within quantum well regions.
    """
    if not show_qw_states or not layers:
        return

    # Check if any layer is designated as a quantum well
    has_well = any("well" in str(l.get("type", "")).lower() or str(l.get("type", "")) == "w" for l in layers)
    if not has_well:
        return

    # Auto-solve QW states if not pre-computed
    if qw_result is None:
        try:
            from aeslibs.quantum_well import solve_quantum_well
            qw_result = solve_quantum_well(
                band_profile=band_data,
                layers=layers,
                temperature_k=temp_k,
                num_electron_states=3,
                num_hole_states=3
            )
        except Exception:
            return

    if qw_result is None or (len(qw_result.electron_energies) == 0 and len(qw_result.hole_energies) == 0):
        return

    z_display = band_data["z"]
    z_min_disp = 100.0
    z_max_disp = 0.0
    well_found = False

    for i, gl in enumerate(geom["layers"]):
        l = layers[i]
        ltype = str(l.get("type", "")).lower()
        if "well" in ltype or ltype == "w":
            well_found = True
            z_min_disp = min(z_min_disp, gl["z0_disp"])
            z_max_disp = max(z_max_disp, gl["z1_disp"])

    if not well_found:
        z_min_disp = 40.0
        z_max_disp = 60.0

    # Slight padding into barrier for wavefunction exponential tails
    span = z_max_disp - z_min_disp
    pad = max(2.0, span * 0.25)
    z_start_disp = max(0.0, z_min_disp - pad)
    z_end_disp = min(100.0, z_max_disp + pad)

    mask = (z_display >= z_start_disp) & (z_display <= z_end_disp)
    if not np.any(mask):
        mask = np.ones(len(z_display), dtype=bool)

    z_sub = z_display[mask]

    c_e = palette.get("qw_e_color", "#F39C12")
    c_h = palette.get("qw_h_color", "#3498DB")
    c_trans = palette.get("annotation_photons", "#F1C40F")

    # 1. Electron Confined States
    if show_electron_states and len(qw_result.electron_energies) > 0:
        lbl_e = np.array(qw_result.electron_energies, dtype=float).copy()
        for k in range(1, len(lbl_e)):
            if abs(lbl_e[k] - lbl_e[k - 1]) < 0.035:
                lbl_e[k] = lbl_e[k - 1] + 0.035

        for idx, e_val in enumerate(qw_result.electron_energies):
            # Eigenenergy level line
            ax.hlines(e_val, z_min_disp - 0.5, z_max_disp + 0.5,
                      colors=c_e, linestyles="-", linewidth=1.8, zorder=5)
            ax.text(z_max_disp + 1.2, lbl_e[idx], f"$e_{idx+1}$ ({e_val:.3f} eV)",
                    color=c_e, fontsize=7.5, va="center", weight="bold", zorder=6)

            # Envelope wavefunction / probability density
            if idx < len(qw_result.electron_wavefunctions):
                psi_raw = qw_result.electron_wavefunctions[idx]
                if len(psi_raw) == len(z_display):
                    psi_sub = psi_raw[mask]
                else:
                    psi_sub = np.interp(z_sub, qw_result.z, psi_raw)

                max_amp = np.max(np.abs(psi_sub)) + 1e-12
                scale_height = 0.06 / max_amp

                if "prob" in qw_display_mode.lower() or "2" in qw_display_mode:
                    prob_sub = psi_sub ** 2
                    scale_p = 0.06 / (np.max(prob_sub) + 1e-12)
                    curve = e_val + prob_sub * scale_p
                    ax.plot(z_sub, curve, color=c_e, linewidth=1.4, zorder=5)
                    ax.fill_between(z_sub, e_val, curve, color=c_e, alpha=0.25, zorder=4)
                else:
                    curve = e_val + psi_sub * scale_height
                    ax.plot(z_sub, curve, color=c_e, linewidth=1.3, linestyle="--", zorder=5)
                    ax.fill_between(z_sub, e_val, curve, color=c_e, alpha=0.15, zorder=4)

    # 2. Hole Confined States
    if show_hole_states and len(qw_result.hole_energies) > 0:
        lbl_h = np.array(qw_result.hole_energies, dtype=float).copy()
        for k in range(1, len(lbl_h)):
            if abs(lbl_h[k] - lbl_h[k - 1]) < 0.035:
                lbl_h[k] = lbl_h[k - 1] - 0.035

        for idx, h_val in enumerate(qw_result.hole_energies):
            htype = qw_result.hole_types[idx] if idx < len(qw_result.hole_types) else f"h{idx+1}"
            ax.hlines(h_val, z_min_disp - 0.5, z_max_disp + 0.5,
                      colors=c_h, linestyles="-", linewidth=1.8, zorder=5)
            ax.text(z_max_disp + 1.2, lbl_h[idx], f"${htype}$ ({h_val:.3f} eV)",
                    color=c_h, fontsize=7.5, va="center", weight="bold", zorder=6)

            if idx < len(qw_result.hole_wavefunctions):
                psi_raw = qw_result.hole_wavefunctions[idx]
                if len(psi_raw) == len(z_display):
                    psi_sub = psi_raw[mask]
                else:
                    psi_sub = np.interp(z_sub, qw_result.z, psi_raw)

                max_amp = np.max(np.abs(psi_sub)) + 1e-12
                scale_height = 0.06 / max_amp

                if "prob" in qw_display_mode.lower() or "2" in qw_display_mode:
                    prob_sub = psi_sub ** 2
                    scale_p = 0.06 / (np.max(prob_sub) + 1e-12)
                    curve = h_val - prob_sub * scale_p
                    ax.plot(z_sub, curve, color=c_h, linewidth=1.4, zorder=5)
                    ax.fill_between(z_sub, curve, h_val, color=c_h, alpha=0.25, zorder=4)
                else:
                    curve = h_val - psi_sub * scale_height
                    ax.plot(z_sub, curve, color=c_h, linewidth=1.3, linestyle="--", zorder=5)
                    ax.fill_between(z_sub, curve, h_val, color=c_h, alpha=0.15, zorder=4)

    # 3. Optical Transition Arrow & Metric Badge
    if show_qw_transitions and qw_result.dominant_transitions:
        top_t = qw_result.dominant_transitions[0]
        e_idx = top_t["e_index"] - 1
        h_idx = top_t["h_index"] - 1

        if e_idx < len(qw_result.electron_energies) and h_idx < len(qw_result.hole_energies):
            e_top = qw_result.electron_energies[e_idx]
            h_bot = qw_result.hole_energies[h_idx]
            z_center = (z_min_disp + z_max_disp) / 2.0

            ax.annotate('', xy=(z_center, e_top), xytext=(z_center, h_bot),
                        arrowprops=dict(arrowstyle="<->", color=c_trans, lw=2.2,
                                       mutation_scale=14, shrinkA=2, shrinkB=2),
                        zorder=7)

            badge_text = (
                f"hν ({top_t['name']})\n"
                f"λ = {top_t['wavelength_nm']:.1f} nm\n"
                f"ΔE = {top_t['energy_ev']:.3f} eV\n"
                f"Γ = {top_t['overlap']:.2f}"
            )
            ax.text(z_center - 1.5, (e_top + h_bot) / 2.0, badge_text,
                    color=c_trans, fontsize=8.0, ha="right", va="center", weight="bold",
                    bbox=dict(boxstyle="round,pad=0.35", facecolor=palette["badge_bg"],
                              edgecolor=c_trans, linewidth=1.2, alpha=0.92),
                    zorder=8)


def _draw_tie_lines(ax_top, ax_bot, geom, palette):
    """Connects top heterostructure layer interfaces to bottom band discontinuities."""
    for gl in geom["layers"]:
        z0 = gl["z0_disp"]
        z1 = gl["z1_disp"]
        ax_top.axvline(z0, color=palette["border"], linestyle=":", linewidth=0.8, alpha=0.4, zorder=1)
        ax_top.axvline(z1, color=palette["border"], linestyle=":", linewidth=0.8, alpha=0.4, zorder=1)
        ax_bot.axvline(z0, color=palette["border"], linestyle=":", linewidth=0.8, alpha=0.4, zorder=1)
        ax_bot.axvline(z1, color=palette["border"], linestyle=":", linewidth=0.8, alpha=0.4, zorder=1)


def _draw_layer_stack_detailed(ax, layers, geom, palette, selected_layer_idx=None, show_annotations=True, device_type="Generic Diode / LED"):
    """
    Renders high-detail vertical cross-sectional layer stack (Substrate at bottom,
    growth upwards, cap at top), standard in semiconductor textbooks.
    """
    ax.set_facecolor(palette["ax_bg"])
    geom_layers = geom["layers"]
    N = len(geom_layers)
    x_min, x_max = 0.2, 0.8
    w_box = x_max - x_min

    # Bottom Contact
    ax.add_patch(patches.Rectangle((x_min, -0.05), w_box, 0.05,
                                   facecolor=palette["contact"], edgecolor=palette["contact_edge"],
                                   hatch="///", linewidth=1.2))
    ax.text(0.5, -0.025, "Back Contact / Cathode Metallization", color=palette["text"],
            fontsize=8, ha="center", va="center", weight="bold")

    # Invert so growth goes bottom-to-top (Substrate at base)
    total_span = 1.0
    for i, gl in enumerate(geom_layers):
        l = layers[i]
        d_phys = gl["d_phys"]
        y0 = gl["z0_disp"] / 100.0
        h = gl["d_disp"] / 100.0

        dtype = str(l.get("doping_type", "n")).lower()
        ltype = str(l.get("type", "barrier")).lower()
        mat = str(l.get("material", ""))

        if "well" in ltype:
            fc = palette["well_fill"]
            ec = palette["well_edge"]
            lw = 1.8
        elif dtype == "p":
            fc = palette["p_fill"]
            ec = palette["p_edge"]
            lw = 1.0
        elif dtype == "n":
            fc = palette["n_fill"]
            ec = palette["n_edge"]
            lw = 1.0
        else:
            fc = palette["i_fill"]
            ec = palette["i_edge"]
            lw = 1.0

        if selected_layer_idx is not None and selected_layer_idx == i:
            ec = palette["highlight"]
            lw = 2.5

        rect = patches.Rectangle((x_min, y0), w_box, h,
                                 facecolor=fc, edgecolor=ec, linewidth=lw,
                                 alpha=0.92, zorder=3)
        ax.add_patch(rect)

        # Labels
        formula = format_chemical_formula(mat, l.get("mole"), l.get("mole_y"), use_mathtext=True)
        th_str = format_thickness_str(d_phys)
        dop_str = format_doping_str(l.get("doping"), dtype)
        role = infer_layer_role(l, i, N, layers)

        ax.text(0.5, y0 + h / 2.0, f"{formula}  ({th_str})\n{role}  |  {dop_str}",
                color=palette["text"], fontsize=8, ha="center", va="center", weight="bold", zorder=4)

    # Top Contact
    ax.add_patch(patches.Rectangle((x_min, 1.0), w_box, 0.05,
                                   facecolor=palette["contact"], edgecolor=palette["contact_edge"],
                                   hatch="///", linewidth=1.2))
    ax.text(0.5, 1.025, "Top Contact / Anode Metallization", color=palette["text"],
            fontsize=8, ha="center", va="center", weight="bold")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.08, 1.12)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color(palette["border"])


def _draw_annotations_horizontal(ax, layers, geom, palette, device_type="Generic Diode / LED"):
    """Draws device-specific physics overlays (photons, carrier injection, optical cavity)."""
    dev = str(device_type).lower()

    # Find active quantum well region or central junction
    well_indices = [i for i, l in enumerate(layers) if "well" in str(l.get("type", "")).lower()]
    if well_indices:
        act_z0 = geom["layers"][well_indices[0]]["z0_disp"]
        act_z1 = geom["layers"][well_indices[-1]]["z1_disp"]
        act_center = (act_z0 + act_z1) / 2.0
    else:
        act_center = 50.0

    # 1. LED Annotations
    if "led" in dev:
        # Radiative recombination marker
        ax.plot([act_center], [0.5], marker="*", markersize=14, color=palette["annotation_photons"], zorder=5)
        # Emission photons (wavy arrows pointing upwards)
        ax.annotate("", xy=(act_center, 1.25), xytext=(act_center, 0.8),
                    arrowprops=dict(arrowstyle="->", color=palette["annotation_photons"], lw=2.0, connectionstyle="arc3,rad=-0.2"), zorder=5)
        ax.text(act_center, 1.28, "$h\\nu$ (Light Emission)", color=palette["annotation_photons"],
                fontsize=8.5, ha="center", va="bottom", weight="bold")

        # Carrier injection
        ax.annotate("e⁻ injection →", xy=(act_center - 15.0, -0.08), xytext=(act_center - 35.0, -0.08),
                    color=palette["annotation_carrier_e"], fontsize=8, weight="bold",
                    arrowprops=dict(arrowstyle="->", color=palette["annotation_carrier_e"], lw=1.5))
        ax.annotate("← h⁺ injection", xy=(act_center + 15.0, -0.08), xytext=(act_center + 35.0, -0.08),
                    color=palette["annotation_carrier_h"], fontsize=8, weight="bold",
                    arrowprops=dict(arrowstyle="->", color=palette["annotation_carrier_h"], lw=1.5))

    # 2. Laser Diode Annotations
    elif "laser" in dev:
        # Cavity axis arrow
        ax.annotate("", xy=(95.0, -0.08), xytext=(5.0, -0.08),
                    arrowprops=dict(arrowstyle="<->", color=palette["annotation_field"], lw=2.0), zorder=5)
        ax.text(50.0, -0.06, "Optical Cavity Longitudinal Axis $\\longleftrightarrow$ ($R_1, R_2$ Facets)",
                color=palette["annotation_field"], fontsize=8, ha="center", va="top", weight="bold")

        # Optical confinement bracket
        ax.annotate("Optical Waveguide / Active Region (Γ)", xy=(act_center, 1.18),
                    color=palette["well_edge"], fontsize=8.5, ha="center", va="bottom", weight="bold")

    # 3. Solar Cell / Photodetector Annotations
    elif "solar" in dev or "detector" in dev or "diode" in dev:
        # Incident sunlight photons
        ax.annotate("", xy=(15.0, 0.95), xytext=(5.0, 1.30),
                    arrowprops=dict(arrowstyle="->", color=palette["annotation_photons"], lw=1.8), zorder=5)
        ax.annotate("", xy=(25.0, 0.95), xytext=(15.0, 1.30),
                    arrowprops=dict(arrowstyle="->", color=palette["annotation_photons"], lw=1.8), zorder=5)
        ax.text(10.0, 1.32, "Sunlight Flux\n(AM1.5G $h\\nu$)", color=palette["annotation_photons"],
                fontsize=8, ha="center", va="bottom", weight="bold")

        # Built-in electric field
        ax.annotate("", xy=(act_center + 10.0, 0.15), xytext=(act_center - 10.0, 0.15),
                    arrowprops=dict(arrowstyle="->", color=palette["annotation_field"], lw=1.5), zorder=5)
        ax.text(act_center, 0.22, "$\\vec{E}_{\\mathrm{bi}}$ (Built-in Field)", color=palette["annotation_field"],
                fontsize=7.5, ha="center", va="bottom", weight="bold")


def _draw_scale_badge(ax, geom, palette):
    """Draws a clean scale mode indicator badge."""
    active_scale = geom.get("active_scale", "AUTO")
    total_phys = geom.get("total_phys_nm", 0.0)
    badge_str = f"Scale: {active_scale}  |  Total: {format_thickness_str(total_phys)}"

    bbox_props = dict(boxstyle="round,pad=0.35", facecolor=palette["badge_bg"],
                      edgecolor=palette["badge_border"], alpha=0.88, linewidth=0.8)
    ax.text(0.985, 0.94, badge_str, transform=ax.transAxes,
            color=palette["subtext"], fontsize=7.5, ha="right", va="top",
            bbox=bbox_props, zorder=6)


# ---------------------------------------------------------------------------
# 9. Interactive Hit Detection
# ---------------------------------------------------------------------------

def find_clicked_layer(event, geom, view_mode="Structure + Bands"):
    """
    Determines which layer index was clicked based on canvas coordinates.
    Returns: int (layer_idx) or None
    """
    if event.xdata is None or event.ydata is None or not geom:
        return None

    x_click = event.xdata
    for gl in geom.get("layers", []):
        z0 = gl["z0_disp"]
        z1 = gl["z1_disp"]
        if z0 <= x_click <= z1:
            return gl["layer_idx"]

    return None


# ---------------------------------------------------------------------------
# 10. Scientific Publication Export
# ---------------------------------------------------------------------------

def export_publication_figure(fig, filepath, dpi=300, format=None):
    """
    Exports publication-grade figure to PNG, SVG, or PDF.
    Ensures transparent or white publication background and crisp vector strokes.
    """
    if not filepath:
        raise ValueError("Invalid target export path.")

    ext = os.path.splitext(filepath)[1].lower().replace(".", "")
    fmt = format if format else (ext if ext else "png")

    # High quality export kwargs
    kwargs = {
        "dpi": dpi,
        "bbox_inches": "tight",
        "pad_inches": 0.05
    }
    if fmt in ("svg", "pdf", "eps"):
        kwargs.pop("dpi", None)

    fig.savefig(filepath, format=fmt, **kwargs)
    return filepath
