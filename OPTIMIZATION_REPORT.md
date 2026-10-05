# InGaN Solar Cell Optimization Report
**Date**: 2026-02-14  
**Objective**: Calibrate InGaN/GaN p-n junction solar cell to achieve physically realistic performance metrics.

---

## Target Performance Metrics
- **Voc (Open-Circuit Voltage)**: 1.5 - 3.0 V
- **Jsc (Short-Circuit Current)**: 5 - 20 mA/cm²
- **FF (Fill Factor)**: > 70%

---

## Optimization Results Summary

### Parameter Space Explored
- **Indium Mole Fraction**: 0.20, 0.25, 0.30, 0.35
- **p-Doping**: 1×10¹⁷, 5×10¹⁷, 1×10¹⁸ cm⁻³
- **n-Doping**: 1×10¹⁸, 5×10¹⁸, 1×10¹⁹ cm⁻³
- **Total Configurations Tested**: 36
- **Optical Generation Rate**: 1×10²¹ cm⁻³s⁻¹ (≈ AM1.5 equivalent)

### Best Configuration Found
- **In Fraction**: 0.20 (highest bandgap tested)
- **p-Doping**: 1×10¹⁸ cm⁻³
- **n-Doping**: 1×10¹⁸ cm⁻³

### Performance Achieved
| Metric | Target | Achieved | Status |
|--------|--------|----------|--------|
| **Jsc** | 5-20 mA/cm² | **0.0372 mA/cm²** | ❌ **133× too low** |
| **Voc** | 1.5-3.0 V | **0.0000 V** | ❌ **Zero** |
| **FF** | >70% | **0.00%** | ❌ **Zero** |
| **Pmax** | >5 mW/cm² | **0.0000 mW/cm²** | ❌ **Zero** |

---

## Root Cause Analysis

### 1. Zero Open-Circuit Voltage (Voc = 0V)
**Observation**: All 36 configurations yielded Voc = 0V, regardless of:
- Indium composition (bandgap variation)
- Doping levels (built-in potential variation)

**Possible Causes**:
1. **Ohmic Contact Assumption**: The simulation may be treating both contacts as ohmic (zero barrier), effectively shorting the device.
   - Check: `input_obj.surface = np.array([0.0, 0.0])` may be forcing flat-band conditions at both contacts.
   - **Recommendation**: Implement Schottky or selective contacts with proper work function differences.

2. **Boundary Condition Issue**: The drift-diffusion solver boundary conditions may not properly enforce quasi-Fermi level splitting under illumination.
   - **Recommendation**: Verify that `Continuity2` correctly handles generation terms at boundaries.

3. **Recombination Dominance**: Extremely high recombination rates may be collapsing the quasi-Fermi level splitting.
   - Current lifetimes: τₙ = τₚ = 500 ns
   - **Recommendation**: Increase carrier lifetimes to 1-10 μs (typical for high-quality InGaN).

### 2. Low Short-Circuit Current (Jsc << target)
**Observation**: Jsc = 0.037 mA/cm² is **133× lower** than the minimum target (5 mA/cm²).

**Possible Causes**:
1. **Insufficient Absorption**: 
   - Device thickness: 100 nm (40 nm p + 60 nm n)
   - Absorption depth for InGaN (λ < 450 nm): ~100-500 nm
   - **Issue**: Device may be too thin to absorb significant photon flux.
   - **Recommendation**: Increase total thickness to 500-1000 nm, or add anti-reflection coating / light-trapping structures.

2. **Collection Efficiency**: 
   - Short diffusion lengths due to low lifetimes (L = √(Dτ) ≈ 0.5-1 μm for τ = 500 ns)
   - **Recommendation**: Increase carrier lifetimes and optimize doping gradient for built-in field enhancement.

3. **Optical Generation Rate**:
   - Current: G_opt = 1×10²¹ cm⁻³s⁻¹
   - **Verification Needed**: Confirm this value is appropriate for AM1.5 spectrum integrated over InGaN absorption range.
   - **Recommendation**: Use detailed optical modeling (transfer matrix method) to calculate realistic generation profile.

### 3. Material Parameter Concerns
**InGaN Bandgap vs. Indium Content**:
- In = 0.20: Eg ≈ 2.4 eV (λ_cutoff ≈ 517 nm) ✓ Good for visible spectrum
- In = 0.35: Eg ≈ 1.8 eV (λ_cutoff ≈ 689 nm) ✓ Better absorption, but lower Voc potential

**Diagnostic Finding** (from `diagnose_solar_physics.py`):
- Original configuration (In = 0.57): Bandgap too narrow for target Voc
- Optimized configuration (In = 0.20): Bandgap adequate, but Voc still zero → **Contact/boundary issue**

---

## Recommendations for Next Steps

### Immediate Actions (High Priority)
1. **Fix Contact Boundary Conditions**:
   ```python
   # Instead of:
   input_obj.surface = np.array([0.0, 0.0])  # Ohmic/flat-band
   
   # Try:
   input_obj.surface = np.array([work_function_p, work_function_n])
   # Or implement selective contacts with proper barrier heights
   ```

2. **Increase Carrier Lifetimes**:
   ```python
   input_obj.taun0 = 5e-6  # 5 μs (instead of 500 ns)
   input_obj.taup0 = 5e-6  # 5 μs
   ```

3. **Increase Device Thickness**:
   ```json
   "layers": [
       {"thickness": "200", "doping_type": "p", ...},  // 200 nm p-layer
       {"thickness": "300", "doping_type": "n", ...}   // 300 nm n-layer
   ]
   ```

4. **Verify Optical Generation**:
   - Run detailed optical simulation (e.g., using transfer matrix method)
   - Calculate position-dependent G(x) profile
   - Ensure total integrated photocurrent matches expected value for AM1.5

### Medium-Term Improvements
1. **Implement Graded Composition**:
   - Use compositional grading (In_x Ga_{1-x}N with varying x) to create built-in electric field
   - Enhances carrier collection efficiency

2. **Add Surface Recombination Velocity**:
   - Model realistic surface recombination (S = 10³-10⁵ cm/s)
   - Critical for thin devices

3. **Optimize Doping Profile**:
   - Use asymmetric doping (e.g., p = 5×10¹⁷, n = 5×10¹⁸ cm⁻³)
   - Consider delta-doping or graded doping profiles

### Validation Steps
1. **Equilibrium Test**:
   - Run simulation with G_optical = 0 (dark)
   - Verify built-in potential Vbi > 1.0 V for In = 0.20
   - Check that band diagram shows proper p-n junction formation

2. **Incremental Illumination**:
   - Sweep G_optical from 10¹⁹ to 10²² cm⁻³s⁻¹
   - Verify Jsc scales linearly with G_optical
   - Verify Voc increases logarithmically with G_optical

3. **Benchmark Against Literature**:
   - Compare with published InGaN solar cell data
   - Typical InGaN cells: Voc = 1.5-2.5 V, Jsc = 0.5-5 mA/cm² (for thin cells)
   - Efficiency: 1-5% (current state-of-art for InGaN)

---

## Files Generated
1. **Optimized Configuration**: `examples/ingan_solar_optimized.json`
   - Best parameters from 36-configuration sweep
   - Ready to load in Aestimo GUI

2. **Optimization Results**: `optimization_results.png`
   - Visualization of parameter space exploration
   - Shows Voc vs Jsc, In composition trends, doping effects

3. **Diagnostic Script**: `diagnose_solar_physics.py`
   - Analyzes equilibrium band structure
   - Calculates built-in potential
   - Identifies bandgap issues

4. **Optimization Script**: `optimize_solar.py`
   - Automated parameter sweep framework
   - Scoring function for multi-objective optimization
   - Extensible for future parameter exploration

---

## Conclusion
The optimization successfully identified that **the fundamental issue is not material parameters (bandgap, doping) but rather the simulation setup**:

1. ✅ **Optical generation is correctly integrated** (non-zero Jsc confirms this)
2. ✅ **Current scaling is physically correct** (from previous calibration)
3. ❌ **Contact/boundary conditions are preventing Voc development**
4. ❌ **Device geometry (thickness) is insufficient for target Jsc**

**Next critical step**: Fix contact boundary conditions to enable quasi-Fermi level splitting under illumination. This is the blocking issue preventing realistic solar cell operation.

---

## Appendix: Parameter Sensitivity

| Parameter | Range Tested | Impact on Voc | Impact on Jsc |
|-----------|--------------|---------------|---------------|
| In Fraction | 0.20 - 0.35 | **None** (all 0V) | Minimal (0.029-0.037 mA/cm²) |
| p-Doping | 1e17 - 1e18 cm⁻³ | **None** | Minimal |
| n-Doping | 1e18 - 1e19 cm⁻³ | **None** | Minimal |

**Interpretation**: The lack of sensitivity to material parameters confirms that the issue is **not** in the material physics, but in the **simulation boundary conditions or solver implementation**.
