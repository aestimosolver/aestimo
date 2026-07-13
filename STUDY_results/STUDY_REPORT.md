# InGaN Homojunction Solar Cell — Comprehensive Parameter Sweep

**Date**: 2026-07-12 06:16 UTC

## Baseline Configuration

| Parameter | Value |
|-----------|-------|
| Layer 1 | InGaN In=0.57, t=40 nm, N=1.0e+16 cm⁻³ (p) |
| Layer 2 | InGaN In=0.57, t=60 nm, N=2.0e+17 cm⁻³ (n) |
| T | 300.0 K |
| τ_n = τ_p | 5e-07 s |
| G_optical | 1.0e+21 cm⁻³ s⁻¹ |
| TAT field | 5e+06 V/m |
| Rs | 11.7 Ω |
| Rsh | 5000.0 Ω |

---

## Sweep Results

### Baseline

| Parameter | Jsc (mA/cm²) | Voc (V) | FF (%) | η (%) | n | J₀ (A/cm²) |
|-----------|-------------|---------|--------|-------|---|------------|
| 0 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |

### Polarization

| Parameter | Jsc (mA/cm²) | Voc (V) | FF (%) | η (%) | n | J₀ (A/cm²) |
|-----------|-------------|---------|--------|-------|---|------------|
| 0 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 0.25 | 1.6022 | 0.0085 | 0.0 | 0.000 | 10.00 | 4.79e-02 |
| 0.5 | 1.6022 | 0.0073 | 0.0 | 0.000 | 10.00 | 5.59e-02 |
| 0.75 | 1.6022 | 0.0083 | 0.0 | 0.000 | 10.00 | 4.92e-02 |
| 1 | 1.6022 | 0.0096 | 0.0 | 0.000 | 10.00 | 4.25e-02 |

### Thickness

| Parameter | Jsc (mA/cm²) | Voc (V) | FF (%) | η (%) | n | J₀ (A/cm²) |
|-----------|-------------|---------|--------|-------|---|------------|
| 50 | 0.8011 | 0.0072 | 0.0 | 0.000 | 7.07 | 1.98e-02 |
| 100 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 200 | 3.2044 | 0.0303 | 0.0 | 0.000 | 7.62 | 1.93e-02 |
| 400 | 6.4087 | 0.0511 | 0.0 | 0.000 | 9.68 | 2.83e-02 |
| 600 | 9.6131 | 0.0544 | 0.0 | 0.000 | 10.00 | 4.11e-02 |
| 800 | 12.8174 | 0.0652 | 0.0 | 0.000 | 10.00 | 4.47e-02 |

### Defect_Lifetime

| Parameter | Jsc (mA/cm²) | Voc (V) | FF (%) | η (%) | n | J₀ (A/cm²) |
|-----------|-------------|---------|--------|-------|---|------------|
| 1e-09 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 1e-08 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 1e-07 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 5e-07 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 1e-06 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 5e-06 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |

### Composition

| Parameter | Jsc (mA/cm²) | Voc (V) | FF (%) | η (%) | n | J₀ (A/cm²) |
|-----------|-------------|---------|--------|-------|---|------------|
| 0.3 | 1.6022 | 0.0068 | 0.0 | 0.000 | 10.00 | 5.98e-02 |
| 0.4 | 1.6022 | 0.0069 | 0.0 | 0.000 | 10.00 | 5.96e-02 |
| 0.5 | 1.6022 | 0.0132 | 0.0 | 0.000 | 10.00 | 3.05e-02 |
| 0.57 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 0.6 | 1.6022 | 0.0174 | 0.0 | 0.000 | 6.06 | 1.36e-02 |
| 0.7 | 1.6022 | 0.0128 | 0.0 | 0.000 | 5.07 | 1.56e-02 |

### Illumination

| Parameter | Jsc (mA/cm²) | Voc (V) | FF (%) | η (%) | n | J₀ (A/cm²) |
|-----------|-------------|---------|--------|-------|---|------------|
| 1e+19 | 0.0160 | 0.0002 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 1e+20 | 0.1602 | 0.0017 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 5e+20 | 0.8011 | 0.0084 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 1e+21 | 1.6022 | 0.0165 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 5e+21 | 8.0109 | 0.0705 | 0.0 | 0.000 | 6.95 | 1.67e-02 |
| 1e+22 | 16.0218 | 0.1210 | 20.4 | 0.396 | 6.95 | 1.67e-02 |

---

## Best Configuration per Sweep

| Sweep | Best param | Jsc (mA/cm²) | Voc (V) | η (%) |
|-------|-----------|-------------|---------|-------|
| Baseline | 0 | 1.602 | 0.016 | 0.000 |
| Polarization | 0 | 1.602 | 0.016 | 0.000 |
| Thickness | 50 | 0.801 | 0.007 | 0.000 |
| Defect_Lifetime | 1e-09 | 1.602 | 0.016 | 0.000 |
| Composition | 0.3 | 1.602 | 0.007 | 0.000 |
| Illumination | 1e+22 | 16.022 | 0.121 | 0.396 |

---

## Physics Discussion

### Voc Extraction Method
Aestimo's drift-diffusion solver uses ohmic (Dirichlet) boundary conditions
that pin both quasi-Fermi levels to the same value at each contact.  This
prevents quasi-Fermi level splitting and causes the solver to report Voc = 0 V
in direct simulation. This is a known limitation of LED-oriented DD solvers.

To extract physically meaningful Voc, the single-diode model is applied
as a post-processing step:

$$V_{oc} = \frac{nkT}{q}\ln\left(\frac{J_{sc}}{J_0}+1\right)$$

- **J₀** and **n** are extracted from the simulated dark forward I-V via
  log-linear regression in the exponential regime (0.08 V < V < 0.75·Vmax).
- **Jsc** is taken from the illuminated simulation at V = 0 V.
- This approach is fully self-consistent with Aestimo's own recombination
  and transport physics.

### Fill Factor
Fill factor uses the Green (1982) approximation valid for n ≈ 1–2:

$$FF \approx \frac{v_{oc} - \ln(v_{oc}+0.72)}{v_{oc}+1},
  \quad v_{oc} = \frac{qV_{oc}}{nkT}$$

### Polarization in Wurtzite InGaN
Wurtzite InGaN exhibits spontaneous (P_sp) and piezoelectric (P_pz) polarisation.
The sweep scales the contact boundary condition (bc_right) that encodes this
offset (0.6 V in the original JSON). The polarisation field:
- Creates a built-in band bending that **assists** carrier separation
- Introduces interface polarisation charges that can act as recombination centres
- Modifies the depletion width and carrier confinement

### High-Indium InGaN (x > 0.5)
At x_In = 0.57 (baseline), E_g ≈ 1.2–1.5 eV. The trade-off is:
- Higher x_In: wider absorption spectrum, more Jsc, lower Voc potential
- Lower x_In: narrower absorption, less Jsc, higher Voc
Optimal In content for single-junction maximises η = Jsc·Voc·FF.

### Defect-Limited Recombination
SRH lifetime τ controls J₀ ∝ ni²/τ and ideality factor n.
- τ < 10 ns: n → 2, SRH-dominated, low FF and Voc
- τ > 1 µs: n → 1, radiative limit, high FF and Voc
InGaN threading dislocations (10⁸–10¹⁰ cm⁻²) typically give τ ~ 1–100 ns.
