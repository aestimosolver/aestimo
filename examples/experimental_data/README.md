# Experimental Data for Si p-n Junction Validation

This directory contains experimental reference data for validating Aestimo drift-diffusion simulations.

## Files

### si_pn_experimental_iv.csv
Experimental I-V characteristics for a Silicon p-n junction diode.

**Device Parameters:**
- Material: Silicon (Si)
- Doping: p-side = 1×10¹⁸ cm⁻³, n-side = 1×10¹⁸ cm⁻³
- Area: 1×10⁻⁴ cm² (100 × 100 μm)
- Temperature: 300 K

**Data Generation:**
This data is synthetically generated using the Shockley diode equation with realistic parameters:
- Saturation current density (Js): 1×10⁻¹² A/cm²
- Ideality factor (n): 1.05 (slight non-ideality from generation-recombination)
- Series resistance (Rs): 5 Ω
- Shunt resistance (Rsh): 10 kΩ

**Format:**
```
# Comments start with #
Voltage (V), Current (A)
-1.0,-1.002e-08
0.0,0.000e+00
0.5,2.340e-07
...
```

### si_pn_experimental_cv.csv
Experimental C-V (capacitance-voltage) characteristics for the same device.

**Data Generation:**
Generated using the abrupt junction C-V relationship:
```
C = A × sqrt(q × εs × Na × Nd / (2 × (Na + Nd) × (Vbi - V)))
```

Where:
- Vbi (built-in voltage): ~0.84 V for Si at 10¹⁸ cm⁻³ doping
- εs (Si permittivity): 11.68 × ε₀
- A: Device area

**Format:**
```
# Comments start with #
Voltage (V), Capacitance (F)
-1.0,4.123e-11
0.0,8.556e-11
...
```

## Usage

To validate your simulation against experimental data, use the example:

```python
# In your input file
enable_experimental_validation = True
experimental_iv_file = "examples/experimental_data/si_pn_experimental_iv.csv"
device_area = 1e-4  # cm²
```

See `examples/sample_pn_with_experimental_validation.py` for a complete example.

## Notes

1. **Synthetic Data**: The experimental data provided here is synthetically generated to match realistic Si p-n junction behavior. For actual research validation, replace with real experimental measurements.

2. **Device Matching**: Ensure your simulation parameters (doping, dimensions, temperature) match the experimental device to get meaningful comparisons.

3. **Units**: 
   - Voltage: Volts (V)
   - Current: Amperes (A)
   - Capacitance: Farads (F)
   - All data uses SI units

4. **Temperature Effects**: The data is for 300 K. For different temperatures, recalibrate or measure at the simulation temperature.

## Adding Your Own Data

To use your own experimental data:

1. Create a CSV file with the same format (voltage, current/capacitance)
2. Add comment lines starting with `#` to document the device parameters
3. Update the path in your simulation input file
4. Ensure device_area matches your experimental device

## References

- Shockley diode equation: W. Shockley, "The theory of p-n junctions in semiconductors and p-n junction transistors," Bell System Technical Journal, 1949
- C-V characteristics: S.M. Sze, "Physics of Semiconductor Devices", 3rd ed., Wiley, 2007
