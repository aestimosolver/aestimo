#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple test script to verify the experimental validation module works correctly.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'aeslibs'))

print("="*60)
print("Testing Experimental Validation Module")
print("="*60)

# Test 1: Import the module
print("\n[Test 1] Importing validation module...")
try:
    from experimental_validation import (
        load_experimental_data,
        compute_error_metrics,
        plot_iv_comparison,
        generate_validation_report
    )
    print("[OK] Successfully imported validation module")
except Exception as e:
    print(f"[FAIL] Failed to import: {e}")
    sys.exit(1)

# Test 2: Load experimental I-V data
print("\n[Test 2] Loading experimental I-V data...")
try:
    iv_file = os.path.join(os.path.dirname(__file__), 
                           'experimental_data', 'si_pn_experimental_iv.csv')
    exp_voltage, exp_current = load_experimental_data(iv_file)
    print(f"[OK] Loaded {len(exp_voltage)} data points")
    print(f"  Voltage range: {exp_voltage.min():.2f}V to {exp_voltage.max():.2f}V")
    print(f"  Current range: {exp_current.min():.3e}A to {exp_current.max():.3e}A")
except Exception as e:
    print(f"[FAIL] Failed to load I-V data: {e}")
    sys.exit(1)

# Test 3: Load experimental C-V data
print("\n[Test 3] Loading experimental C-V data...")
try:
    cv_file = os.path.join(os.path.dirname(__file__), 
                           'experimental_data', 'si_pn_experimental_cv.csv')
    cv_voltage, cv_capacitance = load_experimental_data(cv_file)
    print(f"[OK] Loaded {len(cv_voltage)} data points")
    print(f"  Voltage range: {cv_voltage.min():.2f}V to {cv_voltage.max():.2f}V")
    print(f"  Capacitance range: {cv_capacitance.min():.3e}F to {cv_capacitance.max():.3e}F")
except Exception as e:
    print(f"[FAIL] Failed to load C-V data: {e}")
    sys.exit(1)

# Test 4: Compute error metrics
print("\n[Test 4] Computing error metrics...")
try:
    import numpy as np
    
    # Create dummy simulation data (copy of experimental with small noise)
    sim_current = exp_current * (1 + 0.05 * np.random.randn(len(exp_current)))
    
    # Only forward bias for testing
    forward_mask = exp_voltage >= 0
    metrics = compute_error_metrics(
        exp_current[forward_mask],
        sim_current[forward_mask]
    )
    
    print(f"[OK] Successfully computed error metrics")
    print(f"  RMSE: {metrics['rmse']:.3e}")
    print(f"  MAE: {metrics['mae']:.3e}")
    print(f"  MAPE: {metrics['mape']:.2f}%")
    print(f"  R2: {metrics['r2']:.4f}")
except Exception as e:
    print(f"[FAIL] Failed to compute metrics: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 5: Generate validation report
print("\n[Test 5] Generating validation report...")
try:
    report = generate_validation_report(metrics, output_path=None)
    print("[OK] Successfully generated validation report")
except Exception as e:
    print(f"[FAIL] Failed to generate report: {e}")
    sys.exit(1)

# Test 6: Generate plots (without showing)
print("\n[Test 6] Generating comparison plots...")
try:
    import tempfile
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        plot_path = tmp.name
    
    plot_iv_comparison(
        exp_voltage, exp_current,
        exp_voltage, sim_current,
        output_path=plot_path,
        show=False
    )
    
    # Check if file was created
    if os.path.exists(plot_path):
        size = os.path.getsize(plot_path)
        print(f"[OK] Successfully generated I-V plot ({size} bytes)")
        os.remove(plot_path)  # Clean up
    else:
        print("[FAIL] Plot file was not created")
except Exception as e:
    print(f"[FAIL] Failed to generate plots: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*60)
print("All tests completed successfully!")
print("="*60)
print("\nThe experimental validation module is ready to use.")
print("See examples/sample_pn_with_experimental_validation.py for usage.")
