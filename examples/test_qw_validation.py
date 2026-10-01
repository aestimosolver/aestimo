#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit and Integration Test Suite for Quantum-Well (QW) Experimental Validation
Tests:
  - CSV dataset loading and metadata parsing
  - Statistical error metrics calculation
  - Miller et al. (1984) QCSE benchmark execution
  - Dingle (1975) quantum confinement benchmark execution
  - 6-panel validation suite plotting
  - Markdown validation report generation
"""

import os
import sys
import unittest
import numpy as np

# Ensure repository root is on sys.path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "examples" else SCRIPT_DIR
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use('Agg', force=True)

import aeslibs.quantum_well as qw
import aeslibs.qw_validation as qw_val


class TestQWValidation(unittest.TestCase):
    """Test suite covering quantum-well experimental validation capabilities."""

    def setUp(self):
        self.exp_dir = os.path.join(REPO_ROOT, "examples", "experimental_data")

    def test_load_qw_experimental_csv(self):
        """Test loading and metadata extraction for all three QW benchmark CSVs."""
        # 1. Miller 1984 QCSE
        miller_csv = os.path.join(self.exp_dir, "qw_miller1984_qcse_stark_shift.csv")
        self.assertTrue(os.path.exists(miller_csv), "Miller CSV must exist")
        data_m, meta_m = qw_val.load_qw_experimental_csv(miller_csv)
        self.assertIn("electric_field_kv_cm", data_m)
        self.assertIn("stark_shift_mev", data_m)
        self.assertEqual(len(data_m["electric_field_kv_cm"]), 12)
        self.assertEqual(data_m["electric_field_kv_cm"][0], 0.0)
        self.assertEqual(data_m["electric_field_kv_cm"][-1], 110.0)

        # 2. Dingle 1975 Confinement
        dingle_csv = os.path.join(self.exp_dir, "qw_dingle1975_energy_vs_width.csv")
        self.assertTrue(os.path.exists(dingle_csv), "Dingle CSV must exist")
        data_d, meta_d = qw_val.load_qw_experimental_csv(dingle_csv)
        self.assertIn("well_width_nm", data_d)
        self.assertIn("e1_hh1_transition_ev", data_d)
        self.assertEqual(len(data_d["well_width_nm"]), 11)

        # 3. Tsang 1981 Spectrum
        tsang_csv = os.path.join(self.exp_dir, "qw_tsang1981_sqw_photoluminescence.csv")
        self.assertTrue(os.path.exists(tsang_csv), "Tsang CSV must exist")
        data_t, meta_t = qw_val.load_qw_experimental_csv(tsang_csv)
        self.assertIn("wavelength_nm", data_t)
        self.assertIn("normalized_intensity", data_t)
        self.assertGreater(len(data_t["wavelength_nm"]), 10)

    def test_compute_qw_error_metrics(self):
        """Test mathematical accuracy of error metrics calculation."""
        # Exact match
        y_true = np.array([1.45, 1.48, 1.52, 1.58, 1.65])
        m_perf = qw_val.compute_qw_error_metrics(y_true, y_true, "test_perfect")
        self.assertAlmostEqual(m_perf["rmse"], 0.0, places=7)
        self.assertAlmostEqual(m_perf["mae"], 0.0, places=7)
        self.assertAlmostEqual(m_perf["r_squared"], 1.0, places=7)
        self.assertAlmostEqual(m_perf["pearson_r2"], 1.0, places=7)

        # Known offset (10 meV = 0.010 eV)
        y_offset = y_true + 0.010
        m_off = qw_val.compute_qw_error_metrics(y_true, y_offset, "test_offset")
        self.assertAlmostEqual(m_off["rmse"], 0.010, places=5)
        self.assertAlmostEqual(m_off["mae"], 0.010, places=5)
        self.assertAlmostEqual(m_off["pearson_r2"], 1.0, places=5)

    def test_miller1984_qcse_benchmark_execution(self):
        """Test full execution of Miller 1984 QCSE benchmark script."""
        import examples.qw_miller1984_qcse_benchmark as miller_bench
        res = miller_bench.run_benchmark()
        self.assertEqual(res["status"], "EXPERIMENTALLY VALIDATED")
        self.assertGreater(res["m_hh"]["pearson_r2"], 0.98)
        self.assertGreater(res["m_lh"]["pearson_r2"], 0.98)
        self.assertTrue(os.path.exists(res["figure_path"]), "Validation PNG must be generated")
        self.assertTrue(os.path.exists(res["report_path"]), "Validation report MD must be generated")

    def test_dingle1975_confinement_benchmark_execution(self):
        """Test full execution of Dingle 1975 confinement benchmark script."""
        import examples.qw_dingle1975_confinement_benchmark as dingle_bench
        res = dingle_bench.run_benchmark()
        self.assertEqual(res["status"], "EXPERIMENTALLY VALIDATED")
        self.assertGreater(res["m_hh1"]["pearson_r2"], 0.98)
        self.assertGreater(res["m_lh1"]["pearson_r2"], 0.97)
        self.assertGreater(res["m_hh2"]["pearson_r2"], 0.98)
        self.assertTrue(os.path.exists(res["figure_path"]), "Validation PNG must be generated")
        self.assertTrue(os.path.exists(res["report_path"]), "Validation report MD must be generated")


if __name__ == "__main__":
    unittest.main()
