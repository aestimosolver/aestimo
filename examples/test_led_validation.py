# -*- coding: utf-8 -*-
"""
Automated Unit and Integration Test Suite for LED Validation, Characterization,
Provenance Traceability, and GUI Figure Builders in Aestimo 1D.
"""

from __future__ import annotations
import os
import sys
import json
import unittest
from pathlib import Path
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from aeslibs.characterize_led import (
    analyze_led_iv_curve,
    generate_electroluminescence_spectrum,
    compute_led_efficiency_droop_curve,
    HC_EV_NM,
)
from aeslibs.led_validation import (
    load_led_experimental_csv,
    compute_quantitative_error_metrics,
    plot_standardized_led_suite,
    generate_led_validation_report_markdown,
    LEDTraceabilityRecord,
)
from aestimo_gui import AestimoGUI
from examples.gui_test_support import DeviceFigureFixtures


class TestLEDCharacterization(unittest.TestCase):
    """Unit tests for core physical LED characterization routines."""

    def test_experimental_csv_loader(self):
        csv_files = [
            "led_nakamura1995_blue_sqw_iv.csv",
            "led_nakamura1995_blue_sqw_el_spectrum.csv",
            "led_meyaard2013_mqw_iv.csv",
            "led_meyaard2013_mqw_droop_iqe.csv",
            "led_schubert2006_algaas_dh_iv.csv",
            "led_schubert2006_algaas_dh_el_spectrum.csv",
        ]
        for fname in csv_files:
            fpath = REPO_ROOT / "examples" / "experimental_data" / fname
            self.assertTrue(fpath.exists(), f"Missing file: {fpath}")
            data, comments = load_led_experimental_csv(str(fpath))
            self.assertIsInstance(data, dict)
            self.assertGreater(len(data), 0)
            self.assertGreater(len(comments), 0)
            # Verify arrays are non-empty and finite
            for col_name, arr in data.items():
                self.assertIsInstance(arr, np.ndarray)
                self.assertGreater(len(arr), 0)
                self.assertTrue(np.all(np.isfinite(arr)))

    def test_error_metrics_calculation(self):
        y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        # Exact match
        metrics_exact = compute_quantitative_error_metrics(y_true, y_true)
        self.assertAlmostEqual(metrics_exact['rmse'], 0.0, places=6)
        self.assertAlmostEqual(metrics_exact['r2'], 1.0, places=6)
        self.assertAlmostEqual(metrics_exact['mape'], 0.0, places=6)

        # 10% offset
        y_pred = y_true * 1.1
        metrics_dev = compute_quantitative_error_metrics(y_true, y_pred)
        self.assertAlmostEqual(metrics_dev['mape'], 10.0, places=4)
        self.assertGreater(metrics_dev['rmse'], 0.0)

    def test_el_spectrum_generator(self):
        spec = generate_electroluminescence_spectrum(
            peak_wavelength_nm=450.0,
            fwhm_nm=25.0,
            temperature_k=300.0,
            wavelength_range_nm=(400.0, 500.0),
            num_points=101,
        )
        self.assertIn('wavelength_nm', spec)
        self.assertIn('intensity_norm', spec)
        self.assertAlmostEqual(spec['peak_wavelength_nm'], 450.0, delta=1.5)
        self.assertAlmostEqual(np.max(spec['intensity_norm']), 1.0, places=4)
        self.assertAlmostEqual(spec['fwhm_nm'], 25.0, delta=4.0)

    def test_iv_curve_analysis(self):
        # Synthetic diode I-V
        v = np.linspace(0.0, 3.5, 100)
        # Ideality n=1.5, Is=1e-12 A, Rs=10 Ohm, Area=1e-3 cm^2
        vt = 0.02585 * 1.5
        i_ma = 1e-9 * (np.exp(np.clip(v / vt, 0, 40)) - 1.0) * 1e3
        metrics = analyze_led_iv_curve(v, i_ma, device_area_cm2=1e-3, nominal_current_ma=20.0)
        self.assertIn('turn_on_voltage_v', metrics)
        self.assertIn('forward_voltage_at_nominal_v', metrics)
        self.assertGreater(metrics['turn_on_voltage_v'], 0.5)

    def test_efficiency_droop_abc(self):
        j_sweep = np.logspace(-1, 2.5, 80)
        droop = compute_led_efficiency_droop_curve(
            current_density_a_cm2=j_sweep,
            A_srh_s=1.0e7,
            B_rad_cm3_s=2.0e-11,
            C_auger_cm6_s=1.5e-30,
            d_active_cm=3.0e-7,
        )
        self.assertIn('current_density_a_cm2', droop)
        self.assertIn('normalized_iqe', droop)
        self.assertAlmostEqual(np.max(droop['normalized_iqe']), 1.0, places=4)
        # Peak density should be in the physically expected 5 - 25 A/cm^2 range
        self.assertGreater(droop['j_peak_a_cm2'], 3.0)
        self.assertLess(droop['j_peak_a_cm2'], 40.0)


class TestTraceabilityAndReports(unittest.TestCase):
    """Unit tests for provenance tracking, serialization, and report rendering."""

    def test_traceability_record_serialization(self):
        record = LEDTraceabilityRecord(
            device_id="TEST-LED-001",
            device_name="Test InGaN LED",
            validation_status="EXPERIMENTALLY VALIDATED",
        )
        record.set_bibliographic_reference({
            "title": "Test Title",
            "authors": "A. Tester et al.",
            "journal": "APL",
            "year": 2026,
            "doi": "10.1063/test",
        })
        record.set_experimental_structure({
            "mesa_area_cm2": 1.0e-3,
            "temperature_k": 300.0,
        })
        record.add_parameter_provenance(
            parameter_name="quantum_well_thickness",
            classification="EXPERIMENTAL",
            source="HRTEM calibration",
            value=3.0,
            unit="nm",
        )
        record.set_simulation_results({
            "turn_on_voltage_v": 2.7,
            "test_array": np.array([1.0, 2.0, 3.0]),
        })
        record.add_assumption("Light Extraction", "Planar escape cone")

        d = record.to_dict()
        self.assertEqual(d["device_id"], "TEST-LED-001")
        self.assertEqual(d["validation_status"], "EXPERIMENTALLY VALIDATED")
        self.assertIn("quantum_well_thickness", d["parameter_provenance"])

        # Test Markdown report generation
        md = generate_led_validation_report_markdown(record)
        self.assertIn("# Experimental Validation Report: Test InGaN LED", md)
        self.assertIn("`EXPERIMENTALLY VALIDATED`", md)
        self.assertIn("HRTEM calibration", md)


class TestGUIIntegration(DeviceFigureFixtures, unittest.TestCase):
    """Integration tests verifying GUI figure building for all 3 LED presets with synthetic current fixtures."""

    def test_build_led_figures_nakamura(self):
        out_dir, cfg = self.current_fixture("led_nakamura1995_blue_sqw")

        figs, titles, metrics, val_fig, val_report = AestimoGUI.build_led_figures(self.gui, str(out_dir), cfg)
        self.assertEqual(len(figs), 4)
        line = figs[0].axes[0].lines[1]
        np.testing.assert_allclose(line.get_ydata(), np.maximum(self.current_ma, 1e-7))
        np.testing.assert_allclose(
            line.get_xdata(), self.voltage_v + self.current_ma * 1e-3 * float(cfg.get("rs", 0.0))
        )
        self.assertIsNotNone(val_fig)
        self.assertIsNotNone(val_report)
        self.assertEqual(metrics["validation_status"], cfg["validation_status"])
        self.assertAlmostEqual(metrics["peak_wavelength_nm"], 450.0, delta=1.0)

    def test_build_led_figures_meyaard(self):
        out_dir, cfg = self.current_fixture("led_meyaard2013_blue_mqw")

        figs, titles, metrics, val_fig, val_report = AestimoGUI.build_led_figures(self.gui, str(out_dir), cfg)
        self.assertEqual(len(figs), 4)
        line = figs[0].axes[0].lines[1]
        np.testing.assert_allclose(line.get_ydata(), np.maximum(self.current_ma, 1e-7))
        np.testing.assert_allclose(
            line.get_xdata(), self.voltage_v + self.current_ma * 1e-3 * float(cfg.get("rs", 0.0))
        )
        self.assertIsNotNone(val_fig)
        self.assertIsNotNone(val_report)
        self.assertEqual(metrics["validation_status"], cfg["validation_status"])
        self.assertAlmostEqual(metrics["peak_wavelength_nm"], 445.0, delta=2.0)

    def test_build_led_figures_schubert(self):
        out_dir, cfg = self.current_fixture("led_schubert2006_algaas_dh")

        figs, titles, metrics, val_fig, val_report = AestimoGUI.build_led_figures(self.gui, str(out_dir), cfg)
        self.assertEqual(len(figs), 4)
        line = figs[0].axes[0].lines[1]
        np.testing.assert_allclose(line.get_ydata(), np.maximum(self.current_ma, 1e-7))
        np.testing.assert_allclose(
            line.get_xdata(), self.voltage_v + self.current_ma * 1e-3 * float(cfg.get("rs", 0.0))
        )
        self.assertIsNotNone(val_fig)
        self.assertIsNotNone(val_report)
        self.assertEqual(metrics["validation_status"], cfg["validation_status"])
        self.assertAlmostEqual(metrics["peak_wavelength_nm"], 870.0, delta=2.0)


    def test_missing_output_does_not_read_unrelated_working_directory(self):
        """An unrelated CWD result must never be shown as this device's result."""
        from unittest.mock import patch
        unrelated = self.fixture_root / "output"
        unrelated.mkdir()
        np.savetxt(unrelated / "av_curr.dat", [[0.0, 0.0], [1.0, 100.0]])
        missing = self.fixture_root / "missing_device"
        with patch("os.getcwd", return_value=str(self.fixture_root)):
            result = AestimoGUI.build_led_figures(self.gui, str(missing), {})
        self.assertEqual(result, ([], [], None, None, None))


class TestRepositoryAudit(unittest.TestCase):
    """Verifies that all examples in the repository have valid audit classifications."""

    def test_all_examples_audited(self):
        audit_file = REPO_ROOT / "examples" / "device_examples_audit.json"
        self.assertTrue(audit_file.exists(), "Missing device_examples_audit.json")
        with open(audit_file, "r", encoding="utf-8") as f:
            records = json.load(f)

        self.assertGreaterEqual(len(records), 40)
        valid_statuses = {
            "EXPERIMENTALLY VALIDATED",
            "PARTIALLY VALIDATED",
            "FITTED TO EXPERIMENT",
            "MODEL-BASED / NOT EXPERIMENTALLY VALIDATED",
        }
        for rec in records:
            self.assertIn(rec["validation_status"], valid_statuses, f"Invalid status in {rec['filename']}")
            self.assertTrue(len(rec["recommended_usage"]) > 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
