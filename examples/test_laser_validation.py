# -*- coding: utf-8 -*-
"""
Automated Unit and Integration Test Suite for Laser Diode Validation, Rate Equations,
Optical Cavity Physics, Provenance Traceability, and GUI Figure Builders in Aestimo 1D.
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

from aeslibs.characterize_laser import (
    calculate_cavity_optical_losses,
    calculate_threshold_gain_and_carrier_density,
    solve_laser_rate_equations_steady_state,
    extract_laser_figures_of_merit,
    generate_laser_fp_spectrum,
    compute_laser_temperature_series,
    HC_EV_NM,
)
from aeslibs.laser_validation import (
    load_laser_experimental_csv,
    compute_laser_error_metrics,
    plot_standardized_laser_suite,
    generate_laser_validation_report_markdown,
    LaserTraceabilityRecord,
)
from aestimo_gui import AestimoGUI
from examples.gui_test_support import DeviceFigureFixtures


class TestLaserCharacterization(unittest.TestCase):
    """Unit tests for physical laser cavity optics and rate-equation solver routines."""

    def test_experimental_csv_loader(self):
        csv_files = [
            "laser_tsang1981_gaas_sqw_li.csv",
            "laser_tsang1981_gaas_sqw_iv.csv",
            "laser_tsang1981_gaas_sqw_spectrum.csv",
            "laser_zah1994_1550nm_mqw_li.csv",
            "laser_zah1994_1550nm_mqw_iv.csv",
            "laser_zah1994_1550nm_mqw_spectrum.csv",
            "laser_nakamura1996_blue_mqw_li.csv",
            "laser_nakamura1996_blue_mqw_iv.csv",
            "laser_nakamura1996_blue_mqw_spectrum.csv",
        ]
        for fname in csv_files:
            fpath = REPO_ROOT / "examples" / "experimental_data" / fname
            self.assertTrue(fpath.exists(), f"Missing experimental dataset: {fpath}")
            data, comments = load_laser_experimental_csv(str(fpath))
            self.assertIsInstance(data, dict)
            self.assertGreater(len(data), 0)
            self.assertGreater(len(comments), 0)
            for col_name, arr in data.items():
                self.assertIsInstance(arr, np.ndarray)
                self.assertGreater(len(arr), 0)
                self.assertTrue(np.all(np.isfinite(arr)))

    def test_cavity_optical_losses(self):
        # 500 um GaAs cavity, R1=R2=0.32, alpha_i=10 cm^-1, n_g=3.6
        res = calculate_cavity_optical_losses(
            cavity_length_um=500.0,
            r1=0.32,
            r2=0.32,
            alpha_internal_cm1=10.0,
            group_index=3.6,
        )
        self.assertAlmostEqual(res['alpha_m_cm1'], 22.79, delta=0.5)
        self.assertAlmostEqual(res['alpha_tot_cm1'], 32.79, delta=0.5)
        self.assertGreater(res['photon_lifetime_s'], 1.0e-13)
        self.assertLess(res['photon_lifetime_s'], 1.0e-11)
        self.assertAlmostEqual(res['v_g_cm_s'], 2.99792458e10 / 3.6, delta=1e6)

    def test_rate_equations_steady_state(self):
        i_sweep = np.linspace(0.0, 80.0, 81)
        vol = 500e-4 * 20e-4 * 20e-7  # 500um x 20um x 20nm
        rate_res = solve_laser_rate_equations_steady_state(
            current_array_ma=i_sweep,
            active_volume_cm3=vol,
            A_srh_s=1e8,
            B_rad_cm3_s=1.5e-10,
            C_auger_cm6_s=3e-30,
            eta_i=0.95,
            photon_lifetime_s=3.6e-12,
            confinement_factor=0.038,
            g0_cm1=1800.0,
            n_tr_cm3=1.8e18,
            v_g_cm_s=8.3e9,
            alpha_m_cm1=22.8,
            beta_sp=1e-4,
            peak_wavelength_nm=845.0,
        )
        self.assertIn('power_single_facet_mw', rate_res)
        self.assertIn('carrier_density_cm3', rate_res)
        p_out = rate_res['power_single_facet_mw']
        # Monotonicity check above threshold
        self.assertGreater(p_out[-1], p_out[0])
        self.assertGreater(p_out[-1], 1.0)
        # Carrier density should clamp around threshold
        n_arr = rate_res['carrier_density_cm3']
        n_ratio = n_arr[-1] / n_arr[-2]
        self.assertAlmostEqual(n_ratio, 1.0, delta=0.02)

    def test_extract_figures_of_merit(self):
        i_arr = np.linspace(0.0, 50.0, 101)
        # Step-like above Ith=20 mA with SE=0.4 W/A
        p_arr = np.maximum(0.01, 0.40 * (i_arr - 20.0))
        p_arr[:40] = 0.01 * (i_arr[:40] / 20.0)
        v_arr = 1.40 + 0.01 * i_arr
        fom = extract_laser_figures_of_merit(
            current_ma=i_arr,
            power_mw=p_arr,
            voltage_v=v_arr,
            peak_wavelength_nm=845.0,
            area_cm2=1e-4,
            ith_calc_ma=20.0,
        )
        self.assertAlmostEqual(fom['threshold_current_ma'], 20.0, delta=1.0)
        self.assertAlmostEqual(fom['slope_efficiency_mw_per_ma'], 0.40, delta=0.05)
        self.assertGreater(fom['differential_quantum_efficiency'], 0.20)
        self.assertGreater(fom['max_optical_power_mw'], 10.0)

    def test_fabry_perot_spectrum(self):
        spec = generate_laser_fp_spectrum(
            peak_wavelength_nm=845.0,
            cavity_length_um=500.0,
            group_index=3.6,
            fwhm_sp_nm=25.0,
            fwhm_lasing_nm=1.8,
            current_ratio_i_over_ith=1.25,
            num_modes=31,
        )
        self.assertIn('wavelength_nm', spec)
        self.assertIn('intensity_norm', spec)
        self.assertAlmostEqual(spec['mode_spacing_nm'], 0.198, delta=0.02)
        self.assertAlmostEqual(np.max(spec['intensity_norm']), 1.0, places=3)
        self.assertEqual(len(spec['wavelength_nm']), 31)

    def test_temperature_series(self):
        i_arr = np.linspace(0.0, 60.0, 61)
        series = compute_laser_temperature_series(
            current_array_ma=i_arr,
            temp_array_k=[280.0, 300.0, 320.0],
            ith_ref_ma=20.0,
            t_ref_k=300.0,
            t0_k=150.0,
            slope_eff_ref=0.40,
        )
        self.assertIn("280K", series)
        self.assertIn("300K", series)
        self.assertIn("320K", series)
        # Ith must increase monotonically with T
        self.assertLess(series["280K"]["threshold_current_ma"], series["300K"]["threshold_current_ma"])
        self.assertLess(series["300K"]["threshold_current_ma"], series["320K"]["threshold_current_ma"])

    def test_error_metrics_calculation(self):
        y_true = np.array([1.0, 5.0, 10.0, 20.0])
        metrics_exact = compute_laser_error_metrics(y_true, y_true, target_name="P_opt")
        self.assertAlmostEqual(metrics_exact['rmse'], 0.0, places=6)
        self.assertAlmostEqual(metrics_exact['r2'], 1.0, places=6)


class TestLaserTraceabilityAndReports(unittest.TestCase):
    """Unit tests for laser provenance records, serialization, and report rendering."""

    def test_laser_record_serialization(self):
        record = LaserTraceabilityRecord(
            device_id="TEST-LASER-001",
            device_name="Test GaAs SQW Laser",
            laser_architecture="Ridge Waveguide SQW",
            material_system="Zincblende",
            validation_status="EXPERIMENTALLY VALIDATED",
        )
        record.set_bibliographic_reference({
            "title": "Extremely low threshold laser",
            "authors": "W. T. Tsang",
            "journal": "APL",
            "year": 1981,
            "doi": "10.1063/1.92690",
        })
        record.set_cavity_parameters({
            "cavity_length_um": 500.0,
            "stripe_width_um": 20.0,
            "r1": 0.32,
            "r2": 0.32,
            "alpha_i_cm1": 10.0,
            "alpha_m_cm1": 22.8,
            "confinement_factor": 0.038,
            "t0_k": 180.0,
        })
        record.set_experimental_structure({
            "active_area_cm2": 1e-4,
            "temperature_k": 300.0,
            "layers": [{"material": "GaAs", "thickness": 20.0, "type": "well"}],
        })
        record.add_parameter_provenance(
            parameter_name="cavity_length_um",
            classification="EXPERIMENTAL",
            source="Cleaved laser bar",
            value=500.0,
            unit="um",
        )
        record.set_simulation_results({
            "threshold_current_ma": 20.2,
            "slope_efficiency_mw_per_ma": 0.41,
            "current_ma": np.linspace(0, 40, 41),
            "power_single_facet_mw": np.linspace(0, 8, 41),
        })

        d = record.to_dict()
        self.assertEqual(d["device_id"], "TEST-LASER-001")
        self.assertEqual(d["validation_status"], "REFERENCE COMPARISON / PROVENANCE UNVERIFIED")
        self.assertIn("cavity_parameters", d)

        # Markdown report generation
        md = generate_laser_validation_report_markdown(record)
        self.assertIn("# Experimental Validation Report: Test GaAs SQW Laser", md)
        self.assertIn("`REFERENCE COMPARISON / PROVENANCE UNVERIFIED`", md)
        self.assertIn("W. T. Tsang", md)


class TestLaserGUIIntegration(DeviceFigureFixtures, unittest.TestCase):
    """Integration tests verifying GUI figure building for all 3 laser presets with synthetic current fixtures."""

    def test_build_laser_figures_tsang1981(self):
        out_dir, cfg = self.current_fixture("laser_tsang1981_gaas_sqw")

        figs, titles, metrics, val_fig, val_report = AestimoGUI.build_laser_figures(self.gui, str(out_dir), cfg)
        self.assertEqual(len(figs), 4)
        line = figs[1].axes[0].lines[0]
        np.testing.assert_allclose(line.get_ydata(), self.current_ma)
        np.testing.assert_allclose(
            line.get_xdata(), self.voltage_v + self.current_ma * 1e-3 * float(cfg.get("rs", 0.0))
        )
        self.assertIsNotNone(val_fig)
        self.assertIsNotNone(val_report)
        self.assertEqual(metrics["validation_status"], cfg["validation_status"])
        self.assertGreater(metrics["threshold_current_ma"], 0.0)
        self.assertLess(metrics["threshold_current_ma"], self.current_ma[-1])
        self.assertGreater(metrics["max_optical_power_mw"], 0.0)
        self.assertAlmostEqual(metrics["peak_wavelength_nm"], 845.0, delta=1.0)

    def test_build_laser_figures_zah1994(self):
        out_dir, cfg = self.current_fixture("laser_zah1994_1550nm_mqw")

        figs, titles, metrics, val_fig, val_report = AestimoGUI.build_laser_figures(self.gui, str(out_dir), cfg)
        self.assertEqual(len(figs), 4)
        line = figs[1].axes[0].lines[0]
        np.testing.assert_allclose(line.get_ydata(), self.current_ma)
        np.testing.assert_allclose(
            line.get_xdata(), self.voltage_v + self.current_ma * 1e-3 * float(cfg.get("rs", 0.0))
        )
        self.assertIsNotNone(val_fig)
        self.assertIsNotNone(val_report)
        self.assertEqual(metrics["validation_status"], cfg["validation_status"])
        self.assertGreater(metrics["threshold_current_ma"], 0.0)
        self.assertLess(metrics["threshold_current_ma"], self.current_ma[-1])
        self.assertGreater(metrics["max_optical_power_mw"], 0.0)
        self.assertAlmostEqual(metrics["peak_wavelength_nm"], 1550.0, delta=1.0)

    def test_build_laser_figures_nakamura1996(self):
        out_dir, cfg = self.current_fixture("laser_nakamura1996_blue_mqw")

        figs, titles, metrics, val_fig, val_report = AestimoGUI.build_laser_figures(self.gui, str(out_dir), cfg)
        self.assertEqual(len(figs), 4)
        line = figs[1].axes[0].lines[0]
        np.testing.assert_allclose(line.get_ydata(), self.current_ma)
        np.testing.assert_allclose(
            line.get_xdata(), self.voltage_v + self.current_ma * 1e-3 * float(cfg.get("rs", 0.0))
        )
        self.assertIsNotNone(val_fig)
        self.assertIsNotNone(val_report)
        self.assertEqual(metrics["validation_status"], cfg["validation_status"])
        self.assertGreater(metrics["threshold_current_ma"], 0.0)
        self.assertLess(metrics["threshold_current_ma"], self.current_ma[-1])
        self.assertGreater(metrics["max_optical_power_mw"], 0.0)
        self.assertAlmostEqual(metrics["peak_wavelength_nm"], 405.0, delta=1.0)


    def test_missing_output_does_not_read_unrelated_working_directory(self):
        """An unrelated CWD result must never be shown as this device's result."""
        from unittest.mock import patch
        unrelated = self.fixture_root / "output"
        unrelated.mkdir()
        np.savetxt(unrelated / "av_curr.dat", [[0.0, 0.0], [1.0, 100.0]])
        missing = self.fixture_root / "missing_device"
        with patch("os.getcwd", return_value=str(self.fixture_root)):
            result = AestimoGUI.build_laser_figures(self.gui, str(missing), {})
        self.assertEqual(result, ([], [], None, None, None))


class TestLaserRepositoryAudit(unittest.TestCase):
    """Verifies that all 3 laser configurations are validated in device_examples_audit.json."""

    def test_laser_examples_audited(self):
        audit_file = REPO_ROOT / "examples" / "device_examples_audit.json"
        self.assertTrue(audit_file.exists(), "Missing device_examples_audit.json")
        with open(audit_file, "r", encoding="utf-8") as f:
            records = json.load(f)

        laser_files = {
            "laser_tsang1981_gaas_sqw.json",
            "laser_zah1994_1550nm_mqw.json",
            "laser_nakamura1996_blue_mqw.json",
        }
        found_lasers = {}
        for rec in records:
            if rec["filename"] in laser_files:
                found_lasers[rec["filename"]] = rec

        self.assertEqual(len(found_lasers), 3, "Not all 3 laser devices found in audit!")
        for fname, rec in found_lasers.items():
            self.assertEqual(rec["validation_status"], "REFERENCE COMPARISON / PROVENANCE UNVERIFIED", f"{fname} not marked EXPERIMENTALLY VALIDATED")
            self.assertIn("Laser", rec["recommended_usage"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
