#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit and Integration Tests for Publication-Quality Visual Layer-Stack Diagram
examples/test_structure_diagram.py
"""

import os
import sys
import unittest
from tempfile import TemporaryDirectory
import json
import numpy as np
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

# Ensure project root is in sys.path
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if base_dir not in sys.path:
    sys.path.insert(0, base_dir)

import aeslibs.structure_diagram as sd
import database

class TestStructureDiagram(unittest.TestCase):
    """Test suite for aeslibs/structure_diagram.py"""

    def setUp(self):
        self.examples_dir = os.path.join(base_dir, "examples")
        self.sample_layers = [
            {"material": "GaAs", "thickness": 200.0, "doping": 2e18, "doping_type": "p", "type": "barrier", "mole": 0.0},
            {"material": "AlGaAs", "thickness": 1000.0, "doping": 5e17, "doping_type": "p", "type": "barrier", "mole": 0.60},
            {"material": "AlGaAs", "thickness": 150.0, "doping": 0.0, "doping_type": "i", "type": "barrier", "mole": 0.25},
            {"material": "GaAs", "thickness": 8.0, "doping": 0.0, "doping_type": "i", "type": "well", "mole": 0.0},
            {"material": "AlGaAs", "thickness": 150.0, "doping": 0.0, "doping_type": "i", "type": "barrier", "mole": 0.25},
            {"material": "AlGaAs", "thickness": 1000.0, "doping": 5e17, "doping_type": "n", "type": "barrier", "mole": 0.60},
            {"material": "GaAs", "thickness": 2000.0, "doping": 2e18, "doping_type": "n", "type": "barrier", "mole": 0.0}
        ]

    def test_parse_and_validate_structure_valid(self):
        """Test parsing and property extraction for valid heterostructure."""
        valid, errors, warnings, parsed = sd.parse_and_validate_structure(self.sample_layers)
        self.assertTrue(valid)
        self.assertEqual(len(errors), 0)
        self.assertEqual(len(parsed), len(self.sample_layers))

        # Check property extraction for GaAs well (index 3)
        qw_props = parsed[3]
        self.assertAlmostEqual(qw_props["Eg"], 1.424, places=2)
        self.assertAlmostEqual(qw_props["thickness"], 8.0)
        self.assertEqual(qw_props["type"], "well")
        self.assertIn("Quantum Well", qw_props["role"])
        self.assertTrue(qw_props["n_refractive"] > 3.0)

    def test_parse_and_validate_structure_invalid(self):
        """Test consistency checker catches physical discrepancies."""
        # Empty structure
        v1, errs1, _, _ = sd.parse_and_validate_structure([])
        self.assertFalse(v1)
        self.assertIn("empty", errs1[0].lower())

        # Negative thickness
        bad_layers = [
            {"material": "GaAs", "thickness": -10.0, "doping": 1e18, "doping_type": "n", "type": "barrier", "mole": 0.0}
        ]
        v2, errs2, _, _ = sd.parse_and_validate_structure(bad_layers)
        self.assertFalse(v2)
        self.assertTrue(any("Thickness must be > 0" in e for e in errs2))

        # Invalid material
        bad_mat = [
            {"material": "Unobtainium", "thickness": 50.0, "doping": 1e18, "doping_type": "n", "type": "barrier", "mole": 0.0}
        ]
        v3, errs3, _, _ = sd.parse_and_validate_structure(bad_mat)
        self.assertFalse(v3)
        self.assertTrue(any("Unknown material" in e for e in errs3))

        # Mole fraction out of bounds
        bad_mole = [
            {"material": "AlGaAs", "thickness": 50.0, "doping": 1e18, "doping_type": "n", "type": "barrier", "mole": 1.5}
        ]
        v4, errs4, _, _ = sd.parse_and_validate_structure(bad_mole)
        self.assertFalse(v4)
        self.assertTrue(any("Mole fraction x must be in [0, 1]" in e for e in errs4))

    def test_dynamic_scaling_modes(self):
        """Test TRUE, SCHEMATIC, and AUTO scaling algorithms."""
        # AUTO should choose SCHEMATIC because ratio is 2000.0 / 8.0 = 250.0 > 15
        geom_auto = sd.calculate_display_geometry(self.sample_layers, scale_mode="AUTO")
        self.assertEqual(geom_auto["active_scale"], "SCHEMATIC SCALE")
        self.assertGreater(geom_auto["ratio"], 15.0)

        # In schematic mode, thin QW (8 nm) should have readable display thickness
        qw_geom = geom_auto["layers"][3]
        self.assertGreater(qw_geom["d_disp"], 2.0)  # Has visible width

        # Test TRUE SCALE
        geom_true = sd.calculate_display_geometry(self.sample_layers, scale_mode="TRUE SCALE")
        self.assertEqual(geom_true["active_scale"], "TRUE SCALE")
        # In true scale, total display should equal total span (100.0)
        self.assertAlmostEqual(geom_true["layers"][-1]["z1_disp"], 100.0, places=3)

    def test_mqw_detection(self):
        """Test periodic multi-quantum-well sequence detector."""
        # 5-period InGaN/GaN MQW
        mqw_layers = [
            {"material": "GaN", "thickness": 1000.0, "doping": 5e18, "doping_type": "n", "type": "barrier", "mole": 0.0}
        ]
        for _ in range(5):
            mqw_layers.append({"material": "InGaN", "thickness": 3.5, "doping": 0.0, "doping_type": "i", "type": "well", "mole": 0.15})
            mqw_layers.append({"material": "GaN", "thickness": 7.0, "doping": 0.0, "doping_type": "i", "type": "barrier", "mole": 0.0})
        mqw_layers.append({"material": "GaN", "thickness": 200.0, "doping": 2e18, "doping_type": "p", "type": "barrier", "mole": 0.0})

        regions = sd.detect_mqw_regions(mqw_layers)
        self.assertEqual(len(regions), 1)
        mqw = regions[0]
        self.assertEqual(mqw["period_count"], 5)
        self.assertEqual(mqw["well_material"], "InGaN")
        self.assertEqual(mqw["barrier_material"], "GaN")
        self.assertAlmostEqual(mqw["well_thick"], 3.5)
        self.assertAlmostEqual(mqw["barrier_thick"], 7.0)
        self.assertAlmostEqual(mqw["total_thickness"], 5 * (3.5 + 7.0))

    def test_chemical_formula_formatting(self):
        """Test publication formula formatting for binary, ternary, and quaternary."""
        f_gaas = sd.format_chemical_formula("GaAs", use_mathtext=False)
        self.assertEqual(f_gaas, "GaAs")

        f_algaas = sd.format_chemical_formula("AlGaAs", mole=0.60, use_mathtext=True)
        self.assertIn("\\mathrm{Al}_{0.60}", f_algaas)
        self.assertIn("\\mathrm{Ga}_{0.40}", f_algaas)

        f_ingan = sd.format_chemical_formula("InGaN", mole=0.15, use_mathtext=True)
        self.assertIn("\\mathrm{In}_{0.15}", f_ingan)

        f_quaternary = sd.format_chemical_formula("InGaAsP", mole=0.70, mole_y=0.82, use_mathtext=True)
        self.assertIn("\\mathrm{In}_{0.70}", f_quaternary)
        self.assertIn("\\mathrm{As}_{0.82}", f_quaternary)

    def test_hit_testing(self):
        """Test interactive click detection on canvas."""
        geom = sd.calculate_display_geometry(self.sample_layers, scale_mode="TRUE SCALE")

        class MockEvent:
            def __init__(self, x, y):
                self.xdata = x
                self.ydata = y

        # Click inside layer 0
        l0 = geom["layers"][0]
        ev0 = MockEvent(x=(l0["z0_disp"] + l0["z1_disp"]) / 2.0, y=0.5)
        self.assertEqual(sd.find_clicked_layer(ev0, geom), 0)

        # Click inside layer 3 (QW)
        l3 = geom["layers"][3]
        ev3 = MockEvent(x=(l3["z0_disp"] + l3["z1_disp"]) / 2.0, y=0.5)
        self.assertEqual(sd.find_clicked_layer(ev3, geom), 3)

        # Click out of bounds
        ev_out = MockEvent(x=999.0, y=0.5)
        self.assertIsNone(sd.find_clicked_layer(ev_out, geom))

    def test_render_all_benchmark_devices(self):
        """Test rendering across benchmark devices in all 3 view modes."""
        benchmark_files = [
            "gaas_tobin1990_benchmark.json",
            "laser_tsang1981_gaas_sqw.json",
            "laser_zah1994_1550nm_mqw.json",
            "laser_nakamura1996_blue_mqw.json",
            "led_meyaard2013_blue_mqw.json",
            "ingan_solar_cell.json"
        ]

        for bfile in benchmark_files:
            fpath = os.path.join(self.examples_dir, bfile)
            if not os.path.exists(fpath):
                continue
            with open(fpath, "r", encoding="utf-8") as f:
                proj_data = json.load(f)

            layers = proj_data.get("layers", [])
            dev_type = proj_data.get("device_type", "Generic Diode / LED")
            prov = proj_data.get("parameter_provenance", None)

            # Test all 3 view modes
            for vmode in ["Structure + Bands", "Layer Structure", "Band Diagram"]:
                fig = plt.figure(figsize=(8, 5))
                geom = sd.render_device_diagram(
                    fig,
                    layers,
                    view_mode=vmode,
                    scale_mode="AUTO",
                    show_annotations=True,
                    theme="Dark (GUI)",
                    device_type=dev_type,
                    provenance_dict=prov
                )
                self.assertIsNotNone(geom)
                self.assertEqual(len(geom["layers"]), len(layers))
                plt.close(fig)

    def test_publication_export_formats(self):
        """Test exporting publication figure to PNG (300 DPI), SVG, and PDF."""
        fig = plt.figure(figsize=(8, 5))
        sd.render_device_diagram(fig, self.sample_layers, view_mode="Structure + Bands", theme="Publication (Light)")

        temporary = TemporaryDirectory(prefix="aestimo-export-test-")
        self.addCleanup(temporary.cleanup)
        export_dir = temporary.name
        png_path = os.path.join(export_dir, "test_pub_export.png")
        svg_path = os.path.join(export_dir, "test_pub_export.svg")
        pdf_path = os.path.join(export_dir, "test_pub_export.pdf")

        try:
            sd.export_publication_figure(fig, png_path, dpi=300)
            self.assertTrue(os.path.exists(png_path))
            self.assertGreater(os.path.getsize(png_path), 5000)

            sd.export_publication_figure(fig, svg_path)
            self.assertTrue(os.path.exists(svg_path))
            self.assertGreater(os.path.getsize(svg_path), 5000)

            sd.export_publication_figure(fig, pdf_path)
            self.assertTrue(os.path.exists(pdf_path))
            self.assertGreater(os.path.getsize(pdf_path), 5000)
        finally:
            plt.close(fig)
            for p in (png_path, svg_path, pdf_path):
                if os.path.exists(p):
                    os.remove(p)

    def test_provenance_badge_integration(self):
        """Test 4-tier provenance badge assignment."""
        prov_dict = {
            "general": "EXPERIMENTAL",
            "layers": {
                "GaAs": "EXPERIMENTAL",
                "AlGaAs": "DERIVED"
            }
        }
        l1 = {"material": "GaAs", "thickness": 10.0, "doping": 1e18, "doping_type": "n"}
        l2 = {"material": "AlGaAs", "thickness": 20.0, "doping": 1e17, "doping_type": "p", "mole": 0.3}

        p1 = sd.extract_layer_properties(l1, provenance_dict=prov_dict)
        self.assertEqual(p1["provenance"], "EXPERIMENTAL")

        p2 = sd.extract_layer_properties(l2, provenance_dict=prov_dict)
        self.assertEqual(p2["provenance"], "DERIVED")


if __name__ == "__main__":
    unittest.main()
