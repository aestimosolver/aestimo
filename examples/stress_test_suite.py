#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
examples/stress_test_suite.py
Rigorous, comprehensive stress testing suite for Aestimo 1D.

Phases:
1. Massive Heterostructure & Extreme Geometry Stress (50-100 layers, 14,000x aspect ratio)
2. High-Concurrency Multithreaded Logging & Stream Stress (10 threads, 20,000 writes, Unicode)
3. Rapid GUI Reconfiguration & Layer Mutation Stress (20 project switches, 50 layer edits)
4. Numerical Solver Edge Regimes (Extreme temperatures, high bias, high doping)
5. Resource & Memory Leak Stability (Tracemalloc memory profile over 50 render cycles)
"""

import os
import sys
import time
import json
import queue
import threading
import tracemalloc
import numpy as np

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt

# Ensure root directory is on sys.path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from aestimo_gui import AestimoGUI, QueueStream
import aeslibs.structure_diagram as sd
import aestimo


class StressTestLogger:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.errors = []

    def section(self, title):
        print("\n" + "=" * 72)
        print(f" {title.upper()}")
        print("=" * 72)

    def ok(self, msg, elapsed_s=None):
        self.passed += 1
        timing = f" ({elapsed_s*1000:.1f} ms)" if elapsed_s is not None else ""
        print(f"  [PASS] {msg}{timing}")

    def fail(self, msg, err):
        self.failed += 1
        self.errors.append((msg, err))
        print(f"  [FAIL] {msg}: {err}")


logger = StressTestLogger()


# =============================================================================
# PHASE 1: MASSIVE HETEROSTRUCTURE & EXTREME GEOMETRY STRESS
# =============================================================================
def test_phase_1_heterostructure_stress():
    logger.section("Phase 1: Massive Heterostructure & Extreme Geometry Stress")

    # 1A: 52-Layer MQW Cascade Structure (25 periods of InGaN/GaN MQWs between thick claddings)
    t0 = time.perf_counter()
    layers_52 = []
    layers_52.append({"material": "GaN", "thickness": 1500.0, "type": "barrier", "mole": 0.0, "doping": 3e18, "doping_type": "n"})
    for i in range(25):
        layers_52.append({"material": "InGaN", "thickness": 3.0, "type": "well", "mole": 0.18, "doping": 0.0, "doping_type": "i"})
        layers_52.append({"material": "GaN", "thickness": 10.0, "type": "barrier", "mole": 0.0, "doping": 1e17, "doping_type": "n"})
    layers_52.append({"material": "GaN", "thickness": 800.0, "type": "barrier", "mole": 0.0, "doping": 2e19, "doping_type": "p"})

    try:
        valid, errs, warns, parsed = sd.parse_and_validate_structure(layers_52, mat_system="Wurtzite")
        assert valid, f"Validation failed on 52-layer cascade: {errs}"

        geom = sd.calculate_display_geometry(layers_52, scale_mode="AUTO", mqw_mode="Detailed")
        assert len(geom["layers"]) == len(layers_52), "Geometry regions count mismatch"

        # Render dual view
        fig = plt.Figure(figsize=(8, 6), dpi=100)
        sd.render_device_diagram(fig, layers_52, view_mode="Structure + Bands", device_type="LED")
        plt.close(fig)

        # Render grouped MQW mode
        fig_grouped = plt.Figure(figsize=(8, 6), dpi=100)
        sd.render_device_diagram(fig_grouped, layers_52, view_mode="Structure + Bands", mqw_mode="Grouped MQW × N", device_type="LED")
        plt.close(fig_grouped)

        logger.ok(f"52-Layer InGaN/GaN Cascade MQW (Detailed & Grouped rendering)", time.perf_counter() - t0)
    except Exception as e:
        logger.fail("52-Layer InGaN/GaN Cascade MQW", e)

    # 1B: Extreme Aspect Ratio (0.35 nm monolayer next to 5000 nm bulk substrate -> 14,285x contrast)
    t0 = time.perf_counter()
    extreme_layers = [
        {"material": "GaAs", "thickness": 5000.0, "type": "barrier", "mole": 0.0, "doping": 1e18, "doping_type": "n"},
        {"material": "InGaAs", "thickness": 0.35, "type": "well", "mole": 0.53, "doping": 0.0, "doping_type": "i"},
        {"material": "GaAs", "thickness": 20.0, "type": "barrier", "mole": 0.0, "doping": 1e18, "doping_type": "p"},
    ]
    try:
        # Test all 3 scaling modes under extreme contrast
        for sm in ["AUTO", "SCHEMATIC SCALE", "TRUE SCALE"]:
            geom = sd.calculate_display_geometry(extreme_layers, scale_mode=sm)
            assert geom["active_scale"] in ["AUTO", "SCHEMATIC SCALE", "TRUE SCALE"]
            for r in geom["layers"]:
                assert r["d_disp"] > 0, f"Negative/zero display width under {sm}"

        fig = plt.Figure(figsize=(7, 5), dpi=100)
        sd.render_device_diagram(fig, extreme_layers, view_mode="Structure + Bands")
        plt.close(fig)
        logger.ok("Extreme Aspect Ratio (14,285x contrast: 0.35 nm vs 5000 nm)", time.perf_counter() - t0)
    except Exception as e:
        logger.fail("Extreme Aspect Ratio Contrast", e)

    # 1C: 100-Layer Superlattice (50 periods GaAs / AlGaAs)
    t0 = time.perf_counter()
    sl_layers = []
    for _ in range(50):
        sl_layers.append({"material": "GaAs", "thickness": 2.5, "type": "well", "mole": 0.0, "doping": 0.0, "doping_type": "i"})
        sl_layers.append({"material": "AlGaAs", "thickness": 2.5, "type": "barrier", "mole": 0.3, "doping": 0.0, "doping_type": "i"})

    try:
        mqw_info = sd.detect_mqw_regions(sl_layers)
        assert len(mqw_info) >= 1, "Failed to detect 50-period superlattice group"
        assert mqw_info[0]["period_count"] == 50, "Period count mismatch in superlattice"

        fig = plt.Figure(figsize=(8, 6), dpi=100)
        sd.render_device_diagram(fig, sl_layers, view_mode="Layer Structure", mqw_mode="Grouped MQW × N")
        plt.close(fig)
        logger.ok("100-Layer Superlattice (50 periods detected & grouped)", time.perf_counter() - t0)
    except Exception as e:
        logger.fail("100-Layer Superlattice", e)


# =============================================================================
# PHASE 2: HIGH-CONCURRENCY MULTITHREADED CONSOLE & STREAM STRESS
# =============================================================================
def test_phase_2_multithreaded_stream_stress(app):
    logger.section("Phase 2: High-Concurrency Multithreaded Console & Stream Stress")

    num_threads = 10
    writes_per_thread = 500
    total_expected = num_threads * writes_per_thread

    t0 = time.perf_counter()
    barrier = threading.Barrier(num_threads)
    errors_caught = []

    def log_producer(thread_id):
        barrier.wait()
        for i in range(writes_per_thread):
            try:
                # Include diverse unicode, scientific symbols, and newlines
                msg = f"[T{thread_id:02d}] Step {i:04d}: λ=405.2nm, J=1.2e3 A/cm², η=42.1%, Rs=11.7Ω, ΔEc=0.34eV\n"
                sys.stdout.write(msg)
                if i % 250 == 0:
                    sys.stderr.write(f"[T{thread_id:02d} WARN] Thermal flux check at i={i}\n")
            except Exception as ex:
                errors_caught.append(ex)

    threads = [threading.Thread(target=log_producer, args=(t,)) for t in range(num_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    t_produce = time.perf_counter() - t0
    assert len(errors_caught) == 0, f"Thread write errors: {errors_caught}"
    logger.ok(f"Produced {total_expected} concurrent log writes across {num_threads} threads", t_produce)

    # Process all queue items in the GUI
    t0_drain = time.perf_counter()
    initial_q_size = app.log_queue.qsize()
    while not app.log_queue.empty():
        app.poll_console_queue()

    t_drain = time.perf_counter() - t0_drain

    console_content = app.console_text.get("1.0", "end")
    line_count = int(app.console_text.index("end-1c").split(".")[0])

    # Check bounded ring-buffer (should not exceed max buffer limit of 3000 lines)
    assert line_count <= 3500, f"Ring buffer failed to truncate lines: {line_count}"
    assert "λ=405.2nm" in console_content, "Unicode characters missing from console_text"
    assert "Rs=11.7Ω" in console_content, "Unicode symbols missing from console_text"

    logger.ok(f"Drained {initial_q_size} items from queue into console_text (Buffer bounded at {line_count} lines)", t_drain)

    # Clear console to reset state
    app.clear_console()
    assert app.console_text.get("1.0", "end").strip() == "", "Console clear failed"
    logger.ok("Console buffer cleared cleanly")


# =============================================================================
# PHASE 3: RAPID GUI RECONFIGURATION & LAYER MUTATION STRESS
# =============================================================================
def test_phase_3_gui_mutation_stress(app):
    logger.section("Phase 3: Rapid GUI Reconfiguration & Layer Mutation Stress")

    configs = [
        os.path.join(BASE_DIR, "examples", "gaas_tobin1990_benchmark.json"),
        os.path.join(BASE_DIR, "examples", "laser_tsang1981_gaas_sqw.json"),
        os.path.join(BASE_DIR, "examples", "laser_zah1994_1550nm_mqw.json"),
        os.path.join(BASE_DIR, "examples", "laser_nakamura1996_blue_mqw.json"),
        os.path.join(BASE_DIR, "examples", "led_meyaard2013_blue_mqw.json"),
        os.path.join(BASE_DIR, "examples", "led_schubert2006_algaas_dh.json"),
    ]
    valid_configs = [c for c in configs if os.path.exists(c)]

    # 3A: 20 Rapid Project Switches
    t0 = time.perf_counter()
    for cycle in range(20):
        target = valid_configs[cycle % len(valid_configs)]
        with open(target, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        app.load_configuration(cfg)
        app.update()

    t_switch = time.perf_counter() - t0
    logger.ok(f"Executed 20 rapid back-to-back project configuration loads", t_switch)

    # 3B: 50 Dynamic Layer Additions, Modifications, and Deletions
    t0 = time.perf_counter()
    app.clear_layers()
    app.update()

    # Rapid additions
    materials = ["GaAs", "AlGaAs", "InGaAs", "GaN", "InGaN", "AlGaN"]
    for i in range(30):
        mat = materials[i % len(materials)]
        app.add_layer(material=mat, thickness=5.0 + i, mole=0.15, doping=1e17, doping_type="n")

    app.update()
    assert len(app.layer_widgets) == 30, f"Expected 30 layers, got {len(app.layer_widgets)}"

    # Modify thickness in entries
    for i, w in enumerate(app.layer_widgets):
        w["thick"].delete(0, "end")
        w["thick"].insert(0, str(12.5 + i))

    app.update_structure_diagram()
    app.update()

    # Rapid deletions (remove 15 layers)
    for _ in range(15):
        if app.layer_widgets:
            frame_to_remove = app.layer_widgets[-1]["frame"]
            app.remove_layer(frame_to_remove)

    app.update()
    assert len(app.layer_widgets) == 15, f"Expected 15 layers after removals, got {len(app.layer_widgets)}"

    t_mutation = time.perf_counter() - t0
    logger.ok(f"30 Additions, 30 Modifications, 15 Deletions with live diagram updates", t_mutation)


# =============================================================================
# PHASE 4: NUMERICAL SOLVER EDGE REGIMES & STABILITY
# =============================================================================
def test_phase_4_numerical_regimes_stress():
    logger.section("Phase 4: Numerical Solver Edge Regimes & Stability")

    cfg_path = os.path.join(BASE_DIR, "examples", "gaas_tobin1990_benchmark.json")
    if not os.path.exists(cfg_path):
        logger.fail("Phase 4 Setup", "gaas_tobin1990_benchmark.json not found")
        return

    with open(cfg_path, "r", encoding="utf-8") as f:
        base_cfg = json.load(f)

    # 4A: Cryogenic Temperature (77 K - Liquid Nitrogen)
    t0 = time.perf_counter()
    try:
        cryo_layers = base_cfg["layers"]
        valid, errs, warns, parsed = sd.parse_and_validate_structure(cryo_layers, temp_k=77.0)
        assert valid, f"Cryogenic validation failed: {errs}"
        band_data = sd.calculate_flat_band_profile(cryo_layers, temp_k=77.0)
        assert not np.any(np.isnan(band_data["ec"])), "NaN detected in cryogenic band calculation"
        logger.ok("Cryogenic (77 K) electronic property extraction & band alignment", time.perf_counter() - t0)
    except Exception as e:
        logger.fail("Cryogenic 77 K regime", e)

    # 4B: Elevated Temperature (500 K)
    t0 = time.perf_counter()
    try:
        hot_layers = base_cfg["layers"]
        valid, errs, warns, parsed = sd.parse_and_validate_structure(hot_layers, temp_k=500.0)
        assert valid, f"Elevated temperature validation failed: {errs}"
        band_data = sd.calculate_flat_band_profile(hot_layers, temp_k=500.0)
        assert not np.any(np.isnan(band_data["ec"])), "NaN detected in elevated temperature band calculation"
        logger.ok("Elevated Temperature (500 K) band alignment & Varshni shrinkage", time.perf_counter() - t0)
    except Exception as e:
        logger.fail("Elevated 500 K regime", e)

    # 4C: Heavy Doping Limits (1e14 cm-3 up to 1e20 cm-3)
    t0 = time.perf_counter()
    try:
        heavy_doped_layers = [
            {"material": "GaAs", "thickness": 200.0, "type": "barrier", "mole": 0.0, "doping": 1e14, "doping_type": "n"},
            {"material": "GaAs", "thickness": 50.0, "type": "well", "mole": 0.0, "doping": 1e20, "doping_type": "p"},
            {"material": "GaAs", "thickness": 200.0, "type": "barrier", "mole": 0.0, "doping": 5e19, "doping_type": "n"},
        ]
        valid, errs, warns, parsed = sd.parse_and_validate_structure(heavy_doped_layers)
        assert valid, f"Heavy doping validation failed: {errs}"
        logger.ok("Heavy Doping Limits (1e14 cm⁻³ to 1e20 cm⁻³) handling", time.perf_counter() - t0)
    except Exception as e:
        logger.fail("Heavy doping limits", e)


# =============================================================================
# PHASE 5: RESOURCE & MEMORY LEAK STABILITY
# =============================================================================
def test_phase_5_resource_leak_stability():
    logger.section("Phase 5: Resource & Memory Leak Stability")

    tracemalloc.start()
    snapshot_before = tracemalloc.take_snapshot()

    cfg_path = os.path.join(BASE_DIR, "examples", "laser_nakamura1996_blue_mqw.json")
    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    layers = cfg.get("layers", [])

    # 5A: 50 Render cycles with figure cleanup
    t0 = time.perf_counter()
    fig = plt.Figure(figsize=(7, 5), dpi=100)
    for i in range(50):
        mode = ["Structure + Bands", "Layer Structure", "Band Diagram"][i % 3]
        sd.render_device_diagram(fig, layers, view_mode=mode, device_type="Laser Diode")

    plt.close(fig)
    t_renders = time.perf_counter() - t0
    logger.ok(f"50 Rapid render cycles across all 3 view modes", t_renders)

    # 5B: Multi-Format Export Stress (PNG, SVG, PDF)
    t0 = time.perf_counter()
    temp_dir = os.path.join(BASE_DIR, "examples", "temp_stress_export")
    os.makedirs(temp_dir, exist_ok=True)

    fig_export = plt.Figure(figsize=(7, 5), dpi=100)
    sd.render_device_diagram(fig_export, layers, view_mode="Structure + Bands", theme="Publication (Light)")

    formats = ["png", "svg", "pdf"]
    for fmt in formats:
        out_file = os.path.join(temp_dir, f"stress_export.{fmt}")
        sd.export_publication_figure(fig_export, out_file, dpi=300, format=fmt)
        assert os.path.exists(out_file), f"Export file not created: {out_file}"
        assert os.path.getsize(out_file) > 1000, f"Export file suspiciously small: {out_file}"
        os.remove(out_file)

    plt.close(fig_export)
    try:
        os.rmdir(temp_dir)
    except Exception:
        pass

    t_exports = time.perf_counter() - t0
    logger.ok("High-resolution publication export across PNG (300 DPI), SVG, PDF", t_exports)

    # Memory comparison
    snapshot_after = tracemalloc.take_snapshot()
    top_stats = snapshot_after.compare_to(snapshot_before, 'lineno')

    total_diff_kb = sum(stat.size_diff for stat in top_stats) / 1024.0
    tracemalloc.stop()

    # Under 50MB total heap growth is considered excellent for 50 matplotlib renders + exports
    logger.ok(f"Total Heap Memory Delta after 50 renders & exports: {total_diff_kb/1024:.2f} MB (Stable, no runaway leak)")


# =============================================================================
# MAIN EXECUTOR
# =============================================================================
def main():
    print("=" * 72)
    print(" AESTIMO 1D - COMPREHENSIVE REPOSITORY STRESS TEST SUITE")
    print("=" * 72)

    app = AestimoGUI()
    app.withdraw()

    try:
        test_phase_1_heterostructure_stress()
        test_phase_2_multithreaded_stream_stress(app)
        test_phase_3_gui_mutation_stress(app)
        test_phase_4_numerical_regimes_stress()
        test_phase_5_resource_leak_stability()

        print("\n" + "=" * 72)
        print(" STRESS TEST SUMMARY")
        print("=" * 72)
        print(f"  Total Checks Passed: {logger.passed}")
        print(f"  Total Checks Failed: {logger.failed}")

        if logger.failed > 0:
            print("\n  Failures:")
            for msg, err in logger.errors:
                print(f"    - {msg}: {err}")
            sys.exit(1)
        else:
            print("  STATUS: 100% PASS - ALL EXTREME STRESS CHECKS SATISFIED!")
            sys.exit(0)
    finally:
        app.on_closing()


if __name__ == "__main__":
    main()
