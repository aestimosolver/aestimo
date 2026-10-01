#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
examples/benchmark_gui_speed.py
Comprehensive GUI performance and responsiveness benchmarking suite for Aestimo 1D.

Evaluates:
1. Project configuration loading latency across multiple device complexities
2. Structure diagram rendering and axes reuse latency
3. Result figure display latency (draw_idle vs synchronous draw)
4. Log poller CPU overhead during idle state
5. Material property caching performance
"""

import os
import sys
import time
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Ensure root directory is on sys.path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from aestimo_gui import AestimoGUI
import aeslibs.structure_diagram as sd


def benchmark_project_loading(app):
    print("\n" + "=" * 65)
    print(" 1. PROJECT CONFIGURATION LOADING BENCHMARK")
    print("=" * 65)

    test_configs = [
        ("Tobin 1990 GaAs (4 layers)", os.path.join(BASE_DIR, "examples", "gaas_tobin1990_benchmark.json")),
        ("Tsang 1981 SQW Laser (7 layers)", os.path.join(BASE_DIR, "examples", "laser_tsang1981_gaas_sqw.json")),
        ("Zah 1994 5-QW Laser (17 layers)", os.path.join(BASE_DIR, "examples", "laser_zah1994_1550nm_mqw.json")),
        ("Nakamura 1996 MQW Laser (14 layers)", os.path.join(BASE_DIR, "examples", "laser_nakamura1996_blue_mqw.json")),
        ("Meyaard 2013 MQW LED (14 layers)", os.path.join(BASE_DIR, "examples", "led_meyaard2013_blue_mqw.json")),
    ]

    results = {}
    for name, path in test_configs:
        if not os.path.exists(path):
            print(f"  [SKIP] {name}: file not found at {path}")
            continue

        with open(path, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        num_layers = len(cfg.get("layers", []))

        # Warmup
        app.load_configuration(cfg)
        app.update()

        # Timed runs
        times = []
        for _ in range(3):
            t0 = time.perf_counter()
            app.load_configuration(cfg)
            app.update()
            t1 = time.perf_counter()
            times.append(t1 - t0)

        avg_time = sum(times) / len(times)
        min_time = min(times)
        results[name] = {"layers": num_layers, "avg_s": avg_time, "min_s": min_time}
        print(f"  {name:<38} | {num_layers:>2} layers | Avg: {avg_time*1000:>6.1f} ms | Min: {min_time*1000:>6.1f} ms")

    return results


def benchmark_structure_diagram_rendering(app):
    print("\n" + "=" * 65)
    print(" 2. STRUCTURE DIAGRAM RENDERING & AXES REUSE BENCHMARK")
    print("=" * 65)

    fig = plt.Figure(figsize=(7, 5), dpi=100)
    cfg_path = os.path.join(BASE_DIR, "examples", "laser_nakamura1996_blue_mqw.json")
    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    layers = cfg.get("layers", [])

    # Initial render (Cold)
    t0 = time.perf_counter()
    sd.render_device_diagram(fig, layers, view_mode="Structure + Bands", device_type="Laser Diode")
    t_cold = time.perf_counter() - t0
    print(f"  Cold Render (Gridspec creation):       {t_cold*1000:>6.1f} ms")

    # Subsequent renders with axes reuse (Warm)
    warm_times = []
    for _ in range(10):
        t0 = time.perf_counter()
        sd.render_device_diagram(fig, layers, view_mode="Structure + Bands", device_type="Laser Diode")
        t1 = time.perf_counter()
        warm_times.append(t1 - t0)

    avg_warm = sum(warm_times) / len(warm_times)
    min_warm = min(warm_times)
    print(f"  Warm Render with Axes Reuse (10 runs): Avg: {avg_warm*1000:>6.1f} ms | Min: {min_warm*1000:>6.1f} ms")
    print(f"  Speedup factor from axes reuse:        {t_cold / avg_warm:>6.2f}x")

    return {"cold_s": t_cold, "warm_avg_s": avg_warm, "warm_min_s": min_warm}


def benchmark_figure_display(app):
    print("\n" + "=" * 65)
    print(" 3. MULTI-FIGURE DISPLAY PERFORMANCE (DRAW_IDLE)")
    print("=" * 65)

    # Generate 5 test figures
    test_figures = []
    for i in range(5):
        fig, ax = plt.subplots(figsize=(6, 4))
        x = [0.1 * k for k in range(100)]
        ax.plot(x, [k**0.5 for k in x], label="Curve A")
        ax.plot(x, [k**0.7 for k in x], label="Curve B")
        ax.set_title(f"Test Benchmark Figure {i+1}")
        test_figures.append(fig)

    times = []
    for _ in range(3):
        t0 = time.perf_counter()
        app.display_figures(test_figures)
        app.update_idletasks()
        t1 = time.perf_counter()
        times.append(t1 - t0)

    avg_time = sum(times) / len(times)
    print(f"  5 Figures Display & Tk Embedding:      Avg: {avg_time*1000:>6.1f} ms | Min: {min(times)*1000:>6.1f} ms")

    # Cleanup test figures
    for f in test_figures:
        plt.close(f)

    return {"display_5_figs_avg_s": avg_time}


def benchmark_log_poller(app):
    print("\n" + "=" * 65)
    print(" 4. ZERO-I/O LOG POLLER BENCHMARK")
    print("=" * 65)

    # Ensure log file exists with some content
    log_path = "aestimo.log"
    if not os.path.exists(log_path):
        with open(log_path, "w", encoding="utf-8") as f:
            for i in range(200):
                f.write(f"Log line {i}: initialization and step information\n")

    t0 = time.perf_counter()
    iterations = 100
    for _ in range(iterations):
        app.poll_log_file()
    t1 = time.perf_counter()

    elapsed = t1 - t0
    per_iter = (elapsed / iterations) * 1e6
    print(f"  100 Log Polls (No size change):        Total: {elapsed*1000:>6.2f} ms | {per_iter:>6.1f} µs/poll")

    return {"total_ms": elapsed * 1000, "us_per_poll": per_iter}


def benchmark_property_cache():
    print("\n" + "=" * 65)
    print(" 5. MATERIAL PROPERTY CACHE BENCHMARK")
    print("=" * 65)

    layer = {
        "material": "AlGaAs",
        "mole": 0.3,
        "mole_y": 0.0,
        "thickness": 15.0,
        "doping": 1e18,
        "doping_type": "n",
        "type": "barrier"
    }

    # Clear cache to measure cold
    sd._LAYER_PROPERTY_CACHE.clear()
    t0 = time.perf_counter()
    p1 = sd.extract_layer_properties(layer)
    t_cold = time.perf_counter() - t0

    # Warm cached lookups
    t0 = time.perf_counter()
    n_lookups = 10000
    for _ in range(n_lookups):
        p2 = sd.extract_layer_properties(layer)
    t_warm = time.perf_counter() - t0

    per_lookup = (t_warm / n_lookups) * 1e6
    print(f"  Cold Property Lookup:                  {t_cold*1e6:>6.1f} µs")
    print(f"  10,000 Cached Lookups:                 Total: {t_warm*1000:>6.2f} ms | {per_lookup:>6.2f} µs/lookup")
    print(f"  Speedup factor:                        {(t_cold * n_lookups) / t_warm:>6.1f}x")

    return {"cold_us": t_cold * 1e6, "cached_us": per_lookup}


def main():
    print("Starting Aestimo 1D GUI Performance Benchmark Suite...")
    app = AestimoGUI()
    app.withdraw()  # Run headless without opening visible window

    try:
        loading_res = benchmark_project_loading(app)
        rendering_res = benchmark_structure_diagram_rendering(app)
        display_res = benchmark_figure_display(app)
        poller_res = benchmark_log_poller(app)
        cache_res = benchmark_property_cache()

        print("\n" + "=" * 65)
        print(" ALL BENCHMARKS COMPLETED SUCCESSFULLY")
        print("=" * 65)
    finally:
        app.destroy()


if __name__ == "__main__":
    main()
