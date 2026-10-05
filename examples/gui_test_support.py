"""Synthetic current-sweep fixtures for GUI rendering tests, not physics validation."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np


class DeviceFigureFixtures:
    def setUp(self):
        super().setUp()
        temporary = TemporaryDirectory(prefix="aestimo-figures-")
        self.addCleanup(temporary.cleanup)
        self.addCleanup(plt.close, "all")
        self.fixture_root = Path(temporary.name)
        self.examples_dir = Path(__file__).resolve().parent
        self.gui = SimpleNamespace(examples_dir=str(self.examples_dir))

    def current_fixture(self, preset):
        """Write a known mA/cm² sweep independently of any solver or experimental CSV."""
        with (self.examples_dir / (preset + ".json")).open(encoding="utf-8") as stream:
            cfg = json.load(stream)
        cfg["validation_status"] = "MODEL-BASED / NOT EXPERIMENTALLY VALIDATED"
        output_dir = self.fixture_root / (preset + "_output")
        output_dir.mkdir()
        self.current_ma = np.linspace(0.0, 200.0, 201)
        self.voltage_v = np.linspace(0.0, 6.0, 201)
        area_cm2 = float(cfg["area"])
        np.savetxt(
            output_dir / "av_curr.dat",
            np.column_stack((self.voltage_v, self.current_ma / area_cm2)),
            header="Synthetic GUI test fixture: voltage_V current_density_mA_per_cm2",
        )
        return output_dir, cfg
