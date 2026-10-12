"""State exports label relative effective masses consistently with their values."""
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import sys
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aeslibs import plotting


class TestStateExportUnits(unittest.TestCase):
    def test_equilibrium_and_bias_exports_identify_relative_mass(self):
        mass = np.array([[0.0], [0.067 * plotting.m_e], [0.0]])
        result = SimpleNamespace(
            xaxis=np.array([0.0, 1e-9]), N_wells_virtual=3,
            Well_boundary=np.zeros((3, 2)), meff_state_general=mass,
            meff_statec_general=mass, E_state_general=np.zeros((3, 1)),
            E_statec_general=np.zeros((3, 1)), N_state_general=np.zeros((3, 1)),
            N_statec_general=np.zeros((3, 1)), Total_Steps=1,
            Va_t=np.array([0.0]), av_curr=np.array([0.0]),
        )
        model = SimpleNamespace(n_max=2, subnumber_h=1, subnumber_e=1)
        controls = SimpleNamespace(
            sigma_out=False, electricfield_out=False, potential_out=False,
            states_out=True, probability_out=False, Drift_Diffusion_out=True,
        )
        with patch.object(plotting, 'config', controls), patch.object(
            plotting, 'amort_wave', return_value=(0, 2, 0, 2)
        ), TemporaryDirectory() as tmp:
            for exporter in (plotting.save_and_plot, plotting.save_and_plot2):
                with self.subTest(exporter=exporter.__name__):
                    target = Path(tmp) / exporter.__name__
                    exporter(result, model, str(target), drawFigures=False, show=False)
                    files = list(target.glob('states_*.dat'))
                    self.assertEqual(len(files), 2)
                    for path in files:
                        header = path.read_text().splitlines()[0]
                        self.assertIn('Subband m* (m_e)', header)
                        self.assertNotIn('(kg)', header)
                        self.assertAlmostEqual(np.loadtxt(path)[-1], 0.067)
