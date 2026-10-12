"""Plotting small devices must not allocate fixed 800 MB helper arrays."""
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import sys
import unittest
from unittest.mock import patch
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aeslibs import plotting


class TestPlotMemory(unittest.TestCase):
    def check_plot(self, exporter, expected_levels):
        n = 4
        profile = np.linspace(0.0, 1.0, n)
        result = SimpleNamespace(
            xaxis=profile * 1e-9, N_wells_virtual=3,
            Well_boundary=np.zeros((3, 2)),
            E_statec_general=np.array([[0.0], [7.0], [0.0]]),
            E_state_general=np.array([[0.0], [-3.0], [0.0]]),
            wfe_general=np.zeros((3, 1, n)), wfh_general=np.zeros((3, 1, n)),
            EF=5.0, fitot=profile * plotting.q, fitotc=profile * plotting.q,
            ro_result=profile, el_field1_result=profile, el_field2_result=profile,
            fi_result=profile, nf_result=np.full(n, 1e16), pf_result=np.full(n, 1e16),
            Total_Steps=1, Va_t=np.array([0.0]), av_curr=np.array([0.0]),
            Ec_result=profile, Ev_result=profile, Ei_result=profile,
            Efn_result=profile, Efp_result=profile,
        )
        model = SimpleNamespace(n_max=n, N_wells_virtual=3)
        controls = SimpleNamespace(
            sigma_out=False, electricfield_out=False, potential_out=False,
            states_out=False, probability_out=False, Drift_Diffusion_out=False,
            wavefunction_scalefactor=400,
        )
        original_ones = np.ones

        def bounded_ones(shape, *args, **kwargs):
            size = int(np.prod(shape)) if isinstance(shape, tuple) else int(shape)
            self.assertLess(size, 100_000_000, 'Fixed-size plot allocation regressed')
            return original_ones(shape, *args, **kwargs)

        try:
            with TemporaryDirectory() as tmp, patch.object(plotting, 'config', controls), \
                 patch.object(plotting, 'amort_wave', return_value=(0, n, 0, n)), \
                 patch.object(plotting.np, 'ones', side_effect=bounded_ones):
                figures = exporter(result, model, tmp, drawFigures=True, show=False)
                self.assertIsNotNone(figures[0])
                levels = [line for line in figures[0].axes[0].lines
                          if line.get_linestyle() == '--']
                self.assertEqual(len(levels), len(expected_levels))
                for line, expected in zip(levels, expected_levels):
                    np.testing.assert_array_equal(line.get_ydata(), np.full(n, expected))
        finally:
            plt.close('all')

    def test_equilibrium_plot_preserves_levels_without_large_allocation(self):
        self.check_plot(plotting.save_and_plot, [7.0, -3.0, 5.0])

    def test_transport_plot_preserves_levels_without_large_allocation(self):
        self.check_plot(plotting.save_and_plot2, [0.007, -0.003])
