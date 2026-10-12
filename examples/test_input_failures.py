"""Rejected inputs fail visibly instead of exiting a worker or reporting success."""
import json
import os
from pathlib import Path
import queue
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import numpy as np
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import aestimo
import database


class TestInputFailures(unittest.TestCase):
    def test_unsupported_and_unknown_schemes_fail_before_structure_arrays(self):
        for scheme, error in [(1, NotImplementedError), (3, NotImplementedError),
                              (4, NotImplementedError), (5, NotImplementedError),
                              (6, NotImplementedError), (99, ValueError)]:
            with self.subTest(scheme=scheme), patch.object(aestimo.Structure, 'create_structure_arrays') as arrays:
                with self.assertRaises(error):
                    aestimo.StructureFrom(SimpleNamespace(comp_scheme=scheme), database)
                arrays.assert_not_called()

    def test_grid_limit_raises_a_catchable_error(self):
        source = SimpleNamespace(comp_scheme=2, gridfactor=1.0, maxgridpoints=1,
                                 material=[[10.0, 'GaAs', 0.0, 0.0, 0.0, 'n', 'w']])
        with self.assertRaisesRegex(ValueError, 'exceeding maxgridpoints=1'):
            aestimo.StructureFrom(source, database)

    def test_carrier_energy_guard_raises_instead_of_exiting_the_caller(self):
        from aeslibs.aestimo_poisson1d import equi_np_fi2
        zeros = np.zeros(1)
        ones = np.ones(1)
        with self.assertRaisesRegex(ValueError, 'Carrier energy offset at grid point 0'):
            equi_np_fi2(zeros, zeros, ones, ones, np.array([-1.0]), zeros,
                        1, zeros, 0.025, None, None, SimpleNamespace(),
                        None, None, None, None, zeros, zeros, 1, ones)

    def run_cli(self, directory, text, extension='.py'):
        source = Path(directory) / ('rejected' + extension)
        source.write_text(text, encoding='utf-8')
        return subprocess.run([sys.executable, '-m', 'aestimo', '-i', str(source)],
                              cwd=directory, env=dict(os.environ, PYTHONPATH=str(ROOT), MPLBACKEND='Agg'),
                              capture_output=True, text=True)

    def test_cli_rejected_calculations_return_failure_without_data(self):
        structure = "gridfactor = 1.0\nmaterial = [[10.0, 'GaAs', 0.0, 0.0, 0.0, 'n', 'w']]\n"
        cases = [('comp_scheme = 1\n' + structure, 'nonparabolicity'),
                 ('comp_scheme = 99\n' + structure, 'Unknown computation scheme'),
                 ("comp_scheme = 2\ngridfactor = 1.0\nmaxgridpoints = 1\n"
                  "material = [[10.0, 'GaAs', 0.0, 0.0, 0.0, 'n', 'w']]\n", 'maxgridpoints')]
        for text, message in cases:
            with self.subTest(message=message), TemporaryDirectory() as tmp:
                result = self.run_cli(tmp, text)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn(message, result.stderr)
                self.assertNotIn('Traceback', result.stderr)
                output = Path(tmp) / 'rejected_output'
                self.assertTrue((output / 'aestimo.log').is_file())
                self.assertEqual(list(output.glob('*.dat')), [])
                self.assertNotIn('Simulation is finished', (output / 'aestimo.log').read_text())

    def test_cli_malformed_inputs_are_usage_errors_without_outputs(self):
        for text, extension in [('material = [', '.py'), ('{', '.json'), ('[]', '.json')]:
            with self.subTest(extension=extension, text=text), TemporaryDirectory() as tmp:
                result = self.run_cli(tmp, text, extension)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn('Could not load input', result.stderr)
                self.assertNotIn('Traceback', result.stderr)
                self.assertFalse((Path(tmp) / 'rejected_output').exists())

    def test_gui_worker_reports_unsupported_solver_as_finished_failure(self):
        import aestimo_gui
        config = json.loads((ROOT / 'examples/sample_2qw_InGaAS_GaAs.json').read_text())
        config['solver'] = '1: Schrodinger + Non-parabolicity'
        saved = list(aestimo.logger.handlers)
        level = aestimo.logger.level
        propagate = aestimo.logger.propagate
        for handler in saved:
            aestimo.logger.removeHandler(handler)
        try:
            with TemporaryDirectory() as tmp, patch.object(aestimo, 'output_directory', tmp):
                worker = SimpleNamespace(examples_dir=tmp, project_name='rejected', sim_queue=queue.Queue())
                try:
                    aestimo_gui.AestimoGUI.run_simulation_worker(worker, config)
                    event, success, payload = worker.sim_queue.get_nowait()
                    self.assertEqual(event, 'finish')
                    self.assertFalse(success)
                    self.assertIn('nonparabolicity', payload[0])
                    self.assertIsNone(payload[1])
                    self.assertEqual(list((Path(tmp) / 'rejected_output').glob('*.dat')), [])
                finally:
                    for handler in list(aestimo.logger.handlers):
                        aestimo.logger.removeHandler(handler)
                        handler.close()
        finally:
            for handler in saved:
                aestimo.logger.addHandler(handler)
            aestimo.logger.setLevel(level)
            aestimo.logger.propagate = propagate
