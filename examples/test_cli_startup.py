"""Release regressions for CLI startup, plot selection and calculation logging."""
import logging
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import aestimo


class TestCLIStartup(unittest.TestCase):
    def setUp(self):
        self.handlers = list(aestimo.logger.handlers)
        self.level = aestimo.logger.level
        self.propagate = aestimo.logger.propagate
        for handler in self.handlers:
            aestimo.logger.removeHandler(handler)

    def tearDown(self):
        for handler in list(aestimo.logger.handlers):
            aestimo.logger.removeHandler(handler)
            handler.close()
        for handler in self.handlers:
            aestimo.logger.addHandler(handler)
        aestimo.logger.setLevel(self.level)
        aestimo.logger.propagate = self.propagate

    def test_import_help_and_version_do_not_create_output_files(self):
        env = dict(os.environ, PYTHONPATH=str(ROOT), MPLBACKEND='Agg')
        commands = [
            ['-c', 'import aestimo; import aestimo_gui'],
            ['-m', 'aestimo', '--help'],
            ['-m', 'aestimo', '--version'],
        ]
        for command in commands:
            with self.subTest(command=command), TemporaryDirectory() as tmp:
                result = subprocess.run([sys.executable, *command], cwd=tmp,
                                        env=env, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_missing_input_is_an_error_without_output_files(self):
        env = dict(os.environ, PYTHONPATH=str(ROOT), MPLBACKEND='Agg')
        for args in ([], ['-i', 'missing.py']):
            with self.subTest(args=args), TemporaryDirectory() as tmp:
                result = subprocess.run([sys.executable, '-m', 'aestimo', *args],
                                        cwd=tmp, env=env, capture_output=True, text=True)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_draw_flag_is_passed_to_the_solver_on_each_invocation(self):
        with TemporaryDirectory() as tmp, patch.object(aestimo, 'output_directory', 'output'):
            source = Path(tmp) / 'cli_startup_input.py'
            source.write_text('T = 300.0\n', encoding='utf-8')
            with patch.object(aestimo, 'run_aestimo') as run:
                self.assertEqual(aestimo.main(['-d', '-i', str(source)]), 0)
                self.assertTrue(run.call_args.kwargs['drawFigures'])
                self.assertEqual(aestimo.main(['-i', str(source)]), 0)
                self.assertFalse(run.call_args.kwargs['drawFigures'])

    def test_logger_retargets_and_preserves_application_handlers(self):
        with TemporaryDirectory() as tmp, patch.object(aestimo, 'output_directory', tmp):
            foreign = logging.NullHandler()
            aestimo.logger.addHandler(foreign)
            aestimo.initialize_logger()
            aestimo.initialize_logger()
            aestimo.logger.info('first-run-marker')
            first = Path(tmp) / 'aestimo.log'
            self.assertEqual(first.read_text().count('first-run-marker'), 1)
            aestimo.output_directory = str(Path(tmp) / 'second')
            aestimo.initialize_logger()
            aestimo.logger.info('second-run-marker')
            self.assertNotIn('second-run-marker', first.read_text())
            self.assertEqual((Path(tmp) / 'second/aestimo.log').read_text().count('second-run-marker'), 1)
            self.assertIn(foreign, aestimo.logger.handlers)
