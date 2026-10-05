"""GUI preset resources remain available and writable after wheel installation."""
import os
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from aeslibs.gui_resources import (
    bundled_examples_directory, prepare_gui_examples, resolve_gui_reference,
)


class TestGUIResources(unittest.TestCase):
    def test_checkout_presets_and_references_are_located(self):
        directory = bundled_examples_directory()
        self.assertTrue((directory / 'gaas_tobin1990_benchmark.json').is_file())
        self.assertTrue((directory / 'experimental_data/qw_dingle1975_energy_vs_width.csv').is_file())
        self.assertEqual(prepare_gui_examples(), directory)

    def test_copy_preserves_existing_user_files_and_reference_bytes(self):
        bundled = bundled_examples_directory()
        with TemporaryDirectory() as temporary:
            target = Path(temporary) / 'examples'
            prepare_gui_examples(target)
            for source in bundled.glob('experimental_data/*.csv'):
                self.assertEqual((target / 'experimental_data' / source.name).read_bytes(),
                                 source.read_bytes())
            project = target / 'gaas_tobin1990_benchmark.json'
            project.write_text('{"user_edited": true}\n')
            prepare_gui_examples(target)
            self.assertEqual(project.read_text(), '{"user_edited": true}\n')
            self.assertNotEqual(project.read_bytes(), (bundled / project.name).read_bytes())
            self.assertFalse((target / 'sample_pn.py').exists())

    def test_preset_reference_paths_resolve_outside_checkout(self):
        with TemporaryDirectory() as temporary:
            target = prepare_gui_examples(Path(temporary) / 'examples')
            name = 'qw_dingle1975_energy_vs_width.csv'
            expected = str(target / 'experimental_data' / name)
            previous = Path.cwd()
            try:
                os.chdir(temporary)
                for reference in (name, 'experimental_data/' + name,
                                  'examples/experimental_data/' + name):
                    self.assertEqual(str(Path(resolve_gui_reference(reference, target)).resolve()), expected)
                self.assertEqual(resolve_gui_reference('missing.csv', target), 'missing.csv')
            finally:
                os.chdir(previous)

    def test_explicit_absolute_reference_is_not_replaced_by_same_named_resource(self):
        with TemporaryDirectory() as temporary:
            target = prepare_gui_examples(Path(temporary) / 'examples')
            path = Path(temporary) / 'custom/qw_dingle1975_energy_vs_width.csv'
            self.assertEqual(resolve_gui_reference(path, target), str(path))
            path.parent.mkdir()
            path.write_text('custom data')
            self.assertEqual(resolve_gui_reference(path, target), str(path))
