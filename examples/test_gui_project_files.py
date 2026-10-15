"""Check GUI file handling without creating a Tk window."""
import builtins
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aestimo_gui import AestimoGUI


class GUIProjectFilesTests(unittest.TestCase):
    def test_default_project_preserves_existing_user_file(self):
        with TemporaryDirectory() as directory:
            project = Path(directory) / 'untitled_project.json'
            original = b'{"user_edited": true}\n'
            project.write_bytes(original)
            app = SimpleNamespace(examples_dir=directory,
                                  get_current_configuration=Mock())
            AestimoGUI.auto_save_default_project(app)
            self.assertEqual(project.read_bytes(), original)
            app.get_current_configuration.assert_not_called()

    def test_default_project_creates_missing_file(self):
        with TemporaryDirectory() as directory:
            config = {'description': 'InGaN quantum well'}
            app = SimpleNamespace(examples_dir=str(Path(directory) / 'examples'),
                                  get_current_configuration=lambda: config)
            AestimoGUI.auto_save_default_project(app)
            project = Path(app.examples_dir) / 'untitled_project.json'
            self.assertEqual(json.loads(project.read_text(encoding='utf-8')), config)

    def test_default_project_handles_concurrent_creation(self):
        with TemporaryDirectory() as directory:
            project = Path(directory) / 'untitled_project.json'
            def configuration():
                project.write_text('{"other_launch": true}', encoding='utf-8')
                return {'replacement': True}
            app = SimpleNamespace(examples_dir=directory,
                                  get_current_configuration=configuration)
            AestimoGUI.auto_save_default_project(app)
            self.assertEqual(json.loads(project.read_text()), {'other_launch': True})

    def test_utf8_project_load_is_independent_of_system_encoding(self):
        with TemporaryDirectory() as directory:
            project = Path(directory) / 'quantum_well.json'
            config = {'description': 'Türkçe – ψ electron state', 'material': []}
            project.write_text(json.dumps(config, ensure_ascii=False), encoding='utf-8')
            app = SimpleNamespace(load_configuration=Mock(), status_label=Mock(),
                                  load_previous_results_for_project=Mock(return_value=None))
            real_open = builtins.open
            def locale_open(file, mode='r', *args, **kwargs):
                if 'b' not in mode and 'encoding' not in kwargs:
                    kwargs['encoding'] = 'cp1252'
                return real_open(file, mode, *args, **kwargs)
            with patch('builtins.open', side_effect=locale_open), \
                    patch('tkinter.messagebox.showerror') as error:
                AestimoGUI.load_project(app, str(project))
            app.load_configuration.assert_called_once_with(config)
            error.assert_not_called()


if __name__ == '__main__':
    unittest.main()
