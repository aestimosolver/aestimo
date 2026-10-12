"""Check an installed Aestimo distribution outside its source checkout.

Run with the Python interpreter of the fresh environment containing the wheel.
No live GUI window is opened. All calculation outputs use a temporary directory.
"""
import argparse
from importlib import metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
from tempfile import TemporaryDirectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True,
                        help='Source checkout used to build the installed wheel')
    parser.add_argument('--expected-version', default='4.0.0')
    args = parser.parse_args()
    source = args.source.resolve()
    example = source / 'examples/sample_1qw_barrierdope_ingaas.py'
    if not example.is_file():
        parser.error(f'Missing CLI example: {example}')
    distribution = metadata.distribution('aestimo')
    assert distribution.version == args.expected_version, distribution.version
    command = Path(sysconfig.get_path('scripts')) / (
        'aestimo.exe' if os.name == 'nt' else 'aestimo'
    )
    assert command.is_file(), command
    env = dict(os.environ, MPLBACKEND='Agg')
    previous = Path.cwd()
    with TemporaryDirectory(prefix='aestimo-release-') as temporary:
        work = Path(temporary)
        os.chdir(work)
        try:
            import aestimo
            import aestimo_gui
            from aeslibs.gui_resources import prepare_gui_examples, resolve_gui_reference

            assert not list(work.iterdir()), 'Import created working-directory files'
            for module in (aestimo, aestimo_gui):
                assert not Path(module.__file__).resolve().is_relative_to(source), module.__file__
            assert aestimo.__version__ == distribution.version
            entry_points = {ep.name: ep.value for ep in distribution.entry_points}
            assert entry_points['aestimo'] == 'aestimo:main'
            assert entry_points['aestimo-gui'] == 'aestimo_gui:main'
            assert entry_points['aestimo-gui-win'] == 'aestimo_gui:main'
            for option in ('--help', '--version', '-v'):
                result = subprocess.run([str(command), option], cwd=work, env=env,
                                        text=True, capture_output=True, check=True)
                if option != '--help':
                    assert result.stdout.strip() == f'Aestimo 1D {distribution.version}'

            for arguments in ([], ['-i', 'missing.py']):
                result = subprocess.run([str(command), *arguments], cwd=work, env=env,
                                        text=True, capture_output=True)
                assert result.returncode == 2, result.stderr
            assert not list(work.iterdir()), 'CLI startup created output files'

            examples = prepare_gui_examples(work / 'examples')
            for pattern in ('*.json', 'experimental_data/*.csv'):
                originals = list((source / 'examples').glob(pattern))
                assert originals, pattern
                for original in originals:
                    copied = examples / original.relative_to(source / 'examples')
                    assert copied.read_bytes() == original.read_bytes(), original
            for preset in examples.glob('*.json'):
                data = json.loads(preset.read_text(encoding='utf-8'))
                if isinstance(data, dict) and data.get('exp_file'):
                    assert Path(resolve_gui_reference(data['exp_file'], examples)).is_file(), preset
            user_project = examples / 'gaas_tobin1990_benchmark.json'
            user_project.write_text('{"user_edited": true}\n', encoding='utf-8')
            prepare_gui_examples(examples)
            assert user_project.read_text(encoding='utf-8') == '{"user_edited": true}\n'

            result = subprocess.run([str(command), '-i', str(example)], cwd=work,
                                    env=env, text=True, capture_output=True)
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            output = work / f'{example.stem}_output'
            assert (output / 'aestimo.log').is_file(), output
            assert (output / 'aestimo.log').read_text().count('Simulation is finished') == 1
            assert (output / 'potn_eh_equi_cond.dat').is_file(), output
        finally:
            os.chdir(previous)
    print(f'Aestimo {distribution.version}: installed CLI calculation, entry points, '
          'GUI imports and bundled resources passed.')


if __name__ == '__main__':
    main()
