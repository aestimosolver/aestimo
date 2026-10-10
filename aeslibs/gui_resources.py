"""Locate bundled GUI examples and prepare writable installed-app copies."""
import os
from pathlib import Path
import shutil


def bundled_examples_directory():
    """Locate checkout resources or the examples package in an installed wheel."""
    source = Path(__file__).resolve().parent.parent / 'examples'
    if (source / 'gaas_tobin1990_benchmark.json').is_file():
        return source
    import aestimo_examples
    return Path(aestimo_examples.__file__).resolve().parent


def prepare_gui_examples(destination=None):
    """Give the installed GUI writable presets without altering package files.

    Preserve existing user files. Source checkouts retain their existing example
    directory behavior; an explicit destination can be used in either environment.
    """
    bundled = bundled_examples_directory()
    source = Path(__file__).resolve().parent.parent / 'examples'
    if destination is None and bundled == source:
        return bundled
    if destination is None:
        root = Path(os.environ.get('AESTIMO_WORKSPACE', Path.home() / 'Aestimo')).expanduser()
        destination = root / 'examples'
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    for pattern in ('*.json', 'experimental_data/*.csv', 'experimental_data/README.md'):
        for item in sorted(bundled.glob(pattern)):
            target = destination / item.relative_to(bundled)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                # Exclusive creation prevents overwriting files even on a second launch.
                try:
                    with target.open('xb') as out, item.open('rb') as source_file:
                        shutil.copyfileobj(source_file, out)
                except FileExistsError:
                    pass
    return destination


def resolve_gui_reference(filename, examples_directory):
    """Honor explicit paths, then resolve bundled preset-relative reference paths."""
    path = Path(filename)
    if path.is_file() or path.is_absolute():
        return str(path)
    relative = path
    if relative.parts and relative.parts[0] == 'examples':
        relative = Path(*relative.parts[1:])
    candidate = Path(examples_directory) / relative
    if candidate.is_file():
        return str(candidate)
    candidate = Path(examples_directory) / 'experimental_data' / path.name
    return str(candidate) if candidate.is_file() else str(path)
