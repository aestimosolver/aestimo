"""Sphinx configuration without importing the solver or initializing Tk."""
from pathlib import Path
import tomllib

project = "Aestimo 1D"
author = "Sefer Bora Lisesivdin and Aestimo contributors"
copyright = "2026, Aestimo contributors"
metadata = tomllib.loads(
    (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(encoding="utf-8")
)
release = metadata["project"]["version"]
version = release
extensions = ["myst_parser"]
myst_enable_extensions = ["dollarmath"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
html_title = "Aestimo documentation"
