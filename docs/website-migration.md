# Original website and documentation migration

This user guide is based on the original
[Aestimo website repository](https://github.com/aestimosolver/aestimosolver.github.io),
with technical instructions checked against the 4.0.0 preparation source.

Source revision: `8274c65f188c3d630acada352d47f9a92fa5ada7` on `gh-pages`.
The old site remains unchanged by this migration.

| Original page | Destination / treatment |
| --- | --- |
| `index.md` | Overview, project citation, current installation links and project-history summary. Time-sensitive “latest release” labels are not copied into current guidance. |
| `about.md` | [About](about.md), [models/materials](models.md) and the historical [publication list](publications.md). |
| `authors.md` | [Contributor credits](authors.md), with original release-by-release source linked. Historical emails/affiliations are not presented as current contacts. |
| `documentation.md` | [Installation](installation.md), [quick start](quickstart.md), [inputs](inputs.md), [outputs](outputs.md) and [troubleshooting](troubleshooting.md). |
| `tutorials.md` | Four original notebooks and author credit in [Tutorials](tutorials.md). |
| `development.md` | Contribution guidance, project aims and historical reading list in [Development](development.md). |
| `releasenotes.md` | [Historical release notes](history.md), separate from current 4.0.0 preparation notes. |

## Corrections against the current code

- Use pip-created launch commands instead of modifying `site-packages`, setting
  executable bits or manually adding the source directory to `PATH`.
- Remove the separate `pylab` installation instruction; current dependencies are
  declared in `pyproject.toml`.
- Use a virtual environment or an explicitly configured conda environment.
  Python/Tk availability depends on the installed distribution.
- Explain that Python statements execute in order. Inputs can include arrays,
  lists and expressions, not only four scalar types.
- Replace the incomplete parameter placeholders with current field/unit tables
  and distinguish Python inputs from GUI JSON projects.
- Document current `*_output` directories and exporter names.
- Replace universal mesh/doping recommendations with convergence investigation
  that keeps the intended device model.
- Separate historical solver labels from routes accepted by the current code.
- Keep earlier version support and release dates in the historical section.

## Attribution and notices

Original manual notice: **Copyright (C) 2013–2020 Aestimo group**.
The site's manual declares reuse under the **GNU Free Documentation License 1.3**;
see the [original notice](https://github.com/aestimosolver/aestimosolver.github.io/blob/8274c65f188c3d630acada352d47f9a92fa5ada7/documentation.md)
and the [license text](https://www.gnu.org/licenses/fdl-1.3.html).
The simulator code remains GNU GPL v3 or later, as stated in its
[license](https://github.com/aestimosolver/aestimo/blob/master/COPYING.md).

Project and contributor attribution is retained. Historical bibliography and
release records describe their original periods; current instructions have been
rewritten for the present package and checked against its source.
