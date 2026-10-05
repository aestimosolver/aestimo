"""Reproducible dark Si pn mesh/tolerance audit; no experimental claim."""
import argparse
import json
import runpy
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import aestimo
from aeslibs.newton_raphson import CoupledNewtonSolver, NewtonConvergenceError


def run_case(spacing_nm, tolerance, max_voltage=0.4):
    cfg = {k: v for k, v in runpy.run_path(str(ROOT / 'examples/sample_pn.py')).items()
           if not k.startswith('__')}
    nodes = int(2000 / spacing_nm + 0.5)
    cfg.update(gridfactor=spacing_nm, dop_profile=np.zeros(nodes),
               vmax=max_voltage, dd_residual_tolerance=tolerance,
               dd_max_iterations=100)
    steps = []
    original = CoupledNewtonSolver.solve_step

    def record(solver, *args, **kwargs):
        result = original(solver, *args, **kwargs)
        step = dict(solver.last_diagnostics, voltage_V=float(kwargs['Va']),
                    accepted=bool(result[-1]))
        if result[-1]:
            current = solver.Jtot
            step.update(current_mA_cm2=float(solver.last_Jtot),
                        current_span_mA_cm2=float(np.ptp(current)))
        steps.append(step)
        return result

    failure = None
    with TemporaryDirectory(prefix='aestimo-mesh-') as temporary, \
         patch.object(aestimo, 'output_directory', temporary), \
         patch.object(CoupledNewtonSolver, 'solve_step', record):
        try:
            aestimo.run_aestimo(cfg, drawFigures=False, show=False)
        except NewtonConvergenceError as exc:
            failure = str(exc)
    return dict(spacing_nm=spacing_nm, nodes=nodes, residual_tolerance=tolerance,
                max_iterations=100, failure=failure, steps=steps)


def compare_cases(cases):
    """Use the finest mesh at each tolerance, matching voltages without extrapolation."""
    comparisons = []
    for tolerance in sorted({c['residual_tolerance'] for c in cases}):
        group = sorted((c for c in cases if c['residual_tolerance'] == tolerance),
                       key=lambda c: c['spacing_nm'])
        reference = {round(s['voltage_V'], 10): s for s in group[0]['steps']
                     if s['accepted']}
        for case in group[1:]:
            for step in case['steps']:
                key = round(step['voltage_V'], 10)
                if not step['accepted'] or key not in reference or key == 0:
                    continue
                fine = reference[key]['current_mA_cm2']
                relative = abs(step['current_mA_cm2'] - fine) / abs(fine) if fine else None
                comparisons.append(dict(residual_tolerance=tolerance,
                    spacing_nm=case['spacing_nm'], reference_spacing_nm=group[0]['spacing_nm'],
                    voltage_V=step['voltage_V'], relative_current_difference=relative))
    return comparisons


def run_audit():
    cases = []
    for tolerance in (0.02, 1e-7, 1e-8):
        for spacing in (10., 5., 2.5):
            print(f'Mesh audit: dx={spacing:g} nm, tolerance={tolerance:g}', flush=True)
            cases.append(run_case(spacing, tolerance))
    return dict(description='Numerical dark Si pn audit, not experimental validation',
                preset='examples/sample_pn.py', voltage_step_V=0.02,
                current_unit='mA/cm^2',
                cases=cases, comparisons=compare_cases(cases))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2, allow_nan=False) + '\n')
