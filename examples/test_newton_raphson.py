"""Residual-based Mode 10 acceptance, Jacobian and failed-result regressions."""

import sys
import runpy
from tempfile import TemporaryDirectory
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from aeslibs.newton_raphson import CoupledNewtonSolver, NewtonConvergenceError


def device_solver(nodes=7):
    """Small SI device with normalized equilibrium carriers and no net generation."""
    z = np.linspace(0, 1, nodes)
    n = np.exp(z)
    p = 1 / n
    fi = 0.2 * z
    model = SimpleNamespace(T=300.0, eps=np.full(nodes, 11.7 * 8.8541878128e-12),
                            ni_phys=np.full(nodes, 1e18), dd_max_iterations=25)
    solver = CoupledNewtonSolver(
        model, nodes, 1e-8, np.full(nodes, 1e18), np.zeros(nodes),
        np.ones(nodes), np.zeros(nodes), np.zeros(nodes), np.ones(nodes),
        np.ones(nodes), fi, n_stat=n, p_stat=p,
    )
    solver._init_equilibrium_state(n, p)
    solver.mun_mid = np.full(nodes - 1, 0.1)
    solver.mup_mid = np.full(nodes - 1, 0.05)
    solver.TAUN0 = solver.TAUP0 = np.full(nodes, 1e-7)
    solver.G_norm = np.zeros(nodes)
    return solver, fi, n, p


class TestNewtonSolver(unittest.TestCase):
    def test_equilibrium_needs_no_factorization_and_has_zero_current(self):
        solver, fi, n, p = device_solver()
        with patch('aeslibs.newton_raphson.spla.spsolve') as linear:
            *_, ok = solver.solve_step(fi, n, p, max_iter=0, tol=1e-10)
        self.assertTrue(ok)
        linear.assert_not_called()
        np.testing.assert_allclose(solver.Jtot, 0, atol=1e-10)
        self.assertEqual(solver.last_diagnostics['iterations'], 0)

    def test_analytic_jacobian_matches_central_differences(self):
        solver, fi, n, p = device_solver()
        fi = fi + np.array([0.0, 0.05, -0.03, 0.04, -0.02, 0.01, 0.0])
        n = n * np.linspace(0.95, 1.1, len(n))
        p = p * np.linspace(1.05, 0.9, len(p))
        state = np.column_stack((fi, n, p)).ravel()
        _, analytic = solver.compute_residual_and_jacobian(fi, n, p, Va=0.003)
        numeric = np.empty(analytic.shape)
        for col in range(len(state)):
            delta = 1e-5 * max(1.0, abs(state[col]))
            plus, minus = state.copy(), state.copy()
            plus[col] += delta
            minus[col] -= delta
            a, b = plus.reshape(-1, 3), minus.reshape(-1, 3)
            numeric[:, col] = (
                solver.compute_residual_only(*a.T, Va=0.003)
                - solver.compute_residual_only(*b.T, Va=0.003)
            ) / (2 * delta)
        np.testing.assert_allclose(analytic.toarray(), numeric, rtol=2e-6, atol=2e-7)

    def test_small_bias_converges_and_total_current_is_conserved(self):
        solver, fi, n, p = device_solver()
        fi, n, p, ok = solver.solve_step(fi, n, p, Va=0.001, max_iter=40, tol=1e-10)
        self.assertTrue(ok, solver.last_diagnostics)
        self.assertLess(np.max(np.abs(solver.compute_residual_only(fi, n, p, Va=0.001))), 1e-10)
        np.testing.assert_allclose(solver.Jtot, solver.Jtot[0], rtol=1e-6, atol=1e-8)

    def test_small_potential_update_cannot_hide_continuity_failure(self):
        solver, fi, n, p = device_solver()
        residual = np.zeros(3 * len(fi))
        residual[1::3] = 1.0
        with patch.object(solver, 'compute_residual_and_jacobian',
                          return_value=(residual, sp.eye(len(residual), format='csc'))):
            *_, ok = solver.solve_step(fi, n, p, max_iter=2)
        self.assertFalse(ok)
        self.assertEqual(solver.last_diagnostics['reason'], 'iteration limit')
        self.assertEqual(solver.last_diagnostics['electron_residual'], 1.0)
        self.assertIsNone(solver.Jtot)

    def test_final_allowed_update_is_checked(self):
        solver, fi, n, p = device_solver()
        initial = fi.copy()
        def linear_residual(phi, electrons, holes, Va=0):
            residual = np.zeros(3 * len(phi))
            residual[0::3] = phi - initial - 0.1
            return residual, sp.eye(len(residual), format='csc')
        with patch.object(solver, 'compute_residual_and_jacobian', side_effect=linear_residual):
            *_, ok = solver.solve_step(fi, n, p, max_iter=1, tol=1e-10)
        self.assertTrue(ok)
        self.assertEqual(solver.last_diagnostics['iterations'], 1)

    def test_nonfinite_state_residual_and_correction_are_rejected(self):
        for failure in ('state', 'residual', 'correction'):
            with self.subTest(failure=failure):
                solver, fi, n, p = device_solver()
                if failure == 'state':
                    fi[2] = np.nan
                    *_, ok = solver.solve_step(fi, n, p)
                else:
                    residual = np.ones(3 * len(fi))
                    if failure == 'residual':
                        residual[3] = np.inf
                    with patch.object(solver, 'compute_residual_and_jacobian',
                                      return_value=(residual, sp.eye(len(residual), format='csc'))), \
                         patch('aeslibs.newton_raphson.spla.spsolve',
                               return_value=np.full(len(residual), np.nan)):
                        *_, ok = solver.solve_step(fi, n, p)
                self.assertFalse(ok)
                self.assertIsNone(solver.Jtot)
                self.assertIn('non-finite', solver.last_diagnostics['reason'])

    def test_singular_linear_system_is_a_failed_step(self):
        solver, fi, n, p = device_solver()
        size = 3 * len(fi)
        with patch.object(solver, 'compute_residual_and_jacobian',
                          return_value=(np.ones(size), sp.csc_matrix((size, size)))):
            *_, ok = solver.solve_step(fi, n, p)
        self.assertFalse(ok)
        self.assertEqual(solver.last_diagnostics['reason'], 'linear solve failed')

    def test_failed_step_clears_previous_current_and_raises_with_voltage(self):
        solver, fi, n, p = device_solver()
        *_, ok = solver.solve_step(fi, n, p, max_iter=0)
        self.assertTrue(ok)
        self.assertIsNotNone(solver.Jtot)
        fi, n, p, ok = solver.solve_step(fi, n, p, Va=0.2, max_iter=0)
        self.assertFalse(ok)
        self.assertIsNone(solver.Jtot)
        self.assertTrue(np.isnan(solver.last_Jtot))
        with self.assertRaises(NewtonConvergenceError) as error:
            solver.require_convergence(ok, fi, n, p, 0.2)
        self.assertEqual(error.exception.voltage, 0.2)
        self.assertIn('iteration limit', str(error.exception))

    def test_solve_honors_model_iteration_budget(self):
        solver, fi, n, p = device_solver()
        solver.model.dd_max_iterations = 0
        fi, n, p, ok = solver.solve(fi, n, p, np.full(len(fi), 0.1),
                                   np.full(len(fi), 0.05), solver.TAUN0,
                                   solver.TAUP0, 0.0, 0.0, 0.0, Va=0.2)
        self.assertFalse(ok)
        self.assertEqual(solver.last_diagnostics['iterations'], 0)

    def test_guard_rejects_nonfinite_potential_even_with_success_flag(self):
        solver, fi, n, p = device_solver()
        fi[2] = np.nan
        with self.assertRaises(NewtonConvergenceError):
            solver.require_convergence(True, fi, n, p, 0.1)

    def test_solve_honors_model_residual_tolerance(self):
        solver, fi, n, p = device_solver()
        solver.model.dd_residual_tolerance = 1e-9
        with patch.object(solver, 'solve_step', return_value=(fi, n, p, True)) as step:
            solver.solve(fi, n, p, np.full(len(fi), 0.1),
                         np.full(len(fi), 0.05), solver.TAUN0,
                         solver.TAUP0, 0., 0., 0.)
        self.assertEqual(step.call_args.kwargs['tol'], 1e-9)
        solver.model.dd_residual_tolerance = float('nan')
        with self.assertRaises(ValueError):
            solver.solve(fi, n, p, np.full(len(fi), 0.1),
                         np.full(len(fi), 0.05), solver.TAUN0,
                         solver.TAUP0, 0., 0., 0.)

    def test_equilibrium_flux_is_exactly_zero_with_large_carrier_contrast(self):
        solver, fi, _, _ = device_solver()
        n = np.geomspace(1e-12, 1e12, len(fi))
        p = 1 / n
        solver.n_eq, solver.p_eq = n.copy(), p.copy()
        solver._init_equilibrium_state(n, p)
        solver.compute_currents(fi, n, p)
        np.testing.assert_array_equal(solver.Jtot, 0.)

    def test_tiny_potential_gradient_retains_small_flux(self):
        from aeslibs.newton_raphson import Ber
        solver, fi, _, _ = device_solver()
        solver.fi_eq = np.zeros_like(fi)
        n = np.full(len(fi), 1e8)
        p = np.full(len(fi), 1e-8)
        solver.n_eq, solver.p_eq = n.copy(), p.copy()
        solver._init_equilibrium_state(n, p)
        phi = np.arange(len(fi)) * 1e-17
        jn, jp = solver._edge_fluxes(phi, n, p)
        delta = np.diff(phi)
        np.testing.assert_allclose(jn, -solver.mun_mid * n[1:] * Ber(delta)
                                   * np.expm1(delta), rtol=1e-14, atol=0)
        self.assertTrue(np.all(jn != 0))
        self.assertTrue(np.all(jp != 0))

    def test_stable_flux_matches_direct_sg_away_from_equilibrium(self):
        from aeslibs.newton_raphson import Ber
        solver, fi, n, p = device_solver()
        phi = fi + np.arange(len(fi)) * .8
        n = n * np.linspace(.2, 2., len(n))
        p = p * np.linspace(2., .2, len(p))
        jn, jp = solver._edge_fluxes(phi, n, p)
        delta = np.diff(phi - solver.fi_eq)
        psi_n = delta + solver.d_psi_n0
        psi_p = delta + solver.d_psi_p0
        np.testing.assert_allclose(jn, solver.mun_mid *
            (n[1:] * Ber(psi_n) - n[:-1] * Ber(-psi_n)), rtol=1e-12)
        np.testing.assert_allclose(jp, solver.mup_mid *
            (p[:-1] * Ber(psi_p) - p[1:] * Ber(-psi_p)), rtol=1e-12)

    def test_residual_acceptance_cannot_bypass_current_conservation(self):
        solver, fi, n, p = device_solver()
        def bad_current(*args):
            solver.Jtot = np.arange(len(fi) - 1, dtype=float)
            solver.last_Jtot = 1.
        with patch.object(solver, 'compute_currents', side_effect=bad_current):
            *_, ok = solver.solve_step(fi, n, p, max_iter=0)
        self.assertFalse(ok)
        self.assertFalse(solver.last_diagnostics['current_conserved'])
        self.assertEqual(solver.last_diagnostics['reason'], 'current conservation not met')
        self.assertIsNone(solver.Jtot)
        with self.assertRaises(NewtonConvergenceError):
            solver.require_convergence(ok, fi, n, p, 0.)

    def test_invalid_current_conservation_controls_are_rejected(self):
        for name, values in [('dd_current_atol', [0., -1., np.nan, np.inf]),
                             ('dd_current_rtol', [-1., np.nan, np.inf])]:
            for value in values:
                solver, fi, n, p = device_solver()
                setattr(solver.model, name, value)
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    solver.solve_step(fi, n, p)

    def test_invalid_iteration_budget_and_tolerance_are_rejected(self):
        solver, fi, n, p = device_solver()
        for budget in (-1, 0.5, None):
            with self.subTest(budget=budget), self.assertRaises(ValueError):
                solver.solve_step(fi, n, p, max_iter=budget)
        for tolerance in (0., -1., np.nan, np.inf):
            with self.subTest(tolerance=tolerance), self.assertRaises(ValueError):
                solver.solve_step(fi, n, p, tol=tolerance)


class TestMode10Workflow(unittest.TestCase):
    """Real core routing with a small Si p-n structure and isolated output."""

    def setUp(self):
        import aestimo
        self.aestimo = aestimo
        temporary = TemporaryDirectory(prefix='aestimo-mode10-')
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)
        directory = patch.object(aestimo, 'output_directory', str(self.output))
        directory.start()
        self.addCleanup(directory.stop)
        preset = Path(__file__).resolve().parent / 'sample_pn.py'
        cfg = {k: v for k, v in runpy.run_path(str(preset)).items() if not k.startswith('__')}
        cfg.update(material=[[10., 'Si', 0., 0., 1e15, 'p', 'b'],
                             [10., 'Si', 0., 0., 1e15, 'n', 'b']],
                   gridfactor=1., vmin=0., vmax=0.001, Each_Step=0.001,
                   subnumber_h=1, subnumber_e=1, Quantum_Regions=False,
                   Quantum_Regions_boundary=np.zeros((1, 2)))
        self.input = SimpleNamespace(**cfg)

    def test_real_small_bias_sweep_exports_finite_accepted_steps(self):
        _, model, _, _ = self.aestimo.run_aestimo(self.input, drawFigures=False, show=False)
        data = np.loadtxt(self.output / 'av_curr.dat')
        np.testing.assert_allclose(data[:, 0], [0., 0.001])
        self.assertTrue(np.all(np.isfinite(data)))
        self.assertEqual(model.newton_solver.last_diagnostics['reason'], 'converged')
        self.assertLess(model.newton_solver.last_diagnostics['residual_norm'], 0.02)

    def test_failed_second_voltage_step_prevents_partial_sweep_export(self):
        original = CoupledNewtonSolver.solve_step
        def reject_second(solver, fi, n, p, **kwargs):
            if kwargs['Va'] == 0:
                return original(solver, fi, n, p, **kwargs)
            solver.last_diagnostics = {'reason': 'forced continuity failure',
                                       'residual_norm': 3., 'iterations': 25}
            return fi, n, p, False
        with patch.object(CoupledNewtonSolver, 'solve_step', reject_second), \
             patch.object(self.aestimo, 'save_and_plot2') as export:
            with self.assertRaises(NewtonConvergenceError) as error:
                self.aestimo.run_aestimo(self.input, drawFigures=False, show=False)
        self.assertEqual(error.exception.voltage, 0.001)
        export.assert_not_called()
        self.assertFalse((self.output / 'av_curr.dat').exists())

    def test_single_bias_is_solved_and_iteration_budget_cannot_bypass_guard(self):
        self.input.vmin = self.input.vmax = 0.1
        self.input.dd_max_iterations = 0
        with self.assertRaises(NewtonConvergenceError) as error:
            self.aestimo.run_aestimo(self.input, drawFigures=False, show=False)
        self.assertEqual(error.exception.voltage, 0.1)
        self.assertEqual(error.exception.diagnostics['iterations'], 0)
        self.assertFalse((self.output / 'av_curr.dat').exists())


if __name__ == '__main__':
    unittest.main()
