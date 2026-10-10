"""Numerical agreement must not invent calibration or experimental provenance."""

import json
import sys
import unittest
from pathlib import Path

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from aeslibs.qw_validation import (
    QWTraceabilityRecord, compute_qw_error_metrics, assess_qw_comparison,
    curve_comparison_metrics, generate_qw_validation_report_markdown,
    plot_standardized_qw_validation_suite,
)
from aeslibs.validation_policy import (
    UNVERIFIED_REFERENCE, SYNTHETIC_REFERENCE, MODEL_BASED,
    reference_status, reviewed_status,
)
from examples.audit_reference_data import inventory


def record():
    return QWTraceabilityRecord(
        benchmark_name='Synthetic comparison fixture', paper_title='Fixture',
        authors='Test', journal='None', year=2026, doi='',
        material_system='GaAs', well_material='GaAs', barrier_material='AlGaAs',
        well_width_nm=10., barrier_width_nm=20., temperature_k=300.,
    )


class TestValidationPolicy(unittest.TestCase):
    def test_offset_with_perfect_correlation_fails_agreement(self):
        ref = np.array([1., 2., 3.])
        metrics = compute_qw_error_metrics(ref, ref + 1., 'energy_ev')
        self.assertAlmostEqual(metrics['pearson_r2'], 1.)
        self.assertAlmostEqual(metrics['r_squared'], -0.5)
        self.assertEqual(assess_qw_comparison(metrics), 'METRICS FAILED')
        report = generate_qw_validation_report_markdown(record(), [metrics])
        self.assertIn('METRICS FAILED', report)
        self.assertNotIn('`CALIBRATED`', report)
        self.assertNotIn('`EXPERIMENTALLY VALIDATED`', report)

    def test_sign_reversal_with_perfect_squared_correlation_fails(self):
        ref = np.array([1., 2., 3.])
        metrics = compute_qw_error_metrics(ref, -ref)
        self.assertAlmostEqual(metrics['pearson_r2'], 1.)
        self.assertEqual(assess_qw_comparison(metrics), 'METRICS FAILED')

    def test_exact_variable_match_passes_metrics_without_verifying_source(self):
        ref = np.array([1., 2., 3.])
        metrics = compute_qw_error_metrics(ref, ref)
        self.assertEqual(assess_qw_comparison(metrics), 'METRICS PASSED')
        self.assertEqual(record().to_dict()['provenance_classification'], UNVERIFIED_REFERENCE)

    def test_invalid_constant_and_unpaired_data_are_not_assessable(self):
        for ref, sim in [([], []), ([1.], [1.]), ([1., 1.], [1., 1.]),
                         ([1., 2.], [1.]), ([1., np.nan], [1., 2.]),
                         ([[1., 2.]], [[1., 2.]])]:
            with self.subTest(reference=ref):
                self.assertEqual(assess_qw_comparison(compute_qw_error_metrics(ref, sim)),
                                 'NOT ASSESSABLE')
        report = generate_qw_validation_report_markdown(record(), [])
        self.assertIn('NOT ASSESSED', report)

    def test_curve_comparison_does_not_extrapolate(self):
        metrics = curve_comparison_metrics([0., 1., 2., 3.], [9., 1., 2., 9.],
                                           [1., 2.], [1., 2.], 'energy_ev')
        self.assertEqual(metrics['n_points'], 2)
        self.assertEqual(metrics['rmse'], 0.)
        no_overlap = curve_comparison_metrics([0., 1.], [1., 2.],
                                              [2., 3.], [1., 2.], 'energy_ev')
        self.assertEqual(assess_qw_comparison(no_overlap), 'NOT ASSESSABLE')

    def test_plot_recomputes_statistics_instead_of_trusting_supplied_scores(self):
        data = {'dingle_data': {'exp_lw': [1., 2., 3.], 'exp_e1_hh1': [1., 2., 3.],
                                'sim_lw': [1., 2., 3.], 'sim_e1_hh1': [2., 3., 4.],
                                'r_squared': 0.998, 'rmse_mev': 2.1}}
        fig = plot_standardized_qw_validation_suite(data, record())
        self.addCleanup(plt.close, fig)
        text = '\n'.join(item.get_text() for item in fig.axes[1].texts)
        self.assertIn('-0.5000', text)
        self.assertIn('1000.00', text)
        self.assertNotIn('0.998', text)

    def test_synthetic_and_claimed_reference_statuses_stay_distinct(self):
        self.assertEqual(reference_status(['si_pn_experimental_iv.csv']), SYNTHETIC_REFERENCE)
        self.assertEqual(reference_status(['ingaas_pn_experimental_iv.csv']), SYNTHETIC_REFERENCE)
        self.assertEqual(reference_status(['qw_miller1984_qcse_stark_shift.csv']), UNVERIFIED_REFERENCE)
        self.assertEqual(reference_status([]), MODEL_BASED)
        self.assertEqual(reviewed_status('EXPERIMENTALLY VALIDATED'), UNVERIFIED_REFERENCE)
        self.assertIn('CALIBRATED MODEL', reviewed_status('FITTED TO EXPERIMENT'))

    def test_reference_inventory_matches_current_files_and_has_no_verified_claims(self):
        actual = inventory(ROOT / 'examples' / 'experimental_data')
        saved = json.loads((ROOT / 'docs' / 'reference-data-audit.json').read_text())
        self.assertEqual(saved, actual)
        self.assertEqual(len(actual), len(list((ROOT / 'examples' / 'experimental_data').glob('*.csv'))))
        self.assertEqual(sum(r['reference_kind'] == 'SYNTHETIC_OR_MODEL_REFERENCE' for r in actual), 3)
        self.assertTrue(all(r['original_source_review'] == 'NOT PERFORMED' for r in actual))


if __name__ == '__main__':
    unittest.main()
