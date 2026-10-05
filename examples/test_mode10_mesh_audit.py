"""Audit comparisons must not fabricate overlap or divide by zero."""
import unittest
from mode10_mesh_audit import compare_cases


class TestMeshComparison(unittest.TestCase):
    def test_only_matching_accepted_nonzero_bias_is_compared(self):
        fine = dict(residual_tolerance=1e-7, spacing_nm=2.5,
                    steps=[dict(voltage_V=0., accepted=True, current_mA_cm2=0.),
                           dict(voltage_V=.2, accepted=True, current_mA_cm2=2.)])
        coarse = dict(residual_tolerance=1e-7, spacing_nm=5.,
                      steps=[dict(voltage_V=0., accepted=True, current_mA_cm2=0.),
                             dict(voltage_V=.1, accepted=True, current_mA_cm2=1.),
                             dict(voltage_V=.2, accepted=True, current_mA_cm2=3.),
                             dict(voltage_V=.4, accepted=False)])
        result = compare_cases([coarse, fine])
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]['relative_current_difference'], .5)
        self.assertEqual(result[0]['reference_spacing_nm'], 2.5)

    def test_failed_finest_mesh_cannot_be_replaced_by_coarse_reference(self):
        cases = [dict(residual_tolerance=1e-8, spacing_nm=spacing,
                      steps=[dict(voltage_V=.2, accepted=accepted, current_mA_cm2=1.)])
                 for spacing, accepted in [(10., True), (5., True), (2.5, False)]]
        self.assertEqual(compare_cases(cases), [])
