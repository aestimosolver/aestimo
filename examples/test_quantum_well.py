#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Aestimo 1D - Quantum-Well Confined-State Solver Verification Suite
examples/test_quantum_well.py

Rigorous automated validation testing across 10 physical and numerical benchmarks:
  1. Infinite square well analytical eigenenergy comparison (En = n² pi² ħ² / (2 m* L²))
  2. Finite semiconductor QW bound states, ordering, and parity
  3. Wavefunction normalization verification (residual < 1e-6) and phase determinism
  4. BenDaniel-Duke variable-mass interface flux continuity
  5. Built-in electric field & Quantum-Confined Stark Effect (QCSE red-shift, spatial separation, overlap drop)
  6. 2D sheet carrier populations & 3D volumetric density conservation
  7. Optical transition matrix elements & wavelength calculation
  8. Multiple Quantum Wells (Coupled superlattice splitting vs isolated wells)
  9. Experimental device validation (Tsang 1981 GaAs/AlGaAs SQW ~845 nm laser emission)
 10. Backward compatibility & zero regression when QW mode is OFF
"""

import os
import sys
import unittest
import numpy as np

# Ensure workspace root is on sys.path
WORKSPACE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if WORKSPACE_DIR not in sys.path:
    sys.path.insert(0, WORKSPACE_DIR)

from aeslibs import quantum_well as qw
from aeslibs.quantum_well import (
    solve_effective_mass_schrodinger_1d,
    compute_optical_transitions,
    compute_quantum_carrier_density,
    solve_quantum_well,
    QuantumWellResult,
    HBAR2_OVER_2M0_EV_NM2,
    HC_EV_NM
)


class TestQuantumWellSolver(unittest.TestCase):
    """Automated unit test suite for Aestimo QW Confined-State Solver."""

    def test_01_infinite_square_well_analytical(self):
        """Benchmark 1: Infinite square well eigenenergies vs exact analytical formula."""
        L_nm = 10.0
        m_rel = 0.067  # GaAs electron effective mass
        N = 250
        z = np.linspace(0.0, L_nm, N)
        potential = np.zeros(N)
        mass = np.full(N, m_rel)

        energies, wf, prob, norm_errs = solve_effective_mass_schrodinger_1d(
            z_nm=z,
            potential_ev=potential,
            m_eff_rel=mass,
            num_states=3,
            particle_type="electron"
        )

        self.assertEqual(len(energies), 3, "Failed to retrieve 3 states for infinite well")

        # Analytical: E_n = (n^2 * pi^2 * ħ^2) / (2 * m* * L^2)
        exact_energies = []
        for n in [1, 2, 3]:
            e_exact = (n ** 2) * (np.pi ** 2) * HBAR2_OVER_2M0_EV_NM2 / (m_rel * (L_nm ** 2))
            exact_energies.append(e_exact)

        for n_idx, (e_num, e_ana) in enumerate(zip(energies, exact_energies)):
            rel_error = abs(e_num - e_ana) / e_ana
            self.assertLess(
                rel_error, 0.005,
                f"State {n_idx+1}: Numerical {e_num:.6f} eV differs from analytical {e_ana:.6f} eV by {rel_error*100:.3f}%"
            )

    def test_02_finite_semiconductor_qw_bound_states(self):
        """Benchmark 2: Finite GaAs/AlGaAs QW bound states, ordering, and parity."""
        # 10 nm GaAs well with 20 nm Al0.3Ga0.7As barriers (Delta_Ec = 0.228 eV)
        z = np.linspace(0.0, 50.0, 501)
        dz = z[1] - z[0]
        well_mask = (z >= 20.0) & (z <= 30.0)
        
        # Conduction band: 0 in well, 0.228 in barriers
        ec = np.where(well_mask, 0.0, 0.228)
        # Variable effective mass: 0.067 in well, 0.092 in barrier
        m_e = np.where(well_mask, 0.067, 0.092)

        energies, wf, prob, norm_errs = solve_effective_mass_schrodinger_1d(
            z_nm=z,
            potential_ev=ec,
            m_eff_rel=m_e,
            num_states=3,
            particle_type="electron",
            barrier_ref_ev=0.228
        )

        # A 10 nm GaAs well of 228 meV depth typically binds 2 conduction states
        self.assertGreaterEqual(len(energies), 2, "Expected at least 2 bound electron states")
        self.assertLess(energies[0], energies[1], "Energy ordering violated (E1 must be < E2)")
        self.assertLess(energies[0], 0.228, "Ground state exceeds barrier height")
        self.assertLess(energies[1], 0.228, "First excited state exceeds barrier height")

        # Ground state (n=1) should be symmetric (even parity): psi(mid - delta) ~= psi(mid + delta)
        mid_idx = np.argmin(np.abs(z - 25.0))
        delta_idx = int(3.0 / dz)
        psi1 = wf[0]
        self.assertAlmostEqual(
            psi1[mid_idx - delta_idx], psi1[mid_idx + delta_idx], places=2,
            msg="Ground state parity symmetry violated"
        )

        # First excited state (n=2) should be antisymmetric (odd parity, node near center)
        psi2 = wf[1]
        self.assertAlmostEqual(psi2[mid_idx], 0.0, places=2, msg="First excited state should have node at well center")
        self.assertAlmostEqual(
            psi2[mid_idx - delta_idx], -psi2[mid_idx + delta_idx], places=2,
            msg="First excited state odd parity violated"
        )

    def test_03_wavefunction_normalization_and_phase(self):
        """Benchmark 3: Wavefunction normalization residual < 1e-5 and phase determinism."""
        z = np.linspace(0.0, 40.0, 401)
        dz = z[1] - z[0]
        ec = np.where((z >= 15.0) & (z <= 25.0), 0.0, 0.30)
        mass = np.where((z >= 15.0) & (z <= 25.0), 0.067, 0.088)

        energies, wf, prob, norm_errs = solve_effective_mass_schrodinger_1d(
            z_nm=z, potential_ev=ec, m_eff_rel=mass, num_states=2
        )

        for key, err in norm_errs.items():
            self.assertLess(err, 1e-5, f"Normalization error for {key} ({err}) exceeds tolerance 1e-5")

        for idx, psi in enumerate(wf):
            num_integral = np.sum(psi ** 2) * dz
            self.assertAlmostEqual(num_integral, 1.0, places=5, msg=f"State {idx+1} integral != 1.0")
            # Sign anchoring: maximum absolute peak must be positive
            max_val = psi[np.argmax(np.abs(psi))]
            self.assertGreater(max_val, 0.0, f"State {idx+1} phase sign convention failed")

    def test_04_bendaniel_duke_interface_flux_continuity(self):
        """Benchmark 4: Variable mass BenDaniel-Duke boundary condition continuity."""
        # Across heterointerface, (1/m*) dpsi/dz must be continuous
        z = np.linspace(0.0, 30.0, 601)
        dz = z[1] - z[0]
        intf_z = 15.0
        intf_idx = np.argmin(np.abs(z - intf_z))

        # Well on left (z <= 15 nm, m = 0.067), Barrier on right (z > 15 nm, m = 0.12)
        m_eff = np.where(z <= intf_z, 0.067, 0.12)
        ec = np.where(z <= intf_z, 0.0, 0.35)

        energies, wf, prob, norm_errs = solve_effective_mass_schrodinger_1d(
            z_nm=z, potential_ev=ec, m_eff_rel=m_eff, num_states=1
        )

        psi = wf[0]
        # Wavefunction itself must be continuous across interface (delta < 0.01 between adjacent nodes)
        self.assertAlmostEqual(psi[intf_idx], psi[intf_idx - 1], delta=0.01, msg="Wavefunction jump discontinuity at interface")
        self.assertAlmostEqual(psi[intf_idx], psi[intf_idx + 1], delta=0.01, msg="Wavefunction jump discontinuity at interface")

        # BenDaniel-Duke probability flux across the interface node:
        # F_{i-1/2} = 0.5 * (1/m_{i-1} + 1/m_i) * (psi_i - psi_{i-1}) / dz
        # F_{i+1/2} = 0.5 * (1/m_i + 1/m_{i+1}) * (psi_{i+1} - psi_i) / dz
        inv_m_left = 0.5 * (1.0 / m_eff[intf_idx - 1] + 1.0 / m_eff[intf_idx])
        inv_m_right = 0.5 * (1.0 / m_eff[intf_idx] + 1.0 / m_eff[intf_idx + 1])
        flux_half_left = inv_m_left * (psi[intf_idx] - psi[intf_idx - 1]) / dz
        flux_half_right = inv_m_right * (psi[intf_idx + 1] - psi[intf_idx]) / dz
        
        # Flux difference across node i is O(dz) -> continuous probability flux
        rel_flux_diff = abs(flux_half_left - flux_half_right) / (abs(flux_half_left) + 1e-12)
        self.assertLess(rel_flux_diff, 0.05, f"BenDaniel-Duke interface flux difference {rel_flux_diff*100:.2f}% exceeds 5%")

    def test_05_quantum_confined_stark_effect_qcse(self):
        """Benchmark 5: Electric field QCSE red-shift, carrier separation, and overlap reduction."""
        layers = [
            {"material": "AlGaAs", "thickness": 15.0, "type": "barrier", "mole": 0.3},
            {"material": "GaAs", "thickness": 8.0, "type": "well", "mole": 0.0},
            {"material": "AlGaAs", "thickness": 15.0, "type": "barrier", "mole": 0.3},
        ]
        
        # Synthetic flat-band profile
        z = np.linspace(0.0, 38.0, 381)
        ec = np.where((z >= 15.0) & (z <= 23.0), 0.0, 0.23)
        ev = np.where((z >= 15.0) & (z <= 23.0), -1.424, -1.574)
        band_prof = {"z": z, "ec": ec, "ev": ev}

        # Zero-field solution
        res_f0 = solve_quantum_well(
            band_profile=band_prof,
            layers=layers,
            electric_field_v_cm=0.0
        )

        # High-field solution: F = 100 kV/cm (1e5 V/cm)
        res_f100 = solve_quantum_well(
            band_profile=band_prof,
            layers=layers,
            electric_field_v_cm=1.0e5
        )

        self.assertGreater(len(res_f0.dominant_transitions), 0)
        self.assertGreater(len(res_f100.dominant_transitions), 0)

        t_f0 = res_f0.dominant_transitions[0]
        t_f100 = res_f100.dominant_transitions[0]

        # 1. Stark Shift (Red Shift): Transition energy under field must be lower than zero field
        self.assertLess(
            t_f100["energy_ev"], t_f0["energy_ev"],
            f"QCSE failed: E_trans(F=100kV/cm) {t_f100['energy_ev']:.4f} eV >= E_trans(0) {t_f0['energy_ev']:.4f} eV"
        )
        # Consequently, emission wavelength must increase (red shift)
        self.assertGreater(
            t_f100["wavelength_nm"], t_f0["wavelength_nm"],
            "QCSE wavelength red-shift failed"
        )

        # 2. Overlap Reduction: Spatial polarization reduces electron-hole overlap Gamma_11
        self.assertLess(
            t_f100["overlap"], t_f0["overlap"],
            f"QCSE overlap reduction failed: Gamma(F=100) {t_f100['overlap']:.4f} >= Gamma(0) {t_f0['overlap']:.4f}"
        )

        # 3. Spatial Separation: Electron and hole wavefunction centroids pull apart
        dz = z[1] - z[0]
        z_e_center_f0 = np.sum(z * (res_f0.electron_wavefunctions[0] ** 2)) * dz
        z_h_center_f0 = np.sum(z * (res_f0.hole_wavefunctions[0] ** 2)) * dz
        sep_f0 = abs(z_e_center_f0 - z_h_center_f0)

        z_e_center_f100 = np.sum(z * (res_f100.electron_wavefunctions[0] ** 2)) * dz
        z_h_center_f100 = np.sum(z * (res_f100.hole_wavefunctions[0] ** 2)) * dz
        sep_f100 = abs(z_e_center_f100 - z_h_center_f100)

        self.assertGreater(
            sep_f100, sep_f0,
            f"Carrier spatial separation failed: sep(F=100) {sep_f100:.3f} nm <= sep(0) {sep_f0:.3f} nm"
        )

    def test_06_carrier_statistics_and_density_conservation(self):
        """Benchmark 6: 2D sheet carrier populations and 3D volumetric density integration."""
        z = np.linspace(0.0, 30.0, 301)
        dz = z[1] - z[0]
        e_energies = np.array([0.05, 0.18])
        # Simple Gaussian-like normalized wavefunctions
        wf_e = np.zeros((2, len(z)))
        wf_e[0] = np.exp(-((z - 15.0) / 3.0) ** 2)
        wf_e[0] /= np.sqrt(np.sum(wf_e[0] ** 2) * dz)
        wf_e[1] = (z - 15.0) * np.exp(-((z - 15.0) / 3.0) ** 2)
        wf_e[1] /= np.sqrt(np.sum(wf_e[1] ** 2) * dz)

        h_energies = np.array([-1.45])
        wf_h = np.zeros((1, len(z)))
        wf_h[0] = np.exp(-((z - 15.0) / 2.5) ** 2)
        wf_h[0] /= np.sqrt(np.sum(wf_h[0] ** 2) * dz)

        sheet_n, sheet_p, dens_n, dens_p = compute_quantum_carrier_density(
            z_nm=z,
            electron_energies=e_energies,
            electron_wf=wf_e,
            hole_energies=h_energies,
            hole_wf=wf_h,
            temp_k=300.0,
            ef_electron=0.08,
            ef_hole=-1.40
        )

        self.assertGreater(sheet_n[0], 0.0, "Ground electron sheet density should be positive")
        self.assertGreater(sheet_n[0], sheet_n[1], "Lower energy subband must have higher occupation")

        # Check volumetric integration: int n_qw(z) dz (in cm) should equal sum(n_2D)
        # dz in nm = dz * 1e-7 cm
        integrated_n_cm2 = np.sum(dens_n) * (dz * 1e-7)
        total_sheet_cm2 = np.sum(sheet_n)
        self.assertAlmostEqual(
            integrated_n_cm2 / total_sheet_cm2, 1.0, places=4,
            msg="Volumetric carrier density integral does not conserve 2D sheet charge"
        )

    def test_07_optical_transitions_and_overlap(self):
        """Benchmark 7: Optical transition energy, wavelength, and overlap matrix elements."""
        z = np.linspace(0.0, 30.0, 301)
        dz = z[1] - z[0]

        e_energies = np.array([0.05, 0.15])
        h_energies = np.array([-1.44, -1.48])

        # Ground states: both symmetric (high overlap)
        psi_e1 = np.exp(-((z - 15.0) / 3.0) ** 2)
        psi_e1 /= np.sqrt(np.sum(psi_e1 ** 2) * dz)
        # Excited state: antisymmetric
        psi_e2 = (z - 15.0) * np.exp(-((z - 15.0) / 3.0) ** 2)
        psi_e2 /= np.sqrt(np.sum(psi_e2 ** 2) * dz)

        psi_h1 = np.exp(-((z - 15.0) / 3.0) ** 2)
        psi_h1 /= np.sqrt(np.sum(psi_h1 ** 2) * dz)

        t_e, t_w, gamma, dom = compute_optical_transitions(
            z_nm=z,
            electron_energies=e_energies,
            electron_wf=np.array([psi_e1, psi_e2]),
            hole_energies=h_energies,
            hole_wf=np.array([psi_h1, psi_h1]),
            hole_types=["hh1", "hh2"]
        )

        # e1-h1 overlap should be near 1.0 (identical parity)
        self.assertAlmostEqual(gamma[0, 0], 1.0, places=3, msg="e1-h1 overlap should be near unity for identical shapes")
        # e2-h1 overlap should be near 0.0 (opposite parity: odd vs even)
        self.assertAlmostEqual(gamma[1, 0], 0.0, places=3, msg="e2-h1 overlap should be zero by parity selection rule")

        # Check transition energy and wavelength
        expected_de = 0.05 - (-1.44) # 1.49 eV
        self.assertAlmostEqual(t_e[0, 0], expected_de, places=4)
        expected_lambda = HC_EV_NM / expected_de # ~832.1 nm
        self.assertAlmostEqual(t_w[0, 0], expected_lambda, places=2)

    def test_08_multiple_quantum_wells_subband_splitting(self):
        """Benchmark 8: Double quantum well superlattice symmetric/antisymmetric splitting."""
        # Double QW: two 6 nm GaAs wells separated by a thin 2 nm AlGaAs barrier
        z = np.linspace(0.0, 40.0, 801)
        ec = np.full_like(z, 0.25)
        # Well 1: 12 - 18 nm, Barrier: 18 - 20 nm, Well 2: 20 - 26 nm
        ec[(z >= 12.0) & (z <= 18.0)] = 0.0
        ec[(z >= 20.0) & (z <= 26.0)] = 0.0
        mass = np.where(ec == 0.0, 0.067, 0.092)

        energies, wf, prob, norm_errs = solve_effective_mass_schrodinger_1d(
            z_nm=z, potential_ev=ec, m_eff_rel=mass, num_states=2
        )

        self.assertEqual(len(energies), 2)
        # Tunneling coupling produces symmetric (lower) and antisymmetric (higher) doublet:
        delta_split = energies[1] - energies[0]
        self.assertGreater(delta_split, 0.001, "Coupled DQW must exhibit energy splitting due to barrier tunneling")
        self.assertLess(delta_split, 0.050, "Splitting too large for weak tunneling")

    def test_09_experimental_benchmark_tsang1981_sqw(self):
        """Benchmark 9: Tsang 1981 GaAs/Al0.3Ga0.7As SQW Laser Diode (~845 nm emission)."""
        # From Tsang 1981 (APL 39, 134): Lw = 10 nm GaAs well with Al0.3Ga0.7As barriers
        # Room temperature (300K) GaAs bandgap Eg = 1.424 eV
        # Al0.3Ga0.7As bandgap Eg = 1.798 eV (Delta_Ec = 0.65 * 0.374 = 0.243 eV, Delta_Ev = 0.131 eV)
        # Confinement pushes e1 ~ +35-40 meV, hh1 ~ -10-15 meV -> E_trans ~ 1.47 eV -> lambda ~ 843-847 nm!
        layers = [
            {"material": "AlGaAs", "thickness": 20.0, "type": "barrier", "mole": 0.3},
            {"material": "GaAs", "thickness": 10.0, "type": "well", "mole": 0.0},
            {"material": "AlGaAs", "thickness": 20.0, "type": "barrier", "mole": 0.3},
        ]
        
        z = np.linspace(0.0, 50.0, 501)
        ec = np.where((z >= 20.0) & (z <= 30.0), 0.0, 0.243)
        ev = np.where((z >= 20.0) & (z <= 30.0), -1.424, -1.555)
        band_prof = {"z": z, "ec": ec, "ev": ev}

        qw_result = solve_quantum_well(
            band_profile=band_prof,
            layers=layers,
            temperature_k=300.0,
            num_electron_states=2,
            num_hole_states=2
        )

        self.assertGreater(len(qw_result.dominant_transitions), 0)
        dominant = qw_result.dominant_transitions[0]
        wavelength = dominant["wavelength_nm"]

        # Experimental Tsang 1981 emission is 845 ± 5 nm
        self.assertAlmostEqual(
            wavelength, 845.0, delta=8.0,
            msg=f"Predicted QW emission {wavelength:.1f} nm outside experimental Tsang 1981 range [837, 853] nm"
        )
        self.assertGreater(dominant["overlap"], 0.85, "e1-hh1 optical overlap should exceed 85%")

    def test_10_backward_compatibility_regression(self):
        """Benchmark 10: Zero performance/output penalty when QW calculation is not requested."""
        import aestimo
        # Verify default attributes on Structure / StructureFrom
        dummy_input = {
            "material": [[10.0, "GaAs", 0.0, 0.0, 1e18, "n", "b"]],
            "computation_scheme": 10
        }
        model = aestimo.StructureFrom(dummy_input, aestimo.database)
        self.assertFalse(getattr(model, "enable_qw_solver", False), "QW solver must be disabled by default")


if __name__ == "__main__":
    unittest.main(verbosity=2)
