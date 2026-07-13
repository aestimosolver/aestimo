# Baseline Calibration Report

Authoritative source: `examples/untitled_project.json`

## Baseline Structure

- Layer stack: 2 layer high-indium InGaN homojunction
- Composition: InGaN with x = 0.57
- Thickness: 100.0 nm total
- Doping: p = 1.000e+16 cm^-3, n = 2.000e+17 cm^-3
- Contacts: left = 0.000 V, right = 0.600 V

## Modified Parameters

- Dark baseline: `G_optical` changed from 1.000e+21 to `0.0 cm^-3 s^-1`
- Leakage diagnostic only: voltage window changed from `0.00 to 1.60 V` to `-1.00 to 0.00 V`

## Comparison Against High-Indium Literature Trend

| Metric | Simulation | Literature-Based Reference | Assessment |
| --- | ---: | ---: | --- |
| Reverse leakage at -0.95 V (A/cm^2) | 6.878e-10 | 9.500e-03 | Too low by 7.14 decades |
| Turn-on at |J| = 0.5 A/cm^2 (V) | >1.60 | 0.616 | Too late / too soft |
| Ideality factor fit (0.10 V to 0.30 V) | 1.605 | 2.617 | Too ideal; not defect-dominated enough |

## Interpretation

- The equilibrium band diagram is internally consistent, but the dark transport does not match the defect-limited behavior expected for high-indium InGaN homojunction devices.
- Reverse leakage is many orders of magnitude lower than the literature-based experimental trend, which points to insufficient defect-assisted leakage and tunneling in the baseline.
- The forward knee is too delayed, and the fitted ideality factor is closer to diffusion-limited behavior than the SRH-dominated regime typically reported for high-indium material.

## Validation Decision

- Baseline physically reasonable: `False`
- Parameter sweeps should remain blocked until the baseline is recalibrated against the high-indium reference regime.

## Reference Basis

- Numerical comparison file: `examples/experimental_data/ingan_pn_experimental_iv.csv`
- This repository reference file encodes a literature-based high-indium InGaN p-n homojunction trend with high defect density, elevated leakage, and ideality above 2.
