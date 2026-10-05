> Source review is pending. CSV filenames and citations are inherited claims, not proof of measurement provenance. No original figures/tables or digitization projects have been independently matched in this review. See `docs/reference-data-audit.json` for every file's origin classification and SHA-256 fingerprint. Si I–V, Si C–V and InGaAs I–V are synthetic/model references; InGaN I–V is literature-based with unverified extraction. Numerical agreement, calibration history and experimental validation are distinct.

# Experimental Datasets for Device Validation in Aestimo 1D

This directory contains reference datasets attributed to literature, plus model-generated datasets and benchmark calibration files for validating electrical, photovoltaic, and optoelectronic simulations in **Aestimo 1D**.

---

## 1. Light-Emitting Diode (LED) Experimental Datasets

### A. Nakamura et al. (1995) Blue InGaN Single Quantum Well LED
- **Dataset Files**:
  - `led_nakamura1995_blue_sqw_iv.csv` (I-V and Optical Power L-I data)
  - `led_nakamura1995_blue_sqw_el_spectrum.csv` (Electroluminescence Emission Spectrum)
- **Bibliographic Reference**:
  - S. Nakamura, M. Senoh, N. Iwasa, and S. Nagahama, *"High-power InGaN single-quantum-well-structure blue and violet light-emitting diodes"*, **Applied Physics Letters**, Vol. 67, No. 13, pp. 1868–1870 (1995).
  - DOI: [10.1063/1.114359](https://doi.org/10.1063/1.114359)
- **Experimental Device Structure**:
  - Substrate: c-face sapphire with 30 nm GaN buffer
  - $n$-type GaN cladding: 4.0 μm, $N_d = 2.0\times 10^{18}\text{ cm}^{-3}$
  - Active SQW: 3.0 nm $\text{In}_{0.20}\text{Ga}_{0.80}\text{N}$ undoped
  - $p$-type $\text{Al}_{0.20}\text{Ga}_{0.80}\text{N}$ EBL: 100 nm, $N_a = 5.0\times 10^{17}\text{ cm}^{-3}$
  - $p$-type GaN contact cap: 500 nm, $N_a = 1.0\times 10^{18}\text{ cm}^{-3}$
  - Mesa active area: $350\,\mu\text{m} \times 350\,\mu\text{m} = 1.225\times 10^{-3}\text{ cm}^2$
- **Reported Experimental Observables at 300 K**:
  - Peak emission wavelength: $\lambda_{\text{peak}} = 450.0\text{ nm}$ ($h\nu = 2.755\text{ eV}$), FWHM $= 25.0\text{ nm}$
  - Operating voltage: $V_f = 3.6\text{ V}$ at $I = 20\text{ mA}$ ($J = 16.3\text{ A/cm}^2$)
  - Turn-on voltage: $V_{\text{on}} \approx 2.7\text{ V}$
  - Optical output power: $P_{\text{opt}} = 5.0\text{ mW}$ at 20 mA ($\text{EQE} = 9.2\%$)
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

---

### B. Meyaard et al. (2013) InGaN/GaN Multiple Quantum Well Droop LED
- **Dataset Files**:
  - `led_meyaard2013_mqw_iv.csv` (Current-Voltage characteristics)
  - `led_meyaard2013_mqw_droop_iqe.csv` (Normalized IQE Droop vs Current Density)
- **Bibliographic Reference**:
  - D. S. Meyaard, G.-B. Lin, J. Cho, E. F. Schubert, H. Shim, S.-H. Han, M.-H. Kim, C. Sone, and Y. S. Kim, *"Identifying the cause of the efficiency droop in GaInN light-emitting diodes by correlating the onset of high injection with the onset of the efficiency droop"*, **Applied Physics Letters**, Vol. 102, No. 25, 251114 (2013).
  - DOI: [10.1063/1.4811558](https://doi.org/10.1063/1.4811558)
- **Experimental Device Structure**:
  - $n$-type GaN layer: $N_d = 3.0\times 10^{18}\text{ cm}^{-3}$
  - Active MQW: 5 pairs of 3.0 nm $\text{In}_{0.15}\text{Ga}_{0.85}\text{N}$ wells separated by 10.0 nm GaN barriers (undoped)
  - $p$-type $\text{Al}_{0.15}\text{Ga}_{0.85}\text{N}$ EBL: 20 nm, $N_a = 3.0\times 10^{17}\text{ cm}^{-3}$
  - $p$-type GaN contact cap: 200 nm, $N_a = 1.0\times 10^{18}\text{ cm}^{-3}$
  - Mesa active area: $300\,\mu\text{m} \times 300\,\mu\text{m} = 9.0\times 10^{-4}\text{ cm}^2$
- **Reported Experimental Observables at 300 K**:
  - Turn-on voltage: $V_{\text{on}} \approx 2.6\text{ V}$, $V_f \approx 3.32\text{ V}$ at 20 mA ($J = 22.2\text{ A/cm}^2$)
  - Peak emission wavelength: $\lambda_{\text{peak}} = 445.0\text{ nm}$ ($h\nu = 2.786\text{ eV}$)
  - Efficiency droop curve: normalized $\text{IQE}(J)$ peaking at $J_{\text{peak}} \approx 12.0\text{ A/cm}^2$, dropping to $67\%$ at $100\text{ A/cm}^2$ and $52\%$ at $200\text{ A/cm}^2$
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

---

### C. Schubert (2006) / Steranka (1988) AlGaAs/GaAs Double-Heterostructure Infrared LED
- **Dataset Files**:
  - `led_schubert2006_algaas_dh_iv.csv` (I-V characteristics)
  - `led_schubert2006_algaas_dh_el_spectrum.csv` (Electroluminescence Emission Spectrum)
- **Bibliographic Reference**:
  - E. F. Schubert, *Light-Emitting Diodes*, 2nd Edition, Cambridge University Press (2006), Chapters 5 & 6.
  - G. E. Steranka et al., *"High-efficiency AlGaAs light-emitting diodes"*, **Hewlett-Packard Journal**, Vol. 39, pp. 84–88 (1988).
- **Experimental Device Structure**:
  - $n$-type $\text{Al}_{0.35}\text{Ga}_{0.65}\text{As}$ cladding: 1000 nm, $N_d = 1.0\times 10^{18}\text{ cm}^{-3}$
  - Active GaAs region: 100 nm undoped GaAs ($E_g = 1.424\text{ eV}$)
  - $p$-type $\text{Al}_{0.35}\text{Ga}_{0.65}\text{As}$ cladding: 1000 nm, $N_a = 1.0\times 10^{18}\text{ cm}^{-3}$
  - $p$-type GaAs contact cap: 50 nm, $N_a = 5.0\times 10^{18}\text{ cm}^{-3}$
  - Mesa active area: $250\,\mu\text{m} \times 250\,\mu\text{m} = 6.25\times 10^{-4}\text{ cm}^2$
- **Reported Experimental Observables at 300 K**:
  - Turn-on voltage: $V_{\text{on}} \approx 1.35\text{ V}$, $V_f \approx 1.55\text{ V}$ at 20 mA ($J = 32.0\text{ A/cm}^2$)
  - Peak emission wavelength: $\lambda_{\text{peak}} = 870.0\text{ nm}$ ($h\nu = 1.425\text{ eV}$), FWHM $= 35.0\text{ nm}$
  - Internal Quantum Efficiency: $\text{IQE} \approx 85-90\%$
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

---

## 2. Semiconductor Laser Diode Benchmark Datasets

### A. Tsang (1981) GaAs/AlGaAs Single Quantum Well GRIN-SCH Laser Diode (~845 nm)
- **Dataset Files**:
  - `laser_tsang1981_gaas_sqw_li.csv` (Output Power vs. Injection Current L-I)
  - `laser_tsang1981_gaas_sqw_iv.csv` (Current-Voltage characteristics)
  - `laser_tsang1981_gaas_sqw_spectrum.csv` (Longitudinal Fabry-Pérot Mode Optical Spectrum)
- **Bibliographic Reference**:
  - W. T. Tsang, *"Extremely low threshold GaAs-AlGaAs graded-index waveguide separate-confinement heterostructure lasers grown by molecular beam epitaxy"*, **Applied Physics Letters**, Vol. 39, No. 2, pp. 134–137 (1981). DOI: [10.1063/1.92690](https://doi.org/10.1063/1.92690)
  - Cross-Reference: L. A. Coldren et al., *Diode Lasers and Photonic Integrated Circuits* (2012), Table 2.1; Derry et al., *IEEE JQE* 28, 2698 (1992).
- **Experimental Device Structure**:
  - Material system: $\text{GaAs/AlGaAs}$ (Zincblende)
  - Waveguide: GRIN-SCH SQW ($d_w = 8.0\text{ nm}$ GaAs well, $\text{Al}_{0.25}\text{Ga}_{0.75}\text{As}$ guide, $\text{Al}_{0.60}\text{Ga}_{0.40}\text{As}$ cladding)
  - Cavity dimensions: Length $L = 500\,\mu\text{m}$, stripe width $w = 20\,\mu\text{m}$ (area $1.0\times 10^{-4}\text{ cm}^2$)
  - As-cleaved uncoated facets: $R_1 = R_2 = 0.32$ ($\alpha_m = 22.8\text{ cm}^{-1}$), $\alpha_i \approx 3.5\text{ cm}^{-1}$
- **Reported Experimental Observables at 300 K**:
  - Threshold current: $I_{\text{th}} \approx 20.0\text{ mA}$ ($J_{\text{th}} \approx 200\text{ A/cm}^2$)
  - Slope efficiency: $SE \approx 0.41\text{ mW/mA}$ per facet ($0.82\text{ mW/mA}$ total), $\eta_d \approx 58\%$
  - Forward threshold voltage: $V_{\text{th}} \approx 1.55\text{ V}$, $R_s \approx 2.5\,\Omega$
  - Emission wavelength: $\lambda_{\text{peak}} = 845.0\text{ nm}$, longitudinal mode spacing $\Delta\lambda \approx 0.198\text{ nm}$
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

---

### B. Zah et al. (1994) InGaAsP/InP 1.55 μm Strained MQW Telecommunication Laser
- **Dataset Files**:
  - `laser_zah1994_1550nm_mqw_li.csv` (Multi-Temperature L-I at 25°C, 50°C, 85°C)
  - `laser_zah1994_1550nm_mqw_iv.csv` (Current-Voltage characteristics)
  - `laser_zah1994_1550nm_mqw_spectrum.csv` (Longitudinal Mode Optical Spectrum)
- **Bibliographic Reference**:
  - C.-E. Zah, R. Bhat, B. N. Pathak, F. Favire, W. Lin, M. C. Wang, N. C. Andreadakis, D. M. Hwang, M. A. Koza, T.-P. Lee, Z. Wang, D. Darby, D. Flanders, and J. J. Hsieh, *"High-performance uncooled 1.3-μm and 1.55-μm strained-layer quantum-well lasers"*, **IEEE Journal of Quantum Electronics**, Vol. 30, No. 2, pp. 511–523 (1994). DOI: [10.1109/3.283799](https://doi.org/10.1109/3.283799)
- **Experimental Device Structure**:
  - Material system: $\text{InGaAsP/InP}$ (Zincblende)
  - Active region: 5 pairs of compressively strained $\text{In}_{0.70}\text{Ga}_{0.30}\text{As}_{0.82}\text{P}_{0.18}$ wells ($d_w = 6.0\text{ nm}$) with tensile-strained $1.25\text{Q}$ barriers ($d_b = 9.0\text{ nm}$)
  - Cavity dimensions: Length $L = 500\,\mu\text{m}$, ridge width $w = 2.5\,\mu\text{m}$ (area $1.25\times 10^{-5}\text{ cm}^2$)
  - As-cleaved facets: $R_1 = R_2 = 0.32$ ($\alpha_m = 22.8\text{ cm}^{-1}$), $\alpha_i \approx 8.0\text{ cm}^{-1}$
- **Reported Experimental Observables**:
  - Threshold current: $I_{\text{th}}(25^\circ\text{C}) \approx 12.5\text{ mA}$, $I_{\text{th}}(50^\circ\text{C}) \approx 20.5\text{ mA}$, $I_{\text{th}}(85^\circ\text{C}) \approx 41.0\text{ mA}$
  - Characteristic temperature: $T_0 \approx 52\text{ K}$ (Auger recombination dominance)
  - Slope efficiency: $SE \approx 0.24\text{ mW/mA}$ per facet ($\eta_d \approx 30\%$)
  - Forward threshold voltage: $V_{\text{th}} \approx 0.95\text{ V}$, $R_s \approx 4.0\,\Omega$
  - Emission wavelength: $\lambda_{\text{peak}} = 1550.0\text{ nm}$, mode spacing $\Delta\lambda \approx 0.677\text{ nm}$
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

---

### C. Nakamura et al. (1996) InGaN/GaN/AlGaN Multi-Quantum-Well Violet-Blue Laser Diode (~405 nm)
- **Dataset Files**:
  - `laser_nakamura1996_blue_mqw_li.csv` (Output Power vs. Injection Current L-I)
  - `laser_nakamura1996_blue_mqw_iv.csv` (Current-Voltage characteristics)
  - `laser_nakamura1996_blue_mqw_spectrum.csv` (Longitudinal Mode Optical Spectrum)
- **Bibliographic Reference**:
  - S. Nakamura, M. Senoh, S. Nagahama, N. Iwasa, T. Yamada, T. Matsushita, H. Kiyoku, and Y. Sugimoto, *"InGaN-based multi-quantum-well-structure laser diodes"*, **Applied Physics Letters**, Vol. 68, Issue 15, pp. 2105–2107 (1996). DOI: [10.1063/1.116084](https://doi.org/10.1063/1.116084)
  - Cross-Reference: S. Nakamura et al., *Appl. Phys. Lett.* 69, 4056 (1996).
- **Experimental Device Structure**:
  - Material system: $\text{InGaN/GaN/AlGaN}$ (Wurtzite)
  - Active region: 5-period $\text{In}_{0.15}\text{Ga}_{0.85}\text{N}$ ($3.5\text{ nm}$) / $\text{In}_{0.02}\text{Ga}_{0.98}\text{N}$ ($7.0\text{ nm}$) MQW with $20\text{ nm}$ $p\text{-Al}_{0.20}\text{Ga}_{0.80}\text{N}$ EBL
  - Cavity dimensions: Length $L = 600\,\mu\text{m}$, stripe width $w = 5.0\,\mu\text{m}$ (area $3.0\times 10^{-5}\text{ cm}^2$)
  - As-cleaved/etched facets: $R_1 = R_2 \approx 0.18$ ($\alpha_m = 28.5\text{ cm}^{-1}$), $\alpha_i \approx 15.0\text{ cm}^{-1}$
- **Reported Experimental Observables at 300 K**:
  - Threshold current: $I_{\text{th}} \approx 80.0\text{ mA}$ ($J_{\text{th}} \approx 2.67\text{ kA/cm}^2$)
  - Slope efficiency: $SE \approx 0.40\text{ mW/mA}$ per facet, $\eta_d \approx 26\%$
  - Forward threshold voltage: $V_{\text{th}} \approx 5.5\text{ V}$, $R_s \approx 28.0\,\Omega$
  - Emission wavelength: $\lambda_{\text{peak}} = 405.0\text{ nm}$, mode spacing $\Delta\lambda \approx 0.0547\text{ nm}$
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

---

## 3. Solar Cell Benchmark Datasets

### Tobin et al. (1990) High-Efficiency GaAs Solar Cell
- **File**: `gaas_tobin1990_experimental_iv.csv`
- **Reference**: S. P. Tobin, S. M. Vernon, C. Bajgar, S. J. Wojtczuk, M. R. Melloch, A. Keshavarzi, T. B. Stellwag, S. Venkatensan, M. S. Lundstrom, and K. A. Emery, *"Assessment of MOCVD- and MBE-grown GaAs for high-efficiency solar cells"*, **IEEE Transactions on Electron Devices**, Vol. 37, No. 2, pp. 469–477 (1990). DOI: [10.1109/16.46369](https://doi.org/10.1109/16.46369)
- **Key Metrics (1-Sun AM1.5G, 300 K)**: $J_{\text{sc}} = 27.80\text{ mA/cm}^2, V_{\text{oc}} = 1.052\text{ V}, \text{FF} = 84.91\%, \eta = 24.83\%$.

---

## 3. Junction Diode Reference Datasets

### A. Silicon $p\text{-}n$ Junction Diode
- **Files**: `si_pn_experimental_iv.csv`, `si_pn_experimental_cv.csv`
- **Physics**: Shockley minority carrier diffusion in high-purity Silicon with realistic contact parasitics ($N_a = 10^{18}\text{ cm}^{-3}, N_d = 10^{18}\text{ cm}^{-3}, R_s = 5.0\,\Omega, R_{sh} = 100\text{ M}\Omega, A = 10^{-4}\text{ cm}^2$).

### B. InGaAs $p\text{-}n$ Junction Diode
- **File**: `ingaas_pn_experimental_iv.csv`
- **Physics**: Narrow-bandgap $\text{In}_{0.53}\text{Ga}_{0.47}\text{As}$ lattice-matched to InP ($E_g = 0.736\text{ eV}, R_s = 10.0\,\Omega, R_{sh} = 100\text{ k}\Omega, A = 10^{-4}\text{ cm}^2$).

### C. InGaN $p\text{-}n$ Junction Diode
- **File**: `ingan_pn_experimental_iv.csv`
- **Physics**: Wurtzite $\text{In}_{0.57}\text{Ga}_{0.43}\text{N}$ homojunction with deep Mg acceptor compensation and defect-assisted shunt mechanisms ($E_g \approx 1.35\text{ eV}, R_s = 50.0\,\Omega, R_{sh} = 150\text{ k}\Omega, A = 5\times 10^{-4}\text{ cm}^2$).

---

## 4. Quantum-Well (QW) Confined-State Experimental Datasets

### A. Miller et al. (1984, 1985) Quantum-Confined Stark Effect (QCSE) in GaAs/AlGaAs
- **Dataset File**: `qw_miller1984_qcse_stark_shift.csv`
- **Bibliographic Reference**:
  - D. A. B. Miller, D. S. Chemla, T. C. Damen, A. C. Gossard, W. Wiegmann, T. H. Wood, and C. A. Burrus, *"Band-Edge Electroabsorption in Quantum Well Structures: The Quantum-Confined Stark Effect"*, **Physical Review Letters**, Vol. 53, No. 22, pp. 2173–2176 (1984). DOI: [10.1103/PhysRevLett.53.2173](https://doi.org/10.1103/PhysRevLett.53.2173)
  - D. A. B. Miller, D. S. Chemla, T. C. Damen, A. C. Gossard, W. Wiegmann, T. H. Wood, and C. A. Burrus, *"Electric field dependence of optical absorption near the band gap of quantum-well structures"*, **Physical Review B**, Vol. 32, No. 2, pp. 1043–1060 (1985). DOI: [10.1103/PhysRevB.32.1043](https://doi.org/10.1103/PhysRevB.32.1043)
- **Experimental System**:
  - $9.5\text{ nm}$ GaAs quantum well embedded between $\text{Al}_{0.32}\text{Ga}_{0.68}\text{As}$ barriers at $T = 300\text{ K}$.
  - Applied perpendicular electric field $\mathcal{E}$ from $0$ to $110\text{ kV/cm}$.
- **Reported Experimental Observables**:
  - Heavy-hole exciton transition energy $E(e_1\text{--}hh_1)$ shifting from $1.4550\text{ eV}$ down to $1.4113\text{ eV}$ (Stark red-shift of $-43.7\text{ meV}$).
  - Light-hole exciton transition energy $E(e_1\text{--}lh_1)$ shifting from $1.4720\text{ eV}$ down to $1.4319\text{ eV}$ (Stark red-shift of $-40.1\text{ meV}$).
  - Progressive reduction of electron-hole envelope wavefunction overlap integral $\Gamma_{11}$ due to spatial field-induced separation.
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

### B. Dingle (1974, 1975) Quantum Confinement Energy vs Well Width in GaAs/AlGaAs
- **Dataset File**: `qw_dingle1975_energy_vs_width.csv`
- **Bibliographic Reference**:
  - R. Dingle, *"Confined Carrier Quantum States in Ultrathin Semiconductor Heterostructures"*, **Festkörperprobleme / Advances in Solid State Physics**, Vol. 15, pp. 21–48 (1975).
  - R. Dingle, W. Wiegmann, and C. H. Henry, *"Quantum States of Confined Carriers in Very Thin $\text{Al}_x\text{Ga}_{1-x}\text{As}\text{-GaAs-Al}_x\text{Ga}_{1-x}\text{As}$ Heterostructures"*, **Physical Review Letters**, Vol. 33, No. 14, pp. 827–830 (1974). DOI: [10.1103/PhysRevLett.33.827](https://doi.org/10.1103/PhysRevLett.33.827)
- **Experimental System**:
  - MBE-grown GaAs quantum wells with thickness $L_w$ systematically varied from $2.5\text{ nm}$ to $25.0\text{ nm}$ between $\text{Al}_{0.30}\text{Ga}_{0.70}\text{As}$ barriers at $T = 300\text{ K}$.
- **Reported Experimental Observables**:
  - Fundamental transition energy $E(e_1\text{--}hh_1)$ shifting from $1.4248\text{ eV}$ ($L_w = 25\text{ nm}$) to $1.6420\text{ eV}$ ($L_w = 2.5\text{ nm}$).
  - Light-hole transition energy $E(e_1\text{--}lh_1)$ shifting from $1.4278\text{ eV}$ ($L_w = 25\text{ nm}$) to $1.6950\text{ eV}$ ($L_w = 2.5\text{ nm}$).
  - Second subband transition energy $E(e_2\text{--}hh_2)$ shifting from $1.4390\text{ eV}$ ($L_w = 25\text{ nm}$) to $1.7850\text{ eV}$ ($L_w = 2.5\text{ nm}$).
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

### C. Tsang (1981) Single Quantum Well Photoluminescence Spectrum
- **Dataset File**: `qw_tsang1981_sqw_photoluminescence.csv`
- **Bibliographic Reference**:
  - W. T. Tsang, *"Extremely low threshold GaAs/AlGaAs graded-index waveguide separate-confinement heterostructure lasers grown by MBE"*, **Applied Physics Letters**, Vol. 39, No. 2, pp. 134–137 (1981). DOI: [10.1063/1.92697](https://doi.org/10.1063/1.92697)
- **Experimental System**:
  - $8.0\text{ nm}$ GaAs Single Quantum Well embedded in a parabolic GRIN-SCH optical cavity at $300\text{ K}$.
- **Reported Experimental Observables**:
  - Room-temperature quantum-well photoluminescence/spontaneous emission peaking at $\lambda = 845.0\text{ nm}$ ($h\nu = 1.467\text{ eV}$) with $\text{FWHM} = 6.0\text{ nm}$.
- **Validation Status**: `REFERENCE COMPARISON / PROVENANCE UNVERIFIED`

---

## 5. Parameter Provenance Classification Standard

Every parameter used in validated simulation models is classified according to four strict tiers:
- **`EXPERIMENTAL`**: Directly reported or measured in peer-reviewed publications.
- **`DERIVED`**: Calculated from experimental measurements via fundamental semiconductor physical relationships (e.g. built-in voltage, subband quantization, dielectric permittivity).
- **`FITTED`**: Determined through rigorous numerical optimization with documented objective functions and physical bounds (e.g. contact series resistance).
- **`ASSUMED`**: Physics-based literature estimates when unmeasured (e.g. planar light extraction efficiency).
