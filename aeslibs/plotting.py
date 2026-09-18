import os
import sys
import numpy as np
import matplotlib.pyplot as pl

# Ensure project root is in sys.path
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

try:
    import config
except ImportError:
    import types
    config = types.ModuleType("config")

# Robust fallback defaults for all config attributes
_config_defaults = {
    'damping': 0.2,
    'Stern_damping': True,
    'max_iterations': 80,
    'convergence_test': 1e-4,
    'predic_correc': True,
    'anti_crossing_length': 0.0001,
    'amort_wave_0': 1.5,
    'amort_wave_1': 1.5,
    'strain': True,
    'piezo': False,
    'piezo1': True,
    'quantum_effect': True,
    'parameters': True,
    'electricfield_out': True,
    'potential_out': True,
    'sigma_out': True,
    'probability_out': True,
    'states_out': True,
    'Drift_Diffusion_out': True,
    'wavefunction_scalefactor': 400.0
}
for _attr, _val in _config_defaults.items():
    if not hasattr(config, _attr):
        setattr(config, _attr, _val)

from aeslibs.aestimo_poisson1d import amort_wave

# Defining constants and material parameters
q = 1.602176e-19  # C
kb = 1.3806504e-23  # J/K
m_e = 9.1093826e-31  # kg
T = 300.0  # Kelvin
Vt = kb * T / q  # [eV]
J2meV = 1e3 / q  # Joules to meV


def save_and_plot2(result, model, output_directory='output', drawFigures=False, show=True):
    xaxis = result.xaxis

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    def saveoutput(fname, datatuple, header=""):
        fname2 = os.path.join(output_directory, fname)
        np.savetxt(
            fname2, np.column_stack(datatuple), fmt="%.6e", delimiter=" ", header=header
        )

    # Plotting results
    if getattr(config, 'Drift_Diffusion_out', True):
         saveoutput("av_curr.dat",(result.Va_t, result.av_curr))
    
    for k in range(0, result.Total_Steps):
        vt = result.Va_t[k]
        if getattr(config, 'Drift_Diffusion_out', True):
            if getattr(config, 'sigma_out', True):
                # Use 2D array if available, else fallback to 1D (for backward compatibility or non-sweep)
                ro_data = result.ro_result_[k,:] if hasattr(result, 'ro_result_') else result.ro_result
                saveoutput("sigma_eh_%.2f.dat" % vt, (xaxis, ro_data))
                
            if config.electricfield_out:
                e1_data = result.el_field1_result_[k,:] if hasattr(result, 'el_field1_result_') else result.el_field1_result
                e2_data = result.el_field2_result_[k,:] if hasattr(result, 'el_field2_result_') else result.el_field2_result
                saveoutput(
                    "efield_eh_%.2f.dat" % vt,
                    (xaxis, e1_data, e2_data),
                )
                
            if config.potential_out:
                ec_data = result.Ec_result_[k,:] if hasattr(result, 'Ec_result_') else result.Ec_result
                ev_data = result.Ev_result_[k,:] if hasattr(result, 'Ev_result_') else result.Ev_result
                nf_data = result.nf_result_[k,:] if hasattr(result, 'nf_result_') else result.nf_result
                pf_data = result.pf_result_[k,:] if hasattr(result, 'pf_result_') else result.pf_result
                
                saveoutput(
                    "potn_eh_%.2f.dat" % vt,
                    (xaxis, ec_data, ev_data),
                )
                saveoutput(
                    "np_data0_%.2f.dat" % vt,
                    (xaxis, nf_data * 1e-6, pf_data * 1e-6),
                )
            
            if config.states_out:
                for j in range(1, result.N_wells_virtual - 1):
                    I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
                    i1 = I1 - I1
                    i2 = I2 - I1
                    rel_meff_state = [
                        meff / m_e for meff in result.meff_state_general[j]
                    ]  # going to report relative effective mass.
                    columns = (
                        range(model.subnumber_h),
                        result.E_state_general[j],
                        result.N_state_general[j],
                        rel_meff_state,
                    )
                    # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
                    header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
                    saveoutput(
                        "states_h_QWR%d_%.2f.dat" % (j, vt), columns, header=header
                    )
                    if config.probability_out:
                        saveoutput(
                            "wavefunctions_h_QWR%d_%.2f.dat" % (j, vt),
                            (xaxis[I1:I2], result.wfh_general[j,:,i1:i2].transpose()),
                        )
    
            if config.states_out:
                for j in range(1, result.N_wells_virtual - 1):
                    I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
                    i1 = I1 - I1
                    i2 = I2 - I1
                    rel_meff_statec = [
                        meff / m_e for meff in result.meff_statec_general[j]
                    ]  # going to report relative effective mass.
                    columns = (
                        range(model.subnumber_e),
                        result.E_statec_general[j],
                        result.N_statec_general[j],
                        rel_meff_statec,
                    )
                    # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
                    header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
                    saveoutput(
                        "states_e_QWR%d_%.2f.dat" % (j, vt), columns, header=header
                    )
                    if config.probability_out:
                        saveoutput(
                            "wavefunctions_e_QWR%d_%.2f.dat" % (j, vt),
                            (xaxis[I1:I2], result.wfe_general[j,:,i1:i2].transpose()),
                        )
    if drawFigures:
        span = np.ones(100000000)
        fig1 = pl.figure(figsize=(10, 8))
        pl.suptitle("Aestimo Results")
        pl.subplot(1, 1, 1)
        pl.plot(
            xaxis * 1e6,
            result.Ec_result,
            xaxis * 1e6,
            result.Ev_result,
            xaxis * 1e6,
            result.Ei_result,
            xaxis * 1e6,
            result.Efn_result,
            "r",
            xaxis * 1e6,
            result.Efp_result,
            "b",
        )
        if model.N_wells_virtual - 2 != 0:
            for j in range(1, result.N_wells_virtual - 1):
                I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
                i1 = I1 - I1
                i2 = I2 - I1
                for levelc, statec in zip(
                    result.E_statec_general[j, :], result.wfe_general[j, :, :]
                ):
                    # pl.axhline(levelc,0.1,0.9, hold=True,color='g',ls='--')
                    pl.plot(
                        xaxis[I1:I2] * 1e6,
                        statec[i1:i2] * config.wavefunction_scalefactor * 1e-3
                        + levelc * 1e-3,
                        "b",
                    )
                    pl.plot(
                        xaxis[I1:I2] * 1e6, levelc * span[I1:I2] * 1e-3, "g", ls="--"
                    )
                for level, state in zip(
                    result.E_state_general[j, :], result.wfh_general[j, :, :]
                ):
                    # pl.axhline(level,0.1,0.9,color='g',ls='--')
                    pl.plot(
                        xaxis[I1:I2] * 1e6,
                        state[i1:i2] * config.wavefunction_scalefactor * 1e-3
                        + level * 1e-3,
                        "b",
                    )
                    pl.plot(
                        xaxis[I1:I2] * 1e6, level * span[I1:I2] * 1e-3, "g", ls="--"
                    )
        pl.xlabel("x [um]")
        pl.ylabel("Energy [eV]")
        pl.title("Quasi Fermi Levels (Efn (red) & Efp (blue)) vs Position", fontsize=12)
        pl.legend(("Ec", "Ev", "Ei", "Efn", "Efp"), loc="best", fontsize=12)
        pl.grid(True)

        fig2 = pl.figure(figsize=(10, 8))
        pl.suptitle(
            "1D Drift Diffusion Model Results - at Applied Bias (%.2f)"
            % vt,
            fontsize=12,
        )
        pl.subplots_adjust(hspace=0.4, wspace=0.4)

        pl.subplot(2, 2, 1)
        pl.plot(xaxis * 1e6, result.ro_result * 1e-6)
        pl.xlabel("x [um]")
        pl.ylabel("Total Charge Density [C/cm^3]")
        pl.title("Total Charge Density vs Position ", fontsize=12)
        pl.legend(("Total Charge"), loc="best", fontsize=12)
        pl.grid(True)
        # Plotting Efield
        # figure(1)
        pl.subplot(2, 2, 2)
        pl.plot(
            xaxis * 1e6,
            result.el_field1_result * 1e-8,
            "r",
            xaxis * 1e6,
            result.el_field2_result * 1e-8,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Electric Field 1(red) & 2 (blue) [MV/cm]")
        pl.title("Field Profile vs Position ", fontsize=12)
        pl.legend(("Electric Field 1", "Electric Field 2"), loc="best", fontsize=12)
        pl.grid(True)
        # Plotting Potential
        # figure(2)
        pl.subplot(2, 2, 3)
        pl.plot(xaxis * 1e6, result.Ec_result)
        pl.xlabel("x [um]")
        pl.ylabel("Conduction Band Energy (eV)")
        pl.title("Conduction Band vs Position ", fontsize=12)
        pl.legend(("Conduction Band"), loc="best", fontsize=12)
        pl.grid(True)
        # Plotting State(s)
        # figure(3)
        pl.subplot(2, 2, 4)
        pl.plot(result.Va_t, result.av_curr, 'b-', linewidth=2.0)
        pl.xlabel("Applied Bias (V)", fontsize=11)
        pl.ylabel("Total Current Density (mA/cm²)", fontsize=11)
        pl.title("J-V Characteristic", fontsize=12, fontweight='bold')
        pl.grid(True, which='both', linestyle='--', alpha=0.4)
        if show:
            pl.show()

        fig3 = pl.figure(figsize=(10, 8))
        pl.suptitle(
            "1D Drift Diffusion Model Results - at Applied Bias (%.2f V)"
            % vt,
            fontsize=12,
        )
        pl.subplots_adjust(hspace=0.4, wspace=0.4)
        pl.subplot(2, 2, 1)
        pl.semilogy(
            xaxis * 1e6,
            np.maximum(result.nf_result * 1e-6, 1.0),
            "r-",
            label="Electron (n)",
            linewidth=1.8,
        )
        pl.semilogy(
            xaxis * 1e6,
            np.maximum(result.pf_result * 1e-6, 1.0),
            "b-",
            label="Hole (p)",
            linewidth=1.8,
        )
        pl.xlabel("Position (µm)", fontsize=11)
        pl.ylabel("Carrier Density (cm⁻³)", fontsize=11)
        pl.title("Carrier Concentrations (Semi-log)", fontsize=12, fontweight='bold')
        pl.legend(loc="best", fontsize=10)
        pl.grid(True, which='both', linestyle='--', alpha=0.4)

        pl.subplot(2, 2, 2)
        pl.plot(
            xaxis * 1e6,
            result.Ec_result,
            xaxis * 1e6,
            result.Ev_result,
            xaxis * 1e6,
            result.Ei_result,
            xaxis * 1e6,
            result.Efn_result,
            "r",
            xaxis * 1e6,
            result.Efp_result,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Energy [eV]")
        pl.title("Quasi Fermi Levels (Efn (red) & Efp (blue)) vs Position", fontsize=12)
        pl.legend(("Ec", "Ev", "Ei", "Efn", "Efp"), loc="best", fontsize=12)
        pl.grid(True)

        pl.subplot(2, 2, 3)
        pl.plot(xaxis * 1e6, Vt * result.fi_result)
        pl.xlabel("x [um]")
        pl.ylabel("Potential [V]")
        pl.title("Potential vs Position", fontsize=12)
        pl.legend(("fi"), loc="best", fontsize=12)
        pl.grid(True)

        pl.subplot(2, 2, 4)
        pl.plot(
            xaxis * 1e6, result.Efn_result, "r", xaxis * 1e6, result.Efp_result, "b"
        )
        pl.xlabel("x [um]")
        pl.ylabel("Energy [eV]")
        pl.title("Quasi Fermi Levels (Efn (red) & Efp (blue)) vs Position", fontsize=12)
        pl.legend(("Efn", "Efp"), loc="best", fontsize=12)
        pl.grid(True)
        if show:
            pl.show()
    
    if drawFigures:
        return [fig1, fig2, fig3]
    else:
        return [None, None, None]


def save_and_plot(result, model, output_directory='output', drawFigures=False, show=True):

    xaxis = result.xaxis
    
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    def saveoutput(fname, datatuple, header=""):
        fname2 = os.path.join(output_directory, fname)
        np.savetxt(
            fname2, np.column_stack(datatuple), fmt="%.6e", delimiter=" ", header=header
        )

    if config.sigma_out:
        saveoutput("sigma_eh_equi_cond.dat", (xaxis, result.ro_result))
    if config.electricfield_out:
        saveoutput(
            "efield_eh_equi_cond.dat", (xaxis, result.el_field1_result, result.el_field2_result)
        )
    if config.potential_out:
        saveoutput("potn_eh_equi_cond.dat", (xaxis, result.fitotc / q, result.fitot / q))
        saveoutput(
            "np_data0_equi_cond.dat",
            (xaxis, result.nf_result * 1e-6, result.pf_result * 1e-6),
        )
    if config.states_out:
        for j in range(1, result.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
            i1 = I1 - I1
            i2 = I2 - I1
            rel_meff_state = [
                meff / m_e for meff in result.meff_state_general[j]
            ]  # going to report relative effective mass.
            columns = (
                range(model.subnumber_h),
                result.E_state_general[j],
                result.N_state_general[j],
                rel_meff_state,
            )
            # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
            header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
            saveoutput("states_h_QWR%d_equi_cond.dat" % j, columns, header=header)
            if config.probability_out:
                saveoutput(
                    "wavefunctions_h_QWR%d_equi_cond.dat" % j,
                    (xaxis[I1:I2], result.wfh_general[j,:,i1:i2].transpose()),
                )
    if config.states_out:
        for j in range(1, result.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
            i1 = I1 - I1
            i2 = I2 - I1
            rel_meff_statec = [
                meff / m_e for meff in result.meff_statec_general[j]
            ]  # going to report relative effective mass.
            columns = (
                range(model.subnumber_e),
                result.E_statec_general[j],
                result.N_statec_general[j],
                rel_meff_statec,
            )
            # header = " ".join([col.ljust(12) for col in ("State No.","Energy (meV)","N (m**-2)","Subband m* (m_e)")])
            header = "State No.    Energy (meV) N (m**-2)    Subband m* (kg)"
            saveoutput("states_e_QWR%d_equi_cond.dat" % j, columns, header=header)
            if config.probability_out:
                saveoutput(
                    "wavefunctions_e_QWR%d_equi_cond.dat" % j,
                    (xaxis[I1:I2], result.wfe_general[j,:,i1:i2].transpose()),
                )
    # Resultviewer
    if drawFigures:
        span = np.ones(100000000)
        fig1 = pl.figure(figsize=(10, 8))
        pl.suptitle("Aestimo Results - at Equilibrium Condition")
        pl.subplot(1, 1, 1)
        pl.plot(xaxis, result.fitot * J2meV, "k", xaxis, result.fitotc * J2meV, "k")
        for j in range(1, result.N_wells_virtual - 1):
            I1, I2, I11, I22 = amort_wave(j, result.Well_boundary, model.n_max)
            i1 = I1 - I1
            i2 = I2 - I1
            for levelc, statec in zip(
                result.E_statec_general[j, :], result.wfe_general[j, :, :]
            ):
                # pl.axhline(levelc,0.1,0.9, hold=True,color='g',ls='--')
                pl.plot(
                    xaxis[I1:I2],
                    statec[i1:i2] * config.wavefunction_scalefactor + levelc,
                    "b",
                )
                pl.plot(xaxis[I1:I2], levelc * span[I1:I2], "g", ls="--")
            for level, state in zip(
                result.E_state_general[j, :], result.wfh_general[j, :, :]
            ):
                # pl.axhline(level,0.1,0.9,color='g',ls='--')
                pl.plot(
                    xaxis[I1:I2],
                    state[i1:i2] * config.wavefunction_scalefactor + level,
                    "b",
                )
                pl.plot(xaxis[I1:I2], level * span[I1:I2], "g", ls="--")
            # pl.plot(xaxis, state**2*1e-9/dx*200.0+level,'b')
        pl.plot(xaxis, result.EF * span[0 : model.n_max], "r", ls="--")
        # pl.axhline(result.E_F,0.1,0.9,color='r',ls='--')
        pl.xlabel("Position (m)")
        pl.ylabel("Energy (meV)")
        pl.grid(True)

        fig2 = pl.figure(figsize=(10, 8))
        pl.suptitle(
            "Aestimo Results - at Equilibrium Condition ",
            fontsize=12,
        )
        pl.subplots_adjust(hspace=0.4, wspace=0.4)

        # Plotting Sigma
        # figure(0)
        pl.subplot(2, 2, 1)
        pl.plot(xaxis * 1e6, result.ro_result * 1e-6)
        pl.xlabel("x [um]")
        pl.ylabel("Total Charge Density [C/cm^3]")
        pl.title("Total Charge Density vs Position ", fontsize=10)
        pl.grid(True)

        # Plotting Efield
        # figure(1)
        pl.subplot(2, 2, 2)
        pl.plot(
            xaxis * 1e6,
            result.el_field1_result * 1e-8,
            "r",
            xaxis * 1e6,
            result.el_field2_result * 1e-8,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Electric Field  [MV/cm]")
        pl.title("Field Profile 1(red) & 2 (blue) vs Position ", fontsize=10)
        pl.grid(True)

        # Plotting Potential
        # figure(2)
        pl.subplot(2, 2, 3)
        pl.plot(xaxis * 1e6, Vt * result.fi_result)
        pl.xlabel("x [um]")
        pl.ylabel("Potential [V]")
        pl.title("Potential vs Position", fontsize=12)
        pl.legend(("fi"), loc="best", fontsize=12)
        pl.grid(True)

        pl.subplot(2, 2, 4)
        pl.plot(
            xaxis * 1e6,
            result.nf_result * 1e-6,
            "r",
            xaxis * 1e6,
            result.pf_result * 1e-6,
            "b",
        )
        pl.xlabel("x [um]")
        pl.ylabel("Electron  & Hole  Densities [1/cm^3]")
        pl.title("Electron (red)& Hole (blue) Densities vs Position ", fontsize=10)
        pl.grid(True)
        if show:
            pl.show()
    if drawFigures:
        return [fig1, fig2, None]
    else:
        return [None, None, None]
