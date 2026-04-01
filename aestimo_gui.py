
import tkinter
import tkinter.messagebox
import tkinter.filedialog
import customtkinter
import os
import sys
import json
import threading
import time
import queue
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# Solar Analysis Import
try:
    from characterize_solar import analyze_iv_curve
except ImportError:
    # Fallback if the file is missing or not in path
    def analyze_iv_curve(*args, **kwargs): return None
from aestimo import run_aestimo
import aestimo
from database import materialproperty, alloyproperty, alloyproperty4
import database
import copy

# Set theme
customtkinter.set_appearance_mode("Dark")
customtkinter.set_default_color_theme("blue")

# Setup generic logging to capture later
logging.basicConfig(filename='aestimo.log', level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

class AestimoGUI(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # --- Window Setup ---
        self.title("Aestimo 1D Semiconductor Solver - Professional Edition")
        self.geometry(f"{1400}x{900}")

        # Grid configuration
        self.grid_columnconfigure(1, weight=1) # Main content area
        self.grid_rowconfigure(0, weight=1) 
        
        # Project State
        self.is_simulating = False
        self.solar_study_active = False # Flag to distinguish full solar characterization
        self.log_queue = queue.Queue()
        self.project_name = "untitled_project"
        self.examples_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "examples"))

        # --- Sidebar (Quick Access & Info) ---
        self.sidebar_frame = customtkinter.CTkFrame(self, width=200, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(6, weight=1)

        self.logo_label = customtkinter.CTkLabel(self.sidebar_frame, text="AESTIMO 1D", font=customtkinter.CTkFont(size=24, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))
        
        self.status_label = customtkinter.CTkLabel(self.sidebar_frame, text="Status: Ready", anchor="w", text_color="gray")
        self.status_label.grid(row=1, column=0, padx=20, pady=(0, 20))

        self.run_button = customtkinter.CTkButton(self.sidebar_frame, text="RUN SIMULATION", 
                                                 command=self.start_simulation_thread, 
                                                 height=50, 
                                                 font=customtkinter.CTkFont(size=14, weight="bold"),
                                                 fg_color="#1F6AA5", hover_color="#144870")
        self.run_button.grid(row=2, column=0, padx=20, pady=10)

        self.val_button = customtkinter.CTkButton(self.sidebar_frame, text="RUN VALIDATION", 
                                                 command=self.run_experimental_validation, 
                                                 height=50, 
                                                 font=customtkinter.CTkFont(size=14, weight="bold"),
                                                 fg_color="#16A085", hover_color="#117A65")
        self.val_button.grid(row=3, column=0, padx=20, pady=10)
        
        self.save_btn = customtkinter.CTkButton(self.sidebar_frame, text="Save Project", command=self.save_project, fg_color="green", hover_color="darkgreen")
        self.save_btn.grid(row=4, column=0, padx=20, pady=5)
        
        self.load_btn = customtkinter.CTkButton(self.sidebar_frame, text="Load Project", command=self.load_project, fg_color="#D35400", hover_color="#A04000")
        self.load_btn.grid(row=5, column=0, padx=20, pady=5)
        
        # Progress Bar
        self.progress_bar = customtkinter.CTkProgressBar(self.sidebar_frame, orientation="horizontal")
        self.progress_bar.grid(row=7, column=0, padx=20, pady=20)
        self.progress_bar.set(0)

        # --- Main Content Area (Tabs) ---
        self.tabview = customtkinter.CTkTabview(self, width=800)
        self.tabview.grid(row=0, column=1, padx=20, pady=10, sticky="nsew")
        
        self.tab_structure = self.tabview.add("Structure")
        self.tab_physics = self.tabview.add("Physics & Environment")
        self.tab_solver = self.tabview.add("Solver & Grid")
        self.tab_results = self.tabview.add("Results")
        self.tab_console = self.tabview.add("Console")
        self.tab_validation = self.tabview.add("Validation")

        self.setup_structure_tab()
        self.setup_physics_tab()
        self.setup_solver_tab()
        self.setup_results_tab()
        self.setup_console_tab()
        self.setup_validation_tab()

        # --- Database Management ---
        # Keep a deep copy of original dictionaries to allow reset
        self.default_material_property = copy.deepcopy(database.materialproperty)
        self.default_alloy_property = copy.deepcopy(database.alloyproperty)
        self.default_alloy_property4 = copy.deepcopy(database.alloyproperty4)
        
        self.tab_database = self.tabview.add("Database")
        self.setup_database_tab()
        
        # Auto-save default project to examples folder
        self.auto_save_default_project()
        
        # Start Log Polling
        self.after(500, self.poll_log_file)

    def setup_structure_tab(self):
        """Setup the Layer Editor"""
        self.tab_structure.grid_columnconfigure(0, weight=1)
        self.tab_structure.grid_rowconfigure(0, weight=1)
        
        # Tools row (Add/Clear)
        tools_frame = customtkinter.CTkFrame(self.tab_structure, height=40)
        tools_frame.grid(row=1, column=0, padx=10, pady=(0,10), sticky="ew")
        
        self.add_layer_btn = customtkinter.CTkButton(tools_frame, text="+ Add Layer", command=self.add_layer, width=100)
        self.add_layer_btn.pack(side="left", padx=10, pady=5)
        
        self.add_subtrate_btn = customtkinter.CTkButton(tools_frame, text="Add Substrate/Buffer", 
                                                       command=lambda: self.add_layer(thickness=500, material="GaAs", type="barrier"), 
                                                       fg_color="gray", width=150)
        self.add_subtrate_btn.pack(side="left", padx=10, pady=5)

        self.preview_btn = customtkinter.CTkButton(tools_frame, text="Preview Bands", 
                                                   command=self.preview_band_structure, 
                                                   fg_color="#D35400", hover_color="#A04000", width=120)
        self.preview_btn.pack(side="left", padx=10, pady=5)

        self.browse_mat_btn = customtkinter.CTkButton(tools_frame, text="Browse Materials", 
                                                      command=self.open_material_browser, 
                                                      fg_color="#16A085", hover_color="#117A65", width=140)
        self.browse_mat_btn.pack(side="left", padx=10, pady=5)

        self.clear_btn = customtkinter.CTkButton(tools_frame, text="Clear All", command=self.clear_layers, fg_color="#D94848", hover_color="#A02020", width=80)
        self.clear_btn.pack(side="right", padx=10, pady=5)

        # Scrollable Area for Layers
        self.layers_frame = customtkinter.CTkScrollableFrame(self.tab_structure, label_text="Heterostructure Definition (Growth Direction ↓)")
        self.layers_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.layers_frame.grid_columnconfigure(0, weight=1)
        
        self.layer_widgets = []
        
        # Initial Demo Structure (InGaN Solar Cell)
        # Layer 1: p-InGaN (x=0.57) barrier
        self.add_layer(material="InGaN", thickness=40.0, type="barrier", mole=0.57, doping=1.0e16, doping_type="p")
        # Layer 2: n-InGaN (x=0.57) barrier (absorber)
        self.add_layer(material="InGaN", thickness=60.0, type="barrier", mole=0.57, doping=2.0e17, doping_type="n")
        
    def setup_physics_tab(self):
        """Environment and Physics Models"""
        # Create a scrollable frame to hold all the settings
        self.physics_scroll_frame = customtkinter.CTkScrollableFrame(self.tab_physics, fg_color="transparent")
        self.physics_scroll_frame.pack(fill="both", expand=True)
        
        frame = self.physics_scroll_frame
        frame.grid_columnconfigure(0, weight=1)
        
        # Group: Environment
        env_frame = customtkinter.CTkFrame(frame)
        env_frame.pack(fill="x", padx=10, pady=10)
        customtkinter.CTkLabel(env_frame, text="Environment Variables", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        self.create_input_row(env_frame, "Temperature (K):", "temp_entry", "300.0")
        self.create_input_row(env_frame, "Applied Electric Field (kV/cm):", "field_entry", "0.0") 
        
        # Group: Boundary Conditions
        bc_frame = customtkinter.CTkFrame(frame)
        bc_frame.pack(fill="x", padx=10, pady=10)
        customtkinter.CTkLabel(bc_frame, text="Boundary Conditions (Potential)", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        self.create_input_row(bc_frame, "Left Boundary (V):", "bc_left_entry", "0.0")
        self.create_input_row(bc_frame, "Right Boundary (V):", "bc_right_entry", "0.6")

        # Group: External Bias Sweep (Optional)
        sweep_frame = customtkinter.CTkFrame(frame)
        sweep_frame.pack(fill="x", padx=10, pady=10)
        customtkinter.CTkLabel(sweep_frame, text="Voltage Sweep (for I-V)", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        self.create_input_row(sweep_frame, "V Min (V):", "vmin_entry", "0.0")
        self.create_input_row(sweep_frame, "V Max (V):", "vmax_entry", "1.6")
        self.create_input_row(sweep_frame, "Step (V):", "vstep_entry", "0.05")
        
        # Group: Experimental Data
        exp_frame = customtkinter.CTkFrame(frame)
        exp_frame.pack(fill="x", padx=10, pady=10)
        customtkinter.CTkLabel(exp_frame, text="Experimental Validation", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        row = customtkinter.CTkFrame(exp_frame, fg_color="transparent")
        row.pack(fill="x", padx=5, pady=2)
        customtkinter.CTkLabel(row, text="Exp. I-V File:", width=150, anchor="w").pack(side="left")
        self.exp_file_entry = customtkinter.CTkEntry(row)
        self.exp_file_entry.pack(side="left", fill="x", expand=True)
        self.exp_file_entry.insert(0, "examples/experimental_data/si_pn_experimental_iv.csv")
        customtkinter.CTkButton(row, text="Browse", width=60, command=self.browse_exp_file).pack(side="left", padx=5)
        
        self.create_input_row(exp_frame, "Device Area (cm²):", "area_entry", "5e-4")
        self.create_input_row(exp_frame, "Series Res. (Ω):", "rs_entry", "11.7")
        
        # Rs Mode Selector
        rs_mode_row = customtkinter.CTkFrame(exp_frame, fg_color="transparent")
        rs_mode_row.pack(fill="x", padx=5, pady=2)
        customtkinter.CTkLabel(rs_mode_row, text="Rs Mode:", width=150, anchor="w").pack(side="left")
        self.rs_mode_combo = customtkinter.CTkComboBox(rs_mode_row, values=["External (Fast)", "Internal (Self-Consistent)"])
        self.rs_mode_combo.pack(side="left", fill="x", expand=True)
        self.rs_mode_combo.set("External (Fast)")
        
        self.create_input_row(exp_frame, "Shunt Res. (Ω):", "rsh_entry", "5000")

        # Group: Junction Physics
        junc_frame = customtkinter.CTkFrame(frame)
        junc_frame.pack(fill="x", padx=10, pady=10)
        customtkinter.CTkLabel(junc_frame, text="Junction Physics", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        self.graded_junc_var = tkinter.BooleanVar(value=False)
        self.graded_junc_cb = customtkinter.CTkCheckBox(junc_frame, text="Enable Graded (Diffused) Junction", variable=self.graded_junc_var)
        self.graded_junc_cb.pack(anchor="w", padx=20, pady=5)
        
        self.create_input_row(junc_frame, "Diffusion Length (nm):", "diffusion_len_entry", "10.0")

        # Group: Device Type & Solar Parameters
        device_frame = customtkinter.CTkFrame(frame)
        device_frame.pack(fill="x", padx=10, pady=10)
        customtkinter.CTkLabel(device_frame, text="Device Type (Solar/Detector)", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        self.device_type_combo = customtkinter.CTkComboBox(device_frame, values=["Generic Diode / LED", "Solar Cell / Photodetector"], command=self.on_device_type_change)
        self.device_type_combo.pack(fill="x", padx=10, pady=5)
        self.device_type_combo.set("Solar Cell / Photodetector")
        
        # Solar Cell / Photodetector Options
        self.solar_options_frame = customtkinter.CTkFrame(frame)
        # We'll pack it in on_device_type_change
        
        customtkinter.CTkLabel(self.solar_options_frame, text="Solar Parameters", font=customtkinter.CTkFont(weight="bold")).pack(pady=5)
        self.create_input_row(self.solar_options_frame, "Optical Gen Rate (cm⁻³s⁻¹):", "g_opt_entry", "1e21")
        
        self.solar_study_btn = customtkinter.CTkButton(self.solar_options_frame, text="RUN FULL SOLAR STUDY", 
                                                      command=self.start_solar_study_thread,
                                                      fg_color="#D35400", hover_color="#A04000")
        self.solar_study_btn.pack(pady=10, padx=20)
        
        # Initial state: hidden if default is not solar
        self.on_device_type_change(self.device_type_combo.get())

    def setup_solver_tab(self):
        """Computational Settings"""
        # Create a scrollable frame to hold all the settings
        self.solver_scroll_frame = customtkinter.CTkScrollableFrame(self.tab_solver, fg_color="transparent")
        self.solver_scroll_frame.pack(fill="both", expand=True)
        
        frame = self.solver_scroll_frame
        
        schema_map = {
            "0: Schrodinger": 0,
            "1: Schrodinger + Non-parabolicity": 1,
            "2: Schrodinger-Poisson": 2,
            "3: Schrodinger-Poisson + Non-parabolicity": 3,
            "4: Schrodinger-Exchange": 4,
            "5: Schrodinger-Poisson + Exchange": 5,
            "6: SP + Exchange + Non-parabolicity": 6,
            "7: SP-Drift Diffusion (Sequential)": 7,
            "8: SP-Drift Diffusion (Simultaneous)": 8,
            "9: SP-DD (Gummel-Newton)": 9
        }
        self.schema_map_rev = {v: k for k, v in schema_map.items()}
        
        # Solver Selection
        customtkinter.CTkLabel(frame, text="Solver Model:").pack(anchor="w", padx=20, pady=(20, 5))
        self.solver_combo = customtkinter.CTkComboBox(frame, values=list(schema_map.keys()), width=300)
        self.solver_combo.pack(anchor="w", padx=20, pady=5)
        self.solver_combo.set("7: SP-Drift Diffusion (Sequential)")

        # Grid Settings
        grid_frame = customtkinter.CTkFrame(frame)
        grid_frame.pack(fill="x", padx=20, pady=20)
        customtkinter.CTkLabel(grid_frame, text="Grid Resolution", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        self.create_input_row(grid_frame, "Grid Step (nm):", "grid_step_entry", "1.0")
        self.create_input_row(grid_frame, "Max Grid Points:", "max_points_entry", "200000")

        # Material System
        mat_frame = customtkinter.CTkFrame(frame)
        mat_frame.pack(fill="x", padx=20, pady=10)
        customtkinter.CTkLabel(mat_frame, text="Material System", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        
        self.mat_system_combo = customtkinter.CTkComboBox(mat_frame, values=["Zincblende", "Wurtzite"])
        self.mat_system_combo.pack(anchor="w", padx=10, pady=10)
        self.mat_system_combo.set("Wurtzite")
        
        # Subbands
        state_frame = customtkinter.CTkFrame(frame)
        state_frame.pack(fill="x", padx=20, pady=10)
        customtkinter.CTkLabel(state_frame, text="Quantum States", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        self.create_input_row(state_frame, "Electron Subbands:", "sub_e_entry", "5")
        self.create_input_row(state_frame, "Hole Subbands:", "sub_h_entry", "5")
        
        # Quantum Regions (Advanced)
        qr_frame = customtkinter.CTkFrame(frame)
        qr_frame.pack(fill="x", padx=20, pady=10)
        customtkinter.CTkLabel(qr_frame, text="Quantum Regions (Advanced)", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)

        self.quantum_regions_var = customtkinter.BooleanVar(value=False)
        self.qr_cb = customtkinter.CTkCheckBox(qr_frame, text="Enable Quantum Corrected DD", variable=self.quantum_regions_var)
        self.qr_cb.pack(fill="x", padx=10, pady=5)
        
        row_frame = customtkinter.CTkFrame(qr_frame, fg_color="transparent")
        row_frame.pack(fill="x", padx=5, pady=2)
        customtkinter.CTkLabel(row_frame, text="QR Boundary (nm, JSON list):", width=150, anchor="w").pack(side="left")
        self.qr_boundary_entry = customtkinter.CTkEntry(row_frame)
        self.qr_boundary_entry.pack(side="left", fill="x", expand=True)
        self.qr_boundary_entry.insert(0, "[[0.0, 0.0]]")
        
        # TAT Field
        tat_frame = customtkinter.CTkFrame(frame)
        tat_frame.pack(fill="x", padx=20, pady=10)
        customtkinter.CTkLabel(tat_frame, text="Trap-Assisted Tunneling (TAT)", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        self.create_input_row(tat_frame, "TAT Field (V/m, 1e10=off):", "tat_field_entry", "5e6")
    def on_device_type_change(self, choice):
        if hasattr(self, 'solar_options_frame'):
            if choice == "Solar Cell / Photodetector":
                self.solar_options_frame.pack(fill="x", padx=10, pady=10)
            else:
                self.solar_options_frame.pack_forget()
    def setup_results_tab(self):
        """Results Display"""
        self.tab_results.grid_columnconfigure(0, weight=1)
        self.tab_results.grid_rowconfigure(0, weight=1)
        
        # Main Container (Vertical Split if solar metrics present)
        self.results_main_frame = customtkinter.CTkFrame(self.tab_results)
        self.results_main_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        self.results_main_frame.grid_columnconfigure(0, weight=3) # Plot area
        self.results_main_frame.grid_columnconfigure(1, weight=1) # Metrics area (right)
        self.results_main_frame.grid_rowconfigure(0, weight=1)

        # A container for the matplotlib canvas
        self.results_container = customtkinter.CTkFrame(self.results_main_frame)
        self.results_container.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        
        # A container for Solar Metrics
        self.solar_metrics_frame = customtkinter.CTkScrollableFrame(self.results_main_frame, label_text="Solar Cell Metrics")
        # Hidden by default, shown only when relevant
        # self.solar_metrics_frame.grid(row=0, column=1, sticky="nsew", padx=2, pady=2)
        
        # Placeholder
        customtkinter.CTkLabel(self.results_container, text="Run a simulation to see results here.").pack(expand=True)

    def setup_console_tab(self):
        self.tab_console.grid_columnconfigure(0, weight=1)
        self.tab_console.grid_rowconfigure(0, weight=1)
        
        self.console_text = customtkinter.CTkTextbox(self.tab_console, font=("Consolas", 12))
        self.console_text.insert("0.0", "--- Aestimo Console ---\n")
        
    def setup_validation_tab(self):
        """Experimental Validation Tab Setup"""
        self.tab_validation.grid_columnconfigure(0, weight=1)
        self.tab_validation.grid_rowconfigure(1, weight=1)
        
        info_frame = customtkinter.CTkFrame(self.tab_validation)
        info_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=10)
        
        customtkinter.CTkLabel(info_frame, text="Experimental Validation Results", font=customtkinter.CTkFont(size=18, weight="bold")).pack(pady=10)
        
        # Results area (Split vertically: Report and Plot) using PanedWindow for flexibility
        # We use standard tkinter.PanedWindow because CTk doesn't have one yet.
        # Styling it to blend with dark theme.
        self.val_paned = tkinter.PanedWindow(self.tab_validation, orient=tkinter.HORIZONTAL, sashwidth=6, bg="#2B2B2B", bd=0)
        self.val_paned.grid(row=1, column=0, sticky="nsew", padx=10, pady=10)
        
        # Report (Left Pane)
        # Wrap in a frame to ensure clean resizing
        self.val_text_container = customtkinter.CTkFrame(self.val_paned)
        self.val_report_text = customtkinter.CTkTextbox(self.val_text_container, font=("Consolas", 12))
        self.val_report_text.pack(fill="both", expand=True)
        self.val_report_text.insert("0.0", "Run validation to see the report here.")
        
        self.val_paned.add(self.val_text_container, minsize=200, stretch="always")
        
        # Plot container (Right Pane)
        self.val_plot_container = customtkinter.CTkFrame(self.val_paned)
        self.val_paned.add(self.val_plot_container, minsize=400, stretch="always")
        
        customtkinter.CTkLabel(self.val_plot_container, text="Validation Plot will appear here.").pack(expand=True)

    def browse_exp_file(self):
        file_path = tkinter.filedialog.askopenfilename(
            filetypes=[("CSV Data", "*.csv"), ("Text File", "*.txt"), ("All Files", "*.*")],
            initialdir=self.examples_dir
        )
        if file_path:
            self.set_entry(self.exp_file_entry, file_path)

    def run_experimental_validation(self):
        """Runs the validation logic using current simulation output"""
        if self.is_simulating:
            tkinter.messagebox.showwarning("Busy", "Please wait for simulation to finish.")
            return
            
        exp_file = self.exp_file_entry.get()
        if not os.path.exists(exp_file):
            # Try relative path
            script_dir = os.path.dirname(os.path.abspath(__file__))
            rel_path = os.path.join(script_dir, exp_file)
            if os.path.exists(rel_path):
                exp_file = rel_path
            else:
                tkinter.messagebox.showerror("Error", f"Experimental file not found: {exp_file}")
                return
        
        try:
            from aeslibs.experimental_validation import (
                load_experimental_data, load_current_from_avcurr,
                compute_error_metrics, plot_iv_comparison, generate_validation_report,
                apply_parasitic_resistances, calculate_ideality_factor
            )
            
            # 1. Load Experimental Data
            exp_voltage, exp_current = load_experimental_data(exp_file)
            
            # 2. Extract Simulation Results (Accurate internal current)
            output_dir = os.path.join(self.examples_dir, self.project_name + "_output")
            if not os.path.exists(output_dir):
                tkinter.messagebox.showerror("Error", f"Simulation output not found at: {output_dir}\nPlease run a simulation first.")
                return
            
            area_cm2 = float(self.area_entry.get())
            calc_v_int, calc_i_int = load_current_from_avcurr(output_dir, device_area_cm2=area_cm2)
            
            # 3. Apply Parasitic Resistances (Only if using External Mode)
            rs = float(self.rs_entry.get())
            rsh = float(self.rsh_entry.get())
            rs_mode = self.rs_mode_combo.get()
            
            if rs_mode == "Internal (Self-Consistent)":
                # Rs already applied in simulation, set applied Rs to 0 for post-processing
                # Rsh is likely not applied internally yet, so keep it? 
                # Usually Rsh is external effect, so we apply it here.
                # But apply_parasitic_resistances calculates V_ext = V_int + I*Rs
                # If V_int already includes Rs drop, then we shouldn't add it again.
                calc_v, calc_i = apply_parasitic_resistances(calc_v_int, calc_i_int, Rs=0.0, Rsh=rsh)
            else:
                 calc_v, calc_i = apply_parasitic_resistances(calc_v_int, calc_i_int, Rs=rs, Rsh=rsh)
            
            # 4. Alignment and Metrics
            # Interpolate simulation results to match experimental voltage points
            sim_current_interp = np.interp(exp_voltage, calc_v, calc_i)
            
            forward_mask = exp_voltage >= 0.1 # Focus metrics on active region
            metrics = compute_error_metrics(exp_current[forward_mask], sim_current_interp[forward_mask])
            
            # Calculate Ideality Factor for both
            v_mid_exp, n_exp = calculate_ideality_factor(exp_voltage, exp_current, temperature=float(self.temp_entry.get()))
            v_mid_sim, n_sim = calculate_ideality_factor(calc_v, calc_i, temperature=float(self.temp_entry.get()))
            
            if len(n_sim) > 0:
                metrics['avg_n'] = np.mean(n_sim)
            
            # 5. Display Results
            report = generate_validation_report(metrics)
            self.val_report_text.delete("1.0", "end")
            self.val_report_text.insert("end", report)
            
            # Display Plots (Linear, Semi-log, and Ideality Factor)
            for w in self.val_plot_container.winfo_children(): w.destroy()
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))
            
            # Linear Plot
            ax1.plot(exp_voltage, exp_current, 'ko', label='Experimental', markersize=4, alpha=0.6)
            ax1.plot(calc_v, calc_i, 'r-', label='Simulation', linewidth=2)
            ax1.set_xlabel("Voltage (V)")
            ax1.set_ylabel("Current (A)")
            ax1.set_title("I-V (Linear)")
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Semi-log Plot
            ax2.semilogy(exp_voltage, np.abs(exp_current), 'ko', label='Experimental', markersize=4, alpha=0.6)
            ax2.semilogy(calc_v, np.abs(calc_i), 'r-', label='Simulation', linewidth=2)
            ax2.set_xlabel("Voltage (V)")
            ax2.set_ylabel("Current (log A)")
            ax2.set_title("I-V (Semi-log)")
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            # --- Fix: Clip all plots to experimental range (X and Y) ---
            v_min, v_max = exp_voltage.min(), exp_voltage.max()
            padding_v = (v_max - v_min) * 0.05
            ax1.set_xlim(v_min - padding_v, v_max + padding_v)
            ax2.set_xlim(v_min - padding_v, v_max + padding_v)
            
            # Y-Axis fitting (Linear)
            i_min, i_max = exp_current.min(), exp_current.max()
            padding_i = (i_max - i_min) * 0.1
            ax1.set_ylim(i_min - padding_i, i_max + padding_i)
            
            # Y-Axis fitting (Log)
            exp_i_pos = np.abs(exp_current[exp_current > 0])
            if len(exp_i_pos) > 0:
                y_min_log = exp_i_pos.min() / 5
                y_max_log = exp_i_pos.max() * 5
                ax2.set_ylim(y_min_log, y_max_log)
            
            # Ideality Factor Plot
            if len(n_exp) > 0:
                ax3.plot(v_mid_exp, n_exp, 'ko', label='Exp. n(V)', markersize=4, alpha=0.6)
            if len(n_sim) > 0:
                ax3.plot(v_mid_sim, n_sim, 'r-', label='Sim. n(V)', linewidth=2)
            ax3.set_xlabel("Voltage (V)")
            ax3.set_ylabel("Ideality Factor n")
            ax3.set_title("Ideality Factor Analysis")
            ax3.set_ylim(0.5, 3.5)
            ax3.set_xlim(v_min - padding_v, v_max + padding_v)
            ax3.legend()
            ax3.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            canvas = FigureCanvasTkAgg(fig, master=self.val_plot_container)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)

            # Add toolbar
            toolbar = NavigationToolbar2Tk(canvas, self.val_plot_container)
            toolbar.update()
            canvas.get_tk_widget().pack(fill="both", expand=True)
            
            self.tabview.set("Validation")
            self.status_label.configure(text="Status: Validation Complete", text_color="green")
            
        except Exception as e:
            tkinter.messagebox.showerror("Validation Error", str(e))
            import traceback
            traceback.print_exc()

    # --- Helper Functions ---
    def create_input_row(self, parent, label_text, var_name, default_val):
        """Creates a labelled entry row"""
        row_frame = customtkinter.CTkFrame(parent, fg_color="transparent")
        row_frame.pack(fill="x", padx=5, pady=2)
        
        lbl = customtkinter.CTkLabel(row_frame, text=label_text, width=150, anchor="w")
        lbl.pack(side="left")
        
        entry = customtkinter.CTkEntry(row_frame)
        entry.pack(side="left", fill="x", expand=True)
        entry.insert(0, default_val)
        
        setattr(self, var_name, entry)

    def add_layer(self, material="GaAs", thickness=10, type="barrier", mole=0.0, mole_y=0.0, doping=0.0, doping_type="n"):
        row = len(self.layer_widgets)
        frame = customtkinter.CTkFrame(self.layers_frame)
        frame.pack(fill="x", padx=5, pady=5)
        
        # Material
        materials = sorted(list(materialproperty.keys()) + list(alloyproperty.keys()) + list(alloyproperty4.keys()))
        mat_option = customtkinter.CTkOptionMenu(frame, values=materials, width=100)
        mat_option.set(material)
        mat_option.grid(row=0, column=0, padx=5, pady=5)
        
        # Mole
        customtkinter.CTkLabel(frame, text="x=").grid(row=0, column=1)
        mole_entry = customtkinter.CTkEntry(frame, width=50)
        mole_entry.insert(0, str(mole))
        mole_entry.grid(row=0, column=2, padx=5)

        # Mole Y (for Quaternary)
        customtkinter.CTkLabel(frame, text="y=").grid(row=0, column=3)
        mole_y_entry = customtkinter.CTkEntry(frame, width=50)
        mole_y_entry.insert(0, str(mole_y))
        mole_y_entry.grid(row=0, column=4, padx=5)

        # Thickness
        tk_entry = customtkinter.CTkEntry(frame, width=70)
        tk_entry.insert(0, str(thickness))
        tk_entry.grid(row=0, column=5, padx=5)
        customtkinter.CTkLabel(frame, text="nm").grid(row=0, column=6)

        # Doping
        dop_entry = customtkinter.CTkEntry(frame, width=80)
        dop_entry.insert(0, f"{doping:.1e}")
        dop_entry.grid(row=0, column=7, padx=5)
        customtkinter.CTkLabel(frame, text="cm⁻³").grid(row=0, column=8)

        # Dop Type
        dop_type_opt = customtkinter.CTkOptionMenu(frame, values=["n", "p", "i"], width=60)
        dop_type_opt.set(doping_type)
        dop_type_opt.grid(row=0, column=9, padx=5)

        # Layer Type
        type_opt = customtkinter.CTkOptionMenu(frame, values=["barrier", "well"], width=90)
        type_opt.set(type)
        type_opt.grid(row=0, column=10, padx=5)
        
        # Remove
        btn = customtkinter.CTkButton(frame, text="×", width=30, fg_color="#C0392B", command=lambda f=frame: self.remove_layer(f))
        btn.grid(row=0, column=11, padx=10)
        
        self.layer_widgets.append({
            "frame": frame,
            "mat": mat_option,
            "mole": mole_entry,
            "mole_y": mole_y_entry,
            "thick": tk_entry,
            "dop": dop_entry,
            "dop_type": dop_type_opt,
            "type": type_opt
        })

    def remove_layer(self, frame):
        frame.destroy()
        self.layer_widgets = [w for w in self.layer_widgets if w["frame"].winfo_exists()]

    def clear_layers(self):
        for w in self.layer_widgets:
            w["frame"].destroy()
        self.layer_widgets = []

    # --- Save / Load Logic ---
    def get_current_configuration(self):
        """Bundles UI state into a dictionary"""
        config = {}
        
        # Layers
        layers_data = []
        for w in self.layer_widgets:
            layers_data.append({
                "material": w["mat"].get(),
                "mole": w["mole"].get(),
                "mole_y": w["mole_y"].get(),
                "thickness": w["thick"].get(),
                "doping": w["dop"].get(),
                "doping_type": w["dop_type"].get(),
                "type": w["type"].get()
            })
        config["layers"] = layers_data
        
        # Physics
        config["temp"] = self.temp_entry.get()
        config["field"] = self.field_entry.get()
        config["bc_left"] = self.bc_left_entry.get()
        config["bc_right"] = self.bc_right_entry.get()
        config["vmin"] = self.vmin_entry.get()
        config["vmax"] = self.vmax_entry.get()
        config["vstep"] = self.vstep_entry.get()
        
        # Solver
        config["solver"] = self.solver_combo.get()
        config["grid_step"] = self.grid_step_entry.get()
        config["max_pts"] = self.max_points_entry.get()
        config["mat_sys"] = self.mat_system_combo.get()
        config["sub_e"] = self.sub_e_entry.get()
        config["sub_h"] = self.sub_h_entry.get()
        
        # Quantum Regions
        config["Quantum_Regions"] = self.quantum_regions_var.get()
        config["Quantum_Regions_boundary"] = self.qr_boundary_entry.get()
        
        # TAT
        config["tat_field"] = self.tat_field_entry.get()
        
        # Validation
        config["exp_file"] = self.exp_file_entry.get()
        config["area"] = self.area_entry.get()
        config["rs"] = self.rs_entry.get()
        config["rs_mode"] = self.rs_mode_combo.get()
        config["rsh"] = self.rsh_entry.get()
        # Solar Cell
        config["device_type"] = self.device_type_combo.get()
        config["G_optical"] = self.g_opt_entry.get()
        
        # Junction Physics
        config["graded_junc"] = self.graded_junc_var.get()
        config["diffusion_len"] = self.diffusion_len_entry.get()
        
        return config

    def load_configuration(self, config):
        """Restores UI state from dictionary"""
        # Layers
        self.clear_layers()
        for l in config.get("layers", []):
            self.add_layer(
                material=l["material"],
                thickness=float(l["thickness"]), # Ensure float cast if needed later, add_layer takes raw usually but we pass strings to entries
                type=l["type"],
                mole=float(l["mole"]),
                mole_y=float(l.get("mole_y", 0.0)),
                doping=float(l["doping"]),
                doping_type=l["doping_type"]
            )
            
        # Physics
        self.set_entry(self.temp_entry, config.get("temp", "300.0"))
        self.set_entry(self.field_entry, config.get("field", "0.0"))
        self.set_entry(self.bc_left_entry, config.get("bc_left", "0.0"))
        self.set_entry(self.bc_right_entry, config.get("bc_right", "0.0"))
        self.set_entry(self.vmin_entry, config.get("vmin", "0.0"))
        self.set_entry(self.vmax_entry, config.get("vmax", "1.0"))
        self.set_entry(self.vstep_entry, config.get("vstep", "0.05"))
        
        # Solver
        self.solver_combo.set(config.get("solver", "2: Schrodinger-Poisson"))
        self.set_entry(self.grid_step_entry, config.get("grid_step", "0.5"))
        self.set_entry(self.max_points_entry, config.get("max_pts", "200000"))
        self.mat_system_combo.set(config.get("mat_sys", "Zincblende"))
        self.set_entry(self.sub_e_entry, config.get("sub_e", "5"))
        self.set_entry(self.sub_h_entry, config.get("sub_h", "5"))
        
        # Quantum Regions
        self.quantum_regions_var.set(config.get("Quantum_Regions", False))
        qr_b = config.get("Quantum_Regions_boundary", "[[0.0, 0.0]]")
        if isinstance(qr_b, list): qr_b = json.dumps(qr_b)
        self.set_entry(self.qr_boundary_entry, qr_b)
        
        # TAT
        tat_f = config.get("tat_field", "1e10")
        if isinstance(tat_f, (int, float)): tat_f = str(tat_f)
        self.set_entry(self.tat_field_entry, tat_f)
        
        # Validation
        self.set_entry(self.exp_file_entry, config.get("exp_file", "examples/experimental_data/si_pn_experimental_iv.csv"))
        self.set_entry(self.area_entry, config.get("area", "1e-4"))
        self.set_entry(self.rs_entry, config.get("rs", "0.0"))
        self.rs_mode_combo.set(config.get("rs_mode", "External (Fast)"))
        self.set_entry(self.rsh_entry, config.get("rsh", "1e12"))
        # Solar Cell
        dev_type = config.get("device_type", "Generic Diode / LED")
        self.device_type_combo.set(dev_type)
        self.set_entry(self.g_opt_entry, config.get("G_optical", "0.0"))
        self.on_device_type_change(dev_type)
        
        # Junction Physics
        self.graded_junc_var.set(config.get("graded_junc", False))
        self.set_entry(self.diffusion_len_entry, config.get("diffusion_len", "10.0"))

    def set_entry(self, entry, value):
        entry.delete(0, "end")
        entry.insert(0, str(value))

    def save_project(self):
        file_path = tkinter.filedialog.asksaveasfilename(
            defaultextension=".json", 
            filetypes=[("Aestimo Project", "*.json")],
            initialdir=self.examples_dir,
            initialfile=f"{self.project_name}.json"
        )
        if file_path:
            try:
                config = self.get_current_configuration()
                with open(file_path, "w") as f:
                    json.dump(config, f, indent=4)
                # Extract project name from filename
                self.project_name = os.path.splitext(os.path.basename(file_path))[0]
                self.status_label.configure(text=f"Saved: {self.project_name}", text_color="green")
            except Exception as e:
                tkinter.messagebox.showerror("Save Error", str(e))

    def load_project(self, file_path=None):
        if file_path is None:
            file_path = tkinter.filedialog.askopenfilename(
                filetypes=[("Aestimo Project", "*.json")],
                initialdir=self.examples_dir
            )
        if file_path:
            try:
                with open(file_path, "r") as f:
                    config = json.load(f)
                self.load_configuration(config)
                # Extract project name from filename
                self.project_name = os.path.splitext(os.path.basename(file_path))[0]
                self.status_label.configure(text=f"Loaded: {self.project_name}", text_color="green")
            except Exception as e:
                tkinter.messagebox.showerror("Load Error", str(e))

    def auto_save_default_project(self):
        """Automatically save the default project configuration to examples folder"""
        try:
            # Create examples directory if it doesn't exist
            if not os.path.isdir(self.examples_dir):
                os.makedirs(self.examples_dir, exist_ok=True)
            
            # Save default configuration
            default_project_path = os.path.join(self.examples_dir, "untitled_project.json")
            config = self.get_current_configuration()
            with open(default_project_path, "w") as f:
                json.dump(config, f, indent=4)
        except Exception as e:
            print(f"Warning: Could not auto-save default project: {e}")

    # --- Async Simulation Logic ---
    def start_simulation_thread(self):
        if self.is_simulating:
            return
            
        self.is_simulating = True
        self.run_button.configure(state="disabled")
        self.progress_bar.configure(mode="indeterminate")
        self.progress_bar.start()
        self.status_label.configure(text="Status: Simulating...", text_color="orange")
        
        # Gather inputs in main thread
        try:
            self.last_config = self.get_current_configuration()
            
            thread = threading.Thread(target=self.run_simulation_worker, args=(self.last_config,))
            thread.daemon = True
            thread.start()
        except Exception as e:
            self.finish_simulation(success=False, error_msg=str(e))

    def start_solar_study_thread(self):
        """Starts the full solar characterization sweep (Dark, Light, Temp)"""
        print("[GUI DEBUG] start_solar_study_thread: Start button clicked.")
        if self.is_simulating:
            return
            
        self.is_simulating = True
        self.solar_study_btn.configure(state="disabled")
        self.run_button.configure(state="disabled")
        self.progress_bar.configure(mode="determinate")
        self.progress_bar.set(0)
        self.status_label.configure(text="Status: Running Solar Study...", text_color="orange")
        
        try:
            config = self.get_current_configuration()
            self.last_config = config
            
            thread = threading.Thread(target=self.run_solar_study_worker, args=(config,))
            thread.daemon = True
            thread.start()
        except Exception as e:
            self.finish_simulation(success=False, error_msg=str(e))

    def run_simulation_worker(self, config):
        try:
            # Set backend to Agg to avoid main thread loop errors when plotting in thread
            import matplotlib
            matplotlib.use('Agg', force=True)
            import matplotlib.pyplot as plt
            
            # Set output directory based on project name
            import aestimo
            output_dir = os.path.join(self.examples_dir, self.project_name + "_output")
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            aestimo.output_directory = output_dir
            
            # Reconstruct InputObject
            
            # Layers
            material_list = []
            for l in config["layers"]:
                th = float(l["thickness"])
                mat = l["material"]
                x = float(l["mole"])
                y = float(l.get("mole_y", 0.0))
                dop = float(l["doping"])
                dtype = l["doping_type"] # user visible type 'n','p'
                ltype = l["type"][0]
                
                # Logic fix for p-type and implicit handling
                 # Note: In previous turn I passed raw dtype. Aestimo usually wants number for doping.
                # If dtype is 'p', we usually make val negative or supply extra column.
                
                if dtype == "i": dtype = "n"
                
                material_list.append([th, mat, x, y, dop, dtype, ltype])

            if not material_list:
                raise ValueError("Structure is empty.")
            
            # Physics / Solver
            scheme_id = int(config["solver"].split(":")[0])
            grid_step = float(config["grid_step"])
            max_pts = int(config["max_pts"])
            sub_e = int(config["sub_e"])
            sub_h = int(config["sub_h"])
            mat_sys = config["mat_sys"]
            
            T = float(config["temp"])
            F_app = float(config["field"]) * 1e5
            
            val_vmin = float(config["vmin"])
            val_vmax = float(config["vmax"])
            val_vstep = float(config["vstep"])
            
            bc_left = float(config["bc_left"])
            bc_right = float(config["bc_right"])

            # Input Object
            class InputObject:
                T_val = T
                computation_scheme = scheme_id
                subnumber_h = sub_h
                subnumber_e = sub_e
                gridfactor = grid_step
                maxgridpoints = max_pts
                mat_type = mat_sys
                dx = grid_step * 1e-9 # m
                
                material = material_list
                
                Fapplied = F_app
                
                # TAT Field
                tat_field = float(config.get("tat_field", 1e10))
                
                # Solar Cell Optical Generation
                G_optical = float(config.get("G_optical", 0.0))
                
                vmax = val_vmax
                vmin = val_vmin
                Each_Step = val_vstep
                
                surface = np.array([bc_left, bc_right])
                
                Quantum_Regions = config.get("Quantum_Regions", False)
                qr_b_str = config.get("Quantum_Regions_boundary", "[[0.0, 0.0]]")
                try:
                    if isinstance(qr_b_str, str):
                        Quantum_Regions_boundary = np.array(json.loads(qr_b_str))
                    else:
                        Quantum_Regions_boundary = np.array(qr_b_str)
                except:
                    Quantum_Regions_boundary = np.zeros((1,2))
                
                dop_profile = None 
                
                # Series Resistance Handling
                # If Internal mode, pass Rs to model
                # If External mode, pass 0.0 to model (applied later in validation)
                rs_val = float(config.get("rs", 0.0))
                rs_mode = config.get("rs_mode", "External (Fast)")
                
                if rs_mode == "Internal (Self-Consistent)":
                    Rs = rs_val
                else:
                    Rs = 0.0
                
                # Also pass device area for current density calc
                device_area_m2 = float(config.get("area", 1e-4)) * 1e-4
                
                # Optimized physical parameters
                photovoltaic_mode = True
                enable_polarization = config.get("polarization", True)
                work_function_left = float(config.get("bc_left", 7.0))
                work_function_right = float(config.get("bc_right", 4.0))
                surface_recomb = (0, 0)

            # Fix annoying class attribute name mismatch if any (T vs T_val)
            InputObject.T = InputObject.T_val # just in case

            # Doping Profile
            tot_thick = sum(row[0] for row in material_list) * 1e-9
            dx_m = grid_step * 1e-9
            n_max = int(tot_thick / dx_m)
            
            is_graded = config.get("graded_junc", False)
            diff_len = float(config.get("diffusion_len", "10.0"))
            
            if is_graded and len(material_list) == 2:
                from scipy.special import erf
                # Calculate graded profile centered at first interface
                junc_pos_nm = material_list[0][0]
                xaxis_nm = np.linspace(0, sum(m[0] for m in material_list), n_max)
                
                d0 = material_list[0][4]
                if material_list[0][5] == 'p': d0 = -d0
                d1 = material_list[1][4]
                if material_list[1][5] == 'p': d1 = -d1
                
                # Scale to m^-3
                d0 *= 1e6
                d1 *= 1e6
                
                dop_arr = (d1 + d0)/2.0 + (d1 - d0)/2.0 * erf((xaxis_nm - junc_pos_nm) / diff_len)
                # Clean up layer doping to avoid double counting in aestimo.py
                # We keep type 'n' or 'p' so Aestimo knows the base type for mobility etc. if needed
                # but set magnitude to 0.1 (effectively 0, but avoiding possible div-by-zero if any)
                for i in range(len(material_list)):
                    material_list[i][4] = 0.0
            else:
                dop_arr = np.zeros(n_max)
                curr = 0
                for row in material_list:
                    th_m = row[0] * 1e-9
                    val = row[4]
                    dtype = row[5]
                    if dtype == 'p': val = -val
                    
                    steps = int(th_m / dx_m)
                    end = min(curr + steps, n_max)
                    dop_arr[curr:end] = val * 1e6
                    curr = end
            
            InputObject.dop_profile = dop_arr
            InputObject.__file__ = os.path.abspath("PRO_GUI_SIM_ASYNC.py")

            # Run
            # We must be careful plotting in a thread. 
            # run_aestimo generates figures. Matplotlib behaves badly in threads sometimes.
            # Best practice: Generate figures, but don't show() them. 
            # We are using FigureCanvasTkAgg in main thread later.
            
            # NOTE: run_aestimo uses plt calls internally. This is risky in non-main thread.
            # However, with Agg backend or careful handling it might work.
            # Ideally we'd refactor aestimo to return data objects, but that's huge work.
            # We try standard run.
            
            input_obj, model, result, figures = run_aestimo(InputObject, drawFigures=True, show=False)
            
            # Post results back to main thread
            self.after(0, self.finish_simulation, True, None, figures)

        except Exception as e:
            self.after(0, self.finish_simulation, False, str(e), None)

    def run_solar_study_worker(self, config):
        try:
            plt.switch_backend('Agg')
            
            # 1. Setup Material List and shared parameters
            material_list = []
            for layer in config["layers"]:
                material_list.append([
                    float(layer["thickness"]),
                    layer["material"],
                    float(layer["mole"]),
                    float(layer["mole_y"]),
                    float(layer["doping"]),
                    layer["doping_type"],
                    layer["type"]
                ])
            
            grid_step = float(config.get("grid_step", 1.0))
            max_pts = int(config.get("max_pts", 1000))
            device_area = float(config.get("area", 1.0)) # cm^2

            # Define helper Inner class for InputObject
            class StudyInputObject:
                def __init__(self, cfg_dict, label, T_override=None, G_override=None):
                    self.T_val = T_override if T_override is not None else float(cfg_dict.get("temp", 300))
                    self.T = self.T_val
                    self.computation_scheme = 7 # Match physical baseline Solver 7 (Sequential)
                    self.comp_scheme = 7
                    self.subnumber_h = int(cfg_dict.get("sub_h", 5))
                    self.subnumber_e = int(cfg_dict.get("sub_e", 5))
                    self.gridfactor = float(cfg_dict.get("grid_step", 1.0))
                    self.dx = self.gridfactor * 1e-9
                    self.maxgridpoints = int(cfg_dict.get("max_pts", 1000))
                    self.mat_type = cfg_dict.get("mat_system", "Wurtzite")
                    self.material = material_list
                    # For Solar Study, we MUST sweep forward bias to see the power quadrant.
                    # 1.1V is usually sufficient for InGaN Voc and more stable for Solver 7.
                    self.vmin = 0.0
                    self.vmax = 1.1
                    self.Each_Step = 0.05  # Finer steps for better gradient and stability
                    self.G_optical = G_override if G_override is not None else float(cfg_dict.get("G_optical", 0.0))
                    self.tat_field = float(cfg_dict.get("tat_field", 1e10))
                    self.device_area_m2 = device_area * 1e-4
                    self.device_area = device_area
                    self.__file__ = f"STUDY_{label}.py"
                    self.inputfilename = f"STUDY_{label}"
                    
                    # Doping Profile (Manual build for p-n)
                    tot_m = sum(row[0] for row in self.material) * 1e-9
                    n_max = int(tot_m / self.dx)
                    dop_arr = np.zeros(n_max)
                    curr = 0
                    for row in self.material:
                        th_m = row[0] * 1e-9
                        val = row[4]
                        if row[5] == 'p': val = -val
                        steps = int(th_m / self.dx)
                        end = min(curr + steps, n_max)
                        dop_arr[curr:end] = val * 1e6
                        curr = end
                    self.dop_profile = dop_arr
                    
                    # Missing physical parameters from optimized baseline
                    self.photovoltaic_mode = True
                    self.enable_polarization = True
                    self.work_function_left = 7.0
                    self.work_function_right = 4.0
                    self.surface_recomb = (0, 0)
                    self.Quantum_Regions = False
                    self.Quantum_Regions_boundary = np.zeros((1, 2))
                    self.Rs = 0.0

            import config as sim_config
            sim_config.Drift_Diffusion_out = True
            
            # --- EXECUTION ---
            total_tasks = 6 # Dark, Light, 4 Temps
            
            # Task 1: Dark
            self.after(0, lambda: self.status_label.configure(text="Solar Study: Running Dark I-V..."))
            dark_in = StudyInputObject(config, "Dark", G_override=0.0)
            aestimo.output_directory = os.path.join(os.getcwd(), "STUDY_Dark_output")
            _, _, res_dark, _ = run_aestimo(dark_in, drawFigures=False, show=False)
            dark_data = np.loadtxt(os.path.join(aestimo.output_directory, "av_curr.dat"))
            dark_metrics = analyze_iv_curve(dark_data[:,0], dark_data[:,1], area_cm2=device_area)
            self.after(0, lambda: self.progress_bar.set(1/total_tasks))

            # Task 2: Light (Standard)
            self.after(0, lambda: self.status_label.configure(text="Solar Study: Running Illuminated I-V..."))
            light_in = StudyInputObject(config, "Light")
            aestimo.output_directory = os.path.join(os.getcwd(), "STUDY_Light_output")
            _, _, res_light, standard_figures = run_aestimo(light_in, drawFigures=True, show=False)
            light_data = np.loadtxt(os.path.join(aestimo.output_directory, "av_curr.dat"))
            light_metrics = analyze_iv_curve(light_data[:,0], light_data[:,1], area_cm2=device_area)
            self.after(0, lambda: self.progress_bar.set(2/total_tasks))

            # Task 3: Temperature Sweep (200K to 500K for detailed sensor/cell characterization)
            temp_results = []
            temps = np.arange(200, 525, 25)
            for i, T in enumerate(temps):
                self.after(0, lambda t=T: self.status_label.configure(text=f"Solar Study: Running {t}K..."))
                t_in = StudyInputObject(config, f"Temp_{T}", T_override=float(T))
                aestimo.output_directory = os.path.join(os.getcwd(), f"STUDY_Temp_{T}_output")
                run_aestimo(t_in, drawFigures=False, show=False)
                t_data = np.loadtxt(os.path.join(aestimo.output_directory, "av_curr.dat"))
                m = analyze_iv_curve(t_data[:,0], t_data[:,1], area_cm2=device_area)
                temp_results.append((T, m))
                self.after(0, lambda idx=i: self.progress_bar.set((3+idx)/total_tasks))

            # Prepare result data
            study_data = {
                "dark_metrics": dark_metrics,
                "light_metrics": light_metrics,
                "temp_results": temp_results,
                "config": config,
                "standard_figures": standard_figures
            }
            
            self.after(0, self.finish_simulation, True, None, study_data)

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.after(0, self.finish_simulation, False, str(e), None)
        finally:
            self.after(0, lambda: self.solar_study_btn.configure(state="normal"))

    def finish_simulation(self, success, error_msg=None, figures=None):
        self.is_simulating = False
        self.progress_bar.stop()
        self.run_button.configure(state="normal")
        
        # Check if this is a Solar Study (passed as data dict)
        is_study_data = isinstance(figures, dict) and "light_metrics" in figures
        
        if success:
            self.status_label.configure(text="Status: Simulation Complete", text_color="green")
            
            final_figures = []
            study_figures = []
            if is_study_data:
                print(f"[GUI DEBUG] finish_simulation: Detected FULL SOLAR STUDY data.")
                self.solar_study_active = True
                study_figures = self.generate_study_figures(figures)
                print(f"[GUI DEBUG] finish_simulation: Generated {len(study_figures)} specialized characterization plots.")
                # Use standard figures from the study if provided
                if "standard_figures" in figures:
                    final_figures = figures["standard_figures"]
                    print(f"[GUI DEBUG] finish_simulation: Found {len(final_figures)} standard physics plots to show.")
            else:
                print(f"[GUI DEBUG] finish_simulation: Detected STANDARD simulation data.")
                self.solar_study_active = False
                final_figures = figures if isinstance(figures, list) else []

            # Check for Solar Analysis (Standard Run Only)
            is_solar = hasattr(self, 'last_config') and self.last_config.get("device_type") == "Solar Cell / Photodetector"
            if is_solar and not is_study_data:
                print("DEBUG: Standard Solar Run. Calculating metrics and appending P-V plot.")
                metrics = self.update_solar_metrics()
                if metrics:
                    try:
                        from matplotlib.figure import Figure
                        fig_pv = Figure(figsize=(6, 4))
                        ax_pv = fig_pv.add_subplot(1, 1, 1)
                        ax_pv.plot(metrics['v'], metrics['p'], 'g-', label='Power Density')
                        ax_pv.axvline(metrics['vmpp'], color='orange', ls=':', label=f"MPP: {metrics['pmpp']:.2f} mW/cm²")
                        ax_pv.axhline(0, color='black', lw=0.5)
                        ax_pv.set_xlabel("Voltage (V)")
                        ax_pv.set_ylabel("Power Density (mW/cm²)")
                        ax_pv.set_title("P-V Characteristic")
                        ax_pv.legend()
                        ax_pv.grid(True, alpha=0.3)
                        final_figures.append(fig_pv)
                    except Exception as e:
                        print(f"Error generating solar plot: {e}")
            elif not is_solar:
                self.solar_metrics_frame.grid_forget()

            print(f"DEBUG: Calling display_figures with {len(final_figures)} std and {len(study_figures)} study figs.")
            self.display_figures(final_figures, study_figures=study_figures)
            
            # Update metrics sidebar with the standard 300K characterization
            metrics_to_show = figures.get('light_metrics') if is_study_data else None
            self.after(0, lambda m=metrics_to_show: self.update_solar_metrics_from_data(m))
            
            # Reset flag AFTER displaying
            self.solar_study_active = False
            tkinter.messagebox.showinfo("Success", "Simulation completed successfully.")
        else:
            self.status_label.configure(text="Status: Error", text_color="red")
            self.solar_study_active = False
            tkinter.messagebox.showerror("Error", f"Simulation Failed:\n{error_msg}")

    def generate_study_figures(self, data):
        """Generates the 4 study figures in the main thread using Figure class"""
        from matplotlib.figure import Figure
        
        dark_m = data["dark_metrics"]
        light_m = data["light_metrics"]
        temp_res = data["temp_results"]
        
        # Fig 1: J-V
        fig_jv = Figure(figsize=(6, 4))
        ax1 = fig_jv.add_subplot(1, 1, 1)
        ax1.plot(dark_m['v'], dark_m['j'], 'k--', label='Dark')
        ax1.plot(light_m['v'], light_m['j'], 'r-', label='Light')
        ax1.set_xlabel('Voltage (V)')
        ax1.set_ylabel('Current Density (mA/cm²)')
        ax1.set_title('Dark vs Illuminated J-V')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Fig 2: P-V
        fig_pv = Figure(figsize=(6, 4))
        ax2 = fig_pv.add_subplot(1, 1, 1)
        ax2.plot(light_m['v'], light_m['p'], 'g-', label='Power')
        ax2.axvline(light_m['vmpp'], color='orange', ls=':', label=f"MPP: {light_m['pmpp']:.2f} mW/cm²")
        ax2.set_xlabel('Voltage (V)')
        ax2.set_ylabel('Power Density (mW/cm²)')
        ax2.set_title('P-V Characteristic')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # Fig 3: Voc vs T
        fig_voct = Figure(figsize=(6, 4))
        ax3 = fig_voct.add_subplot(1, 1, 1)
        ax3.plot([r[0] for r in temp_res], [r[1]['voc'] for r in temp_res], 'bo-')
        ax3.set_xlabel('Temperature (K)')
        ax3.set_ylabel('Voc (V)')
        ax3.set_title('Open-Circuit Voltage vs Temperature')
        ax3.grid(True, alpha=0.3)

        # Fig 4: Eff vs T
        fig_efft = Figure(figsize=(6, 4))
        ax4 = fig_efft.add_subplot(1, 1, 1)
        ax4.plot([r[0] for r in temp_res], [r[1]['eta'] for r in temp_res], 'ro-')
        ax4.set_xlabel('Temperature (K)')
        ax4.set_ylabel('Efficiency (%)')
        ax4.set_title('Efficiency vs Temperature')
        ax4.grid(True, alpha=0.3)
        
        return [fig_jv, fig_pv, fig_voct, fig_efft]

    def update_solar_metrics_from_data(self, metrics=None):
        """Update the metrics sidebar from a metrics dictionary or data file"""
        if metrics is None:
            metrics = self.update_solar_metrics()
        
        if not metrics:
            return
            
        # Show the frame
        self.solar_metrics_frame.grid(row=0, column=1, sticky="nsew", padx=2, pady=2)
        
        # Clear existing
        for widget in self.solar_metrics_frame.winfo_children():
            widget.destroy()
            
        # Populate metrics
        display_map = [
            ("Short-Circuit Current", "jsc", "mA/cm²"),
            ("Open-Circuit Voltage", "voc", "V"),
            ("Fill Factor", "ff", "%"),
            ("Efficiency", "eta", "%"),
            ("Max Power Point (Pmax)", "pmpp", "mW/cm²"),
            ("V at Pmax", "vmpp", "V"),
            ("J at Pmax", "jmpp", "mA/cm²"),
            ("Series Resistance (Rs)", "rs", "Ω·cm²"),
            ("Shunt Resistance (Rsh)", "rsh", "Ω·cm²")
        ]
        
        for label, key, unit in display_map:
            row = customtkinter.CTkFrame(self.solar_metrics_frame, fg_color="transparent")
            row.pack(fill="x", pady=2)
            customtkinter.CTkLabel(row, text=f"{label}:", anchor="w", font=customtkinter.CTkFont(size=12)).pack(side="left", padx=5)
            val = metrics.get(key, 0)
            if isinstance(val, (float, int)):
                val_str = f"{abs(val):.4g}" # Using abs for Jsc/Power for display if helpful, but Voc stays sign-sensitive?
                if key in ["jsc", "pmpp", "jmpp"]: val_str = f"{abs(val):.4g}"
            else:
                val_str = str(val)
            customtkinter.CTkLabel(row, text=f"{val_str} {unit}", font=customtkinter.CTkFont(size=12, weight="bold"), text_color="#1F6AA5").pack(side="right", padx=5)

        return metrics

    def update_solar_metrics(self):
        """Extract and display solar cell metrics"""
        # Try to find av_curr.dat
        possible_paths = [
            os.path.join(os.getcwd(), "PRO_GUI_SIM_ASYNC_output", "av_curr.dat"),
            os.path.join(os.getcwd(), "output", "av_curr.dat"),
            os.path.join(os.getcwd(), "test_solar_cell_output", "av_curr.dat") # as a stretch fallback
        ]
        
        iv_path = None
        for p in possible_paths:
            if os.path.exists(p):
                iv_path = p
                break
        
        if not iv_path:
            return # Silent fail or log
            
        try:
            data = np.loadtxt(iv_path)
            # Get area in cm^2 from config
            area_cm2 = float(self.last_config.get("area", 1.0))
            if area_cm2 <= 0: area_cm2 = 1.0
            
            metrics = analyze_iv_curve(data[:,0], data[:,1], area_cm2=area_cm2)
            if not metrics: return
            
            return self.update_solar_metrics_from_data(metrics)

        except Exception as e:
            print(f"Error in solar analysis: {e}")
            return None
            
        return metrics

    def display_figures(self, figures, study_figures=None):
        # Clear previous
        for widget in self.results_container.winfo_children():
            widget.destroy()

        if not figures and not study_figures:
            customtkinter.CTkLabel(self.results_container, text="No figures generated.").pack()
            return
            
        fig_tabview = customtkinter.CTkTabview(self.results_container)
        fig_tabview.pack(fill="both", expand=True)
        
        # Determine titles dynamically for main figures
        is_solar = hasattr(self, 'last_config') and self.last_config.get("device_type") == "Solar Cell / Photodetector"
        
        if is_solar:
             titles = ["Band Diagram", "Charge Density", "Field", "I-V / Sweep", "Solar Analysis", "Other"]
        else:
            titles = ["Band Diagram", "Charge Density", "Field", "I-V / Sweep", "Other"]
        
        # 1. Add standard figures
        for i, fig in enumerate(figures):
            if fig is None: 
                print(f"[GUI DEBUG] display_figures: Skipping None figure at index {i}")
                continue
                
            name = titles[i] if i < len(titles) else f"Figure {i+1}"
            fig_tabview.add(name)
            
            canvas = FigureCanvasTkAgg(fig, master=fig_tabview.tab(name))
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)
            
            toolbar = NavigationToolbar2Tk(canvas, fig_tabview.tab(name))
            toolbar.update()
            
            # Hover cursor
            for ax in fig.axes:
                annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="w"),
                                arrowprops=dict(arrowstyle="->"))
                annot.set_visible(False)
                ax.annot = annot
            canvas.mpl_connect("motion_notify_event", lambda event, c=canvas: self.on_plot_hover(event, c))

        # 2. Add Solar Study Cluster (Nested)
        if study_figures:
            print(f"[GUI DEBUG] display_figures: Creating 'Solar Study' nested tab with {len(study_figures)} graphs.")
            study_tab_name = "Solar Study"
            fig_tabview.add(study_tab_name)
            
            study_tabview = customtkinter.CTkTabview(fig_tabview.tab(study_tab_name))
            study_tabview.pack(fill="both", expand=True)
            
            study_titles = ["J-V Comparison", "P-V Analysis", "Voc vs T", "Efficiency vs T"]
            for i, fig in enumerate(study_figures):
                 if fig is None: 
                     print(f"[GUI DEBUG] display_figures: Skipping None study figure at index {i}")
                     continue
                 name = study_titles[i] if i < len(study_titles) else f"Study {i+1}"
                 study_tabview.add(name)
                 
                 canvas = FigureCanvasTkAgg(fig, master=study_tabview.tab(name))
                 canvas.draw()
                 canvas.get_tk_widget().pack(fill="both", expand=True)
                 
                 toolbar = NavigationToolbar2Tk(canvas, study_tabview.tab(name))
                 toolbar.update()
                 
                 # Hover cursor
                 for ax in fig.axes:
                     annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
                                     bbox=dict(boxstyle="round", fc="w"),
                                     arrowprops=dict(arrowstyle="->"))
                     annot.set_visible(False)
                     ax.annot = annot
                 canvas.mpl_connect("motion_notify_event", lambda event, c=canvas: self.on_plot_hover(event, c))
        else:
            print(f"[GUI DEBUG] display_figures: No study_figures provided. Not creating Solar Study tab.")

    def on_plot_hover(self, event, canvas):
        if event.inaxes:
            ax = event.inaxes
            # Check lines
            found = False
            for line in ax.get_lines():
                cont, ind = line.contains(event)
                if cont:
                    # Update annotation
                    x, y = line.get_data()
                    # Use the first point found
                    idx = ind["ind"][0]
                    # Check if annot exists
                    if hasattr(ax, "annot"):
                        annot = ax.annot
                        annot.xy = (x[idx], y[idx])
                        text = f"x={x[idx]:.4g}\ny={y[idx]:.4g}"
                        annot.set_text(text)
                        annot.set_visible(True)
                        found = True
                        break # Only show for one line at a time
            
            if hasattr(ax, "annot"):
                # If we didn't find a line, hide the annotation
                if not found and ax.annot.get_visible():
                    ax.annot.set_visible(False)
                    canvas.draw_idle()
                elif found:
                    canvas.draw_idle()


    def preview_band_structure(self):
        """Preview band structure without running full simulation"""
        try:
            config = self.get_current_configuration()
            
            # Reconstruct Material List
            material_list = []
            for l in config["layers"]:
                th = float(l["thickness"])
                mat = l["material"]
                x = float(l["mole"])
                y = float(l.get("mole_y", 0.0))
                dop = float(l["doping"])
                dtype = l["doping_type"]
                ltype = l["type"][0]
                
                if dtype == "i": dtype = "n"
                # For preview, we pass doping but it affects potential via Poisson which isn't solved yet.
                # Structure class calculates band edges (fi_e, fi_h) purely from material properties (Eg, band offset).
                # Doping is used for 'dop' array but doesn't shift bands until Poisson is run.
                # So this shows the "flat band" diagram + built-in potentials if we calculated them, 
                # but Aestimo Structure just initializes arrays. fi_e/fi_h are band offsets relative to... vacuum? or ref?
                # Usually relative to a reference.
                
                material_list.append([th, mat, x, y, dop, dtype, ltype])

            if not material_list:
                raise ValueError("Structure is empty.")

            # Build Input Dictionary (mimics InputObject)
            grid_step = float(config["grid_step"])
            T_val = float(config["temp"])
            
            input_dict = {
                'material': material_list,
                'gridfactor': grid_step,
                'maxgridpoints': int(config["max_pts"]),
                'mat_type': config["mat_sys"],
                'T': T_val,
                'Fapplied': float(config["field"]) * 1e5,
                'subnumber_e': int(config["sub_e"]),
                'subnumber_h': int(config["sub_h"]),
                # Dummy values
                'computation_scheme': 2,
                'vmax': 0.0,
                'vmin': 0.0,
                'Each_Step': 0.1,
                'surface': [float(config["bc_left"]), float(config["bc_right"])],
                'Quantum_Regions': False,
                'dop_profile': np.array([0.0]) # Structure will resize/ignore
            }
            
            # Initialize Structure
            import aestimo
            # We must import database inside or use the one we have
            # aestimo uses 'from database import ...' so it uses the module.
            # since we update the module in-memory, it should use our updates.
            
            struct = aestimo.StructureFrom(input_dict, aestimo.database)
            
            # Plotting
            x = np.linspace(0, struct.x_max * 1e9, struct.n_max)
            # q = 1.602176634e-19
            # aestimo defines q? It does but it's not easily accessible if not in vars.
            # Using standard constant.
            q_const = 1.602176634e-19
            
            # fi_e, fi_h are in Joules. Convert to eV.
            ec = struct.fi_e / q_const
            ev = struct.fi_h / q_const
            
            # Add Applied Field tilt?
            # Structure does NOT add Fapplied field to fi_e/fi_h initially. 
            # It's added during solver loop usually.
            # We can add it manually for preview: V = -F*z
            # F is V/m. z is m.
            if struct.Fapp != 0:
                potential_tilt = -struct.Fapp * (x * 1e-9) # V
                # e*V energy shift? 
                # Energy E = -e * phi. 
                # If potential tilts down, electron energy tilts down.
                # So add potential_tilt (in eV) to bands?
                # Actually F = -grad(V). V = -F*z.
                # Band energy E = -q*V = q*F*z.
                # So we ADD F*z (in eV... wait F is V/m, z is m -> V. )
                # Ec (eV) = Ec_flat - Potential(V).
                # Potential(V) = -F*z.
                # Ec(eV) = Ec_flat + F*z.
                ec += (struct.Fapp * x * 1e-9)
                ev += (struct.Fapp * x * 1e-9)

            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(x, ec, label="Ec", color='red')
            ax.plot(x, ev, label="Ev", color='blue')
            ax.set_xlabel("Position (nm)")
            ax.set_ylabel("Energy (eV)")
            ax.set_title("Band Structure Preview (Flat Band + Field)")
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            self.display_figures([fig])
            self.tabview.set("Results")
            
        except Exception as e:
            tkinter.messagebox.showerror("Preview Error", str(e))

    def open_material_browser(self):
        """Open Material Library Browser window"""
        # Create popup window
        browser = customtkinter.CTkToplevel(self)
        browser.title("Material Library Browser")
        browser.geometry("900x600")
        
        # Left panel - Filters and Material List
        left_panel = customtkinter.CTkFrame(browser, width=400)
        left_panel.pack(side="left", fill="both", expand=False, padx=10, pady=10)
        
        # Right panel - Visualization
        right_panel = customtkinter.CTkFrame(browser)
        right_panel.pack(side="right", fill="both", expand=True, padx=10, pady=10)
        
        # Filter controls
        filter_frame = customtkinter.CTkFrame(left_panel)
        filter_frame.pack(fill="x", padx=10, pady=10)
        customtkinter.CTkLabel(filter_frame, text="Filter by Type:", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=5, pady=5)
        
        filter_var = tkinter.StringVar(value="All")
        
        filters = ["All", "Binary", "Ternary", "Quaternary"]
        for filt in filters:
            customtkinter.CTkRadioButton(filter_frame, text=filt, variable=filter_var, value=filt,
                                        command=lambda: self.update_material_list(materials_frame, filter_var.get(), plot_canvas, plot_ax)).pack(anchor="w", padx=20, pady=2)
        
        # Materials list (scrollable)
        materials_frame = customtkinter.CTkScrollableFrame(left_panel, label_text="Materials")
        materials_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Visualization area
        customtkinter.CTkLabel(right_panel, text="Bandgap vs. Lattice Constant", font=customtkinter.CTkFont(size=14, weight="bold")).pack(pady=5)
        
        # Create matplotlib figure
        fig, ax = plt.subplots(figsize=(5, 4))
        plot_canvas = FigureCanvasTkAgg(fig, master=right_panel)
        plot_canvas.get_tk_widget().pack(fill="both", expand=True)
        plot_ax = ax
        
        # Initial population
        self.update_material_list(materials_frame, "All", plot_canvas, plot_ax)
        
    def update_material_list(self, container, filter_type, plot_canvas, plot_ax):
        """Update material list based on filter"""
        # Clear existing
        for widget in container.winfo_children():
            widget.destroy()
            
        # Gather materials
        materials_data = []
        
        # Binary materials
        if filter_type in ["All", "Binary"]:
            for mat_name, props in database.materialproperty.items():
                materials_data.append({
                    "name": mat_name,
                    "type": "Binary",
                    "Eg": props.get("Eg", 0.0),
                    "a0": props.get("a0", 0.0)
                })
        
        # Ternary alloys (use x=0.5 for visualization)
        if filter_type in ["All", "Ternary"]:
            for alloy_name in database.alloyproperty.keys():
                # Estimate at x=0.5
                materials_data.append({
                    "name": alloy_name,
                    "type": "Ternary",
                    "Eg": "Variable",
                    "a0": "Variable"
                })
        
        # Quaternary alloys
        if filter_type in ["All", "Quaternary"]:
            for alloy_name in database.alloyproperty4.keys():
                materials_data.append({
                    "name": alloy_name,
                    "type": "Quaternary",
                    "Eg": "Variable",
                    "a0": "Variable"
                })
        
        # Display materials
        for mat in materials_data:
            card = customtkinter.CTkFrame(container, fg_color=("gray85", "gray25"))
            card.pack(fill="x", padx=5, pady=5)
            
            customtkinter.CTkLabel(card, text=mat["name"], font=customtkinter.CTkFont(size=12, weight="bold")).grid(row=0, column=0, sticky="w", padx=10, pady=5)
            customtkinter.CTkLabel(card, text=f"Type: {mat['type']}", font=customtkinter.CTkFont(size=10)).grid(row=1, column=0, sticky="w", padx=10)
            
            if isinstance(mat["Eg"], (int, float)):
                customtkinter.CTkLabel(card, text=f"Eg: {mat['Eg']:.3f} eV", font=customtkinter.CTkFont(size=10)).grid(row=2, column=0, sticky="w", padx=10)
                customtkinter.CTkLabel(card, text=f"a0: {mat['a0']:.3f} Å", font=customtkinter.CTkFont(size=10)).grid(row=3, column=0, sticky="w", padx=10)
            else:
                customtkinter.CTkLabel(card, text=f"Eg: {mat['Eg']}", font=customtkinter.CTkFont(size=10)).grid(row=2, column=0, sticky="w", padx=10)
                customtkinter.CTkLabel(card, text=f"a0: {mat['a0']}", font=customtkinter.CTkFont(size=10)).grid(row=3, column=0, sticky="w", padx=10)
        
        # Update plot
        plot_ax.clear()
        
        # Extract plottable data (only binary for now)
        eg_vals = []
        a0_vals = []
        names = []
        
        for mat in materials_data:
            if isinstance(mat["Eg"], (int, float)) and isinstance(mat["a0"], (int, float)):
                if mat["Eg"] > 0 and mat["a0"] > 0:  # Valid values
                    eg_vals.append(mat["Eg"])
                    a0_vals.append(mat["a0"])
                    names.append(mat["name"])
        
        if eg_vals:
            plot_ax.scatter(a0_vals, eg_vals, s=100, alpha=0.6, c='blue')
            # Add labels
            for i, name in enumerate(names):
                plot_ax.annotate(name, (a0_vals[i], eg_vals[i]), fontsize=8, alpha=0.7, 
                               xytext=(5, 5), textcoords='offset points')
            
            plot_ax.set_xlabel("Lattice Constant a0 (Å)")
            plot_ax.set_ylabel("Bandgap Eg (eV)")
            plot_ax.set_title(f"Material Properties ({filter_type})")
            plot_ax.grid(True, alpha=0.3)
        else:
            plot_ax.text(0.5, 0.5, "No plottable data\n(Alloys have variable properties)", 
                       ha='center', va='center', transform=plot_ax.transAxes, fontsize=12)
            plot_ax.set_xlabel("Lattice Constant a0 (Å)")
            plot_ax.set_ylabel("Bandgap Eg (eV)")
        
        plot_canvas.draw()




    def setup_database_tab(self):
        self.tab_database.grid_columnconfigure(0, weight=1)
        self.tab_database.grid_rowconfigure(1, weight=1)

        # Top Bar: Material Selection
        top_frame = customtkinter.CTkFrame(self.tab_database)
        top_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        
        customtkinter.CTkLabel(top_frame, text="Select Material/Alloy:").pack(side="left", padx=10)
        
        # Combine keys from all dictionaries
        self.all_materials = sorted(list(database.materialproperty.keys()) + 
                                   list(database.alloyproperty.keys()) + 
                                   list(database.alloyproperty4.keys()))
        self.db_selector = customtkinter.CTkComboBox(top_frame, values=self.all_materials, command=self.update_material_editor, width=200)
        self.db_selector.set(self.all_materials[0])
        self.db_selector.pack(side="left", padx=10)
        
        customtkinter.CTkButton(top_frame, text="+ New", command=self.add_new_material_dialog, width=80, fg_color="green").pack(side="left", padx=5)

        # Main Editor Area with Scrollbar
        self.db_editor_frame = customtkinter.CTkScrollableFrame(self.tab_database, label_text="Properties Editor")
        self.db_editor_frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
        
        self.db_entries = {} # store entry widgets reference

        # Bottom Bar: Actions
        action_frame = customtkinter.CTkFrame(self.tab_database)
        action_frame.grid(row=2, column=0, padx=10, pady=10, sticky="ew")
        
        customtkinter.CTkButton(action_frame, text="Update In-Memory", command=self.update_db_in_memory, fg_color="#1F6AA5").pack(side="left", padx=10, pady=10)
        customtkinter.CTkButton(action_frame, text="Save to File", command=self.save_database_to_file, fg_color="green").pack(side="left", padx=10, pady=10)
        customtkinter.CTkButton(action_frame, text="Load from File", command=self.load_database_from_file, fg_color="#D35400").pack(side="left", padx=10, pady=10)
        customtkinter.CTkButton(action_frame, text="Reset Defaults", command=self.reset_database_defaults, fg_color="#C0392B").pack(side="right", padx=10, pady=10)
        
        # Initial population
        self.update_material_editor(self.all_materials[0])

    def update_material_editor(self, choice):
        # Clear existing
        for widget in self.db_entries.values():
            widget.destroy() # Actually this destroys the entry, but we need to clear rows inside frame.
        
        for widget in self.db_editor_frame.winfo_children():
            widget.destroy()
            
        self.db_entries = {}
        
        # Determine source dict
        if choice in database.materialproperty:
            data = database.materialproperty[choice]
            self.current_db_type = "material"
        elif choice in database.alloyproperty:
            data = database.alloyproperty[choice]
            self.current_db_type = "alloy"
        elif choice in database.alloyproperty4:
            data = database.alloyproperty4[choice]
            self.current_db_type = "quaternary"
        else:
            return

        # Generate inputs
        for i, (key, value) in enumerate(sorted(data.items())):
            row_frame = customtkinter.CTkFrame(self.db_editor_frame, fg_color="transparent")
            row_frame.pack(fill="x", padx=5, pady=2)
            
            customtkinter.CTkLabel(row_frame, text=key, width=150, anchor="w").pack(side="left")
            entry = customtkinter.CTkEntry(row_frame)
            entry.insert(0, str(value))
            entry.pack(side="left", fill="x", expand=True)
            
            self.db_entries[key] = entry

    def update_db_in_memory(self):
        mat_name = self.db_selector.get()
        if not mat_name: return
        
        try:
            if self.current_db_type == "material":
                target_dict = database.materialproperty[mat_name]
            elif self.current_db_type == "alloy":
                target_dict = database.alloyproperty[mat_name]
            else:
                target_dict = database.alloyproperty4[mat_name]
            
            for key, entry in self.db_entries.items():
                val_str = entry.get()
                # Try to infer type: float, int, str
                # Most are floats.
                try:
                    if "." in val_str or "e" in val_str.lower():
                        val = float(val_str)
                    else:
                        try:
                            val = int(val_str)
                        except ValueError:
                            val = val_str # keep as string
                except ValueError:
                     val = val_str # fallback
                
                target_dict[key] = val
            
            self.status_label.configure(text=f"Updated: {mat_name}", text_color="green")
            tkinter.messagebox.showinfo("Update", f"Updated properties for {mat_name} in memory.")
            
        except Exception as e:
            tkinter.messagebox.showerror("Error", str(e))

    def save_database_to_file(self):
        file_path = tkinter.filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("Aestimo Database", "*.json")],
            initialdir=self.examples_dir,
            title="Save User Database"
        )
        if file_path:
            try:
                full_db = {
                    "materialproperty": database.materialproperty,
                    "alloyproperty": database.alloyproperty,
                    "alloyproperty4": database.alloyproperty4
                }
                with open(file_path, "w") as f:
                    json.dump(full_db, f, indent=4)
                tkinter.messagebox.showinfo("Success", "Database saved successfully.")
            except Exception as e:
                tkinter.messagebox.showerror("Error", str(e))

    def load_database_from_file(self):
        file_path = tkinter.filedialog.askopenfilename(
            filetypes=[("Aestimo Database", "*.json")],
            initialdir=self.examples_dir,
            title="Load User Database"
        )
        if file_path:
            try:
                with open(file_path, "r") as f:
                    full_db = json.load(f)
                
                if "materialproperty" in full_db:
                    database.materialproperty.update(full_db["materialproperty"])
                if "alloyproperty" in full_db:
                    database.alloyproperty.update(full_db["alloyproperty"])
                if "alloyproperty4" in full_db:
                    database.alloyproperty4.update(full_db["alloyproperty4"])
                
                # Refresh current view
                self.update_material_editor(self.db_selector.get())
                tkinter.messagebox.showinfo("Success", "Database loaded and updated in memory.")
            except Exception as e:
                tkinter.messagebox.showerror("Error", str(e))

    def reset_database_defaults(self):
        if not tkinter.messagebox.askyesno("Reset", "Are you sure you want to revert all database changes to default?"):
            return
            
        # Restore from deep copies
        database.materialproperty = copy.deepcopy(self.default_material_property)
        database.alloyproperty = copy.deepcopy(self.default_alloy_property)
        database.alloyproperty4 = copy.deepcopy(self.default_alloy_property4)
        
        # We also need to update the module level variable if possible, 
        # but since we did `import database`, `database.materialproperty = ...` updates the name in that module object.
        # But `from database import materialproperty` creates a local name in this file which won't auto-update if we reassign `database.materialproperty`.
        # However, `aestimo_gui.py` imports them locally too. 
        # WARNING: In `aestimo_gui.py`, we have `from database import materialproperty, alloyproperty`.
        # The `add_layer` method uses `materialproperty` global (local to this file) directly.
        # This means modifying `database.materialproperty` WON'T affect the local `materialproperty` name here unless we update it too.
        # FIX: We should update the local names or use `database.materialproperty` everywhere.
        # Better: Update local names here too.
        
        # Wait, Python modules: `database.materialproperty` is a dict. mutating it works.
        # Re-assigning `database.materialproperty = new_dict` breaks reference for others holding the old dict.
        # Best approach: Clear and update the EXISTING dictionary objects to preserve references.
        
        self._restore_dict(database.materialproperty, self.default_material_property)
        self._restore_dict(database.alloyproperty, self.default_alloy_property)
        self._restore_dict(database.alloyproperty4, self.default_alloy_property4)

        self.update_material_editor(self.db_selector.get())
        tkinter.messagebox.showinfo("Reset", "Database reset to defaults.")

    def _restore_dict(self, target, source):
        target.clear()
        target.update(source)

    def poll_log_file(self):
        """Simple poller to read aestimo.log and update console"""
        try:
            if os.path.exists("aestimo.log"):
                with open("aestimo.log", "r") as f:
                    # Seek to end? No, we want to read new lines.
                    # Simple MVP: read all and keep last N lines, or track position.
                    # We'll just read last 50 lines to keep it simple for now.
                    lines = f.readlines()
                    last_lines = "".join(lines[-50:])
                    
                    # Update text widget if changed
                    # Check if different? roughly
                    current_text = self.console_text.get("0.0", "end")
                    if last_lines.strip() not in current_text:
                        self.console_text.delete("1.0", "end")
                        self.console_text.insert("end", last_lines)
                        self.console_text.see("end")
        except Exception:
            pass
        finally:
            self.after(2000, self.poll_log_file) # Poll every 2s

    def add_new_material_dialog(self):
        """Opens a dialog to add a new material entry"""
        dialog = customtkinter.CTkInputDialog(text="Enter unique name for new material/alloy:", title="New Entry")
        name = dialog.get_input()
        
        if name:
            if name in self.all_materials:
                tkinter.messagebox.showerror("Error", "Material name already exists.")
                return
            
            # Simple Selection for Type
            type_dialog = customtkinter.CTkToplevel(self)
            type_dialog.title("Select Type")
            type_dialog.geometry("300x200")
            type_dialog.after(10, type_dialog.focus_force)
            
            customtkinter.CTkLabel(type_dialog, text=f"Select category for '{name}':").pack(pady=20)
            
            def select_type(t):
                if t == "Material":
                    database.materialproperty[name] = copy.deepcopy(database.materialproperty["GaAs"])
                elif t == "Alloy":
                    database.alloyproperty[name] = copy.deepcopy(database.alloyproperty["AlGaAs"])
                else:
                    database.alloyproperty4[name] = copy.deepcopy(database.alloyproperty4["InGaAsP"])
                
                # Refresh UI
                self.all_materials = sorted(list(database.materialproperty.keys()) + 
                                           list(database.alloyproperty.keys()) + 
                                           list(database.alloyproperty4.keys()))
                self.db_selector.configure(values=self.all_materials)
                self.db_selector.set(name)
                self.update_material_editor(name)
                
                # Update Structure Dropdowns (global ref in this file)
                # Note: add_layer dynamically builds this list, so future adds are okay.
                # However, existing layers won't see the new values in their OptionMenus 
                # unless we refresh them all. For now, new layers will have it.
                
                type_dialog.destroy()
                tkinter.messagebox.showinfo("Success", f"Added '{name}' as {t}. Please edit its properties and Save.")

            customtkinter.CTkButton(type_dialog, text="Material (Base)", command=lambda: select_type("Material")).pack(pady=5)
            customtkinter.CTkButton(type_dialog, text="Alloy (Ternary)", command=lambda: select_type("Alloy")).pack(pady=5)
            customtkinter.CTkButton(type_dialog, text="Quaternary", command=lambda: select_type("Quaternary")).pack(pady=5)

if __name__ == "__main__":
    app = AestimoGUI()
    app.mainloop()
