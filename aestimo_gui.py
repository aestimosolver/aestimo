
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
from aestimo import run_aestimo
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

        # Action Buttons
        self.run_button = customtkinter.CTkButton(self.sidebar_frame, text="RUN SIMULATION", 
                                                 command=self.start_simulation_thread, 
                                                 height=50, 
                                                 font=customtkinter.CTkFont(size=14, weight="bold"),
                                                 fg_color="#1F6AA5", hover_color="#144870")
        self.run_button.grid(row=2, column=0, padx=20, pady=10)
        
        self.save_btn = customtkinter.CTkButton(self.sidebar_frame, text="Save Project", command=self.save_project, fg_color="green", hover_color="darkgreen")
        self.save_btn.grid(row=3, column=0, padx=20, pady=5)
        
        self.load_btn = customtkinter.CTkButton(self.sidebar_frame, text="Load Project", command=self.load_project, fg_color="#D35400", hover_color="#A04000")
        self.load_btn.grid(row=4, column=0, padx=20, pady=5)
        
        # Progress Bar
        self.progress_bar = customtkinter.CTkProgressBar(self.sidebar_frame, orientation="horizontal")
        self.progress_bar.grid(row=5, column=0, padx=20, pady=20)
        self.progress_bar.set(0)

        # --- Main Content Area (Tabs) ---
        self.tabview = customtkinter.CTkTabview(self, width=800)
        self.tabview.grid(row=0, column=1, padx=20, pady=10, sticky="nsew")
        
        self.tab_structure = self.tabview.add("Structure")
        self.tab_physics = self.tabview.add("Physics & Environment")
        self.tab_solver = self.tabview.add("Solver & Grid")
        self.tab_results = self.tabview.add("Results")
        self.tab_console = self.tabview.add("Console")

        self.setup_structure_tab()
        self.setup_physics_tab()
        self.setup_solver_tab()
        self.setup_results_tab()
        self.setup_console_tab()

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

        self.clear_btn = customtkinter.CTkButton(tools_frame, text="Clear All", command=self.clear_layers, fg_color="#D94848", hover_color="#A02020", width=80)
        self.clear_btn.pack(side="right", padx=10, pady=5)

        # Scrollable Area for Layers
        self.layers_frame = customtkinter.CTkScrollableFrame(self.tab_structure, label_text="Heterostructure Definition (Growth Direction ↓)")
        self.layers_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.layers_frame.grid_columnconfigure(0, weight=1)
        
        self.layer_widgets = []
        
        # Initial Demo Structure (Based on sample_1qw_barrierdope_ingaas)
        # Layer 1: p-doped AlGaAs contact/barrier
        self.add_layer(material="AlGaAs", thickness=250.0, type="barrier", mole=0.3, doping=0.8e17, doping_type="p")
        # Layer 2: Intrinsic AlGaAs spacer
        self.add_layer(material="AlGaAs", thickness=50.0, type="barrier", mole=0.3, doping=0.0, doping_type="n")
        # Layer 3: Intrinsic GaAs Quantum Well
        self.add_layer(material="GaAs", thickness=15.0, type="well", mole=0.0, doping=0.0, doping_type="n")
        # Layer 4: Intrinsic AlGaAs spacer
        self.add_layer(material="AlGaAs", thickness=5.0, type="barrier", mole=0.3, doping=0.0, doping_type="n")
        # Layer 5: n-doped AlGaAs barrier
        self.add_layer(material="AlGaAs", thickness=20.0, type="barrier", mole=0.3, doping=0.8e18, doping_type="n")
        # Layer 6: n-doped GaAs contact
        self.add_layer(material="GaAs", thickness=15.0, type="barrier", mole=0.0, doping=0.8e18, doping_type="n")
        
    def setup_physics_tab(self):
        """Environment and Physics Models"""
        frame = self.tab_physics
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
        self.create_input_row(sweep_frame, "V Max (V):", "vmax_entry", "1.0")
        self.create_input_row(sweep_frame, "Step (V):", "vstep_entry", "0.05")

    def setup_solver_tab(self):
        """Computational Settings"""
        frame = self.tab_solver
        
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
        self.solver_combo.set("2: Schrodinger-Poisson")

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
        
        # Subbands
        state_frame = customtkinter.CTkFrame(frame)
        state_frame.pack(fill="x", padx=20, pady=10)
        customtkinter.CTkLabel(state_frame, text="Quantum States", font=customtkinter.CTkFont(weight="bold")).pack(anchor="w", padx=10, pady=5)
        self.create_input_row(state_frame, "Electron Subbands:", "sub_e_entry", "5")
        self.create_input_row(state_frame, "Hole Subbands:", "sub_h_entry", "5")

    def setup_results_tab(self):
        """Results Display"""
        self.tab_results.grid_columnconfigure(0, weight=1)
        self.tab_results.grid_rowconfigure(0, weight=1)
        
        # A container for the matplotlib canvas
        self.results_container = customtkinter.CTkFrame(self.tab_results)
        self.results_container.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        
        # Placeholder
        customtkinter.CTkLabel(self.results_container, text="Run a simulation to see results here.").pack(expand=True)

    def setup_console_tab(self):
        self.tab_console.grid_columnconfigure(0, weight=1)
        self.tab_console.grid_rowconfigure(0, weight=1)
        
        self.console_text = customtkinter.CTkTextbox(self.tab_console, font=("Consolas", 12))
        self.console_text.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        self.console_text.insert("0.0", "--- Aestimo Console ---\n")

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

    def load_project(self):
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
            input_data = self.get_current_configuration()
            # Build InputObject here or inside thread? 
            # Safer to build InputObject structure here to catch validation errors immediately
            # But the object class needs to be passed.
            # Let's pass the config dict to thread and build object there.
            
            thread = threading.Thread(target=self.run_simulation_worker, args=(input_data,))
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
                
                vmax = val_vmax
                vmin = val_vmin
                Each_Step = val_vstep
                
                surface = np.array([bc_left, bc_right])
                
                Quantum_Regions = False
                Quantum_Regions_boundary = np.zeros((2,2)) # Default
                dop_profile = None 

            # Fix annoying class attribute name mismatch if any (T vs T_val)
            InputObject.T = InputObject.T_val # just in case

            # Doping Profile
            tot_thick = sum(row[0] for row in material_list) * 1e-9
            dx_m = grid_step * 1e-9
            n_max = int(tot_thick / dx_m)
            
            dop_arr = np.zeros(n_max)
            curr = 0
            for row in material_list:
                th_m = row[0] * 1e-9
                val = row[4]
                dtype = row[5]
                if dtype == 'p': val = -val
                
                steps = int(th_m / dx_m)
                end = min(curr + steps, n_max)
                dop_arr[curr:end] = val
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

    def finish_simulation(self, success, error_msg=None, figures=None):
        self.is_simulating = False
        self.progress_bar.stop()
        self.run_button.configure(state="normal")
        
        if success:
            self.status_label.configure(text="Status: Simulation Complete", text_color="green")
            self.display_figures(figures)
            tkinter.messagebox.showinfo("Success", "Simulation completed successfully.")
        else:
            self.status_label.configure(text="Status: Error", text_color="red")
            tkinter.messagebox.showerror("Error", f"Simulation Failed:\n{error_msg}")

    def display_figures(self, figures):
        # Clear previous
        for widget in self.results_container.winfo_children():
            widget.destroy()

        if not figures:
            customtkinter.CTkLabel(self.results_container, text="No figures generated.").pack()
            return
            
        fig_tabview = customtkinter.CTkTabview(self.results_container)
        fig_tabview.pack(fill="both", expand=True)
        
        titles = ["Band Diagram", "Charge Density", "Field", "I-V / Sweep", "Other"]
        
        for i, fig in enumerate(figures):
            name = titles[i] if i < len(titles) else f"Figure {i+1}"
            fig_tabview.add(name)
            
            canvas = FigureCanvasTkAgg(fig, master=fig_tabview.tab(name))
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)
            
            toolbar = NavigationToolbar2Tk(canvas, fig_tabview.tab(name))
            toolbar.update()
            canvas.get_tk_widget().pack(fill="both", expand=True)


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
