# -*- coding: utf-8 -*-
"""
Utilize the trained neural network to finish downstream tasks.

Created on Tue 16 Apr 2024 02:22:29 PM CST
@author: Xu Zhu
"""
import sys
import os, re
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cantera as ct
import pandas as pd
import scipy.integrate
import argparse
from tqdm import tqdm
import PyCSP.Functions as csp
import PyCSP.utils as utils
from props_1D import AVAILABLE_CASES

import warnings

_seen_nasa_warning = set()
_original_showwarning = warnings.showwarning  # Save original first!

def nasa_filter(message, category, filename, lineno, file=None, line=None):
    if "NasaPoly2::validate" in str(message):
        key = (str(message), filename, lineno)
        if key in _seen_nasa_warning:
            return  # Suppress repeat
        _seen_nasa_warning.add(key)

    # Call the original function (not the overridden one!)
    _original_showwarning(message, category, filename, lineno, file, line)

warnings.showwarning = nasa_filter



# Parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Run 1D Counterflow Flame simulation")
    parser.add_argument(
        '--case_name',
        type=str,
        choices=AVAILABLE_CASES,
        required=True,
        help=(
            "Specify the case name. "
            f"Valid options: {AVAILABLE_CASES}"
        )
    )
    parser.add_argument(
        '--mech_file',
        type=str,
        required=True,
        help="Full path to the mechanism file (cti or yaml)"
    )
    parser.add_argument(
        '--element_num',
        type=str,
        required=True,
        help="Number of elements in the mechanism (e.g., '5' for 5 elements)"
    )
    
    return parser.parse_args()

args = parse_args()
case_name = args.case_name
mech_file = args.mech_file
element_num = args.element_num
mech_folder = os.path.dirname(mech_file)          # "mechanism/NH3_otomo"
mech_filename = os.path.basename(mech_file)        # "NH3_otomo.cti"
mech_name, _ = os.path.splitext(mech_filename)     # ("NH3_otomo", ".cti")

from props_1D import generate_props  # Import the function to generate props
props_list = generate_props(case_name)
sys.path.append(os.path.abspath(mech_folder))
# from RCCE import RCCE
from subprocess import run
from pathlib import Path
import graphviz
import random
import glob
import matplotlib as mpl
import matplotlib.font_manager as fm
import subprocess
import matplotlib.ticker as mticker

# Path to the Times New Roman font file
font_path = '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf'

# Create a font properties object using the font file
font_prop = fm.FontProperties(fname=font_path)
inch = 2.54

class ReactorOde:
    def __init__(self, gas, main_species):
        """Initialize the ReactorOde object with the gas object."""
        self.gas = gas
        self.P = gas.P
        self.main_species = main_species
        species_names = gas.species_names
        self.mask = np.array([species in self.main_species for species in species_names])

    def __call__(self, t, y):
        """The ODE function, y' = f(t,y)"""
        
        self.gas.set_unnormalized_mass_fractions(y[1:])
        self.gas.TP = y[0], self.P  
        rho = self.gas.density
        HRR = self.gas.net_production_rates
        # HRR[self.mask] = 0.0
        dTdt = np.zeros_like(y[0]) 
        dYdt = HRR * self.gas.molecular_weights / rho 
        dYdt[self.mask] = 0.0
        return np.hstack((dTdt, dYdt)) 

def run_CEQ(gas, fig_dir, main_species_names, initial_fractions, fuel, oxyd, T, P, dt, time_end, full_data, element_num, mech_file, moleFraction):
    """
    Run the combustion simulation, extract data, and save results to a file.

    Parameters:
    - gas: Cantera gas object
    - fig_dir: Directory for saving results
    - main_species_names: List of species to control (main species)
    - initial_mole_fractions: Dictionary of initial mole fractions for controlled species
    - fuel: Fuel species name
    - oxyd: Oxidizer species name
    - T: Temperature
    - P: Pressure
    - dt: Time step for the simulation
    - time_end: End time for the simulation
    - full_data: Additional data for processing
    """
    # Call the core simulation function
    gas = run_CEQ_core(gas, main_species_names, initial_fractions, T, P, dt, time_end, moleFraction)

    # Extract mole fractions and species names
    species_names = gas.species_names
    temp, pressure = gas.T, gas.P
    mole_fractions = gas.X
    
    # Create a DataFrame for the results with modified species names
    data = {}

    # Add "X_" or "Y_" to the species names based on moleFraction value
    for species, y in zip(species_names, mole_fractions):
        if moleFraction:
            new_species_name = f"X_{species}"
        else:
            new_species_name = f"Y_{species}"
        data[new_species_name] = y
    
    updated_data = {}
    for species, value in data.items():
        # Mole‑fraction keys start with the prefix "X_"
        if species.startswith("X_"):
            # Strip the prefix to find the base species name
            base_species = species[2:]
            # Use the initial mole‑fraction (if provided); otherwise keep the original value
            updated_data[species] = initial_fractions.get(base_species, value)
        elif species.startswith("Y_"):
            # Strip the prefix to find the base species name
            base_species = species[2:]
            # Use the initial mole‑fraction (if provided); otherwise keep the original value
            updated_data[species] = initial_fractions.get(base_species, value)
        else:
            # Keys without the prefix are copied unchanged
            updated_data[species] = value
    df = pd.DataFrame(updated_data, index=[0])

    # Perform additional calculations
    # 1. Make sure the Cantera object you pass in **is** a PyCSP wrapper
    #    (run_CEQ_core should create it via csp.CanteraCSP).
    #    gas.jacobiantype should already be 'full'.

    # 2. Solve eigen-system of the chemical-source Jacobian
    #    D is a (ns, ns) diagonal matrix; L, R are left/right eigenvectors.
        #push pressure
    gas_csp = csp.CanteraCSP(mech_file)
    gas_csp.TPY = temp, pressure, mole_fractions
    gas_csp.constP = pressure
    gas_csp.jacobiantype = 'full'
    lam,R,L,f = gas_csp.get_kernel()
    real_lam = lam.real
    sorted_idx = np.argsort(np.abs(real_lam))
    cem_idx = sorted_idx[-1 - element_num]
    cem_real = real_lam[cem_idx]
    TLOG = np.sign(cem_real) * np.log10(1.0 + abs(cem_real))
    # Calculate mixture fraction and heat release rate
    MF = gas.mixture_fraction(fuel, oxyd)
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density

    # Add calculated values to the DataFrame
    df['CEM'] = TLOG
    df['MF'] = MF
    df['HRR'] = HRR
    if 't_res' in full_data:
        df['t_res'] = full_data['t_res']
    if 'P' in full_data:
        df['P'] = full_data['P']
    if 'x' in full_data:
        df['x'] = full_data['x']
    if 'grid' in full_data:
        df['grid'] = full_data['grid']
    if 'T' in full_data:
        df['T'] = full_data['T']
    if 'Normalized X' in full_data:
        df['Normalized X'] = full_data['Normalized X']
    if 'Normalized Y' in full_data:
        df['Normalized Y'] = full_data['Normalized Y']
        
# Save the results to a CSV file
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close()

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header=not file_exists, index=False)
        

def run_CEQ_core(
    gas,
    main_species_names,
    initial_fractions,
    T,
    P,
    dt,
    time_end,
    moleFraction,
    steady_time=0.1,
    steady_tol=1e-12
):
    """
    Core function to run the combustion simulation and return the updated gas object.
    
    Stops early if the solution reaches steady state (no change > steady_tol after steady_time).
    """
    # 1) Identify inert and build composition vector
    inert_specie    = gas.species_names[-1]
    total_controlled = sum(initial_fractions.values())
    inert_value      = max(0.0, 1.0 - total_controlled)

    X = np.zeros(gas.n_species)
    for sp, mf in initial_fractions.items():
        X[gas.species_index(sp)] = mf
    if inert_value > 0.0:
        X[gas.species_index(inert_specie)] = inert_value

    # 2) Apply to gas
    gas.TP = T, P
    if moleFraction:
        gas.set_unnormalized_mole_fractions(X)
    else:
        mw = gas.molecular_weights
        # Y = X * mw / np.sum(X * mw)
        gas.set_unnormalized_mass_fractions(X)

    # 3) Pack initial state for the ODE solver
    y0 = np.hstack((gas.T, gas.Y))
    ode = ReactorOde(gas, main_species_names)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator(
        'vode',
        method='bdf',
        with_jacobian=True,
        # rtol=1e-10,
        # atol=1e-12,
        # max_step=dt
    )
    solver.set_initial_value(y0, 0.0)

    # 4) Integrate, with early‐exit on steady state
    prev_y = y0.copy()
    while solver.successful() and solver.t < time_end:
        solver.integrate(solver.t + dt)
        y_new = solver.y
        gas.TPY = y_new[0], P, y_new[1:]

        # once we're past a small burn-in, test for steady state
        if solver.t >= steady_time:
            if np.allclose(y_new, prev_y, atol=steady_tol, rtol=0.0):
                # nothing changed -> assume steady and break
                break
        prev_y[:] = y_new

    # 5) Zero out inert and return
    Y_final = gas.Y.copy()
    Y_final[gas.species_index(inert_specie)] = 0.
    gas.TPY = gas.T, P, Y_final
    return gas


def run_RCCE_fortran(gas, fig_dir, T, P, modeci,element_num, full_data, props):
    """
    Run the simulation to monitor the concentrations of controlled and target species.

    Parameters:
    - gas: Cantera gas object
    - main_species: List of species to control (main species)
    - initial_mole_fractions: Dictionary of initial mole fractions for controlled species
    - T: Temperature
    - P: Pressure
    - time_step: Time step for the simulation
    - time_end: End time for the simulation
    """
    ##### Define the path to the input and output files
    input_file_path = "src/recon_fortran/fc_xe_fe.op"
    if modeci == "RCCE_fortran" : 
        output_file_path = "src/recon_fortran/rcce_result.op"
    if modeci == "ICE_PIC":
        output_file_path = "src/recon_fortran/icepic_result.op"
        
    fortran_executable = os.path.abspath("src/recon_fortran/recon_species")

    def modify_input_file(fortran_input):
        # Read the input file content
        with open(input_file_path, "r") as file:
            lines = file.readlines()
        
        # Substitute the second line with `fortran_input`
        lines[1] = ' '.join(map(str, fortran_input)) + '\n'
        
        # Write the modified content back to the file
        with open(input_file_path, "w") as file:
            file.writelines(lines)
            
    all_species_names = gas.species_names
    T = gas.T
    mass_fractions = gas.Y  # Mass fractions
    enthalpies_SI = gas.h   # Enthalpies in J/kmol
    enthalpies_CGS = enthalpies_SI / 1e3 * 1e7
    molecular_weights = gas.molecular_weights  # in g/mol
    specific_mole_numbers = mass_fractions / molecular_weights
    enthalpies_CGS_array = np.array([enthalpies_CGS])
    fortran_input = np.concatenate((specific_mole_numbers, enthalpies_CGS_array))

    modify_input_file(fortran_input)
    
    def run_fortran_code():
        subprocess.run([fortran_executable], check=True, cwd="src/recon_fortran", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    run_fortran_code()

    def read_output_file(output_file_path, molecular_weights):
        # Read the output data from the file using pandas
        # df = pd.read_csv(output_file_path, delim_whitespace=True, header=None)
        df = pd.read_csv(output_file_path, sep=r"\s+", header=None)

        # Convert the DataFrame to a NumPy array of floats
        W_recon_RCCE = np.array(df.values.flatten(), dtype=float)

        # Ensure molecular_weights is a NumPy array of floats
        molecular_weights = np.array(molecular_weights, dtype=float)

        # Calculate Y_recon_RCCE
        Y_recon_RCCE = W_recon_RCCE[:-1] * molecular_weights
        T_recon_RCCE = W_recon_RCCE[-1]
        return Y_recon_RCCE, T_recon_RCCE
    
    Y_recon_RCCE, T_recon_RCCE = read_output_file(output_file_path, molecular_weights)
    
    Y_recon_mass_fractions = dict(zip(all_species_names, Y_recon_RCCE))
    
    gas.TPY = T_recon_RCCE, P, Y_recon_mass_fractions
    
    # here we deal with the new gas 
    mole_fractions = gas.X
    data = {species: y for species, y in zip(all_species_names, mole_fractions)}
    df = pd.DataFrame(data, index=[0])
    # mass_fractions, MF = get_mass_fractions_from_1D_flame(f, gas, props)
    
    gas_csp = csp.CanteraCSP(mech_file)
    gas_csp.TPY = T_recon_RCCE, P, mole_fractions
    gas_csp.constP = P
    gas_csp.jacobiantype = 'full'
    lam,R,L,f = gas_csp.get_kernel()
    real_lam = lam.real
    sorted_idx = np.argsort(np.abs(real_lam))
    cem_idx = sorted_idx[-1 - element_num]
    cem_real = real_lam[cem_idx]
    TLOG = np.sign(cem_real) * np.log10(1.0 + abs(cem_real))
    # 7. Append to DataFrame
    df['CEM'] = TLOG
    # Calculate mixture fraction and heat release rate
    MF = gas.mixture_fraction(props['fuel'], props['oxidizer'])
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density

    # Add calculated values to the DataFrame
    df['CEM'] = TLOG
    df['MF'] = MF
    df['HRR'] = HRR
    if 't_res' in full_data:
        df['t_res'] = full_data['t_res']
    if 'x' in full_data:
        df['x'] = full_data['x']
    if 'grid' in full_data:
        df['grid'] = full_data['grid']
    if 'T' in full_data:
        df['T'] = full_data['T']
    if 'Normalized X' in full_data:
        df['Normalized X'] = full_data['Normalized X']
    if 'Normalized Y' in full_data:
        df['Normalized Y'] = full_data['Normalized Y']

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close() 

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header = not file_exists, index=False)
        

def create_reaction_path_diagram(gas, fig_dir, element, cantera_row, pred_cantera_row, mole_fraction=True):
    # Prepare species compositions from cantera row
    all_species_names = gas.species_names
    # if mole fractions 
    if mole_fraction:
        all_species_values = [cantera_row[f'X_{name}'] for name in all_species_names]
    else:
        all_species_values = [cantera_row[f'Y_{name}'] for name in all_species_names]
    all_fractions = dict(zip(all_species_names, all_species_values))
    species_composition = ", ".join([f"{sp}:{mf}" for sp, mf in all_fractions.items()])
    gas.TPX = cantera_row['T'], cantera_row['P'], species_composition
    diagram = ct.ReactionPathDiagram(gas, element)
    diagram.show_details = True
    diagram.title = 'Reaction path diagram following {0}'.format(element)
    diagram.label_threshold = 1e-20

    # Define the output directory and files
    output_dir = fig_dir
    dot_file = f'{element}_rxnpath.dot'
    modified_dot_file = f'{element}_rxnpath_modified.dot'
    img_file = f'{element}_rxnpath.pdf'  # Change to PDF

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    dot_path = os.path.join(output_dir, dot_file)
    img_path = os.path.join(output_dir, img_file)
    diagram.write_dot(dot_path)

    ####IF YOU DO DO NOT NEED TO MODIFY THE DOT FILE, YOU CAN SKIP THIS PART####
    print(f"Wrote graphviz input file to '{dot_path}'.")
    modify_dot_file(dot_path, cantera_row, pred_cantera_row, all_species_names)
    dot_path = os.path.join(output_dir, modified_dot_file)
    ############################################################################
    
    # Change the Graphviz command to produce a PDF instead of PNG
    run(['dot', '-Tpdf', dot_path, '-o', img_path, '-Gdpi=200'], check=True)
    print(f"Wrote graphviz output file to '{img_path}'.")
    
    # Display the reaction path diagram
    graphviz.Source(diagram.get_dot())

def modify_dot_file(dot_file_path, cantera_row, pred_cantera_row, all_species_names):
    # Read .dot
    with open(dot_file_path, 'r') as f:
        dot_content = f.read()

    species_pattern = re.compile(r'label="([A-Za-z0-9()]+)"')

    # Normalize for matching
    all_species_lower = [sp.lower() for sp in all_species_names]

    # Revised compute_error
    def compute_error(species_name):
        key = f"X_{species_name}"
        cantera_value = cantera_row.get(key, 0.0)
        pred_value    = pred_cantera_row.get(key, 0.0)
        if cantera_value != 0.0:
            error = (pred_value - cantera_value) / cantera_value
        else:
            error = 0.0
        print(f"Species: {species_name}, key: {key}, "
              f"Cantera: {cantera_value}, Pred: {pred_value}, Error: {error}")
        return error

    def error_to_color(err):
        pct = min(abs(err) * 100, 100) / 100
        r = int(255 * pct) + 100
        g = int(255 * (1 - pct)) + 100
        return f"#{min(r,255):02x}{min(g,255):02x}00"

    def modify_label_and_color(m):
        sp = m.group(1).strip()
        if sp.lower() in all_species_lower:
            err = compute_error(sp)
            pct = "set" if abs(err) < 1e-3 else f"{err*100:.2f}%"
            col = error_to_color(err)
            new_label = f'{sp} {pct}'
            return f'label="{new_label}", style=filled, fillcolor="{col}"'
        else:
            return m.group(0)

    modified = species_pattern.sub(modify_label_and_color, dot_content)

    out_path = dot_file_path.replace('.dot', '_modified.dot')
    with open(out_path, 'w') as f:
        f.write(modified)
    print(f"Modified dot file saved as {out_path}")

import os, random
import pandas as pd
from typing import Sequence, Optional

# ──────────────────────────────────────────────────────────────────────
#  Main driver
# ──────────────────────────────────────────────────────────────────────
def process_experiments_RCCE(
    gas,
    props: dict,
    mech,
    file_name: str,
    fig_dir: str,
    main_species_names: Sequence[str],
    fuel: str,
    oxidizer: str,
    modeci: str,
    dt: float,
    time_end: float,
    rand: float = 0.0,
    mole_fraction: bool = True,
    assist_species: Optional[Sequence[str]] = ('H2O', 'CO2'),
) -> None:
    """
    Read experimental (and optional *assist*) CSV rows and run CEQ / RCCE / ILDM.

    assist_species  : list of species that, if present, are taken from the
                      *assist* CSV instead of the main CSV.
    """

    # -------------------------------------------------- data & utilities
    all_species_names = gas.species_names
    df_main           = pd.read_csv(file_name)

    # optional assist file
    df_assist = None
    if props.get("assist_path") and assist_species:
        print("assist path is read here!")
        assist_file = os.path.join(props["assist_path"], props["file_name"])
        df_assist   = pd.read_csv(assist_file)
    # -------------------------------------------------- mechanism objects
    if modeci == "RCCE":
        rcce = RCCE(mech, main_species_names, use_element_constraints=True)
  
    # -------------------------------------------------- iterate over rows
    
    # I found that we need to interpolate first 
    for idx, row in df_main.iterrows():
        input_T = row["T"]
        input_P = props["P"]
        # ------------ build input_species_values with assist override ----
        input_species_values = []
        for sp in main_species_names:
            col = f"{'X' if mole_fraction else 'Y'}_{sp}"
            # default from main CSV
            value = row[col]

            # override from assist CSV if requested & available
            if (
                df_assist is not None
                and sp in assist_species
                and col in df_assist.columns
            ):
                value = df_assist.loc[idx, col]
            # ppm → mole fraction for NO
            if sp == "NO":
                value = value / 1e6

            input_species_values.append(value)
            

        # ------------ random perturbation (optional) ---------------------
        perturb = [
            val * rand * random.uniform(-1, 1) for val in input_species_values
        ]
        perturbed = [max(v + dv, 0.0) for v, dv in zip(input_species_values, perturb)]
        initial_mf = dict(zip(main_species_names, perturbed))

        # ------------ complete mixture (for RCCE / ILDM) -----------------
        if modeci != "QSSA":
            all_vals = [
                initial_mf.get(sp, row[sp]) for sp in all_species_names
            ]
            mole_fractions = dict(zip(all_species_names, all_vals))

        # -------------------------------------------------- run solvers
        if modeci == "QSSA":
            run_CEQ(gas, fig_dir, main_species_names, initial_mf,
                fuel, oxidizer, input_T, input_P, mech_file=props['mech_file'],
                full_data=row, dt=dt, time_end=time_end,
                element_num=props["element_num"], moleFraction=mole_fraction)
        elif modeci == "RCCE":
            gas.TPX = input_T, input_P, mole_fractions
            rcce.set_conditions(*gas.TPX)
            run_RCCE(gas, rcce, fig_dir, props)
        else:
            raise ValueError(f"Unknown modeci: {modeci}")


def process_simulation_RCCE(gas,props, mech, file_name, fig_dir, main_species_names, fuel, oxidizer, modeci, dt, time_end,denoted_value_name, rand=0, mole_fraction=True):
    # -------------------------------------------------- data & utilities
    all_species_names = gas.species_names
    df_main           = pd.read_csv(file_name, comment="#")
    # optional assist file
    df_assist = None
    if props.get("assist_path") and assist_species:
        print("assist path is read here!")
        assist_species = props['species_from_exp']
        assist_file = os.path.join(props["assist_path"], props["file_name"])
        df_assist   = pd.read_csv(assist_file, comment="#")
    # -------------------------------------------------- mechanism objects

  # -------------------------------------------------- iterate over rows
    # max_hrr_value,max_hrr_index = 0,0
    max_T_value, max_T_index = 0,0 
    # for idx, row in df_main.iterrows():
    for idx, row in tqdm(df_main.iterrows(),
                     total=len(df_main),
                     desc=f"Running {modeci}",
                     unit="case"):
        
        input_T = row["T"]
        input_P = props["P"]
        denoted_value = row[denoted_value_name]
        
        if denoted_value > max_T_value:
            max_T_value = denoted_value
            max_T_index = idx
        # ------------ build input_species_values with assist override ----
        input_species_values = []
        for sp in main_species_names:
            col = f"{'X' if mole_fraction else 'Y'}_{sp}"
            # default from main CSV
            value = row[col]
            # override from assist CSV if requested & available
            if (
                df_assist is not None
                and sp in assist_species
                and col in df_assist.columns
            ):
                value = df_assist.loc[idx, col]
            input_species_values.append(value)

        perturbation = [value * rand * random.uniform(-1, 1) for value in input_species_values]
        perturbed_species_values = [value + perturb for value, perturb in zip(input_species_values, perturbation)]
        
        initial_model_fractions = dict(zip(main_species_names, perturbed_species_values))
        if modeci != "QSSA":
            if mole_fraction == True:
                all_species_values = [row[f'X_{name}'] if name not in main_species_names else 
                                      initial_model_fractions.get(name, row[f'X_{name}']) 
                                      for name in all_species_names]
            else:
                all_species_values = [row[f'Y_{name}'] if name not in main_species_names else 
                                      initial_model_fractions.get(name, row[f'Y_{name}']) 
                                      for name in all_species_names]
            
            mole_fractions = dict(zip(all_species_names, all_species_values))
            
        if modeci == "RCCE_fortran":
            gas.TPX = input_T, input_P, mole_fractions
            mole_fractions['NNH'] = 0
            run_RCCE_fortran(gas, fig_dir, input_T, input_P, modeci,props['element_num'],row, props)

        if modeci == "QSSA":
            run_CEQ(gas, fig_dir, main_species_names, initial_model_fractions, fuel, oxidizer, input_T, input_P, full_data=row, dt=dt, time_end=time_end, element_num=props['element_num'], mech_file=mech,moleFraction=props["ifMoleFraction"])
            


    print("Index for maximum HRR:", max_T_index)
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    df_pred = pd.read_csv(predicted_file_dir, comment='#')
        
    # Visualize reaction path for largest HRR
    max_T_row = df_main.loc[max_T_index]
    pred_T_row = df_pred.loc[max_T_index]
    print("max_hrr_index=",max_T_index)
    print("file_name=",file_name)
    print("predicted_file_dir=", predicted_file_dir)
    print("max_T_row=",max_T_row)
    print("pred_T_row=",pred_T_row)
    elements = ['O','H','N']
    for element in elements:
        create_reaction_path_diagram(gas,fig_dir, element, max_T_row, pred_T_row,mole_fraction)

        
def create_fig_dir(file_name,mech_name,prefix, extra_name=None):
    base_name, ext = os.path.splitext(file_name)
    # If extra_name is provided, append it to the base name
    if extra_name:
        base_name = f"{base_name}_{extra_name}"
    fig_dir = os.path.join("figs", prefix, base_name,mech_name)
    os.makedirs(fig_dir, exist_ok=True)
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    if os.path.isfile(predicted_file_dir):
        os.remove(predicted_file_dir)
    return fig_dir

def create_noisy_fig_dir(file_name, prefix, noise_level):
    base_name, ext = os.path.splitext(file_name)
    # Round the noise level to 3 decimal places
    noise_level_rounded = round(noise_level, 3)
    # Format the directory path with the rounded noise level
    fig_dir = os.path.join("figs", prefix, f'{base_name}_{noise_level_rounded}')
    os.makedirs(fig_dir, exist_ok=True)
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    if os.path.isfile(predicted_file_dir):
        os.remove(predicted_file_dir)
    return fig_dir



def modify_streams_in_file(main_species_names, modeci):
    input_file = 'src/recon_fortran/streams_ori.in'
    output_file = 'src/recon_fortran/streams.in'

    with open(input_file, 'r') as file:
        lines = file.readlines()

    start_index = None
    end_index = None

    for i, line in enumerate(lines):
        if "REPRESENTED_SPECIES BEGIN" in line:
            start_index = i
        if "REPRESENTED_SPECIES END" in line:
            end_index = i

    if start_index is not None and end_index is not None:
        new_species_line = " ".join(main_species_names) + "\n"
        lines = lines[:start_index + 1] + [new_species_line] + lines[end_index:]

    # Modify the MODECI line
    for i, line in enumerate(lines):
        if line.startswith("MODECI"):
            lines[i] = f"MODECI         {modeci}\n"
            break

    with open(output_file, 'w') as file:
        file.writelines(lines)

    print(f"Modified file saved as '{output_file}'")

    
def main():
    gas = ct.Solution(mech_file, thermo='ideal-gas', kinetics='GasKinetics')
    # gas = ct.Solution(thermo='ideal-gas', kinetics='GasKinetics',
    #                 species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
    #                 reactions=gas.reactions())
    print(f"Number of species: {gas.n_species}")
    print(f"Number of reactions: {gas.n_reactions}")
    print("species_names:", gas.species_names)
    
    for props in props_list: 
        props['mech_file'] = mech_file
        props['mech_path'] = mech_folder
        props['mech_name'] = mech_name
        props['element_num'] = int(element_num)
        props['perturb'] = 0
        
        file_name = props['file_name']
        data_dir = props['path']
        
        file_name = os.path.join(data_dir, file_name)
        print("file_name=", file_name)
        # valid_modes = {"RCCE_fortran", "QSSA"} # RCCE_fortran
        valid_modes = { "QSSA"} # RCCE_fortran
        
        for modeci in valid_modes: 

            # main_species_names = ["NH3","H2", "O2", "H2O", "N2","OH", "NH2", "NNH"]

            main_species_names = ["NH3", "H2", "O2", "H2O", "N2"]
            extra_name = None
            # main_species_names = ["NH3", "H2", "O2", "H2O", "N2"]
            # extra_name = None

            denoted_value_name = "HRR"
            # main_species_names = [ "H2", "O2", "H2O", "N2","NH3", "O"]

            # main_species_names = ['NH3', 'NO2','H2', 'O2', 'H', 'O', 'OH', 'HO2', 'H2O', 'H2O2', 'NH2', 'NH', 'N', 'NNH', 'NH2OH', 'H2NO', 'HNOH', 'HNO', 'HON','HONO', 'HNO2', 'NO3', 'HONO2', 'N2O', 'N2H4', 'N2H3', 'N2H2', 'H2NN', 'N2']

            print("Reconstruction starts with method:", modeci)
            if modeci not in valid_modes:
                print(f"Warning: Invalid modeci '{modeci}'. Using default 'QSSA'.")
            if modeci == "RCCE_fortran":
                modify_streams_in_file(main_species_names, 'RCCE')
            fig_dir = create_fig_dir(file_name, mech_name,modeci,extra_name)
            print("fig_dir=",fig_dir)
            if props['Exp_data'] == True:
                speciesFromSim = None
                if props['if_assist'] == True:
                    speciesFromSim = props['species_from_sim']
                process_experiments_RCCE(gas, props, props['mech_file'], file_name, fig_dir, main_species_names,props['fuel'], props['oxidizer'], modeci, props['dt'], props['time_end'], rand=props['perturb'], mole_fraction=props['ifMoleFraction'], assist_species=speciesFromSim)
            else:
                file_name = os.path.join(props['path'], mech_name, props['file_name'])
                print("file_name=", file_name)  
                process_simulation_RCCE(gas, props, props['mech_file'], file_name, fig_dir, main_species_names,props['fuel'], props['oxidizer'], modeci, props['dt'], props['time_end'], denoted_value_name = denoted_value_name,rand=props['perturb'], mole_fraction=props['ifMoleFraction'])
        
        
if __name__ == "__main__":
    main()