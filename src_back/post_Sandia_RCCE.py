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
from props_Sandia import * 
sys.path.append(os.path.abspath(props_list[0]['mech_path']))
from subprocess import run
from pathlib import Path
import graphviz
from lib_cema import *
import random
import glob
import matplotlib as mpl
import matplotlib.font_manager as fm
import subprocess
from run_1D_CEQ import run_CEQ,run_CEQ_mass
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
        Qdot = self.gas.net_production_rates
        # Qdot[self.mask] = 0.0
        dTdt = np.zeros_like(y[0]) 
        dYdt = Qdot * self.gas.molecular_weights / rho 
        dYdt[self.mask] = 0.0
        return np.hstack((dTdt, dYdt)) 


def run_CEQ_stop(gas, fig_dir, initial_species_names, initial_mole_fractions, T, P, dt, time_end, full_data):
    """
    Run the simulation to monitor the concentrations of controlled and target species.

    Parameters:
    - gas: Cantera gas object
    - main_species_names: List of species to control (main species)
    - initial_mole_fractions: Dictionary of initial mole fractions for controlled species
    - T: Temperature
    - P: Pressure
    - dt: Time step for the simulation
    - time_end: End time for the simulation
    """
    time_points = []  # Store time
    species_names = gas.species_names
    all_species_concentration = {sp: [] for sp in species_names}
    species_composition = ", ".join([f"{sp}:{mf}" for sp, mf in initial_mole_fractions.items()])
    # gas.TPY = T, P, species_composition
    gas.TPX = T, P, species_composition
    
    y0 = np.hstack((gas.T, gas.Y))
    # Initialize the ODE system
    ode = ReactorOde(gas, initial_species_names)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True)
    solver.set_initial_value(y0, 0.0)

    oh_checked = False
    while solver.successful() and solver.t < time_end:
        solver.integrate(solver.t + dt)
        time_points.append(solver.t)
        gas.TPY = solver.y[0], P, solver.y[1:]
        for sp in species_names:
            all_species_concentration[sp].append(gas[sp].Y)
        if not oh_checked and gas['OH'].X > full_data['OH']:
            # print("Condition met: 'OH' concentration exceeded the threshold in full_data. Reinitializing...")
            oh_checked = True
            if 'OH' not in initial_species_names:
                initial_species_names.append('OH')
            ode = ReactorOde(gas, initial_species_names)
            y0 = np.hstack((gas.T, gas.Y))
            solver = scipy.integrate.ode(ode)
            solver.set_integrator('vode', method='bdf', with_jacobian=True)
            solver.set_initial_value(y0, time_points[-1])
            continue  # Continue from the updated state with the modified system
    mole_fractions = gas.X
    data = {species: y for species, y in zip(species_names, mole_fractions)}
    df = pd.DataFrame(data, index=[0])
    
    # mass_fractions, MF = get_mass_fractions_from_1D_flame(f, gas, props)
    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num= props['element_num'])
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(props['fuel'], props['oxidizer'])
    Qdot = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy)
    
    df['CEM'] = TLOG
    df['MF'] = MF
    df['Qdot'] = Qdot

    predicted_file_dir = os.path.join(fig_dir, "predicted_Y.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close()

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header = not file_exists, index=False)



def create_reaction_path_diagram(gas, fig_dir, element, cantera_row, pred_cantera_row):
    # Prepare species compositions from cantera row
    all_species_names = gas.species_names
    all_species_values = [cantera_row[name] for name in all_species_names]
    all_fractions = dict(zip(all_species_names, all_species_values))
    species_composition = ", ".join([f"{sp}:{mf}" for sp, mf in all_fractions.items()])
    gas.TPX = cantera_row['T'], cantera_row['P'], species_composition
    diagram = ct.ReactionPathDiagram(gas, element)
    diagram.title = 'Reaction path diagram following {0}'.format(element)
    diagram.label_threshold = 0.001

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
    # Read the contents of the .dot file
    with open(dot_file_path, 'r') as file:
        dot_content = file.read()

    # Debug: Print the species names from the dot file
    print("Species found in the .dot file:")

    # Simplified regular expression to match only species labels inside "label="
    species_pattern = re.compile(r'label="([A-Za-z0-9()]+)"')

    found_species = re.findall(species_pattern, dot_content)
    for species in found_species:
        print(f"Found species in dot file: {species}")  # Print species name from the dot file

    # Normalize species names in `all_species_names` to lowercase
    all_species_names_normalized = [species.lower().strip() for species in all_species_names]

    # Print all_species_names for comparison
    print("Normalized species in all_species_names:")
    print(all_species_names_normalized)

    # Define a function to compute error
    def compute_error(species_name):
        cantera_value = cantera_row.get(species_name, 0)  # Get species value, default to 0 if not found
        pred_value = pred_cantera_row.get(species_name, 0)
        if cantera_value != 0:  # Avoid division by zero
            error = (pred_value - cantera_value) / cantera_value
        else:
            error = 0
        print(f"Species: {species_name}, Cantera value: {cantera_value}, Predicted value: {pred_value}, Error: {error}")
        return error

    def error_to_color(error):
        # Calculate absolute error percentage
        abs_error_percent = min(abs(error) * 100, 100)
        
        # Normalize error between 0 (green) and 1 (red)
        norm_error = abs_error_percent / 100
        
        # Convert to RGB: green (0,1,0) to red (1,0,0)
        red = int(255 * norm_error)
        green = int(255 * (1 - norm_error))

        # Lighten the colors by adding a constant value
        lightening_factor = 100
        red = min(red + lightening_factor, 255)
        green = min(green + lightening_factor, 255)

        color_code = f"#{red:02x}{green:02x}00"
        return color_code

    # Modify the species blocks in the .dot file
    def modify_label_and_color(match):
        species_label = match.group(1).strip()  # Extract the species label and strip any extra spaces
        species_label_lower = species_label.lower()  # Normalize to lowercase for comparison
        print(f"Checking species: {species_label}")

        # Check if the species name matches (case-insensitive check)
        if species_label_lower in all_species_names_normalized:
            # Compute the error for this species
            error = compute_error(species_label)
            if abs(error) < 0.001:  # Check if the error is less than 0.1%
                error_percentage = "set"
            else:
                error_percentage = f"{error * 100:.2f}%"
            # Map the error to a color
            color = error_to_color(abs(error))  # Use absolute error for color mapping
            # Update the label and color
            new_label = f'{species_label} {error_percentage}'
            # print(f"Modifying species: {species_label}, New label: {new_label}, Color: {color}")
            return f'label="{new_label}", style=filled, fillcolor="{color}"'

        else:
            print(f"Species '{species_label}' not found in all_species_names, skipping modification.")
            return match.group(0)  # If species not found, return the original line unchanged

    # Apply all the modifications in one pass
    modified_dot_content = re.sub(species_pattern, modify_label_and_color, dot_content)

    # Save the modified content back to the .dot file **once** after all modifications
    output_dot_file = dot_file_path.replace('.dot', '_modified.dot')
    with open(output_dot_file, 'w') as file:
        file.write(modified_dot_content)

    print(f"Modified dot file saved as {output_dot_file}")


def process_simulation_RCCE(gas,props, mech, file_name, fig_dir, main_species_names, fuel, oxidizer, modeci, dt, time_end, rand=0):
    all_species_names = gas.species_names
    color_arr = ('k','r','b','m','y','g', 'c')
    print("file_name=", file_name)
    df = pd.read_csv(file_name)
    # Find the index with maximum HRR
    for index, row in df.iterrows():
        # input_T = row['T'] * (1+rand* random.uniform(-1, 1) /5)
        input_T = row['T']
        input_P = props['P']
        input_species_values = [ row[name] for name in main_species_names]
        perturbation = [value * rand * random.uniform(-1, 1) for value in input_species_values]
        perturbed_species_values = [value + perturb for value, perturb in zip(input_species_values, perturbation)]
        initial_model_fractions = dict(zip(main_species_names, perturbed_species_values))
        print("initial_model_fractions",initial_model_fractions)
        print("sum of initial_model_fractions", sum(initial_model_fractions.values()))
        if modeci == "QSSA":
            run_CEQ_mass(gas, fig_dir, main_species_names, initial_model_fractions, fuel, oxidizer, input_T, input_P, full_data=row, dt=dt, time_end=time_end)

    # predicted_file_dir = os.path.join(fig_dir, "predicted_Y.csv")
    # print("predicted_file_dir=", predicted_file_dir)
    # elements = ['O','H','N']
    # for element in elements:
    #     create_reaction_path_diagram(gas,fig_dir, element, max_hrr_row, pred_hrr_row)
        
def create_fig_dir(file_name,prefix):
    base_name, ext = os.path.splitext(file_name)
    fig_dir = os.path.join("figs", prefix, base_name)
    os.makedirs(fig_dir, exist_ok=True)
    predicted_file_dir = os.path.join(fig_dir, "predicted_Y.csv")
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
    predicted_file_dir = os.path.join(fig_dir, "predicted_Y.csv")
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
    props = props_list[0]
    last_spec = 'AR'
    gas = ct.Solution(props['mech_file'])
    specs = gas.species()[:]
    N2_ind = gas.species_index(last_spec)
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                    species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
                    reactions=gas.reactions())
    print(f"Number of species: {gas.n_species}")
    print(f"Number of reactions: {gas.n_reactions}")
    print("species_names:", gas.species_names)
    
    for props in props_list: 
        props['perturb'] = 0
        file_name = props['file_name']
        data_dir = "data/SandiaFlame"
        file_name = os.path.join(data_dir, file_name)
        print("file_name=", file_name)
        valid_modes = {"QSSA"} # RCCE_fortran
        for modeci in valid_modes: 
            main_species_names = ["O2", "H2", "H2O", "CH4", "CO", "N2"]
            print("Reconstruction starts with method:", modeci)
            if modeci not in valid_modes:
                print(f"Warning: Invalid modeci '{modeci}'. Using default 'QSSA'.")
            if modeci == "RCCE_fortran":
                modify_streams_in_file(main_species_names, 'RCCE')
            fig_dir = create_fig_dir(file_name, modeci)
            process_simulation_RCCE(gas, props, props['mech_file'], file_name, fig_dir, main_species_names,props['fuel'], props['oxidizer'], modeci, props['dt'], props['time_end'], rand=props['perturb'])
        
if __name__ == "__main__":
    main()