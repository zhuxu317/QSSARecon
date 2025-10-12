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
# sys.path.append(os.path.abspath('src/gri30/'))
# sys.path.append(os.path.abspath('src/gri12/'))
sys.path.append(os.path.abspath('src/Li_demo/'))

from subprocess import run
from pathlib import Path
import graphviz
from lib_cema import *
import random
import glob
import matplotlib as mpl
import matplotlib.font_manager as fm

# Path to the Times New Roman font file
font_path = '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf'

# Create a font properties object using the font file
font_prop = fm.FontProperties(fname=font_path)
inch = 2.54
props = {
    'element_num':5,
    'fuel': 'CH4:1.0',         # fuel components [molar ratio]
    'oxyd': 'O2:1.0, N2:3.76', # oxydizer components [molar ratio]
}
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


def run_CEQ_core(gas, main_species_names, initial_mole_fractions, T, P, dt, time_end):
    """
    Core function to run the combustion simulation and return the updated gas object.
    
    Parameters:
    - gas: Cantera gas object
    - main_species_names: List of species to control (main species)
    - initial_mole_fractions: Dictionary of initial mole fractions for controlled species
    - T: Temperature
    - P: Pressure
    - dt: Time step for the simulation
    - time_end: End time for the simulation
    
    Returns:
    - gas: Updated Cantera gas object after the simulation with the inert species mole fraction set to 0.
    """
    # Identify the inert species (assumed to be the last species in gas)
    inert_specie = gas.species_names[-1]
    # print("inert_specie =", inert_specie)
    
    # Calculate total mole fraction of controlled species
    total_controlled = sum(initial_mole_fractions.values())
    # print(f"Sum of controlled species mole fractions: {total_controlled}")
    
    # Determine inert species mole fraction as 1 - total_controlled
    # if total_controlled < 1:
    inert_value = 1.0 - total_controlled
    # print(f"Setting {inert_specie} mole fraction to: {inert_value}")
    
    # Create species composition string including the inert species adjustment
    species_composition = ", ".join([f"{sp}:{mf}" for sp, mf in initial_mole_fractions.items()])
    species_composition += f", {inert_specie}:{inert_value}"
    
    # Set initial conditions in the gas object
    gas.TPX = T, P, species_composition
    
    # Initialize the solver
    y0 = np.hstack((gas.T, gas.Y))
    ode = ReactorOde(gas, main_species_names)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True)
    solver.set_initial_value(y0, 0.0)
    
    # Run the simulation
    while solver.successful() and solver.t < time_end:
        solver.integrate(solver.t + dt)
        gas.TPY = solver.y[0], P, solver.y[1:]
    
    # After simulation, set the inert species mole fraction to 0
    inert_index = gas.species_index(inert_specie)
    Y_new = gas.Y.copy()
    Y_new[inert_index] = 0.0
    gas.TPY = gas.T, P, Y_new
    return gas



def run_CEQ_core_mass(gas, main_species_names, initial_mole_fractions, T, P, dt, time_end):
    """
    Core function to run the combustion simulation and return the updated gas object.

    Parameters:
    - gas: Cantera gas object
    - main_species_names: List of species to control (main species)
    - initial_mole_fractions: Dictionary of initial mole fractions for controlled species
    - T: Temperature
    - P: Pressure
    - dt: Time step for the simulation
    - time_end: End time for the simulation

    Returns:
    - gas: Updated Cantera gas object after the simulation
    """
    # Set initial conditions in the gas object
    species_composition = ", ".join([f"{sp}:{mf}" for sp, mf in initial_mole_fractions.items()])
    gas.TPY = T, P, species_composition

    # Initialize the solver
    y0 = np.hstack((gas.T, gas.Y))
    ode = ReactorOde(gas, main_species_names)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True)
    solver.set_initial_value(y0, 0.0)

    # Run the simulation
    while solver.successful() and solver.t < time_end:
        solver.integrate(solver.t + dt)
        gas.TPY = solver.y[0], P, solver.y[1:]

    ## STOP here and return gas
    return gas



def run_CEQ(gas, fig_dir, main_species_names, initial_mole_fractions, fuel, oxyd, T, P, dt, time_end, full_data):
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
    gas = run_CEQ_core(gas, main_species_names, initial_mole_fractions, T, P, dt, time_end)

    # Extract mole fractions and species names
    species_names = gas.species_names
    mole_fractions = gas.X
    
    # Create a DataFrame for the results
    data = {species: y for species, y in zip(species_names, mole_fractions)}
    
    updated_data = {species: initial_mole_fractions.get(species, y) for species, y in data.items()}
    
    df = pd.DataFrame(updated_data, index=[0])

    # Perform additional calculations
    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num=props['element_num'])
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(fuel, oxyd)
    Qdot = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy)

    # Add calculated values to the DataFrame
    df['CEM'] = TLOG
    df['MF'] = MF
    df['Qdot'] = Qdot
    if 't_res' in full_data:
        df['t_res'] = full_data['t_res']

    # Save the results to a CSV file
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close()

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header=not file_exists, index=False)
        
def run_CEQ_mass(gas, fig_dir, main_species_names, initial_mass_fractions, fuel, oxyd, T, P, dt, time_end, full_data):
    """
    Run the combustion simulation, extract data, and save results to a file.

    Parameters:
    - gas: Cantera gas object
    - fig_dir: Directory for saving results
    - main_species_names: List of species to control (main species)
    - initial_mass_fractions: Dictionary of initial mole fractions for controlled species
    - fuel: Fuel species name
    - oxyd: Oxidizer species name
    - T: Temperature
    - P: Pressure
    - dt: Time step for the simulation
    - time_end: End time for the simulation
    - full_data: Additional data for processing
    """
    # Call the core simulation function
    gas = run_CEQ_core_mass(gas, main_species_names, initial_mass_fractions, T, P, dt, time_end)

    # Extract mole fractions and species names
    species_names = gas.species_names
    mass_fractions = gas.Y

    # Create a DataFrame for the results
    data = {species: y for species, y in zip(species_names, mass_fractions)}
    df = pd.DataFrame(data, index=[0])

    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num=props['element_num'])
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(fuel, oxyd)
    Qdot = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy)

    # Add calculated values to the DataFrame
    df['CEM'] = TLOG
    df['MF'] = MF
    df['Qdot'] = Qdot
    if 't_res' in full_data:
        df['t_res'] = full_data['t_res']
    if 'z' in full_data:
        df['z'] = full_data['z']
    if 'r/d' in full_data:
        df['r/d'] = full_data['r/d']
    # Save the results to a CSV file
    predicted_file_dir = os.path.join(fig_dir, "predicted_Y.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close()

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header=not file_exists, index=False)
        


def create_reaction_path_diagram(gas, fig_dir, element, cantera_row, pred_cantera_row):
    # Prepare species compositions from cantera row
    all_species_names = gas.species_names
    all_species_values = [cantera_row[name] for name in all_species_names]
    all_fractions = dict(zip(all_species_names, all_species_values))
    species_composition = ", ".join([f"{sp}:{mf}" for sp, mf in all_fractions.items()])
    gas.TPY = cantera_row['T'], cantera_row['P'], species_composition
    diagram = ct.ReactionPathDiagram(gas, element)
    diagram.title = 'Reaction path diagram following {0}'.format(element)
    diagram.label_threshold = 0.01

    # Define the output directory and files
    output_dir = fig_dir
    dot_file = f'{element}_rxnpath.dot'
    modified_dot_file = f'{element}_rxnpath_modified.dot'
    img_file = f'{element}_rxnpath.png'

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
    run(['dot', '-Tpng', dot_path, '-o', img_path, '-Gdpi=200'], check=True)
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


def process_simulation(gas, file_name, fig_dir, main_species_names, x_range, rand=0):
    color_arr = ('k','r','b','m','y','g', 'c')
    df = pd.read_csv(file_name)
    # Find the index with maximum HRR
    max_hrr_index = None
    max_hrr_value = 0
    for index, row in df.iterrows():
        input_T = row['T']
        input_P = row['P']
        HRR = row['Qdot']
        input_species_values = [row[name] for name in main_species_names]
        # Add 5% random perturbation to input_species_values
        perturbation = [value * rand * random.uniform(-1, 1) for value in input_species_values]
        perturbed_species_values = [value + perturb for value, perturb in zip(input_species_values, perturbation)]
        initial_model_fractions = dict(zip(main_species_names, perturbed_species_values))
        run_CEQ(gas, fig_dir, main_species_names, initial_model_fractions, input_T, input_P, dt=1e-2, time_end=1e-1)
        
        if HRR > max_hrr_value:
            max_hrr_value = HRR
            max_hrr_index = index
    print("Index for maximum HRR:", max_hrr_index)
    predicted_file_dir = os.path.join(fig_dir, "predicted_Y.csv")
    df_pred = pd.read_csv(predicted_file_dir)
    # Visualize reaction path for largest HRR
    max_hrr_row = df.loc[max_hrr_index]
    pred_hrr_row = df_pred.loc[max_hrr_index]
    

    HRR_error_max = (df_pred['Qdot'].max() - df['Qdot'].max()) /  df['Qdot'].max() * 100
    HRR_error_L1 = np.sum(np.abs(df_pred['Qdot'] - df['Qdot'])) / np.sum(np.abs(df['Qdot'])) * 100
    
    elements = ['O','C', 'H']
    for element in elements:
        create_reaction_path_diagram(gas,fig_dir, element, max_hrr_row, pred_hrr_row)
    # Plot data for CH4, O2, H2O
    normalized_grid = df['normalized_grid']
    CH4_df, O2_df, H2O_df = df['CH4'], df['O2'], df['H2O']
    CH4_df_pred, O2_df_pred, H2O_df_pred = df_pred['CH4'], df_pred['O2'], df_pred['H2O']

    plt.figure(figsize=(6, 5))
    plt.scatter(normalized_grid[::5], CH4_df[::5], label='CH4', color=color_arr[0])
    plt.scatter(normalized_grid[::5], O2_df[::5], label='O2', color=color_arr[1])
    plt.scatter(normalized_grid[::5], H2O_df[::5], label='H2O', color=color_arr[2])

    plt.plot(normalized_grid, CH4_df_pred, label='CH4 GFRI', linestyle='-.', color=color_arr[0])
    plt.plot(normalized_grid, O2_df_pred, label='O2 GFRI', linestyle='-.', color=color_arr[1])
    plt.plot(normalized_grid, H2O_df_pred, label='H2O GFRI', linestyle='-.', color=color_arr[2])
    plt.legend()
    plt.xlabel('x', fontsize=14)
    plt.ylabel('Y', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_main_species_Y.png"))

    # Plot data for OH, CH2O, HCO
    OH_df, CH2O_df, HCO_df = df['OH'], df['CH2O'], df['HCO']
    OH_df_pred, CH2O_df_pred, HCO_df_pred = df_pred['OH'], df_pred['CH2O'], df_pred['HCO']

    fig, axs = plt.subplots(3, figsize=(15.24/inch, 12/inch), sharex=True)
    for ax in axs:
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(0.8)

    axs[0].scatter(normalized_grid[::3], OH_df[::3], label='OH', color=color_arr[0])
    axs[0].plot(normalized_grid, OH_df_pred, label='OH predicted', linestyle='-', color=color_arr[0])
    axs[0].legend(fontsize=8)
    axs[0].set_ylabel('Y', fontsize=10)
    axs[0].set_xlim(x_range)

    axs[1].scatter(normalized_grid[::5], CH2O_df[::5], label='CH2O', color=color_arr[1])
    axs[1].plot(normalized_grid, CH2O_df_pred, label='CH2O predicted', linestyle='-', color=color_arr[1])
    axs[1].legend(fontsize=8)
    axs[1].set_ylabel('Y', fontsize=10)
    axs[1].set_xlim(x_range)

    axs[2].scatter(normalized_grid[::5], HCO_df[::5], label='HCO', color=color_arr[2])
    axs[2].plot(normalized_grid, HCO_df_pred, label='HCO predicted', linestyle='-', color=color_arr[2])
    axs[2].legend(fontsize=8)
    axs[2].set_xlabel('x', fontsize=10)
    axs[2].set_ylabel('Y', fontsize=10)
    axs[2].set_xlim(x_range)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_predict_species_Y.png"))

    # Plot for CEM, MF, HRR (Qdot)
    CEM_df, MF_df, Qdot_df = df['CEM'], df['MF'], df['Qdot']
    CEM_df_pred, MF_df_pred, Qdot_df_pred = df_pred['CEM'], df_pred['MF'], df_pred['Qdot']

    fig, axs = plt.subplots(3, figsize=(8/inch, 12/inch), sharex=True)
    for ax in axs:
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(0.8)

    axs[0].scatter(normalized_grid[::2], CEM_df[::2], label='CEM', color=color_arr[1])
    axs[0].plot(normalized_grid, CEM_df_pred, label='CEM predicted', linestyle='-', color=color_arr[1])
    axs[0].legend(fontsize=8)
    axs[0].set_ylabel('CEM', fontsize=10)
    axs[0].set_xlim(x_range)

    axs[1].scatter(normalized_grid[::2], MF_df[::2], label='MF', color=color_arr[2])
    axs[1].plot(normalized_grid, MF_df_pred, label='MF predicted', linestyle='-', color=color_arr[2])
    axs[1].legend(fontsize=8)
    axs[1].set_ylabel('MF', fontsize=10)
    axs[1].set_xlim(x_range)

    axs[2].scatter(normalized_grid[::2], Qdot_df[::2], label='HRR', color=color_arr[3])
    axs[2].plot(normalized_grid, Qdot_df_pred, label='HRR predicted', linestyle='-', color=color_arr[3])
    axs[2].legend(fontsize=8)
    axs[2].set_xlabel('x', fontsize=10)
    axs[2].set_ylabel('HRR', fontsize=10)
    axs[2].set_xlim(x_range)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_OH_CEM_MF.png"))
    
    fig, axs = plt.subplots(3, 2, figsize=(15/inch, 15/inch), sharex=True)
    # Plot for OH, CH2O, HCO
    OH_df, CH2O_df, HCO_df = df['OH'], df['CH2O'], df['HCO']
    OH_df_pred, CH2O_df_pred, HCO_df_pred = df_pred['OH'], df_pred['CH2O'], df_pred['HCO']
    # Plot for CEM, MF, HRR (Qdot)
    CEM_df, MF_df, Qdot_df = df['CEM'], df['MF'], df['Qdot']
    CEM_df_pred, MF_df_pred, Qdot_df_pred = df_pred['CEM'], df_pred['MF'], df_pred['Qdot']
    # Define a list of data for easy iteration
    data_list = [(OH_df, OH_df_pred, 'OH', 0), (CH2O_df, CH2O_df_pred, 'CH2O', 1), 
                (HCO_df, HCO_df_pred, 'HCO', 2), (CEM_df, CEM_df_pred, 'CEM', 0), 
                (MF_df, MF_df_pred, 'MF', 1), (Qdot_df, Qdot_df_pred, 'HRR', 2)]
    for i, (df, df_pred, label, color_index) in enumerate(data_list):
        row = i // 2
        col = i % 2
        ax = axs[row, col]
        ax.scatter(normalized_grid[::2], df[::2], label=label, color=color_arr[color_index])
        ax.plot(normalized_grid, df_pred, label=f'{label} predicted', linestyle='-', color=color_arr[color_index])
        ax.legend(fontsize=8)
        ax.set_ylabel(label, fontsize=10)
        ax.set_xlim(x_range)
        if row == 2:
            ax.set_xlabel('x', fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_combine.png"))

    return HRR_error_max, HRR_error_L1
    
    
def create_fig_dir(file_name):
    base_name, ext = os.path.splitext(file_name)
    fig_dir = os.path.join("figs", base_name)
    os.makedirs(fig_dir, exist_ok=True)
    predicted_file_dir = os.path.join(fig_dir, "predicted_Y.csv")
    if os.path.isfile(predicted_file_dir):
        os.remove(predicted_file_dir)
        
    return fig_dir
    
def main():
    mech_file = "src/gri12/gri12.cti"
    last_spec = 'AR'
    gas = ct.Solution(mech_file)
    specs = gas.species()[:]
    N2_ind = gas.species_index(last_spec)
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                    species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
                    reactions=gas.reactions())
    print(f"Number of species: {gas.n_species}")
    print(f"Number of reactions: {gas.n_reactions}")
    print("species_names:", gas.species_names)
    perturb = 0
    x_range = [0,1]
    # extra_species = [None]
    extra_species = [None, "H", "O", "OH", "HO2", "H2O2", "C", "CH", "CH2", "CH2(S)", "CH3", "HCO", "CH2O", "CH2OH", "CH3O", "CH3OH", "C2H", "C2H2", "C2H3", "C2H4", "C2H5", "C2H6", "HCCO", "CH2CO", "HCCOH"]
    for extra_specie in extra_species:
        main_species_names = ["H2", "O2", "CO", "CO2", "CH4", "H2O", "N2"]
        if extra_specie is None:
            extra_specie = "original"
        else:
            main_species_names.extend([extra_specie])
        # If you want to process the dir csv files
        dir_path = "data/case_CH4_counterflow_premixed/phi1.0"
        csv_files = glob.glob(os.path.join(dir_path, "*.csv"))
        data = []
        for file_name in csv_files:
            csv_name = os.path.basename(file_name)
            fig_dir = create_fig_dir(file_name)
            parts = csv_name.split("_")
            strain_rate = float(parts[2])
            HRR_error_max, HRR_error_L1 = process_simulation(gas, file_name, fig_dir, main_species_names, x_range, rand=perturb)
            data.append([strain_rate, HRR_error_max, HRR_error_L1])
        df = pd.DataFrame(data, columns=['strain_rate', 'HRR_error_max', 'HRR_error_L1'])
        df.to_csv(os.path.join(os.path.join('figs',dir_path), f'{extra_specie}_HRR_data.csv'), index=False)
        
        
    # # If you need only process 1 single file
    # file_name = "data/case_CH4_counterflow_premixed/phi1.0/strain_rate_2.2326_consumption_speed_705.8656.csv"
    # fig_dir = create_fig_dir(file_name)
    # HRR_error_max, HRR_error_L1 = process_simulation(gas, file_name, fig_dir, main_species_names, x_range, rand=perturb)
    # print("HRR_error_max=",HRR_error_max)
    
if __name__ == "__main__":
    main()