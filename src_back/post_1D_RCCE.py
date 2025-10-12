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

from props_1D import AVAILABLE_CASES

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
    parser.add_argument(
        '--inert_specie',
        type=str,
        required=True,
        help="Inert species to be used in the mechanism (e.g., 'AR', 'N2')"
    )
    return parser.parse_args()


# def parse_args():
#     parser = argparse.ArgumentParser(description="Run 1D Counterflow Flame simulation")
#     # Existing case_name argument
#     parser.add_argument(
#         '--case_name', 
#         type=str, 
#         choices=['NH3_CTM', 'NH3_NP_sim', 'NH3_NP_exp','CH4'], 
#         required=True, 
#         help="Specify the case name: 'NH3_CTM' or 'NH3_NP' or 'CH4'"f
#     )
#     return parser.parse_args()


args = parse_args()
case_name = args.case_name
mech_file = args.mech_file
element_num = args.element_num
inert_specie = args.inert_specie    # get folder and full filename
mech_folder = os.path.dirname(mech_file)          # "mechanism/NH3_otomo"
mech_filename = os.path.basename(mech_file)        # "NH3_otomo.cti"
mech_name, _ = os.path.splitext(mech_filename)     # ("NH3_otomo", ".cti")

from props_1D import generate_props  # Import the function to generate props
props_list = generate_props(case_name)
sys.path.append(os.path.abspath(mech_folder))
from RCCE import RCCE
from lib_cema import *

from subprocess import run
from pathlib import Path
import graphviz
import random
import glob
import matplotlib as mpl
import matplotlib.font_manager as fm
import subprocess
# from run_1D_CEQ import run_CEQ
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
    # solver.set_integrator('vode', method='bdf', with_jacobian=True)
    solver.set_integrator(
        'vode',
        method='bdf',
        with_jacobian=True,
        # rtol=1e-6,        # relative tolerance
        # atol=1e-10,       # absolute tolerance
        # max_step=dt       # never take a step bigger than your desired dt
    )

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
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    
    df['CEM'] = TLOG
    df['MF'] = MF
    df['HRR'] = HRR

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close()

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header = not file_exists, index=False)



def run_RCCE(gas, rcce, fig_dir, props):
    """
    Run the simulation to monitor the concentrations of controlled and target species.

    Parameters:
    - time_step: Time step for the simulation
    - time_end: End time for the simulation
    """
            
    (Tsim, Psim, Xsim), History = rcce.equilibrium(method="TP")
    gas.TPX = Tsim, Psim, Xsim
    mole_fractions = gas.X
    data = {species: y for species, y in zip(gas.species_names, mole_fractions)}
    df = pd.DataFrame(data, index=[0])

    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num= props['element_num'])
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(props['fuel'], props['oxidizer'])
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    df['CEM'] = TLOG
    df['MF'] = MF
    df['HRR'] = HRR

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close() 

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header = not file_exists, index=False)


def run_RCCE_fortran(gas, fig_dir, T, P, modeci,props):
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
        df = pd.read_csv(output_file_path, delim_whitespace=True, header=None)

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
    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num= props['element_num'])
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(props['fuel'], props['oxidizer'])
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    
    df['CEM'] = TLOG
    df['MF'] = MF
    df['HRR'] = HRR

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close() 

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header = not file_exists, index=False)


def run_RCEQ(gas, rceq, fig_dir, props):
    """
    Run the simulation to monitor the concentrations of controlled and target species.

    Parameters:
    - time_step: Time step for the simulation
    - time_end: End time for the simulation
    """
            
    (Tsim, Psim, Xsim), History = rceq.equilibrium(method="TP")
    gas.TPX = Tsim, Psim, Xsim
    mole_fractions = gas.X
    data = {species: y for species, y in zip(gas.species_names, mole_fractions)}
    df = pd.DataFrame(data, index=[0])

    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num= props['element_num'])
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(props['fuel'], props['oxidizer'])
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    df['CEM'] = TLOG
    df['MF'] = MF
    df['HRR'] = HRR

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    file_exists = os.path.isfile(predicted_file_dir)

    if not file_exists:
        open(predicted_file_dir, 'w').close() 

    # Append to the file without writing the header if it already exists
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header = not file_exists, index=False)
        
        
def run_CEQ(gas, fig_dir, main_species_names, initial_mole_fractions, fuel, oxyd, T, P, dt, time_end, full_data, element_num, moleFraction=True):
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
    gas = run_CEQ_core(gas, main_species_names, initial_mole_fractions, T, P, dt, time_end, moleFraction=True)

    # Extract mole fractions and species names
    species_names = gas.species_names
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
            updated_data[species] = initial_mole_fractions.get(base_species, value)
        else:
            # Keys without the prefix are copied unchanged
            updated_data[species] = value
    df = pd.DataFrame(updated_data, index=[0])

    # Perform additional calculations
    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num=element_num)
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(fuel, oxyd)
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
    initial_mole_fractions,
    T,
    P,
    dt,
    time_end,
    moleFraction=True,
    steady_time=0.1,
    steady_tol=1e-12
):
    """
    Core function to run the combustion simulation and return the updated gas object.
    
    Stops early if the solution reaches steady state (no change > steady_tol after steady_time).
    """
    # 1) Identify inert and build composition vector
    inert_specie    = gas.species_names[-1]
    total_controlled = sum(initial_mole_fractions.values())
    inert_value      = max(0.0, 1.0 - total_controlled)

    X = np.zeros(gas.n_species)
    for sp, mf in initial_mole_fractions.items():
        X[gas.species_index(sp)] = mf
    if inert_value > 0.0:
        X[gas.species_index(inert_specie)] = inert_value

    # 2) Apply to gas
    gas.TP = T, P
    if moleFraction:
        gas.set_unnormalized_mole_fractions(X)
    else:
        mw = gas.molecular_weights
        Y = X * mw / np.sum(X * mw)
        gas.set_unnormalized_mass_fractions(Y)

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
    Y_final[gas.species_index(inert_specie)] = 0.0
    gas.TPY = gas.T, P, Y_final
    return gas

def run_CEQ_core_old(gas, main_species_names, initial_mole_fractions, T, P, dt, time_end, moleFraction=True):
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
    
    # Calculate total mole fraction of controlled species
    total_controlled = sum(initial_mole_fractions.values())
    # print(f"Sum of controlled species mole fractions: {total_controlled}")
    
    # amount of inert we must add (≥0)
    inert_value = max(0.0, 1.0 - total_controlled)

    # build the composition string
    species_composition = ", ".join(
        f"{sp}:{mf}" for sp, mf in initial_mole_fractions.items()
    )
    if inert_value > 0:
        species_composition += f", {inert_specie}:{inert_value}"

    # ----‑‑‑ diagnostic print ------------------------------------------------
    total_feed = total_controlled + inert_value
    # print(f"Total mole fraction sent to Cantera: {total_feed:.6g}")
    
    X = np.zeros(gas.n_species)        # will hold unnormalised amount
    for sp, mf in initial_mole_fractions.items():
        X[gas.species_index(sp)] = mf
    if inert_value > 0.0:
        X[gas.species_index(inert_specie)] = inert_value

    # (If you start with mass fractions, assemble Y instead of X.)

    # ------------------------------------------------------------------ apply to Cantera – no normalisation
    gas.TP = T, P                      # set thermo state first

    if moleFraction:
        gas.set_unnormalized_mole_fractions(X)
    else:
        mw = gas.molecular_weights
        Y = X * mw / np.sum(X * mw)
        gas.set_unnormalized_mass_fractions(Y)
        
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


def run_ILDM(gas, ildm,Tsim, Psim, X_ini, fig_dir, props):
    """
    Run the simulation to monitor the concentrations of controlled and target species.

    Parameters:
    - time_step: Time step for the simulation
    - time_end: End time for the simulation
    """
    # ildm.gas.equilibrate('HP')
    N_ini = ildm.get_N()
    Hini = ildm.gas.enthalpy_mass

    Xsim = ildm.solve_dynamics_pseudo_time(
                dt= 1e-7,###5.0e-7,     # Adjust as necessary
                N0 = X_ini,
                max_steps=5000,
                tol=1e-3,
                min_steps=2000,        # Ensure at least 50 iterations
                verbose=False
            )
    
    # Xsim = N_ini /N_ini.sum()
    
    # Xsim = ildm.solve_dynamics_pseudo_time_adaptive(
    #     N0=X_ini,
    #     Q_L=None,
    #     dt_initial=1e-2,
    #     max_steps=1000,
    #     tol=1e-6,
    #     min_steps=100,        # Minimum iterations before convergence check
    #     max_dt_factor=2.0,     # Maximum factor to increase dt
    #     min_dt_factor=0.5,     # Maximum factor to decrease dt
    #     verbose=False
    #     )
    
    gas.HPX = Hini, Psim, Xsim
    # gas.TPX = Tsim, Psim, Xsim
    mole_fractions = gas.X
    data = {species: y for species, y in zip(gas.species_names, mole_fractions)}
    df = pd.DataFrame(data, index=[0])

    D, L, R = solve_eig_gas(gas)
    eig_vals = highest_val_excl_0(D, N_val=4, element_num= props['element_num'])
    CEM = eig_vals[-1]
    CEM_real = np.real(CEM)
    sign_CEM = np.sign(CEM_real)
    log_CEM = np.log10(1 + np.abs(CEM_real))
    TLOG = sign_CEM * log_CEM
    MF = gas.mixture_fraction(props['fuel'], props['oxidizer'])
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    df['CEM'] = TLOG
    df['MF'] = MF
    df['HRR'] = HRR

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
    elif modeci == "ILDM":
        ildm = ILDM(mech, main_species_names, use_element_constraints=True)

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
            run_CEQ(
                gas, fig_dir, main_species_names, initial_mf,
                fuel, oxidizer, input_T, input_P,
                full_data=row, dt=dt, time_end=time_end,
                element_num=props["element_num"], moleFraction=mole_fraction
            )

        elif modeci == "RCCE":
            gas.TPX = input_T, input_P, mole_fractions
            rcce.set_conditions(*gas.TPX)
            run_RCCE(gas, rcce, fig_dir, props)

        elif modeci == "RCCE_fortran":
            gas.TPX = input_T, input_P, mole_fractions
            # mole_fractions["NNH"] = 0.0   # explicit zero
            run_RCCE_fortran(gas, fig_dir, input_T, input_P, modeci, props)

        elif modeci == "ILDM":
            gas.TPX = input_T, input_P, mole_fractions
            ildm.set_conditions(*gas.TPX)
            run_ILDM(gas, ildm, fig_dir, props)

        else:
            raise ValueError(f"Unknown modeci: {modeci}")


def process_simulation_RCCE(gas,props, mech, file_name, fig_dir, main_species_names, fuel, oxidizer, modeci, dt, time_end,denoted_value_name, rand=0, mole_fraction=True):
    # -------------------------------------------------- data & utilities
    all_species_names = gas.species_names
    df_main           = pd.read_csv(file_name)
    # optional assist file
    df_assist = None
    if props.get("assist_path") and assist_species:
        print("assist path is read here!")
        assist_species = props['species_from_exp']
        assist_file = os.path.join(props["assist_path"], props["file_name"])
        df_assist   = pd.read_csv(assist_file)
    # -------------------------------------------------- mechanism objects
    if modeci == "RCCE":
        rcce = RCCE(mech, main_species_names, use_element_constraints=True)
    elif modeci == "ILDM":
        ildm = ILDM(mech, main_species_names, use_element_constraints=True)

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
        # HRR = row['HRR']
        
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
            all_species_values = [row[name] if name not in main_species_names else initial_model_fractions.get(name, row[name]) 
                                for name in all_species_names]
            
            mole_fractions = dict(zip(all_species_names, all_species_values))
            
        if modeci == "QSSA":
            run_CEQ(gas, fig_dir, main_species_names, initial_model_fractions, fuel, oxidizer, input_T, input_P, full_data=row, dt=dt, time_end=time_end, element_num=props['element_num'], moleFraction=props["ifMoleFraction"])
            
        if modeci == "RCCE":
            gas.TPX = input_T, input_P, mole_fractions
            input_T, input_P, mole_fractions = gas.TPX
            rcce.set_conditions(input_T, input_P, mole_fractions)
            run_RCCE(gas, rcce, fig_dir, props)
        if modeci == "RCCE_fortran":
            gas.TPX = input_T, input_P, mole_fractions
            mole_fractions['NNH'] = 0
            # print("mole_fractions=", mole_fractions)
            run_RCCE_fortran(gas, fig_dir, input_T, input_P, modeci,props)
        if modeci == "RCEQ":
            all_species_names= []
            lif_species_names = ["OH", "NO"]
            all_species_names = main_species_names + lif_species_names
            # Update mole_fractions: keep values for species in all_species_names, set others to 0
            for species in list(mole_fractions.keys()):
                if species not in all_species_names:
                    mole_fractions[species] = 0
            rceq.set_conditions(input_T, input_P, mole_fractions)
            run_RCEQ(gas, rceq, fig_dir, props)
        if modeci == "ILDM":
                ############ RUN QSSA#########################
            gas = run_CEQ_core(gas, main_species_names, initial_model_fractions, input_T, input_P, dt=1e-4, time_end=1e-3)
            # X_qssa = gas.X  # Mole fractions from QSSA
            X_qssa = np.array(all_species_values)
            try:
                ildm.initialize_conditions(input_T, input_P, mole_fractions)
                ildm.set_constrained_species_fractions(initial_model_fractions)
                run_ILDM(gas, ildm, input_T, input_P, X_qssa, fig_dir, props)
            except Exception as e:
                print(f"ILDM run failed: {e}")
                print("Falling back to CEQ method")
                run_CEQ(gas, fig_dir, main_species_names, initial_model_fractions, fuel, oxidizer, input_T, input_P, full_data=row, dt=dt, time_end=time_end)

    print("Index for maximum HRR:", max_T_index)
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    df_pred = pd.read_csv(predicted_file_dir)

    # # 5) compute the absolute deviation series:
    # if denoted_value_name not in df_pred.columns:
    #     raise KeyError(f"{denoted_value_name} not found in predicted_X.csv")
    # deviations = (df_main[denoted_value_name] - df_pred[denoted_value_name]).abs()

    # # 6) find the index of the largest deviation
    # max_dev_index = deviations.idxmax()
    # print(f"Row with largest |Δ{denoted_value_name}| is {max_dev_index}, deviation = {deviations[max_dev_index]:.3e}")

    # # 7) grab the two rows
    # max_row_sim  = df_main.loc[max_dev_index]
    # max_row_pred = df_pred.loc[max_dev_index]

    # # 8) now build the reaction‐path diagrams
    # elements = ['O','H','N']
    # for element in elements:
    #     create_reaction_path_diagram(
    #         gas, fig_dir, element,
    #         max_row_sim,    # simulation row
    #         max_row_pred    # prediction row
    #     )
        
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

    gas = ct.Solution(mech_file)
    specs = gas.species()[:]
    N2_ind = gas.species_index(inert_specie)
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                    species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
                    reactions=gas.reactions())
    print(f"Number of species: {gas.n_species}")
    print(f"Number of reactions: {gas.n_reactions}")
    print("species_names:", gas.species_names)
    
    for props in props_list: 
        props['mech_file'] = mech_file
        props['mech_path'] = mech_folder
        props['mech_name'] = mech_name
        props['element_num'] = int(element_num)
        props['inert_specie'] = inert_specie
        props['perturb'] = 0
        
        file_name = props['file_name']
        data_dir = props['path']
        
        file_name = os.path.join(data_dir, file_name)
        print("file_name=", file_name)
        # valid_modes = {"RCCE_fortran", "QSSA"} # RCCE_fortran
        valid_modes = { "QSSA"} # RCCE_fortran
        
        for modeci in valid_modes: 
            # main_species_names = ["NH3","H2", "O2", "H2O", "N2","OH", "NH2", "NNH"]
            main_species_names = [ "H2", "O2", "H2O", "N2"]
            extra_name = None
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