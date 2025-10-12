# -*- coding: utf-8 -*-
"""
Utilize the CEQ method to predict PaSR manifold.
Created on Mon 4 November 2024 15:51:29 PM CST
@author: Xu Zhu
"""

import sys
import os
import numpy as np
import cantera as ct
import pandas as pd

from props_PaSR import *
sys.path.append(os.path.abspath(props['mech_path']))
from lib_cema import *
from run_1D_CEQ import run_CEQ
from post_1D_RCCE import modify_streams_in_file, create_fig_dir
from post_1D_RCCE import process_simulation_RCCE
    

inch = 2.4

def npy_to_csv(npy_file, csv_file, species_names):
    # Load the .npy file
    data = np.load(npy_file)
    
    # Extract the row data from -6 to -1 (inclusive)
    rows = data[-5:, :, :]
    
    # Initialize an empty DataFrame
    df = pd.DataFrame()
    
    # Create column headers
    column_headers = ['t_res', 'T', 'P'] + species_names
    
    # Concatenate the rows into the DataFrame
    for row in rows:
        temp_df = pd.DataFrame(row, columns=column_headers)
        df = pd.concat([df, temp_df], ignore_index=True)
    
    # Save the DataFrame to a .csv file
    df.to_csv(csv_file, index=False)
    

def csv_add_CEM_MF_Qdot(gas, csv_file, output_file, fuel, oxidizer):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)

    # Prepare to store CEM, MF, and Qdot values
    CEM_list = []
    MF_list = []
    Qdot_list = []

    for index, row in df.iterrows():
        input_T = row['T']
        input_P = row['P']
        all_species_values = [row[name] for name in gas.species_names]
        mole_fractions = dict(zip(gas.species_names, all_species_values))
        # Set the state of the gas
        gas.TPX = input_T, input_P, mole_fractions
        # Calculate CEM
        D, L, R = solve_eig_gas(gas)
        eig_vals = highest_val_excl_0(D, N_val=4, element_num=5)
        CEM = eig_vals[-1]
        CEM_real = np.real(CEM)
        sign_CEM = np.sign(CEM_real)
        log_CEM = np.log10(1 + np.abs(CEM_real))
        TLOG = sign_CEM * log_CEM
        # Calculate MF and Qdot
        MF = gas.mixture_fraction(fuel, oxidizer)
        Qdot = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy)

        # Append the calculated values to the lists
        CEM_list.append(TLOG)
        MF_list.append(MF)
        Qdot_list.append(Qdot)
    
    # Add new columns to the DataFrame
    df['CEM'] = CEM_list
    df['MF'] = MF_list
    df['Qdot'] = Qdot_list

    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_file, index=False)
    print(f"Updated CSV file saved as: {output_file}")

    
def main():
    target_case = props['target_case']
    npy_file = os.path.join('data/PaSR', target_case + '.npy')
    csv_file = os.path.join('data/PaSR', target_case + '.csv')
    output_file_name = os.path.join('data/PaSR', target_case + 'CEM.csv')
    
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
    npy_to_csv(npy_file, csv_file, gas.species_names)
    csv_add_CEM_MF_Qdot(gas, csv_file, output_file_name, props['fuel'], props['oxidizer'])
    perturb = 0
    
    file_name = output_file_name
    valid_modes = { "QSSA", "RCCE_fortran"}
    # valid_modes = { "RCCE"}
    for modeci in valid_modes: 
        # main_species_names = ["H2", "O2", "CO", "CO2", "CH4", "H2O"] 
        main_species_names = ["NH3", "O2", "H2", "H2O", "N2", "AR"]
        
        print("Reconstruction starts with method:", modeci)
        if modeci not in valid_modes:
            print(f"Warning: Invalid modeci '{modeci}'. Using default 'QSSA'.")
            modeci = "QSSA"
        fig_dir = create_fig_dir(file_name, modeci)
        
        process_simulation_RCCE(gas, props, props['mech_file'], file_name, fig_dir, main_species_names,props['fuel'], props['oxidizer'], modeci, props['dt'], props['time_end'], rand=props['perturb'])
    
        
if __name__ == "__main__":
    main()
    

        