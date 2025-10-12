import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import pandas as pd
import time as t
from props_0D import *
sys.path.append(os.path.abspath(props_list[0]['mech_path']))
from lib_cema import *

# Define the simulation function
def run_simulation(props):
    # Initialize the gas object
    gas = ct.Solution(props['mech_file'])
    
    specs = gas.species()[:]
    N2_ind = gas.species_index(props['last_spec'])
    
    # Modify gas solution by removing the last species and adding it at the end
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                      species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
                      reactions=gas.reactions())

    # Create the list of species names
    EI_keys = [''] * gas.n_species
    EI_keys[0] = 'T'
    for i in range(1, gas.n_species):
        EI_keys[i] = gas.species_name(i - 1)

    # Set initial conditions
    gas.set_equivalence_ratio(props['phi'], props['fuel'], props['oxyd'])
    gas.TP = props['T'], props['P']
    
    # Create a constant pressure reactor
    r = ct.IdealGasConstPressureReactor(gas)

    # Create the reactor network
    sim = ct.ReactorNet([r])

    # Initialize time and data vectors
    time = 0.0
    tim = []
    temp = []
    pressure = []  # To store pressure values
    eigenvalues = []
    CEM = []
    Qdot = []
    X_values = []  # To store gas.Y values
    dt_values = []  # Store time step changes
    previous_time = 0.0

    # Start the simulation
    start = t.time()
    current_time = 0
    while current_time < props['simulation_time']:
        sim.step()
        current_time = sim.time
        tim.append(current_time)
        temp.append(r.T)

        # Calculate the time step (difference between current and previous time)
        dt = current_time - previous_time
        dt_values.append(dt)
        previous_time = current_time

        D, L, R = solve_eig_gas(gas)
        eig_vals = highest_val_excl_0(D, N_val=props['N_eig'], element_num=props['element_num'])
        eigenvalues.append(eig_vals)
        cem_value = eig_vals[-1]
        CEM_real = np.real(cem_value)
        sign_CEM = np.sign(CEM_real)
        log_CEM = np.log10(1 + np.abs(CEM_real))
        TLOG = sign_CEM * log_CEM
        CEM.append(TLOG)
    
        pressure.append(gas.P)  # Append pressure at each step
        Qdot.append(-1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy))
        X_values.append(gas.X)  # Save gas.X values at each timestep

    end = t.time()

    print(f"Simulation time: {end - start} seconds")

    # Convert lists to numpy arrays for easier manipulation and plotting
    tim = np.array(tim)
    temp = np.array(temp)
    pressure = np.array(pressure)  # Convert pressure list to numpy array
    eigenvalues = np.array(eigenvalues).T
    CEM = np.array(CEM)
    Qdot = np.array(Qdot)
    X_values = np.array(X_values)

    # Save the results to a CSV file
    df = pd.DataFrame({
        't': tim,
        'T': temp,
        'P': pressure,  # Store pressure values in the CSV
        'CEM': CEM,
        'Qdot': Qdot
    })

    # Add species mole fractions (gas.X) to DataFrame
    gas_species = [spec.name for spec in gas.species()]
    for i, species in enumerate(gas_species):
        df[species] = X_values[:, i]

    # Save DataFrame to CSV
    df.to_csv(os.path.join(props['save_dir'], props['file_name']), index=False)
    print(f"Simulation results saved with {props['file_name']}")

def main():
    for i, props in enumerate(props_list):
        print(f"Running simulation for configuration {i + 1}")
        run_simulation(props)

if __name__ == "__main__":
    main()
