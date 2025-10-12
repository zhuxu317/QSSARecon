import cantera as ct
import os
# Define the temperatures and equivalence ratios
temperatures = [1500]
phi_values = [0.5, 0.8, 1.0, 1.2, 1.5]

# Define the fuel composition (assuming a single composition for simplicity)
fuel_composition = 'NH3:0.8, H2:0.2'
oxidizer = 'O2:1, N2:3.76'

# Initialize a list to store the input files properties
props_list = []

# Generate the combinations of temperature and equivalence ratio
index = 0
for T in temperatures:
    for phi in phi_values:
        # Create the props dictionary for each combination of T and phi
        props = {
            'T': T,                        # Temperature [K]
            'P': 101325,                   # Pressure [Pa]
            'phi': phi,                    # Equivalence ratio
            'fuel': fuel_composition,      # Fuel species
            'oxidizer': oxidizer,                  # Oxidizer components [molar ratio]
            'perturb': 0,
            'simulation_time': 1e-2,
            'dt': 1e-2,
            'time_end': 1e-1,             # Maximum number of points to record
            'N_eig': 4,                    # Number of eigenvalues to store at each timestep
            'N_EI': 1,                     # Number of explosive index species
            'element_num': 5,              # Number of element numbers in a system
            'mech_file': "/data/ZhuXu/Cantera/CEMA/mechanism/NH3/NH3_otomo.cti",  # Mechanism file
            'mech_path': '/data/ZhuXu/Cantera/CEMA/mechanism/NH3',
            'last_spec': 'AR',             # Last species
            'save_dir':"data/case_NH3_0D/",
            'file_name': f'N_HR_{index + 1}.csv'  # Save path, with unique index
        }
        # Append this combination to the list
        props_list.append(props)
        index += 1
        
# if not os.path.exists(props_list[0]['save_dir']):
#     os.makedirs(props_list['save_dir'])
    


# For CH4
# props = {
#     'phi': 1.0,                   # Equivalence ratio
#     'P': 101325,                  # Pressure [Pa]
#     'T': 1500,                    # Temperature [K]
#     'fuel': 'CH4:1',               # Fuel species
#     'oxidizer': 'O2:1, N2:3.76',     # Oxidizer components [molar ratio]
#     'end_time': 1e-2,             # Maximum number of points to record
#     'N_eig': 4,                   # Number of eigenvalues to store at each timestep
#     'N_EI': 1,                    # Number of explosive index species
#     'element_num': 5,             # Number of element numbers in a system
#     'mech_path': 'src/gri30/',
#     'mech_file': "src/gri30/gri30.cti",  # Mechanism file
#     'last_spec': 'AR',             # Last species
#     'save_path': 'CH4_0Dignition'
# }


