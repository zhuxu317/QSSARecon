import cantera as ct
import os
import yaml
from yaml import Loader

def update_props_from_input(input_file, props):
    with open(input_file, 'r') as f:
        pars = yaml.load(f, Loader=yaml.FullLoader)
    
    # Assuming the input YAML contains keys 'fuel' and 'oxidizer' as dictionaries
    fuel_str = ', '.join([f"{key}:{value}" for key, value in pars['fuel'].items()])
    oxidizer_str = ', '.join([f"{key}:{value}" for key, value in pars['oxidizer'].items()])
    
    # Update the props dictionary with formatted strings
    props['T'] = pars['temperature']
    props['P'] = pars['pressure']
    props['fuel'] = fuel_str
    props['oxidizer'] = oxidizer_str
    props['mech_file'] = pars['mech']
    
# file_name = 'KerM01_NH3_5e-2'
# file_name = 'KerM01_NH3_1e-1'
file_name = 'MC_NH3'

props = {
    'perturb': 0,
    'dt': 1e0,
    'time_end': 1e1,             # Maximum number of points to record
    'N_eig': 4,                    # Number of eigenvalues to store at each timestep
    'N_EI': 1,                     # Number of explosive index species
    'element_num': 5,              # Number of element numbers in a system
    'last_spec': 'AR',             # Last species
    'save_dir':"data/case_NH3_PaSR/",
}
props['target_case'] = file_name

input_file = os.path.join('src/PaSR/inputs', f'{file_name}.yaml')
print("processing", input_file)
update_props_from_input(input_file, props)

props['mech_path'] = os.path.dirname(props['mech_file'])