import cantera as ct
import os
import numpy as np
# Define the function to populate props_list based on the case_name

# ─── Make a single, shared list of all supported case‐names ─────────────────
AVAILABLE_CASES = [
    "NH3_CTM",
    "NH3_PCI_KAUST",
    "NH3_PCI_Matrix",
    "NH3_NP_sim",
    "NH3_NP_exp",
    "NH3_DNS",
    "NH3_premixed_sim",
    "CH4",
]

def generate_props(case_name):
    props_list = []

    if case_name == "NH3_CTM":
        fuel_compositions = [
            'NH3:1, H2:0',
            'NH3:0.80, H2:0.20',
            'NH3:0.60, H2:0.40',
            'NH3:0.50, H2:0.50',
            'NH3:0, H2:1',
            'NH3:0.42, H2:0.42, N2:0.16',
            'NH3:0.35, H2:0.35, N2:0.3',
            'NH3:0.27, H2:0.27, N2:0.46',
            'NH3:0.19, H2:0.19, N2:0.62',
            'NH3:0.12, H2:0.12, N2:0.76'
        ]

        fuel_air = [
            0.2,
            0.2,
            0.2, 
            0.2,
            0.2,
            0.2,
            0.18,
            0.18,
            0.18,
            0.18,
        ]

        for i in range(10):  # Adjusted the range to match the number of fuel compositions
            props = {
                'file_name': f'N_CF_{i}.csv',
                'T': 290.,  
                'P': 5.0 * ct.one_atm,
                'fuel': fuel_compositions[i],
                'oxidizer': 'O2:1.0, N2:3.762',
                'Exp_data': False,
                'ifMoleFraction': True,
                'fuel_v': 0.23,
                'oxyd_v': fuel_air[i],
                'strain_factor': 1.2,
                'width': 0.05,
                'dt': 1e-2,
                'time_end': 1e-1,
                'perturb': 0,
                'mixing': None,
                'mech_file': "/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo/NH3_otomo.cti",
                'mech_path': '/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo',
                'type': "nonpremixed_counterflow",
                'path': "SIM_results/case_NH3_counterflow_KAUST/",
                'Tfuel': 290,
                'Toxyd': 290,
                'Text': 290,
                'element_num': 5,
                'inert_specie': 'AR'
            }
            props_list.append(props)

    if case_name == "NH3_PCI_KAUST":
        
        fuel_compositions = [
            'H2:0.6,  N2:0.4',
            # 'NH3:0.80, H2:0.20',
            # 'NH3:0.60, H2:0.40',
            # 'NH3:0.50, H2:0.50',
            # 'NH3:0.42, H2:0.42, N2:0.16',
            # 'NH3:0.35, H2:0.35, N2:0.3',
            # 'NH3:0.27, H2:0.27, N2:0.46',
            # 'NH3:0.19, H2:0.19, N2:0.62',
            # 'NH3:0.12, H2:0.12, N2:0.76'
        ]

        file_names = [
            "H_CF_1", 
            # "N_CF_1", 
            # "N_CF_2", 
            # "N_CF_3",
            # "N_CF_4", 
            # "N_CF_5",
            # "N_CF_6", 
            # "N_CF_7",
            # "N_CF_8",
        ]

        U_air = [
            0.23, #0.23
            # 0.2, 
            # 0.2,
            # 0.2,
            # 0.2,
            # 0.18,
            # # 0.18,
            # 0.18,
            # 0.18,
        ]


        for i, file_name in enumerate(file_names):
            props = {
                'file_name': f'{file_name}.csv',
                'Exp_data': False,
                'T': 290, #290.,  
                'P': 5e5, #5.0 * ct.one_atm,
                'fuel': fuel_compositions[i],
                'oxidizer': 'O2:1, N2:3.727, AR:0.0446',
                # 'oxidizer': 'O2:1, N2:3.727',
                'fuel_v': 0.23,
                'oxyd_v': U_air[i],
                'strain_factor': 1.2,
                'width': 0.007,#0.005,#0.005,
                'dt': 1e-2,
                'time_end': 1e-0,
                'perturb': 0,
                'mixing': None,
                'ifMoleFraction': True,
                # 'mech_file': "/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo/NH3_otomo.cti",
                # 'mech_path': '/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo',
                'type': "nonpremixed_counterflow",
                'path': "SIM_results/case_NH3_counterflow_KAUST_Validation/",
                # 'assist_path': "EXP_data/NH3_NP_CF_EXP_assist/",
                # 'species_from_exp':  ['NH3'],#['NH3','NO','O2', 'H2O','H2','N2'],
                'Tfuel': 290,
                'Toxyd': 290,
                'Text': 290,
                # 'element_num': 5,
                # 'inert_specie': 'AR'
            }
            props_list.append(props)
                
    if case_name == "NH3_PCI_Matrix":
        # 1) define a linear grid for r_H2 and r_N2
        #    r_H2 runs from 0 → 1 in 11 steps
        H2_frac_list = np.linspace(0.0, 1.0, 11)
        # H2_frac_list = [0]
        
        # N2_frac_list = np.linspace(0.0, 0.6, 9)
        N2_frac_list = [0, 0.2, 0.4, 0.6]
        # N2_frac_list = [0.2]
        print("N2_frac_list=", N2_frac_list)
        

        for r_N2 in N2_frac_list:
            rem = 1.0 - r_N2  # total fraction available for NH3 + H2
            for r_H2 in H2_frac_list:
                H2  = rem * r_H2
                NH3 = rem * (1.0 - r_H2)
                N2  = r_N2

                # Round each to two decimals
                NH3_str = f"{NH3:.2f}"
                H2_str  = f"{H2:.2f}"
                N2_str  = f"{N2:.2f}"

                # Build the composition string; omit zero‐fraction species
                parts = []
                if NH3 > 0.0:
                    parts.append(f"NH3:{NH3_str}")
                if H2 > 0.0:
                    parts.append(f"H2:{H2_str}")
                if N2 > 0.0:
                    parts.append(f"N2:{N2_str}")
                fuel_str = ", ".join(parts)

                # Build filename exactly as requested:
                file_name = f"N_CF_NH3_{NH3_str}_H2_{H2_str}_N2_{N2_str}.csv"
                print(f"Generating props for {file_name} with fuel composition: {fuel_str}")
                props = {
                    'file_name': file_name,
                    'Exp_data': False,
                    'T': 290.0,
                    'P': 5.0 * ct.one_atm,
                    'fuel': fuel_str,
                    'oxidizer': 'O2:1, N2:3.727, AR:0.0446',
                    'fuel_v': 0.23,
                    'oxyd_v': 0.20,           # adjust if you want dependence on r_N2 or r_H2
                    'strain_factor': 1.2,
                    # 'width': 0.01,
                    'width': 0.01,
                    
                    'dt': 1e-2,
                    'time_end': 1e-1,
                    'perturb': 0,
                    'mixing': None,
                    'mech_file': "/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo/NH3_otomo.cti",
                    'mech_path': '/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo',
                    'type': "nonpremixed_counterflow",
                    'path': "SIM_results/case_NH3_counterflow_PCI_matrix/",
                    'ifMoleFraction': True,
                    
                    'Tfuel': 290,
                    'Toxyd': 290,
                    'Text': 290,
                    'element_num': 5,
                    'inert_specie': 'AR'
                }
                props_list.append(props)

    # =========================== NH3 premixed counterflow ==========================
    if case_name == "NH3_premixed_sim":

        # (1) NH3:H2 摩尔配比
        blend_ratios = [
            (1, 0),
            (0.8, 0.20),
        ]

        # (2) 当量比 sweep
        phi_list = [0.85, 1.0, 1.15]
        
        pressure_list = [1, 5, 10]

        U_inlet = 0.032        # m/s, fuel & air 共用
        width   = 0.01        # m
        Text    = 290.0      # K  (灭火判据)

        for nh3, h2 in blend_ratios:
            fuel_str = f"NH3:{nh3}, H2:{h2}"
            case_stub = f"NH3_{nh3:.2f}_{h2:.2f}"

            for phi in phi_list:
                
                for pressure in pressure_list:

                    props = {
                        'case_name'   : f"{case_stub}_phi{phi:.2f}_p{pressure}",
                        'file_name'   : f"{case_stub}_phi{phi:.2f}_p{pressure}.csv",
                        'Exp_data'    : False,
                        'T'           : 290.0,
                        'P'           : pressure * ct.one_atm,
                        'ER'          : phi,
                        'fuel'        : fuel_str,
                        'oxidizer'    : "O2:1, N2:3.727, AR:0.0446",
                        'fuel_v'      : U_inlet,
                        'oxyd_v'      : U_inlet,
                        'strain_factor': 1.2,
                        'width'       : width,
                        'dt'          : 1e-2,
                        'time_end'    : 1.0,
                        'perturb'     : 0,
                        'mixing'      : None,
                        'ifMoleFraction': True,
                        'type'        : "premixed_counterflow",
                        'path'        : "SIM_results/NH3_2P_CF/",
                        'Tfuel'       : 290.0,
                        'Toxyd'       : 290.0,
                        'Text'        : Text,
                        'element_num' : 5,
                        'inert_specie': 'AR',
                    }
                    props_list.append(props)
# ==============================================================================


    if case_name == "NH3_NP_sim":
        
        case_names = [
            'AH82',
            'AH64',
            'AH46',
            'N_0',
            'N_08',
            'N_16',
            'N_24',
            'N_30',
            'AH64_K80',
            'AH64_K120',
            'AH64_K140'
        ]
        
        fuel_compositions = [
            'NH3:0.8, H2:0.2',
            'NH3:0.60, H2:0.40',
            'NH3:0.40, H2:0.60',
            'NH3:0.50, H2:0.50',
            'NH3:0.46, H2:0.46, N2:0.08',
            'NH3:0.42, H2:0.42, N2:0.16',
            'NH3:0.38, H2:0.38, N2:0.24',
            'NH3:0.35, H2:0.35, N2:0.3',
            'NH3:0.6, H2:0.4',
            'NH3:0.6, H2:0.4',
            'NH3:0.6, H2:0.4',
        ]

        
        U_fuel= [
            0.6,
            0.6,
            0.62, 
            0.60,
            0.60,
            0.60,
            0.60,
            0.60,
            0.48,
            0.72,
            0.84,
        ]

        U_air= [
            0.4,
            0.4,
            0.38, 
            0.40,
            0.40,
            0.40,
            0.40,
            0.40,
            0.32,
            0.48,
            0.56,
        ]
        
        K_ratio = [
            80/100,
            80/100,
            75/100,
            75/100,
            75/100,
            75/100,
            75/100,
            75/100,
            65/80,
            85/120,
            95/140,
            ]
        
        for i, case_name in enumerate(case_names):
            props = {
                'case_name':f'{case_name}',
                'Exp_data': False,
                'file_name': f'{case_name}.csv',
                'T': 290.,  
                'P': 1.0 * ct.one_atm,
                'fuel': fuel_compositions[i],
                'oxidizer': 'O2:1, N2:3.727, AR:0.0446',
                'fuel_v': U_fuel[i],
                'oxyd_v':   U_air[i] * K_ratio[i],
                'strain_factor': 1.2,
                'width': 0.01,
                'dt': 1e-2,
                'time_end': 1e-0,
                'perturb': 0,
                'mixing': None,
                'ifMoleFraction': True,
                'type': "nonpremixed_counterflow",
                'path': "SIM_results/NH3_NP_CF_reduce/",
                # 'assist_path': "EXP_data/NH3_NP_CF_EXP_assist/",
                # 'species_from_exp':  ['NH3'],#['NH3','NO','O2', 'H2O','H2','N2'],
                'Tfuel': 290,
                'Toxyd': 290,
                'Text': 290,
                'element_num': 5,
                'inert_specie': 'AR'
            }
            props_list.append(props)
            

    if case_name == "NH3_DNS":
        
        case_names = [
            'DNS_reduced_data',
            # "DNS_line_data",
        ]
        
        fuel_compositions = [
            'H2:1',
        ]
        U_fuel= [
            0.6,
        ]

        U_air= [
            0.4,
        ]
        
        K_ratio = [
            1,
            ]
        
        for i, case_name in enumerate(case_names):
            props = {
                'case_name':f'{case_name}',
                'Exp_data': False,
                'file_name': f'{case_name}.csv',
                'T': 290.,  
                'P': 1.0 * ct.one_atm,
                'fuel': fuel_compositions[i],
                'oxidizer': 'O2:1, N2:3.727, AR:0.0446',
                'fuel_v': U_fuel[i],
                'oxyd_v':   U_air[i] * K_ratio[i],
                'strain_factor': 1.2,
                'width': 0.01,
                'dt': 1e-1,
                'time_end': 1e-0,
                'perturb': 0,
                'mixing': None,
                'type': "premixed_counterflow",
                'ifMoleFraction': False, 
                'path': "SIM_results/DNS_dataset/NH3/plt17000",
                # 'assist_path': "EXP_data/NH3_NP_CF_EXP_assist/",
                # 'species_from_exp':  ['NH3'],#['NH3','NO','O2', 'H2O','H2','N2'],
                'Tfuel': 290,
                'Toxyd': 290,
                'Text': 290,
                'element_num': 5,
                'inert_specie': 'AR'
            }
            props_list.append(props)
            
            
    elif case_name == "NH3_NP_exp":
        
        case_names = [
            'AH82',
            'AH64',
            'AH46',
            'N_0',
            'N_08',
            'N_16',
            'N_24',
            'N_30',
            'AH64_K80',
            'AH64_K120',
            'AH64_K140'
        ]
        
        fuel_compositions = [
            'NH3:0.8, H2:0.2',
            'NH3:0.60, H2:0.40',
            'NH3:0.40, H2:0.60',
            'NH3:0.50, H2:0.50',
            'NH3:0.46, H2:0.46, N2:0.08',
            'NH3:0.42, H2:0.42, N2:0.16',
            'NH3:0.38, H2:0.38, N2:0.24',
            'NH3:0.35, H2:0.35, N2:0.3',
            'NH3:0.6, H2:0.4',
            'NH3:0.6, H2:0.4',
            'NH3:0.6, H2:0.4',
        ]

        
        U_fuel= [
            0.6,
            0.6,
            0.62, 
            0.60,
            0.60,
            0.60,
            0.60,
            0.60,
            0.48,
            0.72,
            0.84,
        ]

        U_air= [
            0.4,
            0.4,
            0.38, 
            0.40,
            0.40,
            0.40,
            0.40,
            0.40,
            0.32,
            0.48,
            0.56,
        ]
        
        K_ratio = [
            80/100,
            80/100,
            75/100,
            75/100,
            75/100,
            75/100,
            75/100,
            75/100,
            65/80,
            85/120,
            95/140,
            ]
        
        for i, case_name in enumerate(case_names):
            props = {
                'Exp_data': True,
                'file_name': f'{case_name}.csv',
                'T': 290.,  
                'P': 1.0 * ct.one_atm,
                'fuel': fuel_compositions[i],
                'oxidizer':  'O2:1, N2:3.727, AR:0.0446',
                'fuel_v': U_fuel[i],
                'oxyd_v':   U_air[i]* K_ratio[i],
                'strain_factor': 1.2,
                'width': 0.01,
                'dt': 1e-2,
                'time_end': 1e-0,
                'perturb': 0,
                'mixing': None,
                'mech_file': "/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo/NH3_otomo.cti",
                'mech_path': '/data/ZhuXu/Cantera/CEMA/mechanism/NH3_otomo',
                'path': "EXP_data/NH3_NP_CF_EXP/",
                'if_assist': True,
                'ifMoleFraction': True,
                'assist_path': "SIM_results/NH3_NP_CF_assist/",
                'species_from_sim':  ['NH3'],#['NH3','NO','O2', 'H2O','H2','N2'],
                'Tfuel': 290,
                'Toxyd': 290,
                'Text': 290,
                'element_num': 5,
                'inert_specie': 'AR'
            }
            props_list.append(props)           
    elif case_name == "CH4":
        velocity = [
            0.05,
            0.1, 
            0.15,
            0.2,
            0.25,
        ]
        
        for i in range(len(velocity)):  # Adjusted the range to match the length of the velocity list
            props = {
                'file_name': f'CH4_NP_{i}.csv',
                'T': 298.15,                # Temperature [K]
                'P': 1.0 * ct.one_atm,     # Pressure [Pa]
                'phi': 1.0,                # Equivalence ratio
                'fuel': 'CH4:1',  # Fuel components [molar ratio]
                'oxidizer': 'O2:0.21, N2:0.78, AR:0.01',   # Oxidizer components [molar ratio]
                'fuel_v': velocity[i],
                'oxyd_v': velocity[i],
                'strain_factor': 2,
                'width': 0.1 * K_ratio[i],             # Domain width for calculating flame [m]
                'mixing': None,
                # 'mech_file': "/data/ZhuXu/Cantera/NeuralRecon/mechanism/gri30/gri30.cti",  # Mechanism file
                # 'mech_path': '/data/ZhuXu/Cantera/NeuralRecon/mechanism/gri30/',
                'mech_file': "/data/ZhuXu/Cantera/NeuralRecon/mechanism/gri12/gri12.cti",  # Mechanism file
                'mech_path': '/data/ZhuXu/Cantera/NeuralRecon/mechanism/gri12/',
                'type': "nonpremixed_counterflow",  # Type of flame
                'path': "./data/case_CH4_counterflow/",  # Data path to save 1D flame data and figures
                'Tfuel': 294,              # Fuel temperature [K]
                'Toxyd': 291,              # Oxidizer temperature [K]
                'Text': 300,               # Extinction temperature [K]
                'element_num': 5,          # Number of element numbers in a system
                'inert_specie': 'AR'       # Inert species
            }
            props_list.append(props)

    return props_list

