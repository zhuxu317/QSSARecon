import sys
import os
import argparse
import cantera as ct
import numpy as np
import pandas as pd
from utils import *
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

    
def get_nonpremixd_counterflow(gas, props):
    """ Get 1D non-premixed counterflow flame """
    f = ct.CounterflowDiffusionFlame(gas, width=props['width'])
    f.P = props['P']

    # Calculate densities
    gas.TPX = props['Tfuel'], props['P'], props['fuel']
    fuel_density = gas.density
    gas.TPX = props['Toxyd'], props['P'], props['oxidizer']
    oxidizer_density = gas.density

    fuel_mdot = fuel_density * props['fuel_v']
    oxidizer_mdot = oxidizer_density * props['oxyd_v']

    # Set inlet conditions
    f.fuel_inlet.mdot = fuel_mdot
    f.fuel_inlet.X = props['fuel']
    f.fuel_inlet.T = props['Tfuel']
    f.oxidizer_inlet.mdot = oxidizer_mdot
    f.oxidizer_inlet.X = props['oxidizer']
    f.oxidizer_inlet.T = props['Toxyd']

    # Refinement and initial guess
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
    f.set_initial_guess()

    enable_soret = False  # or False

    if enable_soret:
        f.transport_model = 'Multi'   # required for Soret
        f.soret_enabled   = True
    else:
        f.transport_model = 'Mix'     # mixture-averaged (faster)
        f.soret_enabled   = False

    # print("Flame transport:", f.transport_model, "| Soret:", f.soret_enabled)
    # Interrupt if flame extinguishes
    Text = props['Text']
    def interrupt_extinction(t):
        if np.max(f.T) < Text:
            raise FlameExtinguished('Flame extinguished')
        return 0
    f.set_interrupt(interrupt_extinction)

    try:
        print('Creating the initial solution')
        f.solve(loglevel=0, auto=True)

        # Write CSV
        file_name = os.path.join(props['path'], props['mech_name'], props['file_name'])
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
        f.write_csv(file_name, quiet=False)

        # Post-process
        mass_fractions, MF, HRR, P = get_mass_fractions_Qdot_P_from_1D_flame(
            f, gas, props["fuel"], props["oxidizer"]
        )
        
        CEM_values = get_cem_from_1D_flame(f, gas, props["P"], element_num=props["element_num"])

        df = pd.read_csv(file_name)
        df['CEM'] = CEM_values
        df['MF'] = MF
        df['HRR'] = HRR
        df['P'] = P

        mf_np = np.array(mass_fractions)
        cols = [f'Y_{sp}' for sp in gas.species_names]
        mf_df = pd.DataFrame(mf_np, columns=cols)
        df = pd.concat([df, mf_df], axis=1)

        df['normalized_grid'] = f.flame.grid / (f.flame.grid[-1] - f.flame.grid[0])
        df.to_csv(file_name, index=False)

    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == "__main__":
    args = parse_args()
    case_name = args.case_name
    mech_file = args.mech_file
    element_num = args.element_num
    inert_specie = args.inert_specie
    

    # get folder and full filename
    mech_folder = os.path.dirname(mech_file)          # "mechanism/NH3_otomo"
    mech_filename = os.path.basename(mech_file)        # "NH3_otomo.cti"

    # split off the extension
    mech_name, _ = os.path.splitext(mech_filename)     # ("NH3_otomo", ".cti")
    print("mech_folder=",mech_folder)
    sys.path.append(os.path.abspath(mech_folder))
    from props_1D import generate_props
    from lib_cema import *

    props_list = generate_props(case_name)
    
    props = props_list[0]
    check_path(props["path"])
    check_path(os.path.join(props["path"], "figs/"))

    # Load gas with the specified mechanism
    gas0 = ct.Solution(mech_file)
    specs = gas0.species()[:]
    N2_ind = gas0.species_index(inert_specie)

    for i, props in enumerate(props_list):
        props['mech_file'] = mech_file
        props['mech_path'] = mech_folder
        props['mech_name'] = mech_name
        props['element_num'] = int(element_num)
        props['inert_specie'] = inert_specie
        
        print(f"Running simulation for configuration {i + 1}")
        gas = ct.Solution(
            thermo='IdealGas', kinetics='GasKinetics', transport_model='Multi',
            species=specs[:N2_ind] + specs[N2_ind+1:] + [specs[N2_ind]],
            reactions=gas0.reactions()
        )
        
        gas.set_equivalence_ratio(1.0, props["fuel"], props["oxidizer"])
        Zst = gas.mixture_fraction(props["fuel"], props["oxidizer"])
        print("Zst =", Zst)

        get_nonpremixd_counterflow(gas, props)

        # terminate worker processes if using MPI
        try:
            from mpi4py import MPI
            rank = MPI.COMM_WORLD.Get_rank()
            if rank != 0:
                sys.exit(0)
        except ImportError:
            pass