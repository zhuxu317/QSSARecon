import sys
import os
import argparse
import cantera as ct
import numpy as np
import pandas as pd
import PyCSP.Functions as csp
import PyCSP.utils as utils
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
    return parser.parse_args()


def get_cem_from_1D_flame_CSP(
    flame,
    mech_file,
    P,
    element_num=0,
    cem_store_num=5,
    fig_path='figs/eigvals_CSP.pdf',
):
    """
    Evaluate CEM along a 1-D flame using PyCSP and store/plot all the first cem_store_num eigenvalues.

    Parameters
    ----------
    flame : Cantera flame object (e.g. ct.FreeFlame)
    mech_file : str               Path to your Cantera mechanism
    P : float                     System pressure [Pa]
    element_num : int             Number of slow (conserved-element) modes to skip when picking CEM
    cem_store_num : int           How many eigenvalues per grid point to store/plot
    fig_path : str                Where to save the eigenvalue scatter plot

    Returns
    -------
    CEM_values : list[float]      Signed-log10(|CEM|) at every grid point
    lam_store   : np.ndarray      Shape (n_points, cem_store_num) real parts of eigenvalues
    """

    # rebuild CSP object
    print("debug0: Creating CSP object with mechanism:", mech_file)
    gas_csp = csp.CanteraCSP(mech_file)
    print("debug0: CSP object created with mechanism:", mech_file)
    n_pts = flame.flame.n_points
    # storage for the real-part of the first cem_store_num eigenvalues:
    lam_store = np.zeros((n_pts, cem_store_num), dtype=float)
    CEM_values = []
    print(f"debug0: n_pts = {n_pts}, cem_store_num = {cem_store_num}, element_num = {element_num}")
    for i in range(n_pts):
        try:
            # set local state
            gas_csp.TPY = (flame.T[i], P, flame.Y[:, i])
            gas_csp.constP = P
            gas_csp.jacobiantype = 'full'

            # eigen-decomposition
            lam, R, L, f = gas_csp.get_kernel()
            real_lam = lam.real

            # store the first cem_store_num modes
            lam_store[i, :] = real_lam[:cem_store_num]
            # pick the “Chemical Explosive Mode”:
            #   sort by |λ|, skip `element_num` slowest (smallest |λ|),
            #   then take the next fastest
            sorted_idx = np.argsort(np.abs(real_lam))
            cem_idx = sorted_idx[-1 - element_num]
            cem_real = real_lam[cem_idx]
            # signed log10 transform
            TLOG = np.sign(cem_real) * np.log10(1.0 + abs(cem_real))

        except Exception as err:
            # 如果某个点出错，就记录一个哑值（例如 0）
            print(f"  Warning: CEM computation failed at point {i}: {err}")
            TLOG = 0.0
        finally:
            # 确保每次都 append，保持长度 == n_pts
            CEM_values.append(TLOG)

    return CEM_values


def get_mass_fractions_Qdot_P_from_1D_flame(f, gas, fuel_component, oxid_component):
    """ Get mass fractions for each grid point in the flame and save to CSV if required """
    mass_fractions = []
    Qdot = []
    P = []
    MF = np.zeros((f.flame.n_points, 1))
    
    # Loop over grid points
    for i in range(f.flame.n_points):
        gas.TPY = f.T[i], f.P, f.Y[:, i]
        mass_fractions.append(gas.Y)  # Append the mass fractions for this grid point
        P.append(gas.P)
        # Calculate Qdot and divide by density
        qdot_value = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy)
        Qdot.append(qdot_value / gas.density)  # Divide by gas density
        MF[i] = gas.mixture_fraction(fuel_component, oxid_component)
    return mass_fractions, MF, Qdot, P


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

    # Interrupt if flame extinguishes
    Text = props['Text']
    def interrupt_extinction(t):
        if np.max(f.T) < Text:
            raise FlameExtinguished('Flame extinguished')
        return 0.
    f.set_interrupt(interrupt_extinction)

    try:
        print('Creating the initial solution')
        f.solve(loglevel=0, auto=True)

        file_name = os.path.join(props['path'],
                                 props['mech_name'],
                                 props['file_name'])
        os.makedirs(os.path.dirname(file_name), exist_ok=True)

        f.save(file_name, basis="mole", overwrite=True)
        # Post-process
        mass_fractions, MF, HRR, P = get_mass_fractions_Qdot_P_from_1D_flame(
            f, gas, props["fuel"], props["oxidizer"]
        )
        mech = props['mech_file']
        CEM_values = get_cem_from_1D_flame_CSP(f, mech,props["P"], element_num=props['element_num'])

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
    # get folder and full filename
    mech_folder = os.path.dirname(mech_file)          # "mechanism/NH3_otomo"
    mech_filename = os.path.basename(mech_file)        # "NH3_otomo.cti"

    # split off the extension
    mech_name, _ = os.path.splitext(mech_filename)     # ("NH3_otomo", ".cti")
    from props_1D import generate_props

    props_list = generate_props(case_name)
    
    props = props_list[0]
    check_path(props["path"])
    check_path(os.path.join(props["path"], "figs/"))

    # Load gas with the specified mechanism
    gas = ct.Solution(mech_file)

    for i, props in enumerate(props_list):
        props['mech_file'] = mech_file
        props['mech_path'] = mech_folder
        props['mech_name'] = mech_name
        props['element_num'] = 0  # Default value for element_num

        print(f"Running simulation for configuration {i + 1}")

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