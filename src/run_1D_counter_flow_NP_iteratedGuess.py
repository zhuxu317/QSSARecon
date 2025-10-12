import sys
import os
import argparse
import cantera as ct
import numpy as np
import pandas as pd
from utils import *
from props_1D import AVAILABLE_CASES

# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(description="Run 1D Counterflow Flame simulation")
    parser.add_argument(
        '--case_name',
        type=str,
        choices=AVAILABLE_CASES,
        required=True,
        help=("Specify the case name. " f"Valid options: {AVAILABLE_CASES}")
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

# ──────────────────────────────────────────────────────────────────────────────
# Flame solve
# ──────────────────────────────────────────────────────────────────────────────
def get_nonpremixd_counterflow(gas, props, prev_profile=None):
    """
    Solve a non-premixed counterflow flame.
    If prev_profile is a non-empty SolutionArray/DataFrame from the *same mechanism*,
    use it as the initial guess; otherwise fall back to default initial guess.
    """
    f = ct.CounterflowDiffusionFlame(gas, width=props['width'])
    f.P = props['P']

    # Inlet mass fluxes from velocities
    gas.TPX = props['Tfuel'], props['P'], props['fuel']
    fuel_density = gas.density
    gas.TPX = props['Toxyd'], props['P'], props['oxidizer']
    oxidizer_density = gas.density

    f.fuel_inlet.mdot     = fuel_density     * props['fuel_v']
    f.fuel_inlet.X        = props['fuel']
    f.fuel_inlet.T        = props['Tfuel']
    f.oxidizer_inlet.mdot = oxidizer_density * props['oxyd_v']
    f.oxidizer_inlet.X    = props['oxidizer']
    f.oxidizer_inlet.T    = props['Toxyd']
    
    f.transport_model = 'multicomponent'
    f.soret_enabled = False
    # Use the gas' configured transport model (don't override here)

    # Refinement (keep as before)
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

    # ---------- Initialization ----------
    if prev_profile is None:
        print("[init] Using default set_initial_guess()")
        f.set_initial_guess()
    else:
        print("[init] Using previous SolutionArray/DataFrame as initial guess")
        init_ok = False
        # 1) Try SolutionArray / DataFrame directly
        try:
            f.set_initial_guess(data=prev_profile)
            init_ok = True
        except Exception as e:
            print(f"[init] Direct data init failed: {e}")
        # 2) If previous is a SolutionArray, try converting to DataFrame
        if not init_ok and hasattr(prev_profile, "to_pandas"):
            try:
                f.set_initial_guess(data=prev_profile.to_pandas())
                init_ok = True
            except Exception as e:
                print(f"[init] DataFrame init failed: {e}")
        # 3) Fallback
        if not init_ok:
            print("[init] Falling back to default initial guess")
            f.set_initial_guess()

    # Extinction interrupt
    Text = props['Text']
    def interrupt_extinction(t):
        if np.max(f.T) < Text:
            raise FlameExtinguished('Flame extinguished')
        return 0
    f.set_interrupt(interrupt_extinction)

    # Solve
    print('Solving...')
    try:
        f.solve(loglevel=0, auto=True)
    except Exception as e:
        print(f"Error occurred during solve: {e}")

    # Save CSV
    file_name = os.path.join(props['path'], props['mech_name'], props['file_name'])
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    try:
        f.write_csv(file_name, quiet=False)
    except Exception as e:
        print(f"Warning: write_csv failed: {e}")

    # Post-process
    try:
        mass_fractions, MF, HRR, P = get_mass_fractions_Qdot_P_from_1D_flame(
            f, gas, props["fuel"], props["oxidizer"]
        )
        CEM_values = get_cem_from_1D_flame(f, gas, props["P"], element_num=props["element_num"])

        df = pd.read_csv(file_name)
        df['CEM'] = CEM_values
        df['MF']  = MF
        df['HRR'] = HRR
        df['P']   = P

        mf_np = np.array(mass_fractions)
        cols = [f'Y_{sp}' for sp in gas.species_names]
        mf_df = pd.DataFrame(mf_np, columns=cols)
        df = pd.concat([df, mf_df], axis=1)

        df['normalized_grid'] = f.flame.grid / (f.flame.grid[-1] - f.flame.grid[0])
        df.to_csv(file_name, index=False)
    except Exception as e:
        print(f"Warning: post-processing failed: {e}")

    # Export a profile suitable for re-use (compatible with Cantera 2.6.0a3)
    next_profile = None

    # 1) Preferred: SolutionArray (no args for this version)
    try:
        sa = f.to_solution_array()
        # Ensure non-empty
        try:
            if sa is not None and len(sa) > 0:
                next_profile = sa
        except Exception:
            pass
    except Exception as e:
        print(f"[init-export] to_solution_array failed: {e}")

    # 2) Fallback: pandas DataFrame
    if next_profile is None:
        try:
            df_prof = f.to_pandas()  # 2.6.0a3 has to_pandas()
            if df_prof is not None and len(df_prof) > 0:
                next_profile = df_prof
        except Exception as e:
            print(f"[init-export] pandas export failed: {e}")

    if next_profile is None:
        print("[init-export] Not updating previous profile (empty export).")

    return f, next_profile

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    args = parse_args()
    case_name     = args.case_name
    mech_file     = args.mech_file
    element_num   = args.element_num
    inert_specie  = args.inert_specie

    mech_folder  = os.path.dirname(mech_file)
    mech_filename = os.path.basename(mech_file)
    mech_name, _  = os.path.splitext(mech_filename)
    print("mech_folder=", mech_folder)
    sys.path.append(os.path.abspath(mech_folder))
    from props_1D import generate_props
    from lib_cema import *

    props_list = generate_props(case_name)

    # Ensure output dirs
    tmp_props = props_list[0]
    check_path(tmp_props["path"])
    check_path(os.path.join(tmp_props["path"], "figs/"))

    # Mechanism base gas
    gas0   = ct.Solution(mech_file)
    specs  = gas0.species()[:]
    N2_ind = gas0.species_index(inert_specie)

    # Keep previous profile **per mechanism**
    prev_profile_by_mech = {}

    total = len(props_list)
    for idx, props in enumerate(props_list):
        # Make the gas for THIS entry:
        gas = ct.Solution(
            thermo='IdealGas', kinetics='GasKinetics',
            transport_model='Multi',  # or 'Mix' if you prefer
            species=specs[:N2_ind] + specs[N2_ind+1:] + [specs[N2_ind]],
            reactions=gas0.reactions()
        )

        # Props metadata
        props['mech_file']   = mech_file
        props['mech_path']   = mech_folder
        props['mech_name']   = mech_name
        props['element_num'] = int(element_num)
        props['inert_specie']= inert_specie

        # Progress with mech/case context
        print(f"\n=== Running ({idx+1}/{total}) mech={props.get('mech_name','?')} "
              f"file={props.get('file_name','?')} ===")

        # Report Zst
        gas.set_equivalence_ratio(1.0, props["fuel"], props["oxidizer"])
        Zst = gas.mixture_fraction(props["fuel"], props["oxidizer"])
        print("Zst =", Zst)

        # Choose previous profile only if from same mechanism
        key = props['mech_name']
        prev_profile = prev_profile_by_mech.get(key, None)

        # Solve
        flame, next_profile = get_nonpremixd_counterflow(gas, props, prev_profile=prev_profile)

        # Store the profile for THIS mechanism only
        if next_profile is not None and (
            (hasattr(next_profile, '__len__') and len(next_profile) > 0) or
            (hasattr(next_profile, 'shape') and next_profile.shape[0] > 0)
        ):
            prev_profile_by_mech[key] = next_profile
        else:
            print("[init-export] Not updating previous profile (empty export).")

        # terminate worker processes if using MPI
        try:
            from mpi4py import MPI
            rank = MPI.COMM_WORLD.Get_rank()
            if rank != 0:
                sys.exit(0)
        except ImportError:
            pass
