import os
import sys
import argparse
import cantera as ct
import numpy as np
import pandas as pd
import shutil


from utils import check_path
from props_loader import load_props_list, discover_cases

# Optional: PyCSP for CEM
import PyCSP.Functions as csp

print(f"[dbg] Python exe: {sys.executable}")
print(f"[dbg] Cantera version: {ct.__version__}")
print(f"[dbg] Cantera module: {ct.__file__}")

# ----------------------------
# CLI
# ----------------------------
# --- imports stay the same ---

def parse_args():
    cases_dir_default = "cases"
    cases = discover_cases(cases_dir_default)
    cases_help = f"Known cases in '{cases_dir_default}': {cases}" if cases else "No cases found yet."

    parser = argparse.ArgumentParser(description="Run 1D Counterflow Diffusion Flame")
    parser.add_argument("--cases_dir", type=str, default=cases_dir_default,
                        help=f"Directory containing case folders (default: {cases_dir_default}). {cases_help}")
    parser.add_argument("--case_name", type=str, required=True,
                        help="Case folder name under --cases_dir (e.g., NH3_PCI_KAUST).")
    parser.add_argument("--mech_file", type=str, required=True,
                        help="Full path to the mechanism file (cti or yaml).")
    parser.add_argument("--element_num", type=int, default=None,
                        help="Override number of elements in the mechanism (optional).")
    parser.add_argument("--inert_specie", type=str, default=None,
                        help="Override inert species (e.g., 'AR' or 'N2').")
    parser.add_argument("--loglevel", type=int, default=0,
                        help="Cantera solve loglevel (default 0).")
    parser.add_argument("--out_root", type=str, default=None,
                        help="Optional output root directory. If set, results go to <out_root>/<MECH_NAME>/")
    parser.add_argument("--props_file", type=str, default=None,
                        help="YAML filename inside the case folder OR an absolute path to the YAML.")
    

    return parser.parse_args()


class FlameExtinguished(RuntimeError):
    pass



# ----------------------------
# Core solver + post-processing
# ----------------------------
def get_nonpremixd_counterflow(gas: ct.Solution, props: dict, out_full_csv: str, loglevel: int = 0):
    """Compute a 1D non-premixed counterflow diffusion flame and write a tidy CSV with per-point fields."""
    f = ct.CounterflowDiffusionFlame(gas, width=props['width'])
    f.P = props['P']

    # inlet densities for mass flux
    gas.TPX = props['Tfuel'], props['P'], props['fuel']
    fuel_density = gas.density
    gas.TPX = props['Toxyd'], props['P'], props['oxidizer']
    oxidizer_density = gas.density

    fuel_mdot = fuel_density * props['fuel_v']
    oxidizer_mdot = oxidizer_density * props['oxyd_v']

    f.fuel_inlet.mdot = fuel_mdot
    f.fuel_inlet.X = props['fuel']
    f.fuel_inlet.T = props['Tfuel']

    f.oxidizer_inlet.mdot = oxidizer_mdot
    f.oxidizer_inlet.X = props['oxidizer']
    f.oxidizer_inlet.T = props['Toxyd']

    # refinement & initial guess
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
    f.set_initial_guess()

    # Transport / Soret from YAML
    enable_soret = bool(props.get('soret_enabled', False))
    if enable_soret:
        # Soret requires multicomponent transport
        f.transport_model = 'multicomponent'
        f.soret_enabled = True
    else:
        f.transport_model = 'mixture-averaged'
        f.soret_enabled = False

    Text = props['Text']

    def interrupt_extinction(_t):
        if np.max(f.T) < Text:
            raise FlameExtinguished('Flame extinguished')
        return 0
    f.set_interrupt(interrupt_extinction)

    f.set_max_grid_points("flame", 10000)

    # solve
    print('Creating the initial solution')
    f.solve(loglevel=loglevel, auto=True)

    # Convert flame to a SolutionArray bound to the same gas object
    # states = f.to_solution_array(gas)  # provides .T, .P, .X, .Y, etc. shape (n_points, ...)
    # npts = states.shape[0]
    grid, T_arr, X_arr, Y_arr = _extract_states(f, gas, props['P'])
    npts = T_arr.shape[0]


    # Precompute reference compositions for Bilger Z
    gas.set_equivalence_ratio(1.0, props["fuel"], props["oxidizer"])  # safe; will be overwritten below
    # Extract reference X vectors (fuel/oxidizer inlets)
    gas.TPX = props['Tfuel'], props['P'], props['fuel']
    X_F = gas.X.copy()
    gas.TPX = props['Toxyd'], props['P'], props['oxidizer']
    X_O = gas.X.copy()

    # Prepare CEM calculator once
    gas_csp = csp.CanteraCSP(props['mech_file'])
    gas_csp.constP = props['P']
    gas_csp.jacobiantype = 'full'

    # Storage
    grid = f.grid
    norm_grid = grid / (grid[-1] - grid[0])
    T = np.empty(npts)
    P = np.empty(npts)
    HRR = np.empty(npts)
    MF = np.empty(npts)
    CEM = np.empty(npts)

    # Species mass fractions (wide table later)
    species_names = gas.species_names
    Xs = X_arr.copy()  # (npts, n_species) 
    Ys = Y_arr.copy()  # (npts, n_species)
    # Ys = np.empty((npts, len(species_names)))

    # Loop over grid points and compute diagnostics
    for i in range(npts):
        Ti = T_arr[i]
        Xi = X_arr[i, :]

        # Ti = states.T[i]
        Pi = props['P']  # f.P is constant; use props['P'] which is Pa
        # Xi = states.X[i, :]

        # set both gas and gas_csp
        gas.TPX = Ti, Pi, Xi
        gas_csp.TPY = Ti, Pi, gas.Y

        # CEM
        lam, R, L, f_kernel = gas_csp.get_kernel()
        real_lam = lam.real
        # pick the (last - element_num) eigenvalue by magnitude, as in your original
        sorted_idx = np.argsort(np.abs(real_lam))
        cem_idx = sorted_idx[-1 - props["element_num"]]
        cem_real = real_lam[cem_idx]
        CEM[i] = np.sign(cem_real) * np.log10(1.0 + abs(cem_real))

        # HRR (specific, consistent with your earlier code)
        HRR[i] = -1.0 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density

        # Bilger mixture fraction at this point
        MF[i] =gas.mixture_fraction(props["fuel"], props["oxidizer"])
        # basic fields + species
        T[i] = Ti
        P[i] = Pi

    # Build tidy DataFrame
    df = pd.DataFrame({
        "grid": grid,
        "normalized_grid": norm_grid,
        "T": T,
        "P": P,
        "MF": MF,
        "HRR": HRR,
        "CEM": CEM,
    })

    # Append species mass fractions as columns Y_<name>
    for j, sp in enumerate(species_names):
        df[f"Y_{sp}"] = Ys[:, j]
        df[f"X_{sp}"] = Xs[:, j]  

    # Ensure directory and write
    os.makedirs(os.path.dirname(out_full_csv), exist_ok=True)
    df.to_csv(out_full_csv, index=False)
    print(f"Wrote {out_full_csv} ({npts} points, {len(species_names)} species)")


def _extract_states(f: ct.CounterflowDiffusionFlame, gas: ct.Solution, P: float):
    """
    Return grid, T, X, Y arrays for the flame without using to_solution_array.
    Shapes:
      grid: (n_pts,)
      T:    (n_pts,)
      X,Y:  (n_pts, n_species)
    """
    grid = f.grid.copy()
    T = f.T.copy()
    n_pts = T.size
    n_sp = gas.n_species

    # f.Y is (n_species, n_pts) on 1D flames
    Y = np.empty((n_pts, n_sp))
    X = np.empty((n_pts, n_sp))
    Y_all = np.asarray(f.Y)

    for i in range(n_pts):
        Yi = Y_all[:, i]
        gas.TPY = T[i], P, Yi
        Y[i, :] = gas.Y
        X[i, :] = gas.X

    return grid, T, X, Y

# ----------------------------
# Main
# ----------------------------
def main():
    args = parse_args()

    case_dir_arg = os.path.join(args.cases_dir, args.case_name)
    if not os.path.isdir(case_dir_arg):
        raise FileNotFoundError(f"Case directory not found: {case_dir_arg}")

    # Allow specifying a YAML by name (inside case folder) or an absolute path
    if args.props_file:
        yaml_path = args.props_file
        if not os.path.isabs(yaml_path):
            yaml_path = os.path.join(case_dir_arg, yaml_path)
        yaml_path = os.path.abspath(yaml_path)
        props_list = load_props_list(yaml_path)
    else:
        props_list = load_props_list(case_dir_arg)

    if not props_list:
        raise RuntimeError(f"No props loaded from: {args.props_file or case_dir_arg}")
    
    # Load base gas with the specified mechanism
    mech_file = os.path.abspath(args.mech_file)
    mech_folder = os.path.dirname(mech_file)
    mech_filename = os.path.basename(mech_file)
    mech_name, _ext = os.path.splitext(mech_filename)

    gas0 = ct.Solution(mech_file)
    specs = gas0.species()[:]

    # inert handling (reorder species to put inert last, as you had)
    default_inert = props_list[0].get("inert_specie")
    inert_specie = args.inert_specie or default_inert
    if inert_specie is None:
        raise ValueError("inert_specie must be set (via YAML defaults or --inert_specie).")

    try:
        inert_idx = gas0.species_index(inert_specie)
    except Exception as e:
        raise ValueError(f"Inert species '{inert_specie}' not in mechanism: {e}")

    # element_num may be overridden by CLI
    default_elem = props_list[0].get("element_num")
    override_elem = args.element_num

    # output root dir
    if args.out_root:  # explicit override
        out_case_dir = os.path.join(os.path.abspath(args.out_root), mech_name)
    else:
        output_root = props_list[0].get("output_root", "outputs")
        out_case_dir = os.path.join(output_root, mech_name)
    check_path(out_case_dir)
    if args.props_file:
        # yaml_path was resolved above when loading props_list
        props_src = yaml_path  # absolute path at this point
        try:
            shutil.copy2(props_src, output_root)
            print(f"[info] Copied props file → {output_root}")
        except Exception as e:
            print(f"[warn] Could not copy props file '{props_src}' → '{output_root}': {e}")
        
    for i, p in enumerate(props_list):
        p['mech_file']   = mech_file
        p['mech_path']   = mech_folder
        p['mech_name']   = mech_name
        p['element_num'] = int(override_elem if override_elem is not None else p.get("element_num", default_elem))
        p['inert_specie']= inert_specie

        print(f"Running simulation [{i+1}/{len(props_list)}] for case={p['case_name']} mech={mech_name}")

        # # reorder so inert is last
        # gas = ct.Solution(
        #     thermo='IdealGas', kinetics='GasKinetics', transport_model='Multi',
        #     species=specs[:inert_idx] + specs[inert_idx+1:] + [specs[inert_idx]],
        #     reactions=gas0.reactions()
        # )
        enable_soret = bool(p.get('soret_enabled', False))
        gas = ct.Solution(
            mech_file,
            transport_model=('multicomponent' if enable_soret else 'mixture-averaged')
        )
        # note: we don't rely on set_equivalence_ratio for states — we set states explicitly later
        try:
            Zst = gas.mixture_fraction(p["fuel"], p["oxidizer"])
            print("Zst =", Zst)
        except Exception:
            pass

        out_csv = os.path.join(out_case_dir, p["file_name"])
        try:
            get_nonpremixd_counterflow(gas, p, out_csv, loglevel=args.loglevel)
        except Exception as e:
            print(f"Error occurred on config {i+1}: {e}")

    # terminate worker processes if using MPI
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
        if rank != 0:
            sys.exit(0)
    except ImportError:
        pass


if __name__ == "__main__":
    main()