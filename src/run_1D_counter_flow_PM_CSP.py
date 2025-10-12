#!/usr/bin/env python3
import os, sys, argparse
import numpy as np
import pandas as pd
import cantera as ct
import shutil

from utils import check_path
from props_loader import load_props_list, discover_cases

# Optional: PyCSP for CEM
import PyCSP.Functions as csp


# ----------------------------
# CLI
# ----------------------------
def parse_args():
    cases_dir_default = "cases"
    cases = discover_cases(cases_dir_default)
    cases_help = f"Known cases in '{cases_dir_default}': {cases}" if cases else "No cases found yet."

    p = argparse.ArgumentParser(description="Run twin premixed counterflow flame (both sides identical premix)")
    p.add_argument("--cases_dir", type=str, default=cases_dir_default,
                   help=f"Directory with case folders (default: {cases_dir_default}). {cases_help}")
    p.add_argument("--case_name", type=str, required=True,
                   help="Case folder name under --cases_dir (e.g., NH3_KAUST_Premixed_1bar).")
    p.add_argument("--mech_file", type=str, required=True,
                   help="Mechanism file (.yaml or .cti).")
    p.add_argument("--element_num", type=int, default=None,
                   help="Override element_num from YAML (optional).")
    p.add_argument("--inert_specie", type=str, default=None,
                   help="Override inert species, e.g. 'AR' or 'N2'.")
    p.add_argument("--loglevel", type=int, default=0,
                   help="Cantera solve loglevel (default 0).")
    p.add_argument("--out_root", type=str, default=None,
                   help="Optional output root directory. If set → <out_root>/<MECH_NAME>/")
    p.add_argument("--props_file", type=str, default=None,
               help="YAML filename (relative to case folder) or absolute path.")

    return p.parse_args()


class FlameExtinguished(RuntimeError):
    pass


# ----------------------------
# Core solver + post-processing
# ----------------------------
def get_premixed_counterflow(gas: ct.Solution, props: dict, out_full_csv: str, loglevel: int = 0, enable_soret: bool = False):
    """
    Twin premixed counterflow using CounterflowTwinPremixedFlame:
      - Both inlets (reactants & products) have the same premixed composition
        built via set_equivalence_ratio(phi, fuel, oxidizer).
      - Mass flux on both sides set from mdot = rho(T,P,premix) * mix_v.
    Outputs grid, normalized_grid, T, P, HRR, CEM, and X_/Y_ species columns.
    """
    # Required properties
    phi       = float(props["phi"])
    fuel_base = props["fuel"]       # e.g. "NH3:0.4, H2:0.45, N2:0.15"
    air_base  = props["oxidizer"]   # e.g. "O2:1, N2:3.727, AR:0.0446"
    T0        = float(props["T"])
    P0        = float(props["P"])
    width     = float(props["width"])
    mix_v     = float(props["mix_v"])  # single inlet velocity for both sides
    Text      = float(props["Text"])

    # Build the premixed composition at target phi
    premix_gas = ct.Solution(props["mech_file"])
    premix_gas.TP = T0, P0
    premix_gas.set_equivalence_ratio(phi, fuel=fuel_base, oxidizer=air_base)
    X_premix = premix_gas.X.copy()

    # Compute density for mass flux
    gas.TPX = T0, P0, X_premix
    rho = gas.density
    mdot = rho * mix_v  # kg/m^2/s

    # Instantiate twin premixed flame
    target_pts = 500  # >= 500
    init_grid = np.linspace(0.0, width, target_pts)
    f = ct.CounterflowTwinPremixedFlame(gas=gas, grid=init_grid)

    # f = ct.CounterflowTwinPremixedFlame(gas=gas, width=width)
    f.set_max_grid_points("flame", 100000)
    f.P = P0

    # # Both inlets identical
    f.reactants.mdot = mdot
    # f.reactants.T    = T0
    # f.reactants.X    = X_premix


    # Mesh / transport / refine
    # f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05, prune=0.001)
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0)

    # ▼ transport & Soret
    enable_soret = bool(props.get('soret_enabled', False))

    if enable_soret:
        f.transport_model = "multicomponent"
        f.soret_enabled   = True
    else:
        f.transport_model = "mixture-averaged"
        f.soret_enabled   = False

    # Extinction guard
    def interrupt_extinction(_t):
        if np.max(f.T) < Text:
            raise FlameExtinguished("Flame extinguished")
        return 0
    f.set_interrupt(interrupt_extinction)

    loglevel = 2
    # Solve
    f.solve(loglevel=loglevel, auto=True)

    # Extract states
    grid, T_arr, X_arr, Y_arr = _extract_states(f, gas, P0)
    npts = T_arr.shape[0]
    species_names = gas.species_names

    # CEM preparation
    gas_csp = csp.CanteraCSP(props["mech_file"])
    gas_csp.constP = P0
    gas_csp.jacobiantype = "full"

    # Allocate and compute diagnostics
    norm_grid = grid / (grid[-1] - grid[0])
    T = np.empty(npts)
    P = np.empty(npts)
    HRR = np.empty(npts)
    CEM = np.empty(npts)

    for i in range(npts):
        Ti = T_arr[i]
        Xi = X_arr[i, :]
        gas.TPX = Ti, P0, Xi
        gas_csp.TPY = Ti, P0, gas.Y

        lam, R, L, f_kernel = gas_csp.get_kernel()
        real_lam = lam.real
        idx = np.argsort(np.abs(real_lam))[-1 - props["element_num"]]
        eig = real_lam[idx]
        CEM[i] = np.sign(eig) * np.log10(1.0 + abs(eig))

        HRR[i] = -1.0 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
        T[i] = Ti
        P[i] = P0

    # Save CSV
    df = pd.DataFrame({
        "grid": grid,
        "normalized_grid": norm_grid,
        "T": T,
        "P": P,
        "HRR": HRR,
        "CEM": CEM,
    })
    for j, sp in enumerate(species_names):
        df[f"Y_{sp}"] = Y_arr[:, j]
        df[f"X_{sp}"] = X_arr[:, j]

    os.makedirs(os.path.dirname(out_full_csv), exist_ok=True)
    df.to_csv(out_full_csv, index=False)
    print(f"Wrote {out_full_csv} ({npts} points, {len(species_names)} species)")


def _extract_states(f: ct.CounterflowTwinPremixedFlame, gas: ct.Solution, P: float):
    """Return grid, T, X, Y arrays for the flame without using to_solution_array."""
    grid = f.grid.copy()
    T = f.T.copy()
    n_pts = T.size
    n_sp = gas.n_species

    Y = np.empty((n_pts, n_sp))
    X = np.empty((n_pts, n_sp))
    Y_all = np.asarray(f.Y)   # (n_species, n_pts)

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

    # Allow absolute path or path relative to the case folder
    if args.props_file:
        yaml_path = args.props_file
        if not os.path.isabs(yaml_path):
            yaml_path = os.path.join(case_dir_arg, yaml_path)
        yaml_path = os.path.abspath(yaml_path)
        props_list = load_props_list(yaml_path)
    else:
        props_list = load_props_list(case_dir_arg)

    mech_file = os.path.abspath(args.mech_file)
    mech_folder = os.path.dirname(mech_file)
    mech_filename = os.path.basename(mech_file)
    mech_name, _ext = os.path.splitext(mech_filename)

    gas0 = ct.Solution(mech_file)

    # inert override
    default_inert = props_list[0].get("inert_specie")
    inert_specie = args.inert_specie or default_inert
    if inert_specie is None:
        raise ValueError("inert_specie must be set (via YAML defaults or --inert_specie).")
    try:
        gas0.species_index(inert_specie)
    except Exception as e:
        raise ValueError(f"Inert species '{inert_specie}' not in mechanism: {e}")

    # element_num override
    default_elem = props_list[0].get("element_num")
    override_elem = args.element_num

    # output root
    if args.out_root:
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
        p["mech_file"]   = mech_file
        p["mech_path"]   = mech_folder
        p["mech_name"]   = mech_name
        p["element_num"] = int(override_elem if override_elem is not None else p.get("element_num", default_elem))
        p["inert_specie"]= inert_specie
        enable_soret = bool(p.get('soret_enabled', False))

        # Required keys for twin premixed:
        required = ("T","P","width","Text","fuel","oxidizer")
        missing = [k for k in required if k not in p]
        if missing:
            raise ValueError(f"Missing keys in props for {p.get('case_name','<unknown>')}: {missing}")

        print(f"Running simulation [{i+1}/{len(props_list)}] for case={p['case_name']} mech={mech_name}")
    
        transport = "multicomponent" if enable_soret else "mixture-averaged"
        gas = ct.Solution(mech_file, transport_model=transport)

        out_csv = os.path.join(out_case_dir, p["file_name"])
        try:
            get_premixed_counterflow(gas, p, out_csv, loglevel=args.loglevel)
        except FlameExtinguished:
            print(f"  → {p['case_name']} extinguished (Text={p['Text']} K)")
        except Exception as e:
            print(f"Error on config {i+1}: {e}")

    # optional: finalize MPI workers
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
        if rank != 0:
            sys.exit(0)
    except ImportError:
        pass


if __name__ == "__main__":
    main()