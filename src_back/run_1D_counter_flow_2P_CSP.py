#!/usr/bin/env python3
import sys, os, argparse
import cantera as ct, numpy as np, pandas as pd
from pathlib import Path

from run_1D_counter_flow_NP_CSP import (
    get_cem_from_1D_flame_CSP,
    get_mass_fractions_Qdot_P_from_1D_flame,
)
from utils import computeStrainRatesTwinPremixedFlame, computeConsumptionSpeed, check_path
from props_1D import generate_props   # <— 关键：统一读取 props

# ------------------- CLI -------------------
def parse_args(): 
    p = argparse.ArgumentParser("NH3/H2 premixed counterflow sweep")
    p.add_argument('--mech_file', required=True,
                   help="Mechanism (cti/yaml)")
    p.add_argument('--case_name', required=True,
                   help="Must match a branch in props_1D.py, e.g. NH3_premixed_sim")
    return p.parse_args()
# -------------------------------------------

# ------------------- 主循环 -----------------
def run_case(mech_file, props):
    gas = ct.Solution(
        mech_file)
    gas.set_equivalence_ratio(props['ER'], props['fuel'], props['oxidizer'])
    gas.TP = props['T'], props['P']

    mdot = gas.density * props['fuel_v']        # kg/m2/s

    f = ct.CounterflowTwinPremixedFlame(gas=gas, width=props['width'])

    # ——— Set *both* inlet mass-fluxes for the symmetric twin flame ———
    f.reactants.mdot = mdot
    f.products.mdot  = mdot

    # ——— Impose the same temperature and composition on both sides ———
    f.reactants.T = props['T']
    f.products.T  = props['T']

    # You can either reuse the gas composition array...
    # f.reactants.X = gas.X
    # f.products.X  = gas.X
    #
    # …or explicitly reconstruct the premixed string:
    f.P = props['P']
    f.set_refine_criteria(ratio=2, slope=0.1, curve=0.2, prune=0.05)

    def interrupt(t):
        if np.max(f.T) < props['Text']:
            raise ct.FlameExtinguished("extinguished")
        return 0.
    f.set_interrupt(interrupt)

    print(f"==> {props['case_name']}  φ={props['ER']}")
    f.solve(loglevel=1, auto=True)

    process_and_save(f, gas, mech_file, props)

    # strain-rate ramp
    sf = props['strain_factor']; n = 0
    try:
        f.solve(loglevel=1)
        process_and_save(f, gas, mech_file, props)
    except ct.FlameExtinguished:
        print("   → extinguished")
    except ct.CanteraError as e:
        print("   → solver error:", e)


# ------------- 结果写文件 -------------------

def process_and_save(f, gas, mech_file, props):
    # --- compute K & Sc just as before ---
    K  = computeStrainRatesTwinPremixedFlame(f, gas.density)
    Sc = computeConsumptionSpeed(f)

    # --- build output path ---
    model_name = Path(mech_file).stem  # "Otomo_32s213r"
    out_dir    = os.path.join(props['path'], model_name)
    check_path(out_dir)
    fname      = f"{props['case_name']}.csv"
    csv_path   = Path(out_dir) / fname

    # --- dump the raw flame data (this creates a header with bare species names) ---
    f.save(csv_path, basis="mole", overwrite=True)

    # --- read once (skipping your metadata comment lines) ---
    df = pd.read_csv(csv_path, comment="#")

    # --- rename the original species columns to X_<species> ---
    rename_map = {sp: f"X_{sp}" for sp in gas.species_names if sp in df.columns}
    df.rename(columns=rename_map, inplace=True)

    # --- compute the extra quantities you want to append ---
    mass_fractions, MF, HRR, P = get_mass_fractions_Qdot_P_from_1D_flame(
        f, gas, props['fuel'], props['oxidizer']
    )
    CEM = get_cem_from_1D_flame_CSP(
        f, mech_file, props['P'], element_num=props['element_num']
    )

    # --- build a DataFrame of Y_<species> columns from your mass_fractions list ---
    y_cols = [f"Y_{sp}" for sp in gas.species_names]
    df_y   = pd.DataFrame(
        np.array(mass_fractions),
        columns=y_cols,
        index=df.index
    )

    # --- concatenate everything in one go ---
    df = pd.concat([df, df_y], axis=1)
    df['CEM']             = CEM
    df['MF']              = MF
    df['HRR']             = HRR
    df['P']               = P
    df['normalized_grid'] = f.flame.grid / (f.flame.grid[-1] - f.flame.grid[0])

    # --- write back out, and re-insert your metadata comment if you like ---
    df.to_csv(csv_path, index=False)
    header = f"# K={K:.4f}, Sc={Sc:.4f}\n"
    with open(csv_path, "r+") as fh:
        content = fh.read()
        fh.seek(0)
        fh.write(header + content)

    print("   saved", csv_path)
    
    
# -------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    for pr in generate_props(args.case_name):
        pr['mech_file'] = args.mech_file
        run_case(args.mech_file, pr)
