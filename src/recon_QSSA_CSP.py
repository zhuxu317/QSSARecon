# -*- coding: utf-8 -*-
"""
QSSA / RCCE reconstruction driven by YAML case configs, with clean separation:

  READ (SIM): <output_root>/<MECH_NAME>/<file_name>
  READ (EXP): <case_dir>/Exp/<file_name>  (or <path> if provided in YAML)
  WRITE    : <recon_output_root>/<MECH_NAME>/<file_basename>/predicted_X.csv

Which mode is used (EXP vs SIM) is controlled by YAML via defaults.Exp_data
and/or per-config override.

Created on Tue 16 Apr 2024 02:22:29 PM CST
@author: Xu Zhu
"""
import sys
import os, re
import json
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.font_manager as fm
import seaborn as snsmae
import shutil
import cantera as ct
import pandas as pd
import scipy.integrate
import argparse
from tqdm import tqdm
import PyCSP.Functions as csp
import PyCSP.utils as utils
from scipy.integrate import solve_ivp

# Load props from YAML case dir or explicit YAML file
from utils import check_path
from props_loader import load_props_list, discover_cases

# --------------------------------------------------------------------
# Warning filter (kept)
# --------------------------------------------------------------------
import warnings
_seen_nasa_warning = set()
_original_showwarning = warnings.showwarning  # Save original first!

def nasa_filter(message, category, filename, lineno, file=None, line=None):
    if "NasaPoly2::validate" in str(message):
        key = (str(message), filename, lineno)
        if key in _seen_nasa_warning:
            return  # Suppress repeat
        _seen_nasa_warning.add(key)
    _original_showwarning(message, category, filename, lineno, file, line)

warnings.showwarning = nasa_filter

# --------------------------------------------------------------------
# CLI — roots come from YAML; add --props_file so we can pick exp vs sim
# --------------------------------------------------------------------
def parse_args():
    cases_dir_default = "cases"
    cases = discover_cases(cases_dir_default)
    cases_help = f"Known cases in '{cases_dir_default}': {cases}" if cases else "No cases found yet."

    parser = argparse.ArgumentParser(
        description="Run QSSA/RCCE reconstruction from case YAMLs (roots from YAML; EXP vs SIM chosen by YAML)."
    )
    parser.add_argument("--cases_dir", type=str, default=cases_dir_default,
                        help=f"Directory containing case folders (default: {cases_dir_default}). {cases_help}")
    parser.add_argument("--case_name", type=str, required=True,
                        help="Case folder name under --cases_dir (e.g., NH3_KAUST_NP_1bar).")
    parser.add_argument("--mech_file", type=str, required=True,
                        help="Full path to the mechanism file (cti or yaml).")
    parser.add_argument("--element_num", type=int, default=None,
                        help="Override number of elements in the mechanism (optional).")
    parser.add_argument("--inert_specie", type=str, default=None,
                        help="Override inert species (e.g., 'AR' or 'N2').")
    parser.add_argument("--loglevel", type=int, default=0,
                        help="Reserved for compatibility; not used here (default 0).")
    # NEW: explicitly choose which YAML to load (e.g., exp_props.yaml vs sim_props.yaml)
    parser.add_argument("--props_file", type=str, default=None,
                        help="Optional YAML filename inside the case directory (e.g., exp_props.yaml or sim_props.yaml). "
                             "If omitted, loader will use the default in the case directory.")
    parser.add_argument(
        "--extra_species",
        nargs="+",
        default=[],
        help="Extra species to include in main_species and to suffix recon_output_root. "
            "Example: --extra_species H O OH HO2 NH NO",
    )

    return parser.parse_args()

# --------------------------------------------------------------------
# (Legacy names kept for compatibility where harmless)
# --------------------------------------------------------------------
from props_1D import AVAILABLE_CASES  # not used directly
from props_1D import generate_props   # not used

# from RCCE import RCCE
from subprocess import run
from pathlib import Path
import graphviz
import random
import glob
import subprocess

# Path to the Times New Roman font file
font_path = '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf'
font_prop = fm.FontProperties(fname=font_path)
inch = 2.54

# ====================================================================
# Core ODE & CEQ
# ====================================================================
class ReactorOde:
    def __init__(self, gas, main_species):
        self.gas = gas
        self.P = gas.P
        self.main_species = main_species
        species_names = gas.species_names
        self.mask = np.array([species in self.main_species for species in species_names])

    def __call__(self, t, y):
        self.gas.set_unnormalized_mass_fractions(y[1:])
        self.gas.TP = y[0], self.P
        rho = self.gas.density
        HRR = self.gas.net_production_rates
        dTdt = np.zeros_like(y[0])
        dYdt = HRR * self.gas.molecular_weights / rho
        dYdt[self.mask] = 0.0
        return np.hstack((dTdt, dYdt))

def run_CEQ_core(
    gas, main_species_names, initial_fractions, T, P, dt, time_end,
    moleFraction, steady_time=0.1, steady_tol=1e-12
):
    # 1) Determine inert specie from the current gas ordering (inert last!)
    inert_specie = gas.species_names[-1]
    total_controlled = sum(initial_fractions.values())
    inert_value = max(0.0, 1.0 - total_controlled)

    X = np.zeros(gas.n_species)
    for sp, mf in initial_fractions.items():
        X[gas.species_index(sp)] = mf
    if inert_value > 0.0:
        X[gas.species_index(inert_specie)] = inert_value

    gas.TP = T, P
    if moleFraction:
        gas.set_unnormalized_mole_fractions(X)
    else:
        gas.set_unnormalized_mass_fractions(X)

    y0 = np.hstack((gas.T, gas.Y))
    ode = ReactorOde(gas, main_species_names)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True)
    solver.set_initial_value(y0, 0.0)

    prev_y = y0.copy()
    while solver.successful() and solver.t < time_end:
        solver.integrate(solver.t + dt)
        y_new = solver.y
        gas.TPY = y_new[0], P, y_new[1:]
        if solver.t >= steady_time:
            if np.allclose(y_new, prev_y, atol=steady_tol, rtol=0.0):
                break
        prev_y[:] = y_new

    Y_final = gas.Y.copy()
    Y_final[gas.species_index(inert_specie)] = 0.
    gas.TPY = gas.T, P, Y_final
    return gas

def run_CEQ(gas, fig_dir, main_species_names, initial_fractions, fuel, oxyd, T, P, dt, time_end, full_data, element_num, mech_file, moleFraction, ifPremixed):
    gas = run_CEQ_core(gas, main_species_names, initial_fractions, T, P, dt, time_end, moleFraction)
    species_names = gas.species_names
    temp, pressure = gas.T, gas.P
    mole_fractions = gas.X

    data = {}
    for species, y in zip(species_names, mole_fractions):
        new_species_name = f"X_{species}" if moleFraction else f"Y_{species}"
        data[new_species_name] = y

    # override with initial provided values where applicable
    updated_data = {}
    for species, value in data.items():
        if species.startswith(("X_", "Y_")):
            base_species = species[2:]
            updated_data[species] = initial_fractions.get(base_species, value)
        else:
            updated_data[species] = value
    df = pd.DataFrame(updated_data, index=[0])

    gas_csp = csp.CanteraCSP(mech_file)
    gas_csp.TPY = temp, pressure, mole_fractions
    gas_csp.constP = pressure
    gas_csp.jacobiantype = 'full'
    lam,R,L,f = gas_csp.get_kernel()
    real_lam = lam.real
    sorted_idx = np.argsort(np.abs(real_lam))
    cem_idx = sorted_idx[-1 - element_num]
    cem_real = real_lam[cem_idx]
    TLOG = np.sign(cem_real) * np.log10(1.0 + abs(cem_real))

    if not ifPremixed:
        MF = gas.mixture_fraction(fuel, oxyd)
        df['MF'] = MF
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    

    df['CEM'] = TLOG
    df['HRR'] = HRR
    for k in ('t_res','P','x','grid','T','Normalized X','Normalized Y', "mixture_fraction", "r_mm", "phi"):
        if k in full_data:
            df[k] = full_data[k]

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    os.makedirs(fig_dir, exist_ok=True)
    file_exists = os.path.isfile(predicted_file_dir)
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header=not file_exists, index=False)

# ====================================================================
# RCCE Fortran path (kept intact; optional)
# ====================================================================
def run_RCCE_fortran(gas, fig_dir, T, P, modeci,element_num, full_data, props):
    input_file_path = "src/recon_fortran/fc_xe_fe.op"
    output_file_path = "src/recon_fortran/rcce_result.op" if modeci == "RCCE_fortran" else "src/recon_fortran/icepic_result.op"
    fortran_executable = os.path.abspath("src/recon_fortran/recon_species")

    def modify_input_file(fortran_input):
        with open(input_file_path, "r") as file:
            lines = file.readlines()
        lines[1] = ' '.join(map(str, fortran_input)) + '\n'
        with open(input_file_path, "w") as file:
            file.writelines(lines)

    all_species_names = gas.species_names
    T = gas.T
    mass_fractions = gas.Y
    enthalpies_SI = gas.h
    enthalpies_CGS = enthalpies_SI / 1e3 * 1e7
    molecular_weights = gas.molecular_weights
    specific_mole_numbers = mass_fractions / molecular_weights
    enthalpies_CGS_array = np.array([enthalpies_CGS])
    fortran_input = np.concatenate((specific_mole_numbers, enthalpies_CGS_array))

    modify_input_file(fortran_input)

    def run_fortran_code():
        subprocess.run([fortran_executable], check=True, cwd="src/recon_fortran",
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    run_fortran_code()

    def read_output_file(output_file_path, molecular_weights):
        df = pd.read_csv(output_file_path, sep=r"\s+", header=None)
        W_recon_RCCE = np.array(df.values.flatten(), dtype=float)
        molecular_weights = np.array(molecular_weights, dtype=float)
        Y_recon_RCCE = W_recon_RCCE[:-1] * molecular_weights
        T_recon_RCCE = W_recon_RCCE[-1]
        return Y_recon_RCCE, T_recon_RCCE

    Y_recon_RCCE, T_recon_RCCE = read_output_file(output_file_path, molecular_weights)
    Y_recon_mass_fractions = dict(zip(all_species_names, Y_recon_RCCE))
    gas.TPY = T_recon_RCCE, P, Y_recon_mass_fractions

    mole_fractions = gas.X
    data = {species: y for species, y in zip(all_species_names, mole_fractions)}
    df = pd.DataFrame(data, index=[0])

    gas_csp = csp.CanteraCSP(props['mech_file'])
    gas_csp.TPY = T_recon_RCCE, P, mole_fractions
    gas_csp.constP = P
    gas_csp.jacobiantype = 'full'
    lam,R,L,f = gas_csp.get_kernel()
    real_lam = lam.real
    sorted_idx = np.argsort(np.abs(real_lam))
    cem_idx = sorted_idx[-1 - element_num]
    cem_real = real_lam[cem_idx]
    TLOG = np.sign(cem_real) * np.log10(1.0 + abs(cem_real))
    df['CEM'] = TLOG

    MF = gas.mixture_fraction(props['fuel'], props['oxidizer'])
    HRR = -1 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    df['MF'] = MF
    df['HRR'] = HRR
    for k in ('t_res','x','grid','T','Normalized X','Normalized Y'):
        if k in full_data:
            df[k] = full_data[k]

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    os.makedirs(fig_dir, exist_ok=True)
    file_exists = os.path.isfile(predicted_file_dir)
    with open(predicted_file_dir, 'a') as f:
        df.to_csv(f, header = not file_exists, index=False)

# ====================================================================
# Reaction path diagram helpers
# ====================================================================
def create_reaction_path_diagram(gas, fig_dir, element, cantera_row, pred_cantera_row, mole_fraction=True):
    all_species_names = gas.species_names
    if mole_fraction:
        all_species_values = [cantera_row[f'X_{name}'] for name in all_species_names]
    else:
        all_species_values = [cantera_row[f'Y_{name}'] for name in all_species_names]
    all_fractions = dict(zip(all_species_names, all_species_values))
    species_composition = ", ".join([f"{sp}:{mf}" for sp, mf in all_fractions.items()])
    gas.TPX = cantera_row['T'], cantera_row['P'], species_composition

    diagram = ct.ReactionPathDiagram(gas, element)
    diagram.show_details = True
    diagram.title = 'Reaction path diagram following {0}'.format(element)
    diagram.label_threshold = 1e-20

    output_dir = fig_dir
    dot_file = f'{element}_rxnpath.dot'
    modified_dot_file = f'{element}_rxnpath_modified.dot'
    img_file = f'{element}_rxnpath.pdf'
    os.makedirs(output_dir, exist_ok=True)
    dot_path = os.path.join(output_dir, dot_file)
    img_path = os.path.join(output_dir, img_file)
    diagram.write_dot(dot_path)

    print(f"Wrote graphviz input file to '{dot_path}'.")
    modify_dot_file(dot_path, cantera_row, pred_cantera_row, all_species_names)
    dot_path = os.path.join(output_dir, modified_dot_file)
    run(['dot', '-Tpdf', dot_path, '-o', img_path, '-Gdpi=200'], check=True)
    print(f"Wrote graphviz output file to '{img_path}'.")
    graphviz.Source(diagram.get_dot())

def modify_dot_file(dot_file_path, cantera_row, pred_cantera_row, all_species_names):
    with open(dot_file_path, 'r') as f:
        dot_content = f.read()

    species_pattern = re.compile(r'label="([A-Za-z0-9()]+)"')
    all_species_lower = [sp.lower() for sp in all_species_names]

    def compute_error(species_name):
        key = f"X_{species_name}"
        cantera_value = cantera_row.get(key, 0.0)
        pred_value    = pred_cantera_row.get(key, 0.0)
        if cantera_value != 0.0:
            error = (pred_value - cantera_value) / cantera_value
        else:
            error = 0.0
        print(f"Species: {species_name}, key: {key}, Cantera: {cantera_value}, Pred: {pred_value}, Error: {error}")
        return error

    def error_to_color(err):
        pct = min(abs(err) * 100, 100) / 100
        r = int(255 * pct) + 100
        g = int(255 * (1 - pct)) + 100
        return f"#{min(r,255):02x}{min(g,255):02x}00"

    def modify_label_and_color(m):
        sp = m.group(1).strip()
        if sp.lower() in all_species_lower:
            err = compute_error(sp)
            pct = "set" if abs(err) < 1e-3 else f"{err*100:.2f}%"
            col = error_to_color(err)
            new_label = f'{sp} {pct}'
            return f'label="{new_label}", style=filled, fillcolor="{col}"'
        else:
            return m.group(0)

    modified = species_pattern.sub(modify_label_and_color, dot_content)
    out_path = dot_file_path.replace('.dot', '_modified.dot')
    with open(out_path, 'w') as f:
        f.write(modified)
    print(f"Modified dot file saved as {out_path}")

# ====================================================================
# Reconstruction runners
# ====================================================================
from typing import Sequence, Optional

def process_experiments_RCCE(
    gas, props: dict, mech, file_name: str, fig_dir: str,
    main_species_names: Sequence[str], fuel: str, oxidizer: str,
    modeci: str, dt: float, time_end: float,
    rand: float = 0.0, mole_fraction: bool = True,
    assist_species: Optional[Sequence[str]] = ('H2O', 'CO2'),
    ifPremixed=False
) -> None:
    df_main = pd.read_csv(file_name)

    # optional assist file
    df_assist = None
    if props.get("assist_path") and assist_species:
        assist_file = os.path.join(props["assist_path"], props["file_name"])
        if os.path.isfile(assist_file):
            print("assist path is read here!")
            df_assist = pd.read_csv(assist_file)

    for idx, row in df_main.iterrows():
        input_T = row["T"]
        input_P = props["P"]

        input_species_values = []
        for sp in main_species_names:
            col = f"{'X' if mole_fraction else 'Y'}_{sp}"
            value = row[col]
            if (df_assist is not None) and (sp in assist_species) and (col in df_assist.columns):
                value = df_assist.loc[idx, col]
            if sp == "NO":
                value = value / 1e6
            input_species_values.append(value)

        perturb = [val * rand * random.uniform(-1, 1) for val in input_species_values]
        perturbed = [max(v + dv, 0.0) for v, dv in zip(input_species_values, perturb)]
        initial_mf = dict(zip(main_species_names, perturbed))
        run_CEQ(
            gas, fig_dir, main_species_names, initial_mf,
            props['fuel'], props['oxidizer'], input_T, input_P,
            mech_file=props['mech_file'], full_data=row,
            dt=props['dt'], time_end=props['time_end'],
            element_num=props["element_num"], moleFraction=props["ifMoleFraction"], ifPremixed=ifPremixed
        )

def process_simulation_RCCE(
    gas, props, mech, file_name, fig_dir, main_species_names,
    fuel, oxidizer, modeci, dt, time_end, denoted_value_name,
    rand=0, mole_fraction=True, ifPremixed=False
):
    all_species_names = gas.species_names
    df_main = pd.read_csv(file_name, comment="#")

    # optional assist file
    df_assist = None
    assist_species = props.get('species_from_exp') or props.get('species_from_sim')
    if props.get("assist_path") and assist_species:
        print("assist path is read here!")
        assist_file = os.path.join(props["assist_path"], props["file_name"])
        if os.path.isfile(assist_file):
            df_assist = pd.read_csv(assist_file, comment="#")

    max_val, max_idx = -np.inf, 0
    for idx, row in tqdm(df_main.iterrows(), total=len(df_main), desc="Running QSSA", unit="case"):
        input_T = row["T"]
        input_P = props["P"]
        denoted_value = row.get(denoted_value_name, 0.0)
        if denoted_value > max_val:
            max_val = denoted_value
            max_idx = idx

        input_species_values = []
        for sp in main_species_names:
            col = f"{'X' if mole_fraction else 'Y'}_{sp}"
            value = row[col]
            if (df_assist is not None) and (sp in (assist_species or [])) and (col in df_assist.columns):
                value = df_assist.loc[idx, col]
            input_species_values.append(value)

        perturbation = [value * rand * random.uniform(-1, 1) for value in input_species_values]
        perturbed_species_values = [max(value + perturb, 0.0) for value, perturb in zip(input_species_values, perturbation)]
        initial_model_fractions = dict(zip(main_species_names, perturbed_species_values))

        run_CEQ(
            gas, fig_dir, main_species_names, initial_model_fractions,
            fuel, oxidizer, input_T, input_P,
            full_data=row, dt=dt, time_end=time_end,
            element_num=props['element_num'], mech_file=mech, moleFraction=props["ifMoleFraction"], ifPremixed=ifPremixed
        )

    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    if os.path.isfile(predicted_file_dir):
        df_pred = pd.read_csv(predicted_file_dir, comment='#')
        if 0 <= max_idx < len(df_pred):
            max_row = df_main.loc[max_idx]
            pred_row = df_pred.loc[max_idx]
            for element in ['O','H','N']:
                create_reaction_path_diagram(gas, fig_dir, element, max_row, pred_row, mole_fraction)

# ====================================================================
# Output directory helpers
# ====================================================================
def create_output_dir(out_mech_dir: str, file_name: str, extra_name: str | None = None):
    base_name = os.path.splitext(os.path.basename(file_name))[0]
    if extra_name:
        base_name = f"{base_name}_{extra_name}"
    fig_dir = os.path.join(out_mech_dir, base_name)
    os.makedirs(fig_dir, exist_ok=True)
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    if os.path.isfile(predicted_file_dir):
        os.remove(predicted_file_dir)
    return fig_dir

def create_noisy_fig_dir(out_mech_dir: str, file_name: str, prefix: str, noise_level: float):
    base_name, _ = os.path.splitext(os.path.basename(file_name))
    noise_level_rounded = round(noise_level, 3)
    fig_dir = os.path.join(out_mech_dir, prefix, f'{base_name}_{noise_level_rounded}')
    os.makedirs(fig_dir, exist_ok=True)
    predicted_file_dir = os.path.join(fig_dir, "predicted_X.csv")
    if os.path.isfile(predicted_file_dir):
        os.remove(predicted_file_dir)
    return fig_dir

def modify_streams_in_file(main_species_names, modeci):
    input_file = 'src/recon_fortran/streams_ori.in'
    output_file = 'src/recon_fortran/streams.in'
    with open(input_file, 'r') as file:
        lines = file.readlines()

    start_index = None
    end_index = None
    for i, line in enumerate(lines):
        if "REPRESENTED_SPECIES BEGIN" in line:
            start_index = i
        if "REPRESENTED_SPECIES END" in line:
            end_index = i

    if start_index is not None and end_index is not None:
        new_species_line = " ".join(main_species_names) + "\n"
        lines = lines[:start_index + 1] + [new_species_line] + lines[end_index:]

    for i, line in enumerate(lines):
        if line.startswith("MODECI"):
            lines[i] = f"MODECI         {modeci}\n"
            break

    with open(output_file, 'w') as file:
        file.writelines(lines)
    print(f"Modified file saved as '{output_file}'")

def _sanitize_species_tag(sp_list):
    # e.g. ["H","O","OH"] -> "withH_O_OH"
    if not sp_list:
        return ""
    return "with" + "_".join(sp_list)


# ====================================================================
# Main
# ====================================================================
def main():
    args = parse_args()

    # Resolve case dir & load YAML props
    case_dir = os.path.join(args.cases_dir, args.case_name)
    if not os.path.isdir(case_dir):
        raise FileNotFoundError(f"Case directory not found: {case_dir}")

    # Allow choosing a specific YAML file (absolute path or path inside the case folder)
    case_dir_arg = os.path.join(args.cases_dir, args.case_name)
    if not os.path.isdir(case_dir_arg):
        raise FileNotFoundError(f"Case directory not found: {case_dir_arg}")

    if getattr(args, "props_file", None):
        yaml_path = args.props_file
        if not os.path.isabs(yaml_path) and os.path.sep not in yaml_path:
            yaml_path = os.path.join(case_dir_arg, yaml_path)
        yaml_path  = os.path.abspath(yaml_path)
        props_list = load_props_list(yaml_path)
    else:
        props_list = load_props_list(case_dir_arg)

    # yaml_case_dir = props_list[0].get("__case_dir__", case_dir_arg)
    if not props_list:
        raise RuntimeError(f"No props loaded from: {getattr(args, 'props_file', None) or case_dir}")

    # Mechanism info
    mech_file = os.path.abspath(args.mech_file)
    mech_folder = os.path.dirname(mech_file)
    mech_filename = os.path.basename(mech_file)
    mech_name, _ = os.path.splitext(mech_filename)
    gas = ct.Solution(mech_file)

    # YAML flags / defaults
    default_inert = props_list[0].get("inert_specie")
    ifPremixed    = props_list[0].get("Premixed", False)
    inert_specie  = args.inert_specie or default_inert

    # ---------------- Roots (from YAML) ----------------
    sim_dir            = props_list[0].get("output_root")          # READ
    exp_dir            = props_list[0].get("exp_root")          # READ
    recon_folder_name  = props_list[0].get("recon_output_root")    # WRITE (folder name or path)

    parent_dir = os.path.dirname(sim_dir) if os.path.isabs(sim_dir) else os.path.dirname(os.path.abspath(sim_dir))
    print("recon_dir (yaml)=", recon_folder_name)

    # If recon_output_root is just a folder name (e.g., "Sim_QSSA"), put it next to output_root.
    recon_dir_root = (recon_folder_name if os.path.isabs(recon_folder_name)
                      else os.path.join(parent_dir, recon_folder_name))
    recon_dir = os.path.join(recon_dir_root, mech_name)
    print("recon_dir (resolved)=", recon_dir)
    check_path(recon_dir)
    output_root = props_list[0].get("output_root", "outputs")

    if args.props_file:
        # yaml_path was resolved above when loading props_list
        props_src = yaml_path  # absolute path at this point
        try:
            shutil.copy2(props_src, output_root)
            print(f"[info] Copied props file → {output_root}")
        except Exception as e:
            print(f"[warn] Could not copy props file '{props_src}' → '{output_root}': {e}")
        

    # --------- NEW: read main species from YAML (with a safe fallback) ----------
    main_species_default = props_list[0].get(
        "main_species",
        ["NH3", "H2", "O2", "H2O", "N2"]  # fallback if not present
    )

    # iterate props
    for i, p in enumerate(props_list):
        p['mech_file']   = mech_file
        p['mech_path']   = mech_folder
        p['mech_name']   = mech_name
        p['element_num'] = int(args.element_num if args.element_num is not None else p.get("element_num", 0))
        p['inert_specie']= inert_specie
        p.setdefault('perturb', 0.0)


        # sanitize assist/path if relative
        if p.get('assist_path') and not os.path.isabs(p['assist_path']):
            p['assist_path'] = os.path.join(case_dir, p['assist_path'])
        if p.get('path') and not os.path.isabs(p['path']):
            p['path'] = os.path.join(case_dir, p['path'])

        base_file_name = p['file_name']
        is_exp = bool(p.get('Exp_data', False))
        if is_exp:
            file_name = os.path.join(exp_dir, base_file_name)
            print("file_nam=", file_name)
        else:
            file_name = os.path.join(sim_dir, mech_name, base_file_name)

        print(f"[{i+1}/{len(props_list)}] case={p.get('case_name', args.case_name)} mech={mech_name}")
        print(f"  mode       : {'EXP' if is_exp else 'SIM'}")
        print("  input CSV  :", file_name)

        # per-config override or fall back to defaults list
        main_species_names = p.get("main_species", main_species_default)

        # method(s)
        valid_modes = {"QSSA"}  # extend to {"QSSA","RCCE_fortran"} if desired
        for modeci in valid_modes:
            extra_name = None
            fig_dir = create_output_dir(recon_dir, base_file_name, extra_name)

            if not os.path.isfile(file_name):
                print(f"  WARNING: input file not found, skip: {file_name}")
                continue

            if is_exp:
                speciesFromSim = p['species_from_sim'] if p.get('if_assist', False) else None
                process_experiments_RCCE(
                    gas, p, p['mech_file'], file_name, fig_dir,
                    main_species_names, p['fuel'], p['oxidizer'],
                    modeci, p['dt'], p['time_end'],
                    rand=p['perturb'], mole_fraction=p['ifMoleFraction'], ifPremixed=ifPremixed
                )
            else:
                process_simulation_RCCE(
                    gas, p, p['mech_file'], file_name, fig_dir,
                    main_species_names, p['fuel'], p['oxidizer'],
                    modeci, p['dt'], p['time_end'],
                    denoted_value_name="HRR",
                    rand=p['perturb'], mole_fraction=p['ifMoleFraction'], ifPremixed=ifPremixed
                )

    # Cleanly exit any MPI workers if present
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
        if rank != 0:
            sys.exit(0)
    except ImportError:
        pass

if __name__ == "__main__":
    main()