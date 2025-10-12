#!/usr/bin/env python3
import sys
import os
import argparse

import cantera as ct
import numpy as np
import pandas as pd

from run_1D_counter_flow_NP_CSP import get_cem_from_1D_flame_CSP,get_mass_fractions_Qdot_P_from_1D_flame
from utils import check_path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run 1D Premixed Free Flame simulation and compute CEM (PP style)")
    parser.add_argument(
        '--mech_file',
        type=str,
        required=True,
        help='Path to mechanism file (.cti or .yaml)')
    parser.add_argument(
        '--width',
        type=float,
        default=0.025,
        help='Domain width [m]')
    parser.add_argument(
        '--ER',
        type=float,
        default=1.0,
        help='Equivalence ratio')
    parser.add_argument(
        '--T_min',
        type=float,
        default=300.0,
        help='Minimum inlet temperature [K]')
    parser.add_argument(
        '--T_max',
        type=float,
        default=600.0,
        help='Maximum inlet temperature [K]')
    parser.add_argument(
        '--nT',
        type=int,
        default=5,
        help='Number of temperature samples')
    parser.add_argument(
        '--element_num',
        type=int,
        default=5,
        help='Index of slow modes to skip when picking CEM')
    parser.add_argument(
        '--output_path',
        type=str,
        default='./data/case_CH4_freeflame_premixed/',
        help='Directory to save CSVs and figs')
    return parser.parse_args()

def run_premixed_free_flame(args):
    # Prepare output directories
    check_path(args.output_path)
    check_path(os.path.join(args.output_path, 'figs'))

    # Load gas
    gas = ct.Solution(args.mech_file)

    # Fuel/oxidizer definitions
    fuel = 'CH4:1.0'
    oxid = 'O2:1.0, N2:3.762'
    P = ct.one_atm

    # Set up sampling
    T_samples = np.linspace(args.T_min, args.T_max, args.nT)
    phi = args.ER

    for Tin in T_samples:
        print(f"--- Solving T={Tin:.1f} K, φ={phi:.2f} ---")
        # 1) Set state & create flame
        gas.set_equivalence_ratio(phi, fuel, oxid)
        gas.TP = Tin, P
        flame = ct.FreeFlame(gas, width=args.width)
        flame.set_refine_criteria(ratio=3, slope=0.05, curve=0.1)

        # 2) Solve
        try:
            flame.solve(loglevel=0)
        except ct.CanteraError as e:
            print(f"  ERROR solving flame at T={Tin}: {e}")
            continue

        # 3) Save raw CSV
        csv_file = os.path.join(
            args.output_path,
            f"flame_T{Tin:.1f}_phi{phi:.2f}.csv"
        )
        flame.write_csv(csv_file, quiet=True)

        # 4) Post‐process: mass fractions & mixture fraction
        mass_fractions, MF = get_mass_fractions_from_1D_flame(flame, gas)

        # 5) Compute CEM via CSP
        CEM = get_cem_from_1D_flame_CSP(
            flame, args.mech_file, P,
            element_num=args.element_num
        )

        # 6) Load & augment DataFrame
        df = pd.read_csv(csv_file)
        # drop mole‐fraction X_… columns
        df = df.loc[:, ~df.columns.str.startswith('X_')]

        # append species‐by‐species mass fractions
        for i, sp in enumerate(gas.species_names):
            df[sp] = [mf[i] for mf in mass_fractions]

        df['CEM']             = CEM
        df['MF']              = MF
        df['normalized_grid'] = flame.grid / (flame.grid[-1] - flame.grid[0])

        # 7) Save augmented CSV
        df.to_csv(csv_file, index=False)
        print(f"  Results saved to {csv_file}")

if __name__ == "__main__":
    args = parse_args()
    run_premixed_free_flame(args)

    # if using MPI, kill workers on non‐root ranks
    try:
        from mpi4py import MPI
        if MPI.COMM_WORLD.Get_rank() != 0:
            sys.exit(0)
    except ImportError:
        pass