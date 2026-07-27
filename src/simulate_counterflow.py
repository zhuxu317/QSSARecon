"""Solve the single NH3/H2 counterflow flame with Cantera only."""
from __future__ import annotations

import argparse
from pathlib import Path

import cantera as ct
import numpy as np
import pandas as pd

from common import load_case, mechanism_path, output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/sim_props.yaml")
    parser.add_argument("--mechanism", help="Override mechanism_file from the YAML.")
    parser.add_argument("--output", help="Override the output CSV path.")
    parser.add_argument("--loglevel", type=int, default=0, help="Cantera solver verbosity.")
    return parser.parse_args()


def solve_counterflow(gas: ct.Solution, cfg: dict, loglevel: int) -> pd.DataFrame:
    flame = ct.CounterflowDiffusionFlame(gas, width=float(cfg["width"]))
    flame.P = float(cfg["P"])

    gas.TPX = float(cfg["Tfuel"]), float(cfg["P"]), cfg["fuel"]
    flame.fuel_inlet.mdot = gas.density * float(cfg["fuel_v"])
    flame.fuel_inlet.T = float(cfg["Tfuel"])
    flame.fuel_inlet.X = cfg["fuel"]

    gas.TPX = float(cfg["Toxyd"]), float(cfg["P"]), cfg["oxidizer"]
    flame.oxidizer_inlet.mdot = gas.density * float(cfg["oxyd_v"])
    flame.oxidizer_inlet.T = float(cfg["Toxyd"])
    flame.oxidizer_inlet.X = cfg["oxidizer"]

    # These are the refinement settings used by the current project case.
    flame.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05, prune=0.03)
    flame.transport_model = cfg.get("transport_model", "mixture-averaged")
    flame.set_initial_guess()
    flame.solve(loglevel=loglevel, auto=True)

    grid = flame.grid
    temperatures = flame.T
    mass_fractions = np.asarray(flame.Y).T
    mole_fractions = np.empty_like(mass_fractions)
    mixture_fraction = np.empty(grid.size)
    heat_release_rate = np.empty(grid.size)

    for row, (temperature, y) in enumerate(zip(temperatures, mass_fractions)):
        gas.TPY = temperature, float(cfg["P"]), y
        mole_fractions[row] = gas.X
        mixture_fraction[row] = gas.mixture_fraction(cfg["fuel"], cfg["oxidizer"])
        # Specific heat-release rate, matching the legacy simulation output.
        heat_release_rate[row] = -np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density

    data: dict[str, np.ndarray] = {
        "grid_m": grid,
        "normalized_grid": grid / (grid[-1] - grid[0]),
        "T_K": temperatures,
        "P_Pa": np.full(grid.size, float(cfg["P"])),
        "mixture_fraction": mixture_fraction,
        "heat_release_rate_W_kg": heat_release_rate,
    }
    for index, species in enumerate(gas.species_names):
        data[f"Y_{species}"] = mass_fractions[:, index]
        data[f"X_{species}"] = mole_fractions[:, index]
    return pd.DataFrame(data)


def main() -> None:
    args = parse_args()
    case, _ = load_case(args.config)
    cfg = case["active_config"]
    mech = mechanism_path(case, args.mechanism)
    destination = Path(args.output).expanduser().resolve() if args.output else output_dir(case) / "simulation.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)

    gas = ct.Solution(mech)
    print(f"Cantera {ct.__version__}; mechanism: {mech}")
    print(f"Solving {case.get('case_name', 'counterflow case')} …")
    result = solve_counterflow(gas, cfg, args.loglevel)
    result.to_csv(destination, index=False)
    print(f"Wrote {destination} ({len(result)} grid points, {gas.n_species} species).")


if __name__ == "__main__":
    main()
