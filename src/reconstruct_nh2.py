"""Reconstruct NH2 with QSSA using the Cantera flame profile as input."""
from __future__ import annotations

import argparse
from pathlib import Path

import cantera as ct
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from common import load_case, mechanism_path, output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/sim_props.yaml")
    parser.add_argument("--mechanism", help="Override mechanism_file from the YAML.")
    parser.add_argument("--input", help="Override outputs/simulation.csv.")
    parser.add_argument("--output-dir", help="Override the configured output directory.")
    parser.add_argument("--max-rows", type=int, help="Process only the first N rows (smoke test).")
    return parser.parse_args()


def validate_species(gas: ct.Solution, species: list[str]) -> np.ndarray:
    missing = [name for name in species if name not in gas.species_names]
    if missing:
        raise ValueError(f"Mechanism does not contain required species: {', '.join(missing)}")
    return np.array([gas.species_index(name) for name in species], dtype=int)


def reconstruct_one(
    gas: ct.Solution,
    *,
    temperature: float,
    pressure: float,
    fixed_mole_fractions: np.ndarray,
    fixed_indices: np.ndarray,
    time_end: float,
    max_step: float,
    rtol: float,
    atol: float,
) -> tuple[np.ndarray, bool, str]:
    """Integrate non-prescribed species to local quasi-steady state.

    The supplied species are kept fixed in mass-fraction form.  Other species
    start at zero, so their final values are a reconstruction rather than a
    copy of the full flame solution.
    """
    x0 = np.zeros(gas.n_species)
    x0[fixed_indices] = fixed_mole_fractions
    if np.any(x0 < -1.0e-14):
        raise ValueError("Input profile contains a negative mole fraction.")
    controlled_total = x0.sum()
    if controlled_total <= 0.0:
        raise ValueError("The prescribed species have zero total mole fraction.")

    # Removing all reconstructed species leaves a sub-normalized state.  Scale
    # the prescribed composition back to one before converting it to Y.
    x0 /= controlled_total
    gas.TPX = temperature, pressure, x0
    y0 = gas.Y.copy()

    fixed_mask = np.zeros(gas.n_species, dtype=bool)
    fixed_mask[fixed_indices] = True

    def rhs(_time: float, y: np.ndarray) -> np.ndarray:
        # Cantera's unnormalised setter permits the BDF method's intermediate
        # states while avoiding artificial renormalisation inside the ODE.
        gas.set_unnormalized_mass_fractions(y)
        gas.TP = temperature, pressure
        d_y = gas.net_production_rates * gas.molecular_weights / gas.density
        d_y[fixed_mask] = 0.0
        return d_y

    solution = solve_ivp(
        rhs,
        (0.0, time_end),
        y0,
        method="BDF",
        max_step=max_step,
        rtol=rtol,
        atol=atol,
    )
    y_final = np.clip(solution.y[:, -1], 0.0, None)
    if y_final.sum() == 0.0:
        raise RuntimeError("QSSA integration produced an empty composition.")
    gas.TPY = temperature, pressure, y_final / y_final.sum()
    return gas.X.copy(), bool(solution.success), solution.message


def reconstruct_variant(
    profile: pd.DataFrame,
    gas: ct.Solution,
    fixed_species: list[str],
    qssa: dict,
    target: str,
) -> pd.DataFrame:
    fixed_indices = validate_species(gas, fixed_species)
    target_index = validate_species(gas, [target])[0]
    source_columns = [f"X_{name}" for name in fixed_species]
    missing = [column for column in source_columns if column not in profile.columns]
    if missing:
        raise ValueError(f"Simulation CSV is missing columns: {', '.join(missing)}")

    rows: list[dict] = []
    for row_index, row in profile.iterrows():
        x, success, message = reconstruct_one(
            gas,
            temperature=float(row["T_K"]),
            pressure=float(row["P_Pa"]),
            fixed_mole_fractions=row[source_columns].to_numpy(dtype=float),
            fixed_indices=fixed_indices,
            time_end=float(qssa["time_end_s"]),
            max_step=float(qssa["max_step_s"]),
            rtol=float(qssa["rtol"]),
            atol=float(qssa["atol"]),
        )
        rows.append(
            {
                "row_index": row_index,
                "grid_m": row["grid_m"],
                "normalized_grid": row["normalized_grid"],
                "T_K": row["T_K"],
                f"X_{target}_qssa": x[target_index],
                f"Y_{target}_qssa": gas.Y[target_index],
                "solver_success": success,
                "solver_message": message,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    case, _ = load_case(args.config)
    mech = mechanism_path(case, args.mechanism)
    destination = output_dir(case, args.output_dir)
    simulation_file = Path(args.input).expanduser().resolve() if args.input else destination / "simulation.csv"
    if not simulation_file.is_file():
        raise FileNotFoundError(f"Simulation CSV not found: {simulation_file}. Run `make simulate` first.")

    profile = pd.read_csv(simulation_file)
    if args.max_rows is not None:
        profile = profile.head(args.max_rows).copy()
    required = {"grid_m", "normalized_grid", "T_K", "P_Pa"}
    if not required.issubset(profile.columns):
        raise ValueError(f"{simulation_file} is not an output from simulate_counterflow.py")

    qssa = case.get("qssa", {})
    target = qssa.get("target_species", "NH2")
    variants = qssa.get("variants", {})
    if not variants:
        raise ValueError("Set qssa.variants in the YAML.")

    for name, options in variants.items():
        gas = ct.Solution(mech)
        fixed_species = options.get("fixed_species", [])
        print(f"Reconstructing {target}: {name} ({', '.join(fixed_species)}) …")
        result = reconstruct_variant(profile, gas, fixed_species, qssa, target)
        path = destination / f"qssa_{name}.csv"
        result.to_csv(path, index=False)
        failed = int((~result["solver_success"]).sum())
        print(f"Wrote {path} ({len(result)} points; failed integrations: {failed}).")


if __name__ == "__main__":
    main()
