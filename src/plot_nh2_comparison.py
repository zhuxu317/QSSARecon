"""Plot and tabulate the baseline vs. baseline + OH NH2 reconstruction."""
from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import load_case, output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/sim_props.yaml")
    parser.add_argument("--output-dir", help="Override the configured output directory.")
    return parser.parse_args()


def metrics(reference: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    error = prediction - reference
    return {
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "max_abs_error": float(np.max(np.abs(error))),
        "peak_X_NH2": float(np.max(prediction)),
    }


def main() -> None:
    args = parse_args()
    case, _ = load_case(args.config)
    directory = output_dir(case, args.output_dir)
    qssa = case["qssa"]
    target = qssa.get("target_species", "NH2")
    variants = qssa["variants"]

    simulation = pd.read_csv(directory / "simulation.csv")
    reference_column = f"X_{target}"
    if reference_column not in simulation:
        raise ValueError(f"simulation.csv does not contain {reference_column}")
    merged = simulation[["grid_m", "normalized_grid", "T_K", reference_column]].copy()

    summary: list[dict] = []
    for name in variants:
        result = pd.read_csv(directory / f"qssa_{name}.csv")
        prediction_column = f"X_{target}_qssa"
        result = result[["row_index", prediction_column]].rename(columns={prediction_column: name})
        merged = merged.merge(result, left_index=True, right_on="row_index", how="inner").drop(columns="row_index")
        summary.append({"variant": name, **metrics(merged[reference_column].to_numpy(), merged[name].to_numpy())})

    merged.to_csv(directory / f"{target.lower()}_comparison.csv", index=False)
    pd.DataFrame(summary).to_csv(directory / f"{target.lower()}_error_summary.csv", index=False)

    fig, axis = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    axis.plot(merged["grid_m"] * 1.0e3, merged[reference_column], color="black", lw=2.2, label="Cantera flame")
    for name in variants:
        axis.plot(merged["grid_m"] * 1.0e3, merged[name], lw=1.5, label=f"QSSA: {name.replace('_', ' ')}")
    axis.set(xlabel="Distance (mm)", ylabel=f"Mole fraction, X({target})")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False)
    fig.savefig(directory / f"{target.lower()}_comparison.png", dpi=220)
    print(f"Wrote {directory / f'{target.lower()}_comparison.csv'}")
    print(f"Wrote {directory / f'{target.lower()}_comparison.png'}")


if __name__ == "__main__":
    main()
