"""Small helpers shared by the three standalone scripts."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str | Path) -> Path:
    """Resolve config paths relative to this open-source project."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def load_case(config_path: str | Path) -> tuple[dict[str, Any], Path]:
    """Load and minimally validate the one supported counterflow case."""
    path = Path(config_path).expanduser().resolve()
    with path.open(encoding="utf-8") as stream:
        case = yaml.safe_load(stream)

    if not isinstance(case, dict):
        raise ValueError(f"Expected a YAML mapping in {path}")
    if len(case.get("configs", [])) != 1:
        raise ValueError("This minimal example supports exactly one config entry.")

    defaults = case.get("defaults", {})
    cfg = {**defaults, **case["configs"][0]}
    required = ("fuel", "oxidizer", "Tfuel", "Toxyd", "width", "Text", "fuel_v", "oxyd_v")
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"Missing required case fields: {', '.join(missing)}")

    if "P" in cfg:
        cfg["P"] = float(cfg["P"])
    elif "P_bar" in cfg:
        cfg["P"] = float(cfg["P_bar"]) * 1.0e5
    elif "P_atm" in cfg:
        cfg["P"] = float(cfg["P_atm"]) * 101325.0
    else:
        raise ValueError("Set one of P, P_bar, or P_atm in defaults.")

    case["active_config"] = cfg
    return case, path


def mechanism_path(case: dict[str, Any], override: str | None) -> Path:
    value = override or case.get("mechanism_file")
    if not value:
        raise ValueError("Set mechanism_file in the YAML or pass --mechanism.")
    path = resolve_path(value)
    if not path.is_file():
        raise FileNotFoundError(
            f"Mechanism not found: {path}\n"
            "See mechanism/README.md for how to supply an authorised mechanism."
        )
    return path


def output_dir(case: dict[str, Any], override: str | None = None) -> Path:
    path = resolve_path(override or case.get("output_dir", "outputs"))
    path.mkdir(parents=True, exist_ok=True)
    return path
