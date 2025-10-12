import os
import yaml

EPS = 1e-12  # for tiny-float comparisons


def discover_cases(cases_root: str):
    """Return case names that have a props.yaml."""
    if not os.path.isdir(cases_root):
        return []
    return sorted(
        d for d in os.listdir(cases_root)
        if os.path.isfile(os.path.join(cases_root, d, "props.yaml"))
    )


def _to_pa(props: dict) -> float:
    """Normalize pressure fields in a config dict into Pascals."""
    if "P" in props:                 # already Pa
        return float(props["P"])
    if "P_bar" in props:             # bar → Pa
        return float(props["P_bar"]) * 1e5
    if "P_atm" in props:             # atm → Pa
        # 1 atm = 101325 Pa
        return float(props["P_atm"]) * 101325.0
    raise ValueError("Pressure not specified (need one of: P, P_bar, or P_atm).")


def _apply_defaults_and_validate(
    data: dict, raw_cfgs: list, case_dir: str
) -> list:
    """
    Apply defaults, convert pressure to Pa, attach bookkeeping fields,
    and validate required keys for all configs.
    """
    case_name = data.get("case_name") or os.path.basename(case_dir)
    defaults = data.get("defaults", {}) or {}
    output_root = data.get("output_root", "outputs")

    processed = []
    required = ["fuel", "oxidizer",  "width", "Text"]

    for cfg in raw_cfgs:
        p = {**defaults, **cfg}

        # Normalize pressure to Pa
        p["P"] = _to_pa(p)

        # Bookkeeping fields used by the runner
        p["case_name"] = case_name
        p["output_root"] = output_root

        # Ensure a filename
        if "file_name" not in p:
            pid = str(p.get("id", "config")).zfill(2)
            p["file_name"] = f"{pid}.csv"

        # Validate required keys
        missing = [k for k in required if k not in p]
        if missing:
            raise ValueError(f"Missing keys in props for {case_name}: {missing}")

        processed.append(p)

    return processed


def _expand_config_list(data: dict) -> list:
    """Handle the simple YAML shape with an explicit `configs:` list."""
    cfgs = data.get("configs", []) or []
    return _apply_defaults_and_validate(data, cfgs, data.get("__case_dir__", ""))

# ---- in props_loader.py ----
import math

def _expand_matrix_case(data: dict) -> list:
    """
    Expand a compact matrix definition into many configs.

    Supported modes:

    1) Premixed matrix (your existing schema):
       matrix:
         blend_ratios:  # list of [NH3, H2]
           - [0.8, 0.2]
         phi:   [0.85, 1.0, 1.15]
         P_atm: [1, 5, 10]
         mix_v: [0.02, 0.04]
         run_id_pattern:    "NH3_{NH3:.2f}_H2_{H2:.2f}_phi{phi:.2f}_P{P_atm}_SR_{strain_rate:.2f}"
         file_name_pattern: "NH3_{NH3:.2f}_H2_{H2:.2f}_phi{phi:.2f}_P{P_atm}_SR_{strain_rate:.2f}.csv"

    2) Non-premixed fuel-grid (new):
       matrix:
         mode: nonpremixed_fuel_grid
         H2_frac_list: [0.0, 0.1, ..., 1.0]
         N2_frac_list: [0.0, 0.2, 0.4, 0.6]
         run_id_pattern:    "N_CF_NH3_{NH3:.2f}_H2_{H2:.2f}_N2_{N2:.2f}"
         file_name_pattern: "N_CF_NH3_{NH3:.2f}_H2_{H2:.2f}_N2_{N2:.2f}.csv"
    """
    defaults = data.get("defaults", {}) or {}
    matrix   = data.get("matrix") or {}
    if not matrix:
        raise ValueError("`matrix` section is empty but matrix mode was selected.")

    mode = (matrix.get("mode") or "").strip().lower()

    # ──────────────────────────────────────────────────────────────────────────
    # Branch A: NON-PREMIXED fuel grid (H2_frac_list × N2_frac_list)
    # ──────────────────────────────────────────────────────────────────────────
    if mode == "nonpremixed_fuel_grid" or (
        "H2_frac_list" in matrix or "N2_frac_list" in matrix
    ):
        H2_list = matrix.get("H2_frac_list")
        N2_list = matrix.get("N2_frac_list")
        if not H2_list or not N2_list:
            raise ValueError("`matrix` requires both H2_frac_list and N2_frac_list for nonpremixed_fuel_grid.")

        run_pat  = matrix.get("run_id_pattern",    "N_CF_NH3_{NH3:.2f}_H2_{H2:.2f}_N2_{N2:.2f}")
        file_pat = matrix.get("file_name_pattern", "N_CF_NH3_{NH3:.2f}_H2_{H2:.2f}_N2_{N2:.2f}.csv")

        cfgs = []
        for r_N2 in N2_list:
            r_N2 = float(r_N2)
            rem  = 1.0 - r_N2  # available to NH3 + H2
            for r_H2 in H2_list:
                r_H2 = float(r_H2)
                H2   = rem * r_H2
                NH3  = rem * (1.0 - r_H2)
                N2   = r_N2

                # Build fuel string (omit zero terms cleanly)
                parts = []
                if NH3 > 0.0: parts.append(f"NH3:{NH3:.2f}")
                if H2  > 0.0: parts.append(f"H2:{H2:.2f}")
                if N2  > 0.0: parts.append(f"N2:{N2:.2f}")
                fuel_str = ", ".join(parts) if parts else "NH3:0.00"

                run_id    = run_pat.format(NH3=NH3, H2=H2, N2=N2)
                file_name = file_pat.format(NH3=NH3, H2=H2, N2=N2)

                cfgs.append({
                    **defaults,
                    "run_id": run_id,
                    "file_name": file_name,
                    "fuel": fuel_str,
                    # Keep pressure in atm; your validator will compute P (Pa) if needed
                    "P_atm": float(defaults.get("P_atm", 1.0)),
                })

        return _apply_defaults_and_validate(data, cfgs, data.get("__case_dir__", ""))

    # ──────────────────────────────────────────────────────────────────────────
    # Branch B: PREMIXED matrix (existing)
    # ──────────────────────────────────────────────────────────────────────────
    blends   = matrix.get("blend_ratios") or []
    phi_list = matrix.get("phi") or []
    p_list   = matrix.get("P_atm") or []
    mv_list  = matrix.get("mix_v") or []

    if not blends:
        raise ValueError("matrix.blend_ratios must be a non-empty list of [NH3, H2].")
    if not phi_list:
        raise ValueError("matrix.phi must be a non-empty list.")
    if not p_list:
        raise ValueError("matrix.P_atm must be a non-empty list.")
    if not mv_list:
        raise ValueError("matrix.mix_v must be a non-empty list.")

    run_pat  = matrix.get("run_id_pattern",
                          "NH3_{NH3:.2f}_H2_{H2:.2f}_phi{phi:.2f}_P{P_atm}_SR_{strain_rate:.2f}")
    file_pat = matrix.get("file_name_pattern",
                          "NH3_{NH3:.2f}_H2_{H2:.2f}_phi{phi:.2f}_P{P_atm}_SR_{strain_rate:.2f}.csv")

    width = float(defaults.get("width", 0.01))
    if width <= 0:
        raise ValueError("defaults.width must be positive to compute strain_rate.")
    cfgs = []

    for br in blends:
        try:
            NH3, H2 = float(br[0]), float(br[1])
        except Exception:
            raise ValueError(f"Invalid blend_ratios entry: {br!r} (expect [NH3, H2]).")

        fuel_str = f"NH3:{NH3}, H2:{H2}"

        for phi in phi_list:
            phi_val = float(phi)
            for p_atm in p_list:
                p_atm_val = float(p_atm)
                for mv in mv_list:
                    mix_v = float(mv)
                    strain_rate = 2.0 * mix_v / width

                    run_id    = run_pat.format(NH3=NH3, H2=H2, phi=phi_val, P_atm=int(p_atm_val),
                                               strain_rate=strain_rate, mix_v=mix_v)
                    file_name = file_pat.format(NH3=NH3, H2=H2, phi=phi_val, P_atm=int(p_atm_val),
                                                strain_rate=strain_rate, mix_v=mix_v)

                    cfgs.append({
                        **defaults,
                        "run_id": run_id,
                        "file_name": file_name,
                        "fuel": fuel_str,
                        "phi": phi_val,
                        "P_atm": p_atm_val,
                        "mix_v": mix_v,
                        "strain_rate": strain_rate,
                    })

    return _apply_defaults_and_validate(data, cfgs, data.get("__case_dir__", ""))

def load_props_list(props_source: str) -> list:
    """
    Load props from either a case directory containing props.yaml,
    or a direct YAML file path (e.g., cases/<case>/exp_props.yaml).

    Returns a list[dict] of expanded configs, with top-level fields
    (output_root, recon_output_root, inert_specie, element_num, etc.)
    propagated into each config, and a __case_dir__ breadcrumb.
    """
    # Resolve the YAML path
    if os.path.isdir(props_source):
        props_path = os.path.join(props_source, "props.yaml")
        case_dir = os.path.abspath(props_source)
    else:
        props_path = os.path.abspath(props_source)
        case_dir = os.path.dirname(props_path)

    if not os.path.isfile(props_path):
        raise FileNotFoundError(f"Props YAML not found at: {props_path}")

    with open(props_path, "r") as f:
        data = yaml.safe_load(f) or {}

    # Breadcrumb for error messages / relative path resolution
    data["__case_dir__"] = case_dir

    # Pull top-level roots (may be absent)
    top_output_root        = data.get("output_root")
    top_recon_output_root  = data.get("recon_output_root")

    # Expand your schema into a list of configs
    if "matrix" in data and data.get("matrix"):
        configs = _expand_matrix_case(data)
    elif "configs" in data and data.get("configs"):
        configs = _expand_config_list(data)
    else:
        raise ValueError(
            "Props YAML must contain either a non-empty `configs:` list or a non-empty `matrix:` section.\n"
            f"  File: {props_path}"
        )

    # Propagate top-level fields into each config (if not already present)
    for p in configs:
        p.setdefault("__case_dir__", case_dir)
        if top_output_root is not None:
            p.setdefault("output_root", top_output_root)
        if top_recon_output_root is not None:
            p.setdefault("recon_output_root", top_recon_output_root)

    return configs