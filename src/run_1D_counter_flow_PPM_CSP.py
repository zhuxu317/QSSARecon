#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Partially-premixed counter-flow flame solved like the simple reference,
but with differential diffusion enabled in the final convergence step.

Workflow:
  1) Build CounterflowDiffusionFlame with Mix transport (stable warm-start)
  2) Warm-start with O2-free (or user-provided) fuel inlet
  3) Linearly ramp fuel-side composition to the target premix
     - H2 and NH3 kept at target values
     - other species (e.g., O2, N2, …) scaled by alpha ∈ (0,1]
     - inlet mdot values kept constant during the ramp
  4) Switch to Multi transport (differential diffusion) and enable Soret
     and re-converge
  5) Output full CSV (grid, T, X/Y, HRR, CEM)
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import cantera as ct
import os
import shutil
from utils import check_path
from props_loader import load_props_list, discover_cases
import PyCSP.Functions as csp  # for CEM

# ------------------------ CLI ------------------------
def parse_args():
    base = "cases"
    known = discover_cases(base)
    p = argparse.ArgumentParser(description="Partially-premixed counterflow (simple ramp → Multi transport)")
    p.add_argument("--cases_dir", default=base,
                   help=f"Root of case folders (default '{base}'). Found: {known}")
    p.add_argument("--case_name", required=True,
                   help="Case folder under --cases_dir")
    p.add_argument("--mech_file", required=True,
                   help="Mechanism file (.yaml/.cti)")
    p.add_argument("--element_num", type=int, default=None,
                   help="Override element_num (CEM eig index)")
    p.add_argument("--loglevel", type=int, default=1,
                   help="Cantera solver loglevel (0–5)")
    p.add_argument("--out_root", default=None,
                   help="Optional output root, else uses YAML 'output_root'")
    # optional: enable/disable Soret on final Multi solve
    p.add_argument("--soret", action="store_true",
                   help="Enable Soret effect (thermal diffusion) in final Multi solve")
    p.add_argument("--props_file", type=str, default=None,
                    help="YAML filename (relative to case folder) or absolute path.")
    return p.parse_args()


# --------------------- helpers -----------------------
def _parse_X(s: str) -> dict[str, float]:
    """'A:0.1, B:0.2' → {'A':0.1, 'B':0.2}"""
    out = {}
    for tok in s.split(','):
        tok = tok.strip()
        if not tok:
            continue
        k, v = tok.split(':')
        out[k.strip()] = float(v)
    return out


def _dict_to_X(d: dict[str, float]) -> str:
    """dict → 'A:.., B:..' (Cantera normalizes automatically)."""
    return ", ".join(f"{k}:{d[k]:.6g}" for k in sorted(d))


def robust_solve(flame: ct.CounterflowDiffusionFlame, loglevel: int = 1):
    """Try auto; else energy OFF then ON, no auto (stable)."""
    try:
        flame.solve(loglevel, refine_grid=True, auto=True)
        return
    except ct.CanteraError as err:
        print(f"[robust_solve] Solve failed with auto=True: {err}")
        pass
    # flame.energy_enabled = False
    # flame.solve(loglevel, refine_grid=True, auto=False)
    # flame.energy_enabled = True
    # flame.solve(loglevel, refine_grid=True, auto=False)


def _extract_states(f: ct.CounterflowDiffusionFlame, gas: ct.Solution, P: float):
    """grid, T, X, Y arrays (no to_solution_array usage)."""
    grid = f.grid.copy()
    T = f.T.copy()
    n_pts, n_sp = T.size, gas.n_species
    Y = np.empty((n_pts, n_sp))
    X = np.empty((n_pts, n_sp))
    Y_all = np.asarray(f.Y)  # (n_sp, n_pts)
    for i in range(n_pts):
        gas.TPY = T[i], P, Y_all[:, i]
        Y[i, :] = gas.Y
        X[i, :] = gas.X
    return grid, T, X, Y


# ------------------- core solver ---------------------
def run_case_simple(gas: ct.Solution, props: dict, out_csv: Path,
                    loglevel: int = 1, n_steps: int | None = None, enable_soret: bool = False):
    """
    Parametrized counterflow solver:
      - Right inlet uses props['oxidizer']
      - Target fuel premix is props['fuel'] (dict)
      - Optional props['fuel_init'] for warm start (no-O2 typical)
      - Keep mdot constant during the ramp
      - COUPLED ramp: fuel-side temperature and composition progress together
        with adaptive backtracking to avoid dropping to the cold/non-reacting branch
      - Transport is MULTICOMPONENT from the start and SORET is ALWAYS ON
      - No final polish step

    Optional YAML keys for the adaptive ramp (with defaults shown):
      Tfuel_warm_start:  props.get("Tfuel_warm_start", T_f)   # e.g., 2000
      n_Tfuel_steps:     props.get("n_Tfuel_steps", 0)         # used to size the ramp
      T_cold_threshold:  props.get("T_cold_threshold", 800.0)  # below this → treat as cold
      min_step:          props.get("min_step", 0.02)           # minimum step in s (0..1)
      alpha_start:       props.get("alpha_start", 0.2)
      alpha_end:         props.get("alpha_end", 1.0)
      tighten_frac:      props.get("tighten_frac", 0.7)        # tighten grid after this s
    """
    # --- basic inputs from YAML
    P0     = float(props["P"])
    width  = float(props["width"])
    T_f    = float(props.get("Tfuel", props["T"]))
    T_o    = float(props.get("Toxyd", props["T"]))
    U_f    = float(props["Ufuel"])
    U_o    = float(props["Uair"])

    # Optional temperature warm-start (fuel side)
    Tfuel_warm_start = float(props.get("Tfuel_warm_start", T_f))   # e.g., 2000 K
    n_Tfuel_steps    = int(props.get("n_Tfuel_steps", 0))          # only used to size the coupled loop

    # Step count for the coupled ramp
    n_steps = int(n_steps if n_steps is not None else props.get("n_steps", 10))
    if n_steps < 1:
        n_steps = 1

    # Adaptive ramp controls (can be set in YAML)
    T_cold_threshold = float(props.get("T_cold_threshold", 800.0))
    min_ds           = float(props.get("min_step", 0.02))
    alpha0           = float(props.get("alpha_start", 0.2))
    alpha1           = float(props.get("alpha_end",   1.0))
    tighten_frac     = float(props.get("tighten_frac", 0.7))

    # Compositions
    comp_o_air = props["oxidizer"]          # e.g., 'O2:0.21, N2:0.79'
    target = _parse_X(props["fuel"])        # target fuel premix (dict)

    # Warm-start fuel composition (no O2 typical)
    fuel_init = props.get("fuel_init", "").strip()
    if fuel_init:
        comp_f_init = fuel_init
    else:
        parts = []
        if "H2" in target and target["H2"] > 0:
            parts.append(f"H2:{target['H2']}")
        if "NH3" in target and target["NH3"] > 0:
            parts.append(f"NH3:{target['NH3']}")
        if not parts:
            parts = ["H2:1"]  # fallback
        comp_f_init = ", ".join(parts)

    # ---------- densities → mdot (kept constant during ramp)
    gas.TPX = T_o, P0, comp_o_air
    rho_o = gas.density
    gas.TPX = T_f, P0, comp_f_init
    rho_f = gas.density

    mdot_o = rho_o * U_o
    mdot_f = rho_f * U_f

    print(f"mdot_o = {mdot_o:.6g} kg/m^2/s")
    print(f"mdot_f = {mdot_f:.6g} kg/m^2/s")

    # ---------- flame object: START DIRECTLY with Multi + Soret ----------

    target_pts = 1200  # >= 500
    init_grid = np.linspace(0.0, width, target_pts)
    # f = ct.CounterflowDiffusionFlame(gas, width=width)
    f = ct.CounterflowDiffusionFlame(gas, grid=init_grid)
    f.P = P0
    enable_soret = bool(props.get('soret_enabled', False))

    if enable_soret:
        # Soret requires multicomponent transport
        f.transport_model = 'multicomponent'
        f.soret_enabled = True
    else:
        f.transport_model = 'mixture-averaged'
        f.soret_enabled = False

    # f.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.062)
    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05, prune=0.001)


    # ---------- inlets
    f.oxidizer_inlet.mdot = mdot_o
    f.oxidizer_inlet.T    = T_o
    f.oxidizer_inlet.X    = comp_o_air

    # Start the fuel at a warm temperature (if specified) and warm-start composition
    T_start = max(T_f, Tfuel_warm_start)
    f.fuel_inlet.mdot = mdot_f
    f.fuel_inlet.T    = T_start
    f.fuel_inlet.X    = comp_f_init

    # Initial robust solve at warm start
    robust_solve(f, loglevel)
    f.set_max_grid_points("flame", 100000)

    # ---------- COUPLED ramp with adaptive backtracking ----------
    n_coupled = max(n_steps, n_Tfuel_steps, 1)
    print(f"[ramp] coupling fuel T {T_start:.0f}→{T_f:.0f} K and alpha {alpha0}→{alpha1} in ~{n_coupled} steps")

    s_prev = 0.0
    while s_prev < 1.0 - 1e-12:
        # propose nominal step (roughly n_coupled steps), gentler near the end
        ds_nom = 1.0 / n_coupled
        if s_prev > 0.7:
            ds_nom *= 0.5
        s_try = min(1.0, s_prev + ds_nom)

        ok = False
        last_chance = False  # becomes True when we try exactly s_prev + min_ds

        while True:
            s = s_try
            a  = alpha0 + s * (alpha1 - alpha0)
            Ti = T_start + s * (T_f   - T_start)

            # update inlets
            f.fuel_inlet.T = float(Ti)
            parts = []
            for sp, val in target.items():
                amt = val if sp in ("H2", "NH3") else a * val
                if amt > 1e-16:
                    parts.append(f"{sp}:{amt}")
            comp_f = ", ".join(parts)
            f.fuel_inlet.X = comp_f

            # optional: tighten refinement late in the ramp
            if s > tighten_frac:
                # f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.10)
                f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05, prune=0.001)


            # try solve from previous state
            try:
                f.solve(loglevel, refine_grid=True, auto=False)
            except ct.CanteraError as err:
                Tmax = float(np.max(getattr(f, "T", np.array([0.]))))
                print(f"[ramp ?/?] solver error at s={s:.3f}: {err}; "
                      f"retrying with smaller step (Tmax={Tmax:.1f} K)")

                # halve step; if that drops below min_ds, do ONE last-chance try at exactly min_ds
                next_try = s_prev + 0.5 * (s - s_prev)
                if (next_try - s_prev) < min_ds:
                    if not last_chance:
                        s_try = s_prev + min_ds
                        last_chance = True
                        print(f"  ↳ step would fall below min_step; trying one last time at s={s_try:.3f}")
                        continue
                    else:
                        print("  ↳ last-chance step at min_step also failed; stopping ramp at last stable state.")
                        ok = False
                        break
                else:
                    s_try = next_try
                    continue

            Tmax = float(f.T.max())
            print(f"[ramp] s={s:.3f} → T_f={Ti:.1f} K, α={a:.3f}, fuel={comp_f}, Tmax={Tmax:.1f} K")

            if Tmax < T_cold_threshold:
                # cold or nearly extinguished → backtrack (reduce step)
                next_try = s_prev + 0.5 * (s - s_prev)
                if (next_try - s_prev) < min_ds:
                    if not last_chance:
                        s_try = s_prev + min_ds
                        last_chance = True
                        print(f"  ↳ near extinction; trying one last time at s={s_try:.3f} (min_step).")
                        continue
                    else:
                        print(f"  ↳ near extinction and min_step attempt failed; "
                              f"freezing progress at s={s_prev:.3f}.")
                        ok = False
                        break
                else:
                    s_try = next_try
                    print(f"  ↳ near extinction; reducing step to s={s_try:.3f} and retrying.")
                    continue
            else:
                ok = True
                break

        if not ok:
            print("[ramp] Stopping early at last sustained burning state.")
            break

        # accept the step
        s_prev = s

        # if we've reached s=1.0 (within tolerance), exit
        if 1.0 - s_prev < 1e-12:
            break

    # (Optional) reset default refinement if you tightened it late:
    # f.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.062)

    # ---------- post-process → CSV (grid, T, X/Y, HRR, CEM)
    grid, T_arr, X_arr, Y_arr = _extract_states(f, gas, P0)
    npts = T_arr.size
    species_names = gas.species_names

    # diagnostics
    gas_csp = csp.CanteraCSP(props["mech_file"])
    gas_csp.constP = P0
    gas_csp.jacobiantype = "full"

    HRR = np.empty(npts)
    CEM = np.empty(npts)
    P   = np.full(npts, P0)

    for i in range(npts):
        gas.TPX = T_arr[i], P0, X_arr[i, :]
        gas_csp.TPY = gas.T, P0, gas.Y
        lam, *_ = gas_csp.get_kernel()
        idx = np.argsort(np.abs(lam.real))[-1 - int(props["element_num"])]
        eig = lam.real[idx]
        CEM[i] = np.sign(eig) * np.log10(1.0 + abs(eig))
        HRR[i] = -gas.net_rates_of_progress.dot(gas.delta_enthalpy) / gas.density

    df = pd.DataFrame({
        "x": grid,
        "normalized_grid": grid / (grid[-1] - grid[0]),
        "T": T_arr,
        "P": P,
        "HRR": HRR,
        "CEM": CEM,
    })
    for j, sp in enumerate(species_names):
        df[f"Y_{sp}"] = Y_arr[:, j]
        df[f"X_{sp}"] = X_arr[:, j]

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"✓ wrote {out_csv}  ({npts} pts, {len(species_names)} species)")

# ------------------------ main ------------------------
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

    mech_path  = Path(args.mech_file).resolve()
    mech_name  = mech_path.stem

    # output root
    output_root = props_list[0].get("output_root", "outputs")
    out_root = Path(args.out_root or props_list[0].get("output_root", "outputs")) / mech_name
    check_path(out_root)

    if args.props_file:
        # yaml_path was resolved above when loading props_list
        props_src = yaml_path  # absolute path at this point
        try:
            shutil.copy2(props_src, output_root)
            print(f"[info] Copied props file → {output_root}")
        except Exception as e:
            print(f"[warn] Could not copy props file '{props_src}' → '{output_root}': {e}")
        

    # element_num override
    elem_override = args.element_num

    for i, p in enumerate(props_list, 1):
        p["mech_file"]   = str(mech_path)
        p["element_num"] = int(elem_override if elem_override is not None else p.get("element_num", 5))

        # required keys
        need = ("T", "P", "width", "fuel", "oxidizer", "Ufuel", "Uair")
        missing = [k for k in need if k not in p]
        if missing:
            raise ValueError(f"Missing keys in props for {p.get('case_name','<unknown>')}: {missing}")

        gas = ct.Solution(str(mech_path))
        out_csv = out_root / p["file_name"]
        try:
            run_case_simple(gas, p, out_csv, loglevel=args.loglevel, enable_soret=args.soret)
        except ct.CanteraError as e:
            print(f"Error on config {i}: {e}")

if __name__ == "__main__":
    main()