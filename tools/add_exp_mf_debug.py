#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, argparse, yaml
import pandas as pd
import numpy as np
import cantera as ct

# ------------------------ helpers ------------------------
def _to_pa(pdict):
    if "P" in pdict:      return float(pdict["P"])
    if "P_bar" in pdict:  return float(pdict["P_bar"]) * 1e5
    if "P_atm" in pdict:  return float(pdict["P_atm"]) * 101325.0
    return 101325.0  # default 1 atm

def _load_run_map(case_dir):
    """
    Read cases/<case>/props.yaml and return:
      defaults: (fuel, oxidizer, P[Pa])
      run_map : {<run_id> -> (fuel, oxidizer, P[Pa]),
                 <file_name_stem> -> (fuel, oxidizer, P[Pa])}
    """
    props_path = os.path.join(case_dir, "props.yaml")
    with open(props_path, "r") as f:
        data = yaml.safe_load(f) or {}

    defaults = data.get("defaults", {}) or {}
    cfgs = data.get("configs", []) or []

    def_fuel = defaults.get("fuel")
    def_ox   = defaults.get("oxidizer")
    def_P    = _to_pa(defaults)

    # if defaults missing, borrow from first config
    if (def_fuel is None or def_ox is None) and cfgs:
        first = cfgs[0] or {}
        def_fuel = def_fuel or first.get("fuel")
        def_ox   = def_ox   or first.get("oxidizer")
        def_P    = _to_pa({**defaults, **first})

    if def_fuel is None or def_ox is None:
        raise ValueError("Could not find default 'fuel'/'oxidizer' in props.yaml")

    run_map = {}
    for cfg in cfgs:
        if not cfg:
            continue
        fuel = cfg.get("fuel", def_fuel)
        ox   = cfg.get("oxidizer", def_ox)
        P    = _to_pa({**defaults, **cfg})

        # map by run_id
        run_id = str(cfg.get("run_id", "")).strip()
        if run_id:
            run_map[run_id] = (fuel, ox, P)

        # also map by file_name stem (e.g., AH64 from AH64.csv)
        fname = str(cfg.get("file_name", "")).strip()
        if fname:
            stem, _ = os.path.splitext(fname)
            if stem:
                run_map[stem] = (fuel, ox, P)

    return (def_fuel, def_ox, def_P), run_map

def _brief_stats(vec):
    if vec.size == 0:
        return "n=0"
    return "n={} sum={:.6g} min={:.3g} max={:.3g}".format(vec.size, vec.sum(), vec.min(), vec.max())

def _beta(g: ct.Solution):
    # Bilger β = 2 Z_C + 0.5 Z_H − Z_O   (mass-based elemental fractions)
    ZC = g.elemental_mass_fraction("C") if "C" in g.element_names else 0.0
    ZH = g.elemental_mass_fraction("H") if "H" in g.element_names else 0.0
    ZO = g.elemental_mass_fraction("O") if "O" in g.element_names else 0.0
    return 2.0*ZC + 0.5*ZH - ZO

# ------------------ core computation ---------------------
def compute_mixture_fraction_for_csv(
    csv_path,
    gas,
    fuel,
    oxidizer,
    P_pa,
    out_path=None,
    head=5,
    verbose=True,
    balance_species="N2"
):
    """
    Read an EXP csv with columns like:
      x, X_NO(ppm), T, X_O2, X_N2, X_NH3, X_H2O, X_H2, ...
    Per row:
      - X_NO /= 1e6 (ppm -> mole fraction)
      - clip negatives to 0
      - optionally balance remainder into `balance_species` (default N2)
      - gas.TPX(T,P,X); MF = gas.mixture_fraction(fuel, oxidizer)
    """
    df = pd.read_csv(csv_path)
    xcols = [c for c in df.columns if c.startswith("X_")]
    if "T" not in df.columns:
        raise ValueError("{}: missing 'T' column.".format(csv_path))
    if not xcols:
        raise ValueError("{}: no X_* columns found.".format(csv_path))

    # Quick header info
    if verbose:
        base = os.path.basename(csv_path)
        print(f"  -> Reading {base}: rows={len(df)}, X-cols={len(xcols)}")
        xv = df[xcols].iloc[0].values.astype(float)
        print("     first-row X-stats (raw): {}".format(_brief_stats(xv)))
        if "X_NO" in df.columns:
            print("     note: treating X_NO as ppm and dividing by 1e6")
        print(f"     using fuel='{fuel}'  oxidizer='{oxidizer}'  P={P_pa/101325.0:.3f} atm")

    # β_F / β_O for this (fuel, oxidizer)
    try:
        gas.TPX = 300.0, P_pa, fuel
        beta_F = _beta(gas)
        gas.TPX = 300.0, P_pa, oxidizer
        beta_O = _beta(gas)
        if verbose:
            print(f"     β_F={beta_F:.6g}, β_O={beta_O:.6g}, Δβ={beta_F-beta_O:.6g}")
    except Exception as e:
        if verbose:
            print(f"     [warn] could not compute β_F/β_O: {e}")

    can_balance = (balance_species is not None) and (balance_species in gas.species_names)
    MF = np.empty(len(df), dtype=float)

    for i, row in df.iterrows():
        T = float(row["T"])

        X_dict = {}
        for c in xcols:
            sp = c[2:]  # strip "X_"
            val = float(row[c]) if pd.notnull(row[c]) else 0.0
            if sp.upper() == "NO":
                val = val / 1e6  # ppm -> mole fraction
            if val < 0.0:
                val = 0.0
            X_dict[sp] = val

        sumX = sum(X_dict.values())

        # (optional) balance small deficit into inert (N2)
        if can_balance:
            rem = 1.0 - sumX
            if rem > 1e-8:
                X_dict[balance_species] = X_dict.get(balance_species, 0.0) + rem
                sumX = sum(X_dict.values())

        try:
            gas.TPX = T, P_pa, X_dict
        except Exception as e:
            raise RuntimeError("{}: row {} -> bad composition for Cantera: {}".format(csv_path, i, e))

        MF[i] = gas.mixture_fraction(fuel, oxidizer)

        if verbose and i < head:
            arr = np.array(list(X_dict.values()))
            print("     [row {:>3}] T={:.1f} K  ΣX={:.6f}  minX={:.3g}  maxX={:.3g}  MF={:.6g}".format(
                i, T, sumX, arr.min() if arr.size else 0.0, arr.max() if arr.size else 0.0, MF[i]
            ))

    df["MF"] = MF
    if out_path is None:
        out_path = csv_path
    df.to_csv(out_path, index=False)
    return out_path, len(df)

# ------------------------- CLI ---------------------------
def main():
    ap = argparse.ArgumentParser(description="Append MF to EXP CSVs using Cantera.mixture_fraction; X_NO in ppm (/1e6).")
    ap.add_argument("--case_dir", required=True, help="cases/<case> (e.g., cases/NH3_KAUST_NP_1bar)")
    ap.add_argument("--mech_file", required=True, help="Mechanism file (.yaml or .cti)")
    ap.add_argument("--exp_dir", default="Exp", help="EXP subdir under case_dir (default: Exp)")
    ap.add_argument("--pattern", default="*.csv", help="Filename pattern (default: *.csv)")
    ap.add_argument("--out_dir", default=None, help="If set, write modified CSVs here; else overwrite in place")
    ap.add_argument("--head", type=int, default=5, help="Debug first N rows per file")
    ap.add_argument("--quiet", action="store_true", help="Less verbose")
    ap.add_argument("--no_balance", action="store_true", help="Disable balancing remainder into N2")
    args = ap.parse_args()

    case_dir = os.path.abspath(args.case_dir)
    exp_dir  = os.path.join(case_dir, args.exp_dir)
    if not os.path.isdir(exp_dir):
        raise FileNotFoundError("EXP directory not found: {}".format(exp_dir))

    (def_fuel, def_ox, def_P), run_map = _load_run_map(case_dir)
    print("[info] case_dir  =", case_dir)
    print("[info] mech_file =", os.path.abspath(args.mech_file))
    print("[info] default fuel     =", def_fuel)
    print("[info] default oxidizer =", def_ox)
    print("[info] default pressure = {:.3f} atm".format(def_P/101325.0))

    gas = ct.Solution(args.mech_file)

    from glob import glob
    files = sorted(glob(os.path.join(exp_dir, args.pattern)))
    if not files:
        raise FileNotFoundError("No CSVs matched {} in {}".format(args.pattern, exp_dir))

    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)

    balance_species = None if args.no_balance else "N2"

    for fpath in files:
        fname = os.path.basename(fpath)
        stem, _ = os.path.splitext(fname)

        # choose per-run fuel/oxidizer/P: prefer run_id==stem; else defaults
        fuel, oxidizer, P_pa = run_map.get(stem, (def_fuel, def_ox, def_P))
        if not args.quiet:
            if stem in run_map:
                print(f"[info] {stem}: matched run_id/file_name → using fuel='{fuel}', oxidizer='{oxidizer}', P={P_pa/101325.0:.3f} atm")
            else:
                print(f"[warn] {stem}: no run_id match; falling back to defaults")

        out_path = os.path.join(args.out_dir, fname) if args.out_dir else None
        target, n = compute_mixture_fraction_for_csv(
            fpath, gas, fuel, oxidizer, P_pa, out_path,
            head=args.head, verbose=(not args.quiet),
            balance_species=balance_species
        )
        print("[done] {}: wrote MF for {} rows -> {}".format(fname, n, target))

if __name__ == "__main__":
    main()
