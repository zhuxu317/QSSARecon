#!/usr/bin/env python3
import os
import re
import sys
import numpy as np
import pandas as pd

# ================== USER CONFIG ==================
MECHANISM_NAME = "H2_pNOX_15_94_TC"
MECH_PATH      = f"mechanism/{MECHANISM_NAME}.yaml"

RAW_ROOT       = f"cases/H2_DNS/DNS2D/raw/{MECHANISM_NAME}"
PROCESSED_ROOT = f"cases/H2_DNS/DNS2D/processed/{MECHANISM_NAME}"

# Only process files matching this suffix (set to "" to process all CSVs)
CSV_SUFFIX_FILTER = ".csv"  # e.g., ".csv" or "_slice_filtered.csv"
# =================================================

# ---- Optional: import Cantera lazily ----
try:
    import cantera as ct
except Exception as e:
    sys.exit(
        "Cantera import failed. Install with `conda install -c conda-forge cantera` "
        "or `pip install cantera`. Error:\n" + str(e)
    )

# ---- Helpers (rename headers to friendly names) ----
unit_suffix_re = re.compile(r"\s*\[[^\]]*\]\s*$")
y_species_re   = re.compile(r"^Y\(([^)]+)\)$")

def normalize_col(name: str) -> str:
    """
    Strip unit suffixes, map 'temp' -> 'T', and 'Y(OH)' -> 'Y_OH'.
    Leaves other names as-is (minus units).
    """
    base = unit_suffix_re.sub("", name.strip())
    if base.lower() == "temp":
        return "T"
    m = y_species_re.match(base)
    if m:
        species = m.group(1)
        species = re.sub(r"[^A-Za-z0-9]+", "_", species).strip("_")
        return f"Y_{species}"
    return base

def dedupe_cols(cols):
    """Avoid duplicate column names by appending __dupN."""
    seen = {}
    final = []
    for c in cols:
        if c not in seen:
            seen[c] = 0
            final.append(c)
        else:
            seen[c] += 1
            final.append(f"{c}__dup{seen[c]}")
    return final

def process_csv(in_csv: str, out_csv: str, gas, species_names):
    """Load one CSV, rename headers, compute X_* and HRR, and save."""
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    df = pd.read_csv(in_csv)

    # Normalize column names (strip units, temp->T, Y(OH)->Y_OH)
    original_cols = list(df.columns)
    new_cols = [normalize_col(c) for c in original_cols]
    df.columns = dedupe_cols(new_cols)

    # Build mapping from CSV Y_* -> mechanism species
    csv_y_cols = [c for c in df.columns if c.startswith("Y_")]
    csv_to_mech = {}
    for col in csv_y_cols:
        sp = col[2:]  # after 'Y_'
        if sp in species_names:
            csv_to_mech[col] = sp
        else:
            cand = sp.replace("__", "_")
            if cand in species_names:
                csv_to_mech[col] = cand
            else:
                # Not fatal; just ignore this Y_ column if not in mechanism
                # (keep in df; not used for TPY)
                pass

    # Prepare outputs
    X_cols = [f"X_{s}" for s in species_names]
    X_data = {col: np.zeros(len(df), dtype=float) for col in X_cols}
    HRR = np.zeros(len(df), dtype=float)

    # Constant pressure (Pa); add/update P column
    P = ct.one_atm
    df["P"] = np.full(len(df), P, dtype=float)

    # Sanity: T available?
    if "T" not in df.columns:
        raise KeyError(
            "Column 'T' not found. If your file has 'temp', it should have been "
            "auto-mapped to 'T' by normalize_col()."
        )

    # Iterate rows (robust, readable; okay for typical slice sizes)
    for i, row in df.iterrows():
        T = float(row["T"])

        # Build Y_dict from available CSV columns that map to mechanism species
        Y_dict = {}
        for col, mech_name in csv_to_mech.items():
            val = row[col]
            if pd.notna(val):
                Y_dict[mech_name] = max(float(val), 0.0)

        if not Y_dict:
            # Nothing to set; leave zeros for X_* and HRR
            continue

        # Normalize mass fractions to 1 (defensive)
        y_sum = sum(Y_dict.values())
        if y_sum > 0.0:
            inv = 1.0 / y_sum
            for k in list(Y_dict.keys()):
                Y_dict[k] *= inv

        try:
            gas.TPY = T, P, Y_dict
        except Exception:
            # If TPY fails (e.g., numerical issues), skip this row
            continue

        # Mole fractions
        X_vec = gas.X
        for j, s in enumerate(species_names):
            X_data[f"X_{s}"][i] = X_vec[j]

        # HRR (W/kg): - Σ(ω̇_r * ΔH_r) / ρ
        try:
            rop = gas.net_rates_of_progress   # kmol/m^3/s
            dH  = gas.delta_enthalpy          # J/kmol
            rho = gas.density                 # kg/m^3
            HRR[i] = -1.0 * float(np.dot(rop, dH)) / float(rho)
        except Exception:
            HRR[i] = 0.0

        if (i + 1) % 2000 == 0:
            print(f"  processed {i+1}/{len(df)} rows…")

    # Append new columns
    for col, arr in X_data.items():
        df[col] = arr
    df["HRR"] = HRR

    # Save
    df.to_csv(out_csv, index=False)
    print(f"[ok] {out_csv}")

def main():
    if not os.path.exists(MECH_PATH):
        sys.exit(f"Mechanism not found: {MECH_PATH}")

    # Prepare Cantera gas
    gas = ct.Solution(MECH_PATH)
    species_names = gas.species_names

    # Walk RAW_ROOT and process matching CSVs
    if not os.path.isdir(RAW_ROOT):
        sys.exit(f"Raw root not found: {RAW_ROOT}")

    n_files = 0
    for root, _, files in os.walk(RAW_ROOT):
        for fname in files:
            if not fname.lower().endswith(".csv"):
                continue
            if CSV_SUFFIX_FILTER and not fname.endswith(CSV_SUFFIX_FILTER):
                continue

            in_csv = os.path.join(root, fname)

            # Mirror relative path under PROCESSED_ROOT
            rel = os.path.relpath(in_csv, RAW_ROOT)  # e.g., "ka1p0/plt09600/DNS_y064_slice_filtered.csv"
            out_csv = os.path.join(PROCESSED_ROOT, rel)

            # Create dest dir and process
            print(f"-> processing: {in_csv}")
            os.makedirs(os.path.dirname(out_csv), exist_ok=True)
            process_csv(in_csv, out_csv, gas, species_names)
            n_files += 1

    if n_files == 0:
        print("[info] No CSV files matched. Check RAW_ROOT / CSV_SUFFIX_FILTER.")
    else:
        print(f"[done] Processed {n_files} file(s).")

if __name__ == "__main__":
    main()
