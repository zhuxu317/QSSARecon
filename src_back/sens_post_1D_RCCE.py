import os
import pandas as pd
import numpy as np
import cantera as ct
from post_1D_RCCE import run_CEQ_core, solve_eig_gas, highest_val_excl_0

def perturb_value(value, pct=0.001):
    """Return value perturbed upward by pct (0.1%)."""
    return value * (1 + pct)

def run_CEQ_and_extract(gas, main_species_names, initial_mole_fractions,
                        fuel, oxyd, T, P, dt, time_end, full_data,
                        element_num, moleFraction=True):
    """
    Run the simulation and return a flat dict of outputs:
      - X_<species> for each gas.species_names
      - CEM, MF, Qdot
    """
    gas = run_CEQ_core(gas, main_species_names,
                       initial_mole_fractions, T, P, dt, time_end,
                       moleFraction)

    out = {}
    # species mole fractions
    for sp, x in zip(gas.species_names, gas.X):
        out[f"X_{sp}"] = x

    # eigen‐value metric
    D, L, R = solve_eig_gas(gas)
    eigs = highest_val_excl_0(D, N_val=4, element_num=element_num)
    CEM = eigs[-1].real
    sign = np.sign(CEM)
    out["CEM"] = sign * np.log10(1 + abs(CEM))

    # mixture fraction & heat release
    out["MF"]   = gas.mixture_fraction(fuel, oxyd)
    out["Qdot"] = -1.0 * np.dot(gas.net_rates_of_progress,
                                gas.delta_enthalpy) / gas.density

    return out

def run_sensitivity_analysis(gas, props, fig_dir, file_name,
                             main_species_names, fuel, oxidizer,
                             dt, time_end,
                             moleFraction=True,
                             perturbation_pct=0.001):
    """
    For each row in the input CSV, run:
      1) baseline
      2) perturbed sims (T, P, each species) by +perturbation_pct
    Compute **relative sensitivity** = ((p - b)/b) / perturbation_pct
    and save to sensitivity_gradients.csv.
    """
    df_in = pd.read_csv(file_name)
    results = []
    os.makedirs(fig_dir, exist_ok=True)

    species_list = gas.species_names

    for idx, row in df_in.iterrows():
        T0 = row["T"]
        P0 = props["P"]
        X0 = {sp: row[f"X_{sp}"] for sp in main_species_names}

        extra = {"T": T0}
        if "grid" in row:
            extra["grid"] = row["grid"]

        # 1) baseline run
        base = run_CEQ_and_extract(
            gas, main_species_names, X0,
            fuel, oxidizer, T0, P0,
            dt, time_end, row,
            props["element_num"], moleFraction
        )

        params = ["T", "P"] + main_species_names

        for param in params:
            # build perturbed inputs
            if param == "T":
                Tp, Pp = perturb_value(T0, perturbation_pct), P0
                Xp = X0.copy()
            elif param == "P":
                Tp, Pp = T0, perturb_value(P0, perturbation_pct)
                Xp = X0.copy()
            else:
                Tp, Pp = T0, P0
                Xp = X0.copy()
                Xp[param] = perturb_value(X0[param], perturbation_pct)

            # perturbed run
            pert = run_CEQ_and_extract(
                gas, main_species_names, Xp,
                fuel, oxidizer, Tp, Pp,
                dt, time_end, row,
                props["element_num"], moleFraction
            )

            grad = {
                "case_index": idx,
                "parameter": param,
                **extra
            }

            # species relative sensitivities
            for sp in species_list:
                key = f"X_{sp}"
                b = base[key]
                p = pert[key]
                if b != 0:
                    raw = (p - b) / b
                    grad[key + "_sens"] = raw / perturbation_pct
                else:
                    grad[key + "_sens"] = np.nan

            # CEM, MF, Qdot relative sensitivities
            for metr in ("CEM", "MF", "Qdot"):
                b = base[metr]
                p = pert[metr]
                if b != 0:
                    raw = (p - b) / b
                    grad[metr + "_sens"] = raw / perturbation_pct
                else:
                    grad[metr + "_sens"] = np.nan

            results.append(grad)

    df_res = pd.DataFrame(results)
    out_path = os.path.join(fig_dir, f"{main_species_names[-1]}_sensitivity_gradients_relative_{perturbation_pct}.csv")
    df_res.to_csv(out_path, index=False)
    print(f"Saved relative sensitivities to {out_path}")

# ----------------------------------
def main():
    from props_1D import generate_props
    props = generate_props("NH3_NP_sim")[0]
    gas = ct.Solution(props["mech_file"])
    # reorder so inert last (as your code does):
    last = props["inert_specie"]
    specs = gas.species()[:]
    i_last = gas.species_index(last)
    gas = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                      species=specs[:i_last]+specs[i_last+1:]+[specs[i_last]],
                      reactions=gas.reactions())

    main_species = ["H2","O2","H2O","N2","NH3"]
    perturbation = 1e-3
    # file_name = os.path.join(props["path"], props["file_name"])
    file_name = "figs/QSSA/SIM_results/NH3_NP_CF_reduce_OH/AH82/predicted_X.csv"
    fig_dir   = "figs/NH3_NP_CF/sensitivity_analysis"

    run_sensitivity_analysis(
        gas, props, fig_dir, file_name,
        main_species, props["fuel"], props["oxidizer"],
        props["dt"], props["time_end"], moleFraction=True, perturbation_pct= perturbation
    )

if __name__ == "__main__":
    main()