import os
import pandas as pd
import numpy as np
import cantera as ct
from post_1D_RCCE import run_CEQ_core, solve_eig_gas, highest_val_excl_0
import warnings
import traceback

def perturb_value(value: float, pct: float) -> float:
    """Perturb a value by ±pct (uniformly)."""
    factor = np.random.uniform(1 - pct, 1 + pct)
    return value * factor

def run_CEQ_and_extract(gas, main_species, init_mf, fuel, oxyd, T, P, dt, t_end, element_num):
    gas = run_CEQ_core(gas, main_species, init_mf, T, P, dt, t_end, moleFraction=True)
    out = {f"X_{sp}": x for sp, x in zip(gas.species_names, gas.X)}
    D, L, R = solve_eig_gas(gas)
    eigs = highest_val_excl_0(D, N_val=4, element_num=element_num)
    CEM = eigs[-1].real
    out["CEM"]  = np.sign(CEM) * np.log10(1 + abs(CEM))
    out["MF"]   = gas.mixture_fraction(fuel, oxyd)
    out["HRR"] = -1.0 * np.dot(gas.net_rates_of_progress, gas.delta_enthalpy) / gas.density
    return out

def run_uncertainty_analysis(
    gas, props, fig_dir, file_name,
    main_species, fuel, oxidizer,
    dt, time_end,
    uncertainty: dict,
    n_iterations: int = 1000
):
    """
    Serial version: loop over iterations in the main process.
    Any UserWarning (e.g. DVODE error test failures) is treated as an exception
    and causes that iteration to be skipped (no CSVs written).
    """
    mech_file  = props["mech_file"]
    inert_name = props["inert_specie"]
    
    df_cases = pd.read_csv(file_name)

    input_dir  = os.path.join(fig_dir, "inputs")
    output_dir = os.path.join(fig_dir, "outputs")
    os.makedirs(input_dir,  exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    for it in range(n_iterations):
        stage = "start iteration"
        try:
            stage = "rebuild gas0"
            with warnings.catch_warnings():
                warnings.simplefilter("error", UserWarning)

                gas0 = ct.Solution(mech_file)
                specs = gas0.species()
                idx   = gas0.species_index(inert_name)
                gas_i = ct.Solution(
                    thermo="IdealGas", kinetics="GasKinetics",
                    species=specs[:idx] + specs[idx+1:] + [specs[idx]],
                    reactions=gas0.reactions()
                )

                iter_inputs  = []
                iter_outputs = []

                for case_idx, row in df_cases.iterrows():
                    stage = f"perturb inputs for case {case_idx}"
                    T0     = row["T"]
                    T_pert = perturb_value(T0, uncertainty.get("T", 0.0))

                    extra = {"case_index": case_idx, "T": T_pert}
                    if "grid" in row:
                        extra["grid"] = row["grid"]

                    stage = f"build perturbed mole fractions for case {case_idx}"
                    X0 = {sp: row[f"X_{sp}"] for sp in main_species}
                    pert_mf = {
                        sp: perturb_value(X0[sp], uncertainty.get(sp, 0.0))
                        for sp in main_species
                    }

                    iter_inputs.append({
                        **extra,
                        "iteration": it,
                        **{f"X_{sp}": pert_mf[sp] for sp in main_species}
                    })

                    stage = f"run_CEQ_and_extract for case {case_idx}"
                    out = run_CEQ_and_extract(
                        gas_i, main_species, pert_mf,
                        fuel, oxidizer, T_pert,
                        props["P"], dt, time_end, props["element_num"]
                    )
                    iter_outputs.append({**extra, "iteration": it, **out})

            stage = "saving CSVs"
            in_df  = pd.DataFrame(iter_inputs)
            out_df = pd.DataFrame(iter_outputs)
            in_df.to_csv(os.path.join(input_dir,  f"uncertainty_inputs_iter{it}.csv"),  index=False)
            out_df.to_csv(os.path.join(output_dir, f"uncertainty_outputs_iter{it}.csv"), index=False)

            print(f"[Iteration {it}] Success")

        except Exception as e:
            print(f"[Iteration {it}] Skipped at stage '{stage}' due to: {e}")
            traceback.print_exc()
            continue

def main():
    from props_1D import generate_props
    props_list = generate_props("NH3_NP_sim")

    # build initial gas (you can keep this or ignore—it isn't used below)
    gas = ct.Solution(props_list[0]["mech_file"])
    inert = props_list[0]["inert_specie"]
    specs = gas.species()[:]
    i_last = gas.species_index(inert)
    gas = ct.Solution(
        thermo="IdealGas", kinetics="GasKinetics",
        species=specs[:i_last] + specs[i_last+1:] + [specs[i_last]],
        reactions=gas.reactions()
    )

    main_species = ["H2","O2","H2O","N2","NH3","NH"]
    
    for props in props_list:
        if props['case_name'] == "N_30":
            file_name = os.path.join(props['path'], props['file_name'])
            fig_dir   = f"figs/NH3_NP_CF/UQ_NH_analysis_noise_10/{props['case_name']}"

            run_uncertainty_analysis(
                gas, props, fig_dir, file_name,
                main_species, props["fuel"], props["oxidizer"],
                props["dt"], props["time_end"],
                uncertainty={
                    "T":   0.01,
                    "H2":  0.02,
                    "O2":  0.01,
                    "H2O": 0.02,
                    "N2":  0.005,
                    "NH3": 0.01,
                    # "OH":  0.1,
                    "NH": 0.1,
                    # "NH2": 0.1,
                },
                n_iterations=1000
            )

if __name__ == "__main__":
    main()