import argparse
import cantera as ct
import numpy as np
import reactorch as rt
import torch

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--yaml', required=True, help="yaml mech file")
args = parser.parse_args()

print(f'Processing mechanism {args.yaml!r}\n')

# build the Solution
sol = rt.Solution(
    mech_yaml= args.yaml,
    device     = torch.device('cpu'),
    vectorize  = True,
    is_clip    = False,
    is_norm    = False,
    is_wdot_vec= False,
)
gas = sol.gas

# set a dummy state so all species exist
gas.TPY = 1200, ct.one_atm, np.ones(gas.n_species)/gas.n_species

# collect everything
molecular_weights       = gas.molecular_weights.tolist()
reactant_stoich_coeffs  = gas.reactant_stoich_coeffs()
product_stoich_coeffs   = gas.product_stoich_coeffs()
Arrhenius_coeffs        = sol.Arrhenius_coeffs
efficiencies_coeffs     = sol.efficiencies_coeffs
reactant_orders         = sol.reactant_orders
is_reversible           = sol.is_reversible

# optional Troe / low-pressure parameters
has_troe = len(sol.list_reaction_type4_Troe) > 0
if has_troe:
    Arr_A0   = sol.Arrhenius_A0
    Arr_b0   = sol.Arrhenius_b0
    Arr_Ea0  = sol.Arrhenius_Ea0
    Troe_A   = sol.Troe_A
    Troe_T1  = sol.Troe_T1
    Troe_T2  = sol.Troe_T2
    Troe_T3  = sol.Troe_T3

# now print them all
print("=== molecular_weights ===")
print(np.array(molecular_weights))

print("\n=== reactant_stoich_coeffs ===")
print(np.array(reactant_stoich_coeffs))

print("\n=== product_stoich_coeffs ===")
print(np.array(product_stoich_coeffs))

print("\n=== reactant_orders ===")
print(np.array(reactant_orders))

print("\n=== is_reversible ===")
print(np.array(is_reversible))

print("\n=== Arrhenius_coeffs (shape {}) ===".format(Arrhenius_coeffs.shape))
print(Arrhenius_coeffs)

print("\n=== efficiencies_coeffs (shape {}) ===".format(efficiencies_coeffs.shape))
print(efficiencies_coeffs)

if has_troe:
    print("\n--- Troe / low-pressure parameters detected ---")
    print("\n=== Arrhenius_A0 ===")
    print(np.array(Arr_A0))
    print("\n=== Arrhenius_b0 ===")
    print(np.array(Arr_b0))
    print("\n=== Arrhenius_Ea0 ===")
    print(np.array(Arr_Ea0))
    print("\n=== Troe_A ===")
    print(np.array(Troe_A))
    print("\n=== Troe_T1 ===")
    print(np.array(Troe_T1))
    print("\n=== Troe_T2 ===")
    print(np.array(Troe_T2))
    print("\n=== Troe_T3 ===")
    print(np.array(Troe_T3))

# finally save
npz_name = args.yaml if args.yaml.endswith('.npz') else args.yaml + '.npz'
print(f"\nSaving to {npz_name!r} …\n")
if not has_troe:
    np.savez(
        npz_name,
        molecular_weights       = molecular_weights,
        reactant_stoich_coeffs  = reactant_stoich_coeffs,
        product_stoich_coeffs   = product_stoich_coeffs,
        reactant_orders         = sol.reactant_orders,
        is_reversible           = sol.is_reversible,
        Arrhenius_coeffs        = Arrhenius_coeffs,
        efficiencies_coeffs     = efficiencies_coeffs,
    )
else:
    np.savez(
        npz_name,
        molecular_weights       = molecular_weights,
        reactant_stoich_coeffs  = reactant_stoich_coeffs,
        product_stoich_coeffs   = product_stoich_coeffs,
        reactant_orders         = sol.reactant_orders,
        is_reversible           = sol.is_reversible,
        Arrhenius_coeffs        = Arrhenius_coeffs,
        efficiencies_coeffs     = efficiencies_coeffs,
        Arrhenius_A0            = Arr_A0,
        Arrhenius_b0            = Arr_b0,
        Arrhenius_Ea0           = Arr_Ea0,
        Troe_A                  = Troe_A,
        Troe_T1                 = Troe_T1,
        Troe_T2                 = Troe_T2,
        Troe_T3                 = Troe_T3,
    )

print("Done.")  