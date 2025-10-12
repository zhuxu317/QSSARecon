import cantera as ct
gas = ct.Solution('san_diego_nitrogen.yaml')
# gas = ct.Solution('san_diego.yaml')
gas.check()
for k, rxn in enumerate(gas.reactions()):
    prob = ct.check_reaction_balance(rxn, quiet=True)
    if prob:
        print(f"[{k}] {gas.reaction_equation(k)}  →  {prob}")