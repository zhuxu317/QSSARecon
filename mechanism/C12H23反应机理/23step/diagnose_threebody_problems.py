#!/usr/bin/env python3
import cantera as ct
import numpy as np

def analyze_reaction_thermodynamics(gas, reaction_index):
    """Analyze thermodynamics of a specific reaction"""
    r = gas.reactions()[reaction_index]
    print(f"\nReaction {reaction_index+1}: {r.equation}")
    
    # Set reasonable conditions
    gas.TPX = 1000, ct.one_atm, 'C12H23:1, O2:17.75, N2:66.74'
    
    # Get reaction properties
    try:
        delta_H = gas.delta_enthalpy[reaction_index] / 1000  # kJ/mol
        delta_S = gas.delta_entropy[reaction_index] / 1000   # kJ/mol/K
        delta_G = gas.delta_gibbs[reaction_index] / 1000     # kJ/mol
        
        print(f"  ΔH = {delta_H:.2f} kJ/mol")
        print(f"  ΔS = {delta_S:.6f} kJ/mol/K")
        print(f"  ΔG = {delta_G:.2f} kJ/mol")
        
        # Check if reaction is thermodynamically feasible
        if delta_G > 100:
            print("  WARNING: Very positive ΔG - thermodynamically unfavorable")
        elif delta_G < -500:
            print("  WARNING: Very negative ΔG - might be too fast")
        
    except Exception as e:
        print(f"  ERROR: Could not calculate thermodynamics: {e}")
    
    # Check mass balance
    reactant_atoms = {}
    product_atoms = {}
    
    for species, coeff in r.reactants.items():
        for element, count in gas.species(species).composition.items():
            reactant_atoms[element] = reactant_atoms.get(element, 0) + coeff * count
    
    for species, coeff in r.products.items():
        for element, count in gas.species(species).composition.items():
            product_atoms[element] = product_atoms.get(element, 0) + coeff * count
    
    print("  Mass balance check:")
    all_elements = set(reactant_atoms.keys()) | set(product_atoms.keys())
    balanced = True
    for element in all_elements:
        r_count = reactant_atoms.get(element, 0)
        p_count = product_atoms.get(element, 0)
        if abs(r_count - p_count) > 1e-10:
            print(f"    {element}: Reactants={r_count}, Products={p_count} - UNBALANCED!")
            balanced = False
        else:
            print(f"    {element}: {r_count} - balanced")
    
    if not balanced:
        print("  ERROR: Reaction is not mass balanced!")
    
    return balanced

def main():
    print("Diagnosing Three-body Pure mechanism problems...")
    
    try:
        gas = ct.Solution('jeta-23steps_threebody_pure.yaml')
        print(f"Loaded mechanism with {len(gas.species_names)} species and {len(gas.reactions())} reactions")
        
        # Focus on the first few problematic reactions
        problematic_reactions = [0, 1, 2]  # First 3 reactions that look suspicious
        
        for i in problematic_reactions:
            balanced = analyze_reaction_thermodynamics(gas, i)
            if not balanced:
                print(f"  >>> CRITICAL: Reaction {i+1} has mass balance issues!")
        
        # Test basic flame initialization
        print("\n" + "="*60)
        print("Testing flame initialization...")
        
        gas.TPX = 300, ct.one_atm, 'C12H23:1, O2:17.75, N2:66.74'
        print(f"Initial mixture: T={gas.T}K, P={gas.P/ct.one_atm:.1f}atm")
        print(f"Composition: {dict(zip(gas.species_names, gas.X))}")
        
        # Try to create flame object
        try:
            flame = ct.FreeFlame(gas, width=0.05)
            print("Flame object created successfully")
            
            # Try initial solution
            flame.set_refine_criteria(ratio=5, slope=0.1, curve=0.2)  # Relaxed criteria
            print("Attempting initial solution...")
            flame.solve(loglevel=1, auto=False)
            print("Initial solution successful!")
            
        except Exception as e:
            print(f"Flame initialization failed: {e}")
        
    except Exception as e:
        print(f"Failed to load mechanism: {e}")

if __name__ == '__main__':
    main()