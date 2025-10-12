#!/usr/bin/env python3
import cantera as ct
import yaml

def fix_threebody_mechanism():
    """Fix the problematic reactions in the three-body mechanism"""
    
    # Load the problematic mechanism
    gas = ct.Solution('jeta-23steps_threebody_pure.yaml')
    
    print("Original three-body mechanism problems identified:")
    print("1. Reaction 1: C12H23 + N2 (+M) <=> 12 CH + 11 H + N2 (+M)")
    print("   - Extremely negative ΔG = -125 million kJ/mol")
    print("   - Should be: N2 + C12H23 => 12 CH + 11 H + N2")
    print("2. Reaction 2: CH + H2 + N2 (+M) <=> CH + 2 NH (+M)")
    print("   - ΔG = -5 million kJ/mol")
    print("   - Should be: CH + H2 + N2 => CH + 2 NH")
    print("3. Reaction 3: HO2 + N2 + O (+M) <=> H + 2 NO + O (+M)")
    print("   - ΔG = -11 million kJ/mol") 
    print("   - Should be: O + N2 + HO2 => 2 NO + H + O")
    
    # Load the working cleaned mechanism as reference
    gas_clean = ct.Solution('jeta-23steps_clean_remove_fixed.yaml')
    
    print("\nComparing with working cleaned mechanism...")
    print("Cleaned mechanism reactions:")
    for i in range(min(5, len(gas_clean.reactions()))):
        print(f"{i+1}: {gas_clean.reactions()[i].equation}")
    
    # Create a corrected YAML by modifying the three-body version
    print("\nCreating corrected three-body mechanism...")
    
    # Load the YAML file directly
    with open('jeta-23steps_threebody_pure.yaml', 'r') as f:
        data = yaml.safe_load(f)
    
    # Fix the problematic reactions by removing the (+M) notation
    # and correcting the reaction equations
    if 'reactions' in data:
        reactions = data['reactions']
        
        # Fix Reaction 1: Remove (+M) and make it irreversible
        reactions[0]['equation'] = 'N2 + C12H23 => 12 CH + 11 H + N2'
        if 'efficiencies' in reactions[0]:
            del reactions[0]['efficiencies']
        if 'type' in reactions[0]:
            del reactions[0]['type']
        
        # Fix Reaction 2: Remove (+M) and correct order
        reactions[1]['equation'] = 'CH + H2 + N2 => CH + 2 NH'
        if 'efficiencies' in reactions[1]:
            del reactions[1]['efficiencies']
        if 'type' in reactions[1]:
            del reactions[1]['type']
        
        # Fix Reaction 3: Remove (+M) and correct equation
        reactions[2]['equation'] = 'O + N2 + HO2 => 2 NO + H + O'
        if 'efficiencies' in reactions[2]:
            del reactions[2]['efficiencies']
        if 'type' in reactions[2]:
            del reactions[2]['type']
    
    # Write the corrected mechanism
    with open('jeta-23steps_threebody_fixed.yaml', 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    
    print("Created corrected mechanism: jeta-23steps_threebody_fixed.yaml")
    
    # Test the corrected mechanism
    try:
        gas_fixed = ct.Solution('jeta-23steps_threebody_fixed.yaml')
        print(f"✓ Corrected mechanism loads successfully")
        print(f"  Species: {len(gas_fixed.species_names)}")
        print(f"  Reactions: {len(gas_fixed.reactions())}")
        
        # Test thermodynamics of first few reactions
        gas_fixed.TPX = 1000, ct.one_atm, 'C12H23:1, O2:17.75, N2:66.74'
        for i in range(3):
            delta_G = gas_fixed.delta_gibbs[i] / 1000  # kJ/mol
            print(f"  Reaction {i+1} ΔG = {delta_G:.2f} kJ/mol")
        
        return True
        
    except Exception as e:
        print(f"✗ Error loading corrected mechanism: {e}")
        return False

if __name__ == '__main__':
    fix_threebody_mechanism()