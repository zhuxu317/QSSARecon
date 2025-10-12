# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np

def calculate_flame_speed(mechanism_file, fuel_species, T=300, P=101325, phi_range=[0.6, 1.4], n_points=9):
    """
    Calculate laminar flame speed for a range of equivalence ratios
    
    Parameters:
    mechanism_file: Path to YAML mechanism file
    fuel_species: Name of fuel species in the mechanism
    T: Temperature in K
    P: Pressure in Pa
    phi_range: Range of equivalence ratios [min, max]
    n_points: Number of points to calculate
    
    Returns:
    phi_values: Array of equivalence ratios
    flame_speeds: Array of laminar flame speeds (m/s)
    """
    
    # Load mechanism
    gas = ct.Solution(mechanism_file)
    
    # Define equivalence ratios
    phi_values = np.linspace(phi_range[0], phi_range[1], n_points)
    flame_speeds = []
    
    print("Calculating flame speeds for mechanism: {}".format(mechanism_file))
    print("Fuel species: {}".format(fuel_species))
    print("Temperature: {} K, Pressure: {} Pa".format(T, P))
    print("-" * 50)
    
    for phi in phi_values:
        try:
            # Set equivalence ratio with air
            # For hydrocarbon fuel: CxHy + (x + y/4)(O2 + 3.76*N2) = x*CO2 + y/2*H2O + 3.76*(x + y/4)*N2
            # Assuming C12H23 fuel
            if fuel_species == 'C12H23':
                # C12H23 + 17.75 O2 + 66.74 N2 = 12 CO2 + 11.5 H2O + 66.74 N2
                fuel_moles = 1.0
                o2_moles = 17.75 / phi  # Divide by phi for equivalence ratio
                n2_moles = o2_moles * 3.76  # Air composition
                
                gas.set_equivalence_ratio(phi, fuel_species, 'O2:1, N2:3.76')
            else:
                # For other fuels, use generic approach
                gas.set_equivalence_ratio(phi, fuel_species, 'O2:1, N2:3.76')
            
            gas.TP = T, P
            
            # Create flame object
            flame = ct.FreeFlame(gas, width=0.05)
            flame.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
            
            # Solve flame
            flame.solve(loglevel=0, auto=True)
            
            # Get flame speed
            Su = flame.velocity[0]  # m/s
            flame_speeds.append(Su)
            
            print("φ = {:.2f}: Su = {:.4f} m/s".format(phi, Su))
            
        except Exception as e:
            print("φ = {:.2f}: Failed - {}".format(phi, str(e)))
            flame_speeds.append(np.nan)
    
    return phi_values, np.array(flame_speeds)

def main():
    # Define mechanisms to compare
    mechanisms = {
        'Three-body Pure (Original)': 'jeta-23steps_threebody_pure.yaml',
        'Three-body Fixed': 'jeta-23steps_threebody_fixed.yaml',
        'Cleaned (Remove FORD)': 'jeta-23steps_clean_remove_fixed.yaml'
    }
    
    # Parameters
    fuel_species = 'C12H23'
    T = 300  # K
    P = 101325  # Pa (1 atm)
    phi_range = [0.8, 1.2]
    n_points = 5
    
    results = {}
    
    for name, mechanism_file in mechanisms.items():
        try:
            print("\n" + "="*60)
            print("MECHANISM: {}".format(name))
            print("="*60)
            
            phi_values, flame_speeds = calculate_flame_speed(
                mechanism_file, fuel_species, T, P, phi_range, n_points
            )
            
            results[name] = {
                'phi': phi_values,
                'Su': flame_speeds,
                'mechanism_file': mechanism_file
            }
            
        except Exception as e:
            print("Error with mechanism {}: {}".format(name, str(e)))
            results[name] = None
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    for name, data in results.items():
        if data is not None:
            valid_speeds = data['Su'][~np.isnan(data['Su'])]
            if len(valid_speeds) > 0:
                print("{}:".format(name))
                print("  Mean flame speed: {:.4f} m/s".format(np.mean(valid_speeds)))
                print("  Max flame speed: {:.4f} m/s".format(np.max(valid_speeds)))
                print("  Min flame speed: {:.4f} m/s".format(np.min(valid_speeds)))
                
                # Find stoichiometric conditions (phi closest to 1.0)
                phi_stoich_idx = np.argmin(np.abs(data['phi'] - 1.0))
                Su_stoich = data['Su'][phi_stoich_idx]
                if not np.isnan(Su_stoich):
                    print("  Stoichiometric flame speed: {:.4f} m/s".format(Su_stoich))
            else:
                print("{}: No valid flame speeds calculated".format(name))
        else:
            print("{}: Failed to calculate".format(name))
    
    return results

if __name__ == '__main__':
    results = main()