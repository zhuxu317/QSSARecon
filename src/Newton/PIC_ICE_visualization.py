import os, sys, time
sys.path.append(os.path.abspath('src'))
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.linalg import lstsq
from functools import lru_cache
import pandas as pd
from run_1D_CEQ import run_CEQ_core
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


epsilon = 1e-32
R = 8314.462175

import cantera as ct
import numpy as np

class PIC:
    def __init__(self, mech, constrained_species, use_element_constraints=True):
        """
        Initialize the PIC class.

        Parameters:
        - mech: Path to the Cantera mechanism file.
        - constrained_species: List of species to constrain.
        - use_element_constraints: If True, uses element constraints (Method 1).
                                    If False, uses Method 2 where Nc = len(constrained_species) + 1.
        """
        self.gas = ct.Solution(mech)
        self.specified_species = constrained_species
        self.use_element_constraints = use_element_constraints

        self.Ns = self.gas.n_species
        self.Ne = self.gas.n_elements
        self.Nr = self.gas.n_reactions

        print(f"Number of species (Ns): {self.Ns}")
        print(f"Number of elements (Ne): {self.Ne}")
        print(f"Number of reactions (Nr): {self.Nr}")
        
        # Select Nc based on the chosen method
        if use_element_constraints:
            self.Nc = self.Ne + len(constrained_species)  # Method 1
        else:
            self.Nc = len(constrained_species) + 1        # Method 2s
        self.sp_MW = self.gas.molecular_weights
        self.sp_names = np.array(self.gas.species_names)
        print("species_names=", self.sp_names)

        self.sp_coeffs = np.zeros((self.gas.n_species, 1 + 7 + 7))

        for i, sp in enumerate(self.gas.species()):
            self.sp_coeffs[i, :] = sp.thermo.coeffs

        self.nu_f = np.zeros((self.gas.n_reactions, self.gas.n_species))
        self.nu_r = np.zeros((self.gas.n_reactions, self.gas.n_species))
        self.B = np.zeros(self.Nc)
        self.Be = np.zeros(self.Ne)
        self.c = np.zeros(self.Nc)
        self.ce = np.zeros(self.Ne)
        
    def set_conditions(self, T, P, X):
        """
        Set the conditions for temperature, pressure, and mole fractions.

        Parameters:
        - T: Temperature.
        - P: Pressure.
        - X: Mole fractions.
        """
        # Set gas state
        self.gas.TPX = T, P, X

        # Prepare species dictionary for specified species
        species_dict = {}
        for species in self.specified_species:
            if species in self.gas.species_names:
                species_dict[species] = self.gas[species].X[0]
            else:
                species_dict[species] = None

        # Initialize mole fractions
        mole_fractions = np.zeros(self.gas.n_species)
        for sp in species_dict:
            index = self.gas.species_index(sp)
            mole_fractions[index] = species_dict[sp] if species_dict[sp] is not None else 0.0

        self.X0 = self.gas.X
        N_tot = 1.0
        N0 = self.X0 * N_tot

        # Compute stoichiometric coefficients for forward and reverse reactions
        for i, r in enumerate(self.gas.reactions()):
            for (key, val) in r.reactants.items():
                self.nu_f[i, self.gas.species_index(key)] = val
            for (key, val) in r.products.items():
                self.nu_r[i, self.gas.species_index(key)] = val

        self.el_names = self.gas.element_names
        self.el_numbers = np.zeros((self.gas.n_species, self.Nc))
        c = np.zeros((self.Nc))

        if self.use_element_constraints:
            # Method 1: Use element constraints
            for i, sp in enumerate(self.gas.species()):
                for (key, val) in sp.composition.items():
                    self.el_numbers[i, self.gas.element_index(key)] = val
            B = self.el_numbers
            c = B.T @ N0
            for i, sp in enumerate(self.specified_species):
                self.el_numbers[self.gas.species_index(sp), i + self.Ne] = 1
                c[i + self.Ne] = self.gas[sp].X
        else:
            # Method 2: Nc = len(constrained_species) + 1
            for i, sp in enumerate(self.specified_species):
                self.el_numbers[self.gas.species_index(sp), i + 1] = 1
                c[i + 1] = self.gas[sp].X
            self.el_numbers[:, 0] = 1
            c[0] = self.X0.sum()

        # self.gas.set_unnormalized_mole_fractions(mole_fractions)
        self.N0 = self.gas.X
        self.B = self.el_numbers
        self.Be = self.B[:, :self.Ne]
        self.c = c
        self.ce = self.c[:self.Ne]
        print("self.Be")
        print(self.Be)
        print('self.ce')
        print(self.ce)
        
        
        

    def reaction_mapping(self, phi, tau):
        """
        Evaluate the reaction mapping R(phi, tau).
        """
        # Create constant pressure reactor
        self.gas.X = phi
        print("before_reacting, gas.X=",self.gas.X)
        r = ct.IdealGasConstPressureReactor(self.gas)
        # Create simulation PSR object
        sim = ct.ReactorNet([r])
        current_time = 0
        while current_time < tau:
            sim.step()
            current_time = sim.time
        print("gas.X=",self.gas.X)
        return self.gas.X


def reaction_mapping(gas, phi, tau, num_points=10):
    """
    Evaluate the reaction mapping R(phi, tau) and record intermediate states.

    Parameters:
    - gas: Cantera gas object.
    - phi: Initial composition (array of mole fractions).
    - tau: Total simulation time.
    - num_points: Number of points to record along the reaction path.

    Returns:
    - times: List of time values (tau) at each recorded step.
    - compositions: List of composition (phi) arrays at each recorded step.
    """
    # Create constant pressure reactor
    gas.X = phi
    print("Starting reaction, initial gas.X =", gas.X)
    r = ct.IdealGasConstPressureReactor(gas)

    # Create simulation PSR object
    sim = ct.ReactorNet([r])

    # Generate logarithmically spaced time points
    log_times = np.logspace(np.log10(1e-15), np.log10(tau), num=num_points)

    # Time tracking
    times = []  # List to store time (tau)
    compositions = []  # List to store composition (phi)

    # Simulate the reaction for each target time in log_times
    for target_time in log_times:
        while sim.time < target_time:
            sim.step()
        # Record the current state
        times.append(sim.time)
        compositions.append(gas.X.copy())
        print(f"Recorded at time {sim.time:.2e}: gas.X = {gas.X}")

    return times, compositions
    
def remove_collinear_points(vertices, tolerance=1e-6):
    """
    Remove collinear points from the set of vertices.

    Parameters:
    - vertices: Array of 3D points (boundary vertices).
    - tolerance: A small value to determine collinearity.

    Returns:
    - cleaned_vertices: Array of vertices with collinear points removed.
    """
    cleaned_vertices = [vertices[0]]  # Start with the first vertex

    for i in range(1, len(vertices) - 1):
        # Vectors between consecutive points
        v1 = vertices[i] - cleaned_vertices[-1]
        v2 = vertices[i + 1] - vertices[i]

        # Calculate the cross-product magnitude to check collinearity
        cross_product = np.cross(v1, v2)
        if np.linalg.norm(cross_product) > tolerance:
            cleaned_vertices.append(vertices[i])  # Add the point if not collinear

    cleaned_vertices.append(vertices[-1])  # Add the last vertex
    return np.array(cleaned_vertices)


def find_boundary_vertices(Be, ce):
    """
    Find the boundary vertices of the feasible region for phi subject to Be * phi = ce.

    Parameters:
    - Be: Constraint matrix (m x n, where m < n).
    - ce: Constraint vector (length m).

    Returns:
    - filtered_projected_vertices: Vertices defining the boundary of the feasible region after removing collinear points.
    - boundary_vertices: Vertices on the convex hull of the feasible region.
    - projected_vertices: All vertices projected into the phi[0], phi[1], phi[2] space.
    """
    num_variables = Be.shape[1]  # Number of variables (phi dimensions)
    vertices = []

    # Solve for the vertices of the polytope by minimizing and maximizing each variable
    for i in range(num_variables):
        # Minimize phi[i]
        c = np.zeros(num_variables)
        c[i] = 1  # Objective is to minimize phi[i]
        res_min = linprog(c, A_eq=Be, b_eq=ce, bounds=(0, None))
        if res_min.success:
            vertices.append(res_min.x)

        # Maximize phi[i] (equivalent to minimizing -phi[i])
        c[i] = -1  # Objective is to maximize phi[i]
        res_max = linprog(c, A_eq=Be, b_eq=ce, bounds=(0, None))
        if res_max.success:
            vertices.append(res_max.x)

    # Remove duplicate vertices
    vertices = np.unique(vertices, axis=0)

    # Project the vertices into the phi[0], phi[1], phi[2] space
    projected_vertices = vertices[:, :3]  # Keep only the first 3 dimensions (phi[0], phi[1], phi[2])

    # Use ConvexHull to filter out redundant vertices
    try:
        hull = ConvexHull(projected_vertices)
        boundary_indices = hull.vertices  # Indices of the vertices on the convex hull
        boundary_vertices = projected_vertices[boundary_indices]  # Keep only boundary vertices
        vertices = vertices[boundary_indices]
        # Remove collinear points
        filtered_projected_vertices = remove_collinear_points(boundary_vertices)
    except Exception as e:
        print("ConvexHull failed:", e)
        boundary_vertices = projected_vertices  # Fallback to using all projected vertices
        filtered_projected_vertices = remove_collinear_points(boundary_vertices)

    return filtered_projected_vertices, projected_vertices, vertices


def plot_boundary_region(vertices):
    """
    Plot the boundary region in 3D with swapped axes:
    - X-axis: phi[1]
    - Y-axis: phi[2]
    - Z-axis: phi[0]

    Parameters:
    - vertices: Array of vertices (phi values) defining the boundary of the feasible region.
    """
    # Compute the convex hull of the vertices
    try:
        hull = ConvexHull(vertices)
    except Exception as e:
        print("ConvexHull failed:", e)
        return

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Plot the vertices with swapped axes
    ax.scatter(vertices[:, 1], vertices[:, 2], vertices[:, 0], c="r", label="Vertices")

    # Plot the boundary faces with the new axis order
    for simplex in hull.simplices:
        # Swap axes for the faces
        face = vertices[simplex]
        poly = Poly3DCollection([face[:, [1, 2, 0]]], alpha=0.3, edgecolor="k")  # Swap axes for face
        ax.add_collection3d(poly)

    # Labels and title
    ax.set_xlabel("H2")
    ax.set_ylabel("O")
    ax.set_zlabel("H2O")
    ax.set_title("Boundary Region of Feasible Space")
    ax.legend()
    plt.savefig("figs/feasible_region.png")
    plt.show()


def plot_boundary_region_with_paths(vertices, reaction_paths):
    """
    Plot the boundary region in 3D with multiple reaction paths.

    Parameters:
    - vertices: Array of vertices (phi values) defining the boundary of the feasible region.
    - reaction_paths: List of arrays, where each array contains points along a reaction path.
    """
    # Compute the convex hull of the vertices
    try:
        hull = ConvexHull(vertices)
    except Exception as e:
        print("ConvexHull failed:", e)
        return

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Plot the vertices with swapped axes
    ax.scatter(vertices[:, 1], vertices[:, 2], vertices[:, 0], c="r", label="Boundary Vertices")

    # Plot the boundary faces with the new axis order
    for simplex in hull.simplices:
        # Swap axes for the faces
        face = vertices[simplex]
        poly = Poly3DCollection([face[:, [1, 2, 0]]], alpha=0.3, edgecolor="k")  # Swap axes for face
        ax.add_collection3d(poly)

    # Plot each reaction path
    colors = ['b', 'g', 'm', 'c', 'y']  # Define colors for different paths
    for i, reaction_path in enumerate(reaction_paths):
        reaction_path = np.array(reaction_path)  # Ensure it's a NumPy array
        color = colors[i % len(colors)]  # Cycle through colors if more paths than colors
        ax.plot(
            reaction_path[:, 1], reaction_path[:, 2], reaction_path[:, 0],
            c=color, label=f"Reaction Path {i + 1}", marker="o"
        )

    # Labels and title
    ax.set_xlabel("H2")
    ax.set_ylabel("O")
    ax.set_zlabel("H2O")
    ax.set_title("Boundary Region with Reaction Paths")
    # ax.legend()
    plt.show()
    plt.savefig("figs/feasible_region_with_paths.png")



def main():
    mech = "src/Li_2003/Li_demo.cti"  # Path to the mechanism file
    Nsum = 0.1234 + 0.0411 + 0.6581
    H =  0.1234 / Nsum
    O = 0.0411 / Nsum
    N2 = 0.6581 / Nsum
    # Print the results
    print(f"Sum of moles (Nsum): {Nsum:.5f}")
    print(f"Mole fraction of H: {H:.5f}")
    print(f"Mole fraction of O: {O:.5f}")
    print(f"Mole fraction of N2: {N2:.5f}")
    T, P, X = 2500, 1 * ct.one_atm, "H:0.15, O:0.05, N2:0.8"
    main_species_names = ["O", "H2", "N2"]
    # Initialize Cantera gas object and equilibrate
    gas = ct.Solution(mech)
    gas.TPX = T, P, X
    gas.equilibrate('HP')
    Tini, Pini = gas.TP
    Xini = gas.X
    print("T =", Tini)
    print("Xini =", Xini)
    # Initialize PIC object
    pic = PIC(mech, main_species_names, use_element_constraints=True)
    input_species_values = [gas[name].X for name in gas.species_names]
    Xini = dict(zip(gas.species_names, input_species_values))
    print("Xini =", Xini)

    pic.set_conditions(Tini, Pini, Xini)
    # Find the feasible region
    filtered_projected_vertices, projected_vertices, vertices = find_boundary_vertices(pic.Be.T, pic.ce)
    # plot_boundary_region(filtered_projected_vertices)
    # Generate reaction path
    input_phi = (vertices[1] + vertices[0])/2
    tau = 1e-4
    gas.TP = Tini, Pini
    #########MOIDIFY STARTING POINTS############
    # 1. 
    # starting_points = vertices  
    # 2. 
    # num_points = 10
    # # Generate points along the line between vertices[0] and vertices[1]
    # line_1 = [
    #     vertices[0] + t * (vertices[1] - vertices[0]) for t in np.linspace(0, 1, num_points)
    # ]

    # # Generate points along the line between vertices[1] and vertices[2]
    # line_2 = [
    #     vertices[1] + t * (vertices[2] - vertices[1]) for t in np.linspace(0, 1, num_points)
    # ]
    # line_3 = [
    #     vertices[2] + t * (vertices[3] - vertices[2]) for t in np.linspace(0, 1, num_points)
    # ]
    # line_4 = [
    #     vertices[3] + t * (vertices[4] - vertices[3]) for t in np.linspace(0, 1, num_points)
    # ]
    
    # Combine all starting points
    # starting_points = line_1 + line_2 + line_3 + line_4
    # 3. 
    # Generate points on a quadrilateral surface
    num_points = 10  # Number of points in each direction
    quadrilateral_surface = []
    for i in range(num_points):
        for j in range(num_points):
            u = i / (num_points - 1)
            v = j / (num_points - 1)
            point = (
                (1 - u) * (1 - v) * vertices[0] +  # Bottom-left corner
                u * (1 - v) * vertices[1] +        # Bottom-right corner
                u * v * vertices[6] +             # Top-right corner
                (1 - u) * v * vertices[7]         # Top-left corner
            )
            quadrilateral_surface.append(point) 

    # Convert to NumPy array
    quadrilateral_surface = np.array(quadrilateral_surface)
    starting_points = quadrilateral_surface  # Or quadrilateral_surface
    ################END of MODIFYING INPUT POINTS###################################
    # Time for the reaction path
    tau = 1e-3
    # Store all reaction paths
    reaction_paths = []  # List to hold all reaction paths

    # Loop through starting points and calculate reaction paths
    for i, input_phi in enumerate(starting_points):
        print(f"Calculating reaction path {i + 1} starting from {input_phi}")
        times, compositions = reaction_mapping(gas, input_phi, tau, num_points=20)
        phi_values = np.array(compositions)
        reaction_paths.append(phi_values)

    # Plot the boundary region with all reaction paths
    plot_boundary_region_with_paths(filtered_projected_vertices, reaction_paths)

if __name__ == '__main__':
    main()
    
    
    