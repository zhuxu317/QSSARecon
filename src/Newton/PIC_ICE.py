import os, sys, time
sys.path.append(os.path.abspath('src'))
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.linalg import lstsq
from functools import lru_cache
import pandas as pd
# from run_1D_CEQ import run_CEQ_core
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.optimize import minimize
import cvxpy as cp
from run_1D_CEQ import run_CEQ_core


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
        print(f"species names: {self.gas.species_names}")
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

    gas.X = phi 
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

    return times, compositions
    
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



def objective_function(phi_0, gas, B, r, tau):
    """
    Compute the objective function J(phi_0) = || B^T R(phi_0, tau) - r ||^2.

    Parameters:
    - phi_0: Initial composition (array of mole fractions).
    - gas: Cantera gas object.
    - B: Matrix B used in the objective function.
    - r: Target vector for the optimization.
    - tau: Simulation time.

    Returns:
    - J: Value of the objective function.
    """
    _, compositions = reaction_mapping(gas, phi_0, tau, num_points=20)  # Simulate reaction mapping
    R_phi_tau = compositions[-1]  # Get the final composition after reaction
    J = np.linalg.norm(np.dot(B.T, R_phi_tau) - r)**2  # Compute the squared norm
    return J


def optimize_phi_0(gas, pic, B, r, tau, boundary_vertices):
    """
    Perform optimization to find the optimal phi_0 on the boundary surface.

    Parameters:
    - gas: Cantera gas object.
    - B: Matrix B used in the objective function.
    - r: Target vector for the optimization.
    - tau: Simulation time.
    - boundary_vertices: Boundary surface points (starting points for optimization).

    Returns:
    - result: Result of the optimization process.
    - optimization_trajectory: List of phi_0 values visited during optimization.
    """
    # Initial guess: Use the centroid of the boundary vertices
    phi_0_initial = np.mean(boundary_vertices, axis=0)

    # Store optimization trajectory
    optimization_trajectory = []

    # Callback function to store intermediate phi_0 values
    def callback(phi_0):
        optimization_trajectory.append(phi_0)

    # Define constraints: phi_0 must be non-negative and sum to 1 (mole fractions)
    constraints = [
        {'type': 'eq', 'fun': lambda phi: np.dot(pic.Be.T, phi) - pic.ce},  # Linear equality constraint
        {'type': 'ineq', 'fun': lambda phi: phi}  # Non-negativity constraint
    ]

    # Perform optimization
    result = minimize(
        objective_function,
        phi_0_initial,
        args=(gas, B, r, tau),
        method='SLSQP',  # Sequential Least Squares Programming
        constraints=constraints,
        callback=callback,
        options={'disp': True}
    )

    return result, optimization_trajectory


def polytope_vertices(Be, ce):
    import cdd
    import pprint
    """
    Compute the vertices of a polytope defined by Be * phi = ce and 0 <= phi <= 1.

    Parameters:
    - Be: Constraint matrix (m x n, where m is the number of equality constraints and n is the number of variables).
    - ce: Constraint vector (length m).

    Returns:
    - vertices: Array of vertices of the polytope (n_vertices x n).
    """
    num_variables = Be.shape[1]
    print("Be")
    print(Be)
    # Inequalities: 0 <= phi <= 1
    G = np.vstack([-np.eye(num_variables), np.eye(num_variables)])  # Identity matrix for phi >= 0 and -I for phi <= 1
    h = np.hstack([np.ones(num_variables), np.zeros(num_variables)])  # Upper bound (1) and lower bound (0)

    #  0 <= 0 + x
    #  0 <= 1  - x 
    A = -Be
    b = ce
    H = np.hstack([h[:, None], G])  # Inequality constraints
    E = np.hstack([b[:, None], A])  # Equality constraints
    H_rep = np.vstack([H, E])  # Combine inequalities and equalities
    mat = cdd.Matrix(H_rep, number_type="float")
    mat.rep_type = cdd.RepType.INEQUALITY
    poly = cdd.Polyhedron(mat)
    vertices = np.array(poly.get_generators())
    print("vertices")
    print(vertices)
    return vertices[:, 1:]  # Skip the first column (homogeneous coordinate)


def polytope_vertices_and_facets(Be, ce):
    import cdd
    """
    Compute the vertices of a polytope and extract its facet equations.

    Parameters:
    - Be: Constraint matrix (m x n, where m is the number of equality constraints and n is the number of variables).
    - ce: Constraint vector (length m).

    Returns:
    - vertices: Array of vertices of the polytope (n_vertices x n).
    - facets: Array of facet equations of the polytope (n_facets x (n+1)).
    """
    num_variables = Be.shape[1]

    # Inequalities: 0 <= phi <= 1
    G = np.vstack([-np.eye(num_variables), np.eye(num_variables)])  # Identity matrix for phi >= 0 and -I for phi <= 1
    h = np.hstack([np.ones(num_variables), np.zeros(num_variables)])  # Upper bound (1) and lower bound (0)

    # Equalities: Be * phi = ce
    A = -Be
    b = ce

    # Construct H-representation
    H = np.hstack([h[:, None], G])  # Inequality constraints
    E = np.hstack([b[:, None], A])  # Equality constraints
    H_rep = np.vstack([H, E])  # Combine inequalities and equalities

    # Convert to cddlib representation
    mat = cdd.Matrix(H_rep, number_type="float")
    mat.rep_type = cdd.RepType.INEQUALITY

    # Create the polyhedron
    poly = cdd.Polyhedron(mat)

    # Extract vertices
    vertices = np.array(poly.get_generators())[:, 1:]  # Skip the first column (homogeneous coordinate)

    # Extract facets (inequalities)
    facets = np.array(poly.get_inequalities())

    return vertices, facets


def plot_point_cloud_from_facets(vertices, facets):
    """
    Generate a point cloud from polytope facets and plot it.

    Parameters:
    - vertices: Array of vertices of the polytope (n_vertices x n).
    - facets: Array of facet equations of the polytope (n_facets x (n+1)).
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the vertices
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c="b", label="Vertices")

    # Generate a random point cloud for each facet
    point_cloud = []
    for facet in facets:
        normal = facet[1:]  # The normal vector of the facet
        offset = facet[0]  # The offset (constant term) of the facet equation

        # Generate random points on the facet plane
        for _ in range(100):  # Number of points per facet
            random_point = np.random.rand(3)  # Random point in 3D
            # Project the random point onto the facet plane
            distance = (np.dot(normal, random_point) + offset) / np.linalg.norm(normal)
            projected_point = random_point - distance * normal
            point_cloud.append(projected_point)

    point_cloud = np.array(point_cloud)

    # Plot the point cloud
    ax.scatter(point_cloud[:, 0], point_cloud[:, 1], point_cloud[:, 2], c="r", s=1, label="Point Cloud")

    # Label axes
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Point Cloud from Polytope Facets")
    plt.legend()
    plt.show()
    
    
def plot_feasible_region(feasible_points):
    """
    Plot the feasible region in 3D using sampled points.

    Parameters:
    - feasible_points: Array of points in the feasible region (n x 3).
    """
    # Select the first 3 columns and reorder them to [1, 2, 0]
    feasible_points = feasible_points[:, :3]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the feasible points
    ax.scatter(feasible_points[:, 0], feasible_points[:, 1], feasible_points[:, 2], c="b", label="Feasible Points")

    # Compute the convex hull for the feasible region
    try:
        hull = ConvexHull(feasible_points)
        # Plot the hull as a collection of triangles
        for simplex in hull.simplices:
            triangle = feasible_points[simplex]
            poly = Poly3DCollection([triangle], alpha=0.4, edgecolor="k")
            ax.add_collection3d(poly)
    except Exception as e:
        print("ConvexHull plotting failed:", e)

    # Label axes
    ax.set_xlabel("H2")
    ax.set_ylabel("O")
    ax.set_zlabel("H2O")
    ax.set_title("Feasible Region")

    plt.legend()
    plt.savefig("figs/feasible_region_cdd.png")
    plt.show()

def plot_feasible_region_with_Xini(feasible_points, Xini):
    """
    Plot the feasible region in 3D using sampled points, scatter Xini,
    and denote the boundary points of the convex hull.

    Parameters:
    - feasible_points: Array of points in the feasible region (n x 3).
    - Xini: Initial point (1D array of length 3 or more).
    """
    # Select the first 3 columns
    feasible_points = feasible_points[:, :3]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the feasible points
    ax.scatter(feasible_points[:, 0], feasible_points[:, 1], feasible_points[:, 2], c="b", label="Feasible Points")

    # Plot Xini as a red scatter point
    ax.scatter(Xini[0], Xini[1], Xini[2], c="r", label="Xini", s=100, marker="o")

    # Compute the convex hull for the feasible region
    try:
        hull = ConvexHull(feasible_points)
        # Plot the hull as a collection of triangles
        for simplex in hull.simplices:
            triangle = feasible_points[simplex]
            poly = Poly3DCollection([triangle], alpha=0.4, edgecolor="k")
            ax.add_collection3d(poly)

        # Identify boundary points (unique vertices on the convex hull)
        boundary_points_indices = np.unique(hull.simplices)
        boundary_points = feasible_points[boundary_points_indices]

        # Annotate the boundary points
        for i, point in zip(boundary_points_indices, boundary_points):
            ax.text(point[0], point[1], point[2], f"{i}", color="black", fontsize=8)
    except Exception as e:
        print("ConvexHull plotting failed:", e)

    # Label axes
    ax.set_xlabel("H2")
    ax.set_ylabel("O")
    ax.set_zlabel("H2O")
    ax.set_title("Feasible Region with Xini and Boundary Points")

    plt.legend()
    plt.savefig("figs/feasible_region_with_Xini_and_boundary_points.png")
    plt.show()

def sample_points_on_surface(vertices, num_points=100):
    from itertools import combinations
    """
    Generate uniformly distributed points on the surface defined by the given vertices.

    Parameters:
    - vertices: Array of vertices defining the surface (n x 3, where n = 2, 3, or 4).
    - num_points: Number of points to sample.

    Returns:
    - points: Array of sampled points (m x 3), where m is the actual number of points.
    """
    num_vertices = len(vertices)
    points = []

    if num_vertices == 2:
        # Case 1: Line (2 vertices)
        # Sample points uniformly along the line
        for t in np.linspace(0, 1, num_points):
            point = (1 - t) * vertices[0] + t * vertices[1]
            points.append(point)

    elif num_vertices == 3:
        # Case 2: Triangle (3 vertices)
        # Uniformly sample using barycentric coordinates
        grid_size = int(np.sqrt(num_points))  # Approximate grid size
        for i in range(grid_size + 1):
            for j in range(grid_size + 1 - i):
                u = i / grid_size
                v = j / grid_size
                w = 1 - u - v
                point = u * vertices[0] + v * vertices[1] + w * vertices[2]
                points.append(point)

    elif num_vertices == 4:
        # Case 3: Quadrilateral (4 vertices)
        # Split the quadrilateral into two triangles and sample points uniformly on each triangle
        triangle1 = vertices[:3]  # First triangle (v0, v1, v2)
        triangle2 = [vertices[0], vertices[2], vertices[3]]  # Second triangle (v0, v2, v3)

        # Distribute points equally between the two triangles
        num_points_triangle = num_points // 2

        # Sample points on the first triangle
        points += list(sample_points_on_surface(triangle1, num_points_triangle))

        # Sample points on the second triangle
        points += list(sample_points_on_surface(triangle2, num_points - num_points_triangle))

    else:
        raise ValueError("The surface must be defined by 2, 3, or 4 vertices.")

    return np.array(points)

def plot_feasible_region_with_Xini_boundary(feasible_points, Xini, starting_points):
    """
    Plot the feasible region in 3D using sampled points, scatter Xini,
    and plot boundary scatter points.

    Parameters:
    - feasible_points: Array of points in the feasible region (n x 3).
    - Xini: Initial point (1D array of length 3 or more).
    - starting_points: Array of sampled points (k x 3) to plot as boundary scatter points.
    """
    # Select the first 3 columns
    feasible_points = feasible_points[:, :3]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the feasible points
    ax.scatter(feasible_points[:, 0], feasible_points[:, 1], feasible_points[:, 2], c="b", label="Feasible Points")

    # Plot Xini as a red scatter point
    ax.scatter(Xini[0], Xini[1], Xini[2], c="r", label="Xini", s=100, marker="o")

    # Compute the convex hull for the feasible region
    try:
        hull = ConvexHull(feasible_points)
        # Plot the hull as a collection of triangles
        for simplex in hull.simplices:
            triangle = feasible_points[simplex]
            poly = Poly3DCollection([triangle], alpha=0.4, edgecolor="k")
            ax.add_collection3d(poly)

        # Scatter-plot the starting points
        ax.scatter(
            starting_points[:, 0], starting_points[:, 1], starting_points[:, 2],
            c="g", label="Boundary Scatter Points", s=50, marker="^"
        )
    except Exception as e:
        print("ConvexHull plotting failed:", e)

    # Label axes
    ax.set_xlabel("H2")
    ax.set_ylabel("O")
    ax.set_zlabel("H2O")
    ax.set_title("Feasible Region with Xini and Boundary Scatter Points")

    plt.legend()
    plt.savefig("figs/feasible_region_with_Xini_boundary.png")
    plt.show()


def plot_feasible_region_with_Xini_boundary_trajectory(feasible_points, Xini, X_qssa, starting_points, reaction_paths, elev=30, azim=45):
    """
    Plot the feasible region in 3D, with sampled boundary points, Xini,
    and reaction paths as trajectories.

    Parameters:
    - feasible_points: Array of points in the feasible region (n x 3).
    - Xini: Initial point (1D array of length 3 or more).
    - starting_points: Array of sampled points (k x 3) to plot as boundary scatter points.
    - reaction_paths: List of reaction paths, where each path is an array of points (m x 3).
    - elev: Elevation angle in degrees (default=30).
    - azim: Azimuthal angle in degrees (default=45).
    """
    # Select the first 3 columns
    feasible_points = feasible_points[:, :3]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the feasible points
    ax.scatter(feasible_points[:, 0], feasible_points[:, 1], feasible_points[:, 2], c="b", label="Feasible Points")

    ax.scatter(Xini[0], Xini[1], Xini[2], c="k", label="Xini", s=50, marker="o")
    
    ax.scatter(X_qssa[0], X_qssa[1], X_qssa[2], c="r", label="Xqssa", s=50, marker="^")

    # Compute the convex hull for the feasible region
    try:
        hull = ConvexHull(feasible_points)
        # Plot the hull as a collection of triangles
        for simplex in hull.simplices:
            triangle = feasible_points[simplex]
            poly = Poly3DCollection([triangle], alpha=0.05, edgecolor="k")
            ax.add_collection3d(poly)

        # Scatter-plot the starting points
        ax.scatter(
            starting_points[:, 0], starting_points[:, 1], starting_points[:, 2],
            c="g", label="Boundary Scatter Points", s=10, marker="o"
        )

        # Identify boundary points (unique vertices on the convex hull)
        boundary_points_indices = np.unique(hull.simplices)
        boundary_points = feasible_points[boundary_points_indices]

        # Annotate the boundary points
        for i, point in zip(boundary_points_indices, boundary_points):
            ax.text(point[0], point[1], point[2], f"{i}", color="black", fontsize=8)
            
        # Plot the reaction paths as trajectories
        for path in reaction_paths:
            path = np.array(path)
            ax.plot(
                path[:, 0], path[:, 1], path[:, 2],
                label="Reaction Path", linestyle="-", linewidth=1, alpha=0.8
            )

    except Exception as e:
        print("ConvexHull plotting failed:", e)

    # Set view angle
    ax.view_init(elev=elev, azim=azim)

    # Reduce the number of ticks to 3 for each axis
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    z_min, z_max = ax.get_zlim()
    ax.set_xticks(np.linspace(0, 0.09, 3))
    ax.set_yticks(np.linspace(0, 0.06, 3))
    ax.set_zticks(np.linspace(0, 0.06, 3))

    # Label axes
    ax.set_xlabel("H2")
    ax.set_ylabel("O")
    ax.set_zlabel("H2O")
    ax.set_title("Feasible Region with Xini, Boundary, and Reaction Paths")
    plt.savefig("figs/feasible_region_with_Xini_boundary_trajectory.png")
    plt.show()



def main():
    mech = "src/Li_2003/Li_demo.cti"  # Path to the mechanism file
    Nsum = 0.1234 + 0.0411 + 0.6581
    H =  0.1234 / Nsum
    O = 0.0411 / Nsum
    N2 = 0.6581 / Nsum
    print(f"Sum of moles (Nsum): {Nsum:.5f}")
    print(f"Mole fraction of H: {H:.5f}")
    print(f"Mole fraction of O: {O:.5f}")
    print(f"Mole fraction of N2: {N2:.5f}")
    T, P, X = 500, 1 * ct.one_atm, "H:0.15, O:0.05, N2:0.8"
    main_species_names = ["O", "H2", "N2"]
    gas = ct.Solution(mech)
    gas.TPX = T, P, X
    gas.equilibrate('HP')
    Tini, Pini = gas.TP
    Xini = gas.X
    print("T =", Tini)
    print("Xini =", Xini)
    
    pic = PIC(mech, main_species_names, use_element_constraints=True)
    pic.set_conditions(Tini, Pini, Xini)

    # Find the feasible region
    print("pic.Be.T")
    print(pic.Be.T)
    print(pic.Be.T[:2,:])
    print("pic.ce")
    print(pic.ce)
    # Sample points in the feasible region
    
    feasible_points = polytope_vertices(pic.Be.T[:2, :-1], pic.ce[:2])

    # plot_feasible_region(feasible_points)
    # plot_feasible_region_with_Xini(feasible_points, Xini)

    # boundary_vertex_indices = [1,0,4,7] #Back
    # boundary_vertex_indices = [1,2,6,7] #Front

    # boundary_vertex_indices = [4,6,7] # left
    # boundary_vertex_indices = [1,0,2] # right
    boundary_vertex_indices = [2,0,4,6] # bottom
    
    
    

    boundary_vertices = feasible_points[boundary_vertex_indices]
    # Step 1: Sample points on the surface
    starting_points = sample_points_on_surface(boundary_vertices, num_points=300)
    # Step 2: Plot the feasible region with Xini and the sampled boundary points
    # plot_feasible_region_with_Xini_boundary(feasible_points, Xini, starting_points)

    tau = 1e2
    N_N2 = 1.70681251/2
    # result, optimization_trajectory = optimize_phi_0(gas, pic, pic.Be.T, pic.c, tau, filtered_projected_vertices)
    reaction_paths = []
    N_new = []
    for input_N in starting_points:
        N_new = np.append(input_N, N_N2)
        N_bar = (input_N.sum() + N_N2)
        phi_new = N_new / N_bar
        _, compositions = reaction_mapping(gas, phi_new, tau, num_points=50)
        reaction_paths.append(np.array(compositions))
    # Plot the QSSA point 
    gas.TPX = Tini, Pini, Xini
    input_species_values = [gas[name].X.item() for name in main_species_names]  # Convert array to scalar
    print("input_species_values")
    print(input_species_values)
    Xini_CEQ = dict(zip(main_species_names, input_species_values))
    gas = run_CEQ_core(gas, main_species_names, Xini_CEQ, Tini, Pini, dt=1e-1, time_end=1e0)
    X_qssa = gas.X  # Mole fractions from QSSA
    
    # Plot the boundary region with reaction paths and optimization trajectory
    plot_feasible_region_with_Xini_boundary_trajectory(feasible_points, Xini,X_qssa, starting_points, reaction_paths)

if __name__ == '__main__':
    main()
    
    
    